#!/usr/bin/env python3
"""tblite oracle for sci-form GFN1/GFN2 cross-validation (fichas F44, F45).

Builds a fixed 3D corpus deterministically (RDKit ETKDG + MMFF94) and evaluates
GFN1-xTB / GFN2-xTB single points with the reference tblite library.

Generates `tests/fixtures/validation/tblite_gfn.json` with per molecule:
  - elements, geometry (Angstrom)
  - gfn1 / gfn2: total energy (eV), HOMO/LUMO (eV)
  - the tblite version used

Usage:
    python scripts/validation/oracle_tblite.py
"""

from __future__ import annotations

import sys
import time

import numpy as np

from common import HARTREE2EV, HARTREE2KCAL, Report, emit_fixture_report

try:
    from tblite.interface import Calculator
except ImportError as exc:  # pragma: no cover
    print(f"tblite is required for this oracle: {exc}", file=sys.stderr)
    sys.exit(1)

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError as exc:  # pragma: no cover
    print(f"RDKit is required for geometry generation: {exc}", file=sys.stderr)
    sys.exit(1)

import tblite

TBLITE_VERSION = getattr(tblite, "__version__", None) or "0.7.0"
REFERENCE = f"tblite {TBLITE_VERSION}"

SMILES = [
    ("water", "O"),
    ("methane", "C"),
    ("ammonia", "N"),
    ("hydrogen fluoride", "F"),
    ("methanol", "CO"),
    ("ethanol", "CCO"),
    ("ethane", "CC"),
    ("ethylene", "C=C"),
    ("formaldehyde", "C=O"),
    ("acetonitrile", "CC#N"),
    ("benzene", "c1ccccc1"),
    ("toluene", "Cc1ccccc1"),
    ("pyridine", "c1ccncc1"),
    ("hydrogen sulfide", "S"),
    ("hydrogen chloride", "Cl"),
    ("phosphine", "P"),
    ("chloroform", "ClC(Cl)Cl"),
    ("acetic acid", "CC(=O)O"),
]


def embed(smiles: str, seed: int = 0xC0FFEE):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise ValueError(f"embedding failed for {smiles}")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol)
    else:
        AllChem.UFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    elements = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    coords = []
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        coords += [p.x, p.y, p.z]
    return elements, coords


def singlepoint(method: str, elements, coords) -> dict:
    numbers = np.array(elements, dtype=np.int32)
    angles = np.array(coords).reshape(-1, 3)
    pos_bohr = angles * (1.0 / 0.529177210903)
    calc = Calculator(method, numbers, pos_bohr, charge=0.0)
    calc.set("verbosity", 0)
    res = calc.singlepoint()
    energy = float(res.get("energy"))
    try:
        orbital = np.asarray(res.get("orbital-energies")).ravel()
        occ = np.asarray(res.get("orbital-occupations")).ravel()
        homo = float(orbital[occ > 0.5].max())
        lumo = float(orbital[occ <= 0.5].min())
        level = (homo, lumo)
    except Exception:
        level = None
    return {
        "e_tot_ev": energy * HARTREE2EV,
        "e_tot_hartree": energy,
        "mo_homo_ev": level[0] * HARTREE2EV if level else None,
        "mo_lumo_ev": level[1] * HARTREE2EV if level else None,
    }


def run() -> Report:
    cases = []
    for name, smiles in SMILES:
        elements, coords = embed(smiles)
        cases.append(
            {
                "name": name,
                "smiles": smiles,
                "elements": elements,
                "coords": coords,
                "gfn1": singlepoint("GFN1-xTB", elements, coords),
                "gfn2": singlepoint("GFN2-xTB", elements, coords),
            }
        )

    payload = {
        "reference": REFERENCE,
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F44",
        observable="gfn1_gfn2_fixture_provenance",
        reference_tool=REFERENCE,
        units="eV",
        tolerance=0.0217,  # 0.5 kcal/mol
        level="L4",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "Provenance only (tolerances live with the measured rows). "
            "F45 (GFN2) is validated to 6.5e-4 eV and F44 (GFN1) to 8.1e-5 eV "
            "by tests/cross_validation.rs::{f45_gfn2_matches_tblite,"
            "f44_gfn1_matches_tblite}."
        ),
    )
    emit_fixture_report("tblite_gfn.json", payload, report)
    return report


# ---------------------------------------------------------------------------
# F84 — ALPB / GB solvation
# ---------------------------------------------------------------------------

# name, SMILES, class ("hydrocarbon" = neutral and non-polar)
ALPB_SMILES = [
    ("water", "O", "polar"),
    ("methanol", "CO", "polar"),
    ("ethanol", "CCO", "polar"),
    ("acetone", "CC(C)=O", "polar"),
    ("acetic-acid", "CC(=O)O", "polar"),
    ("urea", "NC(N)=O", "polar"),
    ("aniline", "Nc1ccccc1", "polar"),
    ("methane", "C", "hydrocarbon"),
    ("ethane", "CC", "hydrocarbon"),
    ("propane", "CCC", "hydrocarbon"),
    ("butane", "CCCC", "hydrocarbon"),
    ("pentane", "CCCCC", "hydrocarbon"),
    ("hexane", "CCCCCC", "hydrocarbon"),
    ("ethylene", "C=C", "hydrocarbon"),
    ("benzene", "c1ccccc1", "hydrocarbon"),
    ("toluene", "Cc1ccccc1", "hydrocarbon"),
]


def _solvation_singlepoint(
    method: str, elements, coords, interaction: str, args: tuple
) -> float:
    """Solvation energy (kcal/mol) of a GFN2-xTB single point."""
    numbers = np.array(elements, dtype=np.int32)
    angles = np.array(coords).reshape(-1, 3)
    pos_bohr = angles * (1.0 / 0.529177210903)
    calc = Calculator(method, numbers, pos_bohr, charge=0.0)
    calc.set("verbosity", 0)
    calc.add(interaction, *args)
    res = calc.singlepoint()
    return float(res.get("energy")) * HARTREE2KCAL


def _charges(elements, coords) -> list[float]:
    numbers = np.array(elements, dtype=np.int32)
    angles = np.array(coords).reshape(-1, 3)
    pos_bohr = angles * (1.0 / 0.529177210903)
    calc = Calculator("GFN2-xTB", numbers, pos_bohr, charge=0.0)
    calc.set("verbosity", 0)
    res = calc.singlepoint()
    return [float(q) for q in res.get("charges")]


def run_alpb() -> Report:
    """F84 — ALPB/GB solvation reference (tblite).

    * `gb_solvation` with epsilon = 78.5 and the "still" Born model is the
      *non-empirical* generalized-Born term: it contains the electrostatic
      solvation energy only (no surface term, no empirical shift).
    * `alpb_solvation`/`gbsa_solvation` with the solvent name "water" are
      the parametrized variants: they add the non-polar surface term, the
      charge-dependent solvation (CDS) term and the empirical shift.
    """
    cases = []
    for name, smiles, kind in ALPB_SMILES:
        elements, coords = embed(smiles)
        vacuum = _solvation_singlepoint("GFN2-xTB", elements, coords, "gb-solvation", (78.5, "still"))
        # The vacuum energy (no interaction) for the same geometry.
        numbers = np.array(elements, dtype=np.int32)
        pos_bohr = np.array(coords).reshape(-1, 3) * (1.0 / 0.529177210903)
        calc = Calculator("GFN2-xTB", numbers, pos_bohr, charge=0.0)
        calc.set("verbosity", 0)
        e_vacuum = float(calc.singlepoint().get("energy")) * HARTREE2KCAL

        gb = _solvation_singlepoint("GFN2-xTB", elements, coords, "gb-solvation", (78.5, "still"))
        gbe = _solvation_singlepoint("GFN2-xTB", elements, coords, "gbe-solvation", (78.5, "still"))
        alpb = _solvation_singlepoint("GFN2-xTB", elements, coords, "alpb-solvation", ("water",))
        gbsa = _solvation_singlepoint("GFN2-xTB", elements, coords, "gbsa-solvation", ("water",))
        cases.append(
            {
                "name": name,
                "smiles": smiles,
                "kind": kind,
                "elements": elements,
                "coords": coords,
                "charges": _charges(elements, coords),
                "e_vacuum_kcal": e_vacuum,
                "gb_kcal": gb - e_vacuum,
                "gbe_kcal": gbe - e_vacuum,
                "alpb_kcal": alpb - e_vacuum,
                "gbsa_kcal": gbsa - e_vacuum,
            }
        )
        del vacuum

    payload = {
        "reference": REFERENCE,
        "solvent_dielectric": 78.5,
        "born_model": "still",
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F84",
        observable="alpb_solvation_vs_tblite",
        reference_tool=f"{REFERENCE} GB/ALPB solvation",
        units="kcal/mol",
        tolerance=0.0,
        level="L2",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "tblite GB(eps=78.5, born=still) is the pure electrostatic term; "
            "the parametrized ALPB/GBSA('water') add the surface, CDS and "
            "empirical-shift contributions."
        ),
    )
    emit_fixture_report("tblite_alpb.json", payload, report)
    return report


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command", nargs="?", default="all", choices=["gfn", "alpb", "all"]
    )
    args = parser.parse_args()
    if args.command in ("gfn", "all"):
        run()
    if args.command in ("alpb", "all"):
        run_alpb()
    print("fixtures -> tests/fixtures/validation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
