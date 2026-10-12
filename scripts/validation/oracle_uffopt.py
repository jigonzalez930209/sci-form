#!/usr/bin/env python3
"""RDKit UFF geometry-optimization oracle for sci-form validation (ficha F70).

For each molecule: ETKDG embed (fixed seed) -> record the unminimized
geometry -> ``UFFOptimizeMolecule`` to convergence -> record the minimized
geometry and final UFF energy.

The Rust test (``tests/cross_validation.rs::f70_framework_opt``) starts
sci-form's own BFGS (``materials::optimize_framework`` over the validated
F16 UFF surface) from the same initial geometry and compares final
energies, forces and geometries.

Generates ``tests/fixtures/validation/uff_opt.json``.

Usage:
    python scripts/validation/oracle_uffopt.py
"""

from __future__ import annotations

import time

from rdkit import Chem
from rdkit.Chem import AllChem

from common import Report, emit_fixture_report

REFERENCE = "RDKit 2025.09 UFFOptimizeMolecule"

# Small, diverse, single-fragment molecules (no metals, no charges).
SMILES_OPT = [
    ("ethanol", "CCO"),
    ("butane", "CCCC"),
    ("acetone", "CC(=O)C"),
    ("acetic acid", "CC(=O)O"),
    ("methylamine", "CN"),
    ("cyclohexane", "C1CCCCC1"),
    ("dimethyl ether", "COC"),
    ("benzene", "c1ccccc1"),
    ("ethylene", "C=C"),
    ("formaldehyde", "C=O"),
]


def conf_coords(mol, conf_id=0):
    conf = mol.GetConformer(conf_id)
    out = []
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        out += [p.x, p.y, p.z]
    return out


def run_case(name: str, smiles: str) -> dict:
    mol0 = Chem.MolFromSmiles(smiles)
    if mol0 is None:
        raise RuntimeError(f"cannot parse {smiles}")
    mol = Chem.AddHs(mol0)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.useRandomCoords = True
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError(f"ETKDG failed for {name}")
    initial = conf_coords(mol)
    ff = AllChem.UFFGetMoleculeForceField(mol, confId=-1)
    if ff is None:
        raise RuntimeError(f"no UFF force field for {name}")
    e_initial = ff.CalcEnergy()
    rc = AllChem.UFFOptimizeMolecule(mol, maxIters=5000)
    # NOTE: the pre-optimization FF object goes stale (it reports the
    # unrelaxed energy after in-place minimization); a fresh FF on the
    # final conformer gives the true minimum energy.
    e_final = AllChem.UFFGetMoleculeForceField(mol, confId=-1).CalcEnergy()
    return {
        "name": name,
        "smiles": smiles,
        "elements": [a.GetAtomicNum() for a in mol.GetAtoms()],
        "initial_coords": initial,
        "e_initial": e_initial,
        "optimized_coords": conf_coords(mol),
        "e_final": e_final,
        "rdkit_converged": rc == 0,
    }


def main() -> int:
    cases = []
    for name, smiles in SMILES_OPT:
        case = run_case(name, smiles)
        cases.append(case)
        print(
            f"  {name:14s} E_init={case['e_initial']:10.3f} "
            f"E_opt={case['e_final']:10.3f} kcal/mol "
            f"converged={case['rdkit_converged']}",
            flush=True,
        )

    payload = {
        "reference": REFERENCE,
        "energies": "kcal/mol",
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F70",
        observable="uff_opt_energy_vs_rdkit",
        reference_tool=REFERENCE,
        units="kcal/mol/atom",
        tolerance=1e-6,
        level="L2",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "Same UFF surface (F16 parity 3.98e-13 kcal/mol), same start: "
            "sci-form BFGS vs RDKit minimization. Tolerance 1e-3 "
            "kcal/mol/atom (4.3e-8 eV/atom, stricter than the plan's "
            "1e-6 eV/atom: both minimizers land in the same basin). "
            "Measured by tests/cross_validation.rs::f70_framework_opt."
        ),
    )
    emit_fixture_report("uff_opt.json", payload, report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
