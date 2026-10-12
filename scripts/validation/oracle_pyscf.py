#!/usr/bin/env python3
"""PySCF oracle for sci-form SCF cross-validation (fichas F36, F37, F23, F38, F39).

Generates fixtures under `tests/fixtures/validation/`:

  - `pyscf_scf.json` (F36/F37/F23): one-electron matrices S, T, V,
    two-electron integrals (μν|λσ) for a small subset, RHF/STO-3G total
    energy, orbital energies and Mulliken charges.
  - `pyscf_uhf.json` (F38): UHF/ROHF energies and <S^2> for open shells.
  - `pyscf_ao2mo.json` (F39): MO coefficients and AO/MO ERIs for the
    4-index AO -> MO transformation.

The Rust side (`tests/cross_validation.rs::{f36_scf, f37_eri, f38_uhf_rohf,
f39_ao2mo}`) recomputes the same quantities with `sci_form::{scf, hf}` and
asserts the tolerances from the validation plan.

Usage:
    python scripts/validation/oracle_pyscf.py [scf|uhf|ao2mo|all]
"""

from __future__ import annotations

import argparse
import sys
import time

import numpy as np

from common import DEBYE_EANG, HARTREE2EV, Report, emit_fixture_report

try:
    from pyscf import gto, scf as pyscf_scf
except ImportError as exc:  # pragma: no cover
    print(f"PySCF is required for this oracle: {exc}", file=sys.stderr)
    sys.exit(1)

PYSCF_VERSION = getattr(__import__("pyscf"), "__version__", "unknown")
REFERENCE = f"PySCF {PYSCF_VERSION}"

# Geometries in Angstrom (standard/experimental values; the fixture carries the
# exact coordinates so both sides evaluate the same system).
MOLECULES = [
    {
        "name": "H2",
        "elements": [1, 1],
        "coords": [0.0, 0.0, 0.0, 0.0, 0.0, 0.74],
        "eri": True,
    },
    {
        "name": "H2O",
        "elements": [8, 1, 1],
        "coords": [
            0.0,
            0.0,
            0.1173,
            0.0,
            0.7572,
            -0.4692,
            0.0,
            -0.7572,
            -0.4692,
        ],
        "eri": True,
    },
    {
        "name": "CH4",
        "elements": [6, 1, 1, 1, 1],
        "coords": [
            0.0,
            0.0,
            0.0,
            0.6276,
            0.6276,
            0.6276,
            -0.6276,
            -0.6276,
            0.6276,
            -0.6276,
            0.6276,
            -0.6276,
            0.6276,
            -0.6276,
            -0.6276,
        ],
        "eri": True,
    },
    {
        "name": "NH3",
        "elements": [7, 1, 1, 1],
        "coords": [
            0.0,
            0.0,
            0.1173,
            0.0,
            0.9377,
            -0.2739,
            0.8121,
            -0.4689,
            -0.2739,
            -0.8121,
            -0.4689,
            -0.2739,
        ],
        "eri": True,
    },
    {
        "name": "HF",
        "elements": [9, 1],
        "coords": [0.0, 0.0, 0.0, 0.0, 0.0, 0.917],
        "eri": False,
    },
    {
        "name": "CO",
        "elements": [6, 8],
        "coords": [0.0, 0.0, 0.0, 0.0, 0.0, 1.128],
        "eri": False,
    },
    {
        "name": "C2H4",
        "elements": [6, 6, 1, 1, 1, 1],
        "coords": [
            0.0,
            0.0,
            0.6695,
            0.0,
            0.0,
            -0.6695,
            0.0,
            0.9289,
            1.2321,
            0.0,
            -0.9289,
            1.2321,
            0.0,
            0.9289,
            -1.2321,
            0.0,
            -0.9289,
            -1.2321,
        ],
        "eri": False,
    },
    {
        "name": "CH2O",
        "elements": [6, 8, 1, 1],
        "coords": [
            0.0,
            0.0,
            0.5349,
            0.0,
            0.0,
            -0.6701,
            0.0,
            0.9426,
            1.1165,
            0.0,
            -0.9426,
            1.1165,
        ],
        "eri": False,
    },
]


import numpy as np


def xh3_geometry(bond: float, angle_deg: float) -> tuple[list[float], list[float]]:
    """Central atom + three hydrogens (C3v pyramid). Returns (coords, elements)."""
    c = np.cos(np.radians(angle_deg))
    cb2 = (c + 0.5) / 1.5
    cb = np.sqrt(cb2)
    sb = np.sqrt(1.0 - cb2)
    coords = [0.0, 0.0, 0.0]
    for phi in (0.0, 120.0, 240.0):
        p = np.radians(phi)
        coords += [bond * sb * np.cos(p), bond * sb * np.sin(p), bond * cb]
    return coords, [1, 1, 1]


def xh2_geometry(bond: float, angle_deg: float) -> list[float]:
    c = np.cos(np.radians(angle_deg))
    cb2 = (c + 0.5) / 1.5
    cb = np.sqrt(cb2)
    sb = np.sqrt(1.0 - cb2)
    coords = [0.0, 0.0, 0.0]
    for phi in (0.0, 180.0):
        p = np.radians(phi)
        coords += [bond * sb * np.cos(p), bond * sb * np.sin(p), bond * cb]
    return coords


_ph3, _h3 = xh3_geometry(1.415, 93.5)
_h2s = xh2_geometry(1.34, 92.1)

MOLECULES += [
    {
        "name": "PH3",
        "elements": [15, *_h3],
        "coords": _ph3,
        "eri": False,
    },
    {
        "name": "H2S",
        "elements": [16, 1, 1],
        "coords": _h2s,
        "eri": False,
    },
    {
        "name": "HCl",
        "elements": [17, 1],
        "coords": [0.0, 0.0, 0.0, 0.0, 0.0, 1.275],
        "eri": False,
    },
]


def atom_spec(case: dict) -> list:
    atoms = []
    for z, (x, y, z_) in zip(
        case["elements"], np.array(case["coords"]).reshape(-1, 3)
    ):
        atoms.append((int(z), (float(x), float(y), float(z_))))
    return atoms


def run_scf(limit: int | None = None) -> Report:
    cases = []
    for case in MOLECULES[: limit or len(MOLECULES)]:
        mol = gto.M(
            atom=atom_spec(case),
            basis="sto-3g",
            unit="Angstrom",
            charge=0,
            spin=0,
            verbose=0,
        )
        s = mol.intor("int1e_ovlp")
        t = mol.intor("int1e_kin")
        v = mol.intor("int1e_nuc")
        entry = {
            "name": case["name"],
            "elements": case["elements"],
            "coords": case["coords"],
            "n_basis": int(mol.nao_nr()),
            "n_electrons": int(mol.nelectron),
            "overlap": s.flatten().tolist(),
            "kinetic": t.flatten().tolist(),
            "nuclear": v.flatten().tolist(),
        }
        if case.get("eri"):
            eri = mol.intor("int2e")
            entry["eri"] = eri.flatten().tolist()
        mf = pyscf_scf.RHF(mol)
        mf.conv_tol = 1e-12
        e_tot = mf.kernel()
        entry["e_tot"] = float(e_tot)
        entry["mo_energies"] = sorted(float(x) for x in mf.mo_energy)
        pop, charges = mf.mulliken_pop(verbose=0)
        entry["mulliken_charges"] = [float(q) for q in charges]
        # Exact RHF dipole (density integrals) and the Mulliken
        # point-charge dipole with the same convention as
        # sci-form's `dipole::compute_dipole` (F26).
        dip_exact = [float(v) for v in mf.dip_moment(unit="Debye")]
        entry["dipole_exact"] = dip_exact
        mu = (
            np.asarray(charges)[:, None]
            * np.array(case["coords"]).reshape(-1, 3)
        ).sum(axis=0) * DEBYE_EANG
        entry["dipole_mulliken"] = [float(v) for v in mu]
        cases.append(entry)

    payload = {
        "reference": REFERENCE,
        "basis": "sto-3g",
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F36",
        observable="scf_matrices_energy_mulliken",
        reference_tool=REFERENCE,
        units="hartree / a.u.",
        tolerance=1e-7,
        level="L2",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "S/T/V tolerance 1e-9 (Boys precision); ERI 1e-10 (F37); E_tot 1e-7 Ha; "
            "mo_energies 1e-6 Ha; Mulliken 0.005 e. Measured by "
            "tests/cross_validation.rs::{f36_scf, f37_eri}."
        ),
    )
    emit_fixture_report("pyscf_scf.json", payload, report)
    return report


# ---------------------------------------------------------------------------
# F38 — open-shell SCF (UHF/ROHF) reference energies and spin expectation
# ---------------------------------------------------------------------------

# (name, elements, coords, charge, multiplicity, method)
UHF_MOLECULES = [
    ("H2 singlet", [1, 1], [0.0, 0.0, 0.0, 0.0, 0.0, 0.74], 0, 1, "uhf"),
    ("H2 triplet", [1, 1], [0.0, 0.0, 0.0, 0.0, 0.0, 0.74], 0, 3, "uhf"),
    ("H2+ doublet", [1, 1], [0.0, 0.0, 0.0, 0.0, 0.0, 0.74], 1, 2, "rohf"),
    ("OH radical", [8, 1], [0.0, 0.0, 0.0, 0.0, 0.0, 0.969], 0, 2, "rohf"),
    ("CH3 radical", [6, 1, 1, 1], [0.0, 0.0, 0.0, 1.079, 0.0, 0.0, -0.5395, 0.9344, 0.0, -0.5395, -0.9344, 0.0], 0, 2, "rohf"),
    ("NH2 radical", [7, 1, 1], [0.0, 0.0, 0.144, 0.0, 0.791, -0.504, 0.0, -0.791, -0.504], 0, 2, "rohf"),
    ("O2 triplet", [8, 8], [0.0, 0.0, 0.0, 0.0, 0.0, 1.208], 0, 3, "uhf"),
    ("CH2 triplet", [6, 1, 1], [0.0, 0.0, 0.0, 0.0, 0.929, 0.588, 0.0, -0.929, 0.588], 0, 3, "uhf"),
]


def run_uhf() -> Report:
    cases = []
    for name, elements, coords, charge, multiplicity, method in UHF_MOLECULES:
        mol = gto.M(
            atom=atom_spec({"elements": elements, "coords": coords}),
            basis="sto-3g",
            unit="Angstrom",
            charge=charge,
            spin=multiplicity - 1,
            verbose=0,
        )

        def uhf(init_guess: str) -> dict:
            mf = pyscf_scf.UHF(mol)
            mf.conv_tol = 1e-12
            mf.init_guess = init_guess
            mf.max_cycle = 512
            e = mf.kernel()
            return {"e_tot": float(e), "s2": float(mf.spin_square()[0])}

        def rohf() -> dict:
            mf = pyscf_scf.ROHF(mol)
            mf.conv_tol = 1e-12
            mf.max_cycle = 512
            e = mf.kernel()
            return {"e_tot": float(e), "s2": float(mf.spin_square()[0])}

        entry = {
            "name": name,
            "elements": elements,
            "coords": coords,
            "charge": charge,
            "multiplicity": multiplicity,
            "method": method,
            # Default (atomic superposition) guess: the physical ground state.
            "uhf_default": uhf("atom"),
            # H_core guess: same starting point as sci-form's UHF solver.
            "uhf_hcore": uhf("hcore"),
            "rohf": rohf(),
            "n_electrons": int(mol.nelectron),
        }
        cases.append(entry)

    payload = {
        "reference": REFERENCE,
        "basis": "sto-3g",
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F38",
        observable="uhf_rohf_energy_s2",
        reference_tool=REFERENCE,
        units="hartree",
        tolerance=1e-6,
        level="L2",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "sci-form's UHF starts from the H_core guess, so the cross-check "
            "uses the matching PySCF UHF run (uhf_hcore). Tolerance 1e-6 Ha; "
            "<S^2> 1e-6. The default-guess (atom) references are stored too "
            "and document the initial-guess dependence for NH2/O2. "
            "Measured by tests/cross_validation.rs::f38_uhf_rohf."
        ),
    )
    emit_fixture_report("pyscf_uhf.json", payload, report)
    return report


# ---------------------------------------------------------------------------
# F39 — AO -> MO 4-index transform (hf::mo_transform)
# ---------------------------------------------------------------------------

AO2MO_MOLECULES = [m for m in MOLECULES if m["name"] in ("H2O", "CH4", "NH3")]


def run_ao2mo() -> Report:
    from pyscf import ao2mo as pyscf_ao2mo

    cases = []
    for case in AO2MO_MOLECULES:
        mol = gto.M(
            atom=atom_spec(case),
            basis="sto-3g",
            unit="Angstrom",
            charge=0,
            spin=0,
            verbose=0,
        )
        mf = pyscf_scf.RHF(mol)
        mf.conv_tol = 1e-12
        mf.kernel()
        c = np.asarray(mf.mo_coeff)
        ao_eri = mol.intor("int2e")
        mo_eri = pyscf_ao2mo.full(mol, c, compact=False)
        cases.append(
            {
                "name": case["name"],
                "elements": case["elements"],
                "coords": case["coords"],
                "n_ao": int(mol.nao_nr()),
                "n_mo": int(c.shape[1]),
                "mo_coeff": c.flatten().tolist(),
                "ao_eri": ao_eri.flatten().tolist(),
                "mo_eri": mo_eri.flatten().tolist(),
            }
        )

    payload = {
        "reference": REFERENCE,
        "basis": "sto-3g",
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F39",
        observable="ao_to_mo_transform_max_abs",
        reference_tool=f"{REFERENCE} ao2mo.full",
        units="hartree (chemist notation (pq|rs))",
        tolerance=1e-10,
        level="L2",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "Same C and same AO ERIs on both sides; tolerance 1e-10 abs per "
            "element. Measured by tests/cross_validation.rs::f39_ao2mo."
        ),
    )
    emit_fixture_report("pyscf_ao2mo.json", payload, report)
    return report


# ---------------------------------------------------------------------------
# F41 — CIS/TDA singlet excitation energies (hf::cis)
# ---------------------------------------------------------------------------

TDA_MOLECULES = [m for m in MOLECULES if m["name"] in ("H2", "H2O", "CH4", "NH3", "HF", "CO", "C2H4", "CH2O")]


def run_tda(nstates: int = 4) -> Report:
    from pyscf import tdscf

    cases = []
    for case in TDA_MOLECULES:
        mol = gto.M(
            atom=atom_spec(case),
            basis="sto-3g",
            unit="Angstrom",
            charge=0,
            spin=0,
            verbose=0,
        )
        mf = pyscf_scf.RHF(mol)
        mf.conv_tol = 1e-12
        e_hf = float(mf.kernel())
        td = tdscf.TDA(mf)
        td.nstates = nstates
        td.singlet = True
        td.kernel()
        energies_ev = [float(e) * HARTREE2EV for e in td.e]
        try:
            osc = [float(x) for x in td.oscillator_strength()]
        except Exception:
            osc = []
        cases.append(
            {
                "name": case["name"],
                "elements": case["elements"],
                "coords": case["coords"],
                "n_basis": int(mol.nao_nr()),
                "n_occupied": int(mol.nelectron) // 2,
                "e_tot": e_hf,
                "excitation_ev": energies_ev,
                "oscillator": osc,
            }
        )

    payload = {
        "reference": REFERENCE,
        "basis": "sto-3g",
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F41",
        observable="cis_singlet_excitations_ev",
        reference_tool=f"{REFERENCE} tdscf.TDA (singlets)",
        units="eV",
        tolerance=0.01,
        level="L2",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "CIS/TDA singlet excitation energies, tolerance 0.01 eV per the "
            "validation plan. Measured by tests/cross_validation.rs::f41_cis."
        ),
    )
    emit_fixture_report("pyscf_tda.json", payload, report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode",
        nargs="?",
        default="all",
        choices=["scf", "uhf", "ao2mo", "tda", "all"],
        help="which fixtures to generate (default: all)",
    )
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()
    if args.mode in ("scf", "all"):
        run_scf(args.limit)
    if args.mode in ("uhf", "all"):
        run_uhf()
    if args.mode in ("ao2mo", "all"):
        run_ao2mo()
    if args.mode in ("tda", "all"):
        run_tda()
    print("fixtures -> tests/fixtures/validation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
