#!/usr/bin/env python3
"""PySCF oracle for sci-form's IR intensities (ficha F79b).

Generates `tests/fixtures/validation/pyscf_ir.json` with, for a small
molecule (H2O, CH4) in a closed-shell RHF/STO-3G reference:

  * `frequencies_cm1` — the vibrational frequencies from PySCF's analytic
    Hessian (`pyscf.hessian`), converted with the same convention as
    sci-form's `spectroscopy::hessian` (Hartree/(amu·Bohr²) mass-weighted
    eigenvalues, 219474.63 cm⁻¹ per Hartree).
  * `intensity_exact_km_mol` — the IR intensity computed from the
    **density** dipole derivative along each normal mode.
  * `intensity_mulliken_km_mol` — the same, but with the Mulliken dipole
    (the model sci-form uses), so that the two error sources (dipole model
    and numerical Hessian) can be separated.

IR intensities use the standard convention

    I[km/mol] = 42.2561 * (dmu/dQ)^2      with dmu/dQ in Debye/(Angstrom*sqrt(amu))

which is the unit conversion published by Gaussian
("1 Debye^2-angstrom^-2-amu^-1 = 42.2561 km-mol^-1").

PySCF 2.14 has no IR-intensity API (`pyscf.hessian.thermo` carries a
"TODO: IR intensity"), so the intensities are computed here from PySCF's
analytic Hessian and its SCF dipoles.

Usage:
    python scripts/validation/oracle_pyscf_ir.py [all]
"""

from __future__ import annotations

import os

# Pin the BLAS/OpenMP thread count *before* numpy is imported: PySCF's
# diagonalisation otherwise sums in a thread-dependent order, which changes
# the last digits of the reference values and makes the oracle
# non-reproducible (the fixtures must be byte-identical on every run).
for _thread_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_thread_var] = "1"

import argparse
import sys
import time

import numpy as np

from common import Report, emit_fixture_report

try:
    from pyscf import gto, scf as pyscf_scf
    from pyscf.hessian import rhf as hess_rhf, thermo
except ImportError as exc:  # pragma: no cover
    print(f"PySCF is required for this oracle: {exc}", file=sys.stderr)
    sys.exit(1)

PYSCF_VERSION = getattr(__import__("pyscf"), "__version__", "unknown")
REFERENCE = f"PySCF {PYSCF_VERSION}"

# Atomic masses (amu) as used by sci-form's `spectroscopy::hessian::atomic_mass`.
MASS = {1: 1.00794, 6: 12.011, 8: 15.999}

# Displacement step (Bohr) used for the numerical dipole derivative.
DELTA = 0.005

# Gaussian's IR intensity unit conversion.
IR_UNIT = 42.2561
DEBYE_PER_EA0 = 2.541746
ANG_PER_BOHR = 0.529177

SYSTEMS = [
    {
        "name": "H2O",
        "elements": [8, 1, 1],
        "coords": [0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692],
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
    },
]


def dipole_exact(mf, mol) -> np.ndarray:
    """Density dipole moment (a.u.)."""
    dm = mf.make_rdm1()
    with mol.with_common_orig(np.zeros(3)):
        ao = mol.intor_symmetric("int1e_r", comp=3)
    el = np.einsum("xij,ji->x", ao, dm).real
    nucl = np.einsum("i,ix->x", mol.atom_charges(), mol.atom_coords())
    return nucl - el


def dipole_mulliken(mf, mol) -> np.ndarray:
    """Mulliken point-charge dipole (a.u.)."""
    _, charges = mf.mulliken_pop(verbose=0)
    return (np.asarray(charges)[:, None] * mol.atom_coords()).sum(axis=0)


def intensity_km_mol(dmu_ddelta: np.ndarray) -> float:
    """IR intensity (km/mol) from d(mu)/d(delta) in e*a0 per Bohr."""
    # d(mu)/d(Q) in Debye/(Angstrom*sqrt(amu)):
    #   dmu/ddelta [e*a0/bohr] * (2.541746 D / e*a0) / (0.529177 A / bohr)
    x = np.linalg.norm(dmu_ddelta) * DEBYE_PER_EA0 / ANG_PER_BOHR
    return IR_UNIT * x * x


def run() -> Report:
    cases = []
    for case in SYSTEMS:
        atoms = []
        for z, (x, y, z_) in zip(
            case["elements"], np.array(case["coords"]).reshape(-1, 3)
        ):
            atoms.append((int(z), (float(x), float(y), float(z_))))
        mol = gto.M(
            atom=atoms, basis="sto-3g", unit="Angstrom", charge=0, spin=0, verbose=0
        )
        mf = pyscf_scf.RHF(mol)
        mf.conv_tol = 1e-12
        mf.max_cycle = 512
        e_tot = float(mf.kernel())

        hess = hess_rhf.Hessian(mf).kernel()
        analysis = thermo.harmonic_analysis(mol, hess)
        freqs = np.asarray(analysis["freq_wavenumber"], dtype=float)
        # Cartesian displacement of each normal mode (per unit of the
        # mass-weighted coordinate, in Bohr, with the mass in amu).
        norm_mode = np.asarray(analysis["norm_mode"], dtype=float)

        # `harmonic_analysis` projects out the translations/rotations, so it
        # only reports the vibrational modes.
        n_modes = len(freqs)
        intensities_exact = []
        intensities_mulliken = []
        dmu_exact = []
        dmu_mulliken = []
        for k in range(n_modes):
            freq = float(freqs[k])
            if freq < 1.0:
                continue  # translation / rotation
            # `norm_mode` has the shape (n_vib, natm, 3)
            disp = norm_mode[k]  # (atom, xyz)
            res_exact = []
            res_mull = []
            for sign in (1.0, -1.0):
                new_coords = mol.atom_coords() + sign * DELTA * disp
                mol2 = gto.M(
                    atom=[
                        (a[0], tuple(c))
                        for a, c in zip(
                            [
                                (z, (0.0, 0.0, 0.0))
                                for z in case["elements"]
                            ],
                            new_coords,
                        )
                    ],
                    basis="sto-3g",
                    unit="Bohr",
                    verbose=0,
                )
                mf2 = pyscf_scf.RHF(mol2)
                mf2.conv_tol = 1e-12
                mf2.max_cycle = 512
                mf2.kernel()
                res_exact.append(dipole_exact(mf2, mol2))
                res_mull.append(dipole_mulliken(mf2, mol2))
            dd_exact = (res_exact[0] - res_exact[1]) / (2.0 * DELTA)
            dd_mull = (res_mull[0] - res_mull[1]) / (2.0 * DELTA)
            dmu_exact.append([float(x) for x in dd_exact])
            dmu_mulliken.append([float(x) for x in dd_mull])
            intensities_exact.append(intensity_km_mol(dd_exact))
            intensities_mulliken.append(intensity_km_mol(dd_mull))

        cases.append(
            {
                "name": case["name"],
                "elements": case["elements"],
                "coords": case["coords"],
                "e_tot": e_tot,
                "frequencies_cm1": [float(f) for f in freqs],
                "intensity_exact_km_mol": intensities_exact,
                "intensity_mulliken_km_mol": intensities_mulliken,
                "dmu_exact": dmu_exact,
                "dmu_mulliken": dmu_mulliken,
                "delta_bohr": DELTA,
            }
        )

    payload = {
        "reference": REFERENCE,
        "basis": "sto-3g",
        "generated": time.strftime("%Y-%m-%d"),
        "config": {
            "delta_bohr": DELTA,
            "ir_unit": IR_UNIT,
            "masses": {str(k): v for k, v in MASS.items()},
        },
        "cases": cases,
    }
    report = Report(
        ficha="F79b",
        observable="ir_intensity_vs_mulliken_reference_max_abs",
        reference_tool=f"{REFERENCE} analytic Hessian + SCF dipole derivatives",
        units="km/mol",
        tolerance=0.5,
        level="L2",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "IR intensities (km/mol) from PySCF's analytic Hessian and its "
            "SCF dipoles; PySCF 2.14 has no IR-intensity API, so the "
            "standard I = 42.2561 (dmu/dQ)^2 conversion is applied here. "
            "Measured by tests/pyscf_parity.rs::f79b_ir_intensities."
        ),
    )
    emit_fixture_report("pyscf_ir.json", payload, report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode",
        nargs="?",
        default="all",
        choices=["all"],
        help="which fixtures to generate (default: all)",
    )
    parser.parse_args()
    run()
    print("fixtures -> tests/fixtures/validation/pyscf_ir.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
