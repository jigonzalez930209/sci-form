#!/usr/bin/env python3
"""PySCF oracle for sci-form's electrostatic-potential grid (ficha F78).

Generates `tests/fixtures/validation/pyscf_esp.json` with, for a neutral
(H2O) and a charged (OH-) test system:

  * the grid layout (origin in Ångström, spacing, dimensions) that
    `esp::compute_esp_grid` builds from the molecular bounding box, and
  * the **density-based** electrostatic potential on that grid, evaluated
    from PySCF's RHF/STO-3G density matrix with exact two-electron
    integrals (`df.incore.aux_e2` + `gto.fakemol_for_charges`, the same
    machinery as `pyscf.tools.cubegen.mep`).

sci-form's ESP grid uses the Mulliken point-charge model

    Phi(r) = sum_A q_A / max(|r - R_A|, d_min)

which is a *different* physical model from the density-based potential, so
the two can only be compared in the far field, where both are dominated by
the same multipole expansion.  The fixture therefore also stores the
minimum distance to a nucleus for every grid point, so that the Rust side
can select the far-field subset and gate on the agreement there.

Usage:
    python scripts/validation/oracle_pyscf_esp.py [all]
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
import math
import sys
import time

import numpy as np

from common import Report, emit_fixture_report

try:
    from pyscf import df, gto, scf as pyscf_scf
except ImportError as exc:  # pragma: no cover
    print(f"PySCF is required for this oracle: {exc}", file=sys.stderr)
    sys.exit(1)

PYSCF_VERSION = getattr(__import__("pyscf"), "__version__", "unknown")
REFERENCE = f"PySCF {PYSCF_VERSION}"

# sci-form's `EspGridConfig` used by the gate.
SPACING_ANG = 0.5
PADDING_ANG = 6.0
MIN_DISTANCE_ANG = 0.1

BOHR_TO_ANG = 0.529177210903

SYSTEMS = [
    {
        "name": "H2O",
        "elements": [8, 1, 1],
        "coords": [0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692],
        "charge": 0,
        "multiplicity": 1,
    },
    {
        "name": "OH-",
        "elements": [8, 1],
        "coords": [0.0, 0.0, 0.0, 0.0, 0.0, 0.97],
        "charge": -1,
        "multiplicity": 1,
    },
]


def grid_bounds(positions: np.ndarray, spacing: float, padding: float):
    """Mirror of `esp::esp::grid_bounds` (origin, dimensions)."""
    lo = positions.min(axis=0) - padding
    hi = positions.max(axis=0) + padding
    origin = [float(lo[k]) for k in range(3)]
    dims = [
        int(math.ceil((hi[k] - lo[k]) / spacing)) + 1 for k in range(3)
    ]
    return origin, dims


def density_esp(mol, dm: np.ndarray, points_bohr: np.ndarray) -> np.ndarray:
    """Electrostatic potential of the electron density + nuclei (a.u.).

    Uses the exact Coulomb integrals between the AO basis and the "fake"
    charge centres, exactly like `pyscf.tools.cubegen.mep`.
    """
    v_nuc = np.zeros(len(points_bohr))
    for i in range(mol.natm):
        r = mol.atom_coord(i)
        rp = r - points_bohr
        with np.errstate(divide="ignore"):
            v_nuc += mol.atom_charge(i) / np.einsum("xi,xi->x", rp, rp) ** 0.5
    v_ele = np.empty_like(v_nuc)
    for start in range(0, len(points_bohr), 400):
        chunk = points_bohr[start : start + 400]
        fakemol = gto.fakemol_for_charges(chunk)
        ints = df.incore.aux_e2(mol, fakemol)
        v_ele[start : start + 400] = np.einsum("ijp,ij->p", ints, dm)
    return v_nuc - v_ele


def point_charge_esp(
    charges: np.ndarray, positions_bohr: np.ndarray, points_bohr: np.ndarray
) -> np.ndarray:
    """Mulliken point-charge ESP with sci-form's distance clamp (a.u.)."""
    out = np.zeros(len(points_bohr))
    min_bohr = MIN_DISTANCE_ANG / BOHR_TO_ANG
    for q, r in zip(charges, positions_bohr):
        d = np.linalg.norm(points_bohr - r, axis=1)
        out += q / np.maximum(d, min_bohr)
    return out


def run() -> Report:
    cases = []
    for case in SYSTEMS:
        atoms = []
        for z, (x, y, z_) in zip(
            case["elements"], np.array(case["coords"]).reshape(-1, 3)
        ):
            atoms.append((int(z), (float(x), float(y), float(z_))))
        mol = gto.M(
            atom=atoms,
            basis="sto-3g",
            unit="Angstrom",
            charge=case["charge"],
            spin=case["multiplicity"] - 1,
            verbose=0,
        )
        mf = pyscf_scf.RHF(mol)
        mf.conv_tol = 1e-12
        mf.max_cycle = 512
        e_tot = float(mf.kernel())
        dm = mf.make_rdm1()
        pop, charges = mf.mulliken_pop(verbose=0)
        charges = np.asarray(charges, dtype=float)

        positions_ang = np.array(case["coords"]).reshape(-1, 3)
        origin, dims = grid_bounds(positions_ang, SPACING_ANG, PADDING_ANG)

        # grid points in the same order as `esp::compute_esp_grid`
        ix = np.arange(dims[0]) * SPACING_ANG + origin[0]
        iy = np.arange(dims[1]) * SPACING_ANG + origin[1]
        iz = np.arange(dims[2]) * SPACING_ANG + origin[2]
        points = np.stack(
            np.meshgrid(ix, iy, iz, indexing="ij"), axis=-1
        ).reshape(-1, 3)
        points_bohr = points / BOHR_TO_ANG

        v_dens = density_esp(mol, dm, points_bohr)

        cases.append(
            {
                "name": case["name"],
                "elements": case["elements"],
                "coords": case["coords"],
                "charge": case["charge"],
                "multiplicity": case["multiplicity"],
                "e_tot": e_tot,
                "mulliken_charges": [float(q) for q in charges],
                "total_charge": float(sum(charges)),
                "origin": origin,
                "spacing": SPACING_ANG,
                "dims": dims,
                "min_distance_ang": MIN_DISTANCE_ANG,
                # `null` marks the grid points that sit exactly on a nucleus,
                # where the density-based potential has the nuclear
                # 1/r singularity (sci-form clamps the distance instead).
                "esp_density": [
                    None if not np.isfinite(v) else float(f"{v:.10e}") for v in v_dens
                ],
            }
        )

    payload = {
        "reference": REFERENCE,
        "basis": "sto-3g",
        "generated": time.strftime("%Y-%m-%d"),
        "config": {
            "spacing_ang": SPACING_ANG,
            "padding_ang": PADDING_ANG,
            "min_distance_ang": MIN_DISTANCE_ANG,
        },
        "cases": cases,
    }
    report = Report(
        ficha="F78",
        observable="esp_far_field_max_abs_dev",
        reference_tool=f"{REFERENCE} RHF/STO-3G density ESP (cubegen scheme)",
        units="hartree/e",
        tolerance=1e-2,
        level="L2",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "sci-form's ESP grid uses Mulliken point charges, PySCF's "
            "reference uses the RHF density: the two models differ by "
            "construction in the near field and are compared in the far "
            "field, where both are described by the same multipole "
            "expansion. Measured by tests/pyscf_parity.rs::f78_esp_grid."
        ),
    )
    emit_fixture_report("pyscf_esp.json", payload, report)
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
    print("fixtures -> tests/fixtures/validation/pyscf_esp.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
