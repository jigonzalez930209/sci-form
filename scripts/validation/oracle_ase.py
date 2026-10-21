#!/usr/bin/env python3
"""ASE oracle for sci-form variable-cell validation (ficha F85).

Builds a small periodic system, computes the stress tensor and the relaxed
cell with ASE (Lennard-Jones calculator + cell filter) and writes the
reference values to `tests/fixtures/validation/ase_stress.json`.

The sci-form side (`src/materials/variable_cell.rs`) is recomputed by
`tests/rdkit_chem_parity.rs` from the committed fixture, so CI does not need
ASE.

Conventions
-----------
* ASE returns the stress in eV/Angstrom^3 with the sign convention
  `sigma_ij = (1/V) dE/deps_ij` (the same convention as sci-form's
  `numerical_stress`); no sign flip is needed.
* The Lennard-Jones calculator of ASE shifts the pair energy by `u(rc)` so
  that it is continuous at the cutoff, and it counts the interaction of an
  atom with its own periodic images with half weight (the pair appears only
  once in the neighbour list while the energy carries a factor 1/2).  The
  Rust model in the test reproduces both conventions.

Usage:
    python scripts/validation/oracle_ase.py
"""

from __future__ import annotations

import sys
import time

import numpy as np

from common import FIXTURE_DIR, REPO, Report, emit_fixture_report

try:
    from ase import Atoms
    from ase.calculators.lj import LennardJones
    from ase.filters import UnitCellFilter
    from ase.optimize import BFGS
except ImportError as exc:  # pragma: no cover
    print(f"ASE is required for this oracle: {exc}", file=sys.stderr)
    sys.exit(1)

import ase

REFERENCE = f"ASE {ase.__version__}"

# Lennard-Jones parameters (eV, Angstrom).
#
# The cutoff is chosen so that exactly the eight body-diagonal neighbours
# of the two-atom cell interact: 0.866*a < rc < a for the whole range of
# cell lengths that the relaxation visits (the equilibrium of the eight
# neighbour interaction is a = 2^(1/6)*sigma/0.866 = 1.2960 A).  No pair
# crosses the cutoff during a strain perturbation, so the numerical and the
# analytic stress are directly comparable.
SIGMA = 1.0
EPSILON = 1.0
RC = 1.21

# bcc-like two-atom cell: one atom at the origin, one at the body centre.
SYSTEMS = [
    ("bcc-a1.25", 1.25),
    ("bcc-a1.30", 1.30),
    ("bcc-a1.35", 1.35),
]


def build(a: float) -> Atoms:
    cell = np.eye(3) * a
    positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]) * a
    atoms = Atoms(symbols="Ar2", positions=positions, cell=cell, pbc=True)
    atoms.calc = LennardJones(sigma=SIGMA, epsilon=EPSILON, rc=RC)
    return atoms


def relaxed_lattice(a: float) -> tuple[list[float], float]:
    """Relax the cell with ASE (BFGS on a unit-cell filter)."""
    atoms = build(a)
    atoms.set_constraint()
    filter_ = UnitCellFilter(atoms, mask=None)
    opt = BFGS(filter_, logfile=None)
    opt.run(fmax=1e-7, steps=500)
    return [float(x) for x in atoms.get_cell().array.ravel()], float(
        atoms.get_potential_energy()
    )


def run() -> Report:
    cases = []
    for name, a in SYSTEMS:
        atoms = build(a)
        energy = float(atoms.get_potential_energy())
        stress = np.array(atoms.get_stress(voigt=False))
        relaxed_cell, relaxed_energy = relaxed_lattice(a)
        cases.append(
            {
                "name": name,
                "lattice": [[float(x) for x in row] for row in atoms.get_cell().array],
                "positions": [[float(x) for x in row] for row in atoms.get_positions()],
                "symbols": list(atoms.get_chemical_symbols()),
                "sigma": SIGMA,
                "epsilon": EPSILON,
                "rc": RC,
                "energy_ev": energy,
                "stress_ev_per_a3": [[float(x) for x in row] for row in stress],
                "relaxed_lattice": relaxed_cell,
                "relaxed_energy_ev": relaxed_energy,
            }
        )

    payload = {
        "reference": REFERENCE,
        "calculator": "ase.calculators.lj.LennardJones",
        "energy_shifted_at_cutoff": True,
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F85",
        observable="variable_cell_stress_vs_ase",
        reference_tool=f"{REFERENCE} LennardJones stress + UnitCellFilter",
        units="eV/Angstrom^3",
        tolerance=0.0,
        level="L2",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "stress convention: sigma_ij = (1/V) dE/deps_ij (identical to "
            "sci-form's numerical_stress); the relaxed cell comes from an "
            "ASE BFGS relaxation with a UnitCellFilter."
        ),
    )
    emit_fixture_report("ase_stress.json", payload, report)
    return report


def main() -> int:
    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    run()
    print(f"fixtures -> {FIXTURE_DIR}")
    print(f"summary  -> {REPO / 'reports' / 'validation' / 'summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
