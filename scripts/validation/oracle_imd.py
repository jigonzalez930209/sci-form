#!/usr/bin/env python3
"""Oracle F85: IMD harmonic steering (E = k/2|r-r0|^2, F = -k(r-r0)).

Independent numpy reference for SteeringForces::accumulate, including
multi-force additivity. Units: kcal/mol, Angstrom.

Usage: oracle_imd.py  (writes tests/fixtures/validation/imd_steering.json)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "tests" / "fixtures" / "validation" / "imd_steering.json"


def main() -> int:
    rng = np.random.default_rng(42)
    cases = []
    for n_forces in (1, 2, 4):
        for _ in range(4):
            n_atoms = 6
            pos = rng.standard_normal(n_atoms * 3).tolist()
            forces = []
            seen = {}
            for _ in range(n_forces):
                idx = int(rng.integers(0, n_atoms))
                # apply() replaces any existing force on the atom: last wins.
                seen[idx] = (rng.standard_normal(3), float(rng.uniform(1.0, 200.0)))
            e_ref, f_ref = 0.0, np.zeros(n_atoms * 3)
            for idx, (tgt, k) in seen.items():
                d = np.array(pos[3 * idx:3 * idx + 3]) - tgt
                e_ref += 0.5 * k * float(d @ d)
                f_ref[3 * idx:3 * idx + 3] -= k * d
                forces.append({"atom": idx, "target": tgt.tolist(), "k": k})
            cases.append({"n_atoms": n_atoms, "positions": pos, "forces": forces,
                          "energy_ref": float(e_ref), "forces_ref": f_ref.tolist()})
    payload = {
        "ficha": "F85",
        "oracle": "oracle_imd.py (harmonic closed form, seed 42)",
        "tolerance": 1e-12,
        "cases": cases,
    }
    OUT.write_text(json.dumps(payload) + "\n")
    print(f"wrote {OUT} ({len(cases)} cases)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
