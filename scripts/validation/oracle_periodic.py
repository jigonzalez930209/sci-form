#!/usr/bin/env python3
"""Oracle F88: MP meshes (ASE), Bloch phases + 1D chain bands (closed form).

- Monkhorst-Pack grids vs ase.dft.kpoints.monkhorst_pack (independent).
- Gamma-centered grids vs textbook i/N folded to [-0.5, 0.5).
- Bloch phases e^{i2pik.R} via numpy.
- 1D tight-binding chain E(k) = e + 2t cos(2 pi k): analytic bands on
  the Gamma->X path for 1- and 2-orbital cells.
- NOTE: build_bloch_hamiltonian keeps the real (cos) part only; the
  oracle references below are real by R/-R symmetry. Complex H(k) for
  non-centrosymmetric cells is unimplemented (documented gap).

Usage: oracle_periodic.py  (writes tests/fixtures/validation/periodic_ref.json)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from ase.dft.kpoints import monkhorst_pack

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "tests" / "fixtures" / "validation" / "periodic_ref.json"


def mp_closed(n: int) -> list[float]:
    return sorted((2 * i + 1 - n) / (2 * n) for i in range(n))


def gamma_closed(n: int) -> list[float]:
    out = []
    for i in range(n):
        k = i / n
        out.append(k - 1.0 if k >= 0.5 else k)
    return sorted(out)


def main() -> int:
    rng = np.random.default_rng(42)
    meshes = []
    for grid in ([1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 2], [2, 3, 4]):
        ase_pts = monkhorst_pack(grid)
        meshes.append({
            "grid": grid, "centering": "MonkhorstPack",
            "ase_ref": sorted([sorted(p.tolist()) for p in ase_pts]),
            "closed_ref_1d": {str(n): mp_closed(n) for n in sorted(set(grid))},
        })
        meshes.append({
            "grid": grid, "centering": "GammaCentered",
            "closed_ref_1d": {str(n): gamma_closed(n) for n in sorted(set(grid))},
        })
    phases = []
    for _ in range(12):
        k = rng.uniform(-0.5, 0.5, 3).tolist()
        R = [int(v) for v in rng.integers(-3, 4, 3)]
        th = 2 * np.pi * float(np.dot(k, R))
        phases.append({"k": k, "R": R, "cos_ref": float(np.cos(th)), "sin_ref": float(np.sin(th))})
    chains = []
    for eps, t in ((0.0, -1.0), (-2.5, -0.8), (1.0, 0.5)):
        ks = np.linspace(0, 0.5, 25).tolist()
        bands = [eps + 2 * t * np.cos(2 * np.pi * k) for k in ks]
        chains.append({"eps": eps, "t": t, "k_path": ks, "bands_ref": bands})
    chains2 = []
    for p in ({"es": -1.0, "ts": -0.5, "ep": 2.0, "tp": -0.2},):
        ks = np.linspace(0, 0.5, 25).tolist()
        b1 = [p["es"] + 2 * p["ts"] * np.cos(2 * np.pi * k) for k in ks]
        b2 = [p["ep"] + 2 * p["tp"] * np.cos(2 * np.pi * k) for k in ks]
        chains2.append({**p, "k_path": ks, "band1_ref": b1, "band2_ref": b2})
    payload = {
        "ficha": "F88",
        "oracle": "oracle_periodic.py (ASE monkhorst_pack + closed forms, seed 42)",
        "note": ("build_bloch_hamiltonian keeps cos only; non-centrosymmetric "
                 "complex H(k) unimplemented (documented gap)."),
        "mesh_tolerance": 1e-12,
        "phase_tolerance": 1e-12,
        "band_tolerance": 1e-9,
        "meshes": meshes,
        "phases": phases,
        "chains": chains,
        "chains2": chains2,
    }
    OUT.write_text(json.dumps(payload) + "\n")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
