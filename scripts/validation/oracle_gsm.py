#!/usr/bin/env python3
"""Oracle F82: Muller-Brown MEP/TS via ASE CI-NEB (external path optimizer).

- Muller-Brown surface + analytic gradient (textbook constants).
- Minima from scipy minimize (high precision, 3 basins).
- ASE NEB (7 images, FIRE, climb=True) with an analytic MB calculator:
  MEP profile + TS energy/coords + barrier. The Rust gate runs OUR
  growing-string implementation on the same surface and must reproduce
  the barrier/TS.

Usage: oracle_gsm.py  (writes tests/fixtures/validation/gsm_mullerbrown.json)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.mep import NEB
from ase.optimize import FIRE
from scipy.optimize import minimize

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "tests" / "fixtures" / "validation" / "gsm_mullerbrown.json"

A = (-200.0, -100.0, -170.0, 15.0)
AA = (-1.0, -1.0, -6.5, 0.7)
BB = (0.0, 0.0, 11.0, 0.6)
X0 = (1.0, 0.0, -0.5, -1.0)
Y0 = (0.0, 0.5, 1.5, 1.0)
CC = (-10.0, -10.0, -6.5, 0.7)


def mb(x: float, y: float) -> float:
    return sum(
        a * np.exp(aa * (x - x0) ** 2 + bb * (x - x0) * (y - y0) + cc * (y - y0) ** 2)
        for a, aa, bb, x0, y0, cc in zip(A, AA, BB, X0, Y0, CC)
    )


def mb_grad(x: float, y: float) -> np.ndarray:
    gx = gy = 0.0
    for a, aa, bb, x0, y0, cc in zip(A, AA, BB, X0, Y0, CC):
        e = np.exp(aa * (x - x0) ** 2 + bb * (x - x0) * (y - y0) + cc * (y - y0) ** 2)
        gx += a * e * (2 * aa * (x - x0) + bb * (y - y0))
        gy += a * e * (bb * (x - x0) + 2 * cc * (y - y0))
    return np.array([gx, gy])


class MB(Calculator):
    implemented_properties = ["energy", "forces"]

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        x, y = atoms.positions[0, :2]
        self.results["energy"] = float(mb(x, y))
        self.results["forces"] = np.array([[-mb_grad(x, y)[0], -mb_grad(x, y)[1], 0.0]])


def main() -> int:
    mins = []
    for guess in [(-0.5, 1.4), (0.6, 0.0), (0.0, 0.5)]:
        # Two-stage: Nelder-Mead (robust on the narrow valley) then BFGS polish.
        r0 = minimize(lambda p: mb(*p), guess, method="Nelder-Mead",
                      options={"xatol": 1e-6, "fatol": 1e-6, "maxiter": 20000})
        r = minimize(lambda p: mb(*p), r0.x, jac=lambda p: mb_grad(*p),
                     method="BFGS", options={"gtol": 1e-7, "maxiter": 10000})
        gnorm = float(np.linalg.norm(mb_grad(*r.x)))
        assert gnorm < 1e-5, f"{guess}: |grad|={gnorm:.2e}"
        mins.append({"coords": r.x.tolist(), "energy": float(r.fun)})
    mins.sort(key=lambda m: m["energy"])
    m1, m2 = mins[0], mins[1]
    print("minima:", [(round(m["energy"], 4), [round(v, 4) for v in m["coords"]]) for m in mins])
    images = [Atoms("H", positions=[[m1["coords"][0], m1["coords"][1], 0.0]])]
    for i in range(1, 6):
        t = i / 6
        p = (1 - t) * np.array(m1["coords"]) + t * np.array(m2["coords"])
        images.append(Atoms("H", positions=[[p[0], p[1], 0.0]]))
    images.append(Atoms("H", positions=[[m2["coords"][0], m2["coords"][1], 0.0]]))
    for im in images:
        im.calc = MB()
    neb = NEB(images, climb=True, method="improvedtangent")
    FIRE(neb, logfile=None).run(fmax=1e-4, steps=2000)
    mep = [{"coords": [float(im.positions[0, 0]), float(im.positions[0, 1])],
            "energy": float(im.get_potential_energy())} for im in images]
    ts = max(mep, key=lambda p: p["energy"])
    barrier = ts["energy"] - m1["energy"]
    print(f"TS: E={ts['energy']:.6f} at {[round(v, 4) for v in ts['coords']]} barrier={barrier:.6f}")
    payload = {
        "ficha": "F82",
        "oracle": "oracle_gsm.py (ASE CI-NEB/FIRE on Muller-Brown)",
        "minima": mins,
        "ts_energy": ts["energy"],
        "ts_coords": ts["coords"],
        "barrier": barrier,
        "mep": mep,
        "barrier_tolerance_rel": 0.02,
        "ts_location_tolerance": 0.1,
    }
    OUT.write_text(json.dumps(payload, indent=1) + "\n")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
