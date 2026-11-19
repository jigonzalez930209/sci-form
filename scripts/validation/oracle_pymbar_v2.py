#!/usr/bin/env python3
"""pymbar oracle for the alchemical free-energy estimators (ficha F109, V2-09).

Draws seeded samples from three 1D harmonic oscillators with reduced
potentials u_k(x) = K_k x^2 / 2 (beta = 1), so the exact answer is
Delta f_0k = ln(K_k / K_0) / 2. The same samples are solved with pymbar
(MBAR, BAR, EXP) and written to ``tests/fixtures/validation/v2_pymbar.json``;
``tests/v2_validation.rs`` runs sci-form on the identical samples.

Usage::

    .venv-validation/bin/python scripts/validation/oracle_pymbar_v2.py
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
import pymbar

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "tests" / "fixtures" / "validation" / "v2_pymbar.json"

SPRINGS = [1.0, 2.0, 4.0]
N_PER_STATE = 500
SEED = 20261119


def main() -> int:
    rng = np.random.default_rng(SEED)
    k = np.array(SPRINGS)
    samples = []
    for state, spring in enumerate(k):
        x = rng.normal(0.0, 1.0 / math.sqrt(spring), N_PER_STATE)
        for value in x:
            samples.append({"origin_state": state, "reduced_potentials": (0.5 * k * value * value).tolist()})
    u_kn = np.array([s["reduced_potentials"] for s in samples]).T
    n_k = np.array([N_PER_STATE] * len(k))

    mbar = pymbar.MBAR(u_kn, n_k, solver_protocol="robust")
    f_mbar = mbar.compute_free_energy_differences()["Delta_f"][0].tolist()

    # Two-state BAR between states 0 and 1 (forward/backward work values).
    w_f = u_kn[1, :N_PER_STATE] - u_kn[0, :N_PER_STATE]
    w_r = u_kn[0, N_PER_STATE : 2 * N_PER_STATE] - u_kn[1, N_PER_STATE : 2 * N_PER_STATE]
    bar = pymbar.other_estimators.bar(w_f, w_r)["Delta_f"]
    # One-directional exponential averaging from state 0 to every state.
    exp_0k = [0.0] + [
        pymbar.other_estimators.exp(u_kn[j, :N_PER_STATE] - u_kn[0, :N_PER_STATE])["Delta_f"]
        for j in range(1, len(k))
    ]

    out = {
        "provenance": {
            "tool": f"pymbar {pymbar.__version__}",
            "script": "scripts/validation/oracle_pymbar_v2.py",
            "seed": SEED,
            "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
        },
        "springs": SPRINGS,
        "exact_delta_f": [0.5 * math.log(s / SPRINGS[0]) for s in SPRINGS],
        "samples": samples,
        "mbar_delta_f": f_mbar,
        "bar_delta_f_01": float(bar),
        "exp_delta_f": [float(v) for v in exp_0k],
    }
    OUT.write_text(json.dumps(out) + "\n")
    print(f"wrote {OUT.relative_to(REPO)}: MBAR {f_mbar}, exact {out['exact_delta_f']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
