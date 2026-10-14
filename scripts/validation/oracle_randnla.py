#!/usr/bin/env python3
"""Oracle F93: randomized Nyström vs deterministic references.

Regimes (all seeded):
- characterization: overlap-like S = I + 0.3 Q diag(exp(-i/5)) Q^T
  (realistic overlap deviations), rank k = N/4. Single-shot sketches
  give ~1e-1 here, NOT 1e-3: resolving clustered eigenvalues needs
  power iterations/oversampling the implementation lacks. The gate
  locks the measured accuracy and documents the gap (F95-style).
- valid regime: near-identity S (amp 0.05, the documented complement~=I
  regime) where the randomized path (fallback OFF) must meet 5e-3.
- production: default config (fallback ON) must match scipy exactly.

References: numpy.linalg.eigh, scipy.linalg.sqrtm, scipy.linalg.eigh(H, S).

Usage: oracle_randnla.py  (writes tests/fixtures/validation/randnla_numpy.json)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh as seigh
from scipy.linalg import sqrtm

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "tests" / "fixtures" / "validation" / "randnla_numpy.json"


def overlap_like(rng: np.random.Generator, n: int, amp: float) -> np.ndarray:
    lam = np.exp(-np.arange(n) / 5.0)
    q, _ = np.linalg.qr(rng.standard_normal((n, n)))
    s = np.eye(n) + amp * (q * lam) @ q.T
    return (s + s.T) / 2


def main() -> int:
    rng = np.random.default_rng(42)
    character, valid, production = [], [], []
    for n in (40, 64):
        s = overlap_like(rng, n, 0.3)
        k = n // 4
        w = np.linalg.eigh(s)[0][::-1][:k]
        character.append({"n": n, "k": k, "matrix": s.tolist(), "top_eigenvalues": w.tolist()})
    for n in (40, 64, 96):
        s = overlap_like(rng, n, 0.05)
        k = n // 4
        w = np.linalg.eigh(s)[0][::-1][:k]
        sq = sqrtm(s)
        a = rng.standard_normal((n, n))
        h = (a + a.T) / 2
        g = sorted(float(v) for v in seigh(h, s, eigvals_only=True))
        valid.append({
            "n": n, "k": k, "matrix": s.tolist(), "top_eigenvalues": w.tolist(),
            "inv_sqrt": np.linalg.inv(sq).tolist(),
            "h": h.tolist(), "energies": g,
        })
    for n in (24, 40, 96):
        s = overlap_like(rng, n, 0.3)
        a = rng.standard_normal((n, n))
        h = (a + a.T) / 2
        g = sorted(float(v) for v in seigh(h, s, eigvals_only=True))
        production.append({"n": n, "h": h.tolist(), "s": s.tolist(), "energies": g})
    payload = {
        "ficha": "F93",
        "oracle": "oracle_randnla.py (numpy eigh + scipy sqrtm/eigh, seed 42)",
        "note": ("single-shot rank-N/4 sketches resolve clustered spectra to ~1e-1, "
                 "not 1e-3; 1e-3 needs power iterations/oversampling (promotion criterion, "
                 "not met). Randomized path gated at 5e-3 in the near-identity regime; "
                 "default config (fallback) gated exact."),
        "characterization_tolerance_rel": 0.5,
        "valid_regime_tolerance": 5e-3,
        "inv_sqrt_tolerance_rel": 1e-2,
        "production_tolerance": 1e-6,
        "characterization": character,
        "valid_regime": valid,
        "production": production,
    }
    OUT.write_text(json.dumps(payload) + "\n")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
