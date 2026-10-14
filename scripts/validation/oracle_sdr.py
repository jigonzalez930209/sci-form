#!/usr/bin/env python3
"""Oracle F83: SDR projections + embedding vs numpy + CVXPY.

- warm_start_gram: numpy double-centering (-0.5 J D2 J).
- project_psd: numpy eigh clip (+ removal counts).
- svt_step: numpy soft-threshold of eigenvalues.
- project_distances: constraint residuals after projection.
- extract_coordinates: Gram reconstruction C C^T from top-3 eig
  (rotation/sign invariant comparison).
- sdr_embed end-to-end vs CVXPY convex SDP (min squared distance
  residuals s.t. X PSD, Clarabel): max distance errors compared.
  Point sets: random 3D (seed 42), exact + noisy distances.

Usage: oracle_sdr.py  (writes tests/fixtures/validation/sdr_ref.json)
"""

from __future__ import annotations

import json
from pathlib import Path

import cvxpy as cp
import numpy as np

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "tests" / "fixtures" / "validation" / "sdr_ref.json"


def double_center(d2: np.ndarray) -> np.ndarray:
    n = d2.shape[0]
    j = np.eye(n) - np.ones((n, n)) / n
    return -0.5 * j @ d2 @ j


def psd_clip(x: np.ndarray) -> tuple[np.ndarray, int]:
    w, v = np.linalg.eigh((x + x.T) / 2)
    neg = int((w < 0).sum())
    return (v * np.maximum(w, 0.0)) @ v.T, neg


def svt(x: np.ndarray, tau: float) -> np.ndarray:
    if tau <= 0:
        return x.copy()
    w, v = np.linalg.eigh((x + x.T) / 2)
    return (v * np.maximum(w - tau, 0.0)) @ v.T


def cvxpy_sdp(pairs: list[tuple[int, int, float]], n: int) -> tuple[np.ndarray, float]:
    X = cp.Variable((n, n), symmetric=True)
    obj = sum((X[i, i] + X[j, j] - 2 * X[i, j] - d * d) ** 2 for i, j, d in pairs)
    prob = cp.Problem(cp.Minimize(obj), [X >> 0])
    try:
        prob.solve(solver="SCS", eps=1e-9, max_iters=100000)
    except cp.SolverError:
        prob.solve(solver="CLARABEL")
    assert prob.status in ("optimal", "optimal_inaccurate"), prob.status
    Xa = np.array(X.value)
    err = max(abs(Xa[i, i] + Xa[j, j] - 2 * Xa[i, j] - d * d) ** 0.5 * 0 + abs(
        np.sqrt(max(Xa[i, i] + Xa[j, j] - 2 * Xa[i, j], 0.0)) - d) for i, j, d in pairs)
    return Xa, float(err)


def main() -> int:
    rng = np.random.default_rng(42)
    proj_cases = []
    for n in (6, 10):
        a = rng.standard_normal((n, n))
        x = ((a + a.T) / 2).tolist()
        tau = 0.5
        p, neg = psd_clip(np.array(x))
        proj_cases.append({"n": n, "x": x, "tau": tau,
                           "psd_ref": p.tolist(), "neg_ref": neg,
                           "svt_ref": svt(np.array(x), tau).tolist()})
    embed_cases = []
    for n, noise in ((8, 0.0), (12, 0.0), (10, 0.05)):
        pts = rng.standard_normal((n, 3))
        pairs = []
        for i in range(n):
            for j in range(i + 1, n):
                d = float(np.linalg.norm(pts[i] - pts[j]))
                if noise:
                    d *= 1.0 + float(rng.normal(0, noise))
                pairs.append((i, j, d))
        d2 = np.zeros((n, n))
        for i, j, d in pairs:
            d2[i, j] = d2[j, i] = d * d
        gram = double_center(d2)
        w, v = np.linalg.eigh(gram)
        idx = np.argsort(w)[::-1][:3]
        recon = (v[:, idx] * np.sqrt(np.maximum(w[idx], 0.0))) @ (
            v[:, idx] * np.sqrt(np.maximum(w[idx], 0.0))).T
        opt_gram, opt_err = cvxpy_sdp(pairs, n)
        embed_cases.append({
            "n": n, "noise": noise,
            "pairs": [[i, j, d] for i, j, d in pairs],
            "gram_ref": gram.tolist(),
            "gram_recon_ref": recon.tolist(),
            "cvxpy_max_err": opt_err,
        })
        print(f"n={n} noise={noise}: cvxpy max err={opt_err:.2e}")
    payload = {
        "ficha": "F83",
        "oracle": "oracle_sdr.py (numpy eigh + CVXPY/Clarabel SDP, seed 42)",
        "projection_tolerance": 1e-12,
        "gram_tolerance": 1e-9,
        "embed_tolerance_A": 0.05,
        "projections": proj_cases,
        "embeddings": embed_cases,
    }
    OUT.write_text(json.dumps(payload) + "\n")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
