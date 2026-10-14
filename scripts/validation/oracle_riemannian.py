#!/usr/bin/env python3
"""Oracle F94: PSD cone ops + Gram-matrix distance geometry vs numpy/scipy.

- psd_projection: numpy eigh clip (independent implementation).
- psd_distance: numpy Frobenius norm. NOTE: our psd_distance is the
  Frobenius proxy, NOT the affine-invariant Riemannian distance
  ||log(A^-1/2 B A^-1/2)||_F (documented gap; geomstats comparison is
  the unmet promotion criterion).
- psd_retraction: P_+(X + xi) via numpy.
- objective/gradient of the Gram-matrix constraint penalty: independent
  numpy implementation of the same formulas.
- optimizer: distance-geometry problems from true 3D point sets with
  5% bound slack; scipy.optimize.minimize (L-BFGS-B over vech) reference.

Usage: oracle_riemannian.py  (writes tests/fixtures/validation/riemannian_numpy.json)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "tests" / "fixtures" / "validation" / "riemannian_numpy.json"


def proj_psd(x: np.ndarray) -> np.ndarray:
    w, v = np.linalg.eigh(x)
    return (v * np.maximum(w, 0.0)) @ v.T


def gram_objective(X: np.ndarray, cons: list[tuple[int, int, float, float]]) -> float:
    f = 0.0
    for i, j, lo, hi in cons:
        d2 = X[i, i] + X[j, j] - 2 * X[i, j]
        if d2 < lo * lo:
            f += (lo * lo - d2) ** 2
        elif d2 > hi * hi:
            f += (d2 - hi * hi) ** 2
    return float(f)


def gram_gradient(X: np.ndarray, cons: list[tuple[int, int, float, float]]) -> np.ndarray:
    g = np.zeros_like(X)
    for i, j, lo, hi in cons:
        d2 = X[i, i] + X[j, j] - 2 * X[i, j]
        if d2 < lo * lo:
            fac = -2.0 * (lo * lo - d2)
        elif d2 > hi * hi:
            fac = 2.0 * (d2 - hi * hi)
        else:
            continue
        g[i, i] += fac
        g[j, j] += fac
        g[i, j] -= 2 * fac
        g[j, i] -= 2 * fac
    return (g + g.T) * 0.5


def vech_index(n: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i, n)]


def main() -> int:
    rng = np.random.default_rng(42)
    proj_cases = []
    for n in (8, 16):
        for _ in range(3):
            a = rng.standard_normal((n, n))
            x = ((a + a.T) / 2).tolist()
            e = rng.standard_normal((n, n)) * 0.5
            xi = ((e + e.T) / 2).tolist()
            xa = np.array(x)
            proj_cases.append({
                "n": n, "x": x, "xi": xi,
                "projection_ref": proj_psd(xa).tolist(),
                "retraction_ref": proj_psd(xa + np.array(xi)).tolist(),
            })
    dist_cases = []
    for n in (8, 16):
        a = rng.standard_normal((n, n))
        x = (a + a.T) / 2
        b = rng.standard_normal((n, n))
        y = (b + b.T) / 2
        dist_cases.append({
            "n": n, "a": x.tolist(), "b": y.tolist(),
            "frobenius_ref": float(np.linalg.norm(x - y)),
        })
    eval_cases = []
    for n in (5, 8):
        pts = rng.standard_normal((n, 3))
        xa = rng.standard_normal((n, n)) * 2.0
        x = (xa + xa.T) / 2
        cons = []
        for i in range(n):
            for j in range(i + 1, n):
                d = float(np.linalg.norm(pts[i] - pts[j]))
                cons.append((i, j, 0.95 * d, 1.05 * d))
        eval_cases.append({
            "n": n, "x": x.tolist(),
            "constraints": [[i, j, lo, hi] for i, j, lo, hi in cons],
            "objective_ref": gram_objective(x, cons),
            "gradient_ref": gram_gradient(x, cons).tolist(),
        })
    opt_cases = []
    for n in (5, 8):
        pts = rng.standard_normal((n, 3))
        g = pts @ pts.T
        gn = g + rng.standard_normal((n, n)) * 0.05
        gn = ((gn + gn.T) / 2).tolist()
        cons = []
        for i in range(n):
            for j in range(i + 1, n):
                d = float(np.linalg.norm(pts[i] - pts[j]))
                cons.append((i, j, 0.95 * d, 1.05 * d))
        idx = vech_index(n)
        x0 = np.array([gn[i][j] for i, j in idx])

        def fun(v: np.ndarray) -> float:
            X = np.zeros((n, n))
            for (i, j), val in zip(idx, v):
                X[i, j] = X[j, i] = val
            return gram_objective(X, cons)
        def jac(v: np.ndarray) -> np.ndarray:
            X = np.zeros((n, n))
            for (i, j), val in zip(idx, v):
                X[i, j] = X[j, i] = val
            G = gram_gradient(X, cons)
            return np.array([G[i, j] * (1.0 if i == j else 2.0) for i, j in idx])
        res = minimize(fun, x0, jac=jac, method="L-BFGS-B",
                       options={"maxiter": 2000, "ftol": 1e-14, "gtol": 1e-12})
        X = np.zeros((n, n))
        for (i, j), val in zip(idx, res.x):
            X[i, j] = X[j, i] = val
        viol = 0.0
        for i, j, lo, hi in cons:
            d = float(np.sqrt(max(X[i, i] + X[j, j] - 2 * X[i, j], 0.0)))
            viol = max(viol, lo - d if d < lo else (d - hi if d > hi else 0.0))
        opt_cases.append({
            "n": n, "initial": gn,
            "constraints": [[i, j, lo, hi] for i, j, lo, hi in cons],
            "scipy_objective": float(res.fun), "scipy_max_violation_A": viol,
        })
        print(f"opt N={n}: scipy obj={res.fun:.3e} viol={viol:.3e}")
    payload = {
        "ficha": "F94",
        "oracle": "oracle_riemannian.py (numpy eigh + scipy L-BFGS-B, seed 42)",
        "note": ("psd_distance is the Frobenius proxy, not the affine-invariant "
                 "geodesic distance (documented gap; geomstats cross-check is the "
                 "unmet promotion criterion)."),
        "projection_tolerance": 1e-10,
        "distance_tolerance": 1e-12,
        "objective_tolerance": 1e-9,
        "optimizer_violation_tolerance_A": 1e-6,
        "projections": proj_cases,
        "distances": dist_cases,
        "eval": eval_cases,
        "optimizer": opt_cases,
    }
    OUT.write_text(json.dumps(payload) + "\n")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
