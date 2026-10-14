#!/usr/bin/env python3
"""Oracle F91: KPM DOS/moments/Mulliken vs numpy eigh (no pybinding needed).

- Matrices: GOE (seed 42), N = 20 and 60.
- Reference eigenvalues from numpy.linalg.eigh (independent of our recurrence).
- Reference DOS: Gaussian broadening sigma = width/order on a fixed eV grid.
- Reference Mulliken (S = I, closed shell, Z = 1): P = V_occ V_occ^T,
  q_i = 1 - P_ii, with nocc = nelec/2.
- The Rust gate recomputes exact Chebyshev moments from these eigenvalues
  (textbook cos(k arccos)) using the expansion's own (a, b), so the moment
  comparison is bounds-independent.

Usage: oracle_kpm.py  (writes tests/fixtures/validation/kpm_numpy.json)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "tests" / "fixtures" / "validation" / "kpm_numpy.json"

ORDER = 500
NPTS = 401  # dx ~= sigma/3 at M=500 so the grid resolves the broadening
T_K = 1160.0  # finite-T gate: kT = 0.1 eV (code takes Kelvin: kT = T*8.617e-5)


def goe(rng: np.random.Generator, n: int) -> np.ndarray:
    a = rng.standard_normal((n, n))
    return (a + a.T) / 2.0


def gershgorin_halfwidth(h: np.ndarray) -> float:
    d = np.diag(h)
    r = np.abs(h).sum(axis=1) - np.abs(d)
    emin, emax = (d - r).min(), (d + r).max()
    m = (emax - emin) * 0.01
    return (emax - emin + 2 * m) / 2.0


def gauss_dos(grid: np.ndarray, w: np.ndarray, sigma: float) -> list[float]:
    g = np.exp(-0.5 * ((grid[:, None] - w[None, :]) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
    return g.sum(axis=1).tolist()


def gershgorin_bounds(h: np.ndarray) -> tuple[float, float]:
    d = np.diag(h)
    r = np.abs(h).sum(axis=1) - np.abs(d)
    emin, emax = (d - r).min(), (d + r).max()
    m = (emax - emin) * 0.01
    return emin - m, emax + m


def jackson_dos(grid: np.ndarray, w: np.ndarray, a: float, b: float, order: int) -> list[float]:
    """Exact DOS from eigenvalues through the Jackson-damped delta expansion.

    delta_M(x - xi) = [1 + 2 sum_{k>=1} g_k T_k(x) T_k(xi)] / (pi sqrt(1-x^2)),
    dos(E) = sum_i delta_M / a. Same textbook formula as dos_at_energy;
    the eigenvalues (and hence all moments) stay external.
    """
    m = order
    mp1 = m + 1
    cot = 1.0 / np.tan(np.pi / mp1)
    kk = np.arange(m)
    gk = ((mp1 - kk) * np.cos(np.pi * kk / mp1) + np.sin(np.pi * kk / mp1) * cot) / mp1
    x = np.clip((grid - b) / a, -1.0 + 1e-12, 1.0 - 1e-12)
    xi = np.clip((w - b) / a, -1.0 + 1e-12, 1.0 - 1e-12)
    Tx = np.cos(kk[:, None] * np.arccos(x)[None, :])
    Ti = np.cos(kk[:, None] * np.arccos(xi)[None, :])
    wgt = np.where(kk == 0, 1.0, 2.0)[:, None, None]
    delta = (gk[:, None, None] * wgt * Tx[:, :, None] * Ti[:, None, :]).sum(axis=0)
    delta /= np.pi * np.sqrt(1 - x**2)[:, None]
    return (delta.sum(axis=1) / a).tolist()


def main() -> int:
    rng = np.random.default_rng(42)
    systems = []
    for n in (20, 60):
        h = goe(rng, n)
        w, vecs = np.linalg.eigh(h)
        lo, hi = float(w[0]), float(w[-1])
        width = hi - lo
        # Jackson ~= Gaussian with sigma = pi*a/order, where a is the
        # rescaling half-width (Gershgorin + 1% margin, mirrored formula).
        a = gershgorin_halfwidth(h)
        sigma = np.pi * a / ORDER
        grid = np.linspace(lo - 0.05 * width, hi + 0.05 * width, NPTS)
        dos = gauss_dos(grid, w, sigma)
        emin_g, emax_g = gershgorin_bounds(h)
        a_g = (emax_g - emin_g) / 2.0
        b_g = (emax_g + emin_g) / 2.0
        dos_jackson = jackson_dos(grid, w, a_g, b_g, ORDER)
        # Mulliken reference: half-filled closed shell, S = I, Z_i = 1.
        nelec = n
        nocc = nelec // 2
        # Closed shell: P = 2 V_occ V_occ^T (spin degeneracy).
        p = 2.0 * vecs[:, :nocc] @ vecs[:, :nocc].T
        charges = (1.0 - np.diag(p)).tolist()
        # Finite-T reference (kT = 0.1 eV): mu bisected for nelec/2
        # spatial states, P = 2 V diag(f) V^T.
        kt = T_K * 8.617333262e-5
        mu_lo, mu_hi = float(w[0]) - 1.0, float(w[-1]) + 1.0
        for _ in range(100):
            mu = (mu_lo + mu_hi) / 2
            if float(np.sum(1.0 / (1.0 + np.exp((w - mu) / kt)))) < nelec / 2:
                mu_lo = mu
            else:
                mu_hi = mu
        mu = (mu_lo + mu_hi) / 2
        f = 1.0 / (1.0 + np.exp((w - mu) / kt))
        p_fd = 2.0 * (vecs * f) @ vecs.T
        charges_fd = (1.0 - np.diag(p_fd)).tolist()
        systems.append({
            "n": n, "matrix": h.tolist(), "eigenvalues": w.tolist(),
            "order": ORDER, "grid_eV": grid.tolist(), "sigma_eV": float(sigma),
            "gershgorin_a_eV": float(a), "gershgorin_b_eV": float(b_g),
            "dos_ref": dos, "dos_jackson_ref": dos_jackson, "nelec": nelec,
            "charges_ref": charges,
            "temperature_K": T_K, "mu_eV": float(mu),
            "charges_fd_ref": charges_fd,
        })
    # Gapped spectrum: T = 0 step sits mid-gap (molecular-like case).
    n = 20
    lam = np.concatenate([np.linspace(-5, -1, 10), np.linspace(1, 5, 10)])
    q, _ = np.linalg.qr(rng.standard_normal((n, n)))
    h = (q * lam) @ q.T
    h = (h + h.T) / 2
    w, vecs = np.linalg.eigh(h)
    p = 2.0 * vecs[:, :10] @ vecs[:, :10].T
    systems.append({
        "n": n, "matrix": h.tolist(), "eigenvalues": w.tolist(),
        "order": ORDER, "grid_eV": [], "sigma_eV": 0.0,
        "gershgorin_a_eV": float(gershgorin_halfwidth(h)),
        "dos_ref": [], "dos_jackson_ref": [], "nelec": n,
        "charges_ref": (1.0 - np.diag(p)).tolist(),
        "temperature_K": 0.0, "mu_eV": 0.0,
        "charges_fd_ref": [],
        "gapped": True,
    })
    payload = {
        "ficha": "F91",
        "oracle": "oracle_kpm.py (numpy eigh + Gaussian broadening, seed 42)",
        "dos_tolerance_l2_rel": 0.05,
        "moment_tolerance": 1e-9,
        "charge_tolerance_abs": 0.05,
        "systems": systems,
    }
    OUT.write_text(json.dumps(payload) + "\n")
    print(f"wrote {OUT} ({len(systems)} systems)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
