#!/usr/bin/env python3
"""Oracle F86: Gouy-Chapman / GCS / Helmholtz EDL vs analytic electrochemistry.

Independent numpy implementation of the textbook closed forms
(Grahame equation, Debye-Huckel kappa, Stern series capacitance):
- GC: phi(x) = (2kT/ze) ln((1+g e^-kx)/(1-g e^-kx)), g = tanh(ze phi0/4kT)
      C_diff = eps0 eps_r kappa cosh(ze phi0/2kT)
- GCS: phi_d = (2kT/ze) asinh(sigma/sqrt(8 c0 eps0 eps_r kT)) [Grahame
      inverted], phi_c = sigma/C_H, C^-1 = C_H^-1 + C_diff(phi_d)^-1
- Helmholtz: C_H = eps0 eps_c/d, phi drop = sigma/C_H
- scan: total C over a sigma grid (must be even, minimum at sigma = 0)

Constants: SI (eps0, kB, e, Na), symmetric 1:1 electrolyte, c0 = I*1000 mol/m^3.

Usage: oracle_edl.py  (writes tests/fixtures/validation/edl_analytic.json)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "tests" / "fixtures" / "validation" / "edl_analytic.json"

EPS0 = 8.8541878128e-12
KB = 1.380649e-23
E = 1.602176634e-19
NA = 6.02214076e23
ANG = 1e-10


def kappa_debye(I_M: float, T: float, eps_b: float) -> float:
    c0 = I_M * 1000.0
    return float(np.sqrt(2 * c0 * NA * E**2 / (EPS0 * eps_b * KB * T)))


def gc_profile(phi0: float, T: float, I_M: float, eps_b: float, n: int, extent_A: float) -> dict:
    z = 1
    Vt = KB * T / E
    kap = kappa_debye(I_M, T, eps_b)
    g = float(np.tanh(E * phi0 / (4 * KB * T)))
    xs = np.linspace(0.0, extent_A, n)
    phi, fld, rho = [], [], []
    c0 = I_M * 1000.0
    for xA in xs:
        xm = xA * ANG
        ge = g * np.exp(-kap * xm)
        p = 2 * Vt * np.log((1 + ge) / (1 - ge)) if abs(ge) < 1 - 1e-12 else phi0
        den = 1 - g * g * np.exp(-2 * kap * xm)
        f = 2 * Vt * kap * 2 * g * np.exp(-kap * xm) / den if abs(den) > 1e-15 else 0.0
        r = -2 * E * c0 * NA * np.sinh(E * p / (KB * T))
        phi.append(float(p))
        fld.append(float(f))
        rho.append(float(r))
    c_diff = EPS0 * eps_b * kap * np.cosh(E * phi0 / (2 * KB * T))
    return {
        "x_A": xs.tolist(), "phi_V": phi, "field_V_per_m": fld, "rho_C_per_m3": rho,
        "c_diffuse_F_per_m2": float(c_diff), "kappa_per_m": kap,
    }


def grahame_sigma(phi_d: float, T: float, I_M: float, eps_b: float) -> float:
    c0 = I_M * 1000.0
    return float(np.sqrt(8 * c0 * NA * EPS0 * eps_b * KB * T) * np.sinh(E * phi_d / (2 * KB * T)))


def gcs_profile(sigma: float, T: float, I_M: float, eps_b: float, eps_c: float,
                stern_A: float, n: int, extent_A: float) -> dict:
    z = 1
    Vt = KB * T / E
    stern_m = stern_A * ANG
    cH = EPS0 * eps_c / stern_m
    phi_c = sigma / cH
    c0 = I_M * 1000.0
    pre = float(np.sqrt(8 * c0 * NA * EPS0 * eps_b * KB * T))
    phi_d = 2 * Vt * np.arcsinh(sigma / pre)
    kap = kappa_debye(I_M, T, eps_b)
    g = float(np.tanh(E * phi_d / (4 * KB * T)))
    e_compact = sigma / (EPS0 * eps_c)
    xs = np.linspace(0.0, extent_A, n)
    phi, fld, rho = [], [], []
    for xA in xs:
        xm = xA * ANG
        if xA <= stern_A:
            phi.append(float(e_compact * xm))
            fld.append(float(e_compact))
            rho.append(float(sigma / stern_m))
        else:
            xd = xm - stern_m
            ge = g * np.exp(-kap * xd)
            p = 2 * Vt * np.log((1 + ge) / (1 - ge)) if abs(ge) < 1 - 1e-12 else phi_d
            den = 1 - g * g * np.exp(-2 * kap * xd)
            f = 2 * Vt * kap * 2 * g * np.exp(-kap * xd) / den if abs(den) > 1e-15 else 0.0
            r = -2 * E * c0 * NA * np.sinh(E * p / (KB * T))
            phi.append(float(phi_c + p))
            fld.append(float(f))
            rho.append(float(r))
    c_diff = EPS0 * eps_b * kap * np.cosh(E * phi_d / (2 * KB * T))
    c_tot = 1.0 / (1.0 / cH + 1.0 / c_diff)
    return {
        "x_A": xs.tolist(), "phi_V": phi, "field_V_per_m": fld, "rho_C_per_m3": rho,
        "phi_compact_V": float(phi_c), "phi_diffuse_V": float(phi_d),
        "c_compact_F_per_m2": float(cH), "c_diffuse_F_per_m2": float(c_diff),
        "c_total_F_per_m2": float(c_tot),
    }


def helmholtz(sigma: float, eps_c: float, stern_A: float) -> dict:
    stern_m = stern_A * ANG
    cH = EPS0 * eps_c / stern_m
    return {"c_compact_F_per_m2": float(cH), "drop_V": float(sigma / cH)}


def main() -> int:
    N, EXT = 65, 50.0
    gc_cases, gcs_cases, h_cases, scans = [], [], [], []
    for T in (283.15, 298.15, 310.0):
        for I_M in (0.01, 0.1, 1.0):
            for phi0 in (-0.25, -0.1, -0.025, 0.025, 0.1, 0.25):
                p = gc_profile(phi0, T, I_M, 78.5, N, EXT)
                gc_cases.append({"T_K": T, "I_M": I_M, "eps_b": 78.5, "phi0_V": phi0,
                                 "n": N, "extent_A": EXT, **p,
                                 "sigma_grahame_C_per_m2": grahame_sigma(phi0, T, I_M, 78.5)})
            for sigma in (-0.2, -0.05, -0.01, 0.01, 0.05, 0.2):
                g = gcs_profile(sigma, T, I_M, 78.5, 6.0, 3.0, N, EXT)
                gcs_cases.append({"T_K": T, "I_M": I_M, "eps_b": 78.5, "eps_c": 6.0,
                                  "stern_A": 3.0, "sigma_C_per_m2": sigma,
                                  "n": N, "extent_A": EXT, **g})
                h_cases.append({"sigma_C_per_m2": sigma, "eps_c": 6.0, "stern_A": 3.0,
                                **helmholtz(sigma, 6.0, 3.0)})
    for T in (298.15,):
        for I_M in (0.1, 1.0):
            grid = np.linspace(-0.2, 0.2, 17).tolist()
            ctots = [gcs_profile(s, T, I_M, 78.5, 6.0, 3.0, N, EXT)["c_total_F_per_m2"] for s in grid]
            scans.append({"T_K": T, "I_M": I_M, "eps_b": 78.5, "eps_c": 6.0, "stern_A": 3.0,
                          "sigma_grid": grid, "c_total_F_per_m2": ctots})
    payload = {
        "ficha": "F86",
        "oracle": "oracle_edl.py (textbook Grahame/Stern closed forms, SI)",
        "tolerance_rel": 1e-9,
        "gc": gc_cases, "gcs": gcs_cases, "helmholtz": h_cases, "scans": scans,
    }
    OUT.write_text(json.dumps(payload) + "\n")
    print(f"wrote {OUT} ({len(gc_cases)} gc, {len(gcs_cases)} gcs, {len(scans)} scans)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
