#!/usr/bin/env python3
"""Oracle F87: HTST/Eyring + microkinetics vs closed form + scipy.

- Eyring k = pref e^(-dG#/kT), pref default kT/h, K = e^(-dG/kT),
  k_rev = k_fwd/K — independent numpy implementation with CODATA
  constants (catches eV/J slips, wrong h/kB).
- A<->B reversible + A->B->C consecutive networks: scipy solve_ivp
  (RK45, rtol=1e-10) mass-action references at shared frame times.
- A<->B equilibrium: exact cB/cA = K (thermodynamic consistency).
- Cantera cross-check is the unmet promotion criterion (needs full
  mechanism/species thermo I/O the implementation lacks).

Usage: oracle_kinetics.py  (writes tests/fixtures/validation/kinetics_ref.json)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "tests" / "fixtures" / "validation" / "kinetics_ref.json"

KB_EV = 8.617333262145e-5
H_EV_S = 4.135667696e-15


def eyring(dg_ts: float, dg_rxn: float, T: float, pref: float | None = None) -> dict:
    p = KB_EV * T / H_EV_S if pref is None else pref
    kf = p * np.exp(-dg_ts / (KB_EV * T))
    K = np.exp(-dg_rxn / (KB_EV * T))
    return {"kf": float(kf), "kr": float(kf / K), "K": float(K), "pref": float(p)}


def main() -> int:
    htst = []
    for T in (300.0, 500.0, 800.0, 1000.0):
        for dg_ts in (0.1, 0.5, 1.0, 1.5, 2.0):
            for dg_rxn in (-0.5, 0.0, 0.5):
                htst.append({"T_K": T, "dg_ts_eV": dg_ts, "dg_rxn_eV": dg_rxn,
                             "pref_override": None, **eyring(dg_ts, dg_rxn, T)})
    for T, dg_ts, dg_rxn, pref in ((300.0, 0.8, -0.2, 1e13), (800.0, 1.2, 0.3, 5e12)):
        htst.append({"T_K": T, "dg_ts_eV": dg_ts, "dg_rxn_eV": dg_rxn,
                     "pref_override": pref, **eyring(dg_ts, dg_rxn, T, pref)})
    networks = []
    # A <-> B reversible (k from Eyring at 500 K).
    T = 500.0
    r = eyring(0.9, -0.2, T)
    kf, kr = r["kf"], r["kr"]
    t_end = 5.0 / min(kf, kr)
    times = np.linspace(0, t_end, 11)
    sol = solve_ivp(lambda t, c: [-kf * c[0] + kr * c[1], kf * c[0] - kr * c[1]],
                    (0, t_end), [1.0, 0.0], t_eval=times, method="RK45",
                    rtol=1e-10, atol=1e-13)
    assert sol.success
    networks.append({
        "id": "AB-reversible", "species": ["A", "B"],
        "steps": [{"step_id": "f", "dg_ts_eV": 0.9, "dg_rxn_eV": -0.2, "pref": None}],
        "stoichiometry": [[(0, -1.0), (1, 1.0)]],
        "initial": [1.0, 0.0], "T_K": T, "t_end_s": t_end,
        "times_s": times.tolist(), "populations": sol.y.T.tolist(),
        "K_ref": r["K"],
    })
    # A -> B -> C consecutive, second step slower.
    r1 = eyring(0.7, -0.4, T)
    r2 = eyring(1.0, -0.3, T)
    t_end = 5.0 / min(r1["kf"], r2["kf"])
    times = np.linspace(0, t_end, 11)
    def rhs(t, c):
        r_a = r1["kf"] * c[0]
        r_b = r2["kf"] * c[1]
        return [-r_a, r_a - r_b, r_b]
    sol = solve_ivp(rhs, (0, t_end), [1.0, 0.0, 0.0], t_eval=times,
                    method="RK45", rtol=1e-10, atol=1e-13)
    assert sol.success
    networks.append({
        "id": "ABC-consecutive", "species": ["A", "B", "C"],
        "steps": [{"step_id": "s1", "dg_ts_eV": 0.7, "dg_rxn_eV": -0.4, "pref": None},
                  {"step_id": "s2", "dg_ts_eV": 1.0, "dg_rxn_eV": -0.3, "pref": None}],
        "stoichiometry": [[(0, -1.0), (1, 1.0)], [(1, -1.0), (2, 1.0)]],
        "initial": [1.0, 0.0, 0.0], "T_K": T, "t_end_s": t_end,
        "times_s": times.tolist(), "populations": sol.y.T.tolist(),
        "K_ref": None,
    })
    payload = {
        "ficha": "F87",
        "oracle": "oracle_kinetics.py (Eyring closed form + scipy solve_ivp)",
        "htst_tolerance_rel": 1e-12,
        # Euler (k*dt<=0.1 substepped) vs RK45 integrator difference dominates.
        "trace_tolerance_abs": 1e-3,
        "htst": htst,
        "networks": networks,
    }
    OUT.write_text(json.dumps(payload) + "\n")
    print(f"wrote {OUT} ({len(htst)} htst, {len(networks)} networks)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
