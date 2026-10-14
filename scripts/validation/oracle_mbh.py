#!/usr/bin/env python3
"""Oracle F92: MBH harmonic frequencies vs SI closed form.

Diatomic E = k/2 (r - r0)^2 (k in kcal/mol/A^2) -> nu = (2 pi c)^-1
sqrt(k_SI / mu_SI), derived from SI constants (also validates the
108.59 conversion factor in the implementation). Expected spectrum:
5 near-zero modes (3 translations + 2 rotations) + 1 stretch.

Usage: oracle_mbh.py  (writes tests/fixtures/validation/mbh_analytic.json)
"""

from __future__ import annotations

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "tests" / "fixtures" / "validation" / "mbh_analytic.json"

KCAL = 4184.0
NA = 6.02214076e23
AMU = 1.66053906660e-27
ANG = 1e-10
C = 2.99792458e10  # cm/s


def stretch_cm(m1_amu: float, m2_amu: float, k_kcal: float) -> float:
    mu = m1_amu * m2_amu / (m1_amu + m2_amu) * AMU
    k_si = k_kcal * KCAL / NA / ANG**2
    return float((k_si / mu) ** 0.5 / (2 * 3.141592653589793 * C))


def main() -> int:
    cases = []
    for name, m1, m2, k in (
        ("H2-like", 1.008, 1.008, 100.0),
        ("N2-like", 14.007, 14.007, 300.0),
        ("CO-like", 12.011, 15.999, 250.0),
        ("HCl-like", 1.008, 35.453, 120.0),
    ):
        nu = stretch_cm(m1, m2, k)
        cases.append({"name": name, "m1_amu": m1, "m2_amu": m2,
                      "k_kcal_per_mol_A2": k, "r0_A": 1.0, "nu_ref_cm": nu})
        print(f"{name}: nu={nu:.4f} cm-1")
    payload = {
        "ficha": "F92",
        "oracle": "oracle_mbh.py (SI closed form; validates the 108.59 factor)",
        "note": ("Analytic harmonic core only. Experimental/literature +-2 cm-1 "
                 "(Ghysels et al.) needs a real PES + anharmonicity: unmet promotion "
                 "criterion, documented gap."),
        "stretch_tolerance_rel": 1e-3,
        # Zero modes: translations are exact; rotations of linear molecules
        # pick up quartic (chord-shortening) leakage into the fd Hessian.
        "zero_tolerance_cm": 5.0,
        "cases": cases,
    }
    OUT.write_text(json.dumps(payload, indent=1) + "\n")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
