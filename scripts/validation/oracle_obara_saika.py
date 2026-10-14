#!/usr/bin/env python3
"""Oracle F80: Boys function + (ss|ss) ERIs + Schwarz bounds.

References (independent implementations):
- Boys F_n(x) = 1/2 * x^-(n+1/2) * gamma(n+1/2, x) via scipy.special
  (regularized gammainc * gamma) — no series, no erf approximation.
- (ss|ss) primitive ERIs via PySCF/libcint: a Mole with 4 distinct
  elements carrying single-primitive s functions, unit Bohr, so
  intor('int2e')[0,1,2,3] is exactly the primitive integral our
  `eri_ssss` must reproduce. Units: bohr / bohr^-2 / Hartree.

Usage: oracle_obara_saika.py  (writes tests/fixtures/validation/obara_saika.json)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from pyscf import gto
from scipy.special import gamma, gammainc

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "tests" / "fixtures" / "validation" / "obara_saika.json"

ELEMS = ["H", "He", "Li", "Be"]


def boys_ref(n: int, x: float) -> float:
    if x <= 0.0:
        return 1.0 / (2 * n + 1)
    a = n + 0.5
    return 0.5 * x**-a * gammainc(a, x) * gamma(a)


def eri_prim(
    alpha: float,
    ca: list[float],
    beta: float,
    cb: list[float],
    gamma_: float,
    cc: list[float],
    delta: float,
    cd: list[float],
) -> float:
    mol = gto.Mole()
    mol.atom = (
        f"{ELEMS[0]} {ca[0]} {ca[1]} {ca[2]}; "
        f"{ELEMS[1]} {cb[0]} {cb[1]} {cb[2]}; "
        f"{ELEMS[2]} {cc[0]} {cc[1]} {cc[2]}; "
        f"{ELEMS[3]} {cd[0]} {cd[1]} {cd[2]}"
    )
    mol.basis = {
        ELEMS[0]: [[0, [alpha, 1.0]]],
        ELEMS[1]: [[0, [beta, 1.0]]],
        ELEMS[2]: [[0, [gamma_, 1.0]]],
        ELEMS[3]: [[0, [delta, 1.0]]],
    }
    mol.unit = "Bohr"
    mol.build()
    raw = float(mol.intor("int2e")[0, 1, 2, 3])
    # PySCF/libcint normalizes contracted basis functions; our `eri_ssss`
    # uses raw primitives. Unnormalize: N(a) = (2a/pi)^(3/4) per primitive.
    norm = float(np.prod([(2.0 * e / np.pi) ** 0.75 for e in (alpha, beta, gamma_, delta)]))
    return raw / norm


def main() -> int:
    rng = np.random.default_rng(42)

    boys = []
    for n in range(7):
        for x in [0.0, 1e-13, 1e-8, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 29.9, 30.0, 30.1, 50.0, 100.0]:
            boys.append({"n": n, "x": x, "ref": boys_ref(n, x)})

    eris = []
    for _ in range(12):
        exps = [float(v) for v in 10.0 ** rng.uniform(np.log10(0.2), np.log10(5.0), 4)]
        ctrs = [list(map(float, rng.uniform(-2.0, 2.0, 3))) for _ in range(4)]
        ref = eri_prim(exps[0], ctrs[0], exps[1], ctrs[1], exps[2], ctrs[2], exps[3], ctrs[3])
        eris.append(
            {
                "alpha": exps[0], "ca": ctrs[0],
                "beta": exps[1], "cb": ctrs[1],
                "gamma": exps[2], "cc": ctrs[2],
                "delta": exps[3], "cd": ctrs[3],
                "ref": ref,
            }
        )

    # Same-center self integral (ij|ij) for Schwarz references.
    schwarz = []
    for _ in range(8):
        a, b = [float(v) for v in 10.0 ** rng.uniform(np.log10(0.2), np.log10(5.0), 2)]
        ca = list(map(float, rng.uniform(-2.0, 2.0, 3)))
        cb = list(map(float, rng.uniform(-2.0, 2.0, 3)))
        hii = eri_prim(a, ca, b, cb, a, ca, b, cb)
        schwarz.append({"alpha": a, "ca": ca, "beta": b, "cb": cb, "bound_ref": float(np.sqrt(hii))})

    payload = {
        "ficha": "F80",
        "units": {"length": "bohr", "exponent": "1/bohr^2", "energy": "hartree"},
        "oracle": "oracle_obara_saika.py (scipy gammainc + PySCF/libcint int2e, seed 42)",
        "boys_tolerance_rel": 1e-12,
        "eri_tolerance_rel": 1e-10,
        "boys": boys,
        "eris": eris,
        "schwarz": schwarz,
    }
    OUT.write_text(json.dumps(payload, indent=1) + "\n")
    print(f"wrote {OUT} ({len(boys)} boys, {len(eris)} eris, {len(schwarz)} schwarz)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
