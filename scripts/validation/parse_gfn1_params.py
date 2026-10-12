#!/usr/bin/env python3
"""Extract GFN1-xTB parameters from tblite source for sci-form validation (F44).

Reads tblite's `src/tblite/xtb/gfn1.f90` (LGPL-3.0-or-later; values reproduced
here for cross-validation against tblite GFN1-xTB) and emits
`tests/fixtures/validation/gfn1_params.json`.

Usage:
    python scripts/validation/parse_gfn1_params.py [--src URL-or-path]
"""

from __future__ import annotations

import json
import re
import sys
import time
import urllib.request
from pathlib import Path

DEFAULT_URL = (
    "https://raw.githubusercontent.com/tblite/tblite/main/src/tblite/xtb/gfn1.f90"
)
LOCAL_FALLBACK = Path("/tmp/opencode/tblite_gfn1.f90")


def load_source(arg: str | None) -> str:
    if arg:
        p = Path(arg)
        if p.exists():
            return p.read_text()
        with urllib.request.urlopen(arg, timeout=60) as r:
            return r.read().decode()
    if LOCAL_FALLBACK.exists():
        return LOCAL_FALLBACK.read_text()
    with urllib.request.urlopen(DEFAULT_URL, timeout=60) as r:
        return r.read().decode()


def floats(block: str) -> list[float]:
    return [float(t) for t in re.findall(r"(-?\d+\.\d+(?:e[+-]?\d+)?)_wp", block)]


# tblite GFN1 uses an older eV -> Hartree conversion "for consistency"
# (`evtoau = 1.0_wp / 27.21138505_wp` in src/tblite/xtb/gfn1.f90).
EVTOAU = 1.0 / 27.21138505


def table1d(src: str, name: str, n: int) -> list[float]:
    m = re.search(name + r"\(max_elem\) = \[&(.*?)\]", src, re.S)
    assert m, name
    vals = floats(m.group(1))
    assert len(vals) == n, (name, len(vals))
    return vals


def table2d(src: str, name: str, rows: int, cols: int, scale: float = 1.0) -> list[list[float]]:
    # Fortran reshape fills column-major: flat[(j)*rows + i].
    # Reshape-close is `],&` + continuation; anchor on shape(name) instead.
    pat = re.escape(name) + r"(.*?)(.*?)shape\(" + re.escape(name) + r"\)"
    m = re.search(pat, src, re.S)
    assert m, name
    blk = m.group(2).split("[&", 1)[1]
    vals = floats(blk)
    assert len(vals) == rows * cols, (name, len(vals))
    return [[vals[j * rows + i] * scale for j in range(cols)] for i in range(rows)]


def scalar(src: str, name: str) -> float:
    m = re.search(r"parameter :: " + name + r" = ([-\d.e+_wp]+)", src)
    assert m, name
    return float(m.group(1).replace("_wp", ""))


def main() -> int:
    arg = sys.argv[1] if len(sys.argv) > 1 else None
    src = load_source(arg)

    out: dict = {
        "source": "tblite main src/tblite/xtb/gfn1.f90 (LGPL-3.0-or-later)",
        "generated": time.strftime("%Y-%m-%d"),
        "globals": {
            "s6": 1.0, "s8": 2.4, "a1": 0.63, "a2": 5.0, "s9": 0.0,
            "gexp": 2.0, "enscale": -7.0e-3, "kdiff": 2.85,
            "kdiag": [1.85, 2.25, 2.0, 2.0, 2.0],
            "halogen_radscale": 1.3, "halogen_damping": 0.44,
            "rep_kexp": 1.5, "rep_rexp": 1.0,
        },
        "elements": {},
    }
    # integer 1-D table (no _wp suffix)
    m = re.search(r"nshell\(max_elem\) = \[&(.*?)\]", src, re.S)
    assert m, "nshell"
    nshell = [int(t) for t in re.findall(r"(\d+),?", m.group(1).replace("&", " ")) if t != ""]
    assert len(nshell) == 86, len(nshell)
    rep_alpha = table1d(src, "rep_alpha", 86)
    rep_zeff = table1d(src, "rep_zeff", 86)
    hubbard = table1d(src, "hubbard_parameter", 86)
    # The Fortran table is declared as `[...] * 0.1_wp`; the trailing scale
    # lives outside the bracketed block the regex captures.
    hubbard_derivs = [0.1 * v for v in table1d(src, "p_hubbard_derivs", 86)]
    m = re.search(r"halogen_bond\(1:max_elem\)\s*=\s*0\.1_wp\s*\*\s*\[\s*&(.*?)\]", src, re.S)
    assert m, "halogen_bond"
    halogen_bond = [0.1 * float(t) for t in re.findall(r"(-?\d+\.\d+(?:e[+-]\d+)?)_wp", m.group(1))]
    assert len(halogen_bond) == 86
    m = re.search(r"shell_hubbard\(0:2, max_elem\)(.*?)shape\(shell_hubbard\)", src, re.S)
    assert m, "shell_hubbard"
    _blk = m.group(1).split("[&", 1)[1]
    _sh = [float(t) for t in re.findall(r"(-?\d+\.\d+(?:e[+-]?\d+)?)_wp", _blk)]
    assert len(_sh) == 258, len(_sh)
    shell_hub = [[1.0 + _sh[j * 3 + i] for j in range(86)] for i in range(3)]
    # integer tables: parse separately (no _wp suffix)
    def inttable(name: str, rows: int, cols: int):
        m = re.search(re.escape(name) + r"(.*?)(.*?)shape\(" + re.escape(name) + r"\)", src, re.S)
        assert m, name
        toks = re.findall(r"(-?\d+),?", m.group(2).replace("&", " "))
        vals = [int(t) for t in toks if t not in ("", "-")]
        assert len(vals) == rows * cols, (name, len(vals))
        return [[vals[j * rows + i] for j in range(cols)] for i in range(rows)]

    ang = inttable("ang_shell", 3, 86)
    pqn = inttable("principal_quantum_number", 3, 86)
    ngauss = inttable("number_of_primitives", 3, 86)
    slater = table2d(src, "slater_exponent", 3, 86)
    selfenergy = table2d(src, "p_selfenergy", 3, 86, scale=EVTOAU)
    kcn = table2d(src, "p_kcn", 3, 86, scale=EVTOAU * 0.01)
    shpoly = table2d(src, "p_shpoly", 3, 86, scale=0.01)
    refocc = table2d(src, "reference_occ", 3, 86)
    for z in range(1, 87):
        i = z - 1
        # tables above are [row][Z]; transpose to per-element shell lists
        out["elements"][str(z)] = {
            "nshell": nshell[i],
            "ang": [ang[r][i] for r in range(3)][: nshell[i]],
            "pqn": [pqn[r][i] for r in range(3)][: nshell[i]],
            "ngauss": [ngauss[r][i] for r in range(3)][: nshell[i]],
            "slater": [slater[r][i] for r in range(3)][: nshell[i]],
            "selfenergy": [selfenergy[r][i] for r in range(3)][: nshell[i]],
            "kcn": [kcn[r][i] for r in range(3)][: nshell[i]],
            "shpoly": [shpoly[r][i] for r in range(3)][: nshell[i]],
            "refocc": [refocc[r][i] for r in range(3)][: nshell[i]],
            "shell_hubbard": [shell_hub[r][i] for r in range(3)][: nshell[i]],
            "rep_alpha": rep_alpha[i],
            "rep_zeff": rep_zeff[i],
            "hubbard": hubbard[i],
            "hubbard_derivs": hubbard_derivs[i],
            "halogen_bond": halogen_bond[i],
        }
    outp = Path("tests/fixtures/validation/gfn1_params.json")
    outp.write_text(json.dumps(out, indent=1) + "\n")
    e = out["elements"]
    print(f"elements 1..86 -> {outp}")
    print("H:", json.dumps(e["1"])[:160])
    print("C:", json.dumps(e["6"])[:160])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
