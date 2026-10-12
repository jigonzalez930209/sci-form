#!/usr/bin/env python3
"""Parse xtb's GFN0 parameter file into JSON for sci-form validation (F43).

Reads the installed xtb parameter file (LGPL-2.1-or-later, (C) Grimme group;
values reproduced here for cross-validation against `xtb --gfn 0`) and
emits `tests/fixtures/validation/gfn0_params.json` plus the global block.

Usage:
    XTB=/path/to/xtb python scripts/validation/parse_gfn0_params.py
    (falls back to $XTBPATH/share/xtb/param_gfn0-xtb.txt)
"""

from __future__ import annotations

import json
import os
import shutil
import sys
import time
from pathlib import Path

KEYS = ["ao", "lev", "exp", "EN", "GAM", "KQAT2", "KCNS", "KCNP", "KCND",
        "REPA", "REPB", "POLYS", "POLYP", "POLYD", "KQS", "KQP", "KQD",
        "XI", "KAPPA", "ALPG", "CXB"]


def find_param_file() -> Path:
    env = os.environ.get("XTB")
    if env:
        cand = Path(env).resolve().parent.parent / "share" / "xtb" / "param_gfn0-xtb.txt"
        if cand.exists():
            return cand
    xp = os.environ.get("XTBPATH")
    if xp and (Path(xp) / "param_gfn0-xtb.txt").exists():
        return Path(xp) / "param_gfn0-xtb.txt"
    raise SystemExit("GFN0 param file not found (set XTB=/path/to/xtb)")


def main() -> int:
    src = find_param_file()
    print(f"reading {src}", flush=True)
    lines = src.read_text().splitlines()
    glob, elements = {}, {}
    cur = None
    for line in lines:
        s = line.strip()
        if s.startswith("$Z="):
            cur = int(s[3:].split()[0])
            elements[cur] = {}
        elif s.startswith("$globpar"):
            cur = "glob"
        elif s.startswith("$end") or s.startswith("$info") or s.startswith("$pairpar"):
            cur = None
        elif cur is not None and "=" in s:
            key, _, val = s.partition("=")
            key, val = key.strip().lower(), val.strip()
            if key == "ao":
                (elements[cur] if isinstance(cur, int) else glob)[key] = val
            else:
                try:
                    nums = [float(x) for x in val.split()]
                except ValueError:
                    continue
                (elements[cur] if isinstance(cur, int) else glob)[key] = nums
        elif cur == "glob" and s and not s.startswith("$"):
            parts = s.split()
            if len(parts) == 2:
                try:
                    glob[parts[0].lower()] = [float(parts[1])]
                except ValueError:
                    pass
    out = Path("tests/fixtures/validation/gfn0_params.json")
    out.write_text(json.dumps(
        {"source": "xtb 6.7.1 param_gfn0-xtb.txt (LGPL-2.1-or-later)",
         "generated": time.strftime("%Y-%m-%d"),
         "globals": glob, "elements": {str(k): v for k, v in elements.items()}},
        indent=1) + "\n")
    zs = sorted(elements)
    print(f"elements Z={zs[0]}..{zs[-1]} (n={len(zs)}) -> {out}")
    # sanity: show H and C records
    print("H:", json.dumps(elements[1])[:200])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
