#!/usr/bin/env python3
"""xtb oracle for sci-form GFN0-xTB cross-validation (ficha F43).

Requires the official `xtb` binary (grimme-lab release, LGPL-3). Point the
`XTB` environment variable at it, or have `xtb` on PATH:

    curl -L -o xtb.tar.xz \\
      https://github.com/grimme-lab/xtb/releases/download/v6.7.1/xtb-6.7.1-linux-x86_64.tar.xz
    tar -xJf xtb.tar.xz
    export XTB=$PWD/xtb-dist/bin/xtb

Uses the deterministic corpus from `oracle_tblite.py` and evaluates
GFN0-xTB single points.

Generates `tests/fixtures/validation/xtb_gfn0.json`.

Usage:
    XTB=/path/to/xtb python scripts/validation/oracle_xtb.py
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from common import EV2KCAL, HARTREE2EV, Report, emit_fixture_report
from oracle_tblite import SMILES, embed

Z2SYM = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 15: "P", 16: "S", 17: "Cl"}


def xtb_path() -> str:
    env = os.environ.get("XTB")
    if env and Path(env).exists():
        return env
    found = shutil.which("xtb")
    if found:
        return found
    print("xtb binary not found; set XTB=/path/to/xtb", file=sys.stderr)
    sys.exit(1)


def xtb_env(binary: str) -> dict:
    """Set XTBPATH so the binary finds its parameter files."""
    env = dict(os.environ)
    if "XTBPATH" not in env:
        prefix = Path(binary).resolve().parent.parent
        share = prefix / "share" / "xtb"
        if share.is_dir():
            env["XTBPATH"] = str(share)
    return env


def run_xtb(binary: str, elements, coords) -> dict:
    lines = [str(len(elements)), "sci-form validation corpus"]
    for z, (x, y, zz) in zip(elements, [coords[i : i + 3] for i in range(0, len(coords), 3)]):
        lines.append(f"{Z2SYM[z]:2s} {x: .10f} {y: .10f} {zz: .10f}")
    with tempfile.TemporaryDirectory() as tmp:
        xyz = Path(tmp) / "mol.xyz"
        xyz.write_text("\n".join(lines) + "\n")
        proc = subprocess.run(
            [binary, "--gfn", "0", "--sp", str(xyz)],
            cwd=tmp,
            capture_output=True,
            text=True,
            timeout=600,
            env=xtb_env(binary),
        )
        out = proc.stdout + proc.stderr
        if proc.returncode != 0:
            raise RuntimeError(f"xtb failed ({proc.returncode}): {out[-400:]}")
        m = re.search(r"TOTAL ENERGY\s+(-?\d+\.\d+)\s*Eh", out)
        if not m:
            raise RuntimeError(f"could not parse xtb energy:\n{out[-400:]}")
        e_hartree = float(m.group(1))
        homo = lumo = None
        mh = re.search(r"(-?\d+\.\d+)\s*\(HOMO\)", out)
        ml = re.search(r"(-?\d+\.\d+)\s*\(LUMO\)", out)
        if mh:
            homo = float(mh.group(1))
        if ml:
            lumo = float(ml.group(1))
        # GFN0 energy components (Hartree) for per-term validation (F43).
        comp = {}
        for key, pat in [
            ("h0", r"H0 energy\s+(-?\d+\.\d+)\s*Eh"),
            ("repulsion", r"repulsion energy\s+(-?\d+\.\d+)\s*Eh"),
            ("electrostatic", r"electrostat energy\s+(-?\d+\.\d+)\s*Eh"),
            ("dispersion", r"dispersion energy\s+(-?\d+\.\d+)\s*Eh"),
        ]:
            m = re.search(pat, out)
            if m:
                comp[key] = float(m.group(1))
        return {
            "e_tot_ev": e_hartree * HARTREE2EV,
            "e_tot_hartree": e_hartree,
            "homo_ev": homo,
            "lumo_ev": lumo,
            "components_hartree": comp,
        }


def main() -> int:
    binary = xtb_path()
    cases = []
    for name, smiles in SMILES:
        elements, coords = embed(smiles)
        entry = {"name": name, "elements": elements, "coords": coords}
        entry["gfn0"] = run_xtb(binary, elements, coords)
        cases.append(entry)

    payload = {
        "reference": "xtb 6.7.1 (--gfn 0 --sp)",
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F43",
        observable="gfn0_single_point_ev",
        reference_tool="xtb GFN0-xTB",
        units="eV",
        tolerance=EV2KCAL,
        level="L4",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "Tolerance 1 kcal/mol (0.0434 eV) per the validation plan. "
            "Measured by tests/cross_validation.rs::f43_gfn0."
        ),
    )
    emit_fixture_report("xtb_gfn0.json", payload, report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
