#!/usr/bin/env python3
"""xtb Hessian oracle for sci-form vibrational validation (ficha F48).

For each molecule: optimize with xtb (`--opt`), then compute the analytical
Hessian (`--hess`) at the optimized geometry. Stores the Hessian matrix
(Hartree/Bohr^2, atom-major XYZ) and the real vibrational frequencies.

This validates sci-form's mass-weighting/projection/frequency-conversion
machinery (`ir::build_vibrational_analysis_from_hessian`) against xtb's,
independent of the energy method: same matrix in, same frequencies out.

Generates `tests/fixtures/validation/xtb_hessian.json`.

Usage:
    XTB=/path/to/xtb python scripts/validation/oracle_hessian.py
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

from common import Report, emit_fixture_report
from oracle_tblite import SMILES, embed

Z2SYM = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 15: "P", 16: "S", 17: "Cl"}

# Small molecules: full Hessian is cheap, frequencies well-defined.
SMILES_HESS = [
    ("water", "O"),
    ("methane", "C"),
    ("ammonia", "N"),
    ("formaldehyde", "C=O"),
    ("ethylene", "C=C"),
    ("hydrogen sulfide", "S"),
]


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
    env = dict(os.environ)
    if "XTBPATH" not in env:
        share = Path(binary).resolve().parent.parent / "share" / "xtb"
        if share.is_dir():
            env["XTBPATH"] = str(share)
    return env


def parse_xyz(text: str):
    lines = text.strip().splitlines()
    elements, coords = [], []
    for line in lines[2:]:
        parts = line.split()
        if len(parts) < 4:
            continue
        sym = parts[0].capitalize()
        z = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "S": 16, "P": 15, "Cl": 17}[sym]
        elements.append(z)
        coords += [float(parts[1]), float(parts[2]), float(parts[3])]
    return elements, coords


def parse_hessian(text: str, n: int) -> list:
    vals = []
    for line in text.splitlines():
        if line.strip().startswith("$"):
            continue
        vals += [float(x) for x in line.split()]
    assert len(vals) == n * n, f"hessian size {len(vals)} != {n*n}"
    return vals


def parse_vibspectrum(text: str) -> list:
    freqs = []
    for line in text.splitlines():
        m = re.match(r"\s*\d+\s+(?:\S+\s+)?(-?\d+\.\d+)\s+(\d+\.\d+)", line)
        if m:
            freqs.append(float(m.group(1)))
    return freqs


def run_case(binary: str, env: dict, name: str, smiles: str) -> dict:
    elements, coords = embed(smiles)
    lines = [str(len(elements)), "sci-form validation corpus"]
    for z, (x, y, zz) in zip(elements, [coords[i : i + 3] for i in range(0, len(coords), 3)]):
        lines.append(f"{Z2SYM[z]:2s} {x: .10f} {y: .10f} {zz: .10f}")
    with tempfile.TemporaryDirectory() as tmp:
        xyz = Path(tmp) / "mol.xyz"
        xyz.write_text("\n".join(lines) + "\n")
        opt = subprocess.run(
            [binary, "--gfn", "2", "--opt", "tight", str(xyz)],
            cwd=tmp, capture_output=True, text=True, timeout=1200, env=env,
        )
        if opt.returncode != 0:
            raise RuntimeError(f"xtb opt failed for {name}:\n{(opt.stdout + opt.stderr)[-500:]}")
        opt_xyz = Path(tmp) / "xtbopt.xyz"
        if not opt_xyz.exists():
            raise RuntimeError(f"no xtbopt.xyz for {name}")
        elements, coords = parse_xyz(opt_xyz.read_text())
        (Path(tmp) / "xtbopt.xyz").rename(Path(tmp) / "mol.xyz")
        hess = subprocess.run(
            [binary, "--gfn", "2", "--hess", str(Path(tmp) / "mol.xyz")],
            cwd=tmp, capture_output=True, text=True, timeout=1200, env=env,
        )
        if hess.returncode != 0:
            raise RuntimeError(f"xtb hess failed for {name}")
        hess_txt = (Path(tmp) / "hessian").read_text()
        vib_txt = (Path(tmp) / "vibspectrum").read_text()
        n = len(elements)
        return {
            "name": name,
            "smiles": smiles,
            "elements": elements,
            "coords": coords,
            "n_atoms": n,
            "hessian_hartree_bohr2": parse_hessian(hess_txt, 3 * n),
            "frequencies_cm1": parse_vibspectrum(vib_txt),
        }


def main() -> int:
    binary = xtb_path()
    env = xtb_env(binary)
    cases = []
    for name, smiles in SMILES_HESS:
        cases.append(run_case(binary, env, name, smiles))
        print(f"  {name}: {len(cases[-1]['frequencies_cm1'])} modes", flush=True)

    payload = {
        "reference": "xtb 6.7.1 GFN2-xTB (opt tight + hess)",
        "hessian_units": "Hartree/Bohr^2, atom-major XYZ",
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F48",
        observable="hessian_diagonalization_vs_xtb",
        reference_tool="xtb GFN2-xTB vibspectrum",
        units="cm^-1",
        tolerance=2.0,
        level="L2",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "Same Hessian matrix in, same frequencies out: validates "
            "mass-weighting, rigid-body projection and cm^-1 conversion. "
            "Tolerance 2 cm^-1 per the validation plan. Measured by "
            "tests/cross_validation.rs::f48_hessian."
        ),
    )
    emit_fixture_report("xtb_hessian.json", payload, report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
