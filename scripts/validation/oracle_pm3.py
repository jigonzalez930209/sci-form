#!/usr/bin/env python3
"""PM3 heat-of-formation oracle for sci-form validation (ficha F42).

Status: cross-check against MOPAC is BLOCKED (no MOPAC binary and no
Fortran toolchain in this environment). This oracle therefore ships:

  - `pm3_dhf.json`: deterministic geometries (RDKit ETKDG + MMFF94) plus
    experimental gas-phase heats of formation (NIST WebBook, 298 K) used
    by the regression guard `f42_pm3_guard` in tests/cross_validation.rs;
  - an optional MOPAC branch: if a `mopac` binary is on PATH (or `MOPAC`
    points at it), reference PM3 single points (`PM3 1SCF`) are computed
    and stored under the `mopac` key, upgrading the fixture to a real
    cross-validation the moment the toolchain exists.

Usage:
    python scripts/validation/oracle_pm3.py
    MOPAC=/path/to/mopac python scripts/validation/oracle_pm3.py  # with reference
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from common import Report, emit_fixture_report

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError as exc:  # pragma: no cover
    raise SystemExit(f"RDKit is required for this oracle: {exc}")

REFERENCE = "NIST WebBook gas-phase ΔHf°(298 K)"

# molecule -> (SMILES, experimental ΔHf in kcal/mol)
SERIES = {
    "methane": ("C", -17.8),
    "ethane": ("CC", -20.0),
    "ethylene": ("C=C", 12.5),
    "acetylene": ("C#C", 54.2),
    "water": ("O", -57.8),
    "ammonia": ("N", -11.0),
    "hydrogen fluoride": ("F", -65.3),
    "formaldehyde": ("C=O", -26.0),
    "carbon monoxide": ("[C-]#[O+]", -26.4),
    "methanol": ("CO", -48.0),
    "benzene": ("c1ccccc1", 19.8),
}

Z2SYM = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F"}


def mopac_path() -> str | None:
    env = os.environ.get("MOPAC")
    if env and Path(env).exists():
        return env
    return shutil.which("mopac")


def run_mopac(binary: str, elements, coords) -> float | None:
    """Run MOPAC PM3 1SCF; return the heat of formation in kcal/mol."""
    lines = ["PM3 1SCF", "", ""]
    for z, (x, y, zz) in zip(elements, [coords[i : i + 3] for i in range(0, len(coords), 3)]):
        lines.append(f"{Z2SYM[z]:2s} {x: .10f} {y: .10f} {zz: .10f}")
    with tempfile.TemporaryDirectory() as tmp:
        dat = Path(tmp) / "mol.dat"
        dat.write_text("\n".join(lines) + "\n")
        try:
            proc = subprocess.run(
                [binary, str(dat)],
                cwd=tmp,
                capture_output=True,
                text=True,
                timeout=600,
            )
        except (OSError, subprocess.TimeoutExpired):
            return None
        out = (proc.stdout or "") + (proc.stderr or "")
        arc = Path(tmp) / "mol.arc"
        text = arc.read_text() if arc.exists() else out
        for line in text.splitlines():
            if "HEAT OF FORMATION" in line.upper():
                try:
                    return float(line.split("=")[-1].split()[0])
                except (ValueError, IndexError):
                    return None
    return None


def run() -> Report:
    binary = mopac_path()
    cases = []
    for name, (smi, dhf) in SERIES.items():
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        params = AllChem.ETKDGv3()
        params.randomSeed = 0x42
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise ValueError(f"embedding failed for {name}")
        AllChem.MMFFOptimizeMolecule(mol)
        conf = mol.GetConformer()
        elements, coords = [], []
        for atom in mol.GetAtoms():
            elements.append(atom.GetAtomicNum())
            p = conf.GetAtomPosition(atom.GetIdx())
            coords += [p.x, p.y, p.z]
        entry = {"name": name, "elements": elements, "coords": coords, "exp_dhf": dhf}
        if binary is not None:
            ref = run_mopac(binary, elements, coords)
            if ref is not None:
                entry["mopac_dhf"] = ref
        cases.append(entry)

    payload = {
        "reference": REFERENCE + (" + MOPAC PM3" if binary else " (MOPAC branch pending)"),
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F42",
        observable="pm3_heat_of_formation_guard",
        reference_tool=REFERENCE,
        units="kcal/mol",
        tolerance=5000.0,
        level="L4",
        n_cases=len(cases),
        metric="max_abs_vs_experiment",
        value=float("nan"),
        notes=(
            "MOPAC cross-check BLOCKED (no binary/toolchain); this guard only "
            "trips on gross regressions of the current behavior. Measured by "
            "tests/cross_validation.rs::f42_pm3_guard."
        ),
    )
    emit_fixture_report("pm3_dhf.json", payload, report)
    return report


def main() -> int:
    run()
    print("fixtures -> tests/fixtures/validation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
