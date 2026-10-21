#!/usr/bin/env python3
"""PM3 heat-of-formation oracle for sci-form validation (ficha F42).

Two fixtures are produced:

  - `pm3_dhf.json`: deterministic geometries (RDKit ETKDG + MMFF94) plus
    experimental gas-phase heats of formation (NIST WebBook, 298 K) used
    by the regression guard `f42_pm3_guard` in tests/cross_validation.rs;
  - `mopac_pm3.json` (mode `mopac-fixture`): the PM3 single-point reference
    values that `tests/cross_validation.rs::f42_mopac_pm3` compares against.
    The reference geometry is **sci-form's own seed-42 embedding**, so the
    fixture has to be regenerated whenever the embedding pipeline changes
    (the heat of formation is geometry dependent).

The MOPAC branch is optional: if a `mopac` binary is on PATH (or `MOPAC`
points at it), reference PM3 single points (`PM3 1SCF`) are computed.

Usage:
    python scripts/validation/oracle_pm3.py
    python scripts/validation/oracle_pm3.py mopac-fixture
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

from common import Report, cli_path, emit_fixture_report

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError as exc:  # pragma: no cover
    raise SystemExit(f"RDKit is required for this oracle: {exc}")

REFERENCE = "NIST WebBook gas-phase ΔHf°(298 K)"

#: The molecules of the `mopac_pm3.json` fixture (name, SMILES).  The
#: reference geometry is sci-form's own seed-42 embedding.
MOPAC_FIXTURE_MOLECULES = [
    ("water", "O"),
    ("methane", "C"),
    ("ethane", "CC"),
    ("ethylene", "C=C"),
    ("acetylene", "C#C"),
    ("benzene", "c1ccccc1"),
    ("methanol", "CO"),
    ("acetic_acid", "CC(=O)O"),
    ("formaldehyde", "C=O"),
    ("acetone", "CC(=O)C"),
    ("ammonia", "N"),
    ("methylamine", "CN"),
    ("pyridine", "c1ccncc1"),
    ("pyrrole", "c1cc[nH]c1"),
    ("thiophene", "c1ccsc1"),
    ("chloromethane", "CCl"),
    ("dichloromethane", "C(Cl)Cl"),
    ("silane", "[SiH4]"),
    ("phosphine", "P"),
    ("hydrogen_sulfide", "S"),
    ("fluoroform", "FC(F)F"),
    ("phenol", "Oc1ccccc1"),
    ("aniline", "Nc1ccccc1"),
    ("toluene", "Cc1ccccc1"),
]

#: Seed of the sci-form embedding used for the MOPAC reference geometries.
MOPAC_FIXTURE_SEED = 42

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

Z2SYM = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 14: "Si", 15: "P", 16: "S", 17: "Cl"}


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
                    # MOPAC prints "HEAT OF FORMATION = <kcal/mol> KCAL/MOL = <kJ/mol> KJ/MOL";
                    # the kcal/mol value is the one after the *first* "=".
                    return float(line.split("=")[1].split()[0])
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


def sciform_embed(smiles: str, seed: int):
    """sci-form's own seed-`seed` embedding (via the sci-form CLI)."""
    import json
    import subprocess

    cmd = [str(cli_path()), "embed", smiles, "--seed", str(seed), "--format", "json"]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if proc.returncode != 0:
        return None
    data = json.loads(proc.stdout)
    return data["coords"], data.get("elements", [])


def run_mopac_fixture() -> Report:
    """Regenerate `mopac_pm3.json` (F42) from sci-form's seed-42 embeddings.

    The reference geometry is sci-form's own embedding, so the fixture has to
    be regenerated whenever the embedding pipeline changes: the PM3 heat of
    formation is geometry dependent.
    """
    binary = mopac_path()
    if binary is None:
        raise SystemExit(
            "a MOPAC binary is required to regenerate mopac_pm3.json "
            "(install it with scripts/validation/install_mopac.sh)"
        )
    cases = []
    max_delta = 0.0
    for name, smi in MOPAC_FIXTURE_MOLECULES:
        embedded = sciform_embed(smi, MOPAC_FIXTURE_SEED)
        if embedded is None:
            raise SystemExit(f"sci-form embedding failed for {name} ({smi})")
        coords, _ = embedded
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        elements = [a.GetAtomicNum() for a in mol.GetAtoms()]
        ref = run_mopac(binary, elements, coords)
        if ref is None:
            raise SystemExit(f"MOPAC failed for {name} ({smi})")
        cases.append(
            {
                "name": name,
                "smiles": smi,
                "mopac_hof_kcal": ref,
            }
        )
    payload = {
        "reference": "MOPAC 23.2.5 PM3 1SCF on the sci-form seed-42 embedding",
        "source": "MOPAC 23.2.5 (PM3 1SCF) at the sci-form embed geometries (seed 42)",
        "notes": [
            "F42 arbiter for the faithful PM3 port (src/pm3/faithful.rs).",
            "The reference geometry is sci-form's own embedding, so this fixture "
            "must be regenerated (python scripts/validation/oracle_pm3.py "
            "mopac-fixture) whenever the embedding pipeline changes.",
        ],
        "seed": MOPAC_FIXTURE_SEED,
        "tolerance": 0.02,
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F42",
        observable="mopac_pm3_max_abs",
        reference_tool="MOPAC 23.2.5 PM3 1SCF",
        units="kcal/mol",
        tolerance=0.02,
        level="L2",
        n_cases=len(cases),
        metric="max_abs",
        # The deviation is measured by the Rust gate
        # (tests/cross_validation.rs::f42_mopac_pm3), which records the real
        # value into target/validation/F42.json; `collect.py` merges it over
        # this placeholder. The oracle cannot evaluate sci-form's PM3 itself,
        # so it must not report a fabricated value.
        value=float("nan"),
        notes=(
            "reference values are recomputed on sci-form's own seed-42 "
            "embedding, so this fixture must be regenerated whenever the "
            "embedding pipeline changes; the measured deviation is recorded "
            "by the Rust gate f42_mopac_pm3"
        ),
    )
    emit_fixture_report("mopac_pm3.json", payload, report)
    return report


def main() -> int:
    if len(sys.argv) > 1 and sys.argv[1] == "mopac-fixture":
        run_mopac_fixture()
        print("fixtures -> tests/fixtures/validation")
        return 0
    run()
    print("fixtures -> tests/fixtures/validation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
