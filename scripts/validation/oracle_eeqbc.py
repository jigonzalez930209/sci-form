#!/usr/bin/env python3
"""Generate `tests/fixtures/validation/eeqbc_charges.json` from multicharge.

Runs the multicharge CLI (`--model eeqbc`, the 2025 bond-capacitor model of
Froitzheim et al., J. Chem. Phys. 162, 214109) on a set of small molecules and
records the resulting partial charges and coordination numbers.

The multicharge CLI is a meson/ninja build of grimme-lab/multicharge (see
the EEQ-BC port notes).  It is not vendored in this repository: point
SCIFORM_REF_DATA at the checkout and the binary is looked up at
`$SCIFORM_REF_DATA/multicharge-src/build/app/multicharge`, or pass
`--binary PATH`.

Usage:
    python scripts/validation/oracle_eeqbc.py [--binary PATH] [--out FILE]
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "tests/fixtures/validation/eeqbc_charges.json"
DEFAULT_BIN = (
    Path(os.environ.get("SCIFORM_REF_DATA", REPO / ".cache" / "refdata"))
    / "multicharge-src"
    / "build"
    / "app"
    / "multicharge"
)

# Atomic numbers for the element symbols used below (mctc/multicharge range).
Z_OF_SYMBOL = {
    "H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "S": 16, "Cl": 17, "Fe": 26,
    "Cu": 29, "Zn": 30, "Sn": 50, "I": 53, "Xe": 54, "Ba": 56, "W": 74, "Au": 79,
    "Hg": 80, "Pb": 82, "Th": 90, "U": 92, "Am": 95,
}

# name -> (charge, [(element, x, y, z)]) with positions in Angstrom.
MOLECULES: dict[str, tuple[float, list[tuple[str, float, float, float]]]] = {
    "water": (
        0.0,
        [("O", 0.000000, 0.000000, 0.117300), ("H", 0.000000, 0.757200, -0.469200),
         ("H", 0.000000, -0.757200, -0.469200)],
    ),
    "ammonia": (
        0.0,
        [("N", 0.000000, 0.000000, 0.127700), ("H", 0.000000, 0.936900, -0.298000),
         ("H", 0.811300, -0.468500, -0.298000), ("H", -0.811300, -0.468500, -0.298000)],
    ),
    "methane": (
        0.0,
        [("C", 0.000000, 0.000000, 0.000000), ("H", 0.627600, 0.627600, 0.627600),
         ("H", -0.627600, -0.627600, 0.627600), ("H", -0.627600, 0.627600, -0.627600),
         ("H", 0.627600, -0.627600, -0.627600)],
    ),
    "methanol": (
        0.0,
        [("C", 0.000000, 0.000000, 0.000000), ("O", 1.400000, 0.000000, 0.000000),
         ("H", -0.380000, 0.950000, 0.000000), ("H", -0.380000, -0.475000, 0.823000),
         ("H", -0.380000, -0.475000, -0.823000), ("H", 1.800000, 0.890000, 0.000000)],
    ),
    "benzene": (
        0.0,
        [("C", 1.396800, 0.000000, 0.000000), ("C", 0.698400, 1.209700, 0.000000),
         ("C", -0.698400, 1.209700, 0.000000), ("C", -1.396800, 0.000000, 0.000000),
         ("C", -0.698400, -1.209700, 0.000000), ("C", 0.698400, -1.209700, 0.000000),
         ("H", 2.481500, 0.000000, 0.000000), ("H", 1.240800, 2.149000, 0.000000),
         ("H", -1.240800, 2.149000, 0.000000), ("H", -2.481500, 0.000000, 0.000000),
         ("H", -1.240800, -2.149000, 0.000000), ("H", 1.240800, -2.149000, 0.000000)],
    ),
    "ammonium": (
        1.0,
        [("N", 0.000000, 0.000000, 0.000000), ("H", 0.630000, 0.630000, 0.630000),
         ("H", -0.630000, -0.630000, 0.630000), ("H", -0.630000, 0.630000, -0.630000),
         ("H", 0.630000, -0.630000, -0.630000)],
    ),
    "hydroxide": (
        -1.0,
        [("O", 0.000000, 0.000000, 0.000000), ("H", 0.000000, 0.000000, 0.960000)],
    ),
    "acetic_acid": (
        0.0,
        [("C", 0.000000, 0.000000, 0.000000), ("C", 1.520000, 0.000000, 0.000000),
         ("O", 2.090000, 1.070000, 0.000000), ("O", 2.090000, -1.100000, 0.000000),
         ("H", -0.390000, 1.010000, 0.000000), ("H", -0.390000, -0.510000, 0.880000),
         ("H", -0.390000, -0.510000, -0.880000), ("H", 3.050000, -1.000000, 0.000000)],
    ),
    "thiophene": (
        0.0,
        [("S", 1.230000, 0.000000, 0.000000), ("C", 0.360000, 1.150000, 0.000000),
         ("C", -1.030000, 0.710000, 0.000000), ("C", -1.030000, -0.710000, 0.000000),
         ("C", 0.360000, -1.150000, 0.000000), ("H", 0.700000, 2.180000, 0.000000),
         ("H", -1.820000, 1.260000, 0.000000), ("H", -1.820000, -1.260000, 0.000000),
         ("H", 0.700000, -2.180000, 0.000000)],
    ),
    "zinc_oxide": (
        0.0,
        [("Zn", 0.000000, 0.000000, 0.000000), ("O", 0.000000, 0.000000, 1.720000)],
    ),
    "copper_chloride": (
        0.0,
        [("Cu", 0.000000, 0.000000, 0.000000), ("Cl", 0.000000, 0.000000, 2.050000),
         ("Cl", 2.050000, 0.000000, 0.000000)],
    ),
    "tin_tetrafluoride": (
        0.0,
        [("Sn", 0.000000, 0.000000, 0.000000), ("F", 1.880000, 0.000000, 0.000000),
         ("F", -1.880000, 0.000000, 0.000000), ("F", 0.000000, 1.880000, 0.000000),
         ("F", 0.000000, -1.880000, 0.000000)],
    ),
    "lead_iodide": (
        0.0,
        [("Pb", 0.000000, 0.000000, 0.000000), ("I", 0.000000, 0.000000, 2.790000),
         ("I", 2.490000, 0.000000, -1.130000)],
    ),
    "mercury_dichloride": (
        0.0,
        [("Hg", 0.000000, 0.000000, 0.000000), ("Cl", 0.000000, 0.000000, 2.250000),
         ("Cl", 0.000000, 0.000000, -2.250000)],
    ),
    "gold_dimer": (
        0.0,
        [("Au", 0.000000, 0.000000, 0.000000), ("Au", 0.000000, 0.000000, 2.470000)],
    ),
    "xenon_difluoride": (
        0.0,
        [("Xe", 0.000000, 0.000000, 0.000000), ("F", 0.000000, 0.000000, 1.980000),
         ("F", 0.000000, 0.000000, -1.980000)],
    ),
    "barium_oxide": (
        0.0,
        [("Ba", 0.000000, 0.000000, 0.000000), ("O", 0.000000, 0.000000, 1.940000)],
    ),
    "thorium_dioxide": (
        0.0,
        [("Th", 0.000000, 0.000000, 0.000000), ("O", 0.000000, 0.000000, 2.090000),
         ("O", 0.000000, 0.000000, -2.090000)],
    ),
    "uranium_hexafluoride": (
        0.0,
        [("U", 0.000000, 0.000000, 0.000000), ("F", 0.000000, 0.000000, 1.996000),
         ("F", 0.000000, 0.000000, -1.996000), ("F", 1.996000, 0.000000, 0.000000),
         ("F", -1.996000, 0.000000, 0.000000), ("F", 0.000000, 1.996000, 0.000000),
         ("F", 0.000000, -1.996000, 0.000000)],
    ),
    "tungsten_trioxide": (
        0.0,
        [("W", 0.000000, 0.000000, 0.000000), ("O", 1.700000, 0.000000, 0.000000),
         ("O", -0.850000, 1.472000, 0.000000), ("O", -0.850000, -1.472000, 0.000000)],
    ),
    "americium_trichloride": (
        0.0,
        [("Am", 0.000000, 0.000000, 0.000000), ("Cl", 0.000000, 0.000000, 2.600000),
         ("Cl", 2.250000, 0.000000, -1.300000), ("Cl", -2.250000, 0.000000, -1.300000)],
    ),
    "aspirin": (
        0.0,
        [("C", -1.210000, 0.390000, 0.000000), ("C", 0.000000, 1.140000, 0.000000),
         ("C", 1.210000, 0.390000, 0.000000), ("C", 1.210000, -1.040000, 0.000000),
         ("C", 0.000000, -1.780000, 0.000000), ("C", -1.210000, -1.040000, 0.000000),
         ("C", 2.550000, 1.170000, 0.000000), ("O", 3.600000, 0.580000, 0.000000),
         ("O", 2.430000, 2.370000, 0.000000), ("O", -2.550000, 1.170000, 0.000000),
         ("C", -3.780000, 0.470000, 0.000000), ("O", -2.480000, 2.370000, 0.000000),
         ("H", 0.000000, 2.220000, 0.000000), ("H", 2.150000, -1.590000, 0.000000),
         ("H", 0.000000, -2.860000, 0.000000), ("H", -2.150000, -1.590000, 0.000000),
         ("H", -4.600000, 1.180000, 0.000000), ("H", -3.950000, -0.140000, 0.000000),
         ("H", -3.680000, -0.180000, 0.900000), ("H", -2.840000, 3.050000, 0.000000),
         ("H", -1.610000, 2.700000, 0.000000)],
    ),
    "ferrocene_fragment": (
        0.0,
        [("Fe", 0.000000, 0.000000, 0.000000), ("C", 1.200000, 0.000000, 1.650000),
         ("C", 0.370000, 1.141000, 1.650000), ("C", -0.970000, 0.706000, 1.650000),
         ("C", -0.970000, -0.706000, 1.650000), ("C", 0.370000, -1.141000, 1.650000),
         ("C", 1.200000, 0.000000, -1.650000), ("C", 0.370000, 1.141000, -1.650000),
         ("C", -0.970000, 0.706000, -1.650000), ("C", -0.970000, -0.706000, -1.650000),
         ("C", 0.370000, -1.141000, -1.650000)],
    ),
}


def run_oracle(binary: Path, name: str, charge: float, atoms: list) -> dict:
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        xyz = tmp / f"{name}.xyz"
        lines = [str(len(atoms)), name]
        for el, x, y, z in atoms:
            lines.append(f"{el} {x:.6f} {y:.6f} {z:.6f}")
        xyz.write_text("\n".join(lines) + "\n")
        cmd = [str(binary), "--model", "eeqbc", "-i", "xyz", "-j", str(xyz)]
        if charge != 0.0:
            cmd = [str(binary), "--model", "eeqbc", "-c", str(charge), "-i", "xyz", "-j", str(xyz)]
        subprocess.run(cmd, cwd=tmp, check=True, capture_output=True, text=True)
        data = json.loads((tmp / "multicharge.json").read_text())
    return {
        "name": name,
        "charge": charge,
        "elements": [Z_OF_SYMBOL[el] for el, *_ in atoms],
        "positions": [[x, y, z] for _, x, y, z in atoms],
        "charges": data["charges"],
        "cn": data["coordination numbers"],
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", type=Path, default=DEFAULT_BIN)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    version = subprocess.run(
        [str(args.binary), "--version"], capture_output=True, text=True, check=True
    ).stdout.strip()

    molecules = []
    for name, (charge, atoms) in MOLECULES.items():
        entry = run_oracle(args.binary, name, charge, atoms)
        net = sum(entry["charges"])
        assert abs(net - charge) < 1e-10, f"{name}: net charge {net} != {charge}"
        molecules.append(entry)
        print(f"  {name:24s} n={len(atoms):3d} Q={charge:+.0f} Σq={net:+.2e}")

    fixture = {
        "source": f"multicharge CLI --model eeqbc ({version})",
        "model": "EEQ-BC 2025 (Froitzheim et al., J. Chem. Phys. 162, 214109)",
        "units": {"positions": "angstrom", "charges": "elementary charge"},
        "molecules": molecules,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(fixture, indent=2) + "\n")
    print(f"wrote {args.out.relative_to(REPO)} ({len(molecules)} molecules)")


if __name__ == "__main__":
    main()
