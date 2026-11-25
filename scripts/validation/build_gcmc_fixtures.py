#!/usr/bin/env python3
"""Build GCMC framework/guest fixtures from a RASPA-style CIF.

The CIF must carry `_atom_site_charge` and a `_symmetry_equiv_pos_as_xyz`
loop. The script expands the asymmetric unit, converts fractional to Cartesian
coordinates and writes the JSON types consumed by `sci-form gcmc`.

Examples:
    python3 scripts/validation/build_gcmc_fixtures.py \
        --cif .../Cu-BTC.cif --output tests/fixtures/validation/gcmc_cu_btc.json
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

# DREIDING/UFF parameters used by RASPA's ExampleMOFsForceField, in K and Å.
LJ_K = {
    "H": (7.64893, 2.84642),
    "C": (47.8562, 3.47299),
    "N": (38.9492, 3.26256),
    "O": (48.1581, 3.03315),
    "F": (36.4834, 3.0932),
    "Cu": (2.5161, 3.11369),
    "Zn": (62.3992, 2.46155),
}
K_TO_KCAL = 0.001_987_204_258_640_83

# TraPPE/Garcia-Sanchez parameters from ExampleMOFsForceField.
GUESTS = {
    "co2": {
        "sites": [
            {"offset": [0.0, 0.0, 1.149], "epsilon_k": 85.671, "sigma": 3.017, "charge": -0.3256},
            {"offset": [0.0, 0.0, 0.0], "epsilon_k": 29.933, "sigma": 2.745, "charge": 0.6512},
            {"offset": [0.0, 0.0, -1.149], "epsilon_k": 85.671, "sigma": 3.017, "charge": -0.3256},
        ],
        "tc": 304.1282,
        "pc_bar": 73.773,
        "omega": 0.22394,
    },
    "n2": {
        "sites": [
            {"offset": [0.0, 0.0, 0.55], "epsilon_k": 38.298, "sigma": 3.306, "charge": -0.405},
            {"offset": [0.0, 0.0, 0.0], "epsilon_k": 0.0, "sigma": 1.0, "charge": 0.810},
            {"offset": [0.0, 0.0, -0.55], "epsilon_k": 38.298, "sigma": 3.306, "charge": -0.405},
        ],
        "tc": 126.192,
        "pc_bar": 33.958,
        "omega": 0.0372,
    },
}


def raspa_type_maps(
    pseudo_atoms: Path | None, mixing_rules: Path | None
) -> tuple[dict[str, float], dict[str, tuple[float, float]]]:
    """Read per-type charges and LJ parameters from RASPA force-field files."""
    charges: dict[str, float] = {}
    if pseudo_atoms is not None:
        for line in pseudo_atoms.read_text().splitlines():
            fields = line.split("#", 1)[0].split()
            if len(fields) >= 6 and fields[0] not in {"#type", "type"}:
                try:
                    charges[fields[0]] = float(fields[6])
                except ValueError:
                    continue
    lj: dict[str, tuple[float, float]] = {}
    if mixing_rules is not None:
        for line in mixing_rules.read_text().splitlines():
            fields = line.split("#", 1)[0].split()
            if len(fields) >= 5 and fields[1] == "lennard-jones":
                try:
                    lj[fields[0]] = (float(fields[2]), float(fields[3]))
                except ValueError:
                    continue
    return charges, lj


def number(text: str, key: str) -> float:
    match = re.search(rf"^{re.escape(key)}\s+([^\s]+)", text, re.M)
    if not match:
        raise SystemExit(f"missing {key}")
    return float(match.group(1).replace("(", "").replace(")", ""))


def symmetry_ops(text: str) -> list[str]:
    block = re.search(
        r"_symmetry_equiv_pos_as_xyz\s*\n(.*?)(?:\n\s*\n|\n_|\Z)", text, re.S
    )
    if not block:
        return ["x,y,z"]
    ops = []
    for line in block.group(1).splitlines():
        line = line.strip().strip("'\"")
        if not line or line.startswith("#"):
            continue
        if line.startswith("_") or line.startswith("loop_") or line.startswith("data_"):
            break
        ops.append(line)
    return ops


def atom_rows(text: str) -> list[dict]:
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.strip() != "_atom_site_label":
            continue
        tags = []
        j = i
        while j < len(lines) and lines[j].strip().startswith("_"):
            tags.append(lines[j].strip())
            j += 1
        rows = []
        while j < len(lines):
            row = lines[j].strip()
            if not row or row.startswith("#") or row.startswith("_") or row.startswith("loop_"):
                break
            parts = row.split()
            if len(parts) < len(tags):
                break
            record = dict(zip(tags, parts))
            rows.append(record)
            j += 1
        return rows
    raise SystemExit("no _atom_site loop")


def evaluate(op: str, x: float, y: float, z: float) -> tuple[float, float, float]:
    # The CIF expressions only use x, y, z, integers, + - * / and parentheses.
    safe = {"__builtins__": {}, "x": x, "y": y, "z": z}
    values = [eval(part.replace("^", "**"), safe) for part in op.split(",")]
    return float(values[0]), float(values[1]), float(values[2])


def lattice(text: str) -> list[list[float]]:
    a = number(text, "_cell_length_a")
    b = number(text, "_cell_length_b")
    c = number(text, "_cell_length_c")
    alpha = math.radians(number(text, "_cell_angle_alpha"))
    beta = math.radians(number(text, "_cell_angle_beta"))
    gamma = math.radians(number(text, "_cell_angle_gamma"))
    sin_gamma = math.sin(gamma)
    c_x = c * math.cos(beta)
    c_y = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / sin_gamma
    c_z = math.sqrt(max(c * c - c_x * c_x - c_y * c_y, 0.0))
    return [
        [a, 0.0, 0.0],
        [b * math.cos(gamma), b * sin_gamma, 0.0],
        [c_x, c_y, c_z],
    ]


def cartesian(lattice_vectors: list[list[float]], frac: tuple[float, float, float]) -> list[float]:
    return [
        sum(frac[k] * lattice_vectors[k][axis] for k in range(3))
        for axis in range(3)
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cif", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--guest", choices=sorted(GUESTS), default="co2")
    parser.add_argument("--temperature-k", type=float, default=323.0)
    parser.add_argument("--pressure-bar", type=float, default=10.0)
    parser.add_argument("--cutoff", type=float, default=12.0)
    parser.add_argument("--pseudo-atoms", type=Path)
    parser.add_argument("--mixing-rules", type=Path)
    args = parser.parse_args()
    type_charges, type_lj = raspa_type_maps(args.pseudo_atoms, args.mixing_rules)

    text = args.cif.read_text()
    cell = lattice(text)
    ops = symmetry_ops(text)
    rows = atom_rows(text)

    seen: set[tuple[str, int, int, int]] = set()
    atoms = []
    for row in rows:
        element = row.get("_atom_site_type_symbol", "").strip()
        label = row.get("_atom_site_label", "").strip()
        use_type = label if label in type_lj else element
        if use_type in type_lj:
            epsilon_k, sigma = type_lj[use_type]
        elif element in LJ_K:
            epsilon_k, sigma = LJ_K[element]
        else:
            raise SystemExit(f"no LJ parameters for {element} / {label}")
        if label in type_charges:
            charge = type_charges[label]
        else:
            charge = float(
                row.get("_atom_site_charge", "0").replace("(", "").replace(")", "")
            )
        frac = (
            float(row["_atom_site_fract_x"]),
            float(row["_atom_site_fract_y"]),
            float(row["_atom_site_fract_z"]),
        )
        for op in ops:
            x, y, z = evaluate(op, *frac)
            wrapped = (x % 1.0, y % 1.0, z % 1.0)
            key = (use_type, *(round(value * 1e5) for value in wrapped))
            if key in seen:
                continue
            seen.add(key)
            atoms.append(
                {
                    "position": cartesian(cell, wrapped),
                    "epsilon_kcal_mol": epsilon_k * K_TO_KCAL,
                    "sigma_angstrom": sigma,
                    "charge": charge,
                }
            )

    # Snap floating-point noise from the trigonometric cell construction so
    # the hot GCMC loop can use the exact-orthogonal minimum-image path.
    cell = [
        [0.0 if abs(value) < 1e-12 else value for value in row]
        for row in cell
    ]
    framework = {
        "cell": {
            "lattice": cell,
        },
        "atoms": atoms,
    }
    guest_spec = GUESTS[args.guest]
    guest = {
        "sites": [
            {
                "offset_angstrom": site["offset"],
                "epsilon_kcal_mol": site["epsilon_k"] * K_TO_KCAL,
                "sigma_angstrom": site["sigma"],
                "charge": site["charge"],
            }
            for site in guest_spec["sites"]
        ],
        "critical_temperature_k": guest_spec["tc"],
        "critical_pressure_bar": guest_spec["pc_bar"],
        "acentric_factor": guest_spec["omega"],
    }
    config = {
        "temperature_k": args.temperature_k,
        "pressures_bar": [args.pressure_bar],
        "equilibration_steps": 100_000,
        "production_steps": 500_000,
        "sample_interval": 100,
        "seed": 42,
        "max_guests": 200,
        "translation_probability": 0.3,
        "rotation_probability": 0.3,
        "cutoff_angstrom": args.cutoff,
        "electrostatics": {"kind": "cutoff"},
        "tail_correction": False,
    }
    payload = {"framework": framework, "guest": guest, "config": config}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"{args.output}: {len(atoms)} framework atoms (sum charge "
        f"{sum(atom['charge'] for atom in atoms):+.4f}), {args.guest}, "
        f"{args.temperature_k} K, {args.pressure_bar} bar"
    )


if __name__ == "__main__":
    main()
