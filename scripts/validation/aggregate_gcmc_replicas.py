#!/usr/bin/env python3
"""Aggregate independent GCMC replicas into one loading estimate.

Each replica is a ``sci-form gcmc`` JSON result with a single pressure point.
The combined mean uses inverse-variance weights on the reported standard
errors.  A between-replica term is added so the reported uncertainty does not
hide run-to-run scatter.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
from pathlib import Path


def load_replicas(pattern: str) -> list[dict]:
    rows = []
    for path in sorted(glob.glob(pattern)):
        try:
            data = json.loads(Path(path).read_text())
            point = data["points"][0]
        except (OSError, ValueError, KeyError, IndexError):
            continue
        rows.append(
            {
                "file": path,
                "loading": float(point["mean_loading"]),
                "standard_error": float(point["standard_error"]),
                "samples": int(point["n_samples"]),
                "energy_per_molecule": float(
                    point.get("mean_energy_per_molecule_kcal_mol", 0.0)
                ),
            }
        )
    return rows


def combine(rows: list[dict]) -> dict:
    if not rows:
        raise SystemExit("no valid replica files found")
    weights = [1.0 / max(row["standard_error"], 1e-9) ** 2 for row in rows]
    total_weight = sum(weights)
    mean = sum(weight * row["loading"] for weight, row in zip(weights, rows)) / total_weight
    within = math.sqrt(1.0 / total_weight)
    if len(rows) > 1:
        variance = sum((row["loading"] - mean) ** 2 for row in rows) / (len(rows) - 1)
        between = math.sqrt(variance / len(rows))
    else:
        between = 0.0
    total = math.sqrt(within * within + between * between)
    return {
        "replicas": len(rows),
        "mean_loading": mean,
        "within_standard_error": within,
        "between_standard_error": between,
        "standard_error": total,
        "mean_energy_per_molecule_kcal_mol": sum(
            row["energy_per_molecule"] for row in rows
        )
        / len(rows),
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--set", action="append", required=True, help="glob or tag pattern")
    parser.add_argument("--reference", type=float, default=116.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    result = {"reference": args.reference, "sets": {}}
    for item in args.set:
        name, _, pattern = item.partition("=") if "=" in item else (item, "", item)
        rows = load_replicas(pattern)
        if not rows:
            continue
        combined = combine(rows)
        combined["relative_error"] = abs(combined["mean_loading"] - args.reference) / args.reference
        result["sets"][name] = combined

    text = json.dumps(result, indent=2)
    if args.output:
        args.output.write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
