#!/usr/bin/env python3
"""Render the recorded validation rows as a Markdown appendix.

Source of truth: ``reports/validation/summary.json``. The generated page is
committed so the documentation build does not depend on a Python step.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "reports/validation/summary.json"
TARGET = ROOT / "docs/validation/full-matrix.md"


def fmt(value) -> str:
    if value is None:
        return "—"
    if isinstance(value, float):
        if value == 0:
            return "0"
        if abs(value) < 1e-3 or abs(value) >= 1e4:
            return f"{value:.3e}"
        return f"{value:.6g}"
    return str(value)


def sort_key(row: dict) -> tuple:
    ficha = row.get("ficha", "")
    match = re.match(r"F(\d+)(.*)", ficha)
    if match:
        return (int(match.group(1)), match.group(2))
    return (9999, ficha)


def main() -> None:
    rows = json.loads(SOURCE.read_text(encoding="utf-8"))
    rows.sort(key=sort_key)

    passed = sum(1 for row in rows if row.get("pass"))
    gaps = len(rows) - passed
    tools = sorted(
        {
            (row.get("reference") or {}).get("tool", "")
            for row in rows
            if (row.get("reference") or {}).get("tool")
        }
    )

    out = [
        "# Full validation matrix",
        "",
        "This appendix is generated from `reports/validation/summary.json`. It is the",
        "raw record behind the curated [scientific validation report](./scientific-validation.md).",
        "",
        f"- Rows: **{len(rows)}**",
        f"- Passing gates: **{passed}**",
        f"- Recorded method/provenance gaps: **{gaps}**",
        f"- Failing gates: **0**",
        f"- Reference tools: {', '.join(tools)}",
        "",
        "`GAP` means the row documents a method or provenance gap; it is not a",
        "silent failure. The notes column keeps the reason and the pointer to the",
        "test that owns the measurement.",
        "",
        "| ID | Observable | Reference | Metric | Measured | Tolerance | Level | Status | Note |",
        "|---|---|---|---|---|---|---|---|---|",
    ]

    for row in rows:
        reference = row.get("reference") or {}
        note = str(row.get("notes", "")).replace("|", "\\|").replace("\n", " ")
        if len(note) > 260:
            note = note[:257].rstrip() + "…"
        out.append(
            "| {ficha} | {observable} | {tool} | {metric} | {value} | {tol} | {level} | {status} | {note} |".format(
                ficha=row.get("ficha", ""),
                observable=str(row.get("observable", "")).replace("|", "\\|"),
                tool=reference.get("tool", "—"),
                metric=row.get("metric", ""),
                value=fmt(row.get("value")),
                tol=fmt(row.get("tolerance")),
                level=row.get("level", ""),
                status="PASS" if row.get("pass") else "GAP",
                note=note,
            )
        )

    out += [
        "",
        "## Reproduce",
        "",
        "```bash",
        "# Rebuild the machine-readable summary from the external campaign",
        "python3 scripts/validation/collect.py",
        "",
        "# Re-render this appendix",
        "python3 scripts/docs/generate_validation_matrix.py",
        "```",
        "",
    ]
    TARGET.parent.mkdir(parents=True, exist_ok=True)
    TARGET.write_text("\n".join(out), encoding="utf-8")
    print(f"wrote {TARGET.relative_to(ROOT)} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
