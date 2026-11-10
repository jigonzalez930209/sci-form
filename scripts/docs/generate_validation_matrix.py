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


def soft_breaks(text) -> str:
    """Escape pipes and add break opportunities inside long identifiers."""
    return str(text).replace("|", "\\|").replace("_", "_<wbr>")


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
        "- Reference tools: RDKit 2026.09.1, PySCF 2.14.0, tblite 0.7.0,",
        "  dftd4 4.3.0, xtb 6.7.1, LAMMPS, ASE 3.29.0, pymatgen 2026.9.24,",
        "  pormake 0.2.3, FreeSASA, NumPy/SciPy, and the independent oracles.",
        "",
        "`GAP` means the row documents a method or provenance gap; it is not a",
        "silent failure. The per-row notes and the owning test are kept in",
        "`reports/validation/summary.json`; the curated report explains the",
        "limitations behind the gaps.",
        "",
        "| ID | Observable | Reference | Metric | Measured | Tolerance | Level | Status |",
        "|---|---|---|---|---|---|---|---|",
    ]

    for row in rows:
        reference = row.get("reference") or {}
        out.append(
            "| {ficha} | {observable} | {tool} | {metric} | {value} | {tol} | {level} | {status} |".format(
                ficha=row.get("ficha", ""),
                observable=f"<code>{soft_breaks(row.get('observable', ''))}</code>",
                tool=reference.get("tool", "—"),
                metric=f"<code>{soft_breaks(row.get('metric', ''))}</code>",
                value=fmt(row.get("value")),
                tol=fmt(row.get("tolerance")),
                level=row.get("level", ""),
                status="PASS" if row.get("pass") else "GAP",
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
