#!/usr/bin/env python3
"""Merge the metrics written by tests/cross_validation.rs into the summary.

The oracle scripts write placeholder rows (fixture provenance); the Rust tests
write `target/validation/{ficha}.json` with the measured values. This script
combines both into `reports/validation/summary.json`.
"""

from __future__ import annotations

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TARGET = REPO / "target" / "validation"
SUMMARY = REPO / "reports" / "validation" / "summary.json"


def main() -> int:
    summary = []
    if SUMMARY.exists():
        try:
            summary = json.loads(SUMMARY.read_text())
        except json.JSONDecodeError:
            summary = []
    by_key = {(r.get("ficha"), r.get("observable")): r for r in summary}

    for path in sorted(TARGET.glob("F*.json")):
        data = json.loads(path.read_text())
        file_rows = data if isinstance(data, list) else [data]
        for row in file_rows:
            key = (row["ficha"], row["observable"])
            base = dict(by_key.get(key, {}))
            base.update(
                {
                    "ficha": row["ficha"],
                    "observable": row["observable"],
                    "metric": row["metric"],
                    "value": row["value"],
                    "tolerance": row["tolerance"],
                    "pass": row["pass"],
                    "n_cases": row["n_cases"],
                    "source": row.get("source", "tests/cross_validation.rs"),
                }
            )
            base.setdefault("reference", {})
            base.setdefault("level", "L2")
            by_key[key] = base

    rows = sorted(by_key.values(), key=lambda r: (r.get("ficha", ""), r.get("observable", "")))
    for row in rows:
        if row.get("value") is None:
            # Fixture provenance rows whose measurement lives in the Rust tests,
            # or strict gates that are currently tracked as ignored tests.
            row["pass"] = None
            row["status"] = "tracked"
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY.write_text(json.dumps(rows, indent=1))

    passed = [r for r in rows if r.get("pass") is True]
    failed = [r for r in rows if r.get("pass") is False]
    pending = [r for r in rows if r.get("pass") is None]
    print(f"validation summary: {len(rows)} reports, {len(passed)} pass, "
          f"{len(failed)} fail, {len(pending)} pending (Rust-side metrics)")
    for r in rows:
        flag = "PASS" if r.get("pass") is True else ("FAIL" if r.get("pass") is False else "----")
        value = r.get("value")
        value_s = "n/a" if value is None else f"{value:.6g}"
        print(
            f"  [{flag}] {r['ficha']:>3} {r.get('observable',''):<48} "
            f"value={value_s:<12} tol={r.get('tolerance')} n={r.get('n_cases')}"
        )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
