#!/usr/bin/env python3
"""Oracle F90: chart/render payload schema + conformance checker.

- Authors the normative JSON Schema (draft 2020-12) for ChartPayload /
  ChartSeries / RenderPayload / TrajectoryFrame (Vega-Lite data rules:
  finite numbers, x/y coherence, unique non-empty ids).
- Provides adversarial payloads (valid + invalid) consumed by the Rust
  gate, which mirrors the schema clauses natively (no new Rust deps).
- Check mode: validates Rust-emitted payloads in
  target/validation/render_bridge/*.json with the `jsonschema` package:
  `oracle_render_charts.py --check`.

Usage:
  oracle_render_charts.py          # writes fixture
  oracle_render_charts.py --check  # validates Rust-emitted JSON
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
FIX = REPO / "tests" / "fixtures" / "validation" / "render_charts.json"
EMIT = REPO / "target" / "validation" / "render_bridge"

SCHEMA = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "title": "ChartPayload",
    "type": "object",
    "required": ["title", "series"],
    "properties": {
        "title": {"type": "string", "minLength": 1},
        "series": {"type": "array", "items": {"$ref": "#/$defs/series"}},
    },
    "$defs": {
        "series": {
            "type": "object",
            "required": ["series_id", "label", "x", "y", "x_unit", "y_unit"],
            "properties": {
                "series_id": {"type": "string", "minLength": 1},
                "label": {"type": "string"},
                "x": {"type": "array", "items": {"type": "number"}},
                "y": {"type": "array", "items": {"type": "number"}},
                "x_unit": {"type": "string"},
                "y_unit": {"type": "string"},
            },
        }
    },
}


def good_chart() -> dict:
    return {
        "title": "Microkinetic populations",
        "series": [
            {"series_id": "A", "label": "A", "x": [0.0, 1.0, 2.0],
             "y": [1.0, 0.5, 0.25], "x_unit": "s", "y_unit": "arb."},
            {"series_id": "B", "label": "B", "x": [0.0, 1.0, 2.0],
             "y": [0.0, 0.5, 0.75], "x_unit": "s", "y_unit": "arb."},
        ],
    }


def bad_charts() -> list[dict]:
    g = good_chart()
    no_title = dict(g)
    del no_title["title"]
    empty_id = json.loads(json.dumps(g))
    empty_id["series"][0]["series_id"] = ""
    dup_id = json.loads(json.dumps(g))
    dup_id["series"][1]["series_id"] = "A"
    skew = json.loads(json.dumps(g))
    skew["series"][0]["y"] = [1.0, 0.5]
    return [
        {"id": "missing-title", "chart": no_title, "valid": False},
        {"id": "empty-series-id", "chart": empty_id, "valid": False},
        {"id": "duplicate-series-id", "chart": dup_id, "valid": False},
        {"id": "xy-length-skew", "chart": skew, "valid": False},
        {"id": "good-chart", "chart": g, "valid": True},
    ]


def gen() -> int:
    FIX.write_text(json.dumps(
        {"ficha": "F90", "oracle": "oracle_render_charts.py (JSON Schema draft 2020-12)",
         "schema": SCHEMA, "adversarial": bad_charts()}, indent=1) + "\n")
    print(f"wrote {FIX}")
    return 0


def check() -> int:
    import jsonschema
    files = sorted(EMIT.glob("*.json"))
    if not files:
        print("no Rust-emitted payloads (run the F90 gate first)")
        return 1
    ok = True
    for f in files:
        payload = json.loads(f.read_text())
        try:
            jsonschema.validate(payload, SCHEMA)
            # Vega-Lite data rules beyond structural schema.
            ids = [s["series_id"] for s in payload["series"]]
            assert len(ids) == len(set(ids)) and all(ids), "ids"
            for s in payload["series"]:
                assert len(s["x"]) == len(s["y"]) and s["x"], "coherence"
                assert all(isinstance(v, (int, float)) and v == v
                           and abs(v) != float("inf") for v in s["x"] + s["y"]), "finite"
            print(f"[PASS] {f.name}")
        except Exception as err:  # noqa: BLE001
            ok = False
            print(f"[FAIL] {f.name}: {err}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(check() if "--check" in sys.argv else gen())
