#!/usr/bin/env python3
"""Generate `src/materials/space_group_table.rs` from the committed reference.

The table holds the 230 ITC space groups with their symmetry operations.
Its provenance is the spglib/pymatgen reference data committed as
`tests/fixtures/validation/space_groups_230.json`; this script transcribes
that fixture into the Rust source so the table can never drift silently
(the audit flagged the missing generator as P2-57).

The generated file is verified against the same reference by the F64 gate
(`tests/validation_gates.rs::f64_space_groups_match_spglib`).

Usage:
    python3 scripts/validation/gen_space_group_table.py            # rewrite the table
    python3 scripts/validation/gen_space_group_table.py --check    # verify it is current
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
FIXTURE = REPO / "tests" / "fixtures" / "validation" / "space_groups_230.json"
OUT = REPO / "src" / "materials" / "space_group_table.rs"

HEADER = """//! Generated space-group table (all 230 ITC groups).
//!
//! Generated from pymatgen/spglib reference data (ITC settings); do not edit
//! by hand. See `scripts/validation/gen_space_group_table.py` and the F64
//! gate in `tests/validation_gates.rs`.

use super::space_groups::LatticeType;

/// Compact symmetry operation: integer rotation (row-major) and
/// translation numerators in twelfths of the lattice vector.
pub(crate) type OpCode = ([[i8; 3]; 3], [i8; 3]);

/// Per-group data: symbol, point group, lattice centering, operations.
pub(crate) struct GroupTableEntry {
    pub symbol: &'static str,
    pub point_group: &'static str,
    pub lattice: LatticeType,
    pub ops: &'static [OpCode],
}
"""


def translation_to_twelfths(translation: list[float]) -> list[int]:
    """Convert a fractional translation to integer twelfths.

    The Rust table stores translations as multiples of 1/12, so 0.5 becomes
    6 and 1/3 becomes 4. A translation that is not a multiple of 1/12 is a
    hard error: the table cannot represent it.
    """
    out = []
    for t in translation:
        twelfths = t * 12.0
        rounded = round(twelfths)
        if abs(twelfths - rounded) > 1e-6:
            raise SystemExit(
                f"translation {t} is not a multiple of 1/12 (got {twelfths})"
            )
        out.append(int(rounded))
    return out


def render_ops(group: dict) -> str:
    rows = []
    for op in group["operations"]:
        # The rotation is a 3x3 matrix: `[row, row, row]`.
        matrix_rows = [
            "[" + ", ".join(str(int(v)) for v in row) + "]" for row in op["rotation"]
        ]
        rot = "[" + ", ".join(matrix_rows) + "]"
        trans = ", ".join(str(v) for v in translation_to_twelfths(op["translation"]))
        rows.append(f"    ({rot}, [{trans}]),")
    body = "\n".join(rows)
    return f"const OPS_{group['number']}: &[OpCode] = &[\n{body}\n];\n"


def render() -> str:
    data = json.loads(FIXTURE.read_text())
    groups = data["groups"]
    if len(groups) != 230:
        raise SystemExit(f"expected 230 groups, got {len(groups)}")

    parts = [HEADER]

    for group in groups:
        parts.append(f"\n// #{group['number']} {group['symbol']}\n")
        parts.append(render_ops(group))

    parts.append(
        "\n/// Look up the generated data for an ITC space-group number.\n"
        "pub(crate) fn group_table(number: u16) -> Option<GroupTableEntry> {\n"
    )
    parts.append("    let entry = match number {\n")
    for group in groups:
        parts.append(
            f"        {group['number']} => GroupTableEntry {{\n"
            f'            symbol: "{group["symbol"]}",\n'
            f'            point_group: "{group["point_group"]}",\n'
            f"            lattice: LatticeType::{group['lattice']},\n"
            f"            ops: OPS_{group['number']},\n"
            f"        }},\n"
        )
    parts.append("        _ => return None,\n")
    parts.append("    };\n    Some(entry)\n}\n")

    return "".join(parts)


def format_rust(text: str) -> str:
    """Format the generated Rust with rustfmt so the output is stable."""
    import subprocess
    import tempfile

    with tempfile.NamedTemporaryFile("w", suffix=".rs", delete=False) as fh:
        fh.write(text)
        tmp = fh.name
    try:
        subprocess.run(
            ["rustfmt", "--edition", "2021", tmp],
            check=True,
            capture_output=True,
        )
        return Path(tmp).read_text()
    except (OSError, subprocess.CalledProcessError) as exc:
        raise SystemExit(f"rustfmt failed: {exc}")
    finally:
        Path(tmp).unlink(missing_ok=True)


def main() -> int:
    text = format_rust(render())
    if "--check" in sys.argv:
        current = OUT.read_text() if OUT.exists() else ""
        if current != text:
            print(
                f"{OUT.relative_to(REPO)} is out of date — "
                f"run scripts/validation/gen_space_group_table.py"
            )
            return 1
        print(f"{OUT.relative_to(REPO)} is up to date")
        return 0

    OUT.write_text(text)
    print(f"wrote {OUT.relative_to(REPO)} ({text.count('GroupTableEntry {')} groups)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
