#!/usr/bin/env python3
"""Provenance check for sci-form's transcribed CSD torsion tables (F76).

`src/smarts/torsion_data.rs` transcribes RDKit's torsion-preference tables:

    Code/GraphMol/ForceFieldHelpers/CrystalFF/torsionData/
        torsionPreferences_v2.in            -> TORSION_PREFS_V2
        torsionPreferences_smallrings.in     -> TORSION_PREFS_SMALL_RINGS
        torsionPreferences_macrocycles.in    -> TORSION_PREFS_MACROCYCLES

RDKit's Python wheel does **not** ship the `.in` files (they are compiled into
the library), so this script

  1. looks for the files inside the installed RDKit package,
  2. falls back to the pinned upstream copy in the RDKit git repository at the
     tag that matches the installed RDKit version, verifying its SHA-256,
  3. parses the upstream tables and `src/smarts/torsion_data.rs`,
  4. reports a per-pattern diff,
  5. writes `tests/fixtures/validation/csd_torsion_table.json` with the
     SHA-256 of the upstream file *and* of the transcribed table.

The Rust side (`tests/rdkit_parity.rs::f76_csd_torsion_table`) pins the
transcribed SHA-256 so that any drift of `torsion_data.rs` is caught in CI
without needing RDKit or network access.

Usage:
    python scripts/validation/check_csd_torsion_table.py [--upstream-dir DIR]
"""

from __future__ import annotations

import argparse
import hashlib
import re
import struct
import sys
import time
import urllib.request
from pathlib import Path

from common import REPO, Report, emit_fixture_report

#: RDKit release tag matching the installed RDKit version.
RDKIT_TAG = "Release_2026_09_1"
UPSTREAM_BASE = (
    "https://raw.githubusercontent.com/rdkit/rdkit/"
    f"{RDKIT_TAG}/Code/GraphMol/ForceFieldHelpers/CrystalFF/torsionData"
)

#: SHA-256 of the three upstream files at RDKIT_TAG (pinned so that a silent
#: upstream change is detected instead of silently accepted).
EXPECTED_UPSTREAM_SHA256 = {
    "torsionPreferences_v2.in": "c1c1ff3d93b9629d1e71502daa02257c9968ca4788478d4970f81e457a37d8e5",
    "torsionPreferences_smallrings.in": "28124c740e9aba1ec835fa9cef8f095258f0b94f5829dda2b995e99bbd1b5efb",
    "torsionPreferences_macrocycles.in": "0a11044bbcdddf842b505e4baa0b4878e1fcc0725fa7c99dd82ff462938eab6e",
}

#: Upstream file -> sci-form table name.
TABLES = {
    "torsionPreferences_v2.in": "TORSION_PREFS_V2",
    "torsionPreferences_smallrings.in": "TORSION_PREFS_SMALL_RINGS",
    "torsionPreferences_macrocycles.in": "TORSION_PREFS_MACROCYCLES",
}

#: Order in which the tables are serialised into the transcribed digest.
TABLE_ORDER = [
    "TORSION_PREFS_V2",
    "TORSION_PREFS_SMALL_RINGS",
    "TORSION_PREFS_EXTENSIONS",
    "TORSION_PREFS_MACROCYCLES",
]

RUST_FILE = REPO / "src" / "smarts" / "torsion_data.rs"


# ---------------------------------------------------------------------------
# Locating the upstream files
# ---------------------------------------------------------------------------


def find_in_installed_rdkit() -> dict[str, Path]:
    """Search the installed RDKit package for the `.in` files."""
    found: dict[str, Path] = {}
    try:
        import rdkit  # noqa: PLC0415
        from rdkit import RDConfig  # noqa: PLC0415
    except ImportError:
        return found
    roots = [Path(RDConfig.RDDataDir), Path(rdkit.__file__).parent]
    for root in roots:
        if not root.exists():
            continue
        for path in root.rglob("torsionPreferences*.in"):
            found[path.name] = path
    return found


def fetch_upstream(name: str, cache_dir: Path) -> Path:
    """Download the pinned upstream file (cached) and verify its SHA-256."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    target = cache_dir / name
    if not target.exists():
        url = f"{UPSTREAM_BASE}/{name}"
        with urllib.request.urlopen(url, timeout=60) as resp:  # noqa: S310
            target.write_bytes(resp.read())
    digest = hashlib.sha256(target.read_bytes()).hexdigest()
    expected = EXPECTED_UPSTREAM_SHA256[name]
    if digest != expected:
        raise SystemExit(
            f"upstream {name} SHA-256 mismatch: got {digest}, expected {expected}\n"
            f"(the pinned tag {RDKIT_TAG} no longer matches; update "
            f"EXPECTED_UPSTREAM_SHA256 deliberately)"
        )
    return target


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

_CPP_STRING = re.compile(r'"((?:[^"\\]|\\.)*)"')


def parse_upstream_in(path: Path) -> list[tuple[str, list[float]]]:
    """Parse an RDKit `torsionPreferences*.in` file.

    Each entry is a C string literal holding `"<SMARTS> s1 V1 ... s6 V6\\n"`.
    """
    text = path.read_text()
    entries: list[tuple[str, list[float]]] = []
    for raw in _CPP_STRING.findall(text):
        line = raw.replace("\\n", "").replace("\\t", " ").strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) != 13:
            raise SystemExit(f"{path.name}: malformed entry {line!r}")
        smarts = parts[0]
        values = [float(p) for p in parts[1:]]
        entries.append((smarts, values))
    return entries


_RUST_TABLE = re.compile(
    r"pub const (\w+): &\[\(&str, \[f32; 12\]\)\] = &\[(.*?)\n\];",
    re.DOTALL,
)
_RUST_ENTRY = re.compile(
    r"\(\s*\"((?:[^\"\\]|\\.)*)\"\s*,\s*\[([^\]]*)\]\s*,?\s*\)", re.DOTALL
)


def parse_rust_tables(path: Path) -> dict[str, list[tuple[str, list[float]]]]:
    """Parse the three transcribed tables out of `src/smarts/torsion_data.rs`."""
    text = path.read_text()
    tables: dict[str, list[tuple[str, list[float]]]] = {}
    for name, body in _RUST_TABLE.findall(text):
        entries: list[tuple[str, list[float]]] = []
        for pattern, nums in _RUST_ENTRY.findall(body):
            values = [float(v.strip()) for v in nums.split(",") if v.strip()]
            if len(values) != 12:
                raise SystemExit(f"{path.name}: {name} entry with {len(values)} values")
            entries.append((pattern, values))
        tables[name] = entries
    return tables


def f32_bits(value: float) -> int:
    """IEEE-754 single precision bit pattern of `value` (round to nearest)."""
    return struct.unpack("<I", struct.pack("<f", value))[0]


def canonical_lines(tables: dict[str, list[tuple[str, list[float]]]]) -> list[str]:
    """Canonical serialisation of the transcribed tables.

    One line per entry: `<table>\\t<smarts>\\t<b0>,<b1>,...,<b11>` where each
    `bi` is the 8-hex-digit IEEE-754 bit pattern of the f32 value.  The Rust
    test reproduces exactly this serialisation, so the SHA-256 of the joined
    lines is a stable fingerprint of the transcribed data.
    """
    lines: list[str] = []
    for table in TABLE_ORDER:
        for pattern, values in tables[table]:
            bits = ",".join(f"{f32_bits(v):08x}" for v in values)
            lines.append(f"{table}\t{pattern}\t{bits}")
    return lines


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--upstream-dir",
        type=Path,
        default=None,
        help="directory holding the upstream torsionPreferences*.in files",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path("/tmp/opencode/rdkit-torsion-cache"),
        help="cache directory for downloaded upstream files",
    )
    args = parser.parse_args()

    installed = find_in_installed_rdkit()
    source = {}
    for name in TABLES:
        if args.upstream_dir is not None:
            path = args.upstream_dir / name
            if not path.exists():
                raise SystemExit(f"missing upstream file {path}")
        elif name in installed:
            path = installed[name]
        else:
            path = fetch_upstream(name, args.cache_dir)
        source[name] = path

    upstream = {name: parse_upstream_in(path) for name, path in source.items()}
    rust_tables = parse_rust_tables(RUST_FILE)

    diffs: list[dict] = []
    n_upstream = 0
    n_transcribed = 0
    n_differing = 0
    for name, table in TABLES.items():
        up = upstream[name]
        mine = rust_tables[table]
        n_upstream += len(up)
        n_transcribed += len(mine)
        up_map: dict[str, list[float]] = {}
        for pattern, values in up:
            up_map.setdefault(pattern, values)
        my_map: dict[str, list[float]] = {}
        for pattern, values in mine:
            my_map.setdefault(pattern, values)
        only_up = sorted(set(up_map) - set(my_map))
        only_mine = sorted(set(my_map) - set(up_map))
        changed = []
        for pattern in sorted(set(up_map) & set(my_map)):
            a = [f32_bits(v) for v in up_map[pattern]]
            b = [f32_bits(v) for v in my_map[pattern]]
            if a != b:
                changed.append(
                    {
                        "smarts": pattern,
                        "upstream": up_map[pattern],
                        "transcribed": my_map[pattern],
                    }
                )
        n_differing += len(changed) + len(only_up) + len(only_mine)
        diffs.append(
            {
                "file": name,
                "table": table,
                "n_upstream": len(up),
                "n_transcribed": len(mine),
                "only_in_upstream": only_up,
                "only_in_transcribed": only_mine,
                "different_parameters": changed,
            }
        )

    # sci-form specific patterns that are deliberately *not* part of RDKit's
    # tables.  They live in their own table so that the three upstream tables
    # stay faithful transcriptions.
    extensions = rust_tables.get("TORSION_PREFS_EXTENSIONS", [])
    n_extensions = len(extensions)

    lines = canonical_lines(rust_tables)
    transcribed_sha = hashlib.sha256("\n".join(lines).encode("utf-8")).hexdigest()
    upstream_sha = {
        name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in source.items()
    }
    for name, digest in upstream_sha.items():
        if digest != EXPECTED_UPSTREAM_SHA256[name]:
            raise SystemExit(f"{name}: unexpected upstream SHA-256 {digest}")

    payload = {
        "reference": f"RDKit {RDKIT_TAG} Code/GraphMol/ForceFieldHelpers/CrystalFF/torsionData",
        "rdkit_tag": RDKIT_TAG,
        "upstream_source": (
            "installed rdkit package"
            if installed
            else f"{UPSTREAM_BASE} (pinned tag; the RDKit wheel does not ship the .in files)"
        ),
        "upstream_sha256": upstream_sha,
        "transcribed_sha256": transcribed_sha,
        "transcribed_serialisation": (
            "one line per entry: '<table>\\t<smarts>\\t<f32 bit patterns, comma separated>'; "
            "the SHA-256 is taken over the lines joined with '\\n'"
        ),
        "n_upstream_entries": n_upstream,
        "n_transcribed_entries": n_transcribed,
        "n_extension_entries": n_extensions,
        "n_differing_entries": n_differing,
        "extension_patterns": [p for p, _ in extensions],
        "tables": diffs,
        "generated": time.strftime("%Y-%m-%d"),
    }
    report = Report(
        ficha="F76",
        observable="csd_torsion_table_drift",
        reference_tool=f"RDKit {RDKIT_TAG} torsionPreferences*.in",
        units="count",
        tolerance=0,
        level="L2",
        n_cases=n_upstream,
        metric="differing_entries",
        value=float(n_differing),
        notes=(
            "the three transcribed CSD torsion tables must match RDKit's upstream "
            "torsionPreferences*.in entry for entry; the sci-form specific "
            "patterns are kept in a separate TORSION_PREFS_EXTENSIONS table; the "
            "SHA-256 of the transcribed tables is pinned in the Rust test"
        ),
    )
    emit_fixture_report("csd_torsion_table.json", payload, report)

    print(f"upstream entries   : {n_upstream}")
    print(f"transcribed entries: {n_transcribed}")
    print(f"extension entries  : {n_extensions}")
    print(f"differing entries  : {n_differing}")
    print(f"transcribed sha256 : {transcribed_sha}")
    for name, digest in upstream_sha.items():
        print(f"  {name:38s} {digest}")
    for d in diffs:
        print(
            f"  {d['table']:26s} upstream={d['n_upstream']:4d} "
            f"transcribed={d['n_transcribed']:4d} "
            f"only_up={len(d['only_in_upstream'])} "
            f"only_mine={len(d['only_in_transcribed'])} "
            f"changed={len(d['different_parameters'])}"
        )
    return 0 if n_differing == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
