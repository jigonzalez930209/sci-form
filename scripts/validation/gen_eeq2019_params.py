#!/usr/bin/env python3
"""Generate `src/charges_eeq/eeq2019_data.rs` from the multicharge sources.

Inputs:
  - `src/multicharge/param/eeq2019.f90` (grimme-lab/multicharge, main):
    `eeq_chi`, `eeq_eta`, `eeq_kcnchi`, `eeq_rad` for Z = 1..103.
  - `src/dftd4/data/covrad.f90` (dftd4 v4.3.0): Pyykko covalent radii in
    Angstrom, used for the EEQ coordination-number cutoff.

The generated values for Z <= 86 are verified against the previously committed
table (same source, same units).

Usage:
    python scripts/validation/gen_eeq2019_params.py [--source-dir DIR] [--dftd4-dir DIR]
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

from refdata import ref_data_dir, skip

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "src" / "charges_eeq" / "eeq2019_data.rs"
MAX_ELEM = 103  # multicharge's eeq2019 coverage


def _floats(text: str) -> list[float]:
    return [float(x) for x in re.findall(r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?", text)]


def parse_multicharge_array(path: Path, name: str) -> list[float]:
    text = path.read_text()
    m = re.search(rf"{name}\s*\(max_elem\)\s*=\s*\[&(.*?)\]", text, re.S)
    if not m:
        raise ValueError(f"{name} not found in {path}")
    vals = _floats(m.group(1))
    if len(vals) != MAX_ELEM:
        raise ValueError(f"{name}: parsed {len(vals)} values, expected {MAX_ELEM}")
    return vals


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--source-dir",
        type=Path,
        default=ref_data_dir() / "multicharge",
        help="grimme-lab/multicharge checkout (default: $SCIFORM_REF_DATA/multicharge)",
    )
    ap.add_argument(
        "--dftd4-dir",
        type=Path,
        default=ref_data_dir() / "dftd4-430",
        help="dftd4 v4.3.0 source tree (default: $SCIFORM_REF_DATA/dftd4-430)",
    )
    ap.add_argument(
        "--check",
        action="store_true",
        help="regenerate into a temporary file, compare byte-for-byte with the "
             "committed table and exit non-zero on a difference",
    )
    args = ap.parse_args()

    missing = [
        p
        for p in (
            args.source_dir / "eeq2019.f90",
            args.dftd4_dir / "covrad.f90",
        )
        if not p.exists()
    ]
    if missing:
        msg = (
            f"the upstream sources are not available ({', '.join(str(p) for p in missing)}); "
            f"set SCIFORM_REF_DATA or pass --source-dir/--dftd4-dir"
        )
        if args.check:
            return skip(msg)
        raise SystemExit(f"error: {msg}")

    chi = parse_multicharge_array(args.source_dir / "eeq2019.f90", "eeq_chi")
    eta = parse_multicharge_array(args.source_dir / "eeq2019.f90", "eeq_eta")
    kcnchi = parse_multicharge_array(args.source_dir / "eeq2019.f90", "eeq_kcnchi")
    rad = parse_multicharge_array(args.source_dir / "eeq2019.f90", "eeq_rad")

    cov_src = (args.dftd4_dir / "covrad.f90").read_text()
    m = re.search(r"covalent_rad_2009\(max_elem\)\s*=\s*aatoau\s*\*\s*\[(.*?)\]", cov_src, re.S)
    covrad = _floats(m.group(1))[:MAX_ELEM]
    assert len(covrad) == MAX_ELEM

    # --- provenance check against the committed table -----------------------
    old = OUT.read_text()

    def committed(name: str) -> list[float] | None:
        m = re.search(rf"{name}: \[f64; 86\] = \[(.*?)\];", old, re.S)
        return _floats(m.group(1)) if m else None

    for name, values in [
        ("EEQ2019_CHI", chi),
        ("EEQ2019_ETA", eta),
        ("EEQ2019_KCNCHI", kcnchi),
        ("EEQ2019_RAD", rad),
        ("EEQ2019_COVRAD_ANG", covrad),
    ]:
        prev = committed(name)
        if prev is None:
            continue
        for z in range(86):
            if abs(prev[z] - values[z]) > 1e-9:
                print(f"[warn] {name}[Z={z + 1}]: committed {prev[z]} vs source {values[z]}")

    def emit(name: str, values: list[float], comment: str, per_line: int = 6) -> str:
        lines = [
            f"/// {comment}",
            "#[rustfmt::skip]",
            f"pub(crate) const {name}: [f64; {len(values)}] = [",
        ]
        for i in range(0, len(values), per_line):
            chunk = ", ".join(f"{v:.8f}" for v in values[i : i + per_line])
            lines.append(f"    {chunk},")
        lines.append("];")
        return "\n".join(lines) + "\n"

    out = (
        "//! EEQ2019 charge-model parameters (generated from multicharge).\n"
        "//!\n"
        "//! Source: grimme-lab/multicharge `src/multicharge/param/eeq2019.f90`\n"
        "//! (Apache-2.0) plus the Pyykko covalent radii from dftd4 v4.3.0.\n"
        "//! Coverage: Z = 1..103 (multicharge's eeq2019 range).\n"
        "//!\n"
        "//! Regenerate with `scripts/validation/gen_eeq2019_params.py`; do not\n"
        "//! edit by hand.\n\n"
    )
    out += emit("EEQ2019_CHI", chi, "Element-specific electronegativity (atomic units).")
    out += "\n" + emit("EEQ2019_ETA", eta, "Chemical hardness (atomic units).")
    out += "\n" + emit("EEQ2019_RAD", rad, "EEQ radius (Bohr).")
    out += "\n" + emit("EEQ2019_KCNCHI", kcnchi, "Coordination-number scaling of chi.")
    out += "\n" + emit(
        "EEQ2019_COVRAD_ANG",
        covrad,
        "Pyykko covalent radii (Angstrom) for the CN cutoff.",
    )
    if args.check:
        current = OUT.read_text() if OUT.exists() else ""
        if current != out:
            print(
                f"{OUT.relative_to(REPO)} is out of date — "
                f"run scripts/validation/gen_eeq2019_params.py"
            )
            return 1
        print(f"{OUT.relative_to(REPO)} is up to date")
        return 0

    OUT.write_text(out)
    print(f"wrote {OUT} (MAX_ELEM={MAX_ELEM})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
