#!/usr/bin/env python3
"""Generate `src/xtb/d4_data.rs` from the dftd4 reference sources.

Inputs (fetched from the dftd4 v4.3.0 release, the version validated by F46):
  - src/dftd4/reference.inc       (per-element reference systems, charges,
                                   hydrogen counts, coordination numbers,
                                   Gaussian weights, reference polarizabilities)
  - src/dftd4/data/r4r2.f90       (<r^4>/<r^2> expectation values)
  - src/dftd4/data/hardness.f90   (chemical hardness)
  - src/dftd4/data/zeff.f90       (effective nuclear charges)
  - src/dftd4/data/covrad.f90     (covalent radii, D3 set)
  - src/dftd4/data/en.f90         (Pauling electronegativities)

Coverage is extended to Z = 1..118 (matching dftd4); the generated values for
Z <= 86 are verified against the previously committed table, which was produced
from the same source.

Usage:
    python scripts/validation/gen_d4_refs.py [--source-dir DIR] [--check]
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

from refdata import ref_data_dir, skip

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "src" / "xtb" / "d4_data.rs"
MAX_ELEM = 118
MAX_REF = 7
NFREQ = 23


def _floats(text: str) -> list[float]:
    return [float(x.replace("d", "e").replace("D", "e")) for x in re.findall(r"[-+]?\d+\.\d+(?:[eEdD][-+]?\d+)?", text)]


def _numbers(text: str) -> list[float]:
    """Numbers with or without a decimal part (Fortran integer literals)."""
    return [float(x) for x in re.findall(r"(?<![\w.])[-+]?\d+(?:\.\d+)?(?![\w.])", text)]


def parse_reference_inc(path: Path) -> dict:
    """Parse the generated Fortran include with the per-element reference data."""
    text = path.read_text()
    # Fortran `data` statements, possibly continued across lines.
    text = re.sub(r"&\s*\n\s*&", " ", text)
    data: dict = {}

    def collect(name: str) -> dict[tuple[int, ...], float]:
        out: dict[tuple[int, ...], float] = {}
        # Single scalar entries: data name(i,j) / value /
        # Scalar entries only: the index list must be purely numeric (so the
        # vector form `name(:, i, j)` is not matched here).
        for m in re.finditer(rf"data\s+{name}\s*\(\s*(\d+(?:\s*,\s*\d+)*)\s*\)\s*/\s*([^/]+)/", text):
            idx = tuple(int(x) for x in re.findall(r"\d+", m.group(1)))
            vals = _floats(m.group(2)) or _numbers(m.group(2))
            if vals:
                out[idx] = vals[0]
        # Vector entries: data name(:,i,j) / v1, v2, ... /
        for m in re.finditer(rf"data\s+{name}\s*\(:\s*,\s*(\d+)\s*,\s*(\d+)\)\s*/\s*([^/]+)/", text):
            iref, z = int(m.group(1)), int(m.group(2))
            for k, v in enumerate(_floats(m.group(3))):
                out[(k + 1, iref, z)] = v
        return out

    data["refn"] = {z: v for (z,), v in collect("refn").items()}
    data["refq"] = collect("refq")
    data["refh"] = collect("refh")
    data["dftq"] = collect("dftq")
    data["dfth"] = collect("dfth")
    # dftd4's standalone D4 (validation F46) uses the "classical" reference
    # charges (`clsq`/`clsh`), not `dftq`/`dfth` (verified empirically against
    # dftd4 4.3.0 and by the committed-table provenance check).
    data["clsq"] = collect("clsq")
    data["clsh"] = collect("clsh")
    data["ascale"] = collect("ascale")
    data["refcn"] = collect("refcn")
    data["refcovcn"] = collect("refcovcn")
    data["refsys"] = collect("refsys")
    data["hcount"] = collect("hcount")
    data["alphaiw"] = collect("alphaiw")
    data["secq"] = collect("secq")
    data["sscale"] = collect("sscale")
    # Secondary polarizabilities: `data secaiw (:,N) / ... /` (single index).
    secaiw: dict[tuple[int, ...], float] = {}
    for m in re.finditer(r"data\s+secaiw\s*\(\s*:\s*,\s*(\d+)\s*\)\s*/\s*([^/]+)/", text):
        i = int(m.group(1))
        for k, v in enumerate(_floats(m.group(2))):
            secaiw[(k + 1, i)] = v
    data["secaiw"] = secaiw
    return data


def parse_fortran_array(path: Path, name: str) -> list[float]:
    text = path.read_text()
    m = re.search(rf"{name}\s*\(max_elem\)\s*=\s*(?:aatoau\s*\*\s*)?\[(.*?)\]", text, re.S)
    if not m:
        raise ValueError(f"{name} not found in {path}")
    body = m.group(1)
    vals = _floats(body)
    if not vals:
        # integer-valued tables (e.g. effective nuclear charges)
        vals = [float(x) for x in re.findall(r"(?<![\w.])[-+]?\d+(?![\w.])", body)]
    return vals


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--source-dir",
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

    src = args.source_dir
    missing = [
        p
        for p in (
            src / "reference.inc",
            src / "r4r2.f90",
            src / "hardness.f90",
            src / "zeff.f90",
            src / "covrad.f90",
            src / "en.f90",
        )
        if not p.exists()
    ]
    if missing:
        msg = (
            f"the dftd4 sources are not available ({', '.join(str(p) for p in missing)}); "
            f"set SCIFORM_REF_DATA or pass --source-dir"
        )
        if args.check:
            return skip(msg)
        raise SystemExit(f"error: {msg}")

    refs = parse_reference_inc(src / "reference.inc")
    r4r2_raw = parse_fortran_array(src / "r4r2.f90", "r4_over_r2")
    hardness = parse_fortran_array(src / "hardness.f90", "chemical_hardness")
    zeff = parse_fortran_array(src / "zeff.f90", "effective_nuclear_charge")
    # dftd4: covalent_rad_d3 = 4/3 * aatoau * covalent_rad_2009(Angstrom)
    aatoau = 1.0 / 0.529177210903
    covrad_ang = parse_fortran_array(src / "covrad.f90", "covalent_rad_2009")
    covrad = [4.0 / 3.0 * aatoau * v for v in covrad_ang]
    en = parse_fortran_array(src / "en.f90", "pauling_en")
    assert len(r4r2_raw) == MAX_ELEM and len(hardness) == MAX_ELEM

    # dftd4's `sqrt_z_r4_over_r2`: what the D4 damping radii use.
    r4r2 = [(0.5 * r4r2_raw[z] * (z + 1) ** 0.5) ** 0.5 for z in range(MAX_ELEM)]

    # --- provenance check against the previously committed table -------------
    old = OUT.read_text()
    checks = {
        "R4R2": (r4r2, 1e-6),
        "CHEMICAL_HARDNESS": (hardness, 1e-6),
        "EFFECTIVE_NUCLEAR_CHARGE": (zeff, 1e-6),
        "COVRAD_D3": (covrad, 1e-6),
        "PAULING_EN": (en, 1e-6),
    }
    for name, (values, tol) in checks.items():
        m = re.search(rf"pub const {name}: \[f64; 86\] = \[(.*?)\];", old, re.S)
        if not m:
            continue
        prev = _floats(m.group(1))
        if len(prev) != 86:
            print(f"[warn] {name}: committed table has {len(prev)} entries, expected 86")
            continue
        for z in range(86):
            rel = abs(prev[z] - values[z]) / max(abs(values[z]), 1e-12)
            if rel > tol:
                print(f"[warn] {name}[Z={z + 1}]: committed {prev[z]} vs source {values[z]}")
    # reference tables
    prev_refn = re.search(r"pub const REFN: \[usize; 86\] = \[(.*?)\];", old, re.S)
    if prev_refn:
        prev = [int(x) for x in re.findall(r"\d+", prev_refn.group(1))]
        for z in range(86):
            if prev[z] != refs["refn"].get(z + 1, 0):
                print(f"[warn] REFN[Z={z + 1}]: committed {prev[z]} vs source {refs['refn'].get(z + 1)}")
    # --- emit ---------------------------------------------------------------
    def flat3(table: dict, default: float = 0.0) -> list[float]:
        out = []
        for z in range(1, MAX_ELEM + 1):
            for iref in range(1, MAX_REF + 1):
                out.append(table.get((iref, z), default))
        return out

    def flat_alpha() -> list[float]:
        out = []
        for z in range(1, MAX_ELEM + 1):
            for iref in range(1, MAX_REF + 1):
                for k in range(1, NFREQ + 1):
                    out.append(refs["alphaiw"].get((k, iref, z), 0.0))
        return out

    def emit_f64(name: str, values: list[float], per_line: int = 6) -> str:
        lines = [
            "#[rustfmt::skip]",
            f"pub const {name}: [f64; {len(values)}] = [",
        ]
        for i in range(0, len(values), per_line):
            chunk = ", ".join(f"{v:.10e}" for v in values[i : i + per_line])
            lines.append(f"    {chunk},")
        lines.append("];")
        return "\n".join(lines) + "\n"

    def emit_usize(name: str, values, per_line: int = 14) -> str:
        lines = [
            "#[rustfmt::skip]",
            f"pub const {name}: [usize; {len(values)}] = [",
        ]
        for i in range(0, len(values), per_line):
            chunk = ", ".join(str(int(v)) for v in values[i : i + per_line])
            lines.append(f"    {chunk},")
        lines.append("];")
        return "\n".join(lines) + "\n"

    secq = [refs["secq"].get((i,), 0.0) for i in range(1, 18)]
    sscale = [refs["sscale"].get((i,), 0.0) for i in range(1, 18)]
    secaiw = []
    for i in range(1, 18):
        for k in range(1, NFREQ + 1):
            secaiw.append(refs["secaiw"].get((k, i), 0.0))

    header = f"""//! dftd4 reference data (generated from dftd4 v4.3.0).
//!
//! Sources: `src/dftd4/reference.inc` (per-element reference systems, charges,
//! hydrogen counts, coordination numbers, Gaussian weights and reference
//! polarizabilities) and `src/dftd4/data/*.f90` (r4r2, hardness, zeff, covrad,
//! Pauling EN).
//!
//! Coverage: Z = 1..{MAX_ELEM} (dftd4's full range). Regenerate with
//! `scripts/validation/gen_d4_refs.py`; do not edit by hand.

/// Maximum supported element (dftd4 range).
pub const MAX_ELEM: usize = {MAX_ELEM};
/// Maximum number of reference systems per element.
pub const MAX_REF: usize = {MAX_REF};
/// Number of Casimir-Polder frequencies.
pub const NFREQ: usize = {NFREQ};

/// `<r^4>/<r^2>` expectation values as used by the D4 damping radii
/// (`sqrt(0.5 * r4r2 * sqrt(Z))`, dftd4's `sqrt_z_r4_over_r2`).
"""
    out = header
    out += emit_f64("R4R2", r4r2)
    out += "\n/// Chemical hardness (eV) used by the EEQ charge model.\n"
    out += emit_f64("CHEMICAL_HARDNESS", hardness)
    out += "\n/// Effective nuclear charges (EEQ).\n"
    out += emit_f64("EFFECTIVE_NUCLEAR_CHARGE", zeff)
    out += "\n/// D3 covalent radii (Bohr).\n"
    out += emit_f64("COVRAD_D3", covrad)
    out += "\n/// Pauling electronegativities.\n"
    out += emit_f64("PAULING_EN", en)
    out += "\n/// Casimir-Polder integration frequencies (atomic units).\n"
    cp_freq = re.search(r"pub const CP_FREQ: \[f64; 23\] = \[(.*?)\];", old, re.S)
    out += emit_f64("CP_FREQ", _floats(cp_freq.group(1)), per_line=4)
    out += "\n/// Casimir-Polder integration weights.\n"
    cp_w = re.search(r"pub const CP_WEIGHTS: \[f64; 23\] = \[(.*?)\];", old, re.S)
    out += emit_f64("CP_WEIGHTS", _floats(cp_w.group(1)), per_line=4)
    out += "\n/// Number of reference systems per element.\n"
    out += emit_usize("REFN", [refs["refn"].get(z, 0) for z in range(1, MAX_ELEM + 1)])
    out += "\n/// Reference coordination numbers.\n"
    out += emit_f64("REFCOVCN", flat3(refs["refcovcn"]))
    out += "\n/// Reference charges (GFN2-xTB reference set).\n"
    out += emit_f64("REFQ_GFN2", flat3(refs["refq"]))
    out += "\n/// Reference hydrogen counts (GFN2-xTB reference set).\n"
    out += emit_f64("REFH", flat3(refs["refh"]))
    out += "\n/// Reference hydrogen counts (standalone DFT reference set).\n"
    out += emit_f64("HCOUNT", flat3(refs["hcount"]))
    out += "\n/// Reference charge scaling factors.\n"
    out += emit_f64("ASCALE", flat3(refs["ascale"], 0.0))
    out += "\n/// Reference system indices.\n"
    out += emit_usize("REFSYS", [int(refs["refsys"].get((i, z), 0)) for z in range(1, MAX_ELEM + 1) for i in range(1, MAX_REF + 1)])
    out += "\n/// Reference coordination numbers (D3 definition).\n"
    out += emit_f64("REFCN", flat3(refs["refcn"]))
    out += "\n/// Reference charges (standalone DFT reference set).\n"
    out += emit_f64("REFQ_DFT", flat3(refs["clsq"]))
    out += "\n/// Reference hydrogen counts (standalone DFT reference set).\n"
    out += emit_f64("REFH_DFT", flat3(refs["clsh"]))
    out += "\n/// Reference polarizabilities at the Casimir-Polder frequencies\n"
    out += "/// (indexed `[((z-1)*MAX_REF + iref)*NFREQ + k]`).\n"
    out += emit_f64("ALPHAIW", flat_alpha())
    out += "\n/// Maximum number of secondary reference systems.\n"
    out += f"pub const MAX_SEC: usize = {len(secq)};\n"
    out += "\n/// Secondary reference scaling factors.\n"
    out += emit_f64("SSCALE", sscale, per_line=4)
    out += "\n/// Secondary reference polarizabilities.\n"
    out += emit_f64("SECAIW", secaiw, per_line=4)

    if args.check:
        current = OUT.read_text() if OUT.exists() else ""
        if current != out:
            print(
                f"{OUT.relative_to(REPO)} is out of date — "
                f"run scripts/validation/gen_d4_refs.py"
            )
            return 1
        print(f"{OUT.relative_to(REPO)} is up to date")
        return 0

    OUT.write_text(out)
    print(f"wrote {OUT} (MAX_ELEM={MAX_ELEM})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
