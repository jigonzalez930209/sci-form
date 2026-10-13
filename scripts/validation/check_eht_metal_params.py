#!/usr/bin/env python3
"""Benchmark the transition-metal EHT parameters against the YAeHMOP table.

`src/eht/params.rs` cites Alvarez, *Tables of Parameters for Extended Hückel
Calculations* (1993) for the 3d/4d/5d metals. The YAeHMOP distribution
(SourceForge, `yaehmop/tightbind/eht_parms.dat`) compiles the same literature
values into a machine-readable table, so it is used as the reference here.

The script parses both tables and reports, per element/orbital:
  - VSIP (eV) and Slater exponent ζ, our value vs the reference;
  - flags rows where |Δ| exceeds the tolerances below.

Rows where the reference is incomplete (zero placeholders: Y, Hf 6p/5d) are
skipped and reported separately. Our table keeps the *first* exponent of the
reference's double-zeta d functions (single-zeta basis, documented).

Usage:
    python scripts/validation/check_eht_metal_params.py [--ref PATH] [--patch]
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
PARAMS = REPO / "src" / "eht" / "params.rs"
DEFAULT_REF = Path("/tmp/opencode/yaehmop/tightbind/eht_parms.dat")

TOL_VSIP = 0.02  # eV
TOL_ZETA = 0.02

# Elements whose reference rows are incomplete (zero placeholders): their
# entries are reported as skipped instead of mismatches.
INCOMPLETE_REFERENCE = {39, 47, 72}

TM_BLOCKS = {
    21: "SC", 22: "TI", 23: "V", 24: "CR", 25: "MN", 26: "FE", 27: "CO",
    28: "NI", 29: "CU", 30: "ZN", 39: "Y", 40: "ZR", 41: "NB", 42: "MO",
    43: "TC", 44: "RU", 45: "RH", 46: "PD", 47: "AG", 48: "CD", 72: "HF",
    73: "TA", 74: "W", 75: "RE", 76: "OS", 77: "IR", 78: "PT", 79: "AU",
    80: "HG",
}


def parse_reference(path: Path) -> dict[int, dict[str, tuple[float, float]]]:
    ref: dict[int, dict[str, tuple[float, float]]] = {}
    for line in path.read_text().splitlines():
        if line.startswith(";") or not line.strip():
            continue
        f = line.split()
        if len(f) < 8:
            continue
        try:
            z = int(f[1])
        except ValueError:
            continue
        if z in TM_BLOCKS:
            ref.setdefault(z, {})[f[5]] = (float(f[6]), float(f[7]))
    return ref


def parse_ours() -> dict[str, list[tuple[str, float, float]]]:
    src = PARAMS.read_text()
    ours: dict[str, list[tuple[str, float, float]]] = {}
    for m in re.finditer(r"static ([A-Z]+)_ORBITALS: \[OrbitalDef; (\d+)\] = \[(.*?)\n\];", src, re.S):
        orbitals = re.findall(
            r'label: "([^"]+)",\s*vsip: ([-\d.]+),\s*zeta: ([-\d.]+)', m.group(3)
        )
        ours[m.group(1)] = [(c, float(d), float(e)) for c, d, e in orbitals]
    return ours


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", type=Path, default=DEFAULT_REF)
    ap.add_argument("--patch", action="store_true", help="print corrected Rust blocks")
    ap.add_argument(
        "--emit",
        type=Path,
        default=None,
        help="write the reference table as a validation fixture",
    )
    args = ap.parse_args()

    ref = parse_reference(args.ref)
    ours = parse_ours()

    n_ok = n_diff = n_skip = 0
    fixes: dict[str, dict[str, tuple[float, float]]] = {}
    for z, sym in sorted(TM_BLOCKS.items()):
        if sym not in ours or z not in ref:
            print(f"{sym:3s} Z={z}: sin bloque propio o sin referencia")
            continue
        if z in INCOMPLETE_REFERENCE:
            print(f"{sym:3s} Z={z}: referencia incompleta en YAeHMOP, se mantiene nuestra tabla")
            n_skip += len(ours[sym])
            continue
        for label, vsip, zeta in ours[sym]:
            ang = label[-1]
            if ang not in ref[z]:
                continue
            r_vsip, r_zeta = ref[z][ang]
            if r_vsip == 0.0 and r_zeta == 0.0:
                print(f"{sym:3s} {label:4s}: referencia incompleta (placeholder 0), se mantiene")
                n_skip += 1
                continue
            d_v, d_z = vsip - r_vsip, zeta - r_zeta
            ok = abs(d_v) <= TOL_VSIP and abs(d_z) <= TOL_ZETA
            if ok:
                n_ok += 1
                continue
            n_diff += 1
            print(
                f"{sym:3s} {label:4s}  ours vsip={vsip:9.3f} zeta={zeta:7.4f} | "
                f"ref vsip={r_vsip:9.3f} zeta={r_zeta:7.4f} | "
                f"Δvsip={d_v:+7.3f} Δzeta={d_z:+7.4f}"
            )
            fixes.setdefault(sym, {})[label] = (r_vsip, r_zeta)

    print(f"\ncoinciden: {n_ok}; discrepancias: {n_diff}; saltadas (ref incompleta): {n_skip}")

    if args.emit is not None:
        import json

        elements = []
        for z, sym in sorted(TM_BLOCKS.items()):
            if z not in ref:
                continue
            rows = []
            for label, vsip, zeta in ours.get(sym, []):
                ang = label[-1]
                if ang not in ref[z]:
                    continue
                r_vsip, r_zeta = ref[z][ang]
                rows.append(
                    {
                        "orbital": label,
                        "reference_vsip_ev": r_vsip,
                        "reference_zeta": r_zeta,
                        "reference_complete": r_vsip != 0.0
                        and r_zeta != 0.0
                        and z not in INCOMPLETE_REFERENCE,
                    }
                )
            elements.append({"z": z, "symbol": sym, "orbitals": rows})
        fixture = {
            "source": "YAeHMOP 3.0.3 tightbind/eht_parms.dat (SourceForge), "
            "compilation of Alvarez, Tables of Parameters for Extended Huckel "
            "Calculations (1993) and the Hoffmann-group papers",
            "units": {"vsip": "eV", "zeta": "bohr^-1"},
            "tolerance": {"vsip_ev": TOL_VSIP, "zeta": TOL_ZETA},
            "note": "Our table keeps the first exponent of the reference's "
            "double-zeta d functions (single-zeta basis). Elements whose "
            "reference rows are incomplete keep our values.",
            "elements": elements,
        }
        args.emit.parent.mkdir(parents=True, exist_ok=True)
        args.emit.write_text(json.dumps(fixture, indent=1) + "\n")
        print(f"fixture escrito en {args.emit}")

    if args.patch and fixes:
        print("\n// valores corregidos (literatura YAeHMOP/Alvarez):")
        for sym, orbitals in sorted(fixes.items()):
            for label, (vsip, zeta) in sorted(orbitals.items()):
                print(f"//   {sym} {label}: vsip={vsip}, zeta={zeta}")


if __name__ == "__main__":
    main()
