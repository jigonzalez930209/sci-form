#!/usr/bin/env python3
"""Generate `src/charges_eeq/eeqbc_data.rs` from the multicharge sources.

Inputs (grimme-lab/multicharge @ main):
  - `src/multicharge/param/eeqbc2025.f90`: `eeqbc_chi`, `eeqbc_eta`,
    `eeqbc_rad`, `eeqbc_kcnchi`, `eeqbc_kqchi`, `eeqbc_kqeta`, `eeqbc_cap`,
    `eeqbc_cov_radii` (declared as `0.5_wp * [...]`), `eeqbc_avg_cn`
    for Z = 1..103.

Inputs (grimme-lab/mctc-lib @ main, used by `new_eeqbc2025_model`):
  - `src/mctc/data/vdwrad.f90`: packed pairwise van-der-Waals radii (raw table
    values, see the note in `main`), lower triangle with
    `idx(num1, num2) = min + max*(max-1)/2` for `min <= max` (1-based).
  - `src/mctc/data/paulingen.f90`: Pauling electronegativities (118 elements),
    normalized to fluorine by dividing by 3.98 inside the model constructor.

Model constants (`multicharge/param.f90::new_eeqbc2025_model`):
    kcnrad = 0.14, kbc = 0.60, cutoff = 25.0, cn_exp = 2.0, norm_exp = 0.75

Usage:
    python scripts/validation/gen_eeqbc_params.py \
        [--multicharge-dir DIR] [--mctc-dir DIR]
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "src" / "charges_eeq" / "eeqbc_data.rs"
MAX_ELEM = 103  # multicharge's eeqbc2025 coverage

def _floats(text: str) -> list[float]:
    return [float(x) for x in re.findall(r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?", text)]


def parse_array(path: Path, name: str, expected: int, scale: float = 1.0) -> list[float]:
    text = path.read_text()
    m = re.search(rf"{name}\s*\([^=]*\)\s*=\s*(?:[A-Za-z0-9._]+\s*\*\s*)?\[\s*&(.*?)\]", text, re.S)
    if not m:
        raise ValueError(f"{name} not found in {path}")
    vals = [v * scale for v in _floats(m.group(1))]
    if len(vals) != expected:
        raise ValueError(f"{name}: parsed {len(vals)} values, expected {expected}")
    return vals


def rust_array(name: str, values: list[float], doc: str) -> str:
    lines = [f"/// {doc}", f"pub const {name}: [f64; {len(values)}] = ["]
    for i in range(0, len(values), 4):
        # repr(): shortest decimal that round-trips to the same f64, so the
        # generated table keeps the full precision of the Fortran literals.
        chunk = ", ".join(repr(v) for v in values[i : i + 4])
        lines.append(f"    {chunk},")
    lines.append("];")
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--multicharge-dir", type=Path, default=Path("/tmp/opencode/multicharge-src"))
    ap.add_argument("--mctc-dir", type=Path, default=Path("/tmp/opencode/mctc-lib-src"))
    args = ap.parse_args()

    param = args.multicharge_dir / "src/multicharge/param/eeqbc2025.f90"
    chi = parse_array(param, "eeqbc_chi", MAX_ELEM)
    eta = parse_array(param, "eeqbc_eta", MAX_ELEM)
    rad = parse_array(param, "eeqbc_rad", MAX_ELEM)
    kcnchi = parse_array(param, "eeqbc_kcnchi", MAX_ELEM)
    kqchi = parse_array(param, "eeqbc_kqchi", MAX_ELEM)
    kqeta = parse_array(param, "eeqbc_kqeta", MAX_ELEM)
    cap = parse_array(param, "eeqbc_cap", MAX_ELEM)
    cov = parse_array(param, "eeqbc_cov_radii", MAX_ELEM, scale=0.5)
    avg_cn = parse_array(param, "eeqbc_avg_cn", MAX_ELEM)

    n_pairs = MAX_ELEM * (MAX_ELEM + 1) // 2
    # NOTE: mctc-lib declares the table as `aatoau * [...]` (Angstrom -> Bohr)
    # but `new_eeqbc2025_model` multiplies it by `autoaa` again before handing it
    # to the model, so the two factors cancel and the *raw* table value (an
    # Angstrom number) is what enters the bond-capacitor switch. Reproduced
    # verbatim: `scale=1.0`, see the probe in the EEQ-BC validation notes.
    vdw = parse_array(args.mctc_dir / "src/mctc/data/vdwrad.f90", "vdwrad", n_pairs)
    en = parse_array(args.mctc_dir / "src/mctc/data/paulingen.f90", "pauling_en", 118)

    parts = [
        "//! EEQ-BC 2025 parameters (generated; do not edit by hand).",
        "//!",
        "//! Source: grimme-lab/multicharge `src/multicharge/param/eeqbc2025.f90`",
        "//! (Froitzheim, Mueller, Hansen, Grimme, J. Chem. Phys. 162, 214109 (2025),",
        "//! DOI 10.1063/5.0268978), plus mctc-lib `vdwrad.f90` / `paulingen.f90`",
        "//! for the bond-capacitor reference distances and EN weights.",
        "//!",
        "//! Regenerate with `scripts/validation/gen_eeqbc_params.py`.",
        "//! Units: Bohr for lengths, Hartree for chi/eta.",
        "",
        "#[rustfmt::skip]",
        rust_array("EEQBC_CHI", chi, "Element electronegativities (Hartree), Z = 1..103."),
        "",
        "#[rustfmt::skip]",
        rust_array("EEQBC_ETA", eta, "Element chemical hardnesses (Hartree), Z = 1..103."),
        "",
        "#[rustfmt::skip]",
        rust_array("EEQBC_RAD", rad, "Gaussian charge widths (Bohr), Z = 1..103."),
        "",
        "#[rustfmt::skip]",
        rust_array("EEQBC_KCNCHI", kcnchi, "CN scaling of electronegativity, Z = 1..103."),
        "",
        "#[rustfmt::skip]",
        rust_array("EEQBC_KQCHI", kqchi, "Local-charge scaling of electronegativity, Z = 1..103."),
        "",
        "#[rustfmt::skip]",
        rust_array("EEQBC_KQETA", kqeta, "Local-charge scaling of hardness, Z = 1..103."),
        "",
        "#[rustfmt::skip]",
        rust_array("EEQBC_CAP", cap, "Bond capacitances, Z = 1..103."),
        "",
        "#[rustfmt::skip]",
        rust_array("EEQBC_COVRAD", cov, "Covalent radii for the CN counting function (Bohr)."),
        "",
        "#[rustfmt::skip]",
        rust_array("EEQBC_AVG_CN", avg_cn, "Average coordination numbers, Z = 1..103."),
        "",
        "#[rustfmt::skip]",
        rust_array(
            "EEQBC_VDW_PAIRS",
            vdw,
            "Packed pairwise van-der-Waals distances as passed by\n/// `new_eeqbc2025_model` (raw mctc-lib table value); index\n/// `min + max*(max-1)/2` (1-based) for `min <= max`.",
        ),
        "",
        "#[rustfmt::skip]",
        rust_array("EEQBC_PAULING_EN", en, "Pauling electronegativities (F = 3.98), Z = 1..118."),
        "",
    ]
    OUT.write_text("\n".join(parts))
    print(f"wrote {OUT.relative_to(REPO)} ({len(chi)} elements, {n_pairs} vdW pairs, {len(en)} EN)")


if __name__ == "__main__":
    main()
