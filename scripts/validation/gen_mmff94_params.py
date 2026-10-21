#!/usr/bin/env python3
"""Generate src/forcefield/mmff94_params.rs from the published MMFF94 tables.

Deterministic code generator (re-run after any upstream refresh):

    python scripts/validation/gen_mmff94_params.py

Upstream source
---------------
The full MMFF94 parameter set is embedded in the RDKit source tree as C++
string literals:

    https://raw.githubusercontent.com/rdkit/rdkit/master/Code/ForceField/MMFF/Params.cpp

RDKit is distributed under the BSD-3-Clause license.  The parameter data it
carries is Copyright (c) Merck and Co., Inc., 1994, 1995, 1996 and is the
published MMFF94 set of Halgren, T., *J. Comput. Chem.* **1996**, *17*,
490-519 (plus the MMFF94s static variant of Halgren's MMFF94s paper).

The file is downloaded once into ``.cache/mmff94/Params.cpp`` (override with
``--params-cpp PATH``) and the SHA-256 of the parsed input is recorded in the
generated Rust header so the provenance of every number is traceable.

What is emitted
---------------
``src/forcefield/mmff94_params.rs`` containing

* ``MMFF94_DEFS``       - atom-type equivalence levels (MMFF.I Table III)
* ``MMFF94_PROPS``      - atom-type properties (MMFF.II Table I)
* ``MMFF94_AROMATIC_TYPES`` - aromatic atom types
* ``MMFF94_BOND_PARAMS``   - bond stretching parameters (kb, r0)
* ``MMFF94_ANGLE_PARAMS``  - angle bending parameters (ka, theta0)
* ``MMFF94_TOR_PARAMS`` / ``MMFF94S_TOR_PARAMS`` - torsion parameters
* ``MMFF94_VDW_PARAMS``    - van der Waals atom-type parameters
* ``MMFF94_BNDK_PARAMS``   - default bond length / force constant rule
* ``MMFF94_HERSCHBACH_LAURIE`` - Badger's rule parameters
* ``MMFF94_COV_RAD_PAU_ELE``   - covalent radius / Pauling electronegativity
"""

from __future__ import annotations

import argparse
import hashlib
import re
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
CACHE_DIR = ROOT / ".cache" / "mmff94"
OUT = ROOT / "src" / "forcefield" / "mmff94_params.rs"

PARAMS_URL = (
    "https://raw.githubusercontent.com/rdkit/rdkit/master/Code/ForceField/MMFF/Params.cpp"
)

# Maximum MMFF94 atom type number (99 = MG+2).
MAX_TYPE = 99


# ─── upstream retrieval ──────────────────────────────────────────────────────


def fetch_params_cpp(cache_dir: Path, params_cpp: Path | None) -> tuple[str, str]:
    """Return (source text, sha256 hex) of the RDKit Params.cpp."""
    if params_cpp is not None:
        text = params_cpp.read_text()
        return text, hashlib.sha256(text.encode()).hexdigest()
    cache_dir.mkdir(parents=True, exist_ok=True)
    target = cache_dir / "Params.cpp"
    if not target.exists():
        print(f"downloading {PARAMS_URL} -> {target}", file=sys.stderr)
        with urllib.request.urlopen(PARAMS_URL) as fh:
            target.write_bytes(fh.read())
    text = target.read_text()
    return text, hashlib.sha256(text.encode()).hexdigest()


# ─── C++ string literal extraction ───────────────────────────────────────────

_ESCAPES = {
    "n": "\n",
    "t": "\t",
    "r": "\r",
    "\\": "\\",
    '"': '"',
    "'": "'",
    "0": "\0",
    "a": "\a",
    "b": "\b",
    "f": "\f",
    "v": "\v",
}


def _decode_literal(lit: str) -> str:
    out: list[str] = []
    i = 0
    while i < len(lit):
        ch = lit[i]
        if ch != "\\":
            out.append(ch)
            i += 1
            continue
        i += 1
        if i >= len(lit):
            break
        esc = lit[i]
        if esc in _ESCAPES:
            out.append(_ESCAPES[esc])
            i += 1
        elif esc == "x":
            j = i + 1
            hex_digits = ""
            while j < len(lit) and len(hex_digits) < 2 and lit[j] in "0123456789abcdefABCDEF":
                hex_digits += lit[j]
                j += 1
            out.append(chr(int(hex_digits, 16)))
            i = j
        else:
            out.append(esc)
            i += 1
    return "".join(out)


def extract_cpp_string(src: str, name: str) -> str:
    """Decode the C++ string literal(s) bound to ``name`` into one text blob.

    Handles both ``const std::string name = "..." "...";`` and the array form
    ``const std::string name[] = { "...", "...", };``.
    """
    m = re.search(r"const\s+std::string\s+" + re.escape(name) + r"\s*(\[\s*\])?\s*=\s*", src)
    if not m:
        raise KeyError(f"{name}: declaration not found")
    is_array = bool(m.group(1))
    i = m.end()
    if is_array:
        if src[i] != "{":
            raise ValueError(f"{name}: expected '{{'")
        i += 1
        terminator = "}"
    else:
        terminator = ";"

    chunks: list[str] = []
    while i < len(src):
        ch = src[i]
        if ch == '"':
            j = i + 1
            while j < len(src):
                if src[j] == "\\":
                    j += 2
                    continue
                if src[j] == '"':
                    break
                j += 1
            chunks.append(_decode_literal(src[i + 1 : j]))
            i = j + 1
            continue
        if ch == terminator:
            break
        i += 1
    return "".join(chunks)


def data_lines(blob: str) -> list[list[str]]:
    """Split a decoded parameter blob into tab-separated data records."""
    rows: list[list[str]] = []
    for line in blob.split("\n"):
        if not line or line.startswith("*"):
            continue
        fields = line.split("\t")
        if not fields or not fields[0].strip():
            continue
        rows.append([f.strip() for f in fields])
    return rows


# ─── table parsers ───────────────────────────────────────────────────────────


def parse_defs(blob: str) -> dict[int, tuple[int, int, int, int]]:
    """MMFFDEF.PAR: symbol, type, four equivalence levels (first row per type wins)."""
    out: dict[int, tuple[int, int, int, int]] = {}
    for row in data_lines(blob):
        if len(row) < 6:
            continue
        atype = int(row[1])
        if atype in out:
            continue
        out[atype] = (int(row[2]), int(row[3]), int(row[4]), int(row[5]))
    return out


def parse_props(blob: str) -> dict[int, tuple[int, ...]]:
    """MMFFPROP.PAR: atype, aspec, crd, val, pilp, mltb, arom, lin, sbmb."""
    out: dict[int, tuple[int, ...]] = {}
    for row in data_lines(blob):
        if len(row) < 9:
            continue
        out[int(row[0])] = tuple(int(v) for v in row[1:9])
    return out


def parse_bonds(blob: str) -> list[tuple[int, int, int, float, float]]:
    """MMFFBND.PAR: bond type, i type, j type, kb (md/A), r0 (A)."""
    out = []
    for row in data_lines(blob):
        if len(row) < 5:
            continue
        out.append((int(row[0]), int(row[1]), int(row[2]), float(row[3]), float(row[4])))
    return out


def parse_bndk(blob: str) -> list[tuple[int, int, float, float]]:
    """MMFFBNDK.PAR: atomic number i, atomic number j, r0 (A), kb (md/A)."""
    out = []
    for row in data_lines(blob):
        if len(row) < 4:
            continue
        out.append((int(row[0]), int(row[1]), float(row[2]), float(row[3])))
    return out


def parse_herschbach_laurie(blob: str) -> list[tuple[int, int, float, float, float]]:
    """Badger's rule table: periodic row i, periodic row j, a_ij, d_ij, dp_ij."""
    out = []
    for row in data_lines(blob):
        if len(row) < 5:
            continue
        out.append((int(row[0]), int(row[1]), float(row[2]), float(row[3]), float(row[4])))
    return out


def parse_cov_rad_pau_ele(blob: str) -> list[tuple[int, float, float]]:
    """Covalent radius (A) and Pauling electronegativity per atomic number."""
    out = []
    for row in data_lines(blob):
        if len(row) < 3:
            continue
        out.append((int(row[0]), float(row[1]), float(row[2])))
    return out


def parse_angles(blob: str) -> list[tuple[int, int, int, int, float, float]]:
    """MMFFANG.PAR: angle type, i, j, k types, ka (md.A/rad^2), theta0 (deg)."""
    out = []
    for row in data_lines(blob):
        if len(row) < 6:
            continue
        out.append(
            (int(row[0]), int(row[1]), int(row[2]), int(row[3]), float(row[4]), float(row[5]))
        )
    return out


def parse_torsions(blob: str) -> list[tuple[int, int, int, int, int, float, float, float]]:
    """MMFFTOR.PAR: torsion type, i, j, k, l types, V1, V2, V3 (kcal/mol)."""
    out = []
    for row in data_lines(blob):
        if len(row) < 8:
            continue
        out.append(
            (
                int(row[0]),
                int(row[1]),
                int(row[2]),
                int(row[3]),
                int(row[4]),
                float(row[5]),
                float(row[6]),
                float(row[7]),
            )
        )
    return out


def parse_vdw(blob: str) -> tuple[tuple[float, float, float, float, float], list]:
    """MMFFVDW.PAR: globals (power, B, Beta, DARAD, DAEPS) plus per-type rows."""
    rows = data_lines(blob)
    globals_row = (float(rows[0][0]), float(rows[0][1]), float(rows[0][2]),
                   float(rows[0][3]), float(rows[0][4]))
    types = []
    for row in rows[1:]:
        if len(row) < 6:
            continue
        da = {"-": 0, "A": 1, "D": 2}[row[5]]
        types.append((int(row[0]), float(row[1]), float(row[2]), float(row[3]),
                      float(row[4]), da))
    return globals_row, types


def parse_aromatic_types(src: str) -> list[int]:
    m = re.search(r"defaultMMFFArom\s*=\s*\{(.*?)\}", src, re.S)
    if not m:
        raise KeyError("defaultMMFFArom not found")
    return [int(v) for v in re.findall(r"\d+", m.group(1))]


# ─── Rust emission helpers ───────────────────────────────────────────────────


def rust_f64(x: float) -> str:
    """Shortest round-trip decimal literal for a float."""
    text = repr(float(x))
    if text in ("inf", "-inf", "nan"):
        raise ValueError(f"non-finite value {x}")
    if "." not in text and "e" not in text:
        text += ".0"
    return text


def rust_rows(name: str, ty: str, rows: list[str]) -> list[str]:
    # `#[rustfmt::skip]` keeps the generated tables in their compact
    # one-row-per-parameter form (the same pattern as `src/xtb/d4_data.rs`).
    out = ["#[rustfmt::skip]", f"pub const {name}: &[{ty}] = &["]
    for row in rows:
        out.append(f"    {row},")
    out.append("];")
    return out


# ─── main ─────────────────────────────────────────────────────────────────────


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--params-cpp", type=Path, default=None,
                    help="local copy of RDKit Code/ForceField/MMFF/Params.cpp")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    src, sha = fetch_params_cpp(CACHE_DIR, args.params_cpp)

    defs = parse_defs(extract_cpp_string(src, "defaultMMFFDef"))
    props = parse_props(extract_cpp_string(src, "defaultMMFFProp"))
    arom = parse_aromatic_types(src)
    bonds = parse_bonds(extract_cpp_string(src, "defaultMMFFBond"))
    bndk = parse_bndk(extract_cpp_string(src, "defaultMMFFBndk"))
    hl = parse_herschbach_laurie(extract_cpp_string(src, "defaultMMFFHerschbachLaurie"))
    cov = parse_cov_rad_pau_ele(extract_cpp_string(src, "defaultMMFFCovRadPauEle"))
    angles = parse_angles(extract_cpp_string(src, "defaultMMFFAngleData"))
    tors = parse_torsions(extract_cpp_string(src, "defaultMMFFTor"))
    tors_s = parse_torsions(extract_cpp_string(src, "defaultMMFFsTor"))
    vdw_globals, vdw_types = parse_vdw(extract_cpp_string(src, "defaultMMFFVdW"))

    out: list[str] = []
    out.append("//! MMFF94 parameter tables (Halgren, *J. Comput. Chem.* 1996, 17, 490-519).")
    out.append("//!")
    out.append("//! Generated by `scripts/validation/gen_mmff94_params.py` - do not edit by hand.")
    out.append("//!")
    out.append("//! Value source: the published MMFF94 parameter set embedded in the RDKit")
    out.append("//! source tree at `Code/ForceField/MMFF/Params.cpp`")
    out.append("//! (<https://github.com/rdkit/rdkit>), which is distributed under the")
    out.append("//! BSD-3-Clause license.  The parameter data itself is")
    out.append("//! Copyright (c) Merck and Co., Inc., 1994, 1995, 1996.")
    out.append("//!")
    out.append(f"//! Upstream SHA-256: `{sha}`")
    out.append("//!")
    out.append("//! Units (as documented by RDKit's `Params.h` and Halgren 1996):")
    out.append("//! bond `kb` in md/A and `r0` in A; angle `ka` in md.A/rad^2 and")
    out.append("//! `theta0` in degrees; torsion `V1/V2/V3` in kcal/mol; van der Waals")
    out.append("//! `alpha_i` in A^3, `N_i` dimensionless, `A_i` in A^4, `G_i` in")
    out.append("//! kcal.A^4/mol.")
    out.append("")

    # ── atom type equivalence levels ──
    out.append("/// Atom-type equivalence levels (MMFF.I Table III, `MMFFDEF.PAR`).")
    out.append("///")
    out.append("/// `eq_level[0]` is the type itself (level 2 in Halgren's numbering);")
    out.append("/// `eq_level[1..3]` are the successively more general \"default types\" used")
    out.append("/// to resolve angle and torsion parameters that are not tabulated for the")
    out.append("/// full type quadruple.")
    out.append("#[derive(Debug, Clone, Copy)]")
    out.append("pub struct Mmff94Def {")
    out.append("    pub eq_level: [u8; 4],")
    out.append("}")
    out.append("")
    out.append("/// Equivalence levels indexed by MMFF94 atom type (index 0 is unused).")
    out.append("#[rustfmt::skip]")
    out.append(f"pub const MMFF94_DEFS: [Mmff94Def; {MAX_TYPE + 1}] = [")
    for t in range(MAX_TYPE + 1):
        eq = defs.get(t)
        if eq is None:
            out.append("    Mmff94Def { eq_level: [0, 0, 0, 0] },")
        else:
            out.append("    Mmff94Def { eq_level: [%d, %d, %d, %d] }," % eq)
    out.append("];")
    out.append("")

    # ── atom type properties ──
    out.append("/// Atom-type properties (MMFF.II Table I, `MMFFPROP.PAR`).")
    out.append("///")
    out.append("/// * `atno`  - atomic number")
    out.append("/// * `crd`  - number of connections (coordination number)")
    out.append("/// * `val`  - number of valence electrons")
    out.append("/// * `pilp` - 1 if the atom has a p-typed lone pair")
    out.append("/// * `mltb` - number of lone pairs on the atom")
    out.append("/// * `arom` - 1 if the atom is aromatic")
    out.append("/// * `linh` - 1 if the atom is linear (sp hybridised with 2 connections)")
    out.append("/// * `sbmb` - 1 if single bonds to this atom are \"special\" (conjugated)")
    out.append("#[derive(Debug, Clone, Copy)]")
    out.append("pub struct Mmff94Prop {")
    out.append("    pub atno: u8,")
    out.append("    pub crd: u8,")
    out.append("    pub val: u8,")
    out.append("    pub pilp: u8,")
    out.append("    pub mltb: u8,")
    out.append("    pub arom: u8,")
    out.append("    pub linh: u8,")
    out.append("    pub sbmb: u8,")
    out.append("}")
    out.append("")
    out.append("/// Atom-type properties indexed by MMFF94 atom type (index 0 is unused).")
    out.append("#[rustfmt::skip]")
    out.append(f"pub const MMFF94_PROPS: [Mmff94Prop; {MAX_TYPE + 1}] = [")
    for t in range(MAX_TYPE + 1):
        p = props.get(t)
        if p is None:
            out.append("    Mmff94Prop { atno: 0, crd: 0, val: 0, pilp: 0, mltb: 0, arom: 0, linh: 0, sbmb: 0 },")
        else:
            out.append(
                "    Mmff94Prop { atno: %d, crd: %d, val: %d, pilp: %d, mltb: %d, arom: %d, linh: %d, sbmb: %d },"
                % p
            )
    out.append("];")
    out.append("")

    out.append("/// Aromatic MMFF94 atom types (`MMFFArom` collection).")
    out.append("pub const MMFF94_AROMATIC_TYPES: &[u8] = &[")
    out.append("    " + ", ".join(str(v) for v in arom) + ",")
    out.append("];")
    out.append("")

    # ── bond stretching ──
    out.append("/// Bond stretching parameters (MMFF.II Table II, `MMFFBND.PAR`).")
    out.append("///")
    out.append("/// `kb` is the force constant in md/A, `r0` the reference length in A.")
    out.append("/// Rows are sorted by `(i_type, j_type, bond_type)` with `i_type <= j_type`")
    out.append("/// so that they can be binary searched.")
    out.append("#[derive(Debug, Clone, Copy)]")
    out.append("pub struct Mmff94BondParam {")
    out.append("    pub bond_type: u8,")
    out.append("    pub i_type: u8,")
    out.append("    pub j_type: u8,")
    out.append("    pub kb: f64,")
    out.append("    pub r0: f64,")
    out.append("}")
    out.append("")
    bonds_sorted = sorted(bonds, key=lambda r: (r[1], r[2], r[0]))
    out.extend(rust_rows(
        "MMFF94_BOND_PARAMS",
        "Mmff94BondParam",
        ["Mmff94BondParam { bond_type: %d, i_type: %d, j_type: %d, kb: %s, r0: %s }"
         % (bt, i, j, rust_f64(kb), rust_f64(r0))
         for (bt, i, j, kb, r0) in bonds_sorted],
    ))
    out.append("")

    # ── angle bending ──
    out.append("/// Angle bending parameters (MMFF.III, `MMFFANG.PAR`).")
    out.append("///")
    out.append("/// `ka` is the force constant in md.A/rad^2 and `theta0` the reference")
    out.append("/// angle in degrees.  Rows are sorted by `(angle_type, j_type, i_type,")
    out.append("/// k_type)` so that they can be binary searched.")
    out.append("#[derive(Debug, Clone, Copy)]")
    out.append("pub struct Mmff94AngleParam {")
    out.append("    pub angle_type: u8,")
    out.append("    pub i_type: u8,")
    out.append("    pub j_type: u8,")
    out.append("    pub k_type: u8,")
    out.append("    pub ka: f64,")
    out.append("    pub theta0: f64,")
    out.append("}")
    out.append("")
    angles_sorted = sorted(angles, key=lambda r: (r[0], r[2], r[1], r[3]))
    out.extend(rust_rows(
        "MMFF94_ANGLE_PARAMS",
        "Mmff94AngleParam",
        ["Mmff94AngleParam { angle_type: %d, i_type: %d, j_type: %d, k_type: %d, ka: %s, theta0: %s }"
         % (at, i, j, k, rust_f64(ka), rust_f64(th))
         for (at, i, j, k, ka, th) in angles_sorted],
    ))
    out.append("")

    # ── torsions ──
    out.append("/// Torsion parameters (MMFF.IV, `MMFFTOR.PAR`).")
    out.append("///")
    out.append("/// `v1`, `v2`, `v3` are in kcal/mol.  Rows are sorted by")
    out.append("/// `(tor_type, i_type, j_type, k_type, l_type)` so that they can be binary")
    out.append("/// searched.")
    out.append("#[derive(Debug, Clone, Copy)]")
    out.append("pub struct Mmff94TorParam {")
    out.append("    pub tor_type: u8,")
    out.append("    pub i_type: u8,")
    out.append("    pub j_type: u8,")
    out.append("    pub k_type: u8,")
    out.append("    pub l_type: u8,")
    out.append("    pub v1: f64,")
    out.append("    pub v2: f64,")
    out.append("    pub v3: f64,")
    out.append("}")
    out.append("")

    def emit_torsion_table(const_name: str, rows) -> None:
        rows_sorted = sorted(rows, key=lambda r: (r[0], r[1], r[2], r[3], r[4]))
        out.extend(rust_rows(
            const_name,
            "Mmff94TorParam",
            ["Mmff94TorParam { tor_type: %d, i_type: %d, j_type: %d, k_type: %d, l_type: %d, v1: %s, v2: %s, v3: %s }"
             % (tt, i, j, k, l, rust_f64(v1), rust_f64(v2), rust_f64(v3))
             for (tt, i, j, k, l, v1, v2, v3) in rows_sorted],
        ))
        out.append("")

    out.append("/// MMFF94 (dynamics) torsion parameters.")
    emit_torsion_table("MMFF94_TOR_PARAMS", tors)
    out.append("/// MMFF94s (static) torsion parameters.")
    emit_torsion_table("MMFF94S_TOR_PARAMS", tors_s)

    # ── van der Waals ──
    out.append("/// Global van der Waals parameters (MMFF.V Table I, `MMFFVDW.PAR`).")
    out.append("#[derive(Debug, Clone, Copy)]")
    out.append("pub struct Mmff94VdWGlobals {")
    out.append("    /// Exponent of `alpha_i` in `R*_i = A_i * alpha_i^power`.")
    out.append("    pub power: f64,")
    out.append("    /// Buffering coefficient `B` of the pair minimum rule.")
    out.append("    pub b: f64,")
    out.append("    /// Exponential damping coefficient `Beta`.")
    out.append("    pub beta: f64,")
    out.append("    /// Radius scaling applied to donor/acceptor pairs.")
    out.append("    pub darad: f64,")
    out.append("    /// Well depth scaling applied to donor/acceptor pairs.")
    out.append("    pub daeps: f64,")
    out.append("}")
    out.append("")
    g = vdw_globals
    out.append("pub const MMFF94_VDW_GLOBALS: Mmff94VdWGlobals = Mmff94VdWGlobals {")
    out.append("    power: %s," % rust_f64(g[0]))
    out.append("    b: %s," % rust_f64(g[1]))
    out.append("    beta: %s," % rust_f64(g[2]))
    out.append("    darad: %s," % rust_f64(g[3]))
    out.append("    daeps: %s," % rust_f64(g[4]))
    out.append("};")
    out.append("")
    out.append("/// Per-atom-type van der Waals parameters (`MMFFVDW.PAR`).")
    out.append("///")
    out.append("/// `da` encodes the MMFF donor/acceptor flag: 0 = none, 1 = acceptor,")
    out.append("/// 2 = donor.")
    out.append("#[derive(Debug, Clone, Copy)]")
    out.append("pub struct Mmff94VdWParam {")
    out.append("    pub atom_type: u8,")
    out.append("    pub alpha_i: f64,")
    out.append("    pub n_i: f64,")
    out.append("    pub a_i: f64,")
    out.append("    pub g_i: f64,")
    out.append("    pub da: u8,")
    out.append("}")
    out.append("")
    out.extend(rust_rows(
        "MMFF94_VDW_PARAMS",
        "Mmff94VdWParam",
        ["Mmff94VdWParam { atom_type: %d, alpha_i: %s, n_i: %s, a_i: %s, g_i: %s, da: %d }"
         % (t, rust_f64(al), rust_f64(n), rust_f64(a), rust_f64(gi), da)
         for (t, al, n, a, gi, da) in vdw_types],
    ))
    out.append("")

    # ── default bond rule ──
    out.append("/// Default bond length / force constant rule (`MMFFBNDK.PAR`).")
    out.append("///")
    out.append("/// `r0` in A, `kb` in md/A; keyed by the pair of atomic numbers.")
    out.append("#[derive(Debug, Clone, Copy)]")
    out.append("pub struct Mmff94BndkParam {")
    out.append("    pub z_i: u8,")
    out.append("    pub z_j: u8,")
    out.append("    pub r0: f64,")
    out.append("    pub kb: f64,")
    out.append("}")
    out.append("")
    bndk_sorted = sorted(bndk, key=lambda r: (min(r[0], r[1]), max(r[0], r[1])))
    out.extend(rust_rows(
        "MMFF94_BNDK_PARAMS",
        "Mmff94BndkParam",
        ["Mmff94BndkParam { z_i: %d, z_j: %d, r0: %s, kb: %s }"
         % (zi, zj, rust_f64(r0), rust_f64(kb))
         for (zi, zj, r0, kb) in bndk_sorted],
    ))
    out.append("")

    # ── Herschbach-Laurie (Badger's rule) ──
    out.append("/// Herschbach-Laurie parameters for Badger's rule (MMFF.V, Table I).")
    out.append("///")
    out.append("/// Indexed by the pair of periodic-table rows (Halgren's `HL` rows).")
    out.append("#[derive(Debug, Clone, Copy)]")
    out.append("pub struct Mmff94HerschbachLaurie {")
    out.append("    pub row_i: u8,")
    out.append("    pub row_j: u8,")
    out.append("    pub a_ij: f64,")
    out.append("    pub d_ij: f64,")
    out.append("    pub dp_ij: f64,")
    out.append("}")
    out.append("")
    hl_sorted = sorted(hl, key=lambda r: (min(r[0], r[1]), max(r[0], r[1])))
    out.extend(rust_rows(
        "MMFF94_HERSCHBACH_LAURIE",
        "Mmff94HerschbachLaurie",
        ["Mmff94HerschbachLaurie { row_i: %d, row_j: %d, a_ij: %s, d_ij: %s, dp_ij: %s }"
         % (ri, rj, rust_f64(a), rust_f64(d), rust_f64(dp))
         for (ri, rj, a, d, dp) in hl_sorted],
    ))
    out.append("")

    # ── covalent radius / Pauling electronegativity ──
    out.append("/// Covalent radius (A) and Pauling electronegativity per atomic number.")
    out.append("#[derive(Debug, Clone, Copy)]")
    out.append("pub struct Mmff94CovRadPauEle {")
    out.append("    pub z: u8,")
    out.append("    pub r0: f64,")
    out.append("    pub chi: f64,")
    out.append("}")
    out.append("")
    out.extend(rust_rows(
        "MMFF94_COV_RAD_PAU_ELE",
        "Mmff94CovRadPauEle",
        ["Mmff94CovRadPauEle { z: %d, r0: %s, chi: %s }" % (z, rust_f64(r), rust_f64(c))
         for (z, r, c) in cov],
    ))
    out.append("")

    args.out.write_text("\n".join(out))
    print(
        "wrote %s (%d bytes): defs=%d props=%d arom=%d bonds=%d angles=%d "
        "tor=%d tor_s=%d vdw=%d bndk=%d hl=%d covrad=%d"
        % (
            args.out,
            args.out.stat().st_size,
            len(defs),
            len(props),
            len(arom),
            len(bonds),
            len(angles),
            len(tors),
            len(tors_s),
            len(vdw_types),
            len(bndk),
            len(hl),
            len(cov),
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
