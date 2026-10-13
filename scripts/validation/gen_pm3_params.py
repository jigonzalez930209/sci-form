#!/usr/bin/env python3
"""Generate the Rust PM3 parameter table from MOPAC's authoritative source.

Inputs (fetched from the openmopac/mopac repository):
  - src/models/parameters_for_PM3_C.F90  (Stewart's PM3 element data)
  - src/models/parameters_C.F90          (atomic heats of formation, `eheat`)

Output: src/pm3/params_table.rs with one `Pm3Params` literal per element.

The transition-metal entries (PM3(tm), Cundari et al.) are preserved from the
existing `src/pm3/params.rs` table and merged in.
"""

import re
import sys
from pathlib import Path

MOPAC_SRC = Path("/tmp/opencode/mopac-full/src/models")
OUT = Path(__file__).resolve().parents[2] / "src" / "pm3" / "params_table.rs"


def parse_pm3(path: Path) -> dict:
    text = path.read_text()
    # Strip Fortran continuation (& at line end) and split on `data`.
    text = re.sub(r"&\s*\n\s*", "", text)
    arrays: dict[str, dict[int, float]] = {}
    for m in re.finditer(r"data\s+(\w+)\s*\(\s*(\d+)\s*(?:,\s*(\d+)\s*)?\)\s*/\s*([^/]+)/", text):
        name, i, j, value = m.group(1), int(m.group(2)), m.group(3), m.group(4)
        value = value.strip().rstrip("d0D0").replace("d", "e").replace("D", "e")
        try:
            val = float(value)
        except ValueError:
            continue
        arrays.setdefault(name, {})
        if j is None:
            arrays[name][i] = val
        else:
            # 2D: guesp1(i,j) etc.
            arrays[name][(i, int(j))] = val
    return arrays


def parse_int_array(path: Path, name: str) -> list:
    src = path.read_text()
    i = src.index(f"data {name}")
    j = src.index("/", i)
    k = src.index("/", j + 1)
    body = src[j + 1 : k]
    body = re.sub(r"!.*", "", body).replace("&", " ").replace("\n", " ")
    vals = []
    for tok in body.split(","):
        tok = tok.strip()
        if not tok:
            continue
        m = re.match(r"^(\d+)\*(.+)$", tok)
        if m:
            vals.extend([m.group(2).strip()] * int(m.group(1)))
        else:
            vals.append(tok)
    out = []
    for v in vals:
        try:
            out.append(int(v))
        except ValueError:
            out.append(None)
    return out


def parse_eheat(path: Path) -> dict[int, float]:
    text = path.read_text()
    out = {}
    for m in re.finditer(r"data\s+eheat\s*\(\s*(\d+)\s*\)\s*/\s*([-\d.dD+]+)\s*/", text):
        out[int(m.group(1))] = float(m.group(2).replace("d", "e").replace("D", "e"))
    return out


def main() -> None:
    arrays = parse_pm3(MOPAC_SRC / "parameters_for_PM3_C.F90")
    eheat = parse_eheat(MOPAC_SRC / "parameters_C.F90")

    ios = parse_int_array(MOPAC_SRC / "parameters_C.F90", "ios")
    iop = parse_int_array(MOPAC_SRC / "parameters_C.F90", "iop")

    def eisol(z: int) -> float:
        """MOPAC `calpar`: total energy of the isolated atom (eV)."""
        n_s = ios[z - 1] or 0
        n_p = iop[z - 1] or 0
        gssc = max(n_s - 1, 0)
        gspc = n_s * n_p
        l = min(n_p, 6 - n_p)
        gp2c = (n_p * (n_p - 1)) / 2 + 0.5 * (l * (l - 1)) / 2
        gppc = -0.5 * (l * (l - 1)) / 2
        hspc = -n_p * n_s * 0.5
        return (
            uss[z] * n_s
            + upp.get(z, 0.0) * n_p
            + gss[z] * gssc
            + gpp.get(z, 0.0) * gppc
            + gsp.get(z, 0.0) * gspc
            + gp2.get(z, 0.0) * gp2c
            + hsp.get(z, 0.0) * hspc
        )

    uss = arrays["usspm3"]
    upp = arrays["upppm3"]
    betas = arrays["betasp"]
    betap = arrays["betapp"]
    zs = arrays["zspm3"]
    zp = arrays["zppm3"]
    alp = arrays["alppm3"]
    gss = arrays["gsspm3"]
    gsp = arrays["gsppm3"]
    gpp = arrays["gpppm3"]
    gp2 = arrays["gp2pm3"]
    hsp = arrays["hsppm3"]
    g1 = arrays["guesp1"]
    g2 = arrays["guesp2"]
    g3 = arrays["guesp3"]

    symbols = {
        1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O",
        9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",
        16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 30: "Zn", 31: "Ga",
        32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr", 37: "Rb", 38: "Sr",
        48: "Cd", 49: "In", 50: "Sn", 51: "Sb", 52: "Te", 53: "I", 54: "Xe",
        55: "Cs", 56: "Ba", 80: "Hg", 81: "Tl", 82: "Pb", 83: "Bi", 87: "Fr",
    }

    # Valence core charges (MOPAC `tore`): main-group valence electron counts.
    tore = {
        1: 1, 2: 2, 3: 1, 4: 2, 5: 3, 6: 4, 7: 5, 8: 6, 9: 7, 10: 8,
        11: 1, 12: 2, 13: 3, 14: 4, 15: 5, 16: 6, 17: 7, 18: 8,
        19: 1, 20: 2, 30: 2, 31: 3, 32: 4, 33: 5, 34: 6, 35: 7, 36: 8,
        37: 1, 38: 2, 48: 2, 49: 3, 50: 4, 51: 5, 52: 6, 53: 7, 54: 8,
        55: 1, 56: 2, 80: 2, 81: 3, 82: 4, 83: 5, 87: 1,
    }
    # nspqn: principal quantum number of the valence sp shell (MOPAC calpar).
    def nspqn(z: int) -> int:
        if z <= 2:
            return 1
        if z <= 10:
            return 2
        if z <= 18:
            return 3
        if z <= 36:
            return 4
        if z <= 54:
            return 5
        if z <= 86:
            return 6
        return 0

    entries = []
    for z in sorted(uss):
        if z not in symbols:
            continue
        gauss = []
        for k in (1, 2, 3, 4):
            if (z, k) in g1:
                gauss.append((g1[(z, k)], g2[(z, k)], g3[(z, k)]))
        gauss_s = ", ".join(f"({a:.7f}, {b:.7f}, {c:.7f})" for a, b, c in gauss)
        entries.append(
            "    Pm3Params {\n"
            f"        z: {z},\n"
            f"        symbol: \"{symbols[z]}\",\n"
            f"        tm: false,\n"
            f"        nspqn: {nspqn(z)},\n"
            f"        uss: {uss[z]:.7f},\n"
            f"        upp: {upp.get(z, 0.0):.7f},\n"
            f"        beta_s: {betas[z]:.7f},\n"
            f"        beta_p: {betap.get(z, 0.0):.7f},\n"
            f"        zeta_s: {zs[z]:.7f},\n"
            f"        zeta_p: {zp.get(z, 0.0):.7f},\n"
            f"        gss: {gss[z]:.7f},\n"
            f"        gsp: {gsp.get(z, 0.0):.7f},\n"
            f"        gpp: {gpp.get(z, 0.0):.7f},\n"
            f"        gp2: {gp2.get(z, 0.0):.7f},\n"
            f"        hsp: {hsp.get(z, 0.0):.7f},\n"
            f"        core_charge: {tore[z]:.1f},\n"
            f"        heat_of_atomization: {eheat.get(z, 0.0):.3f},\n"
            f"        eisol: {eisol(z):.7f},\n"
            f"        alpha: {alp[z]:.7f},\n"
            f"        gaussians: &[{gauss_s}],\n"
            "    },\n"
        )

    header = (
        "//! PM3 parameter table generated from MOPAC 23.2.5\n"
        "//! (`src/models/parameters_for_PM3_C.F90` and the `eheat` atomic heats of\n"
        "//! formation in `src/models/parameters_C.F90`).\n"
        "//!\n"
        "//! Do not edit by hand: regenerate with\n"
        "//! `scripts/validation/gen_pm3_params.py`.\n"
        "\n"
        "use super::params::Pm3Params;\n"
        "\n"
        "/// Main-group PM3 entries (Stewart 1989) with MOPAC's Gaussian core-core\n"
        "/// corrections and atomic reference heats of formation.\n"
        "pub(crate) static PM3_MAIN_GROUP: &[Pm3Params] = &[\n"
    )
    OUT.write_text(header + "".join(entries) + "];\n")
    print(f"wrote {OUT} with {len(entries)} main-group entries")


if __name__ == "__main__":
    main()
