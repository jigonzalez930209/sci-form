#!/usr/bin/env python3
"""Independent numpy oracle for sci-form's Extended Hückel Theory (ficha F86-EHT).

The oracle implements Extended Hückel Theory from scratch in numpy:

  * STO-nG expansion of the Slater valence orbitals (HSP 1969 / Stewart 1970
    contraction tables, transcribed below),
  * the analytic overlap matrix `S` of the contracted Gaussian basis
    (Hermite/moment expansion, an independent route to the standard
    Cartesian Gaussian overlap integral),
  * the Wolfsberg-Helmholtz Hamiltonian
    `H_ij = (K/2) S_ij (H_ii + H_jj)` with `K = 1.75` and `H_ii = VSIP`,
  * the generalized eigenproblem `H C = S C E` (Löwdin orthogonalisation is
    avoided by using scipy's symmetric generalized solver),
  * the valence-electron count and the HOMO/LUMO index and gap.

The *published* tables that the oracle needs are handled in two different
ways on purpose:

  * the STO-nG contraction coefficients are **transcribed into this file**
    (published values, Hehre/Stewart/Pople JCP 51, 2657 (1969) and Stewart
    JCP 54, 70 (1970)), and the oracle *verifies* that sci-form's
    `src/eht/basis.rs` still uses the same tables; any divergence is written
    into the fixture as `sto3g_table_matches_rust` so the Rust gate can
    report it,
  * the EHT *parameters* (VSIP and Slater exponents) and the hydrogenic
    metal ionisation potentials are **read from the Rust sources**
    (`src/eht/params.rs`, `src/eht/metal_fallback.rs`), because the parameter
    tables themselves are already gated by F31 against the YAeHMOP
    reference.  This oracle gates the *energy*, not the table.

Provenance of the parameters: sci-form's `src/eht/params.rs`, itself
transcribed from Hoffmann, R. J. Chem. Phys. 39, 1397 (1963) and the
Alvarez (1993) compilation of Extended-Hückel parameters.

Usage:
    python scripts/validation/oracle_eht.py
"""

from __future__ import annotations

import os

# Pin the BLAS/OpenMP thread count *before* numpy is imported: the
# generalized eigenproblem otherwise sums in a thread-dependent order,
# which changes the last digits of the reference values and makes the
# oracle non-reproducible (the fixtures must be byte-identical on every run).
for _thread_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_thread_var] = "1"

import argparse
import json
import math
import re
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from common import Report, emit_fixture_report

REPO = Path(__file__).resolve().parents[2]

# ---------------------------------------------------------------------------
# Constants / conventions
# ---------------------------------------------------------------------------

#: Wolfsberg-Helmholtz constant (sci-form: `eht::params::DEFAULT_K`).
WOLFSBERG_HELMHOLTZ_K = 1.75

#: Angstrom -> Bohr conversion used by sci-form's `eht::basis::build_basis`.
ANG_TO_BOHR = 1.0 / 0.529177249

#: Hartree in eV (only used for the derived-metal exponents).
HARTREE_EV = 13.605693

# ---------------------------------------------------------------------------
# STO-nG contraction tables (published values, transcribed here on purpose)
# ---------------------------------------------------------------------------
#
# Hehre, Stewart & Pople, J. Chem. Phys. 51, 2657 (1969)  -> 1s, 2sp, 3sp
# Stewart, J. Chem. Phys. 54, 70 (1970)                    -> 4sp, 4d, 5d
#
# Each entry is `(contraction coefficient, primitive exponent)` for a
# Slater orbital with zeta = 1; the exponents are scaled by zeta^2.

STO3G_1S = (
    (0.154329, 2.227660),
    (0.535328, 0.405771),
    (0.444635, 0.109818),
)
STO3G_2SP_INNER = (
    (-0.099967, 2.581058),
    (0.399513, 0.532013),
    (0.700115, 0.168856),
)
STO3G_2SP_OUTER = (
    (0.155916, 2.581058),
    (0.607684, 0.532013),
    (0.391957, 0.168856),
)
STO3G_3SP_INNER = (
    (-0.219620, 0.994203),
    (0.225595, 0.231031),
    (0.900398, 0.075142),
)
STO3G_3SP_OUTER = (
    (0.010588, 0.994203),
    (0.595167, 0.231031),
    (0.462001, 0.075142),
)
STO3G_4SP_INNER = (
    (-0.310438, 0.494718),
    (0.015515, 0.120193),
    (1.023468, 0.040301),
)
STO3G_4SP_OUTER = (
    (-0.063311, 0.494718),
    (0.574459, 0.120193),
    (0.499768, 0.040301),
)
STO3G_5SP_INNER = (
    (-0.383204, 0.289515),
    (-0.159438, 0.073780),
    (1.143456, 0.025335),
)
STO3G_5SP_OUTER = (
    (-0.094810, 0.289515),
    (0.570520, 0.073780),
    (0.497340, 0.025335),
)
STO3G_3D = (
    (0.219767, 2.488296),
    (0.655547, 0.798105),
    (0.286573, 0.331966),
)
STO3G_4D = (
    (0.2308552718, 0.1777717219),
    (0.6042409177, 0.0804064735),
    (0.2595768926, 0.03949855551),
)
STO3G_5D = (
    (-0.02010175008, 0.4913352950),
    (0.5899370608, 0.07329090601),
    (0.4658445960, 0.03594209290),
)

#: (n, l) -> contraction table (same lookup as sci-form's `sto3g_expansion`).
STO3G_TABLES = {
    (1, 0): STO3G_1S,
    (2, 0): STO3G_2SP_INNER,
    (2, 1): STO3G_2SP_OUTER,
    (3, 0): STO3G_3SP_INNER,
    (3, 1): STO3G_3SP_OUTER,
    (3, 2): STO3G_3D,
    (4, 0): STO3G_4SP_INNER,
    (4, 1): STO3G_4SP_OUTER,
    (4, 2): STO3G_4D,
    (5, 0): STO3G_5SP_INNER,
    (5, 1): STO3G_5SP_OUTER,
    (5, 2): STO3G_5D,
    (6, 0): STO3G_5SP_INNER,
    (6, 1): STO3G_5SP_OUTER,
    (6, 2): STO3G_5D,
    (7, 0): STO3G_5SP_INNER,
    (7, 1): STO3G_5SP_OUTER,
    (7, 2): STO3G_5D,
}

#: Rust source of the contraction tables (used for the divergence check).
RUST_BASIS_RS = REPO / "src" / "eht" / "basis.rs"

# ---------------------------------------------------------------------------
# Molecular geometries (Angstrom)
# ---------------------------------------------------------------------------
#
# Every geometry is built from standard structural parameters so that the
# oracle stays deterministic and free of external dependencies.


def _tetrahedral(center, ref_dir, bond, h_angles_deg, azimuth_deg):
    """Three substituents around `center` at `h_angles_deg` from `ref_dir`."""
    ref = np.asarray(ref_dir, dtype=float)
    ref = ref / np.linalg.norm(ref)
    # build a local frame
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(tmp, ref))) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    e1 = np.cross(ref, tmp)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(ref, e1)
    out = []
    for az in azimuth_deg:
        a = math.radians(h_angles_deg)
        phi = math.radians(az)
        direction = math.cos(a) * ref + math.sin(a) * (
            math.cos(phi) * e1 + math.sin(phi) * e2
        )
        out.append(center + bond * direction)
    return out


def _ferrocene():
    """Staggered (D5d) ferrocene: Fe-Cp centroid 1.65 A, ring radius 1.20 A."""
    atoms = []
    ring_r = 1.20
    ring_z = 1.65
    ch = 1.06
    for ring, (z, phase) in enumerate(((ring_z, 0.0), (-ring_z, math.pi / 5.0))):
        for k in range(5):
            ang = phase + 2.0 * math.pi * k / 5.0
            c = np.array([ring_r * math.cos(ang), ring_r * math.sin(ang), z])
            atoms.append((6, c))
            h = c + ch * np.array([math.cos(ang), math.sin(ang), 0.0])
            atoms.append((1, h))
    atoms.append((26, np.zeros(3)))
    return atoms


def _ethanol():
    """Ethanol from standard internal coordinates (r / angles in the usual
    experimental values: rCC 1.5126, rCO 1.4312, rOH 0.9709)."""
    c1 = np.zeros(3)
    c2 = np.array([1.5126, 0.0, 0.0])
    atoms = [(6, c1), (6, c2)]

    # CH3 hydrogens: 110.4 deg from the C-C axis, tetrahedral azimuths.
    for p in _tetrahedral(c1, c2 - c1, 1.0962, 110.4, (0.0, 120.0, 240.0)):
        atoms.append((1, p))

    # O: 107.99 deg C-C-O angle.
    ang_cco = math.radians(107.99)
    o = c2 + 1.4312 * np.array(
        [math.cos(math.pi - ang_cco), math.sin(math.pi - ang_cco), 0.0]
    )
    atoms.append((8, o))

    # CH2 hydrogens: symmetric about the C1-C2-O plane.
    a1 = (c1 - c2) / np.linalg.norm(c1 - c2)
    a2 = (o - c2) / np.linalg.norm(o - c2)
    m = a1 + a2
    m /= np.linalg.norm(m)
    b = np.cross(a1, a2)
    b /= np.linalg.norm(b)
    cos_alpha = math.cos(math.radians(109.5)) / math.cos(math.radians(54.0))
    alpha = math.acos(cos_alpha)
    for sign in (+1.0, -1.0):
        direction = math.cos(alpha) * m + sign * math.sin(alpha) * b
        atoms.append((1, c2 + 1.0934 * direction))

    # O-H: 108.0 deg C-O-H, in the C1-C2-O plane.
    o_to_c2 = (c2 - o) / np.linalg.norm(c2 - o)
    rot = np.array(
        [
            [
                math.cos(math.radians(108.0)),
                -math.sin(math.radians(108.0)),
                0.0,
            ],
            [
                math.sin(math.radians(108.0)),
                math.cos(math.radians(108.0)),
                0.0,
            ],
            [0.0, 0.0, 1.0],
        ]
    )
    atoms.append((1, o + 0.9709 * (rot @ o_to_c2)))
    return atoms


def _benzene():
    """Planar benzene: ring radius 1.3971 A, C-H 1.0878 A radially out."""
    atoms = []
    ring_r = 1.3971
    ch = 1.0878
    for k in range(6):
        ang = 2.0 * math.pi * k / 6.0
        c = np.array([ring_r * math.cos(ang), ring_r * math.sin(ang), 0.0])
        atoms.append((6, c))
        atoms.append((1, c + ch * np.array([math.cos(ang), math.sin(ang), 0.0])))
    return atoms


def _acetylene():
    return [
        (6, np.array([0.0, 0.0, 0.0])),
        (6, np.array([1.2033, 0.0, 0.0])),
        (1, np.array([-1.0605, 0.0, 0.0])),
        (1, np.array([2.2638, 0.0, 0.0])),
    ]


def _ammonia():
    r = 1.0124
    theta = math.radians(106.67)
    atoms = [(7, np.zeros(3))]
    for k in range(3):
        phi = 2.0 * math.pi * k / 3.0
        atoms.append(
            (
                1,
                r * np.array(
                    [
                        math.sin(theta) * math.cos(phi),
                        math.sin(theta) * math.sin(phi),
                        math.cos(theta),
                    ]
                ),
            )
        )
    return atoms


def _methane():
    r = 1.0870
    atoms = [(6, np.zeros(3))]
    # tetrahedral CH4: the four H point away from the C along the four
    # tetrahedral directions (the "ref" direction is the H opposite to +x).
    for p in _tetrahedral(
        np.zeros(3), np.array([1.0, 1.0, 1.0]), r, 109.471, (0.0, 120.0, 240.0)
    ):
        atoms.append((1, p))
    atoms.append(
        (
            1,
            -r * np.array([1.0, 1.0, 1.0]) / math.sqrt(3.0),
        )
    )
    return atoms


def _water():
    r = 0.9584
    theta = math.radians(104.45)
    atoms = [(8, np.zeros(3))]
    for sign in (+1.0, -1.0):
        atoms.append(
            (
                1,
                r * np.array([math.sin(theta / 2.0) * sign, 0.0, math.cos(theta / 2.0)]),
            )
        )
    return atoms


def _hydrogen():
    return [(1, np.zeros(3)), (1, np.array([0.7414, 0.0, 0.0]))]


def _uo2():
    return [
        (8, np.array([0.0, 0.0, 1.7500])),
        (92, np.zeros(3)),
        (8, np.array([0.0, 0.0, -1.7500])),
    ]


MOLECULES = [
    ("H2", _hydrogen()),
    ("H2O", _water()),
    ("CH4", _methane()),
    ("NH3", _ammonia()),
    ("C2H2", _acetylene()),
    ("benzene", _benzene()),
    ("ethanol", _ethanol()),
    ("ferrocene", _ferrocene()),
    ("UO2", _uo2()),
]


# ---------------------------------------------------------------------------
# Reading the EHT parameters from the Rust sources
# ---------------------------------------------------------------------------


def _read(path: Path) -> str:
    return path.read_text()


def parse_params_rs() -> dict[int, dict]:
    """Parse `src/eht/params.rs` into {Z: {'symbol': str, 'orbitals': [...]}}."""
    text = _read(REPO / "src" / "eht" / "params.rs")
    blocks = re.findall(
        r"static\s+(\w+)_ORBITALS\s*:\s*\[OrbitalDef;\s*\d+\]\s*=\s*\[(.*?)\];",
        text,
        re.S,
    )
    by_name: dict[str, list[dict]] = {}
    for name, body in blocks:
        orbitals = []
        for m in re.finditer(
            r"n:\s*(\d+)\s*,\s*l:\s*(\d+)\s*,\s*label:\s*\"([^\"]+)\"\s*,\s*"
            r"vsip:\s*(-?[\d.]+)\s*,\s*zeta:\s*(-?[\d.]+)",
            body,
        ):
            orbitals.append(
                {
                    "n": int(m.group(1)),
                    "l": int(m.group(2)),
                    "label": m.group(3),
                    "vsip": float(m.group(4)),
                    "zeta": float(m.group(5)),
                }
            )
        by_name[name] = orbitals

    out: dict[int, dict] = {}
    for m in re.finditer(
        r"EhtParams\s*\{\s*z:\s*(\d+)\s*,\s*symbol:\s*\"([^\"]+)\"\s*,\s*"
        r"orbitals:\s*&(\w+)_ORBITALS\s*,",
        text,
    ):
        z = int(m.group(1))
        out[z] = {"symbol": m.group(2), "orbitals": by_name[m.group(3)]}
    if not out:
        raise RuntimeError("could not parse src/eht/params.rs")
    return out


def parse_first_ip_ev() -> dict[int, float]:
    """Parse the hydrogenic first ionisation potentials from metal_fallback.rs."""
    text = _read(REPO / "src" / "eht" / "metal_fallback.rs")
    body = re.search(r"fn first_ip_ev\(z: u8\) -> f64 \{(.*?)\n\}", text, re.S)
    if not body:
        raise RuntimeError("could not parse first_ip_ev()")
    out: dict[int, float] = {}
    for m in re.finditer(r"(\d+)\s*=>\s*(-?[\d.]+)\s*,", body.group(1)):
        out[int(m.group(1))] = float(m.group(2))
    return out


def parse_symbol_map() -> dict[int, str]:
    """Parse the element symbols of the derived (hydrogenic) parameters."""
    text = _read(REPO / "src" / "eht" / "metal_fallback.rs")
    body = re.search(r"fn symbol\(z: u8\) -> &'static str \{(.*?)\n\}", text, re.S)
    if not body:
        raise RuntimeError("could not parse symbol()")
    out: dict[int, str] = {}
    for m in re.finditer(r"(\d+)\s*=>\s*\"([A-Za-z]+)\"", body.group(1)):
        out[int(m.group(1))] = m.group(2)
    return out


def parse_derived_orbital_rules() -> dict:
    """Extract the *rules* that build the hydrogenic (derived) valence shells."""
    text = _read(REPO / "src" / "eht" / "metal_fallback.rs")
    rules: dict[str, str] = {}

    body = re.search(r"fn valence_n\(z: u8\) -> u8 \{(.*?)\n\}", text, re.S)
    rules["valence_n_body"] = body.group(1) if body else ""

    body = re.search(r"fn uses_d_valence\(z: u8\) -> bool \{(.*?)\n\}", text, re.S)
    rules["uses_d_body"] = body.group(1) if body else ""

    body = re.search(r"fn build_orbitals\(z: u8\) -> Vec<OrbitalDef> \{(.*?)\n\}", text, re.S)
    rules["build_orbitals_body"] = body.group(1) if body else ""

    body = re.search(r"fn valence_electrons\(z: u8\) -> usize \{(.*?)\n\}", text, re.S)
    rules["valence_electrons_body"] = body.group(1) if body else ""

    body = re.search(r"fn zeta_from_ip\(ip_ev: f64\) -> f64 \{(.*?)\n\}", text, re.S)
    rules["zeta_from_ip_body"] = body.group(1) if body else ""
    return rules


def valence_n_rust(z: int, rules: dict) -> int:
    """Re-implementation of `metal_fallback::valence_n`."""
    body = rules["valence_n_body"]
    if re.search(r"3\s*\|\s*4\s*=>", body) and (z == 3 or z == 4):
        return 2
    if re.search(r"11\.\.=13", body) and 11 <= z <= 13:
        return 3
    if re.search(r"19\s*\|\s*20\s*\|\s*31", body) and z in (19, 20, 31):
        return 4
    if re.search(r"37\s*\|\s*38\s*\|\s*49\s*\|\s*50", body) and z in (37, 38, 49, 50):
        return 5
    if re.search(r"55\.\.=83", body) and 55 <= z <= 83:
        return 6
    return 7


def uses_d_valence_rust(z: int, rules: dict) -> bool:
    body = rules["uses_d_body"]
    m = re.search(r"matches!\(z,\s*([0-9]+)\.\.=([0-9]+)\s*\|\s*([0-9]+)\.\.=([0-9]+)\s*\)", body)
    if not m:
        raise RuntimeError("could not parse uses_d_valence()")
    lo1, hi1, lo2, hi2 = (int(m.group(i)) for i in range(1, 5))
    return (lo1 <= z <= hi1) or (lo2 <= z <= hi2)


def with_p_rust(z: int, rules: dict) -> bool:
    """`with_p` in `build_orbitals` (transcribed from the Rust source)."""
    body = rules["build_orbitals_body"]
    m = re.search(
        r"let with_p = matches!\(z,\s*([0-9]+)\s*\|\s*([0-9]+)\s*\|\s*([0-9]+)\s*\|\s*"
        r"([0-9]+)\s*\|\s*([0-9]+)\s*\|\s*([0-9]+)\s*\|\s*([0-9]+)\)\s*\|\|\s*uses_d_valence\(z\);",
        body,
    )
    if not m:
        raise RuntimeError("could not parse with_p rule")
    listed = {int(m.group(i)) for i in range(1, 8)}
    return z in listed or uses_d_valence_rust(z, rules)


def zeta_from_ip_rust(ip_ev: float) -> float:
    """`metal_fallback::zeta_from_ip` (zeta = sqrt(IP / Hartree), clamped)."""
    return min(max(math.sqrt(max(ip_ev, 1.0) / HARTREE_EV), 0.35), 6.0)


def derived_params_rust(z: int, ips: dict[int, float], rules: dict) -> list[dict]:
    """Re-implementation of `metal_fallback::build_orbitals`."""
    n = valence_n_rust(z, rules)
    ip = ips.get(z, 6.0)
    orbitals = [
        {
            "n": n,
            "l": 0,
            "label": _shell_label_rust(n, 0, rules),
            "vsip": -ip,
            "zeta": zeta_from_ip_rust(ip),
        }
    ]
    if with_p_rust(z, rules):
        ip_p = max(ip - 3.5, 2.0)
        orbitals.append(
            {
                "n": n,
                "l": 1,
                "label": _shell_label_rust(n, 1, rules),
                "vsip": -ip_p,
                "zeta": zeta_from_ip_rust(ip_p),
            }
        )
    if uses_d_valence_rust(z, rules):
        ip_d = ip + 2.0
        orbitals.append(
            {
                "n": n - 1,
                "l": 2,
                "label": _shell_label_rust(n, 2, rules),
                "vsip": -ip_d,
                "zeta": zeta_from_ip_rust(ip_d),
            }
        )
    return orbitals


def _shell_label_rust(n: int, l: int, rules: dict) -> str:
    """`metal_fallback::shell_label`."""
    table = {
        (2, 0): "2s",
        (3, 0): "3s",
        (3, 1): "3p",
        (4, 0): "4s",
        (4, 1): "4p",
        (4, 2): "3d",
        (5, 0): "5s",
        (5, 1): "5p",
        (5, 2): "4d",
        (6, 0): "6s",
        (6, 1): "6p",
        (6, 2): "5d",
        (7, 0): "7s",
        (7, 1): "7p",
        (7, 2): "6d",
    }
    return table.get((n, l), "val")


def valence_electrons_rust(z: int, rules: dict) -> int:
    """`eht::solver::count_valence_electrons` (per element)."""
    body = rules["valence_electrons_body"]
    fixed = {
        1: 1,
        5: 3,
        6: 4,
        7: 5,
        8: 6,
        9: 7,
        14: 4,
        15: 5,
        16: 6,
        17: 7,
        35: 7,
        53: 7,
        21: 3,
        39: 3,
        22: 4,
        40: 4,
        72: 4,
        23: 5,
        41: 5,
        73: 5,
        24: 6,
        42: 6,
        74: 6,
        25: 7,
        43: 7,
        75: 7,
        26: 8,
        44: 8,
        76: 8,
        27: 9,
        45: 9,
        77: 9,
        28: 10,
        46: 10,
        78: 10,
        29: 11,
        47: 11,
        79: 11,
        30: 12,
        48: 12,
        80: 12,
    }
    if z in fixed:
        return fixed[z]
    # the fallback table lives in `metal_fallback::valence_electrons`
    fallback = {
        3: 1, 11: 1, 19: 1, 37: 1, 55: 1, 87: 1,
        4: 2, 12: 2, 20: 2, 38: 2, 56: 2, 88: 2,
        13: 3, 31: 3, 49: 3, 81: 3,
        50: 4, 82: 4,
        83: 5,
        63: 2, 70: 2,
        89: 3,
        90: 4,
        91: 5,
        92: 6,
        93: 7,
        94: 8,
        95: 9,
        96: 10,
        97: 11,
        98: 12,
        99: 13,
        100: 14,
        101: 15,
        102: 16,
        103: 3,
    }
    if 57 <= z <= 71:
        return 3
    return fallback.get(z, 0)


def count_valence_electrons_rust(elements: list[int], rules: dict) -> int:
    return sum(valence_electrons_rust(z, rules) for z in elements)


# ---------------------------------------------------------------------------
# Gaussian overlap machinery (independent implementation)
# ---------------------------------------------------------------------------


def odd_double_factorial(n: int) -> float:
    if n <= 0:
        return 1.0
    acc = 1.0
    k = n
    while k > 0:
        acc *= float(k)
        k -= 2
    return acc


def cartesian_norm(alpha: float, lx: int, ly: int, lz: int) -> float:
    """Normalisation of the Cartesian Gaussian x^lx y^ly z^lz exp(-alpha r^2)."""
    lsum = lx + ly + lz
    pref = (2.0 * alpha / math.pi) ** 0.75
    ang = (4.0 * alpha) ** (lsum / 2.0)
    denom = math.sqrt(
        odd_double_factorial(2 * lx - 1)
        * odd_double_factorial(2 * ly - 1)
        * odd_double_factorial(2 * lz - 1)
    )
    return pref * ang / denom


def gaussian_moment(m: int, p: float) -> float:
    """Integral of y^m exp(-p y^2) over the real line."""
    if m % 2:
        return 0.0
    k = m // 2
    val = math.sqrt(math.pi / p)
    for j in range(1, k + 1):
        val *= (2.0 * j - 1.0) / (2.0 * p)
    return val


def overlap_1d(i: int, j: int, a: float, b: float, alpha: float, beta: float) -> float:
    """Raw 1D overlap integral of (x-a)^i (x-b)^j exp(-alpha (x-a)^2 - beta (x-b)^2).

    The Gaussian product theorem is used:

        alpha (x-a)^2 + beta (x-b)^2 = p (x-P)^2 + mu (a-b)^2,
        p = alpha + beta,  P = (alpha a + beta b) / p,  mu = alpha beta / p,

    so the integral factorises into the "prefactor" exp(-mu (a-b)^2) and the
    moment integrals of the shifted Gaussian.
    """
    p = alpha + beta
    center = (alpha * a + beta * b) / p
    total = 0.0
    for t in range(i + 1):
        for u in range(j + 1):
            total += (
                math.comb(i, t)
                * math.comb(j, u)
                * (center - a) ** (i - t)
                * (center - b) ** (j - u)
                * gaussian_moment(t + u, p)
            )
    mu = alpha * beta / p
    total *= math.exp(-mu * (a - b) ** 2)
    return total


def cartesian_terms(l: int, m: int) -> list[tuple[float, tuple[int, int, int]]]:
    """Real-harmonic Cartesian terms; identical convention to sci-form."""
    if l == 0:
        return [(1.0, (0, 0, 0))]
    if l == 1:
        return {
            -1: [(1.0, (1, 0, 0))],
            0: [(1.0, (0, 1, 0))],
            1: [(1.0, (0, 0, 1))],
        }[m]
    if l == 2:
        s2 = math.sqrt(2.0)
        s6 = math.sqrt(6.0)
        return {
            -2: [(1.0, (1, 1, 0))],
            -1: [(1.0, (1, 0, 1))],
            0: [(1.0, (0, 1, 1))],
            1: [(1.0 / s2, (2, 0, 0)), (-1.0 / s2, (0, 2, 0))],
            2: [
                (-1.0 / s6, (2, 0, 0)),
                (-1.0 / s6, (0, 2, 0)),
                (2.0 / s6, (0, 0, 2)),
            ],
        }[m]
    raise ValueError(f"unsupported l = {l}")


class Ao:
    """One atomic orbital of the molecular basis (Gaussian expansion)."""

    __slots__ = ("atom", "center", "n", "l", "m", "vsip", "zeta", "label", "prims")

    def __init__(self, atom, center, n, l, m, vsip, zeta, label, prims):
        self.atom = atom
        self.center = center
        self.n = n
        self.l = l
        self.m = m
        self.vsip = vsip
        self.zeta = zeta
        self.label = label
        self.prims = prims


def sto3g_expansion(n: int, l: int, zeta: float) -> list[tuple[float, float]]:
    """STO-nG primitive expansion (coefficient, exponent) for a Slater orbital."""
    table = STO3G_TABLES[(n, l)]
    zeta_sq = zeta * zeta
    radial_scale = 1.0
    if n in (6, 7):
        ratio = 5.0 / float(n)
        radial_scale = ratio * ratio
    return [(c, a * zeta_sq * radial_scale) for (c, a) in table]


def build_basis(
    elements: list[int],
    positions_ang: np.ndarray,
    params: dict[int, dict],
    ips: dict[int, float],
    rules: dict,
    symbols: dict[int, str] | None = None,
) -> list[Ao]:
    """Build the molecular EHT basis (same AO ordering convention as sci-form)."""
    symbols = symbols or {}
    basis: list[Ao] = []
    for idx, (z, pos) in enumerate(zip(elements, positions_ang)):
        entry = params.get(z)
        if entry is not None:
            symbol = entry["symbol"]
            orbitals = entry["orbitals"]
        else:
            derived = derived_params_rust(z, ips, rules)
            if not derived:
                continue
            symbol = symbols.get(z, "M")
            orbitals = derived
        center = np.asarray(pos, dtype=float) * ANG_TO_BOHR
        for orb in orbitals:
            n, l = orb["n"], orb["l"]
            prims = sto3g_expansion(n, l, orb["zeta"])
            if l == 0:
                basis.append(
                    Ao(idx, center, n, 0, 0, orb["vsip"], orb["zeta"], f"{symbol}_{orb['label']}", prims)
                )
            elif l == 1:
                for m, tag in ((-1, "x"), (0, "y"), (1, "z")):
                    basis.append(
                        Ao(
                            idx,
                            center,
                            n,
                            1,
                            m,
                            orb["vsip"],
                            orb["zeta"],
                            f"{symbol}_{orb['label']}{tag}",
                            prims,
                        )
                    )
            elif l == 2:
                for m, tag in ((-2, "xy"), (-1, "xz"), (0, "yz"), (1, "x2-y2"), (2, "z2")):
                    basis.append(
                        Ao(
                            idx,
                            center,
                            n,
                            2,
                            m,
                            orb["vsip"],
                            orb["zeta"],
                            f"{symbol}_{orb['label']}{tag}",
                            prims,
                        )
                    )
    return basis


def ao_overlap(a: Ao, b: Ao) -> float:
    """Overlap of two contracted Gaussians (analytic, all terms)."""
    terms_a = cartesian_terms(a.l, a.m)
    terms_b = cartesian_terms(b.l, b.m)
    total = 0.0
    for (ca, alpha) in a.prims:
        for (cb, beta) in b.prims:
            prim = 0.0
            for coef_a, (lax, lay, laz) in terms_a:
                norm_a = cartesian_norm(alpha, lax, lay, laz)
                for coef_b, (lbx, lby, lbz) in terms_b:
                    norm_b = cartesian_norm(beta, lbx, lby, lbz)
                    ix = overlap_1d(lax, lbx, a.center[0], b.center[0], alpha, beta)
                    iy = overlap_1d(lay, lby, a.center[1], b.center[1], alpha, beta)
                    iz = overlap_1d(laz, lbz, a.center[2], b.center[2], alpha, beta)
                    prim += coef_a * coef_b * norm_a * norm_b * ix * iy * iz
            total += ca * cb * prim
    return total


def overlap_matrices(basis: list[Ao]) -> tuple[np.ndarray, np.ndarray]:
    """Return (S with the EHT convention S_ii = 1, S with the true diagonal)."""
    n = len(basis)
    s_true = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1):
            val = ao_overlap(basis[i], basis[j])
            s_true[i, j] = val
            s_true[j, i] = val
    s_norm = s_true.copy()
    np.fill_diagonal(s_norm, 1.0)
    return s_norm, s_true


def hamiltonian(basis: list[Ao], s: np.ndarray, k: float) -> np.ndarray:
    """Wolfsberg-Helmholtz Hamiltonian in eV."""
    n = len(basis)
    h = np.zeros((n, n))
    for i in range(n):
        h[i, i] = basis[i].vsip
    for i in range(n):
        for j in range(i + 1, n):
            hij = 0.5 * k * s[i, j] * (h[i, i] + h[j, j])
            h[i, j] = hij
            h[j, i] = hij
    return h


def solve_generalized(h: np.ndarray, s: np.ndarray) -> np.ndarray:
    """Eigenvalues of H C = S C E (ascending)."""
    w, _ = eigh(h, s)
    return np.asarray(w, dtype=float)


# ---------------------------------------------------------------------------
# STO-nG table cross-check against the Rust source
# ---------------------------------------------------------------------------

_RUST_TABLE_NAMES = {
    (1, 0): "STO3G_1S",
    (2, 0): "STO3G_2SP_INNER",
    (2, 1): "STO3G_2SP_OUTER",
    (3, 0): "STO3G_3SP_INNER",
    (3, 1): "STO3G_3SP_OUTER",
    (3, 2): "STO3G_3D",
    (4, 0): "STO3G_4SP_INNER",
    (4, 1): "STO3G_4SP_OUTER",
    (4, 2): "STO3G_4D",
    (5, 0): "STO3G_5SP_INNER",
    (5, 1): "STO3G_5SP_OUTER",
    (5, 2): "STO3G_5D",
}


def check_rust_sto3g_tables() -> dict:
    """Compare the transcribed tables with `src/eht/basis.rs`."""
    text = _read(RUST_BASIS_RS)
    out = {}
    for key, name in _RUST_TABLE_NAMES.items():
        m = re.search(
            rf"static\s+{name}\s*:\s*\[\(f64, f64\);\s*3\]\s*=\s*\[(.*?)\];", text, re.S
        )
        if not m:
            out[f"{key[0]}{key[1]}"] = {"found": False}
            continue
        rows = re.findall(r"\(\s*(-?[\d.]+)\s*,\s*(-?[\d.]+)\s*\)", m.group(1))
        rust = [(float(a), float(b)) for a, b in rows]
        mine = list(STO3G_TABLES[key])
        max_diff = max(
            max(abs(r[0] - o[0]), abs(r[1] - o[1])) for r, o in zip(rust, mine)
        )
        out[f"{key[0]}{key[1]}"] = {
            "found": True,
            "rust": [[c, a] for c, a in rust],
            "oracle": [[c, a] for c, a in mine],
            "max_abs_diff": max_diff,
            "identical": rust == mine,
        }
    return out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def _round_matrix(mat: np.ndarray, digits: int = 12) -> list[list[float]]:
    return [[float(round(v, digits)) for v in row] for row in mat]


def run() -> Report:
    params = parse_params_rs()
    ips = parse_first_ip_ev()
    rules = parse_derived_orbital_rules()
    symbols = parse_symbol_map()

    sto3g_check = check_rust_sto3g_tables()
    sto3g_identical = all(v.get("identical") for v in sto3g_check.values())

    cases = []
    for name, atoms in MOLECULES:
        elements = [z for z, _ in atoms]
        positions = np.array([p for _, p in atoms], dtype=float)
        basis = build_basis(elements, positions, params, ips, rules, symbols)
        s_norm, s_true = overlap_matrices(basis)
        h = hamiltonian(basis, s_norm, WOLFSBERG_HELMHOLTZ_K)

        energies = solve_generalized(h, s_norm)
        energies_true_diag = solve_generalized(h, s_true)

        n_electrons = count_valence_electrons_rust(elements, rules)
        n_orbitals = len(basis)
        n_occupied = -(-n_electrons // 2)  # ceil division
        homo_index = n_occupied - 1 if 0 < n_occupied <= n_orbitals else 0
        lumo_index = n_occupied if n_occupied < n_orbitals else homo_index
        homo = float(energies[homo_index])
        lumo = float(energies[lumo_index])

        diag_true = np.diag(s_true).copy()
        cases.append(
            {
                "name": name,
                "elements": elements,
                "coords": [round(float(v), 12) for v in positions.reshape(-1)],
                "n_basis": n_orbitals,
                "ao_labels": [ao.label for ao in basis],
                "n_electrons": n_electrons,
                "overlap": _round_matrix(s_norm),
                "overlap_true_diagonal": [float(round(v, 12)) for v in diag_true],
                "hamiltonian": _round_matrix(h),
                "energies_ev": [float(round(v, 12)) for v in energies],
                "energies_ev_true_diagonal": [
                    float(round(v, 12)) for v in energies_true_diag
                ],
                "homo_index": homo_index,
                "lumo_index": lumo_index,
                "homo_ev": float(round(homo, 12)),
                "lumo_ev": float(round(lumo, 12)),
                "gap_ev": float(round(lumo - homo, 12)),
            }
        )

    # --- derived (hydrogenic) metal parameters ---------------------------
    derived_cases = []
    for z in sorted(ips):
        orbitals = derived_params_rust(z, ips, rules)
        derived_cases.append(
            {
                "z": z,
                "orbitals": [
                    {
                        "n": o["n"],
                        "l": o["l"],
                        "label": o["label"],
                        "vsip": o["vsip"],
                        "zeta": o["zeta"],
                    }
                    for o in orbitals
                ],
            }
        )

    payload = {
        "ficha": "F86-EHT",
        "oracle": "oracle_eht.py (independent numpy Extended Hückel Theory)",
        "reference": (
            "independent numpy implementation of EHT (STO-nG overlap, "
            "Wolfsberg-Helmholtz Hamiltonian, generalized eigenproblem); "
            "parameters read from sci-form's src/eht/params.rs "
            "(transcribed from Hoffmann 1960 / Alvarez 1993)"
        ),
        "generated": time.strftime("%Y-%m-%d"),
        "conventions": {
            "k_wolfsberg_helmholtz": WOLFSBERG_HELMHOLTZ_K,
            "angstrom_to_bohr": ANG_TO_BOHR,
            "hartree_ev": HARTREE_EV,
            "overlap_diagonal": (
                "S_ii = 1 (Extended-Huckel convention: the valence orbitals are "
                "normalised); 'overlap_true_diagonal' reports the integral of the "
                "unnormalised STO-nG contraction for reference"
            ),
        },
        "sto3g_tables": sto3g_check,
        "sto3g_tables_identical_to_rust": sto3g_identical,
        "cases": cases,
        "derived_params": derived_cases,
    }

    report = Report(
        ficha="F86-EHT",
        observable="eht_overlap_max_abs",
        reference_tool="oracle_eht.py (independent numpy EHT)",
        units="dimensionless (S), eV (H, energies)",
        tolerance=1e-10,
        level="L2",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "Extended Huckel Theory: STO-nG overlap, Wolfsberg-Helmholtz "
            "Hamiltonian (K=1.75) and the generalized eigenproblem H C = S C E, "
            "compared against an independent numpy implementation. The STO-nG "
            "contraction tables (HSP 1969 / Stewart 1970) are transcribed in the "
            "oracle and cross-checked against sci-form's src/eht/basis.rs."
        ),
    )
    emit_fixture_report("eht_parity.json", payload, report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode",
        nargs="?",
        default="all",
        choices=["all"],
        help="which fixtures to generate (default: all)",
    )
    parser.parse_args()
    run()
    print("fixtures -> tests/fixtures/validation/eht_parity.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
