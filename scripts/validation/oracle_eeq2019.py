#!/usr/bin/env python3
"""Independent numpy oracle for sci-form's EEQ2019 charge model (ficha F87-EEQ2019).

The oracle re-implements the EEQ2019 model of Grimme et al. (JCTC 2019)
directly from the published equations, in numpy/scipy only:

  * fractional coordination numbers with an erf counting function and a
    smooth log cut-off,

        CN_i   = sum_{j != i} 1/2 [ 1 + erf( -kcn (r_ij - rc_ij) / rc_ij ) ]
        rc_ij  = (R_i + R_j) * 4/3                       (D3 radii)
        CN_i <- ln(1 + e^{CN_max}) - ln(1 + e^{CN_max - CN_i})

  * the Gaussian-width Coulomb matrix (atomic units),

        A_ij   = erf(r_ij gamma_ij) / r_ij ,  gamma_ij^2 = 1/(rad_i^2 + rad_j^2)
        A_ii   = eta_i + sqrt(2/pi) / rad_i

  * the CN-dependent electronegativity

        x_i    = -chi_i^0 + kcnchi_i * CN_i / sqrt(CN_i)

  * the electronegativity-equilibration system with the total-charge
    constraint

        [ A   1 ] [ q ]   [ x  ]
        [ 1^T 0 ] [ l ] = [ Q  ] .

The element parameters (chi, eta, rad, kcnchi and the Pyykko covalent radii)
are read from `src/charges_eeq/eeq2019_data.rs`, which is generated from
grimme-lab/multicharge (`src/multicharge/param/eeq2019.f90`) plus the dftd4
Pyykko radii by `scripts/validation/gen_eeq2019_params.py`.

Unlike sci-form's implementation the oracle uses the exact `scipy.special.erf`
(sci-form uses the Abramowitz & Stegun 7.1.26 approximation), so the two
sides differ by the accuracy of that rational approximation only.

Usage:
    python scripts/validation/oracle_eeq2019.py
"""

from __future__ import annotations

import os

# Pin the BLAS/OpenMP thread count *before* numpy is imported: the linear
# solver otherwise sums in a thread-dependent order, which changes the last
# digits of the reference values and makes the oracle non-reproducible.
for _thread_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_thread_var] = "1"

import argparse
import math
import re
import sys
import time
from pathlib import Path

import numpy as np
from scipy.special import erf

from common import Report, emit_fixture_report

REPO = Path(__file__).resolve().parents[2]
RUST_DATA = REPO / "src" / "charges_eeq" / "eeq2019_data.rs"

# ---------------------------------------------------------------------------
# Model constants (multicharge defaults)
# ---------------------------------------------------------------------------

#: Steepness of the erf coordination-number counting function.
CN_KCN = 7.5
#: Smooth CN cut-off (log cut).
CN_MAX = 8.0
#: Real-space CN cut-off (Bohr).
CN_CUTOFF = 25.0
#: Bohr per Angstrom (the value sci-form uses for the radius conversion).
BOHR_PER_ANGSTROM = 1.8897259886
#: D3 radii are 4/3 of the Pyykko (2009) set.
D3_RADIUS_SCALE = 4.0 / 3.0


# ---------------------------------------------------------------------------
# Parameter tables (read from the Rust source)
# ---------------------------------------------------------------------------


def parse_table(name: str) -> list[float]:
    text = RUST_DATA.read_text()
    m = re.search(
        rf"pub\(crate\)\s+const\s+{name}\s*:\s*\[f64;\s*(\d+)\]\s*=\s*\[(.*?)\];",
        text,
        re.S,
    )
    if not m:
        raise RuntimeError(f"could not parse {name} from {RUST_DATA}")
    values = [float(v) for v in re.findall(r"-?\d+\.\d+", m.group(2))]
    if len(values) != int(m.group(1)):
        raise RuntimeError(f"{name}: expected {m.group(1)} values, got {len(values)}")
    return values


def load_parameters() -> dict:
    return {
        "chi": parse_table("EEQ2019_CHI"),
        "eta": parse_table("EEQ2019_ETA"),
        "rad": parse_table("EEQ2019_RAD"),
        "kcnchi": parse_table("EEQ2019_KCNCHI"),
        "covrad": parse_table("EEQ2019_COVRAD_ANG"),
    }


# ---------------------------------------------------------------------------
# Molecular geometries (Angstrom)
# ---------------------------------------------------------------------------


def _tetrahedral(center, ref_dir, bond, h_angle_deg, azimuth_deg):
    ref = np.asarray(ref_dir, dtype=float)
    ref = ref / np.linalg.norm(ref)
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(tmp, ref))) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    e1 = np.cross(ref, tmp)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(ref, e1)
    out = []
    for az in azimuth_deg:
        a = math.radians(h_angle_deg)
        phi = math.radians(az)
        direction = math.cos(a) * ref + math.sin(a) * (
            math.cos(phi) * e1 + math.sin(phi) * e2
        )
        out.append(np.asarray(center, dtype=float) + bond * direction)
    return out


def _water():
    return [
        (8, np.array([0.0, 0.0, 0.0])),
        (1, np.array([0.7572, 0.5865, 0.0])),
        (1, np.array([-0.7572, 0.5865, 0.0])),
    ]


def _methane():
    atoms = [(6, np.zeros(3))]
    for p in _tetrahedral(
        np.zeros(3), np.array([1.0, 1.0, 1.0]), 1.0870, 109.471, (0.0, 120.0, 240.0)
    ):
        atoms.append((1, p))
    atoms.append((1, -1.0870 * np.array([1.0, 1.0, 1.0]) / math.sqrt(3.0)))
    return atoms


def _ammonia():
    r = 1.0124
    theta = math.radians(106.67)
    atoms = [(7, np.zeros(3))]
    for k in range(3):
        phi = 2.0 * math.pi * k / 3.0
        atoms.append(
            (
                1,
                r
                * np.array(
                    [
                        math.sin(theta) * math.cos(phi),
                        math.sin(theta) * math.sin(phi),
                        math.cos(theta),
                    ]
                ),
            )
        )
    return atoms


def _ammonium():
    atoms = [(7, np.zeros(3))]
    for p in _tetrahedral(
        np.zeros(3), np.array([1.0, 1.0, 1.0]), 1.0300, 109.471, (0.0, 120.0, 240.0)
    ):
        atoms.append((1, p))
    atoms.append((1, -1.0300 * np.array([1.0, 1.0, 1.0]) / math.sqrt(3.0)))
    return atoms


def _hydrogen_fluoride():
    return [(9, np.zeros(3)), (1, np.array([0.0, 0.0, 0.9170]))]


def _hydroxide():
    return [(8, np.zeros(3)), (1, np.array([0.0, 0.0, 0.9709]))]


def _hydronium():
    r = 1.0120
    theta = math.radians(113.6)
    atoms = [(8, np.zeros(3))]
    for k in range(3):
        phi = 2.0 * math.pi * k / 3.0
        atoms.append(
            (
                1,
                r
                * np.array(
                    [
                        math.sin(theta) * math.cos(phi),
                        math.sin(theta) * math.sin(phi),
                        math.cos(theta),
                    ]
                ),
            )
        )
    return atoms


def _ethane():
    """Staggered ethane: r(CC) = 1.5390, r(CH) = 1.0940."""
    c1 = np.zeros(3)
    c2 = np.array([1.5390, 0.0, 0.0])
    atoms = [(6, c1), (6, c2)]
    # the CH3 groups are staggered: the second group's azimuths are rotated
    # by 60 degrees with respect to the first one.
    for p in _tetrahedral(c1, c2 - c1, 1.0940, 110.2, (0.0, 120.0, 240.0)):
        atoms.append((1, p))
    for p in _tetrahedral(c2, c1 - c2, 1.0940, 110.2, (60.0, 180.0, 300.0)):
        atoms.append((1, p))
    return atoms


def _benzene():
    atoms = []
    ring_r = 1.3971
    ch = 1.0878
    for k in range(6):
        ang = 2.0 * math.pi * k / 6.0
        c = np.array([ring_r * math.cos(ang), ring_r * math.sin(ang), 0.0])
        atoms.append((6, c))
        atoms.append((1, c + ch * np.array([math.cos(ang), math.sin(ang), 0.0])))
    return atoms


def _ethanol():
    c1 = np.zeros(3)
    c2 = np.array([1.5126, 0.0, 0.0])
    atoms = [(6, c1), (6, c2)]
    for p in _tetrahedral(c1, c2 - c1, 1.0962, 110.4, (0.0, 120.0, 240.0)):
        atoms.append((1, p))
    ang_cco = math.radians(107.99)
    o = c2 + 1.4312 * np.array(
        [math.cos(math.pi - ang_cco), math.sin(math.pi - ang_cco), 0.0]
    )
    atoms.append((8, o))
    a1 = (c1 - c2) / np.linalg.norm(c1 - c2)
    a2 = (o - c2) / np.linalg.norm(o - c2)
    m = a1 + a2
    m /= np.linalg.norm(m)
    b = np.cross(a1, a2)
    b /= np.linalg.norm(b)
    alpha = math.acos(math.cos(math.radians(109.5)) / math.cos(math.radians(54.0)))
    for sign in (+1.0, -1.0):
        atoms.append((1, c2 + 1.0934 * (math.cos(alpha) * m + sign * math.sin(alpha) * b)))
    o_to_c2 = (c2 - o) / np.linalg.norm(c2 - o)
    rot = np.array(
        [
            [math.cos(math.radians(108.0)), -math.sin(math.radians(108.0)), 0.0],
            [math.sin(math.radians(108.0)), math.cos(math.radians(108.0)), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    atoms.append((1, o + 0.9709 * (rot @ o_to_c2)))
    return atoms


def _acetate():
    """Acetate anion CH3COO- (planar COO-, staggered CH3)."""
    c1 = np.array([0.0, 0.0, 0.0])
    c2 = np.array([1.2700, 0.0, 0.0])
    atoms = [(6, c1), (6, c2)]
    for p in _tetrahedral(c1, c2 - c1, 1.0900, 110.0, (0.0, 120.0, 240.0)):
        atoms.append((1, p))
    for sign in (+1.0, -1.0):
        o = c2 + 1.2500 * np.array([math.cos(math.radians(58.0)), sign * math.sin(math.radians(58.0)), 0.0])
        atoms.append((8, o))
    return atoms


def _sodium_chloride():
    return [(11, np.zeros(3)), (17, np.array([0.0, 0.0, 2.3607]))]


def _diiodine():
    return [(53, np.zeros(3)), (53, np.array([0.0, 0.0, 2.6660]))]


def _ferrohexaaqua():
    """[Fe(H2O)6]2+ (octahedral, Fe-O = 2.16 A)."""
    atoms = [(26, np.zeros(3))]
    dirs = [
        (1.0, 0.0, 0.0),
        (-1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, -1.0, 0.0),
        (0.0, 0.0, 1.0),
        (0.0, 0.0, -1.0),
    ]
    for d in dirs:
        d = np.asarray(d, dtype=float)
        o = 2.1600 * d
        atoms.append((8, o))
        # the two water hydrogens point away from the metal, in the plane
        # perpendicular to the Fe-O axis
        perp = np.cross(d, np.array([0.0, 0.0, 1.0]))
        if np.linalg.norm(perp) < 1e-8:
            perp = np.cross(d, np.array([0.0, 1.0, 0.0]))
        perp /= np.linalg.norm(perp)
        perp2 = np.cross(d, perp)
        for sign in (+1.0, -1.0):
            direction = 0.5372 * d + 0.8434 * sign * perp
            atoms.append((1, o + 0.9572 * direction))
    return atoms


MOLECULES = [
    ("H2O", 0.0, _water()),
    ("CH4", 0.0, _methane()),
    ("NH3", 0.0, _ammonia()),
    ("HF", 0.0, _hydrogen_fluoride()),
    ("C2H6", 0.0, _ethane()),
    ("benzene", 0.0, _benzene()),
    ("ethanol", 0.0, _ethanol()),
    ("NaCl", 0.0, _sodium_chloride()),
    ("I2", 0.0, _diiodine()),
    ("OH-", -1.0, _hydroxide()),
    ("NH4+", 1.0, _ammonium()),
    ("H3O+", 1.0, _hydronium()),
    ("acetate", -1.0, _acetate()),
    ("Fe(H2O)6_2+", 2.0, _ferrohexaaqua()),
]


# ---------------------------------------------------------------------------
# The EEQ2019 model
# ---------------------------------------------------------------------------


def log_cn_cut(cn: np.ndarray, cn_max: float) -> np.ndarray:
    """Smooth log cut-off of the coordination numbers."""
    return np.log1p(np.exp(cn_max)) - np.log1p(np.exp(cn_max - cn))


def coordination_numbers(
    elements: list[int], positions_bohr: np.ndarray, covrad: list[float]
) -> np.ndarray:
    """Fractional coordination numbers (erf counting function + log cut)."""
    n = len(elements)
    cn = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d = positions_bohr[i] - positions_bohr[j]
            r2 = float(d @ d)
            if r2 < 1e-12 or r2 > CN_CUTOFF * CN_CUTOFF:
                continue
            r = math.sqrt(r2)
            rc = (
                (covrad[elements[i] - 1] + covrad[elements[j] - 1])
                * BOHR_PER_ANGSTROM
                * D3_RADIUS_SCALE
            )
            cn[i] += 0.5 * (1.0 + erf(-CN_KCN * (r - rc) / rc))
    return log_cn_cut(cn, CN_MAX)


def eeq2019_charges(
    elements: list[int],
    positions_bohr: np.ndarray,
    total_charge: float,
    params: dict,
) -> tuple[np.ndarray, np.ndarray]:
    """EEQ2019 partial charges; returns (charges, coordination numbers)."""
    n = len(elements)
    chi = np.array([params["chi"][z - 1] for z in elements])
    eta = np.array([params["eta"][z - 1] for z in elements])
    rad = np.array([params["rad"][z - 1] for z in elements])
    kcnchi = np.array([params["kcnchi"][z - 1] for z in elements])

    cn = coordination_numbers(elements, positions_bohr, params["covrad"])

    amat = np.zeros((n + 1, n + 1))
    for i in range(n):
        for j in range(i):
            d = positions_bohr[j] - positions_bohr[i]
            r2 = float(d @ d)
            gam = 1.0 / (rad[i] * rad[i] + rad[j] * rad[j])
            val = erf(math.sqrt(r2 * gam)) / math.sqrt(r2)
            amat[i, j] += val
            amat[j, i] += val
        amat[i, i] += eta[i] + math.sqrt(2.0 / math.pi) / rad[i]
    amat[n, :] = 1.0
    amat[:, n] = 1.0
    amat[n, n] = 0.0

    xvec = np.zeros(n + 1)
    for i in range(n):
        c = max(cn[i], 1e-14)
        xvec[i] = -chi[i] + kcnchi[i] * cn[i] / math.sqrt(c)
    xvec[n] = total_charge

    sol = np.linalg.solve(amat, xvec)
    return sol[:n], cn


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def run() -> Report:
    params = load_parameters()

    cases = []
    for name, total_charge, atoms in MOLECULES:
        elements = [z for z, _ in atoms]
        positions = np.array([p for _, p in atoms], dtype=float)
        pos_bohr = positions * BOHR_PER_ANGSTROM
        charges, cn = eeq2019_charges(elements, pos_bohr, total_charge, params)
        cases.append(
            {
                "name": name,
                "elements": elements,
                "coords": [round(float(v), 12) for v in positions.reshape(-1)],
                "total_charge": total_charge,
                "charges": [float(round(v, 12)) for v in charges],
                "coordination_numbers": [float(round(v, 12)) for v in cn],
                "charge_sum": float(round(float(charges.sum()), 12)),
            }
        )

    payload = {
        "ficha": "F87-EEQ2019",
        "oracle": "oracle_eeq2019.py (independent numpy implementation of EEQ2019)",
        "reference": (
            "independent numpy/scipy implementation of the published EEQ2019 "
            "equations (Grimme et al., JCTC 2019); element parameters read from "
            "sci-form's src/charges_eeq/eeq2019_data.rs (generated from "
            "grimme-lab/multicharge eeq2019.f90 + dftd4 Pyykko radii)"
        ),
        "generated": time.strftime("%Y-%m-%d"),
        "conventions": {
            "cn_kcn": CN_KCN,
            "cn_max": CN_MAX,
            "cn_cutoff_bohr": CN_CUTOFF,
            "bohr_per_angstrom": BOHR_PER_ANGSTROM,
            "d3_radius_scale": D3_RADIUS_SCALE,
            "units": "atomic units (positions converted to bohr, charges in e)",
        },
        "cases": cases,
    }

    report = Report(
        ficha="F87-EEQ2019",
        observable="eeq2019_charge_max_abs",
        reference_tool="oracle_eeq2019.py (independent numpy EEQ2019)",
        units="e",
        tolerance=1e-5,
        level="L2",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "EEQ2019 partial charges: erf coordination numbers with a log cut at "
            "CN_max = 8, Gaussian-width Coulomb matrix and CN-dependent "
            "electronegativity, solved with the total-charge constraint. The "
            "oracle uses the exact scipy erf; sci-form uses the Abramowitz & "
            "Stegun 7.1.26 approximation."
        ),
    )
    emit_fixture_report("eeq2019_charges.json", payload, report)
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
    print("fixtures -> tests/fixtures/validation/eeq2019_charges.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
