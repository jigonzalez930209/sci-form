#!/usr/bin/env python3
"""Oracle F81: CGA motors vs Rodrigues/quaternion + clifford library.

- Rodrigues rotation matrix (closed form) + quaternion composition:
  reference transformed points for random rigid motions.
- `clifford` (independent GA library, Cl(4,1)): rotor/translator/motor
  sandwich on conformal points — validates our bitmask blade encoding.
- Dihedral rotation about a bond axis via Rodrigues; dihedral angle
  via the standard atan2 formula.

Usage: oracle_cga.py  (writes tests/fixtures/validation/cga_clifford.json)
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
from clifford import Cl  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "tests" / "fixtures" / "validation" / "cga_clifford.json"

layout, blades = Cl(4, 1)
e1, e2, e3, e4, e5 = (blades[f"e{i}"] for i in range(1, 6))
EP, EM = e4, e5
NI = EP + EM
NO = (EM - EP) / 2


def cpoint(p: np.ndarray):
    p = np.asarray(p, float)
    return p[0] * e1 + p[1] * e2 + p[2] * e3 + 0.5 * (p @ p) * NI + NO


def cextract(X) -> list[float]:
    d = float((X | NI).value[0])
    return [float(-(X | e).value[0] / d) for e in (e1, e2, e3)]


def clifford_motor(axis: np.ndarray, angle: float, t: np.ndarray):
    B = axis[0] * (e2 ^ e3) + axis[1] * (e3 ^ e1) + axis[2] * (e1 ^ e2)
    R = np.cos(angle / 2) - np.sin(angle / 2) * B
    T = 1 - 0.5 * ((t[0] * e1 + t[1] * e2 + t[2] * e3) ^ NI)
    return T * R


def rodrigues(axis: np.ndarray, angle: float) -> np.ndarray:
    k = axis / np.linalg.norm(axis)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)


def dihedral_np(p1, p2, p3, p4) -> float:
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    u = b2 / np.linalg.norm(b2)
    return float(np.arctan2(np.dot(np.cross(n1, n2), u), np.dot(n1, n2)))


def main() -> int:
    rng = np.random.default_rng(42)
    motions = []
    for _ in range(16):
        v = rng.standard_normal(3)
        axis = (v / np.linalg.norm(v)).tolist()
        angle = float(rng.uniform(-np.pi, np.pi))
        t = rng.standard_normal(3).tolist()
        pts = rng.standard_normal((4, 3)).tolist()
        R = rodrigues(np.array(axis), angle)
        rod = [(R @ np.array(p) + np.array(t)).tolist() for p in pts]
        M = clifford_motor(np.array(axis), angle, np.array(t))
        cli = [cextract(M * cpoint(np.array(p)) * (~M)) for p in pts]
        motions.append({"axis": axis, "angle": angle, "t": t, "points": pts,
                        "rodrigues_ref": rod, "clifford_ref": cli})
    dihedrals = []
    for _ in range(8):
        a = rng.standard_normal(3)
        axis = a / np.linalg.norm(a)
        b = rng.standard_normal(3)
        angle = float(rng.uniform(-np.pi, np.pi))
        pts = rng.standard_normal((5, 3)).tolist()
        R = rodrigues(axis, angle)
        rot = [(R @ (np.array(p) - np.array(b)) + np.array(b)).tolist() for p in pts]
        p1, p2, p3, p4 = (rng.standard_normal(3) for _ in range(4))
        dihedrals.append({
            "axis_point": b.tolist(), "axis_dir": axis.tolist(), "angle": angle,
            "points": pts, "rotated_ref": rot,
            "dihedral_points": [p1.tolist(), p2.tolist(), p3.tolist(), p4.tolist()],
            "dihedral_ref": dihedral_np(p1, p2, p3, p4),
        })
    payload = {
        "ficha": "F81",
        "oracle": "oracle_cga.py (Rodrigues + clifford Cl(4,1), seed 42)",
        "tolerance": 1e-12,
        "motions": motions,
        "dihedrals": dihedrals,
    }
    OUT.write_text(json.dumps(payload) + "\n")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
