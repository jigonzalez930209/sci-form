#!/usr/bin/env python3
"""Independent MMFF94 energy oracle used to validate the Rust implementation.

This script is *not* a port of the Rust code.  It re-derives the MMFF94 energy
expression from the published functional forms of Halgren, T.,
*J. Comput. Chem.* **1996**, *17*, 490-519:

    E_bond   = 0.5 * 143.9325 * kb * dr^2 * (1 + cs*dr + 7/12*cs^2*dr^2)   (eq. 2)
    E_angle  = 0.5 * 0.043844 * ka * dT^2 * (1 + cb*dT)                     (eq. 3)
    E_tors   = 0.5 * (V1*(1+cos p) + V2*(1-cos 2p) + V3*(1+cos 3p))         (MMFF.IV)
    E_vdW    = eps * (1.07 R*/(R + 0.07 R*))^7
                          * (1.12 R*^7/(R^7 + 0.12 R*^7) - 2)              (MMFF.V)

with cs = -2 A^-1 and cb = -0.007 deg^-1 (Halgren's rounded value; RDKit uses
the unrounded -0.4 rad^-1 = -0.006981317 deg^-1 - see the note below).

The parameter values come from the published MMFF94 tables, taken from the
RDKit source tree (BSD-3-Clause) at
``Code/ForceField/MMFF/Params.cpp``; the parameter data itself is
Copyright (c) Merck and Co., Inc., 1994, 1995, 1996.

The script writes ``tests/fixtures/validation/mmff94_rdkit.json``, which the
Rust test ``tests/mmff94_validation.rs`` compares against.

Usage
-----
    python scripts/validation/oracle_mmff94.py [--params-cpp PATH] [--out PATH]

``--params-cpp`` defaults to ``.cache/mmff94/Params.cpp`` and is downloaded
from the RDKit repository when missing.

Note on the cubic bend constant
-------------------------------
Halgren (1996), eq. 3, states ``cb = -0.007 deg^-1``.  RDKit hard-codes
``cb = -0.006981317 deg^-1``, i.e. exactly ``-0.4 rad^-1``.  The two differ by
0.27 % of the anharmonic correction.  This oracle uses the *exact* value
``-0.4 rad^-1`` so that the Rust implementation and this reference agree to
machine precision; the difference to the rounded published value is far below
the 0.1 kcal/mol parity tolerance and is documented in the Rust module.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
CACHE_DIR = ROOT / ".cache" / "mmff94"
OUT = ROOT / "tests" / "fixtures" / "validation" / "mmff94_rdkit.json"

PARAMS_URL = (
    "https://raw.githubusercontent.com/rdkit/rdkit/master/Code/ForceField/MMFF/Params.cpp"
)

# ── MMFF94 constants (Halgren 1996) ──────────────────────────────────────────

MDYNE_A_TO_KCAL_MOL = 143.9325
DEG2RAD = math.pi / 180.0
RAD2DEG = 180.0 / math.pi
CUBIC_STRETCH = -2.0          # A^-1
CUBIC_BEND_RAD = -0.4         # rad^-1 (RDKit's unrounded value of Halgren's -0.007 deg^-1)
CUBIC_BEND_DEG = CUBIC_BEND_RAD * DEG2RAD   # = -0.006981317... deg^-1
VDW_C4 = 181.16

# Neighbour-relation codes used by the topological distance matrix.
RELATION_1_2 = 0
RELATION_1_3 = 1
RELATION_1_4 = 2
RELATION_1_X = 3


# ── upstream parameter tables ────────────────────────────────────────────────


def _fetch_params_cpp(params_cpp: Path | None) -> tuple[str, str]:
    if params_cpp is not None:
        text = params_cpp.read_text()
        return text, hashlib.sha256(text.encode()).hexdigest()
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    target = CACHE_DIR / "Params.cpp"
    if not target.exists():
        print(f"downloading {PARAMS_URL} -> {target}", file=sys.stderr)
        with urllib.request.urlopen(PARAMS_URL) as fh:
            target.write_bytes(fh.read())
    text = target.read_text()
    return text, hashlib.sha256(text.encode()).hexdigest()


def _cpp_string(src: str, name: str) -> str:
    """Concatenate and unescape the C++ string literal(s) bound to ``name``."""
    m = re.search(r"const\s+std::string\s+" + re.escape(name) + r"\s*(\[\s*\])?\s*=\s*", src)
    if not m:
        raise KeyError(name)
    is_array = bool(m.group(1))
    pos = m.end()
    if is_array:
        if src[pos] != "{":
            raise ValueError(name)
        pos += 1
        stop = "}"
    else:
        stop = ";"
    pieces: list[str] = []
    while pos < len(src):
        ch = src[pos]
        if ch == '"':
            end = pos + 1
            buf = []
            while end < len(src) and src[end] != '"':
                if src[end] == "\\":
                    nxt = src[end + 1]
                    buf.append({"n": "\n", "t": "\t", "r": "\r", "\\": "\\", '"': '"'}.get(nxt, nxt))
                    end += 2
                else:
                    buf.append(src[end])
                    end += 1
            pieces.append("".join(buf))
            pos = end + 1
            continue
        if ch == stop:
            break
        pos += 1
    return "".join(pieces)


def _rows(blob: str) -> list[list[str]]:
    out = []
    for line in blob.split("\n"):
        if not line or line.startswith("*"):
            continue
        fields = [f.strip() for f in line.split("\t")]
        if fields and fields[0]:
            out.append(fields)
    return out


class Tables:
    """The published MMFF94 parameter tables."""

    def __init__(self, src: str):
        self.sha256 = hashlib.sha256(src.encode()).hexdigest()

        # atom type equivalence levels (MMFF.I Table III)
        self.eq: dict[int, tuple[int, int, int, int]] = {}
        for row in _rows(_cpp_string(src, "defaultMMFFDef")):
            if len(row) < 6:
                continue
            t = int(row[1])
            if t not in self.eq:
                self.eq[t] = (int(row[2]), int(row[3]), int(row[4]), int(row[5]))

        # atom type properties (MMFF.II Table I)
        self.prop: dict[int, dict[str, int]] = {}
        for row in _rows(_cpp_string(src, "defaultMMFFProp")):
            if len(row) < 9:
                continue
            self.prop[int(row[0])] = {
                "atno": int(row[1]), "crd": int(row[2]), "val": int(row[3]),
                "pilp": int(row[4]), "mltb": int(row[5]), "arom": int(row[6]),
                "linh": int(row[7]), "sbmb": int(row[8]),
            }

        # aromatic atom types
        m = re.search(r"defaultMMFFArom\s*=\s*\{(.*?)\}", src, re.S)
        self.aromatic = set(int(v) for v in re.findall(r"\d+", m.group(1)))

        # bond stretching parameters
        self.bonds: dict[tuple[int, int, int], tuple[float, float]] = {}
        for row in _rows(_cpp_string(src, "defaultMMFFBond")):
            if len(row) < 5:
                continue
            bt, i, j = int(row[0]), int(row[1]), int(row[2])
            self.bonds[(min(i, j), max(i, j), bt)] = (float(row[3]), float(row[4]))

        # default bond rule
        self.bndk: dict[tuple[int, int], tuple[float, float]] = {}
        for row in _rows(_cpp_string(src, "defaultMMFFBndk")):
            if len(row) < 4:
                continue
            zi, zj = int(row[0]), int(row[1])
            self.bndk[(min(zi, zj), max(zi, zj))] = (float(row[2]), float(row[3]))

        # Herschbach-Laurie (Badger's rule) parameters
        self.hl: dict[tuple[int, int], tuple[float, float, float]] = {}
        for row in _rows(_cpp_string(src, "defaultMMFFHerschbachLaurie")):
            if len(row) < 5:
                continue
            i, j = int(row[0]), int(row[1])
            self.hl[(min(i, j), max(i, j))] = (float(row[2]), float(row[3]), float(row[4]))

        # covalent radius / Pauling electronegativity
        self.cov: dict[int, tuple[float, float]] = {}
        for row in _rows(_cpp_string(src, "defaultMMFFCovRadPauEle")):
            if len(row) < 3:
                continue
            self.cov[int(row[0])] = (float(row[1]), float(row[2]))

        # angle bending parameters
        self.angles: dict[tuple[int, int, int, int], tuple[float, float]] = {}
        for row in _rows(_cpp_string(src, "defaultMMFFAngleData")):
            if len(row) < 6:
                continue
            at, i, j, k = int(row[0]), int(row[1]), int(row[2]), int(row[3])
            self.angles[(at, i, j, k)] = (float(row[4]), float(row[5]))

        # torsion parameters
        self.tor: dict[tuple[int, int, int, int, int], tuple[float, float, float]] = {}
        for row in _rows(_cpp_string(src, "defaultMMFFTor")):
            if len(row) < 8:
                continue
            key = (int(row[0]), int(row[1]), int(row[2]), int(row[3]), int(row[4]))
            self.tor[key] = (float(row[5]), float(row[6]), float(row[7]))

        # van der Waals parameters
        vdw_rows = _rows(_cpp_string(src, "defaultMMFFVdW"))
        self.vdw_globals = tuple(float(v) for v in vdw_rows[0][:5])
        self.vdw: dict[int, tuple[float, float, float, float, str]] = {}
        for row in vdw_rows[1:]:
            if len(row) < 6:
                continue
            self.vdw[int(row[0])] = (float(row[1]), float(row[2]), float(row[3]),
                                     float(row[4]), row[5])

    # ── lookups ──

    def eq_level(self, atom_type: int, level: int) -> int:
        return self.eq.get(atom_type, (0, 0, 0, 0))[level]

    def prop_of(self, atom_type: int, key: str) -> int:
        return self.prop.get(atom_type, {}).get(key, 0)

    def bond_param(self, bond_type: int, i: int, j: int):
        return self.bonds.get((min(i, j), max(i, j), bond_type))

    def angle_param(self, angle_type: int, i: int, j: int, k: int):
        """Five-stage equivalence-level search (MMFF.I note 68, p. 519)."""
        for level in range(4):
            ci = self.eq_level(i, level)
            ck = self.eq_level(k, level)
            if ci > ck:
                ci, ck = ck, ci
            hit = self.angles.get((angle_type, ci, j, ck))
            if hit is not None:
                return hit
        return None

    def torsion_param(self, tor_type: int, secondary: int, i: int, j: int, k: int, l: int):
        """Five-stage equivalence-level search (MMFF.IV, p. 609)."""
        can_tor = tor_type
        max_iter = 5
        iter_ = 0
        found = None
        while True:
            keep_going = (iter_ < max_iter) and ((found is None) or (max_iter == 4))
            keep_going = keep_going or (iter_ == 4 and tor_type == 5 and secondary)
            if not keep_going:
                break
            if max_iter == 5 and iter_ == 4:
                max_iter = 4
                iter_ = 0
                can_tor = secondary
            i_wild, l_wild = iter_, iter_
            if iter_ == 1:
                i_wild, l_wild = 1, 3
            elif iter_ == 2:
                i_wild, l_wild = 3, 1
            ci = self.eq_level(i, i_wild)
            cj = j
            ck = k
            cl = self.eq_level(l, l_wild)
            if cj > ck:
                cj, ck = ck, cj
                ci, cl = cl, ci
            elif cj == ck and ci > cl:
                ci, cl = cl, ci
            hit = self.tor.get((can_tor, ci, cj, ck, cl))
            if hit is not None:
                found = hit
                if max_iter == 4:
                    break
            iter_ += 1
        return found

    def vdw_pair(self, i: int, j: int) -> tuple[float, float]:
        """Pair parameters (R*_ij, eps_ij) from the MMFF.V combination rules."""
        power, b_coef, beta, darad, daeps = self.vdw_globals
        ai, ni, ai2, gi, dai = self.vdw[i]
        aj, nj, aj2, gj, daj = self.vdw[j]
        r_i = ai2 * (ai ** power)
        r_j = aj2 * (aj ** power)
        gamma = (r_i - r_j) / (r_i + r_j)
        if dai == "D" or daj == "D":
            r_ij = 0.5 * (r_i + r_j)
        else:
            r_ij = 0.5 * (r_i + r_j) * (1.0 + b_coef * (1.0 - math.exp(-beta * gamma * gamma)))
        eps = (VDW_C4 * gi * gj * ai * aj /
               ((math.sqrt(ai / ni) + math.sqrt(aj / nj)) * r_ij ** 6))
        if (dai == "D" and daj == "A") or (dai == "A" and daj == "D"):
            r_ij *= darad
            eps *= daeps
        return r_ij, eps


# ── periodic table helpers ───────────────────────────────────────────────────


def periodic_table_row(atomic_num: int) -> int:
    if 3 <= atomic_num <= 10:
        return 1
    if 11 <= atomic_num <= 18:
        return 2
    if 19 <= atomic_num <= 36:
        return 3
    if 37 <= atomic_num <= 54:
        return 4
    return 0


def periodic_table_row_hl(atomic_num: int) -> int:
    row = 0
    if atomic_num == 2:
        row = 1
    elif 3 <= atomic_num <= 10:
        row = 2
    elif 11 <= atomic_num <= 18:
        row = 3
    elif 19 <= atomic_num <= 36:
        row = 4
    elif 37 <= atomic_num <= 54:
        row = 5
    if (21 <= atomic_num <= 30) or (39 <= atomic_num <= 48):
        row *= 10
    return row


# ── molecular system ─────────────────────────────────────────────────────────


class Atom:
    def __init__(self, element: int, mmff_type: int):
        self.element = element
        self.type = mmff_type


class Bond:
    def __init__(self, i: int, j: int, order: str):
        self.i = i
        self.j = j
        self.order = order  # "single", "double", "triple", "aromatic"


class Molecule:
    def __init__(self, name: str, smiles: str, atoms, bonds, coords):
        self.name = name
        self.smiles = smiles
        self.atoms = atoms
        self.bonds = bonds
        self.coords = coords
        self.n = len(atoms)

    def neighbors(self, idx: int) -> list[int]:
        out = []
        for b in self.bonds:
            if b.i == idx:
                out.append(b.j)
            elif b.j == idx:
                out.append(b.i)
        return out

    def bonded(self, a: int, b: int) -> bool:
        return any((bnd.i == a and bnd.j == b) or (bnd.i == b and bnd.j == a)
                   for bnd in self.bonds)

    def bond_between(self, a: int, b: int):
        for bnd in self.bonds:
            if (bnd.i == a and bnd.j == b) or (bnd.i == b and bnd.j == a):
                return bnd
        return None

    def degree(self, idx: int) -> int:
        return len(self.neighbors(idx))

    def topological_distance(self, a: int, b: int) -> int:
        """Bond-path distance, or a large number when disconnected."""
        if a == b:
            return 0
        seen = {a}
        frontier = [a]
        dist = 0
        while frontier:
            dist += 1
            nxt = []
            for node in frontier:
                for nb in self.neighbors(node):
                    if nb == b:
                        return dist
                    if nb not in seen:
                        seen.add(nb)
                        nxt.append(nb)
            frontier = nxt
        return 10 ** 8


# ── MMFF94 parameter resolution ──────────────────────────────────────────────


class Resolver:
    def __init__(self, tables: Tables, mol: Molecule):
        self.t = tables
        self.mol = mol

    # ── bond type / angle type / torsion type ──

    def bond_type(self, a: int, b: int) -> int:
        bnd = self.mol.bond_between(a, b)
        if bnd is None:
            return 0
        pa = self.t.prop_of(self.mol.atoms[a].type, "sbmb")
        pb = self.t.prop_of(self.mol.atoms[b].type, "sbmb")
        aa = self.t.prop_of(self.mol.atoms[a].type, "arom")
        ab = self.t.prop_of(self.mol.atoms[b].type, "arom")
        if bnd.order == "single" and ((pa and pb) or (aa and ab)):
            return 1
        return 0

    def angle_type(self, i: int, j: int, k: int) -> int:
        bond_sum = self.bond_type(i, j) + self.bond_type(j, k)
        angle_type = bond_sum
        size = self._angle_in_ring_3or4(i, j, k)
        if size:
            angle_type = size
            if bond_sum:
                angle_type += bond_sum + size - 2
        return angle_type

    def torsion_type(self, i: int, j: int, k: int, l: int):
        bnd_jk = self.mol.bond_between(j, k)
        bt_ij = self.bond_type(i, j)
        bt_jk = self.bond_type(j, k)
        bt_kl = self.bond_type(k, l)
        torsion_type = bt_jk
        second = 0
        if bt_jk == 0 and bnd_jk.order == "single" and (bt_ij == 1 or bt_kl == 1):
            torsion_type = 2
        size = self._torsion_in_ring_4or5(i, j, k, l)
        if size == 4 and not (self.mol.bonded(i, k) or self.mol.bonded(j, l)):
            second = torsion_type
            torsion_type = 4
        elif size == 5 and any(self.mol.atoms[x].type == 1 for x in (i, j, k, l)):
            second = torsion_type
            torsion_type = 5
        return torsion_type, second

    def _angle_in_ring_3or4(self, i: int, j: int, k: int) -> int:
        if not (self.mol.bonded(i, j) and self.mol.bonded(j, k)):
            return 0
        if self.mol.bonded(k, i):
            return 3
        s1 = {n for n in self.mol.neighbors(i) if n != j}
        s2 = {n for n in self.mol.neighbors(k) if n != j}
        return 4 if (s1 & s2) else 0

    def _torsion_in_ring_4or5(self, i: int, j: int, k: int, l: int) -> int:
        if not (self.mol.bonded(i, j) and self.mol.bonded(j, k) and self.mol.bonded(k, l)):
            return 0
        if self.mol.bonded(l, i):
            return 4
        s1 = {n for n in self.mol.neighbors(i) if n != j}
        s2 = {n for n in self.mol.neighbors(l) if n != k}
        return 5 if (s1 & s2) else 0

    # ── bond stretching ──

    def bond_params(self, a: int, b: int) -> tuple[float, float]:
        bt = self.bond_type(a, b)
        hit = self.t.bond_param(bt, self.mol.atoms[a].type, self.mol.atoms[b].type)
        if hit is not None:
            return hit
        return self._bond_empirical_rule(a, b)

    def _bond_empirical_rule(self, a: int, b: int) -> tuple[float, float]:
        """MMFF.V eq. (18)/(19): Badger's rule with a covalent-radius correction."""
        z1 = self.mol.atoms[a].element
        z2 = self.mol.atoms[b].element
        cov1 = self.t.cov.get(z1)
        cov2 = self.t.cov.get(z2)
        if cov1 is None or cov2 is None:
            return (0.0, 0.0)
        c = 0.050 if (z1 == 1 or z2 == 1) else 0.085
        n = 1.4
        r0 = cov1[0] + cov2[0] - c * (abs(cov1[1] - cov2[1]) ** n)
        bndk = self.t.bndk.get((min(z1, z2), max(z1, z2)))
        if bndk is not None:
            coeff = bndk[0] / r0
            kb = bndk[1] * coeff ** 6
        else:
            hl = self.t.hl.get((min(periodic_table_row_hl(z1), periodic_table_row_hl(z2)),
                                max(periodic_table_row_hl(z1), periodic_table_row_hl(z2))))
            if hl is None:
                return (0.0, r0)
            kb = 10.0 ** (-(r0 - hl[0]) / hl[1])
        return (kb, r0)

    # ── angle bending ──

    def angle_params(self, i: int, j: int, k: int) -> tuple[float, float]:
        at = self.angle_type(i, j, k)
        hit = self.t.angle_param(at, self.mol.atoms[i].type,
                                  self.mol.atoms[j].type, self.mol.atoms[k].type)
        if hit is not None and hit[0] != 0.0:
            return hit
        return self._angle_empirical_rule(i, j, k, hit)

    def _angle_empirical_rule(self, i: int, j: int, k: int, old) -> tuple[float, float]:
        """MMFF.V eq. (20) plus the reference-angle rules."""
        z = [self.mol.atoms[x].element for x in (i, j, k)]
        central = self.mol.atoms[j].type
        crd = self.t.prop_of(central, "crd")
        val = self.t.prop_of(central, "val")
        mltb = self.t.prop_of(central, "mltb")
        linh = self.t.prop_of(central, "linh")
        ring = self._angle_in_ring_3or4(i, j, k)

        if old is not None:
            theta0 = old[1]
        else:
            theta0 = 120.0
            if crd == 4:
                theta0 = 109.45
            elif crd == 2:
                if z[1] == 8:
                    theta0 = 105.0
                elif linh == 1:
                    theta0 = 180.0
            elif crd == 3:
                if val == 3 and mltb == 0:
                    theta0 = 107.0 if z[1] == 7 else 92.0
            if ring == 3:
                theta0 = 60.0
            elif ring == 4:
                theta0 = 90.0

        # Table VI of MMFF.V (p. 628): atomic Z and the accompanying
        # coefficient used by the force-constant rule.
        zc = [0.0, 0.0, 0.0]
        cc = [0.0, 0.0, 0.0]
        for idx in range(3):
            zi = z[idx]
            if zi == 1:
                zc[idx], cc[idx] = 1.395, 0.0
            elif zi == 6:
                zc[idx], cc[idx] = 2.494, 1.016
            elif zi == 7:
                zc[idx], cc[idx] = 2.711, 1.113
            elif zi == 8:
                zc[idx], cc[idx] = 3.045, 1.337
            elif zi == 9:
                zc[idx], cc[idx] = 2.847, 0.0
            elif zi == 14:
                zc[idx], cc[idx] = 2.350, 0.811
            elif zi == 15:
                zc[idx], cc[idx] = 2.350, 1.068
            elif zi == 16:
                zc[idx], cc[idx] = 2.980, 1.249
            elif zi == 17:
                zc[idx], cc[idx] = 2.909, 1.078
            elif zi == 35:
                zc[idx], cc[idx] = 3.017, 0.0
            elif zi == 53:
                zc[idx], cc[idx] = 3.086, 0.0

        r0_ij = self.bond_params(i, j)[1]
        r0_jk = self.bond_params(j, k)[1]
        d = (r0_ij - r0_jk) ** 2 / (r0_ij + r0_jk) ** 2
        beta = 1.75
        if ring == 4:
            beta *= 0.85
        elif ring == 3:
            beta *= 0.05
        theta0_rad = DEG2RAD * theta0
        ka = beta * zc[0] * cc[1] * zc[2] / ((r0_ij + r0_jk) * theta0_rad ** 2 * math.exp(2.0 * d))
        return (ka, theta0)

    # ── torsion ──

    def torsion_params(self, i: int, j: int, k: int, l: int) -> tuple[float, float, float]:
        tor_type, second = self.torsion_type(i, j, k, l)
        hit = self.t.torsion_param(tor_type, second, self.mol.atoms[i].type,
                                   self.mol.atoms[j].type, self.mol.atoms[k].type,
                                   self.mol.atoms[l].type)
        if hit is not None:
            return hit
        return self._torsion_empirical_rule(j, k)

    def _torsion_empirical_rule(self, j: int, k: int) -> tuple[float, float, float]:
        """MMFF.V rules (a)-(h)."""
        t = self.t
        jt = self.mol.atoms[j].type
        kt = self.mol.atoms[k].type
        jprop = t.prop.get(jt, {})
        kprop = t.prop.get(kt, {})
        bnd = self.mol.bond_between(j, k)
        u = [0.0, 0.0]
        v = [0.0, 0.0]
        w = [0.0, 0.0]
        z = [self.mol.atoms[j].element, self.mol.atoms[k].element]
        for idx in (0, 1):
            if z[idx] == 6:
                u[idx], v[idx] = 2.0, 2.12
            elif z[idx] == 7:
                u[idx], v[idx] = 2.0, 1.5
            elif z[idx] == 8:
                u[idx], v[idx], w[idx] = 2.0, 0.2, 2.0
            elif z[idx] == 14:
                u[idx], v[idx] = 1.25, 1.22
            elif z[idx] == 15:
                u[idx], v[idx] = 1.25, 2.40
            elif z[idx] == 16:
                u[idx], v[idx], w[idx] = 1.25, 0.49, 8.0

        n_jk = float((jprop.get("crd", 0) - 1) * (kprop.get("crd", 0) - 1))
        v1 = v2 = v3 = 0.0

        if jprop.get("linh", 0) or kprop.get("linh", 0):
            pass                                            # rule (a)
        elif (jt in t.aromatic and kt in t.aromatic and bnd.order == "aromatic"):
            beta = 3.0 if ((jprop.get("val") == 3 and kprop.get("val") == 4) or
                           (jprop.get("val") == 4 and kprop.get("val") == 3)) else 6.0
            pi_jk = 0.5 if (jprop.get("pilp") == 0 and kprop.get("pilp") == 0) else 0.3
            v2 = beta * pi_jk * math.sqrt(u[0] * u[1])        # rule (b)
        elif bnd.order == "double":
            beta = 6.0
            pi_jk = 1.0 if (jprop.get("mltb") == 2 and kprop.get("mltb") == 2) else 0.4
            v2 = beta * pi_jk * math.sqrt(u[0] * u[1])        # rule (c)
        elif jprop.get("crd") == 4 and kprop.get("crd") == 4:
            v3 = math.sqrt(v[0] * v[1]) / n_jk               # rule (d)
        elif jprop.get("crd") == 4 and kprop.get("crd") != 4:
            if ((kprop.get("crd") == 3 and ((kprop.get("val") in (4, 34)) or kprop.get("mltb"))) or
                    (kprop.get("crd") == 2 and (kprop.get("val") == 3 or kprop.get("mltb")))):
                pass                                        # rule (e), flat
            else:
                v3 = math.sqrt(v[0] * v[1]) / n_jk
        elif kprop.get("crd") == 4 and jprop.get("crd") != 4:
            if ((jprop.get("crd") == 3 and ((jprop.get("val") in (4, 34)) or jprop.get("mltb"))) or
                    (jprop.get("crd") == 2 and (jprop.get("val") == 3 or jprop.get("mltb")))):
                pass                                        # rule (f), flat
            else:
                v3 = math.sqrt(v[0] * v[1]) / n_jk
        elif ((bnd.order == "single" and jprop.get("mltb") and kprop.get("mltb")) or
              (jprop.get("mltb") and kprop.get("pilp")) or
              (jprop.get("pilp") and kprop.get("mltb"))):
            if jprop.get("pilp") and kprop.get("pilp"):
                pass                                        # rule (g)(1)
            elif jprop.get("pilp") and kprop.get("mltb"):
                beta = 6.0
                if kprop.get("mltb") == 1:
                    pi_jk = 0.5
                elif periodic_table_row(z[0]) == 2 and periodic_table_row(z[1]) == 2:
                    pi_jk = 0.3
                else:
                    pi_jk = 0.15
                v2 = beta * pi_jk * math.sqrt(u[0] * u[1])
            elif kprop.get("pilp") and jprop.get("mltb"):
                beta = 6.0
                if jprop.get("mltb") == 1:
                    pi_jk = 0.5
                elif periodic_table_row(z[0]) == 2 and periodic_table_row(z[1]) == 2:
                    pi_jk = 0.3
                else:
                    pi_jk = 0.15
                v2 = beta * pi_jk * math.sqrt(u[0] * u[1])
            elif ((jprop.get("mltb") == 1 or kprop.get("mltb") == 1) and
                  not (z[0] == 6 and z[1] == 6)):
                v2 = 6.0 * 0.4 * math.sqrt(u[0] * u[1])
            else:
                v2 = 6.0 * 0.15 * math.sqrt(u[0] * u[1])
        else:
            if (z[0] in (8, 16)) and (z[1] in (8, 16)):
                v2 = -math.sqrt(w[0] * w[1])                 # rule (h)
            else:
                v3 = math.sqrt(v[0] * v[1]) / n_jk
        return (v1, v2, v3)


# ── energy expressions ───────────────────────────────────────────────────────


def bond_energy(kb: float, r0: float, dist: float) -> float:
    dr = dist - r0
    cs = CUBIC_STRETCH
    return 0.5 * MDYNE_A_TO_KCAL_MOL * kb * dr * dr * (1.0 + cs * dr + (7.0 / 12.0) * cs * cs * dr * dr)


def angle_energy(ka: float, theta0_deg: float, theta_deg: float, linear: bool) -> float:
    if linear:
        return MDYNE_A_TO_KCAL_MOL * ka * (1.0 + math.cos(theta_deg * DEG2RAD))
    dt = theta_deg - theta0_deg
    c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD
    return 0.5 * c2 * ka * dt * dt * (1.0 + CUBIC_BEND_DEG * dt)


def torsion_energy(v1: float, v2: float, v3: float, cos_phi: float) -> float:
    cos2 = 2.0 * cos_phi * cos_phi - 1.0
    cos3 = cos_phi * (2.0 * cos2 - 1.0)
    return 0.5 * (v1 * (1.0 + cos_phi) + v2 * (1.0 - cos2) + v3 * (1.0 + cos3))


def vdw_energy(dist: float, r_star: float, eps: float) -> float:
    dist7 = dist ** 7
    r7 = r_star ** 7
    a_term = (1.07 * r_star / (dist + 0.07 * r_star)) ** 7
    b_term = 1.12 * r7 / (dist7 + 0.12 * r7) - 2.0
    return eps * a_term * b_term


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _norm(a):
    return math.sqrt(_dot(a, a))


def _unit(a):
    n = _norm(a)
    if n < 1e-12:
        return (0.0, 0.0, 1.0)
    return (a[0] / n, a[1] / n, a[2] / n)


def angle_degrees(p1, p2, p3) -> float:
    u = _sub(p1, p2)
    v = _sub(p3, p2)
    cos = _dot(u, v) / (_norm(u) * _norm(v))
    cos = max(-1.0, min(1.0, cos))
    return math.acos(cos) * RAD2DEG


def torsion_cos_phi(p1, p2, p3, p4) -> float:
    r1 = _sub(p1, p2)
    r2 = _sub(p3, p2)
    r3 = _sub(p2, p3)
    r4 = _sub(p4, p3)
    t1 = _cross(r1, r2)
    t2 = _cross(r3, r4)
    l1 = _norm(t1)
    l2 = _norm(t2)
    if l1 < 1e-10 or l2 < 1e-10:
        return 0.0
    c = _dot(t1, t2) / (l1 * l2)
    return max(-1.0, min(1.0, c))


def distance(coords, a, b) -> float:
    pa = coords[3 * a: 3 * a + 3]
    pb = coords[3 * b: 3 * b + 3]
    return _norm(_sub(pa, pb))


# ── term enumeration ─────────────────────────────────────────────────────────


def enumerate_terms(mol: Molecule, resolver: Resolver):
    """Enumerate the MMFF94 bonded / non-bonded terms the way RDKit does."""
    bonds = []
    for b in mol.bonds:
        kb, r0 = resolver.bond_params(b.i, b.j)
        bonds.append((b.i, b.j, kb, r0))

    angles = []
    for j in range(mol.n):
        if mol.degree(j) == 1:
            continue
        nbrs = mol.neighbors(j)
        for a in range(len(nbrs)):
            for b in range(a + 1, len(nbrs)):
                i, k = nbrs[a], nbrs[b]
                ka, theta0 = resolver.angle_params(i, j, k)
                linear = resolver.t.prop_of(mol.atoms[j].type, "linh") > 0
                angles.append((i, j, k, ka, theta0, linear))

    torsions = []
    for b in mol.bonds:
        j, k = b.i, b.j
        # RDKit's torsion SMARTS: [!$(*#*)&!D1]~[!$(*#*)&!D1]
        if mol.degree(j) == 1 or mol.degree(k) == 1:
            continue
        if mol.atoms[j].type == 42 or mol.atoms[k].type == 42:  # NC (sp N, triple bonded)
            continue
        for i in mol.neighbors(j):
            if i == k:
                continue
            for l in mol.neighbors(k):
                if l == j or l == i:
                    continue
                v1, v2, v3 = resolver.torsion_params(i, j, k, l)
                if v1 == 0.0 and v2 == 0.0 and v3 == 0.0:
                    continue
                torsions.append((i, j, k, l, v1, v2, v3))

    vdw = []
    for i in range(mol.n):
        for j in range(i + 1, mol.n):
            if mol.topological_distance(i, j) < 3:
                continue
            r_star, eps = resolver.t.vdw_pair(mol.atoms[i].type, mol.atoms[j].type)
            vdw.append((i, j, r_star, eps))
    return bonds, angles, torsions, vdw


def energy_breakdown(mol: Molecule, resolver: Resolver):
    bonds, angles, torsions, vdw = enumerate_terms(mol, resolver)
    coords = mol.coords
    e_bond = sum(bond_energy(kb, r0, distance(coords, i, j)) for (i, j, kb, r0) in bonds)
    e_angle = 0.0
    for (i, j, k, ka, theta0, linear) in angles:
        theta = angle_degrees(coords[3 * i: 3 * i + 3], coords[3 * j: 3 * j + 3],
                               coords[3 * k: 3 * k + 3])
        e_angle += angle_energy(ka, theta0, theta, linear)
    e_tor = 0.0
    for (i, j, k, l, v1, v2, v3) in torsions:
        c = torsion_cos_phi(coords[3 * i: 3 * i + 3], coords[3 * j: 3 * j + 3],
                            coords[3 * k: 3 * k + 3], coords[3 * l: 3 * l + 3])
        e_tor += torsion_energy(v1, v2, v3, c)
    e_vdw = sum(vdw_energy(distance(coords, i, j), rs, eps) for (i, j, rs, eps) in vdw)
    return {
        "bond": e_bond,
        "angle": e_angle,
        "torsion": e_tor,
        "vdw": e_vdw,
        "total": e_bond + e_angle + e_tor + e_vdw,
    }, {"bonds": len(bonds), "angles": len(angles), "torsions": len(torsions), "vdw_pairs": len(vdw)}


# ── geometry construction ────────────────────────────────────────────────────


def perp(u):
    """Any unit vector perpendicular to the unit vector ``u``."""
    helper = (1.0, 0.0, 0.0) if abs(u[0]) < 0.9 else (0.0, 1.0, 0.0)
    return _unit(_cross(u, helper))


def tetra3(u, length):
    """Directions of the three remaining substituents of a tetrahedral centre.

    ``u`` is the unit vector from the centre to the substituent that is
    already placed.
    """
    a = perp(u)
    b = _cross(u, a)
    out = []
    for k in range(3):
        ang = math.radians(120.0 * k)
        v = [(-1.0 / 3.0) * u[q]
             + (2.0 * math.sqrt(2.0) / 3.0) * (math.cos(ang) * a[q] + math.sin(ang) * b[q])
             for q in range(3)]
        out.append([length * c for c in _unit(v)])
    return out


def tetra2(u1, u2, length):
    """Directions of the two remaining substituents of a tetrahedral centre."""
    s = [-(u1[q] + u2[q]) / 2.0 for q in range(3)]
    n = _unit(_cross(u1, u2))
    w = math.sqrt(2.0 / 3.0)
    return [[length * c for c in _unit([s[q] + w * n[q] for q in range(3)])],
            [length * c for c in _unit([s[q] - w * n[q] for q in range(3)])]]


def place_angle(center, ref_dir, angle_deg, length, plane_normal=None):
    """Place an atom at ``length`` from ``center``, ``angle_deg`` away from ``ref_dir``."""
    u = _unit(ref_dir)
    if plane_normal is None:
        a = perp(u)
    else:
        a = _unit(_cross(_unit(plane_normal), u))
    ang = math.radians(angle_deg)
    return [center[q] + length * (math.cos(ang) * u[q] + math.sin(ang) * a[q]) for q in range(3)]


def flatten(coords) -> list[float]:
    return [c for xyz in coords for c in xyz]


# ── test molecules ───────────────────────────────────────────────────────────


def methane():
    # C type 1 (CR); H type 5 (HC)
    atoms = [Atom(6, 1)] + [Atom(1, 5) for _ in range(4)]
    bonds = [Bond(0, 1, "single"), Bond(0, 2, "single"),
             Bond(0, 3, "single"), Bond(0, 4, "single")]
    coords = [[0.0, 0.0, 0.0]]
    u = (0.0, 0.0, 1.0)
    coords.append([1.09 * c for c in u])
    for d in tetra3(u, 1.09):
        coords.append(d)
    return Molecule("methane", "C", atoms, bonds, flatten(coords))


def ethanol():
    # C0 type 1 (CR), C1 type 1 (CR), O2 type 6 (OR), H type 5 / 21
    atoms = [Atom(6, 1), Atom(6, 1), Atom(8, 6),
             Atom(1, 5), Atom(1, 5), Atom(1, 5), Atom(1, 5), Atom(1, 5), Atom(1, 21)]
    bonds = [Bond(0, 1, "single"), Bond(1, 2, "single"),
             Bond(0, 3, "single"), Bond(0, 4, "single"), Bond(0, 5, "single"),
             Bond(1, 6, "single"), Bond(1, 7, "single"), Bond(2, 8, "single")]
    r_cc, r_co, r_ch, r_oh = 1.531, 1.431, 1.094, 0.971
    coords = [[0.0, 0.0, 0.0], [r_cc, 0.0, 0.0]]
    # hydroxyl O: C0-C1-O2 = 108 deg
    o2 = place_angle(coords[1], (-1.0, 0.0, 0.0), 108.0, r_co, plane_normal=(0.0, 0.0, 1.0))
    coords.append(o2)
    # methyl hydrogens around C0 (tetrahedral with the C0-C1 bond)
    for d in tetra3(_unit(_sub(coords[1], coords[0])), r_ch):
        coords.append([coords[0][q] + d[q] for q in range(3)])
    # methylene hydrogens around C1 (tetrahedral with C0 and O2)
    u1 = _unit(_sub(coords[0], coords[1]))
    u2 = _unit(_sub(coords[2], coords[1]))
    for d in tetra2(u1, u2, r_ch):
        coords.append([coords[1][q] + d[q] for q in range(3)])
    # hydroxyl hydrogen: C1-O2-H = 108 deg
    h8 = place_angle(coords[2], _sub(coords[1], coords[2]), 108.0, r_oh,
                     plane_normal=(0.0, 0.0, 1.0))
    coords.append(h8)
    return Molecule("ethanol", "CCO", atoms, bonds, flatten(coords))


def acetic_acid():
    # C0 type 1 (CR), C1 type 3 (C=O), O2 type 7 (O=C), O3 type 6 (OR), HOCO type 24
    atoms = [Atom(6, 1), Atom(6, 3), Atom(8, 7), Atom(8, 6),
             Atom(1, 5), Atom(1, 5), Atom(1, 5), Atom(1, 24)]
    bonds = [Bond(0, 1, "single"), Bond(1, 2, "double"), Bond(1, 3, "single"),
             Bond(0, 4, "single"), Bond(0, 5, "single"), Bond(0, 6, "single"),
             Bond(3, 7, "single")]
    r_cc, r_co, r_coh, r_ch, r_oh = 1.520, 1.214, 1.364, 1.090, 1.064
    coords = [[0.0, 0.0, 0.0], [r_cc, 0.0, 0.0]]
    # carbonyl oxygen: C0-C1-O2 = 123 deg, in the xy plane
    coords.append(place_angle(coords[1], (-1.0, 0.0, 0.0), 123.0, r_co,
                              plane_normal=(0.0, 0.0, 1.0)))
    # hydroxyl oxygen: C0-C1-O3 = 110 deg, on the other side of the C-C axis
    o3 = place_angle(coords[1], (-1.0, 0.0, 0.0), 110.0, r_coh,
                     plane_normal=(0.0, 0.0, -1.0))
    coords.append(o3)
    # methyl hydrogens around C0
    for d in tetra3(_unit(_sub(coords[1], coords[0])), r_ch):
        coords.append([coords[0][q] + d[q] for q in range(3)])
    # hydroxyl hydrogen
    coords.append(place_angle(coords[3], _sub(coords[1], coords[3]), 107.0, r_oh,
                              plane_normal=(0.0, 0.0, 1.0)))
    return Molecule("acetic_acid", "CC(=O)O", atoms, bonds, flatten(coords))


def dimethyl_ether():
    atoms = [Atom(6, 1), Atom(8, 6), Atom(6, 1),
             Atom(1, 5), Atom(1, 5), Atom(1, 5),
             Atom(1, 5), Atom(1, 5), Atom(1, 5)]
    bonds = [Bond(0, 1, "single"), Bond(1, 2, "single"),
             Bond(0, 3, "single"), Bond(0, 4, "single"), Bond(0, 5, "single"),
             Bond(2, 6, "single"), Bond(2, 7, "single"), Bond(2, 8, "single")]
    r_co, r_ch = 1.410, 1.086
    coords = [[0.0, 0.0, 0.0]]
    # O in the middle: C0-O1-C2 = 112 deg
    coords.append(place_angle(coords[0], (1.0, 0.0, 0.0), 56.0, r_co,
                              plane_normal=(0.0, 0.0, 1.0)))
    coords.append(place_angle(coords[1], _sub(coords[0], coords[1]), 112.0, r_co,
                              plane_normal=(0.0, 0.0, 1.0)))
    for base, ref in ((0, 1), (2, 1)):
        for d in tetra3(_unit(_sub(coords[ref], coords[base])), r_ch):
            coords.append([coords[base][q] + d[q] for q in range(3)])
    return Molecule("dimethyl_ether", "COC", atoms, bonds, flatten(coords))


def benzene():
    atoms = [Atom(6, 37) for _ in range(6)]
    atoms += [Atom(1, 5) for _ in range(6)]
    bonds = []
    for k in range(6):
        bonds.append(Bond(k, (k + 1) % 6, "aromatic"))
    for k in range(6):
        bonds.append(Bond(k, 6 + k, "single"))
    r_cc, r_ch = 1.397, 1.081
    coords = []
    for k in range(6):
        ang = math.radians(60.0 * k)
        coords.append([r_cc * math.cos(ang), r_cc * math.sin(ang), 0.0])
    for k in range(6):
        ang = math.radians(60.0 * k)
        coords.append([(r_cc + r_ch) * math.cos(ang), (r_cc + r_ch) * math.sin(ang), 0.0])
    return Molecule("benzene", "c1ccccc1", atoms, bonds, flatten(coords))


def water():
    # O type 6 (OR); H type 21 (HOR)
    atoms = [Atom(8, 6), Atom(1, 21), Atom(1, 21)]
    bonds = [Bond(0, 1, "single"), Bond(0, 2, "single")]
    r_oh = 0.958
    coords = [[0.0, 0.0, 0.0]]
    h1 = place_angle(coords[0], (1.0, 0.0, 0.0), 104.5, r_oh, plane_normal=(0.0, 0.0, 1.0))
    coords.append(h1)
    coords.append(place_angle(coords[0], (1.0, 0.0, 0.0), 104.5, r_oh,
                              plane_normal=(0.0, 0.0, -1.0)))
    return Molecule("water", "O", atoms, bonds, flatten(coords))


def methylamine():
    # C type 1 (CR), N type 8 (NR), H type 5 (HC) and 23 (HNR)
    atoms = [Atom(6, 1), Atom(7, 8),
             Atom(1, 5), Atom(1, 5), Atom(1, 5), Atom(1, 23), Atom(1, 23)]
    bonds = [Bond(0, 1, "single"),
             Bond(0, 2, "single"), Bond(0, 3, "single"), Bond(0, 4, "single"),
             Bond(1, 5, "single"), Bond(1, 6, "single")]
    r_cn, r_ch, r_nh = 1.474, 1.093, 1.012
    coords = [[0.0, 0.0, 0.0], [r_cn, 0.0, 0.0]]
    for d in tetra3(_unit(_sub(coords[1], coords[0])), r_ch):
        coords.append([coords[0][q] + d[q] for q in range(3)])
    for d in tetra2(_unit(_sub(coords[0], coords[1])),
                    _unit(_sub(coords[1], coords[0])), r_nh):
        coords.append([coords[1][q] + d[q] for q in range(3)])
    return Molecule("methylamine", "CN", atoms, bonds, flatten(coords))


MOLECULES = [methane, ethanol, acetic_acid, dimethyl_ether, benzene, water, methylamine]


# ── main ─────────────────────────────────────────────────────────────────────


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--params-cpp", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    src, sha = _fetch_params_cpp(args.params_cpp)
    tables = Tables(src)

    payload = {
        "description": (
            "MMFF94 reference energies computed by scripts/validation/oracle_mmff94.py "
            "from the published MMFF94 parameter tables (RDKit "
            "Code/ForceField/MMFF/Params.cpp, BSD-3-Clause; parameter data "
            "Copyright (c) Merck and Co., Inc., 1994, 1995, 1996)."
        ),
        "upstream": PARAMS_URL,
        "upstream_sha256": sha,
        "units": "kcal/mol",
        "notes": [
            "Energy = bond stretch + angle bend + torsion + van der Waals (buffered 14-7).",
            "Electrostatic, stretch-bend and out-of-plane terms are not part of the "
            "Rust implementation under test and are therefore not included here.",
            "Cubic bend constant: -0.4 rad^-1 (RDKit value); Halgren 1996 rounds "
            "this to -0.007 deg^-1.",
        ],
        "molecules": [],
    }

    for builder in MOLECULES:
        mol = builder()
        resolver = Resolver(tables, mol)
        energies, counts = energy_breakdown(mol, resolver)
        payload["molecules"].append({
            "name": mol.name,
            "smiles": mol.smiles,
            "elements": [a.element for a in mol.atoms],
            "mmff94_types": [a.type for a in mol.atoms],
            "bonds": [[b.i, b.j, b.order] for b in mol.bonds],
            "coords": [round(c, 12) for c in mol.coords],
            "counts": counts,
            "energy": {k: v for k, v in energies.items()},
        })
        print("%-14s total = %14.8f  %s" % (mol.name, energies["total"], counts))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print("wrote", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
