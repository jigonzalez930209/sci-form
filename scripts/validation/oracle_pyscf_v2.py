#!/usr/bin/env python3
"""PySCF oracle for the v2.x candidates (fichas F101–F115, docs/v2-roadmap.md).

Writes ``tests/fixtures/validation/v2_pyscf.json``; ``tests/v2_validation.rs``
recomputes the sci-form side, so CI checks the parity without PySCF.

Every reference below is the *same method* as the sci-form prototype
(STO-3G, all electrons, no frozen core) so the comparison isolates
implementation errors (gate level L2). Where the roadmap asks for a stronger
physical reference the script also writes an L3 value (coupled response,
electrostatic embedding) that the prototype is expected to miss.

| ficha | V2 id | reference |
|-------|-------|-----------|
| F101  | V2-01 | RKS PBE / PBE0 total energy (same basis, PySCF grid level 4) |
| F102  | V2-02 | RHF + MP2 correlation energy |
| F103  | V2-03 | mechanical-embedding QM/MM energy (vacuum RHF + Mulliken) and, as L3, electrostatic embedding |
| F104  | V2-04 | UHF spin density at the nuclei -> isotropic hyperfine (MHz) |
| F105  | V2-05 | uncoupled (UCHF) polarizability, checked by finite differences; coupled (CPHF) as L3 |
| F106  | V2-06 | Koopmans 1s binding energy (-epsilon of the atom's core MO, UHF) |

Usage::

    .venv-validation/bin/python scripts/validation/oracle_pyscf_v2.py

Dependencies: pyscf (>= 2.4), numpy.
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from pyscf import dft, gto, mp, qmmm, scf

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "tests" / "fixtures" / "validation" / "v2_pyscf.json"

HARTREE_TO_EV = 27.211386245988
HARTREE_TO_KCAL = 627.5094740631
BOHR3_TO_ANG3 = 0.148184711
COULOMB_KCAL_ANG = 332.06371
# Same anchor as src/epr.rs: one unpaired electron in an exact hydrogen 1s
# orbital (spin density 1/pi bohr^-3) gives the 21-cm line, 1420.405751 MHz.
PROTON_HYPERFINE_MHZ = 1420.405751
PROTON_G = 5.5856946893
# Nuclear g-factors of the magnetic isotopes used for A_iso (CODATA / NIST).
NUCLEAR_G = {1: 5.5856946893, 6: 1.40482, 7: 0.40376100, 8: -0.757516, 9: 5.257736}


# ---------------------------------------------------------------------------
# Geometries (experimental, Å). Built from bond lengths/angles so the source
# is auditable: NIST CCCBDB experimental geometries.
# ---------------------------------------------------------------------------


def _bent(center: int, ligand: int, r: float, angle_deg: float):
    half = math.radians(angle_deg) / 2.0
    return [
        (center, (0.0, 0.0, 0.0)),
        (ligand, (r * math.sin(half), 0.0, r * math.cos(half))),
        (ligand, (-r * math.sin(half), 0.0, r * math.cos(half))),
    ]


def _pyramidal(center: int, r: float, angle_deg: float):
    # Three equivalent ligands with H-X-H angle `angle_deg`.
    theta = math.radians(angle_deg)
    # distance of each H from the C3 axis
    rho = r * math.sqrt((2.0 - 2.0 * math.cos(theta)) / 3.0)
    h = math.sqrt(r * r - rho * rho)
    atoms = [(center, (0.0, 0.0, 0.0))]
    for k in range(3):
        phi = 2.0 * math.pi * k / 3.0
        atoms.append((1, (rho * math.cos(phi), rho * math.sin(phi), -h)))
    return atoms


def _tetrahedral(center: int, r: float):
    s = r / math.sqrt(3.0)
    return [(center, (0.0, 0.0, 0.0))] + [
        (1, (s, s, s)),
        (1, (-s, -s, s)),
        (1, (-s, s, -s)),
        (1, (s, -s, -s)),
    ]


def _planar_ch3(r: float):
    atoms = [(6, (0.0, 0.0, 0.0))]
    for k in range(3):
        phi = 2.0 * math.pi * k / 3.0
        atoms.append((1, (r * math.cos(phi), r * math.sin(phi), 0.0)))
    return atoms


def _formaldehyde():
    r_co, r_ch, hch = 1.205, 1.111, 116.1
    half = math.radians(hch) / 2.0
    return [
        (6, (0.0, 0.0, 0.0)),
        (8, (0.0, 0.0, r_co)),
        (1, (r_ch * math.sin(half), 0.0, -r_ch * math.cos(half))),
        (1, (-r_ch * math.sin(half), 0.0, -r_ch * math.cos(half))),
    ]


CLOSED_SHELL = {
    "H2O": _bent(8, 1, 0.9572, 104.52),
    "NH3": _pyramidal(7, 1.012, 106.67),
    "CH4": _tetrahedral(6, 1.087),
    "HF": [(9, (0.0, 0.0, 0.0)), (1, (0.0, 0.0, 0.917))],
    "CO": [(6, (0.0, 0.0, 0.0)), (8, (0.0, 0.0, 1.128))],
    "HCN": [(1, (0.0, 0.0, -1.064)), (6, (0.0, 0.0, 0.0)), (7, (0.0, 0.0, 1.156))],
    "H2CO": _formaldehyde(),
}

RADICALS = {
    "H": [(1, (0.0, 0.0, 0.0))],
    "CH3": _planar_ch3(1.079),
    "NH2": _bent(7, 1, 1.024, 103.4),
    "OH": [(8, (0.0, 0.0, 0.0)), (1, (0.0, 0.0, 0.970))],
}

SYMBOL = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F"}


def build(atoms, spin=0):
    mol = gto.M(
        atom=[(SYMBOL[z], xyz) for z, xyz in atoms],
        basis="sto-3g",
        unit="Angstrom",
        spin=spin,
        cart=True,
        verbose=0,
    )
    return mol


def as_record(atoms):
    return {
        "elements": [z for z, _ in atoms],
        "positions_angstrom": [list(xyz) for _, xyz in atoms],
    }


# ---------------------------------------------------------------------------
# F101 hybrid DFT, F102 MP2
# ---------------------------------------------------------------------------


def dft_energy(mol, xc):
    mf = dft.RKS(mol)
    mf.xc = xc
    mf.grids.level = 4
    mf.conv_tol = 1e-10
    return float(mf.kernel())


def mp2_case(mol):
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-11
    e_hf = mf.kernel()
    pt = mp.MP2(mf)  # all electrons: sci-form has no frozen core
    e_corr, _ = pt.kernel()
    return float(e_hf), float(e_corr)


# ---------------------------------------------------------------------------
# F103 QM/MM
# ---------------------------------------------------------------------------


def qmmm_case():
    qm = CLOSED_SHELL["H2O"]
    mm_positions = [[0.0, 0.0, 3.5], [2.8, 0.0, -1.0]]
    mm_charges = [-0.834, 0.417]
    mol = build(qm)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-11
    e_qm = float(mf.kernel())
    pop = mf.mulliken_pop(verbose=0)[1]
    qm_mm = 0.0
    for (z, xyz), q in zip(qm, pop):
        for pos, q_mm in zip(mm_positions, mm_charges):
            r = math.dist(xyz, pos)
            qm_mm += COULOMB_KCAL_ANG * q * q_mm / r
    mm_mm = COULOMB_KCAL_ANG * mm_charges[0] * mm_charges[1] / math.dist(*mm_positions)
    mechanical = e_qm * HARTREE_TO_KCAL + qm_mm + mm_mm

    # L3: electrostatic embedding, the MM charges polarise the QM density.
    mf_ee = qmmm.mm_charge(scf.RHF(mol), mm_positions, mm_charges, unit="Angstrom")
    mf_ee.conv_tol = 1e-11
    e_ee = float(mf_ee.kernel())  # includes QM-MM electrostatics (not MM-MM)
    electrostatic = e_ee * HARTREE_TO_KCAL + mm_mm
    return {
        "qm": as_record(qm),
        "mm_positions": mm_positions,
        "mm_charges": mm_charges,
        "mechanical_total_kcal_mol": mechanical,
        "qm_vacuum_kcal_mol": e_qm * HARTREE_TO_KCAL,
        "mulliken_charges": [float(q) for q in pop],
        "electrostatic_total_kcal_mol": electrostatic,
    }


# ---------------------------------------------------------------------------
# F104 EPR
# ---------------------------------------------------------------------------


def epr_case(atoms):
    mol = build(atoms, spin=1)
    mf = scf.UHF(mol)
    mf.conv_tol = 1e-11
    mf.kernel()
    dm_a, dm_b = mf.make_rdm1()
    coords = mol.atom_coords()  # Bohr
    ao = mol.eval_gto("GTOval_cart", coords)
    spin = np.einsum("pi,ij,pj->p", ao, dm_a - dm_b, ao)
    entries = []
    for index, (z, _) in enumerate(atoms):
        g = NUCLEAR_G[z]
        a_iso = PROTON_HYPERFINE_MHZ * (g / PROTON_G) * spin[index] * math.pi
        entries.append({"atom_index": index, "element": z, "spin_density_au": float(spin[index]), "a_iso_mhz": float(a_iso)})
    return {"molecule": as_record(atoms), "s2": float(mf.spin_square()[0]), "entries": entries}


# ---------------------------------------------------------------------------
# F105 polarizability
# ---------------------------------------------------------------------------


def uncoupled_polarizability(mf, mol):
    """UCHF alpha = 2 sum_spin sum_ab <a|mu_i|b><b|mu_j|a>/(e_b - e_a)  (bohr^3)."""
    with mol.with_common_orig((0.0, 0.0, 0.0)):
        dip = mol.intor("int1e_r")  # <mu|r|nu>; the electron dipole is -r
    c, e, occ = mf.mo_coeff, mf.mo_energy, mf.mo_occ
    o, v = occ > 0, occ == 0
    d = np.einsum("xpq,pi,qj->xij", dip, c, c)
    de = e[v][None, :] - e[o][:, None]
    per_spin = np.einsum("xab,yab,ab->xy", d[:, o][:, :, v], d[:, o][:, :, v], 1.0 / de)
    return 2.0 * 2.0 * per_spin  # two spins x the factor 2 of -d2E/dF2


def uncoupled_by_finite_difference(mf, mol, h=1e-4):
    """Independent check: -d2/dl2 of sum_spin sum_occ eig(F0 + l*mu)."""
    with mol.with_common_orig((0.0, 0.0, 0.0)):
        dip = mol.intor("int1e_r")
    s = mol.intor("int1e_ovlp")
    f0 = mf.get_fock()
    nocc = int((mf.mo_occ > 0).sum())
    from scipy.linalg import eigh

    def occ_sum(lam, i):
        w = eigh(f0 + lam * dip[i], s, eigvals_only=True)
        return 2.0 * w[:nocc].sum()  # both spins

    diag = []
    for i in range(3):
        second = (occ_sum(h, i) - 2.0 * occ_sum(0.0, i) + occ_sum(-h, i)) / (h * h)
        diag.append(-second)
    return diag


def coupled_polarizability(mol, h=1e-3):
    """CPHF alpha by finite field on the full RHF SCF (bohr^3)."""
    alpha = np.zeros((3, 3))
    with mol.with_common_orig((0.0, 0.0, 0.0)):
        dip = mol.intor("int1e_r")

    def dipole(field):
        mf = scf.RHF(mol)
        h0 = mf.get_hcore()
        mf.get_hcore = lambda *args: h0 + np.einsum("x,xij->ij", field, dip)
        mf.conv_tol = 1e-12
        mf.kernel()
        return mf.dip_moment(unit="AU", verbose=0)

    for i in range(3):
        f = np.zeros(3)
        f[i] = h
        alpha[:, i] = (dipole(f) - dipole(-f)) / (2.0 * h)
    return alpha


def polarizability_case(atoms):
    mol = build(atoms)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.kernel()
    unc = uncoupled_polarizability(mf, mol)
    fd = uncoupled_by_finite_difference(mf, mol)
    for i in range(3):
        if abs(unc[i, i] - fd[i]) > 1e-4 * max(1.0, abs(fd[i])):
            raise SystemExit(f"UCHF formula disagrees with its finite difference: {unc[i, i]} vs {fd[i]}")
    cph = coupled_polarizability(mol)
    return {
        "molecule": as_record(atoms),
        "uncoupled_isotropic_ang3": float(np.trace(unc) / 3.0 * BOHR3_TO_ANG3),
        "uncoupled_tensor_ang3": (unc * BOHR3_TO_ANG3).tolist(),
        "coupled_isotropic_ang3": float(np.trace(cph) / 3.0 * BOHR3_TO_ANG3),
    }


# ---------------------------------------------------------------------------
# F106 Koopmans XPS
# ---------------------------------------------------------------------------


def xps_case(atoms):
    mol = build(atoms)
    mf = scf.UHF(mol)
    mf.conv_tol = 1e-11
    mf.kernel()
    c = mf.mo_coeff[0]
    e = mf.mo_energy[0]
    s = mol.intor("int1e_ovlp")
    heavy = [i for i, (z, _) in enumerate(atoms) if z > 2]
    labels = mol.ao_labels(fmt=False)
    bindings = []
    for k in range(len(heavy)):  # the lowest len(heavy) MOs are the 1s cores
        mulliken = np.einsum("p,pq,q->p", c[:, k], s, c[:, k])
        per_atom = {}
        for p, (atom, *_rest) in enumerate(labels):
            per_atom[atom] = per_atom.get(atom, 0.0) + mulliken[p]
        atom = max(heavy, key=lambda i: per_atom.get(i, 0.0))
        bindings.append({"atom_index": atom, "element": atoms[atom][0], "binding_energy_ev": float(-e[k] * HARTREE_TO_EV)})
    return {"molecule": as_record(atoms), "bindings": sorted(bindings, key=lambda b: b["atom_index"])}


def main() -> int:
    import pyscf

    t0 = time.time()
    out = {
        "provenance": {
            "tool": f"PySCF {pyscf.__version__}",
            "basis": "STO-3G (cartesian), all electrons",
            "script": "scripts/validation/oracle_pyscf_v2.py",
            "geometries": "NIST CCCBDB experimental bond lengths/angles",
            "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
        },
        "dft": [],
        "mp2": [],
        "epr": [],
        "polarizability": [],
        "xps": [],
    }
    for name, atoms in CLOSED_SHELL.items():
        mol = build(atoms)
        e_hf, e_corr = mp2_case(mol)
        out["mp2"].append({"name": name, **as_record(atoms), "hf_energy_hartree": e_hf, "correlation_energy_hartree": e_corr})
        out["dft"].append({"name": name, **as_record(atoms), "pbe_hartree": dft_energy(mol, "pbe,pbe"), "pbe0_hartree": dft_energy(mol, "pbe0")})
        out["xps"].append({"name": name, **xps_case(atoms)})
        if name in ("H2O", "NH3", "CH4", "HF", "CO"):
            out["polarizability"].append({"name": name, **polarizability_case(atoms)})
        print(f"  {name:5s} done ({time.time() - t0:.1f}s)")
    for name, atoms in RADICALS.items():
        out["epr"].append({"name": name, **epr_case(atoms)})
    out["qmmm"] = qmmm_case()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {OUT.relative_to(REPO)} in {time.time() - t0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
