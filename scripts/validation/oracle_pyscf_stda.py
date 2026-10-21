#!/usr/bin/env python3
"""PySCF oracle for sci-form's sTDA UV-Vis spectrum (ficha F77).

Generates `tests/fixtures/validation/pyscf_stda.json` with, for every test
molecule (closed-shell RHF reference, STO-3G basis):

  * `tda_*`      — PySCF `tdscf.TDA` singlet excitation energies (eV) and
                   oscillator strengths (the cross-library reference).
  * `tdhf_*`     — PySCF `tdscf.TDHF` (= `tdscf.TDDFT` for an RHF reference,
                   because a Hartree-Fock reference carries no density
                   functional) energies and oscillator strengths.
  * `tda_active_*` — the same TDA restricted to sci-form's *active space*
                   (occupied/virtual window of `StdaConfig`), which is the
                   like-for-like reference for the sTDA matrix.
  * `tda_noex_active_*` — the same active-space TDA **without the exchange
                   term**.  sci-form's sTDA has no exchange contribution at
                   all (`StdaConfig::ax` is declared but never used), so this
                   reference separates the two error sources: the missing
                   exchange term and the monopole approximation of the
                   Coulomb interaction.
  * `stda_ref_*` — a faithful Python port of sci-form's sTDA algorithm
                   (monopole Coulomb, no exchange) evaluated on PySCF's
                   orbitals.  This is the *method control* that pins the
                   sci-form implementation; the residual against PySCF's
                   TDA is then attributable to the method, not to a coding
                   error.

Usage:
    python scripts/validation/oracle_pyscf_stda.py [all]
"""

from __future__ import annotations

import os

# Pin the BLAS/OpenMP thread count *before* numpy is imported: PySCF's
# diagonalisation otherwise sums in a thread-dependent order, which changes
# the last digits of the reference values and makes the oracle
# non-reproducible (the fixtures must be byte-identical on every run).
for _thread_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_thread_var] = "1"

import argparse
import sys
import time

import numpy as np

from common import HARTREE2EV, Report, emit_fixture_report

try:
    from pyscf import ao2mo, gto, scf as pyscf_scf, tdscf
except ImportError as exc:  # pragma: no cover
    print(f"PySCF is required for this oracle: {exc}", file=sys.stderr)
    sys.exit(1)

PYSCF_VERSION = getattr(__import__("pyscf"), "__version__", "unknown")
REFERENCE = f"PySCF {PYSCF_VERSION}"

# sci-form's `StdaConfig` used by the gate (documented in the fixture).
OCC_WINDOW_EV = 20.0
VIRT_WINDOW_EV = 20.0
N_ROOTS = 4

# sci-form's sTDA parameters (see `spectroscopy::stda_uvvis`).
HARTREE_TO_EV = HARTREE2EV
DAMPING_ETA = 0.3
PRESCREEN_THRESHOLD = 1e-6
CORE_FLOOR_EV = -20.0

MOLECULES = [
    {
        "name": "H2",
        "elements": [1, 1],
        "coords": [0.0, 0.0, 0.0, 0.0, 0.0, 0.74],
    },
    {
        "name": "HF",
        "elements": [9, 1],
        "coords": [0.0, 0.0, 0.0, 0.0, 0.0, 0.917],
    },
    {
        "name": "H2O",
        "elements": [8, 1, 1],
        "coords": [
            0.0,
            0.0,
            0.1173,
            0.0,
            0.7572,
            -0.4692,
            0.0,
            -0.7572,
            -0.4692,
        ],
    },
    {
        "name": "NH3",
        "elements": [7, 1, 1, 1],
        "coords": [
            0.0,
            0.0,
            0.1173,
            0.0,
            0.9377,
            -0.2739,
            0.8121,
            -0.4689,
            -0.2739,
            -0.8121,
            -0.4689,
            -0.2739,
        ],
    },
    {
        "name": "CH4",
        "elements": [6, 1, 1, 1, 1],
        "coords": [
            0.0,
            0.0,
            0.0,
            0.6276,
            0.6276,
            0.6276,
            -0.6276,
            0.6276,
            0.6276,
            0.6276,
            -0.6276,
            0.6276,
            -0.6276,
            -0.6276,
            -0.6276,
        ],
    },
    {
        "name": "CH2O",
        "elements": [6, 8, 1, 1],
        "coords": [
            0.0,
            0.0,
            0.5349,
            0.0,
            0.0,
            -0.6701,
            0.0,
            0.9426,
            1.1165,
            0.0,
            -0.9426,
            1.1165,
        ],
    },
    {
        "name": "C2H4",
        "elements": [6, 6, 1, 1, 1, 1],
        "coords": [
            0.0,
            0.0,
            0.6695,
            0.0,
            0.0,
            -0.6695,
            0.0,
            0.9289,
            1.2321,
            0.0,
            -0.9289,
            1.2321,
            0.0,
            0.9289,
            -1.2321,
            0.0,
            -0.9289,
            -1.2321,
        ],
    },
]


def atom_spec(case: dict) -> list:
    atoms = []
    for z, (x, y, z_) in zip(
        case["elements"], np.array(case["coords"]).reshape(-1, 3)
    ):
        atoms.append((int(z), (float(x), float(y), float(z_))))
    return atoms


def select_active_space(mol, mo_energy: np.ndarray) -> tuple[list[int], list[int]]:
    """Mirror of `spectroscopy::stda_uvvis::select_active_space`."""
    n_occ = mol.nelectron // 2
    nao = mol.nao_nr()
    homo = mo_energy[n_occ - 1]
    lumo = mo_energy[n_occ] if n_occ < nao else homo + 1.0

    occ_cutoff = max(
        homo - OCC_WINDOW_EV / HARTREE_TO_EV, CORE_FLOOR_EV / HARTREE_TO_EV
    )
    virt_cutoff = lumo + VIRT_WINDOW_EV / HARTREE_TO_EV

    occ_idx = [i for i in range(n_occ) if mo_energy[i] >= occ_cutoff]
    if not occ_idx:
        n_include = min(max(n_occ, 1), 3)
        occ_idx = list(range(n_occ - n_include, n_occ))
    virt_idx = [a for a in range(n_occ, nao) if mo_energy[a] <= virt_cutoff]
    return occ_idx, virt_idx


def mo_dipole_matrices(mol, coeff: np.ndarray) -> np.ndarray:
    """D^x_{pq} = <p| r_x |q> in atomic units (origin at the coordinates)."""
    int1e_r = mol.intor("int1e_r")
    return np.einsum("xij,ip,jq->xpq", int1e_r, coeff, coeff)


def transition_dipoles_singlet(
    x: np.ndarray, occ_idx: list[int], virt_idx: list[int], dip_mo: np.ndarray
) -> np.ndarray:
    """Transition dipole of a singlet TDA state.

    `x` has shape (n_occ_active, n_virt_active); the singlet transition
    density is T_pq = 2 sum_ia X_ia C_pa C_qi - sum_ia X_ia C_pi C_qa.
    """
    mol_dip = 0.0
    dip = np.zeros(3)
    for i_idx, i in enumerate(occ_idx):
        for a_idx, a in enumerate(virt_idx):
            xi = x[i_idx, a_idx]
            if xi == 0.0:
                continue
            for comp in range(3):
                dip[comp] += xi * (2.0 * dip_mo[comp, a, i] - dip_mo[comp, i, a])
    del mol_dip
    return dip


def tda_matrix(
    mo_eri: np.ndarray,
    mo_energy: np.ndarray,
    occ_idx: list[int],
    virt_idx: list[int],
    exchange: bool,
) -> np.ndarray:
    """A matrix of the (Tamm-Dancoff) singles response in the MO basis.

    A_{ia,jb} = delta_ij delta_ab (eps_a - eps_i) + 2 (ia|jb) - (ij|ab)

    with the chemist-notation MO ERIs; the last (exchange) term is dropped
    when `exchange` is False.
    """
    n_occ, n_virt = len(occ_idx), len(virt_idx)
    dim = n_occ * n_virt
    amat = np.zeros((dim, dim))
    for i1, i in enumerate(occ_idx):
        for a1, a in enumerate(virt_idx):
            k1 = i1 * n_virt + a1
            amat[k1, k1] = mo_energy[a] - mo_energy[i]
            for i2, j in enumerate(occ_idx):
                for a2, b in enumerate(virt_idx):
                    k2 = i2 * n_virt + a2
                    val = 2.0 * mo_eri[i, a, j, b]
                    if exchange:
                        val -= mo_eri[i, j, a, b]
                    amat[k1, k2] += val
    return amat


def solve_singles(amat: np.ndarray, n_occ: int, n_virt: int, n_roots: int):
    """Lowest roots of the singles matrix -> (energies, X vectors)."""
    values, vectors = np.linalg.eigh(amat)
    order = np.argsort(values)
    energies = []
    xvecs = []
    for root in order[: min(n_roots, len(values))]:
        energies.append(values[root])
        xvecs.append(vectors[:, root].reshape(n_occ, n_virt))
    return np.array(energies), xvecs


def stda_reference(
    mol, mf, mo_energy: np.ndarray, coeff: np.ndarray, overlap: np.ndarray,
    mo_eri: np.ndarray | None = None, exchange: bool = False,
):
    """Faithful port of sci-form's `compute_stda` (monopole Coulomb).

    With `exchange=True` the Tamm-Dancoff exchange term `-2 ax (ij|ab)`
    (ax = 0.5, i.e. `-(ij|ab)`) is added to the A matrix; this is the
    variant implemented by sci-form's `compute_stda_with_exchange`.

    Returns (energies in eV, oscillator strengths, active space).
    """
    n_occ = mol.nelectron // 2
    nao = mol.nao_nr()
    occ_idx, virt_idx = select_active_space(mol, mo_energy)
    n_active_occ, n_active_virt = len(occ_idx), len(virt_idx)

    # atom -> AO index map (`ao_labels` carries the atom index of every AO;
    # `mol.ao_loc` is a per-shell offset in PySCF 2.x and must not be used
    # as a per-atom start index).
    basis_to_atom = np.zeros(nao, dtype=int)
    for mu, label in enumerate(mol.ao_labels()):
        basis_to_atom[mu] = int(label.split()[0])

    sc = overlap @ coeff
    # transition charges q[i][a][atom]
    q = np.zeros((n_active_occ, n_active_virt, mol.natm))
    for i_l, i in enumerate(occ_idx):
        for a_l, a in enumerate(virt_idx):
            for mu in range(nao):
                q[i_l, a_l, basis_to_atom[mu]] += coeff[mu, i] * sc[mu, a]

    # damped monopole Coulomb kernel
    positions = mol.atom_coords()
    gamma = np.zeros((mol.natm, mol.natm))
    for a in range(mol.natm):
        for b in range(mol.natm):
            if a == b:
                gamma[a, b] = DAMPING_ETA
            else:
                r2 = float(np.sum((positions[a] - positions[b]) ** 2))
                avg = 1.0 / (2.0 * DAMPING_ETA) + 1.0 / (2.0 * DAMPING_ETA)
                gamma[a, b] = 1.0 / np.sqrt(r2 + avg * avg)

    pairs = [(i, a) for i in occ_idx for a in virt_idx]
    dim = len(pairs)
    amat = np.zeros((dim, dim))
    norms = np.linalg.norm(q.reshape(dim, mol.natm), axis=1)
    for k1, (i1, a1) in enumerate(pairs):
        il1, al1 = occ_idx.index(i1), virt_idx.index(a1)
        amat[k1, k1] = mo_energy[a1] - mo_energy[i1]
        for k2, (i2, a2) in enumerate(pairs):
            # Cauchy-Schwarz prescreening of the *monopole* term only; the
            # exchange integrals are not screened.
            if norms[k1] * norms[k2] >= PRESCREEN_THRESHOLD:
                il2, al2 = occ_idx.index(i2), virt_idx.index(a2)
                jval = 0.0
                for atom_a in range(mol.natm):
                    qia = q[il1, al1, atom_a]
                    if abs(qia) < PRESCREEN_THRESHOLD:
                        continue
                    for atom_b in range(mol.natm):
                        jval += qia * gamma[atom_a, atom_b] * q[il2, al2, atom_b]
                amat[k1, k2] += 2.0 * jval
            if exchange:
                amat[k1, k2] -= mo_eri[i1, i2, a1, a2]

    values, vectors = np.linalg.eigh(amat)
    order = np.argsort(values)
    energies = []
    strengths = []
    for root in order[: min(N_ROOTS, dim)]:
        energy = values[root]
        if energy < 0.0:
            continue
        xvec = vectors[:, root]
        tdm = np.zeros(3)
        for k, (i, a) in enumerate(pairs):
            x = xvec[k]
            if abs(x) < 1e-10:
                continue
            il, al = occ_idx.index(i), virt_idx.index(a)
            for atom in range(mol.natm):
                tdm += x * q[il, al, atom] * positions[atom]
        tdm *= np.sqrt(2.0)
        energies.append(energy * HARTREE_TO_EV)
        strengths.append((2.0 / 3.0) * energy * float(np.dot(tdm, tdm)))
    del n_occ
    return np.array(energies), np.array(strengths), occ_idx, virt_idx


def run() -> Report:
    cases = []
    for case in MOLECULES:
        mol = gto.M(
            atom=atom_spec(case),
            basis="sto-3g",
            unit="Angstrom",
            charge=0,
            spin=0,
            verbose=0,
        )
        mf = pyscf_scf.RHF(mol)
        mf.conv_tol = 1e-12
        mf.max_cycle = 512
        e_tot = float(mf.kernel())
        mo_energy = np.asarray(mf.mo_energy)
        coeff = np.asarray(mf.mo_coeff)
        overlap = np.asarray(mf.get_ovlp())
        nao = mol.nao_nr()

        # --- PySCF references (full active space) -------------------------
        td = tdscf.TDA(mf)
        td.nstates = N_ROOTS
        td.kernel()
        tda_ev = [float(e) * HARTREE2EV for e in td.e]
        tda_osc = [float(x) for x in td.oscillator_strength()]

        tdhf = tdscf.TDHF(mf)
        tdhf.nstates = N_ROOTS
        tdhf.kernel()
        tdhf_ev = [float(e) * HARTREE2EV for e in tdhf.e]
        tdhf_osc = [float(x) for x in tdhf.oscillator_strength()]

        # --- active space (sci-form's window selection) -------------------
        occ_idx, virt_idx = select_active_space(mol, mo_energy)
        frozen = [
            i
            for i in range(mol.nelectron // 2 + mol.nao_nr() - mol.nelectron // 2)
            if i not in occ_idx and i < mol.nelectron // 2
        ]
        frozen += [
            a
            for a in range(mol.nelectron // 2, nao)
            if a not in virt_idx
        ]

        td_act = tdscf.TDA(mf)
        td_act.frozen = frozen
        td_act.nstates = N_ROOTS
        td_act.kernel()
        tda_act_ev = [float(e) * HARTREE2EV for e in td_act.e]
        tda_act_osc = [float(x) for x in td_act.oscillator_strength()]

        # --- MO ERIs and the exchange-free active-space TDA ---------------
        mo_eri = ao2mo.full(mol, coeff, compact=False).reshape(nao, nao, nao, nao)
        dip_mo = mo_dipole_matrices(mol, coeff)

        amat_noex = tda_matrix(mo_eri, mo_energy, occ_idx, virt_idx, exchange=False)
        ev_noex, xvecs_noex = solve_singles(
            amat_noex, len(occ_idx), len(virt_idx), N_ROOTS
        )
        noex_ev = [float(e) * HARTREE2EV for e in ev_noex]
        noex_osc = []
        for xvec in xvecs_noex:
            dip = transition_dipoles_singlet(xvec, occ_idx, virt_idx, dip_mo)
            noex_osc.append((2.0 / 3.0) * float(ev_noex[len(noex_osc)]) * float(np.dot(dip, dip)))

        # --- sTDA reference port ------------------------------------------
        stda_ev, stda_osc, occ_ref, virt_ref = stda_reference(
            mol, mf, mo_energy, coeff, overlap
        )
        assert occ_ref == occ_idx and virt_ref == virt_idx

        # --- sTDA reference port with the exchange term -------------------
        stda_ex_ev, stda_ex_osc, _, _ = stda_reference(
            mol, mf, mo_energy, coeff, overlap, mo_eri=mo_eri, exchange=True
        )

        cases.append(
            {
                "name": case["name"],
                "elements": case["elements"],
                "coords": case["coords"],
                "n_basis": int(nao),
                "n_electrons": int(mol.nelectron),
                "e_tot": e_tot,
                "mo_energy": [float(x) for x in mo_energy],
                "occ_active": [int(i) for i in occ_idx],
                "virt_active": [int(a) for a in virt_idx],
                "tda_ev": tda_ev,
                "tda_oscillator": tda_osc,
                "tdhf_ev": tdhf_ev,
                "tdhf_oscillator": tdhf_osc,
                "tda_active_ev": tda_act_ev,
                "tda_active_oscillator": tda_act_osc,
                "tda_noex_active_ev": noex_ev,
                "tda_noex_active_oscillator": noex_osc,
                "stda_ref_ev": [float(x) for x in stda_ev],
                "stda_ref_oscillator": [float(x) for x in stda_osc],
                "stda_ex_ref_ev": [float(x) for x in stda_ex_ev],
                "stda_ex_ref_oscillator": [float(x) for x in stda_ex_osc],
            }
        )

    payload = {
        "reference": REFERENCE,
        "basis": "sto-3g",
        "generated": time.strftime("%Y-%m-%d"),
        "config": {
            "occ_window_ev": OCC_WINDOW_EV,
            "virt_window_ev": VIRT_WINDOW_EV,
            "n_roots": N_ROOTS,
            "damping_eta": DAMPING_ETA,
            "prescreen_threshold": PRESCREEN_THRESHOLD,
            "core_floor_ev": CORE_FLOOR_EV,
        },
        "cases": cases,
    }
    report = Report(
        ficha="F77",
        observable="stda_vs_pyscf_tda_max_abs_ev",
        reference_tool=f"{REFERENCE} tdscf.TDA / tdscf.TDHF",
        units="eV / dimensionless (oscillator strength)",
        tolerance=0.01,
        level="L2",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "sTDA (monopole Coulomb) vs PySCF TD-SCF, plus the exchange-enabled "
            "variant (A_x = -2 ax (ij|ab), ax = 0.5) that sci-form implements in "
            "`compute_stda_with_exchange`. The implementation-parity gate (sci-form "
            "sTDA vs the Python port stda_ref / stda_ex_ref) is tight; the "
            "cross-library deviation is dominated by the monopole approximation "
            "of the Coulomb interaction. Measured by "
            "tests/pyscf_parity.rs::f77_stda_uvvis."
        ),
    )
    emit_fixture_report("pyscf_stda.json", payload, report)
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
    print("fixtures -> tests/fixtures/validation/pyscf_stda.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
