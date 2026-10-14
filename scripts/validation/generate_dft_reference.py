#!/usr/bin/env python3
"""Oracle F77: PySCF RKS LDA(VWN)/PBE STO-3G references for our KS-DFT gate.

Same Å geometries as the Rust gate; basis STO-3G (matches our minimal
basis on H/C/N/O); xc 'lda,vwn' for SVWN, 'pbe,pbe' for PBE.
Reports total energy (Ha), HOMO/LUMO/gap (eV), and convergence.

Usage: generate_dft_reference.py  (writes tests/fixtures/validation/dft_pyscf.json)
"""

from __future__ import annotations

import json
from pathlib import Path

from pyscf import dft, gto

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "tests" / "fixtures" / "validation" / "dft_pyscf.json"

HA_TO_EV = 27.211386

GEOMS = {
    "H2": ([1, 1], [[0, 0, 0], [0.74, 0, 0]]),
    "H2O": ([8, 1, 1], [[0, 0, 0], [0.757, 0.586, 0], [-0.757, 0.586, 0]]),
    "CH4": ([6, 1, 1, 1, 1], [[0, 0, 0], [0.629, 0.629, 0.629], [0.629, -0.629, -0.629],
                              [-0.629, 0.629, -0.629], [-0.629, -0.629, 0.629]]),
    "NH3": ([7, 1, 1, 1], [[0, 0, 0.117], [0, 0.939, -0.313], [0.813, -0.47, -0.313],
                           [-0.813, -0.47, -0.313]]),
    "N2": ([7, 7], [[0, 0, 0], [1.10, 0, 0]]),
    "CO": ([6, 8], [[0, 0, 0], [1.13, 0, 0]]),
}

XCS = {"svwn": "lda,vwn", "pbe": "pbe,pbe"}


def run_case(name: str, elements: list[int], coords: list, xc: str) -> dict:
    mol = gto.Mole()
    mol.atom = [[e, *c] for e, c in zip(elements, coords)]
    mol.unit = "Angstrom"
    mol.basis = "sto-3g"
    mol.build()
    mf = dft.RKS(mol)
    mf.xc = xc
    mf.grids.level = 4
    e_tot = float(mf.kernel())
    assert mf.converged, f"{name}/{xc} did not converge"
    mo_e = mf.mo_energy * HA_TO_EV
    nocc = mol.nelectron // 2
    homo = float(sorted(mo_e[:nocc])[-1])
    lumo = float(sorted(mo_e[nocc:])[0])
    return {
        "name": name, "elements": elements, "coords_A": coords, "xc": xc,
        "total_energy_Ha": e_tot, "homo_eV": homo, "lumo_eV": lumo,
        "gap_eV": lumo - homo, "converged": True,
    }


def main() -> int:
    cases = []
    for name, (elements, coords) in GEOMS.items():
        for key, xc in XCS.items():
            case = run_case(name, elements, coords, xc)
            case["xc"] = key
            cases.append(case)
            print(f"{name}/{key}: E={case['total_energy_Ha']:.8f} gap={case['gap_eV']:.3f}")
    payload = {
        "ficha": "F77",
        "oracle": "generate_dft_reference.py (PySCF RKS STO-3G, grids.level=4)",
        "basis": "sto-3g",
        "energy_tolerance_Ha": 1e-3,
        "gap_tolerance_eV": 0.05,
        "cases": cases,
    }
    OUT.write_text(json.dumps(payload, indent=1) + "\n")
    print(f"wrote {OUT} ({len(cases)} cases)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
