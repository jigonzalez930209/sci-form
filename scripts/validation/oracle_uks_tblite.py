#!/usr/bin/env python3
"""Generate `tests/fixtures/validation/uks_tblite.json` from tblite 0.7.0.

Runs the tblite `Calculator` (`GFN2-xTB`, `uhf=N`) on a set of open-shell
molecules and records total energies. NOTE: tblite's Python interface takes
positions in BOHR (mctc convention); geometries below are in Angstrom and are
converted. NOTE 2: tblite 0.7.0 wires no spin-polarization container for GFN2
(`uhf` only splits alpha/beta occupations, the SCF cycles are identical to
restricted for doublets); our UKS adds the on-site spin interaction of
tblite's own `spin.f90`/`data/spin.f90`, so the gate compares
`E_uks - E_spin` vs these references (relaxation is second order, <0.05 eV).

Usage:
    ~/.venvs/sci-form-validation/bin/python scripts/validation/oracle_uks_tblite.py [--out FILE]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "tests/fixtures/validation/uks_tblite.json"
ANGSTROM_TO_BOHR = 1.8897261246257702

# name -> (uhf, [(Z, x, y, z)]) with positions in Angstrom.
SYSTEMS: dict[str, tuple[int, list[tuple[int, float, float, float]]]] = {
    "OH": (1, [(8, 0.0, 0.0, 0.0), (1, 0.0, 0.0, 0.97)]),
    "CH3": (
        1,
        [
            (6, 0.0, 0.0, 0.0),
            (1, 1.09, 0.0, 0.0),
            (1, -0.545, 0.944, 0.0),
            (1, -0.545, -0.944, 0.0),
        ],
    ),
    "NO": (1, [(7, 0.0, 0.0, 0.0), (8, 0.0, 0.0, 1.15)]),
    "O2": (2, [(8, 0.0, 0.0, 0.0), (8, 0.0, 0.0, 1.21)]),
    # Closed-shell control for the uhf=0 equivalence check.
    "H2O": (
        0,
        [(8, 0.0, 0.0, 0.1173), (1, 0.0, 0.7572, -0.4692), (1, 0.0, -0.7572, -0.4692)],
    ),
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    from tblite.interface import Calculator

    cases = []
    for name, (uhf, atoms) in SYSTEMS.items():
        numbers = np.array([z for z, _, _, _ in atoms])
        pos_a = np.array([[x, y, z] for _, x, y, z in atoms])
        entry = {
            "name": name,
            "uhf": uhf,
            "elements": [int(z) for z, _, _, _ in atoms],
            "coords": [c for _, x, y, z in atoms for c in (x, y, z)],
        }
        for method in ("GFN2-xTB", "GFN1-xTB"):
            calc = Calculator(
                method=method,
                numbers=numbers,
                positions=pos_a * ANGSTROM_TO_BOHR,
                charge=0.0,
                uhf=uhf if uhf > 0 else None,
            )
            res = calc.singlepoint()
            e_ha = float(res.get("energy"))
            key = "tblite_energy_hartree" if method == "GFN2-xTB" else "tblite_gfn1_energy_hartree"
            entry[key] = e_ha
            print(f"{method} {name}: uhf={uhf} E={e_ha:.9f} Ha")
        cases.append(entry)

    payload = {
        "reference": "tblite 0.7.0 Calculator GFN2-xTB uhf=N (positions converted A->bohr)",
        "generated": "scripts/validation/oracle_uks_tblite.py",
        "cases": cases,
    }
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print("wrote", args.out)


if __name__ == "__main__":
    main()
