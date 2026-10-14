#!/usr/bin/env python3
"""Oracle F84: RDKit UFF gradients for the dynamics force path.

Replicates the F16 construction (AddHs, fixture coords) and additionally
records ff.CalcGrad() (dE/dx, kcal/mol/A). Our LiveMolecularSystem with
MdBackend::Uff must reproduce -grad as forces (F16 already gates the
energy to 0.0; this gates the force path used by Verlet).

Usage: oracle_dynamics.py [N]  (writes tests/fixtures/validation/uff_forces_rdkit.json)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "tests" / "fixtures" / "validation" / "uff_rdkit.json"
OUT = REPO / "tests" / "fixtures" / "validation" / "uff_forces_rdkit.json"


def main() -> int:
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    fixture = json.loads(SRC.read_text())
    cases = []
    for rec in fixture["cases"]:
        if len(cases) >= limit:
            break
        smi = rec["smiles"]
        mol0 = Chem.MolFromSmiles(smi)
        if mol0 is None:
            continue
        mol = Chem.AddHs(mol0)
        if [a.GetAtomicNum() for a in mol.GetAtoms()] != rec["elements"]:
            continue
        conf = Chem.Conformer(mol.GetNumAtoms())
        xyz = rec["coords"]
        for i in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(i, (xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]))
        mol.RemoveAllConformers()
        mol.AddConformer(conf, assignId=True)
        try:
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=0)
            if ff is None:
                continue
            grad = list(ff.CalcGrad())
        except Exception:
            continue
        cases.append({"smiles": smi, "elements": rec["elements"], "coords": xyz,
                      "reference_energy": rec["reference_energy"],
                      "rdkit_forces": [-g for g in grad]})
    payload = {
        "ficha": "F84",
        "oracle": "oracle_dynamics.py (RDKit UFF CalcGrad, F16 construction)",
        "force_tolerance": 1e-3,
        "cases": cases,
    }
    OUT.write_text(json.dumps(payload) + "\n")
    print(f"wrote {OUT} ({len(cases)} cases)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
