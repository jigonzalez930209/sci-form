#!/usr/bin/env python3
"""Reference MMFF94 energies and partial charges from RDKit itself.

The independent oracle (``oracle_mmff94.py``) checks each energy term against
the published functional forms.  This companion checks the *full* energy and
the partial charges against the actual RDKit implementation, using the same
coordinates as the oracle fixture.

It reads ``tests/fixtures/validation/mmff94_rdkit.json`` and writes
``tests/fixtures/validation/mmff94_rdkit_lib.json``.

Usage:
    python scripts/validation/oracle_mmff94_rdkit.py

The RDKit version used for the committed fixture is recorded in
``validation/refs.lock.json``.
"""

from __future__ import annotations

import json
from pathlib import Path

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parent.parent.parent
SOURCE = ROOT / "tests" / "fixtures" / "validation" / "mmff94_rdkit.json"
OUT = ROOT / "tests" / "fixtures" / "validation" / "mmff94_rdkit_lib.json"


def main() -> int:
    source = json.loads(SOURCE.read_text())
    payload = {
        "description": (
            "Full MMFF94 energies and partial charges computed by RDKit "
            "MMFFGetMoleculeForceField/MMFFGetMoleculeProperties on the "
            "coordinates of mmff94_rdkit.json."
        ),
        "rdkit_version": rdkit.__version__,
        "variant": "MMFF94",
        "molecules": [],
    }

    for molecule in source["molecules"]:
        mol = Chem.AddHs(Chem.MolFromSmiles(molecule["smiles"]))
        if mol.GetNumAtoms() != len(molecule["elements"]):
            raise SystemExit(f"{molecule['name']}: atom count mismatch")
        conformer = Chem.Conformer(mol.GetNumAtoms())
        for index in range(mol.GetNumAtoms()):
            conformer.SetAtomPosition(
                index, molecule["coords"][3 * index : 3 * index + 3]
            )
        mol.RemoveAllConformers()
        mol.AddConformer(conformer)

        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
        if props is None:
            raise SystemExit(f"{molecule['name']}: no MMFF94 properties")
        force_field = AllChem.MMFFGetMoleculeForceField(mol, props, confId=0)
        if force_field is None:
            raise SystemExit(f"{molecule['name']}: no MMFF94 force field")

        types = [props.GetMMFFAtomType(i) for i in range(mol.GetNumAtoms())]
        if types != molecule["mmff94_types"]:
            raise SystemExit(
                f"{molecule['name']}: RDKit typing {types} differs from the "
                f"oracle fixture {molecule['mmff94_types']}"
            )
        charges = [
            props.GetMMFFPartialCharge(i) for i in range(mol.GetNumAtoms())
        ]
        energy = force_field.CalcEnergy()
        payload["molecules"].append(
            {
                "name": molecule["name"],
                "energy": energy,
                "partial_charges": charges,
                "mmff94_types": types,
            }
        )
        print(f"{molecule['name']:<16} E = {energy: .10f}  ({rdkit.__version__})")

    OUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("wrote", OUT.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
