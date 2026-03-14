import json
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

def generate_3d_coordinates(smiles_list):
    results = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        mol = Chem.AddHs(mol)
        res = AllChem.EmbedMolecule(mol, randomSeed=42)
        if res != 0:
            continue
        # NO MMFF optimization - raw ETKDG output
        conf = mol.GetConformer()
        atoms = []
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            atoms.append({
                "element": atom.GetAtomicNum(),
                "x": pos.x,
                "y": pos.y,
                "z": pos.z,
                "formal_charge": atom.GetFormalCharge(),
                "hybridization": str(atom.GetHybridization())
            })
        bonds = []
        for bond in mol.GetBonds():
            bonds.append({
                "start": bond.GetBeginAtomIdx(),
                "end": bond.GetEndAtomIdx(),
                "order": str(bond.GetBondType())
            })
        results.append({
            "smiles": smiles,
            "atoms": atoms,
            "bonds": bonds
        })
    return results

if __name__ == "__main__":
    input_file = sys.argv[1]
    with open(input_file, 'r') as f:
        smiles = [line.strip() for line in f if line.strip()]
    print(f"Processing {len(smiles)} SMILES strings...")
    data = generate_3d_coordinates(smiles)
    output_file = "tests/fixtures/reference_coords_noff.json"
    import os
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Generated {len(data)} molecules")
