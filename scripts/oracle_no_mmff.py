"""Generate reference coordinates WITHOUT MMFF optimization.
This tests pure ETKDG output to isolate MMFF's contribution to the gap."""
import json
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
        # NO MMFF optimization — pure ETKDG output
        conf = mol.GetConformer()
        atoms = []
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            atoms.append({
                "element": atom.GetAtomicNum(),
                "x": pos.x, "y": pos.y, "z": pos.z,
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
        results.append({"smiles": smiles, "atoms": atoms, "bonds": bonds})
    return results

if __name__ == "__main__":
    with open("scripts/10k_smiles.smi", 'r') as f:
        smiles = [line.strip() for line in f if line.strip()]
    print(f"Processing {len(smiles)} SMILES (no MMFF)...")
    data = generate_3d_coordinates(smiles)
    print(f"Generated {len(data)} molecules")
    with open("tests/fixtures/reference_coords_no_mmff.json", "w") as f:
        json.dump(data, f)
