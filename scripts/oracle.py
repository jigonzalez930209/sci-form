import json
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

def generate_3d_coordinates(smiles_list):
    results = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Failed to parse SMILES: {smiles}")
            continue
        
        # Add hydrogens (standard for 3D embedding)
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates using standard ETKDGv3
        res = AllChem.EmbedMolecule(mol, randomSeed=42)
        if res != 0:
            print(f"Failed to embed SMILES: {smiles}")
            continue
            
        # Force field optimization (MMFF94)
        AllChem.MMFFOptimizeMolecule(mol)
        
        # Extract atoms and coordinates
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
            
        # Extract bonds
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
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        with open(input_file, 'r') as f:
            smiles = [line.strip() for line in f if line.strip()]
    else:
        # Default mini test set
        smiles = [
            "C",          # Methane
            "CCO",        # Ethanol
            "c1ccccc1",   # Benzene
            "CC(=O)O",    # Acetic Acid
            "N#C",        # Hydrogen cyanide
        ]
        
    print(f"Processing {len(smiles)} SMILES strings...")
    data = generate_3d_coordinates(smiles)
    
    output_file = "tests/fixtures/reference_coords.json"
    import os
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    with open(output_file, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Saved {len(data)} reference structures to {output_file}")
