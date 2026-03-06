import json
import time
import random
import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

def run_benchmark():
    smiles_file = "scripts/10k_smiles.smi"
    if not os.path.exists(smiles_file):
        print("SMILES file not found.")
        sys.exit(1)
        
    with open(smiles_file, 'r') as f:
        smiles_list = [line.strip() for line in f if line.strip()]
        
    random.seed(42)
    selected_smiles = random.sample(smiles_list, 100)
    
    results = []
    rdkit_total_time = 0.0
    
    for smiles in selected_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
            
        mol = Chem.AddHs(mol)
        
        # We only measure embed time (ETKDG raw generation)
        start_time = time.perf_counter()
        res = AllChem.EmbedMolecule(mol, randomSeed=42)
        end_time = time.perf_counter()
        
        if res != 0:
            continue
            
        rdkit_total_time += (end_time - start_time)
        
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
        
    output_data = {
        "rdkit_time_ms": rdkit_total_time * 1000.0,
        "count": len(results),
        "molecules": results
    }
    
    os.makedirs("tests/fixtures", exist_ok=True)
    with open("tests/fixtures/benchmark_data.json", "w") as f:
        json.dump(output_data, f, indent=2)
        
    print(f"RDKit Benchmark Done. Generated {len(results)} valid conformers in {rdkit_total_time*1000.0:.2f} ms.")

if __name__ == "__main__":
    run_benchmark()
