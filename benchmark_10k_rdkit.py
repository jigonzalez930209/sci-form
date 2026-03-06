import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import time
import json
import os

def run_10k_rdkit_benchmark():
    smi_file = "scripts/10k_smiles.smi"
    if not os.path.exists(smi_file):
        print(f"Error: {smi_file} not found")
        return

    with open(smi_file, "r") as f:
        smiles_list = [line.strip() for line in f if line.strip()]

    print(f"Starting RDKit benchmark for {len(smiles_list)} molecules...")
    
    count = 0
    start_time = time.perf_counter()
    
    # We will save coordinates for a subset (first 500) to keep reference file manageable
    reference_data = []
    
    for i, smi in enumerate(smiles_list):
        try:
            m = Chem.MolFromSmiles(smi)
            if m is None: 
                print(f"Failed to parse SMILES: {smi}")
                continue
            m = Chem.AddHs(m)
            
            # Use ETKDGv3
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            res = AllChem.EmbedMolecule(m, params)
            
            if res == -1:
                # Try fallback if ETKDG fails for very small molecules
                res = AllChem.EmbedMolecule(m, AllChem.ETKDG())
            
            if res != -1:
                # try:
                #     AllChem.MMFFOptimizeMolecule(m)
                # except:
                #     pass 
                
                if count < 500:
                    try:
                        conf = m.GetConformer()
                        atoms = []
                        for k, atom in enumerate(m.GetAtoms()):
                            pos = conf.GetAtomPosition(k)
                            atoms.append({
                                "element": atom.GetAtomicNum(),
                                "x": pos.x,
                                "y": pos.y,
                                "z": pos.z,
                                "formal_charge": atom.GetFormalCharge(),
                                "hybridization": str(atom.GetHybridization())
                            })
                        
                        bonds = []
                        for bond in m.GetBonds():
                            bonds.append({
                                "start": bond.GetBeginAtomIdx(),
                                "end": bond.GetEndAtomIdx(),
                                "order": str(bond.GetBondType())
                            })
                            
                        reference_data.append({
                            "smiles": smi,
                            "atoms": atoms,
                            "bonds": bonds
                        })
                    except Exception as e:
                        if i < 20:
                            print(f"Error extracting coords for {smi}: {e}")
                        continue
                
                count += 1
                if count % 1000 == 0:
                    print(f"Processed {count} molecules...")
        except Exception as e:
            if i < 10: # Print only first 10 errors to avoid spam
                print(f"Error processing {smi}: {e}")
            continue

    end_time = time.perf_counter()
    total_ms = (end_time - start_time) * 1000.0
    
    print(f"\n=== RDKit 10k Benchmark Results ===")
    print(f"Molecules Successfully Embedded: {count}")
    print(f"Total Time: {total_ms:.2f} ms")
    print(f"Average Time: {total_ms/count:.2f} ms/mol")
    
    # Save reference data
    os.makedirs("tests/fixtures", exist_ok=True)
    with open("tests/fixtures/rdkit_10k_reference.json", "w") as f:
        json.dump(reference_data, f)
    print(f"Saved reference coordinates for first 500 molecules to tests/fixtures/rdkit_10k_reference.json")

if __name__ == "__main__":
    run_10k_rdkit_benchmark()
