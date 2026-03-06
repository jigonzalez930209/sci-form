import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import json
import time
import numpy as np

def run_benchmark():
    with open("tests/fixtures/benchmark_data.json", "r") as f:
        data = json.load(f)
    
    molecules = data['molecules']
    count = len(molecules)
    
    # 1. Measure RDKit RAW ETKDG (Matches what's in the JSON mostly)
    start_raw = time.perf_counter()
    for item in molecules:
        m = Chem.MolFromSmiles(item['smiles'])
        m = Chem.AddHs(m)
        AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
    end_raw = time.perf_counter()
    raw_time_ms = (end_raw - start_raw) * 1000.0
    
    # 2. Measure RDKit FULL (ETKDG + MMFF Minimization)
    start_full = time.perf_counter()
    for item in molecules:
        m = Chem.MolFromSmiles(item['smiles'])
        m = Chem.AddHs(m)
        AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
        # MMFF94 is the standard RDKit minimizer
        AllChem.MMFFOptimizeMolecule(m)
    end_full = time.perf_counter()
    full_time_ms = (end_full - start_full) * 1000.0
    
    print(f"=== RDKit (Python) Performance Baseline ===")
    print(f"Molecules Processed: {count}")
    print(f"RDKit Raw (ETKDG) Time: {raw_time_ms:.2f} ms ({raw_time_ms/count:.2f} ms/mol)")
    print(f"RDKit Full (ETKDG + MMFF) Time: {full_time_ms:.2f} ms ({full_time_ms/count:.2f} ms/mol)")

if __name__ == "__main__":
    run_benchmark()
