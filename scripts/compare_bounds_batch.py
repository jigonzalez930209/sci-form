"""Compare RDKit bounds matrices with sci-form bounds for test molecules.
Reads the same reference data and generates comparison statistics."""
import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom

# Load reference data  
with open("tests/fixtures/reference_coords.json") as f:
    all_mols = json.load(f)

# Same shuffling as test
import random
rng = random.Random()
rng.seed(123)
# Python's shuffle uses a different algo than Rust's, so let's match by using numpy
import numpy as np
np_rng = np.random.RandomState(123)
indices = np_rng.permutation(len(all_mols))[:100]
molecules = [all_mols[i] for i in indices]

# For each molecule, generate RDKit bounds and dump stats
results = []
for idx, mol_data in enumerate(molecules):
    smiles = mol_data["smiles"]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        continue
    mol = Chem.AddHs(mol)
    n = mol.GetNumAtoms()
    
    try:
        bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
    except:
        continue
    
    # Count tight pairs (1-2, 1-3 with range < 0.2)
    tight_pairs = 0
    medium_pairs = 0  
    wide_pairs = 0
    for i in range(n):
        for j in range(i+1, n):
            rng = bm[i][j] - bm[j][i]
            if rng < 0.2:
                tight_pairs += 1
            elif rng < 1.0:
                medium_pairs += 1
            else:
                wide_pairs += 1
    
    results.append({
        "idx": idx,
        "smiles": smiles,
        "n_atoms": n,
        "tight": tight_pairs,
        "medium": medium_pairs,
        "wide": wide_pairs,
    })

# Print bounds summary for first 10 molecules
print("Mol# | SMILES | Atoms | Tight | Medium | Wide")
for r in results[:20]:
    print(f"  {r['idx']:3d} | {r['smiles'][:30]:30s} | {r['n_atoms']:5d} | {r['tight']:5d} | {r['medium']:6d} | {r['wide']:4d}")

# Detailed comparison for specific molecules we know are problematic
problem_smiles = [
    "N=COC1=COC=N1",       # Mol 91
    "CC(=NO)CCOC",          # Mol 65
    "CC(=O)N(C)CC=O",       # Mol 92
    "CCOC(=O)N(C)C",        # Mol 28
    "CC(C)C1CC1C",          # Mol 85
]

for ps in problem_smiles:
    mol = Chem.MolFromSmiles(ps)
    mol = Chem.AddHs(mol)
    n = mol.GetNumAtoms()
    bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
    
    print(f"\n=== {ps} ({n} atoms) ===")
    # Print atom info
    for i, atom in enumerate(mol.GetAtoms()):
        print(f"  {i}: {atom.GetSymbol()} hyb={atom.GetHybridization()}")
    
    # Print all tight bounds (1-2 and 1-3 pairs) 
    print("  Tight bounds (range < 0.15):")
    for i in range(n):
        for j in range(i+1, n):
            lower = bm[j][i]
            upper = bm[i][j]
            rng = upper - lower
            if rng < 0.15:
                elem_i = mol.GetAtomWithIdx(i).GetAtomicNum()
                elem_j = mol.GetAtomWithIdx(j).GetAtomicNum()
                print(f"    ({i},{j}) e({elem_i},{elem_j}) [{lower:.4f}, {upper:.4f}]")
