"""
Compare bounds matrices between RDKit and our code for the 2 failing molecules.
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDistGeom
import json

smis = ['CC1(O)CN2CC12', 'CC1(O)C2OCC12C']

for smi in smis:
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    n = mol.GetNumAtoms()
    
    # Get RDKit bounds matrix
    bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
    
    # Print bounds as JSON for comparison
    bounds = {}
    for i in range(n):
        for j in range(i+1, n):
            # bm[i][j] = upper bound, bm[j][i] = lower bound (for i<j)
            bounds[f"{i},{j}"] = [float(bm[j][i]), float(bm[i][j])]
    
    with open(f'/tmp/rdkit_bounds_{n}.json', 'w') as f:
        json.dump(bounds, f, indent=2)
    
    # Print heavy-atom bounds
    heavy = [i for i in range(n) if mol.GetAtomWithIdx(i).GetAtomicNum() > 1]
    print(f"\n=== {smi} (n={n}) ===")
    print(f"Heavy atoms: {heavy}")
    for i in heavy:
        for j in heavy:
            if j > i:
                lb = bm[j][i]
                ub = bm[i][j]
                print(f"  ({i:2d},{j:2d}): [{lb:.4f}, {ub:.4f}] (range={ub-lb:.4f})")
