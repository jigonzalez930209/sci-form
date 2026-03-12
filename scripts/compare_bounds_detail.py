"""
Compare our bounds matrix vs RDKit for CC1(O)C2OCC12C
"""
import json, subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDistGeom

smi = 'CC1(O)C2OCC12C'
mol = Chem.AddHs(Chem.MolFromSmiles(smi))
n = mol.GetNumAtoms()

# RDKit bounds
bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)

print(f"RDKit bounds for {smi} (n={n}):")
heavy = [i for i in range(n) if mol.GetAtomWithIdx(i).GetAtomicNum() > 1]
for i in range(n):
    for j in range(i+1, n):
        lb = bm[j][i]  # lower
        ub = bm[i][j]  # upper
        if i in heavy and j in heavy:
            tag = "HEAVY"
        else:
            tag = ""
        print(f"  ({i:2},{j:2}): [{lb:.6f}, {ub:.6f}] {tag}")
