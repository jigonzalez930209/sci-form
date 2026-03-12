"""Compare RDKit bounds vs sci-form bounds for specific SMILES.
Generates RDKit bounds and dumps them for side-by-side comparison."""
import sys
import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDistGeom

smiles = sys.argv[1] if len(sys.argv) > 1 else "CC(=NO)CCOC"

mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
n = mol.GetNumAtoms()
bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)

# Get topological distance matrix
from rdkit.Chem import rdmolops
dm = rdmolops.GetDistanceMatrix(mol)

print(f"SMILES: {smiles}, {n} atoms")
for i, atom in enumerate(mol.GetAtoms()):
    nbs = [x.GetIdx() for x in atom.GetNeighbors()]
    print(f"  {i}: {atom.GetSymbol()} hyb={atom.GetHybridization()} neighbors={nbs}")

# Print bonds
print("Bonds:")
for bond in mol.GetBonds():
    i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx() 
    print(f"  {i}-{j}: {bond.GetBondTypeAsDouble()}")

# Print bounds grouped by topological distance
for td in [1, 2, 3, 4, 5]:
    print(f"\nTopological distance {td} pairs:")
    for i in range(n):
        for j in range(i+1, n):
            if dm[i][j] == td:
                lower = bm[j][i]
                upper = bm[i][j]
                rng = upper - lower
                print(f"  ({i},{j}) [{lower:.4f}, {upper:.4f}] range={rng:.4f}")

# Remaining pairs
print(f"\nTopological distance >= 6:")
for i in range(n):
    for j in range(i+1, n):
        if dm[i][j] >= 6:
            lower = bm[j][i]
            upper = bm[i][j]
            rng = upper - lower
            print(f"  ({i},{j}) [{lower:.4f}, {upper:.4f}] range={rng:.4f}")
