"""Compare RDKit bounds vs sci-form bounds for a specific molecule."""
import sys
import json
import subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom

smiles = sys.argv[1] if len(sys.argv) > 1 else "N=COC1=COC=N1"

mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
n = mol.GetNumAtoms()
bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)

# Also embed to get reference conformer
AllChem.EmbedMolecule(mol, randomSeed=42)
conf = mol.GetConformer()

print(f"SMILES: {smiles}, {n} atoms")
print(f"\nAtoms:")
for i, atom in enumerate(mol.GetAtoms()):
    pos = conf.GetAtomPosition(i)
    print(f"  {i}: {atom.GetSymbol()} hyb={atom.GetHybridization()}")

# Print bounds with comparison tracking
print(f"\nRDKit bounds (pairs with range < 0.2 = tightly constrained):")
for i in range(n):
    for j in range(i+1, n):
        lower = bm[j][i]
        upper = bm[i][j]
        # Compute distance in conformer  
        pi = conf.GetAtomPosition(i)
        pj = conf.GetAtomPosition(j)
        dist = pi.Distance(pj)
        rng = upper - lower
        if rng < 0.2:
            viol = max(0, lower - dist, dist - upper) 
            mark = " ***" if viol > 0.01 else ""
            print(f"  ({i},{j}) [{lower:.3f}, {upper:.3f}] ref={dist:.3f}{mark}")

# Print all 1-5+ pair bounds (wide range) 
print(f"\nWide-range pairs (range > 0.5):")
for i in range(n):
    for j in range(i+1, n):
        lower = bm[j][i]
        upper = bm[i][j]
        rng = upper - lower
        if rng > 0.5:
            print(f"  ({i},{j}) [{lower:.3f}, {upper:.3f}] range={rng:.3f}")
