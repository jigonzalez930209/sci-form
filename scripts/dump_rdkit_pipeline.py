#!/usr/bin/env python3
"""Dump RDKit's full ETKDG pipeline intermediate state for a single molecule."""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom

SMILES = "C#CCOC(C)CC1CC2C3CCC(C)C(O)(C3)C2O1"

mol = Chem.MolFromSmiles(SMILES)
mol = Chem.AddHs(mol)

# Get bounds matrix
bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
n = mol.GetNumAtoms()
print(f"Molecule: {SMILES}")
print(f"Atoms: {n}")

# Save bounds
np.save("/tmp/rdkit_bounds_fail.npy", bm)
print(f"Bounds saved to /tmp/rdkit_bounds_fail.npy")
print(f"  UB sum: {sum(bm[i][j] for i in range(n) for j in range(i+1, n)):.10f}")
print(f"  LB sum: {sum(bm[j][i] for i in range(n) for j in range(i+1, n)):.10f}")

# Generate conformer with seed=42
params = AllChem.ETKDGv3()
params.randomSeed = 42
params.useRandomCoords = False
params.numThreads = 1

cid = AllChem.EmbedMolecule(mol, params)
if cid >= 0:
    conf = mol.GetConformer(cid)
    coords = np.array([[conf.GetAtomPosition(i).x,
                         conf.GetAtomPosition(i).y,
                         conf.GetAtomPosition(i).z] for i in range(n)])
    np.save("/tmp/rdkit_coords_fail.npy", coords)
    print(f"\nFinal coordinates saved to /tmp/rdkit_coords_fail.npy")
    print(f"  First 5 atoms:")
    for i in range(min(5, n)):
        print(f"    {i}: ({coords[i][0]:.6f}, {coords[i][1]:.6f}, {coords[i][2]:.6f})")
else:
    print("Embedding FAILED!")
