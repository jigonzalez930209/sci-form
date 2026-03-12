#!/usr/bin/env python3
"""Extract intermediate distance geometry values from RDKit for mol 84 (CCC1CCC1(C)C).
We need initial random distances, metric matrix, eigenvalues, eigenvectors, and initial coords."""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
import numpy as np

smiles = "CCC1CCC1(C)C"
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
n = mol.GetNumAtoms()
print(f"Molecule: {smiles}, {n} atoms")

# Get bounds matrix
bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
print(f"\nBounds matrix ({bm.shape}):")
print(f"  bm[0,1] = {bm[0,1]:.10f} (upper)")
print(f"  bm[1,0] = {bm[1,0]:.10f} (lower)")
print(f"  bm[0,2] = {bm[0,2]:.10f} (upper)")
print(f"  bm[2,0] = {bm[2,0]:.10f} (lower)")

# Try embedding
ps = AllChem.ETKDGv2()
ps.randomSeed = 42
ps.maxIterations = 1  # Single attempt
cid = AllChem.EmbedMolecule(mol, ps)

if cid >= 0:
    conf = mol.GetConformer(cid)
    print(f"\nFinal coordinates:")
    for i in range(n):
        p = conf.GetAtomPosition(i)
        print(f"  atom {i:2d} ({mol.GetAtomWithIdx(i).GetSymbol():2s}): ({p.x:12.8f}, {p.y:12.8f}, {p.z:12.8f})")
    
    # Print all pairwise distances
    print(f"\nPairwise distances:")
    for i in range(n):
        for j in range(i+1, n):
            p1 = conf.GetAtomPosition(i)
            p2 = conf.GetAtomPosition(j)
            d = ((p1.x-p2.x)**2 + (p1.y-p2.y)**2 + (p1.z-p2.z)**2)**0.5
            print(f"  d({i:2d},{j:2d}) = {d:.10f}")
else:
    print("Embedding failed with maxIterations=1, trying with default")
    ps2 = AllChem.ETKDGv2()
    ps2.randomSeed = 42
    cid2 = AllChem.EmbedMolecule(mol, ps2)
    if cid2 >= 0:
        conf = mol.GetConformer(cid2)
        print(f"\nFinal coordinates (with retries):")
        for i in range(n):
            p = conf.GetAtomPosition(i)
            print(f"  atom {i:2d} ({mol.GetAtomWithIdx(i).GetSymbol():2s}): ({p.x:12.8f}, {p.y:12.8f}, {p.z:12.8f})")
