#!/usr/bin/env python3
"""
Compare RDKit's distance geometry output with our Rust implementation
by embedding ethane (CC) with seed=42 and comparing the initial 4D coordinates.

RDKit doesn't directly expose intermediate distance geometry values,
but we can compare the final coordinates as a sanity check.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

mol = Chem.MolFromSmiles("CC")
mol = Chem.AddHs(mol)
n = mol.GetNumAtoms()
print(f"Ethane: {n} atoms")

# Embed with ETKDGv2, seed=42
ps = AllChem.ETKDGv2()
ps.randomSeed = 42
ps.useRandomCoords = False
cid = AllChem.EmbedMolecule(mol, ps)
conf = mol.GetConformer(cid)
print(f"\nRDKit ETKDGv2 coordinates (seed=42):")
for i in range(n):
    p = conf.GetAtomPosition(i)
    print(f"  atom {i} ({mol.GetAtomWithIdx(i).GetSymbol()}): ({p.x:12.8f}, {p.y:12.8f}, {p.z:12.8f})")

# Also try pure DG (no ETKDG) for a cleaner comparison
ps2 = AllChem.EmbedParameters()
ps2.randomSeed = 42
ps2.useExpTorsionAnglePrefs = False
ps2.useBasicKnowledge = False
ps2.useRandomCoords = False
cid2 = AllChem.EmbedMolecule(mol, ps2)
if cid2 >= 0:
    conf2 = mol.GetConformer(cid2)
    print(f"\nRDKit pure DG coordinates (no ETKDG, seed=42):")
    for i in range(n):
        p = conf2.GetAtomPosition(i)
        print(f"  atom {i} ({mol.GetAtomWithIdx(i).GetSymbol()}): ({p.x:12.8f}, {p.y:12.8f}, {p.z:12.8f})")

# Print pairwise distances
print(f"\nRDKit ETKDGv2 pairwise distances (heavy atoms only):")
for i in range(n):
    for j in range(i+1, n):
        if mol.GetAtomWithIdx(i).GetAtomicNum() != 1 and mol.GetAtomWithIdx(j).GetAtomicNum() != 1:
            p1 = conf.GetAtomPosition(i)
            p2 = conf.GetAtomPosition(j)
            d = ((p1.x-p2.x)**2 + (p1.y-p2.y)**2 + (p1.z-p2.z)**2)**0.5
            print(f"  d({i},{j}) = {d:.8f}")

print(f"\nRDKit ETKDGv2 ALL pairwise distances:")
for i in range(n):
    for j in range(i+1, n):
        p1 = conf.GetAtomPosition(i)
        p2 = conf.GetAtomPosition(j)
        d = ((p1.x-p2.x)**2 + (p1.y-p2.y)**2 + (p1.z-p2.z)**2)**0.5
        print(f"  d({i},{j}) = {d:.10f}")
