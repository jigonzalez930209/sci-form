"""Trace RDKit's ETKDG embedding pipeline for a specific molecule.

Dumps coordinates after each stage:
1. Initial 4D coords (from pickRandomDistMat + computeInitialCoords)
2. After 1st bounds FF minimization (chiral_w=1.0, 4d_w=0.1)
3. After 2nd bounds FF minimization (chiral_w=0.2, 4d_w=1.0)
4. After 3D projection
5. After ETKDG 3D minimization
"""
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom

smiles = sys.argv[1] if len(sys.argv) > 1 else "CC(C#N)(C#N)OC"

mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
n = mol.GetNumAtoms()

print(f"SMILES: {smiles}")
print(f"Atoms: {n}")
for i, atom in enumerate(mol.GetAtoms()):
    print(f"  {i}: {atom.GetSymbol()} ({atom.GetAtomicNum()}) hyb={atom.GetHybridization()}")

# We can't easily get intermediate coords from RDKit
# But we can get the final result and see if it differs
ps = AllChem.ETKDGv3()
ps.randomSeed = 42
ps.useBasicKnowledge = True
ps.useExpTorsionAnglePrefs = False  # no CSD torsions (matching our test)
ps.numThreads = 1

res = AllChem.EmbedMolecule(mol, ps)
if res != 0:
    print("Embedding failed!")
    sys.exit(1)

conf = mol.GetConformer()
print("\n=== FINAL COORDS (ETKDGv3, no CSD torsion prefs) ===")
for i in range(n):
    pos = conf.GetAtomPosition(i)
    print(f"  {i} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}")

# Also try with full ETKDG (with CSD torsions) for comparison
mol2 = Chem.MolFromSmiles(smiles)
mol2 = Chem.AddHs(mol2)
ps2 = AllChem.ETKDGv3()
ps2.randomSeed = 42
ps2.numThreads = 1
res2 = AllChem.EmbedMolecule(mol2, ps2)
if res2 == 0:
    conf2 = mol2.GetConformer()
    print("\n=== FINAL COORDS (ETKDGv3, WITH CSD torsion prefs) ===")
    for i in range(n):
        pos = conf2.GetAtomPosition(i)
        print(f"  {i} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}")
