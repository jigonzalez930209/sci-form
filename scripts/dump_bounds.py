"""Dump RDKit's bounds matrix for specific SMILES molecules.

Usage:
  python dump_bounds.py SMILES       # single molecule, text output
  python dump_bounds.py --json S1 S2 # JSON output for multiple molecules 
"""
import sys
import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom


def get_bounds(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, None
    mol = Chem.AddHs(mol)
    bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
    return bm, mol.GetNumAtoms(), mol


if "--json" in sys.argv:
    sys.argv.remove("--json")
    smiles_list = sys.argv[1:]
    results = {}
    for smiles in smiles_list:
        bm, n, mol = get_bounds(smiles)
        if bm is None:
            continue
        # Store as flat list: bounds[i*n+j] for all i,j
        flat = []
        for i in range(n):
            for j in range(n):
                flat.append(round(float(bm[i][j]), 6))
        results[smiles] = {"n": n, "bounds": flat}
    print(json.dumps(results))
    sys.exit(0)

# Single molecule text output
smiles = sys.argv[1] if len(sys.argv) > 1 else "N=COC1=COC=N1"
bm, n, mol = get_bounds(smiles)

print(f"SMILES: {smiles}")
print(f"Atoms: {n}")
for i, atom in enumerate(mol.GetAtoms()):
    print(f"  {i}: {atom.GetSymbol()} ({atom.GetAtomicNum()}) hyb={atom.GetHybridization()}")

print(f"\nBonds:")
for bond in mol.GetBonds():
    i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    print(f"  {i}-{j}: {bond.GetBondTypeAsDouble()}")

print(f"\nBounds matrix (lower triangle = lower, upper triangle = upper):")
print(f"Format: (i,j) lower upper")
for i in range(n):
    for j in range(i+1, n):
        lower = bm[j][i]  # lower bound
        upper = bm[i][j]  # upper bound
        if upper < 999:
            print(f"  ({i},{j}) [{lower:.3f}, {upper:.3f}]  range={upper-lower:.3f}")
        else:
            print(f"  ({i},{j}) [{lower:.3f}, INF]")
