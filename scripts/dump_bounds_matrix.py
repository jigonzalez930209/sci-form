#!/usr/bin/env python3
"""Dump the distance bounds matrix from RDKit for specific SMILES.

Usage: python3 scripts/dump_bounds_matrix.py "SMILES1" "SMILES2" ...

Outputs JSON with upper and lower bounds matrices for each SMILES.
"""
import json
import sys
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom

def get_bounds_matrix(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    n = mol.GetNumAtoms()
    
    # Get the bounds matrix using RDKit's distance geometry
    # This uses set15 (1-5 bounds) by default
    bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
    
    # bm is a numpy array where upper triangle = upper bounds, lower triangle = lower bounds
    # bm[i][j] for i<j is upper bound, bm[j][i] is lower bound
    result = {
        "smiles": smiles,
        "n_atoms": n,
        "atoms": [],
        "bonds": [],
        "upper": [],
        "lower": [],
    }
    
    for i, atom in enumerate(mol.GetAtoms()):
        result["atoms"].append({
            "element": atom.GetAtomicNum(),
            "hybridization": str(atom.GetHybridization()),
            "formal_charge": atom.GetFormalCharge(),
        })
    
    for bond in mol.GetBonds():
        result["bonds"].append({
            "start": bond.GetBeginAtomIdx(),
            "end": bond.GetEndAtomIdx(),
            "order": str(bond.GetBondType()),
        })
    
    for i in range(n):
        row_upper = []
        row_lower = []
        for j in range(n):
            if i < j:
                row_upper.append(float(bm[i][j]))  # upper bound
                row_lower.append(float(bm[j][i]))  # lower bound
            elif i > j:
                row_upper.append(float(bm[j][i]))  # mirror upper
                row_lower.append(float(bm[i][j]))  # mirror lower
            else:
                row_upper.append(0.0)
                row_lower.append(0.0)
        result["upper"].append(row_upper)
        result["lower"].append(row_lower)
    
    return result

if __name__ == "__main__":
    smiles_list = sys.argv[1:]
    if not smiles_list:
        smiles_list = ["C#CCC1COCC2(C)CC(NC23CCOC3)C1=O"]
    
    results = []
    for smi in smiles_list:
        r = get_bounds_matrix(smi)
        if r is not None:
            results.append(r)
        else:
            print(f"Failed to parse: {smi}", file=sys.stderr)
    
    print(json.dumps(results, indent=2))
