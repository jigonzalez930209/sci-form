#!/usr/bin/env python3
"""
Determine which attempt RDKit uses for a specific molecule.
Uses binary search on maxIterations to find the first successful attempt.
"""
import json, sys
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
import numpy as np

def find_rdkit_attempt(smiles, max_check=200):
    """Find which attempt number RDKit succeeds on."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    
    for max_iter in range(1, max_check + 1):
        params = rdDistGeom.ETKDGv3()
        params.randomSeed = 42
        params.maxIterations = max_iter
        res = AllChem.EmbedMolecule(mol, params)
        if res == 0:
            return max_iter - 1  # attempt is 0-indexed
    return -1

def main():
    with open('tests/fixtures/gdb20_reference.json') as f:
        data = json.load(f)
    
    # Check specific molecules
    indices = [int(x) for x in sys.argv[1:]] if len(sys.argv) > 1 else list(range(20))
    
    for idx in indices:
        mol_data = data[idx]
        smiles = mol_data['smiles']
        attempt = find_rdkit_attempt(smiles)
        print(f"mol[{idx}]: smiles={smiles}, rdkit_attempt={attempt}")

if __name__ == '__main__':
    main()
