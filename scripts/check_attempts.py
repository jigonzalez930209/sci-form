#!/usr/bin/env python3
"""Check which attempt RDKit succeeds on for the failing molecules.
If RDKit fails on some attempts and our code doesn't (or vice versa),
the RNG states would diverge, giving different final conformers."""

from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

# Failing molecules
smiles_list = [
    "CC(C#N)(C#N)OC",      # mol 27
    "CCOC(=O)N(C)C",       # mol 28
    "CCCNC(=O)C=O",        # mol 34
    "C1C(=O)CC1(C=O)O",    # mol 39
    "CNC(=[NH2+])C([O-])=O", # mol 45
    "CC(CCO)OC=O",         # mol 48
    "CCCc1cocn1",           # mol 54
    "OC1(CC=CC1)C=O",      # mol 59
    "CN(C)CC(=O)OC",       # mol 62
    "CC(=NO)CCOC",         # mol 65
    "CC(C)C1CC(O)C1",      # mol 67
    "CCC1CCC1(C)C",        # mol 84
    "CC1(CC1)CCOC",        # mol 94
    "CC1(C)CC(=N)O1",      # mol 98
]

for smi in smiles_list:
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    n = mol.GetNumAtoms()
    
    # Try embedding with trackFailures to see which attempt succeeds
    ps = AllChem.ETKDGv2()
    ps.randomSeed = 42
    ps.trackFailures = True
    cid = AllChem.EmbedMolecule(mol, ps)
    
    if cid >= 0:
        # Check failure counts
        failures = ps.GetFailureCounts()
        total_failures = sum(failures)
        print(f"{smi:30s} n={n:2d} cid={cid:2d} total_failures={total_failures:3d} failures={list(failures)}")
    else:
        print(f"{smi:30s} n={n:2d} FAILED")
