"""
Dump intermediate coordinates from RDKit for specific molecules.
Compares distances, eigenvalues, initial coords, and post-4D BFGS coords.
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
from rdkit.Geometry import Point3D
import json

def dump_rdkit_intermediates(smiles, seed=42):
    """Use RDKit's Python API to get as much intermediate data as possible."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    n = mol.GetNumAtoms()
    
    # Get bounds matrix
    bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
    
    # Run embedding with specific seed to get final coords
    ps = rdDistGeom.ETKDGv3()
    ps.randomSeed = seed
    ps.useRandomCoords = False
    ps.numThreads = 1
    ps.pruneRmsThresh = -1.0
    ps.useSmallRingTorsions = True
    ps.useMacrocycleTorsions = True
    ps.useBasicKnowledge = True
    ps.ETversion = 2
    
    AllChem.EmbedMolecule(mol, ps)
    
    conf = mol.GetConformer()
    coords_3d = []
    for i in range(n):
        p = conf.GetAtomPosition(i)
        coords_3d.append([p.x, p.y, p.z])
    
    print(f"\n=== {smiles} (n={n}) ===")
    print(f"Final 3D coords from RDKit:")
    for i, c in enumerate(coords_3d):
        print(f"  atom {i:2d}: ({c[0]:10.6f}, {c[1]:10.6f}, {c[2]:10.6f})")
    
    return {
        'smiles': smiles,
        'n': n,
        'bounds': bm.tolist(),
        'coords_3d': coords_3d,
    }

# Test molecules that fail
failing = [
    "CC(C#N)(C#N)OC",     # mol 27
    "CCCc1cocn1",          # mol 54  
    "CC1(CC1)CCOC",        # mol 94
]

for smi in failing:
    dump_rdkit_intermediates(smi)
