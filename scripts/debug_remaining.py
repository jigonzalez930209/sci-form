"""
Compare bounds matrices and coordinates for the 2 remaining failing molecules.
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
import json

smis = ['CC1(O)CN2CC12', 'CC1(O)C2OCC12C']

for smi in smis:
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    n = mol.GetNumAtoms()
    
    # Get bounds matrix
    bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
    
    # Embed with ETKDGv3, seed=42
    ps = rdDistGeom.ETKDGv3()
    ps.randomSeed = 42
    ps.useRandomCoords = False
    ps.numThreads = 1
    ps.pruneRmsThresh = -1.0
    ps.useSmallRingTorsions = True
    ps.useMacrocycleTorsions = True
    ps.useBasicKnowledge = True
    ps.ETversion = 2
    
    cid = AllChem.EmbedMolecule(mol, ps)
    if cid < 0:
        print(f"\n{smi}: EMBEDDING FAILED")
        continue
    
    conf = mol.GetConformer()
    coords = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, 
                         conf.GetAtomPosition(i).z] for i in range(n)])
    
    print(f"\n=== {smi} (n={n}) ===")
    
    # Print ring info
    ri = mol.GetRingInfo()
    print(f"Rings: {[list(r) for r in ri.AtomRings()]}")
    
    # Print atoms
    for i in range(n):
        a = mol.GetAtomWithIdx(i)
        print(f"  Atom {i:2d}: {a.GetSymbol():2s} hyb={str(a.GetHybridization()):5s} deg={a.GetDegree()} "
              f"coords=({coords[i][0]:8.4f}, {coords[i][1]:8.4f}, {coords[i][2]:8.4f})")
    
    # Print bounds for heavy atoms only
    heavy = [i for i in range(n) if mol.GetAtomWithIdx(i).GetAtomicNum() > 1]
    print(f"\nBounds matrix (heavy atoms only, {len(heavy)} atoms):")
    print("     ", "".join(f"{i:8d}" for i in heavy))
    for i in heavy:
        row = ""
        for j in heavy:
            if i == j:
                row += "    ----"
            elif j > i:
                row += f"{bm[i][j]:8.3f}"  # upper bound
            else:
                row += f"{bm[i][j]:8.3f}"  # lower bound
        print(f"  {i:2d}: {row}")
    
    # Save bounds for comparison with Rust
    bounds_dict = {}
    for i in range(n):
        for j in range(i+1, n):
            bounds_dict[f"{i},{j}"] = {"lower": float(bm[j][i]), "upper": float(bm[i][j])}
    
    with open(f'/tmp/bounds_{smi.replace("(", "").replace(")", "").replace("/", "")}.json', 'w') as f:
        json.dump(bounds_dict, f)
    print(f"Bounds saved to /tmp/bounds_*.json")
