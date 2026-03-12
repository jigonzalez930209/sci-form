"""Compare raw ETKDG output vs MMFF94-optimized output to measure the MMFF effect."""
import json
import numpy as np

with open("tests/fixtures/reference_coords.json") as f:
    mmff = json.load(f)

with open("tests/fixtures/reference_coords_noff.json") as f:
    raw = json.load(f)

# Index by smiles
mmff_idx = {m["smiles"]: m for m in mmff}
raw_idx = {m["smiles"]: m for m in raw}

# Use the same 100 molecules from the test (seed 123)
from random import Random
rng = Random(123)
smiles_list = list(set(mmff_idx.keys()) & set(raw_idx.keys()))
rng.shuffle(smiles_list)
smiles_list = smiles_list[:100]

rmsds = []
for smi in smiles_list:
    m = mmff_idx[smi]
    r = raw_idx[smi]
    n = len(m["atoms"])
    if len(r["atoms"]) != n:
        continue
    
    diffs_sq = []
    for i in range(n):
        for j in range(i+1, n):
            dx_m = m["atoms"][i]["x"] - m["atoms"][j]["x"]
            dy_m = m["atoms"][i]["y"] - m["atoms"][j]["y"]
            dz_m = m["atoms"][i]["z"] - m["atoms"][j]["z"]
            d_m = (dx_m**2 + dy_m**2 + dz_m**2)**0.5
            
            dx_r = r["atoms"][i]["x"] - r["atoms"][j]["x"]
            dy_r = r["atoms"][i]["y"] - r["atoms"][j]["y"]
            dz_r = r["atoms"][i]["z"] - r["atoms"][j]["z"]
            d_r = (dx_r**2 + dy_r**2 + dz_r**2)**0.5
            
            diffs_sq.append((d_m - d_r)**2)
    
    rmsd = (sum(diffs_sq) / len(diffs_sq))**0.5
    rmsds.append((rmsd, smi))

rmsds.sort(key=lambda x: x[0], reverse=True)
print(f"MMFF94 effect: Average pairwise distance RMSD between raw ETKDG and MMFF-optimized")
print(f"Average: {np.mean([r for r,s in rmsds]):.3f} Å")
print(f"Max: {rmsds[0][0]:.3f} Å")
print(f"Above 0.5: {sum(1 for r,s in rmsds if r > 0.5)}")
print(f"\nWorst 10:")
for rmsd, smi in rmsds[:10]:
    print(f"  {rmsd:.3f} {smi}")
print(f"\nBest 10:")
for rmsd, smi in rmsds[-10:]:
    print(f"  {rmsd:.3f} {smi}")
