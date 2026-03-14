"""
Compare random distances picked by our code vs RDKit for mol 27 (CC(C#N)(C#N)OC).
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDistGeom

smiles = "CC(C#N)(C#N)OC"
mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
n = mol.GetNumAtoms()
print(f"Molecule: {smiles}, n={n}")

# Get bounds matrix
bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
print(f"Bounds matrix shape: {bm.shape}")

# Simulate RDKit's MinstdRand
class MinstdRand:
    def __init__(self, seed):
        self.state = seed
        self.A = 48271
        self.M = 2147483647  # 2^31 - 1
    
    def next_int(self):
        self.state = (self.state * self.A) % self.M
        return self.state
    
    def next_double(self):
        val = self.next_int()
        return (val - 1) / (self.M - 2)  # M-2 divisor (our current)

rng = MinstdRand(42)

# Pick random distances exactly matching RDKit's loop
dists = np.zeros((n, n))
for i in range(1, n):
    for j in range(i):
        lb = bm[j][i]  # lower bound at [j][i] for j<i
        ub = bm[i][j]  # upper bound at [i][j] for i>j... no wait
        # RDKit stores: mmat[i][j] = upper for i<j, lower for i>j
        # Python rdDistGeom returns: bm[i][j] = upper for i<j, lower for i>j
        # Actually the Python API returns: bm[i][j] (i<j) = upper, bm[j][i] (j>i) = lower???
        # Let me check
        # RDKit: mmat is upper triangular for upper bounds, lower triangular for lower bounds
        # GetMoleculeBoundsMatrix returns full symmetric with upper in [i][j] where i<j
        # Wait: RDKit docs say bm[i][j] = upper bound distance when i<j
        # When the C++ code does mmat->getUpperBound(i,j), for i<j it returns mmat[i*N+j]
        # When it does mmat->getLowerBound(i,j), for i<j it returns mmat[j*N+i]
        
        # In the pick loop: i=1..n, j=0..i (so j<i always)
        # C++ does: lb = mmat.getLowerBound(j, i) — since j<i, this is mmat[i*N+j] 
        # C++ does: ub = mmat.getUpperBound(j, i) — since j<i, this is mmat[j*N+i]
        # Wait that seems backwards...
        
        # Let me just use the Python bounds matrix as-is
        # bm[i][j] with i<j gives upper bound
        # bm[j][i] with j>i gives lower bound
        # In our loop: i from 1, j from 0..i-1, so j<i always
        # lb should be bm[i][j] and ub should be bm[j][i]??? No...
        
        # Actually in Python rdDistGeom.GetMoleculeBoundsMatrix:
        # For i<j: upper bound is bm[i][j], lower bound is bm[j][i]
        # In our loop i>j, so:
        # upper bound = bm[j][i] (since j<i)
        # lower bound = bm[i][j] (since i>j)
        
        ub_val = bm[j][i]  # j<i -> upper bound
        lb_val = bm[i][j]  # i>j -> lower bound
        
        r = rng.next_double()
        d = lb_val + r * (ub_val - lb_val)
        dists[i][j] = d
        dists[j][i] = d

# Print first 10 distances
print("\nFirst 20 random distances (i, j, lb, ub, rand, dist):")
count = 0
rng2 = MinstdRand(42)
for i in range(1, n):
    for j in range(i):
        ub_val = bm[j][i]
        lb_val = bm[i][j]
        r = rng2.next_double()
        d = lb_val + r * (ub_val - lb_val)
        if count < 20:
            print(f"  ({i:2d},{j:2d}): lb={lb_val:.6f}, ub={ub_val:.6f}, r={r:.10f}, d={d:.6f}")
        count += 1

print(f"\nTotal distance picks: {count}")
print(f"RNG state after all picks: {rng.state}")
