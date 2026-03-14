"""Dump RDKit's sampled distances and initial 4D coords for a molecule.

Matches the RDKit embedding pipeline: pickRandomDistMat → computeInitialCoords.
"""
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom

smiles = sys.argv[1] if len(sys.argv) > 1 else "CC(C#N)(C#N)OC"

mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
n = mol.GetNumAtoms()

# Get bounds matrix
bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)

# Now simulate the RDKit embedding pipeline 
# Use ETKDGv3 with seed=42
ps = AllChem.ETKDGv3()
ps.randomSeed = 42
ps.useRandomCoords = False
ps.maxIterations = 0  # We just want the initial embedding, not optimization

# But we can't easily get intermediate state from RDKit's embed...
# Instead, let's compute manually what RDKit does.

# Step 1: minstd_rand RNG
class MinstdRand:
    A = 48271
    M = 2147483647  # 2^31 - 1
    
    def __init__(self, seed):
        self.state = seed
    
    def next_int(self):
        self.state = (self.state * self.A) % self.M
        return self.state
    
    def next_double(self):
        val = self.next_int()
        return (val - 1.0) / (self.M - 1.0)

rng = MinstdRand(42)

# Step 2: pickRandomDistMat — for i=1..n, for j=0..i
dists = np.zeros((n, n))
print("Sampled distances (first 20 pairs):")
count = 0
for i in range(1, n):
    for j in range(0, i):
        ub = float(bm[j][i]) if j < i else float(bm[i][j])  # upper: bm[min][max]
        lb = float(bm[i][j]) if i > j else float(bm[j][i])  # lower: bm[max][min]
        rval = rng.next_double()
        d = lb + rval * (ub - lb)
        dists[i][j] = d
        dists[j][i] = d
        if count < 20:
            print(f"  d({i},{j}) = {d:.6f}")
        count += 1

# Step 3: computeInitialCoords
# Build squared distance matrix (symmetric packed)
d_size = n * (n + 1) // 2
sq_packed = np.zeros(d_size)
sum_sq_all = 0.0
for i in range(n):
    idx = i * (i + 1) // 2
    for j in range(i + 1):
        d = dists[i][j]
        sq_packed[idx + j] = d * d
        sum_sq_all += d * d
sum_sq_all /= (n * n)

# D0_i computation
d0 = np.zeros(n)
for i in range(n):
    row_sum = 0.0
    for j in range(n):
        if i >= j:
            idx = i * (i + 1) // 2 + j
        else:
            idx = j * (j + 1) // 2 + i
        row_sum += sq_packed[idx]
    d0[i] = row_sum / n - sum_sq_all

# Build metric matrix T
t_packed = np.zeros(d_size)
for i in range(n):
    idx_base = i * (i + 1) // 2
    for j in range(i + 1):
        if i >= j:
            sq_idx = i * (i + 1) // 2 + j
        else:
            sq_idx = j * (j + 1) // 2 + i
        sq_val = sq_packed[sq_idx]
        t_packed[idx_base + j] = 0.5 * (d0[i] + d0[j] - sq_val)

# Power eigensolver
ndim = 4
eigen_seed = int(sum_sq_all * n)

def power_eigen_solver(num_eig, mat, n, seed):
    MAX_ITERATIONS = 1000
    TOLERANCE = 0.001
    HUGE_EIGVAL = 1.0e10
    TINY_EIGVAL = 1.0e-10
    
    eigenvalues = []
    eigenvectors = []
    
    for ei in range(num_eig):
        eig_val = -HUGE_EIGVAL
        seed += ei
        
        # setToRandom with seed
        r = MinstdRand(seed)
        v = np.array([r.next_double() for _ in range(n)])
        norm = np.sqrt(np.sum(v * v))
        if norm > 0:
            v /= norm
        
        converged = False
        for iteration in range(MAX_ITERATIONS):
            # z = mat * v (using packed symmetric)
            z = np.zeros(n)
            for i in range(n):
                idx_base = i * (i + 1) // 2
                for j in range(i + 1):
                    z[i] += mat[idx_base + j] * v[j]
                for j in range(i + 1, n):
                    z[i] += mat[j * (j + 1) // 2 + i] * v[j]
            
            prev_val = eig_val
            eval_id = np.argmax(np.abs(z))
            eig_val = z[eval_id]
            
            if abs(eig_val) < TINY_EIGVAL:
                break
            
            v = z / eig_val
            
            if abs(eig_val - prev_val) < TOLERANCE:
                converged = True
                break
        
        if not converged:
            return None
        
        norm = np.sqrt(np.sum(v * v))
        if norm > 0:
            v /= norm
        
        eigenvalues.append(eig_val)
        eigenvectors.append(v.copy())
        
        # Deflate
        for i in range(n):
            idx_base = i * (i + 1) // 2
            for j in range(i + 1):
                mat[idx_base + j] -= eig_val * v[i] * v[j]
    
    return eigenvalues, eigenvectors

result = power_eigen_solver(ndim, t_packed.copy(), n, eigen_seed)
if result is None:
    print("Eigensolver failed!")
    sys.exit(1)

eigenvalues, eigenvectors = result

# Process eigenvalues
EIGVAL_TOL = 1e-3
eig_sqrt = []
for ev in eigenvalues:
    if ev > EIGVAL_TOL:
        eig_sqrt.append(np.sqrt(ev))
    elif abs(ev) < EIGVAL_TOL:
        eig_sqrt.append(0.0)
    else:
        eig_sqrt.append(ev)  # keep negative marker

# Generate coordinates
print(f"\n4D initial coords (atom x y z w):")
for i in range(n):
    coords = []
    for j in range(ndim):
        if eig_sqrt[j] >= 0:
            coords.append(eig_sqrt[j] * eigenvectors[j][i])
        else:
            coords.append(1.0 - 2.0 * rng.next_double())
    print(f"  {i} {coords[0]:10.6f} {coords[1]:10.6f} {coords[2]:10.6f} {coords[3]:10.6f}")
