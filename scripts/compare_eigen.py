#!/usr/bin/env python3
"""Compute metric matrix and eigenvalues from RDKit bounds to compare with Rust."""
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDistGeom

SMILES = "C#CCOC(C)CC1CC2C3CCC(C)C(O)(C3)C2O1"
mol = Chem.MolFromSmiles(SMILES)
mol = Chem.AddHs(mol)

bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
n = mol.GetNumAtoms()
print(f"Atoms: {n}")

# Reproduce RDKit's MinstdRand
class MinstdRand:
    def __init__(self, seed):
        self.state = seed if seed != 0 else 1

    def next_int(self):
        self.state = (self.state * 16807) % 2147483647
        return self.state

    def next_double(self):
        return self.next_int() / 2147483646.0

rng = MinstdRand(42)

# Pick random distances (RDKit loop order: i=1..n, j=0..i)
dists = np.zeros((n, n))
for i in range(1, n):
    for j in range(i):
        lb = bm[i][j]  # lower bound
        ub = bm[j][i]  # upper bound (stored in upper triangle)
        # Wait, RDKit stores UB in bm[i][j] where i<j, LB in bm[j][i]
        # So for pair (i,j) with i>j: UB = bm[j][i], LB = bm[i][j]
        lb = bm[i][j]
        ub = bm[j][i]
        r = rng.next_double()
        d = lb + r * (ub - lb)
        dists[i][j] = d
        dists[j][i] = d

print("=== Distances (first 10 pairs) ===")
for i in range(min(10, n)):
    for j in range(i+1, min(i+3, n)):
        print(f"  d({i},{j}) = {dists[i][j]:.10f}")

# Squared distance matrix
dsq = dists ** 2

# Metric matrix (RDKit's compute_metric_matrix)
# t[i][j] = 0.5 * (dsq[0][i] + dsq[0][j] - dsq[i][j])
# Packed symmetric format: upper triangle row-major
packed_size = n * (n + 1) // 2
t_packed = np.zeros(packed_size)

sum_sq_all = 0.0
for i in range(n):
    for j in range(n):
        sum_sq_all += dsq[i][j]

print(f"\nsum_sq_all = {sum_sq_all:.10f}")
idx = 0
for i in range(n):
    for j in range(i, n):
        d_val = 0.5 * (dsq[0][i] + dsq[0][j] - dsq[i][j])
        t_packed[idx] = d_val
        idx += 1

# Save for comparison
np.save("/tmp/rdkit_metric_packed.npy", t_packed)
np.save("/tmp/rdkit_distances.npy", dists)
print("Metric matrix and distances saved")

# Print first few metric matrix values
print("\n=== Metric matrix (first 10 packed values) ===")
for i in range(10):
    print(f"  t[{i}] = {t_packed[i]:.12f}")

# Compute eigenvalues using numpy for comparison
# Convert packed to full matrix
full = np.zeros((n, n))
idx = 0
for i in range(n):
    for j in range(i, n):
        full[i][j] = t_packed[idx]
        full[j][i] = t_packed[idx]
        idx += 1

eigenvalues, eigenvectors = np.linalg.eigh(full)
# Sort descending
idx_sort = np.argsort(-eigenvalues)
eigenvalues = eigenvalues[idx_sort]
eigenvectors = eigenvectors[:, idx_sort]

print(f"\n=== Top 5 Eigenvalues (numpy) ===")
for i in range(5):
    print(f"  λ{i} = {eigenvalues[i]:.10f}")

# Now simulate the power iteration to get the SAME eigenvalues/vectors as RDKit
# This is the key comparison point
# RDKit's power iteration uses:
# - eigen_seed = int(sum_sq_all * n)
# - MAX_ITERATIONS = 1000
# - TOLERANCE = 0.001
# - TINY_EIGVAL = 1e-10

eigen_seed = int(sum_sq_all * n)
print(f"\neigen_seed = {eigen_seed}")

# Power iteration
def power_eigen(t_packed, n, ndim, seed):
    """Simulate RDKit's powerEigenSolver."""
    local_rng = MinstdRand(seed)
    size = n * (n + 1) // 2

    def sym_mat_vec_multiply(packed, vec, n):
        """Multiply packed symmetric matrix by vector."""
        result = np.zeros(n)
        idx = 0
        for i in range(n):
            # Diagonal
            result[i] += packed[idx] * vec[i]
            idx += 1
            # Off-diagonal
            for j in range(i+1, n):
                result[i] += packed[idx] * vec[j]
                result[j] += packed[idx] * vec[i]
                idx += 1
        return result

    eigenvalues = []
    eigenvectors = []

    for dim in range(ndim):
        # Initialize random vector
        v = np.array([local_rng.next_double() for _ in range(n)])
        norm = np.sqrt(np.dot(v, v))
        v /= norm

        for iteration in range(1000):
            w = sym_mat_vec_multiply(t_packed, v, n)
            eigenvalue = np.dot(v, w)
            norm = np.sqrt(np.dot(w, w))
            if norm < 1e-10:
                break
            w /= norm

            diff = np.sqrt(np.dot(w - v, w - v))
            v = w.copy()
            if diff < 0.001:
                break

        eigenvalue = np.dot(v, sym_mat_vec_multiply(t_packed, v, n))
        eigenvalues.append(eigenvalue)
        eigenvectors.append(v.copy())

        # Deflate
        idx = 0
        for i in range(n):
            for j in range(i, n):
                t_packed[idx] -= eigenvalue * v[i] * v[j]
                idx += 1

    return eigenvalues, eigenvectors

t_packed_copy = t_packed.copy()
pi_eigenvalues, pi_eigenvectors = power_eigen(t_packed_copy, n, 3, eigen_seed)

print(f"\n=== Top 3 Eigenvalues (power iteration) ===")
for i in range(3):
    print(f"  λ{i} = {pi_eigenvalues[i]:.10f}")

# Compute coordinates from power iteration
print(f"\n=== Initial Coordinates from Power Iteration (first 5) ===")
for i in range(min(5, n)):
    coords = []
    for d in range(3):
        if pi_eigenvalues[d] > 0:
            coords.append(np.sqrt(pi_eigenvalues[d]) * pi_eigenvectors[d][i])
        else:
            coords.append(0.0)
    print(f"  {i}: ({coords[0]:.8f}, {coords[1]:.8f}, {coords[2]:.8f})")
