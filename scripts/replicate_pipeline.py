#!/usr/bin/env python3
"""Replicate the Rust embedding pipeline in Python to find where results diverge from RDKit.

For a given molecule, this script:
1. Computes bounds (same as our Rust code, verified identical to RDKit's)
2. Samples random distances (using minstd_rand with seed 42)
3. Builds metric matrix
4. Runs power iteration eigendecomposition
5. Prints eigenvalues and initial coordinates
6. Compares with RDKit's actual conformer coordinates
"""

import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
from fixture_io import load_json_fixture

class MinstdRand:
    """Match our Rust MinstdRand exactly."""
    A = 48271
    M = 2147483647  # 2^31 - 1

    def __init__(self, seed):
        self.state = int(seed)

    def next_int(self):
        self.state = (self.state * self.A) % self.M
        return self.state

    def next_double(self):
        val = self.next_int()
        return (val - 1.0) / (self.M - 1.0)


def pick_random_distances(rng, bounds):
    """Match our Rust pick_rdkit_distances exactly."""
    n = bounds.shape[0]
    dists = np.zeros((n, n))
    for i in range(1, n):
        for j in range(i):
            ub = bounds[min(i,j), max(i,j)]  # upper in (row<col)
            lb = bounds[max(i,j), min(i,j)]  # lower in (row>col)
            rval = rng.next_double()
            d = lb + rval * (ub - lb)
            dists[i, j] = d
            dists[j, i] = d
    return dists


def build_metric_matrix_packed(dists):
    """Match our Rust compute_initial_coords_rdkit metric matrix exactly."""
    n = dists.shape[0]
    d_size = n * (n + 1) // 2
    sq_packed = np.zeros(d_size)
    sum_sq_all = 0.0

    for i in range(n):
        id_i = i * (i + 1) // 2
        for j in range(i + 1):
            d = dists[i, j]
            sq_packed[id_i + j] = d * d
            sum_sq_all += d * d

    sum_sq_all /= (n * n)

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

    # Build packed metric matrix
    t_packed = np.zeros(d_size)
    for i in range(n):
        id_i = i * (i + 1) // 2
        for j in range(i + 1):
            sq_val = sq_packed[id_i + j] if i >= j else sq_packed[j * (j + 1) // 2 + i]
            t_packed[id_i + j] = 0.5 * (d0[i] + d0[j] - sq_val)

    return t_packed, d0, sum_sq_all


def power_eigen_solver(num_eig, mat, n, seed):
    """Match our Rust power_eigen_solver exactly."""
    MAX_ITERATIONS = 1000
    TOLERANCE = 0.001
    HUGE_EIGVAL = 1.0e10
    TINY_EIGVAL = 1.0e-10

    eigenvalues = []
    eigenvectors = []

    for ei in range(num_eig):
        eig_val = -HUGE_EIGVAL
        s = seed + ei

        # Initialize random vector using minstd_rand
        rng = MinstdRand(s)
        v = np.array([rng.next_double() for _ in range(n)])
        norm = np.sqrt(np.sum(v * v))
        if norm > 0:
            v /= norm

        converged = False
        z = np.zeros(n)

        for iteration in range(MAX_ITERATIONS):
            # z = mat * v (packed symmetric format)
            for i in range(n):
                z[i] = 0.0
                id_i = i * (i + 1) // 2
                for j in range(i + 1):
                    z[i] += mat[id_i + j] * v[j]
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
            id_i = i * (i + 1) // 2
            for j in range(i + 1):
                mat[id_i + j] -= eig_val * v[i] * v[j]

    return eigenvalues, eigenvectors


# Load a specific failing molecule
ref_mols = load_json_fixture("tests/fixtures/gdb20_reference.json")

# Find a molecule with high RMSD from our test output
# Pick one from the worst-RMSD list
TARGET_SMILES = None
for mol in ref_mols[:200]:
    smiles = mol['smiles']
    has_triple = any(b['order'] == 'TRIPLE' for b in mol['bonds'])
    has_aromatic = any(b['order'] == 'AROMATIC' for b in mol['bonds'])
    n_heavy = sum(1 for a in mol['atoms'] if a['element'] != 1)
    if has_aromatic and has_triple and n_heavy <= 16:
        TARGET_SMILES = smiles
        TARGET_MOL = mol
        break

if TARGET_SMILES is None:
    # Fallback: pick first one with triple bond
    for mol in ref_mols[:100]:
        if any(b['order'] == 'TRIPLE' for b in mol['bonds']):
            TARGET_SMILES = mol['smiles']
            TARGET_MOL = mol
            break

print(f"Molecule: {TARGET_SMILES}")
n = len(TARGET_MOL['atoms'])
print(f"Atoms: {n}")

# Get RDKit's bounds matrix for this molecule
rdkit_mol = Chem.MolFromSmiles(TARGET_SMILES)
rdkit_mol = Chem.AddHs(rdkit_mol)
rdkit_bm = rdDistGeom.GetMoleculeBoundsMatrix(rdkit_mol)

print(f"\n=== Step 1: Bounds Matrix ===")
# Convert to our format (lower triangle = lower bound, upper triangle = upper bound)
bounds = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        bounds[i, j] = rdkit_bm[i][j]  # upper
        bounds[j, i] = rdkit_bm[j][i]  # lower
print(f"Bounds computed: {n}x{n}")

# Step 2: Sample distances
print(f"\n=== Step 2: Distance Sampling (seed=42) ===")
rng = MinstdRand(42)
dists = pick_random_distances(rng, bounds)

# Print first few distances
print("First 5 sampled distances:")
count = 0
for i in range(1, n):
    for j in range(i):
        if count < 5:
            print(f"  ({i},{j}): {dists[i,j]:.10f}")
            count += 1

# Step 3: Build metric matrix
print(f"\n=== Step 3: Metric Matrix ===")
t_packed, d0, sum_sq_all = build_metric_matrix_packed(dists)
print(f"sum_sq_all: {sum_sq_all:.10f}")
print(f"D0[0:5]: {d0[:5]}")

EIGVAL_TOL = 1e-3
bad_d0 = any(d0[i] < EIGVAL_TOL for i in range(n)) and n > 3
print(f"Any D0 < tolerance: {bad_d0}")
if bad_d0:
    print("  D0 values:", d0)
    print("  EMBEDDING WOULD FAIL (D0 check)")

# Step 4: Power iteration
print(f"\n=== Step 4: Eigendecomposition ===")
ndim = 3  # or 4 for chiral molecules
eigen_seed = int(sum_sq_all * n)
print(f"Eigen seed: {eigen_seed}")

t_copy = t_packed.copy()
result = power_eigen_solver(ndim, t_copy, n, eigen_seed)
if result is None:
    print("Power iteration FAILED to converge!")
else:
    eigenvalues, eigenvectors = result
    print(f"Eigenvalues: {eigenvalues}")
    for d in range(ndim):
        ev = eigenvectors[d]
        print(f"  Eigvec[{d}] first 5: {ev[:5]}")
        print(f"  Eigvec[{d}] norm: {np.sqrt(np.sum(ev*ev)):.10f}")

    # Step 5: Build initial coordinates
    coords = np.zeros((n, ndim))
    zero_eigs = 0
    for d in range(ndim):
        ev = eigenvalues[d]
        if ev > EIGVAL_TOL:
            coords[:, d] = np.sqrt(ev) * eigenvectors[d]
        elif abs(ev) < EIGVAL_TOL:
            zero_eigs += 1
            for i in range(n):
                coords[i, d] = 1.0 - 2.0 * rng.next_double()
        else:
            for i in range(n):
                coords[i, d] = 1.0 - 2.0 * rng.next_double()

    print(f"\nZero eigenvalues: {zero_eigs}")
    if zero_eigs >= 1 and n > 3:
        print("EMBEDDING WOULD FAIL (zero eigenvalue)")
    else:
        print(f"\n=== Step 5: Initial Coordinates ===")
        print("First 5 atoms:")
        for i in range(min(5, n)):
            print(f"  {i}: ({coords[i,0]:.6f}, {coords[i,1]:.6f}, {coords[i,2]:.6f})")

    # Step 6: Compare with RDKit's conformer
    print(f"\n=== Step 6: RDKit Reference Coordinates ===")
    print("First 5 atoms:")
    for i in range(min(5, n)):
        a = TARGET_MOL['atoms'][i]
        print(f"  {i}: ({a['x']:.6f}, {a['y']:.6f}, {a['z']:.6f})")

    # Compute pairwise distance RMSD
    ref_dists = []
    our_dists = []
    for i in range(n):
        for j in range(i+1, n):
            ai, aj = TARGET_MOL['atoms'][i], TARGET_MOL['atoms'][j]
            ref_d = np.sqrt((ai['x']-aj['x'])**2 + (ai['y']-aj['y'])**2 + (ai['z']-aj['z'])**2)
            ref_dists.append(ref_d)
    print(f"\nPairwise RMSD between initial coords and RDKit reference: N/A (initial coords are pre-optimization)")
