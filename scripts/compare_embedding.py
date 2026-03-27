#!/usr/bin/env python3
"""
Compare our embedding pipeline with RDKit's for specific failing molecules.
Implements the MinstdRand RNG and the full embedding pipeline in Python
to validate against both our Rust code and RDKit.
"""
import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
import sys
from fixture_io import load_json_fixture

class MinstdRand:
    """boost::minstd_rand (a=48271, c=0, m=2^31-1)"""
    A = 48271
    M = 2147483647
    
    def __init__(self, seed):
        self.state = seed if seed != 0 else 1
    
    def next(self):
        self.state = (self.A * self.state) % self.M
        return self.state
    
    def next_double(self):
        val = self.next()
        return (val - 1) / (self.M - 1)

def get_bounds_matrix_from_rdkit(mol):
    """Get the bounds matrix from RDKit after triangle smoothing."""
    bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
    return bm

def pick_distances(rng, bounds):
    """Pick random distances matching RDKit's pickRandomDistMat."""
    n = bounds.shape[0]
    dists_packed = {}
    # RDKit loop: i=1..n-1, j=0..i-1
    for i in range(1, n):
        for j in range(0, i):
            ub = bounds[min(i,j), max(i,j)]  # upper triangle
            lb = bounds[max(i,j), min(i,j)]  # lower triangle
            rval = rng.next_double()
            d = lb + rval * (ub - lb)
            dists_packed[(i, j)] = d
    return dists_packed

def build_metric_matrix(dists_packed, n):
    """Build the metric matrix T from packed distances."""
    # Build packed symmetric array
    d_size = n * (n + 1) // 2
    sq_packed = [0.0] * d_size
    sum_sq = 0.0
    
    for i in range(n):
        idx = i * (i + 1) // 2
        for j in range(i + 1):
            if i == j:
                sq_packed[idx + j] = 0.0
            else:
                d = dists_packed.get((i, j), dists_packed.get((j, i), 0.0))
                sq_packed[idx + j] = d * d
            sum_sq += sq_packed[idx + j]
    
    sum_sq /= (n * n)
    
    # Compute D0 (row averages minus grand mean)
    d0 = [0.0] * n
    for i in range(n):
        row_sum = 0.0
        for j in range(n):
            if i >= j:
                row_sum += sq_packed[i * (i + 1) // 2 + j]
            else:
                row_sum += sq_packed[j * (j + 1) // 2 + i]
        d0[i] = row_sum / n - sum_sq
    
    # Check D0
    EIGVAL_TOL = 0.001
    for i in range(n):
        if d0[i] < EIGVAL_TOL and n > 3:
            return None, sum_sq, d0
    
    # Build T matrix (packed symmetric)
    t_packed = [0.0] * d_size
    for i in range(n):
        idx = i * (i + 1) // 2
        for j in range(i + 1):
            sq_val = sq_packed[i * (i + 1) // 2 + j]
            t_packed[idx + j] = 0.5 * (d0[i] + d0[j] - sq_val)
    
    return t_packed, sum_sq, d0

def power_eigen_solver(n_eig, t_packed, n, seed):
    """Power eigenvalue solver matching RDKit exactly."""
    MAX_ITERATIONS = 1000
    TOLERANCE = 0.001
    HUGE_EIGVAL = 1.0e10
    TINY_EIGVAL = 1.0e-10
    
    mat = list(t_packed)  # copy since we modify
    eigenvalues = []
    eigenvectors = []
    
    for ei in range(n_eig):
        eig_val = -HUGE_EIGVAL
        seed += ei
        
        # Initialize random vector using MinstdRand
        rng = MinstdRand(seed if seed > 0 else 1)
        v = [rng.next_double() for _ in range(n)]
        
        # Normalize
        ns = sum(x * x for x in v)
        norm = ns ** 0.5
        if norm > 0:
            v = [x / norm for x in v]
        
        converged = False
        z = [0.0] * n
        
        for iteration in range(MAX_ITERATIONS):
            # z = mat * v (packed symmetric)
            for i in range(n):
                accum = 0.0
                idx = i * (i + 1) // 2
                for j in range(i + 1):
                    accum += mat[idx + j] * v[j]
                for j in range(i + 1, n):
                    accum += mat[j * (j + 1) // 2 + i] * v[j]
                z[i] = accum
            
            prev_val = eig_val
            # Find largest absolute value index
            max_abs = -1.0
            eval_id = 0
            for i in range(n):
                a = abs(z[i])
                if a > max_abs:
                    max_abs = a
                    eval_id = i
            
            eig_val = z[eval_id]
            
            if abs(eig_val) < TINY_EIGVAL:
                break
            
            # v = z / eigVal
            v = [z[i] / eig_val for i in range(n)]
            
            if abs(eig_val - prev_val) < TOLERANCE:
                converged = True
                break
        
        if not converged:
            return None
        
        # Normalize v
        ns = sum(x * x for x in v)
        norm = ns ** 0.5
        if norm > 0:
            v = [x / norm for x in v]
        
        eigenvalues.append(eig_val)
        eigenvectors.append(list(v))
        
        # Deflate
        for i in range(n):
            idx = i * (i + 1) // 2
            for j in range(i + 1):
                mat[idx + j] -= eig_val * v[i] * v[j]
    
    return eigenvalues, eigenvectors

def compute_initial_coords(rng, dists_packed, n, ndim=3):
    """Full initial coordinate computation matching RDKit."""
    EIGVAL_TOL = 0.001
    
    t_packed, sum_sq, d0 = build_metric_matrix(dists_packed, n)
    if t_packed is None:
        return None, sum_sq, d0, None, None
    
    n_eigs = min(ndim, n)
    eigen_seed = int(sum_sq * n)  # matching (int)(sumSqD2 * N)
    
    result = power_eigen_solver(n_eigs, t_packed, n, eigen_seed)
    if result is None:
        return None, sum_sq, d0, eigen_seed, None
    
    eigenvalues, eigenvectors = result
    
    # Process eigenvalues
    eig_sqrt = []
    zero_eigs = 0
    for i in range(n_eigs):
        if eigenvalues[i] > EIGVAL_TOL:
            eig_sqrt.append(eigenvalues[i] ** 0.5)
        elif abs(eigenvalues[i]) < EIGVAL_TOL:
            eig_sqrt.append(0.0)
            zero_eigs += 1
        else:
            eig_sqrt.append(eigenvalues[i])  # keep negative
    
    if zero_eigs >= 1 and n > 3:
        return None, sum_sq, d0, eigen_seed, eigenvalues
    
    # Generate coordinates
    coords = np.zeros((n, ndim))
    for i in range(n):
        for j in range(ndim):
            if j < n_eigs and eig_sqrt[j] >= 0.0:
                coords[i, j] = eig_sqrt[j] * eigenvectors[j][i]
            else:
                coords[i, j] = 1.0 - 2.0 * rng.next_double()
    
    return coords, sum_sq, d0, eigen_seed, eigenvalues

def pairwise_distance_rmsd(coords1, coords2):
    """Compute alignment-free RMSD based on pairwise distances."""
    n = coords1.shape[0]
    sum_sq = 0.0
    count = 0
    for i in range(n):
        for j in range(i + 1, n):
            d1 = np.linalg.norm(coords1[i] - coords1[j])
            d2 = np.linalg.norm(coords2[i] - coords2[j])
            sum_sq += (d1 - d2) ** 2
            count += 1
    return (sum_sq / count) ** 0.5

def main():
    mol_idx = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    
    # Load reference data
    data = load_json_fixture('tests/fixtures/gdb20_reference.json')
    
    mol_data = data[mol_idx]
    smiles = mol_data['smiles']
    n = len(mol_data['atoms'])
    print(f"Molecule {mol_idx}: {smiles}")
    print(f"  N atoms: {n}")
    
    # Get RDKit reference coords
    ref_coords = np.array([[a['x'], a['y'], a['z']] for a in mol_data['atoms']])
    
    # Build mol WITH hydrogens and get bounds matrix
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
    n_bm = bm.shape[0]
    print(f"  Bounds matrix shape: {bm.shape}")
    print(f"  N atoms in reference: {n}")
    
    # Simulate our embedding pipeline with ALL atoms (including H)
    rng = MinstdRand(42)
    
    # Pick random distances 
    dists = pick_distances(rng, bm)
    
    # Print first few distances for verification
    print(f"\n  First 10 picked distances (i,j -> d):")
    count = 0
    for i in range(1, n_bm):
        for j in range(0, i):
            if count < 10:
                print(f"    ({i},{j}) -> {dists[(i,j)]:.15f}")
            count += 1
    print(f"  RNG state after distance picking: {rng.state}")
    
    # Compute initial coordinates
    coords, sum_sq, d0, eigen_seed, eigenvalues = compute_initial_coords(
        rng, dists, n_bm, ndim=3
    )
    
    print(f"\n  sum_sq_all: {sum_sq:.15f}")
    if eigen_seed is not None:
        print(f"  eigen_seed: {eigen_seed}")
    
    if d0 is not None:
        print(f"  D0 values (first 5): {[f'{x:.15f}' for x in d0[:5]]}")
        neg_d0 = [i for i, x in enumerate(d0) if x < 0.001]
        if neg_d0:
            print(f"  D0 < EIGVAL_TOL at indices: {neg_d0}")
    
    if eigenvalues is not None:
        print(f"\n  Eigenvalues:")
        for i, ev in enumerate(eigenvalues):
            print(f"    [{i}] = {ev:.15f}")
    
    if coords is not None:
        print(f"\n  Initial coords (first 5 atoms):")
        for i in range(min(5, n_bm)):
            print(f"    atom[{i}] = ({coords[i,0]:.15f}, {coords[i,1]:.15f}, {coords[i,2]:.15f})")
        
        # Compare with reference (using matching atom count)
        n_compare = min(n, n_bm)
        rmsd = pairwise_distance_rmsd(coords[:n_compare], ref_coords[:n_compare, :3])
        print(f"\n  Pairwise RMSD (python initial vs RDKit final): {rmsd:.6f}")
    else:
        print("\n  Embedding FAILED (attempt 0)")
        print(f"  RNG state at failure: {rng.state}")
    
    # DUMP DATA for Rust comparison
    print(f"\n=== COMPARISON DATA ===")
    print(f"RNG_STATE_POST_DIST: {rng.state}")
    if eigen_seed is not None:
        print(f"EIGEN_SEED: {eigen_seed}")
    if eigenvalues is not None:
        for i, ev in enumerate(eigenvalues):
            print(f"EIGENVALUE_{i}: {ev:.20f}")
    if coords is not None:
        for i in range(n_bm):
            for j in range(3):
                print(f"COORD_{i}_{j}: {coords[i,j]:.20f}")

if __name__ == '__main__':
    main()
