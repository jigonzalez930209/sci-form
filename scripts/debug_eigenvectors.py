#!/usr/bin/env python3
"""Debug: dump eigenvalues, eigenvectors, and initial coords for a specific molecule.

Patches RDKit's internals to extract intermediate data from computeInitialCoords.
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
from rdkit import RDLogger
import json

RDLogger.logger().setLevel(RDLogger.ERROR)

# Use the first failing molecule
SMILES = "C#CC(C)(COCCC1CCCNCC1)OCC(C)=O"
SEED = 42


def minstd_rand(seed):
    """Boost minstd_rand compatible generator."""
    state = seed
    while True:
        state = (state * 48271) % 2147483647
        yield state / 2147483646.0


def main():
    mol = Chem.MolFromSmiles(SMILES)
    mol = Chem.AddHs(mol)
    n = mol.GetNumAtoms()
    print(f"Molecule: {SMILES}")
    print(f"Atoms: {n}")

    # Get bounds matrix
    ps = rdDistGeom.srETKDGv3()
    ps.randomSeed = SEED
    bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
    print(f"Bounds matrix shape: {bm.shape}")

    # Build distance matrix using our RNG (matching RDKit's pick_rdkit_distances)
    rng = minstd_rand(SEED)
    dists = np.zeros((n, n))
    for i in range(1, n):
        for j in range(i):
            r = next(rng)
            lb = bm[j, i]  # lower bound: below diagonal after smoothing? 
            # Actually, rdDistGeom.GetMoleculeBoundsMatrix returns upper bound at [i,j] where i<j
            # and lower bound at [j,i] where j>i
            # But numpy: bm[i,j] = upper, bm[j,i] = lower when i < j
            ub = bm[i, j] if i > j else bm[j, i]
            lb_val = bm[j, i] if i > j else bm[i, j]
            # RDKit: d = lb + r * (ub - lb) where ub is upper, lb is lower
            # Check bounds matrix convention
            d = lb_val + r * (ub - lb_val)
            dists[i, j] = d
            dists[j, i] = d

    # Build squared distance matrix
    sq_dists = dists * dists

    # Compute sumSqD2 from packed lower triangle
    sum_sq = 0.0
    for i in range(n):
        for j in range(i + 1):
            sum_sq += sq_dists[i, j]
    sum_sq_d2 = sum_sq / (n * n)
    print(f"sumSqD2 = {sum_sq_d2:.15e}")
    print(f"eigen_seed = {int(sum_sq_d2 * n)}")

    # Compute D0
    d0 = np.zeros(n)
    for i in range(n):
        row_sum = 0.0
        for j in range(n):
            row_sum += sq_dists[i, j]  # symmetric access
        d0[i] = row_sum / n - sum_sq_d2
    print(f"D0[0:5] = {d0[:5]}")

    # Build metric matrix T (packed symmetric)
    T_full = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1):
            T_full[i, j] = 0.5 * (d0[i] + d0[j] - sq_dists[i, j])
            T_full[j, i] = T_full[i, j]

    print(f"T[0,0:5] = {T_full[0, :5]}")

    # Power eigendecomposition matching RDKit
    eigen_seed = int(sum_sq_d2 * n)
    nEigs = 3  # 3D embedding

    eigenvalues = []
    eigenvectors = []

    T_work = T_full.copy()
    # Pack into lower triangle for power iteration
    t_packed = []
    for i in range(n):
        for j in range(i + 1):
            t_packed.append(T_work[i, j])
    t_packed = np.array(t_packed, dtype=np.float64)

    def get_packed(i, j):
        if i >= j:
            return t_packed[i * (i + 1) // 2 + j]
        else:
            return t_packed[j * (j + 1) // 2 + i]

    def set_packed(i, j, val):
        if i >= j:
            t_packed[i * (i + 1) // 2 + j] = val
        else:
            t_packed[j * (j + 1) // 2 + i] = val

    for ei in range(nEigs):
        eig_val = -1e10
        eigen_seed += ei

        # Initialize random vector (RDKit's setToRandom + normalize)
        rng2 = minstd_rand(eigen_seed)
        v = np.array([next(rng2) for _ in range(n)])
        norm = np.sqrt(np.sum(v * v))
        v /= norm

        converged = False
        for it in range(1000):
            # z = T * v (using packed symmetric)
            z = np.zeros(n)
            for i in range(n):
                accum = 0.0
                row_start = i * (i + 1) // 2
                for j in range(i + 1):
                    accum += t_packed[row_start + j] * v[j]
                for j in range(i + 1, n):
                    accum += t_packed[j * (j + 1) // 2 + i] * v[j]
                z[i] = accum

            prev_val = eig_val
            eval_id = np.argmax(np.abs(z))
            eig_val = z[eval_id]

            if abs(eig_val) < 1e-10:
                break

            v = z / eig_val

            if abs(eig_val - prev_val) < 0.001:
                converged = True
                break

        if not converged:
            print(f"  Eigenvector {ei}: NOT CONVERGED")
            break

        # Normalize
        norm = np.sqrt(np.sum(v * v))
        v /= norm

        eigenvalues.append(eig_val)
        eigenvectors.append(v.copy())

        print(f"  Eigenvector {ei}: eigenvalue={eig_val:.15e}, v[0:5]={v[:5]}, iters={it+1}")

        # Deflate
        for i in range(n):
            for j in range(i + 1):
                t_packed[i * (i + 1) // 2 + j] -= eig_val * v[i] * v[j]

    # Compute initial coordinates
    print("\nInitial coordinates (first 5 atoms):")
    eig_sqrt = [np.sqrt(ev) if ev > 0.001 else ev for ev in eigenvalues]
    for i in range(min(5, n)):
        coords = [eig_sqrt[j] * eigenvectors[j][i] if eig_sqrt[j] >= 0 else 0.0 for j in range(3)]
        print(f"  Atom {i}: [{coords[0]:.10f}, {coords[1]:.10f}, {coords[2]:.10f}]")

    # Also embed with RDKit for comparison
    mol2 = Chem.MolFromSmiles(SMILES)
    mol2 = Chem.AddHs(mol2)
    res = AllChem.EmbedMolecule(mol2, randomSeed=SEED)
    if res == 0:
        conf = mol2.GetConformer()
        print("\nRDKit final coordinates (first 5 atoms):")
        for i in range(min(5, n)):
            pos = conf.GetAtomPosition(i)
            print(f"  Atom {i}: [{pos.x:.10f}, {pos.y:.10f}, {pos.z:.10f}]")


if __name__ == "__main__":
    main()
