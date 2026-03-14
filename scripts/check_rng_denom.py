#!/usr/bin/env python3
"""Check if the RNG double conversion matches by comparing
the random distances for the first attempt of a known molecule.

We compare by checking the metric matrix and eigen seed.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom, rdMolDescriptors
from rdkit import RDLogger
import numpy as np
import json

RDLogger.logger().setLevel(RDLogger.ERROR)

class MinstdRand:
    """Our Python implementation - must match Rust."""
    A = 48271
    M = 2147483647
    
    def __init__(self, seed):
        self.state = seed
    
    def next_int(self):
        self.state = (self.state * self.A) % self.M
        return self.state
    
    def next_double(self):
        val = self.next_int()
        return (val - 1.0) / (self.M - 1.0)  # denominator = 2147483646

SMILES = "C#CC(C)(COCCC1CCCNCC1)OCC(C)=O"

def main():
    mol = Chem.MolFromSmiles(SMILES)
    mol = Chem.AddHs(mol)
    n = mol.GetNumAtoms()
    
    # Load reference data to get bounds
    with open("tests/fixtures/gdb20_reference.json") as f:
        refs = json.load(f)
    
    # Find this molecule
    ref = None
    for r in refs:
        if r["smiles"] == SMILES:
            ref = r
            break
    
    if ref is None:
        print("Molecule not found in reference data")
        return
    
    print(f"Found molecule: {SMILES}")
    print(f"Atoms: {n}, Torsions: {len(ref['torsions'])}")
    
    # The first attempt uses seed=42
    # pick_rdkit_distances consumes n*(n-1)/2 RNG values
    rng = MinstdRand(42)
    
    n_dist = n * (n - 1) // 2
    print(f"Random distances to generate: {n_dist}")
    
    # Generate first few distances and compute sum_sq_all
    dists = [[0.0] * n for _ in range(n)]
    doubles = []
    for i in range(1, n):
        for j in range(i):
            val = rng.next_double()
            doubles.append(val)
    
    print(f"First 5 random doubles: {doubles[:5]}")
    print(f"Last 5 random doubles: {doubles[-5:]}")
    print(f"RNG state after distances: {rng.state}")
    
    # The power iteration uses a local RNG seeded with (sum_sq_all * n) as i32
    # We need the bounds matrix to compute the actual distances
    # But we can at least check the RNG state

    # After pick_rdkit_distances, the main RNG state is known.
    # The power iteration doesn't consume main RNG.
    # Coordinate fill for negative eigenvalues consumes from main RNG.
    
    # For this molecule with 47 atoms and 3rd eigenvalue negative:
    # Dimensions 0, 1: positive (use eigenvector, no main RNG)
    # Dimension 2: negative (use main RNG, consume 47 values)
    
    print(f"\nAfter pick_rdkit_distances, next few main RNG doubles:")
    for i in range(5):
        print(f"  {i}: {rng.next_double():.18f}")
    
    # Now check alternate denominator: (M - 2) = 2147483645
    rng2_state = 42
    A = 48271
    M = 2147483647
    
    def next_int_alt(state):
        return (state * A) % M
    
    def next_double_alt(val):
        """Alternate formula: (val - 1) / (M - 2)"""
        return (val - 1.0) / (M - 2.0)  # denominator = 2147483645
    
    rng2_state = next_int_alt(42)
    d_ours = (rng2_state - 1.0) / 2147483646.0
    d_alt = (rng2_state - 1.0) / 2147483645.0
    print(f"\nFirst value comparison:")
    print(f"  int: {rng2_state}")
    print(f"  Ours (denom 2147483646): {d_ours:.18f}")
    print(f"  Alt  (denom 2147483645): {d_alt:.18f}")
    print(f"  Difference: {abs(d_ours - d_alt):.2e}")
    print(f"  Ratio difference: {abs(1 - d_ours/d_alt):.2e}")

if __name__ == "__main__":
    main()
