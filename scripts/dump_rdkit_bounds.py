"""
Dump RDKit intermediate embedding data for a specific molecule.
This dumps:
1. The smoothed bounds matrix
2. Topological distances
"""
import sys
import numpy as np

# Pick a failing molecule from the test results
SMILES = "C#CC(C)(COCCC1CCCNCC1)OCC(C)=O"  # worst failing molecule

def main():
    # Must run from outside the rdkit-master directory
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDistGeom
    
    mol = Chem.MolFromSmiles(SMILES)
    mol = Chem.AddHs(mol)
    
    n = mol.GetNumAtoms()
    print(f"Molecule: {SMILES}")
    print(f"Atoms: {n}")
    
    # Get bounds matrix (smoothed)
    bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
    
    print(f"\nBounds matrix shape: {bm.shape}")
    print(f"Bounds matrix dtype: {bm.dtype}")
    
    # Convention: bm[i,j] for i < j is the upper bound
    #             bm[i,j] for i > j is the lower bound
    
    # Save bounds matrix
    np.save("/tmp/rdkit_bounds.npy", bm)
    print("Saved bounds matrix to /tmp/rdkit_bounds.npy")
    
    # Print some key bounds
    print("\n=== 1-2 Bond bounds (first 5) ===")
    bond_count = 0
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i > j:
            i, j = j, i
        ub = bm[i, j]
        lb = bm[j, i]
        print(f"  ({i:2d},{j:2d}): lb={lb:.6f}, ub={ub:.6f}, range={ub-lb:.6f}")
        bond_count += 1
        if bond_count >= 10:
            break
    
    print("\n=== Bounds matrix checksums ===")
    # Sum of upper triangle
    ub_sum = 0.0
    lb_sum = 0.0
    for i in range(n):
        for j in range(i+1, n):
            ub_sum += bm[i, j]
            lb_sum += bm[j, i]
    print(f"  UB sum: {ub_sum:.10f}")
    print(f"  LB sum: {lb_sum:.10f}")
    
    # Print first row
    print(f"\n=== First row of bounds matrix (upper bounds, i=0) ===")
    for j in range(min(n, 20)):
        if j > 0:
            print(f"  bm[0,{j:2d}] = {bm[0,j]:.10f} (UB), bm[{j:2d},0] = {bm[j,0]:.10f} (LB)")
    
    # Atom info
    print(f"\n=== Atom info ===")
    for i in range(n):
        atom = mol.GetAtomWithIdx(i)
        hyb = str(atom.GetHybridization())
        print(f"  {i:2d}: {atom.GetSymbol():2s} hyb={hyb:5s} #H={atom.GetTotalNumHs()}")
        if i > 25:
            print(f"  ... ({n - 26} more atoms)")
            break

if __name__ == "__main__":
    main()
