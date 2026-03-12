#!/usr/bin/env python3
"""Verify MinstdRand double conversion matches RDKit/boost."""

# Our implementation:
# state = seed
# next_int: state = (state * 48271) % 2147483647; return state
# next_double: (next_int() - 1) / (2147483647 - 1)  [= (val-1)/2147483646]

class MinstdRand:
    A = 48271
    M = 2147483647
    
    def __init__(self, seed):
        self.state = seed
    
    def next_int(self):
        self.state = (self.state * self.A) % self.M
        return self.state
    
    def next_double(self):
        val = self.next_int()
        return (val - 1.0) / (self.M - 1.0)

# Generate first 10 doubles with seed=42
rng = MinstdRand(42)
print("Our MinstdRand doubles (seed=42):")
for i in range(10):
    v = rng.next_double()
    print(f"  {i}: int={rng.state}, double={v:.18f}")

# Now use RDKit to generate the same
try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdDistGeom, AllChem
    import ctypes
    print("\nRDKit version check available")
    
    # We can't directly access RDKit's internal RNG.
    # But we can check by generating a KNOWN embedding and comparing coordinates.
    # A simpler check: use a 2-atom molecule (H-H) where the embedding is trivial.
    mol = Chem.MolFromSmiles("[H][H]")
    mol = Chem.AddHs(mol)
    
    ps = rdDistGeom.EmbedParameters()
    ps.randomSeed = 42
    ps.useExpTorsionAnglePrefs = False
    ps.useBasicKnowledge = False
    
    cid = AllChem.EmbedMolecule(mol, ps)
    if cid >= 0:
        conf = mol.GetConformer(cid)
        print(f"\nH-H molecule coordinates:")
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            print(f"  Atom {i}: ({pos.x:.18f}, {pos.y:.18f}, {pos.z:.18f})")
    
    # Try a larger molecule for more RNG values
    mol2 = Chem.MolFromSmiles("C")
    mol2 = Chem.AddHs(mol2)
    cid2 = AllChem.EmbedMolecule(mol2, ps)
    if cid2 >= 0:
        conf2 = mol2.GetConformer(cid2)
        print(f"\nMethane coordinates:")
        for i in range(mol2.GetNumAtoms()):
            pos = conf2.GetAtomPosition(i)
            print(f"  Atom {i}: ({pos.x:.18f}, {pos.y:.18f}, {pos.z:.18f})")
    
except ImportError:
    print("RDKit not available for comparison")
