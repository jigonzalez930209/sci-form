#!/usr/bin/env python3
"""Compare our pipeline with RDKit's step-by-step for the first failing molecules.

For each molecule that fails (RMSD > 0.5 Å), dump:
1. RDKit's initial coordinates (after eigendecomposition, before any FF)
2. RDKit's coordinates after firstMinimization
3. RDKit's coordinates after 3D FF minimization
4. Number of torsion/distance/angle/inversion contributions in each FF
"""

import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)

# Load reference data
with open("tests/fixtures/gdb20_reference.json") as f:
    ref_data = json.load(f)

# Load SMILES
smiles_list = []
with open("scripts/GDB20.50000.smi") as f:
    for line in f:
        parts = line.strip().split()
        if parts:
            smiles_list.append(parts[0])

# Find failing molecules by checking our Rust output
# For now, let's just test the first 100 molecules and find those we suspect fail
# We'll focus on molecules with triple bonds and aromatic rings

def compute_rmsd(coords1, coords2):
    """Compute RMSD between two coordinate matrices."""
    diff = coords1 - coords2
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))

def get_rdkit_intermediate_coords(smiles, seed=42):
    """Get RDKit conformer at various stages using the detailed API."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    
    # Standard ETKDGv3 embedding
    ps = rdDistGeom.ETKDGv3()
    ps.randomSeed = seed
    ps.useSmallRingTorsions = False
    ps.useMacrocycleTorsions = False
    
    cid = AllChem.EmbedMolecule(mol, ps)
    if cid < 0:
        return None
    
    conf = mol.GetConformer(cid)
    n = mol.GetNumAtoms()
    coords = np.array([[conf.GetAtomPosition(i).x,
                         conf.GetAtomPosition(i).y,
                         conf.GetAtomPosition(i).z] for i in range(n)])
    
    return {
        'n_atoms': n,
        'coords': coords,
        'heavy_atoms': sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1),
    }

def analyze_molecule(idx, smiles, ref_coords):
    """Analyze a specific molecule."""
    result = get_rdkit_intermediate_coords(smiles)
    if result is None:
        return None
    
    ref = np.array(ref_coords)
    our = result['coords']
    
    if ref.shape != our.shape:
        return None
    
    rmsd = compute_rmsd(ref, our)
    if rmsd > 1e-4:
        print(f"WARNING: RDKit reference mismatch for mol {idx}: RMSD={rmsd:.6f}")
    
    return {
        'idx': idx,
        'smiles': smiles,
        'n_atoms': result['n_atoms'],
        'rmsd_ref': rmsd,
    }

# Look for molecules with specific properties
count = 0
tested = 0
print("=== Analyzing molecules for failure patterns ===\n")

for idx, smi in enumerate(smiles_list[:200]):
    if str(idx) not in ref_data:
        continue
    
    ref = ref_data[str(idx)]
    ref_coords = np.array(ref['coordinates'])
    
    # Get basic info
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        continue
    mol = Chem.AddHs(mol)
    
    has_triple = any(b.GetBondType() == Chem.rdchem.BondType.TRIPLE for b in mol.GetBonds())
    has_aromatic = any(b.GetIsAromatic() for b in mol.GetBonds())
    
    # Focus on molecules likely to fail (triple bonds, aromatic)
    if not (has_triple or has_aromatic):
        continue
    
    tested += 1
    
    # Get the ETKDG details for this molecule
    ps = rdDistGeom.ETKDGv3()
    ps.randomSeed = 42
    ps.useSmallRingTorsions = False
    ps.useMacrocycleTorsions = False
    
    # Get bounds matrix
    mol2 = Chem.MolFromSmiles(smi)
    mol2 = Chem.AddHs(mol2)
    bm = rdDistGeom.GetMoleculeBoundsMatrix(mol2)
    n = mol2.GetNumAtoms()
    
    # Embed
    cid = AllChem.EmbedMolecule(mol2, ps)
    if cid < 0:
        continue
    
    conf = mol2.GetConformer(cid)
    rdkit_coords = np.array([[conf.GetAtomPosition(i).x,
                               conf.GetAtomPosition(i).y,
                               conf.GetAtomPosition(i).z] for i in range(n)])
    
    rmsd_check = compute_rmsd(ref_coords, rdkit_coords)
    
    # Now we need to know our RMSD. Since we can't easily run Rust here,
    # let's just print the molecule info and the reference coords for manual comparison
    
    if count < 20:
        # Print details for molecules that are likely to fail
        print(f"Mol {idx}: {smi}")
        print(f"  N={n}, heavy={sum(1 for a in mol2.GetAtoms() if a.GetAtomicNum() > 1)}")
        print(f"  Triple={has_triple}, Aromatic={has_aromatic}")
        print(f"  RDKit self-consistency: RMSD={rmsd_check:.6f}")
        
        # Count torsion angles
        from rdkit.Chem import rdMolTransforms
        n_rotatable = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol2)
        print(f"  Rotatable bonds: {n_rotatable}")
        
        # Get chirality info
        chiral_centers = Chem.FindMolChiralCenters(mol2, includeUnassigned=True)
        print(f"  Chiral centers: {len(chiral_centers)}")
        
        # Analyze bounds matrix
        tight_pairs = 0
        loose_pairs = 0
        for i in range(n):
            for j in range(i+1, n):
                gap = bm[i][j] - bm[j][i]  # upper - lower
                if gap < 0.5:
                    tight_pairs += 1
                elif gap > 3.0:
                    loose_pairs += 1
        print(f"  Bounds: {tight_pairs} tight (<0.5Å), {loose_pairs} loose (>3.0Å)")
        print()
        count += 1

print(f"\nTested {tested} molecules with triple/aromatic bonds")
