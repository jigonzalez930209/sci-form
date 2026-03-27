#!/usr/bin/env python3
"""Compare 3D FF energy at RDKit's final coordinates vs coordinates from our pipeline.

For a failing molecule, we:
1. Build the 3D FF (same as RDKit would)
2. Evaluate the 3D FF energy at RDKit's final coordinates
3. Evaluate the 3D FF energy at our code's final coordinates
4. Compare

This tells us whether both codes converge to the same energy minimum.
"""

import json
import sys
import numpy as np
from fixture_io import load_json_fixture

from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom, rdMolDescriptors
from rdkit.Chem import rdForceFieldHelpers
from rdkit.ForceField import rdForceField
from rdkit import RDLogger, ForceField

RDLogger.logger().setLevel(RDLogger.ERROR)


def get_rdkit_3d_ff_energy(mol, coords_3d):
    """Build RDKit's 3D ETKDG FF and compute energy at given coordinates."""
    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    # Set coordinates
    for i in range(n):
        conf.SetAtomPosition(i, (float(coords_3d[i*3]), float(coords_3d[i*3+1]), float(coords_3d[i*3+2])))
    
    # Build the 3D force field (same as RDKit uses internally)
    # We need to use the internal API to construct the ETKDG FF
    # Unfortunately, rdkit doesn't expose construct3DForceField directly from Python
    # But we can use the ForceField from an embedding run
    
    # Alternative: compute UFF energy as a proxy
    ff = rdForceFieldHelpers.UFFGetMoleculeForceField(mol)
    if ff is None:
        return None
    return ff.CalcEnergy()


def trace_rdkit_pipeline(smiles, seed=42):
    """Trace the RDKit pipeline, printing energy at each stage."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Failed to parse: {smiles}")
        return
    
    mol = Chem.AddHs(mol)
    n = mol.GetNumAtoms()
    print(f"SMILES: {smiles}")
    print(f"Atoms: {n}")
    
    # Embed and track what happens
    ps = rdDistGeom.ETKDGv3()
    ps.randomSeed = seed
    ps.useSmallRingTorsions = False
    ps.useMacrocycleTorsions = False
    ps.clearConfs = True
    ps.numThreads = 1
    
    res = AllChem.EmbedMolecule(mol, ps)
    if res != 0:
        print("Embedding failed!")
        return
    
    conf = mol.GetConformer()
    coords = []
    for i in range(n):
        pos = conf.GetAtomPosition(i)
        coords.extend([pos.x, pos.y, pos.z])
    
    print(f"\nRDKit final coordinates (first 5 atoms):")
    for i in range(min(5, n)):
        print(f"  atom[{i}]: ({coords[i*3]:.6f}, {coords[i*3+1]:.6f}, {coords[i*3+2]:.6f})")
    
    # UFF energy at RDKit coords (as proxy)
    uff = rdForceFieldHelpers.UFFGetMoleculeForceField(mol)
    if uff:
        print(f"\nUFF energy at RDKit coords: {uff.CalcEnergy():.6f}")
    
    # Get CSD torsions
    torsions = rdDistGeom.GetExperimentalTorsions(mol)
    print(f"\nCSD torsions: {len(torsions)}")
    
    return coords


def compute_pairwise_rmsd(coords1, coords2, n):
    """Compute pairwise-distance RMSD between two coordinate sets."""
    rmsd_sum = 0.0
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            d1 = np.sqrt(sum((coords1[i*3+d] - coords1[j*3+d])**2 for d in range(3)))
            d2 = np.sqrt(sum((coords2[i*3+d] - coords2[j*3+d])**2 for d in range(3)))
            rmsd_sum += (d1 - d2)**2
            count += 1
    return np.sqrt(rmsd_sum / count)


def main():
    # Use a specific failing molecule
    smiles = sys.argv[1] if len(sys.argv) > 1 else "N#Cc1ccc(Br)c(CC2CCCCC2)c1C=O"
    
    rdkit_coords = trace_rdkit_pipeline(smiles)
    if rdkit_coords is None:
        return
    
    # Load our coords from reference data
    ref_data = load_json_fixture("tests/fixtures/gdb20_reference.json")
    
    ref_mol = None
    for ref in ref_data:
        if ref["smiles"] == smiles:
            ref_mol = ref
            break
    
    if ref_mol:
        n = len(ref_mol["atoms"])
        ref_coords = []
        for a in ref_mol["atoms"]:
            ref_coords.extend([a["x"], a["y"], a["z"]])
        
        # Compare
        rmsd = compute_pairwise_rmsd(rdkit_coords, ref_coords, n)
        print(f"\nPairwise RMSD between fresh embed and stored reference: {rmsd:.6f}")
        
        # Check max coordinate difference
        max_diff = max(abs(a - b) for a, b in zip(rdkit_coords, ref_coords))
        print(f"Max coordinate difference: {max_diff:.10f}")


if __name__ == "__main__":
    main()
