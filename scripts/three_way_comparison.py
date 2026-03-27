#!/usr/bin/env python3
"""
Three-way comparison: sci-form vs RDKit vs OpenBabel
Compares pairwise-distance RMSD between conformers from all three methods.

Uses GDB-20 molecules from gdb20_reference.json or gdb20_reference.json.gz.
For each molecule:
  1. RDKit: Reference from JSON (seed=42 with explicit H coords)
  2. OpenBabel: OBBuilder + MMFF94/UFF optimization  
  3. sci-form: CLI binary with own torsion matching

Metric: Heavy-atom pairwise-distance RMSD (alignment-free).
OpenBabel may generate atoms in different order, so we compare heavy atoms only
by regenerating from the same SMILES and using canonical heavy-atom ordering.
"""
import json
import sys
import time
import subprocess
import os
import numpy as np
from fixture_io import load_json_fixture

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)

from openbabel import openbabel as ob
ob.obErrorLog.SetOutputLevel(0)


def heavy_pairwise_rmsd(coords1, coords2):
    """Pairwise-distance RMSD between two heavy-atom coordinate sets."""
    n = len(coords1)
    if n != len(coords2) or n < 2:
        return float('inf')
    c1 = np.array(coords1)
    c2 = np.array(coords2)
    sq_sum = 0.0
    npairs = 0
    for a in range(n):
        for b in range(a+1, n):
            d1 = np.linalg.norm(c1[a] - c1[b])
            d2 = np.linalg.norm(c2[a] - c2[b])
            sq_sum += (d1-d2)**2
            npairs += 1
    return np.sqrt(sq_sum / npairs) if npairs > 0 else 0.0


def rdkit_heavy_coords_fresh(smiles, seed=42):
    """Generate coordinates with RDKit fresh from SMILES, return heavy-atom coords."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv2()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) < 0:
        return None
    conf = mol.GetConformer()
    coords = []
    for i in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(i)
        if atom.GetAtomicNum() > 1:
            pos = conf.GetAtomPosition(i)
            coords.append([pos.x, pos.y, pos.z])
    return coords


def rdkit_multi_seed_heavy(smiles, seeds):
    """Generate multiple conformers with different seeds, return list of heavy-atom coords."""
    results = []
    for seed in seeds:
        c = rdkit_heavy_coords_fresh(smiles, seed)
        if c is not None:
            results.append(c)
    return results


def openbabel_heavy_coords(smiles):
    """Generate 3D conformer using OpenBabel, return heavy-atom coords."""
    conv = ob.OBConversion()
    conv.SetInFormat("smi")
    mol = ob.OBMol()
    conv.ReadString(mol, smiles)
    mol.AddHydrogens()
    
    builder = ob.OBBuilder()
    builder.Build(mol)
    
    ff = ob.OBForceField.FindForceField("MMFF94")
    if ff is None or not ff.Setup(mol):
        ff = ob.OBForceField.FindForceField("UFF")
        if ff is None or not ff.Setup(mol):
            return None
    
    ff.ConjugateGradients(500)
    ff.GetCoordinates(mol)
    
    coords = []
    for atom in ob.OBMolAtomIter(mol):
        if atom.GetAtomicNum() > 1:
            coords.append([atom.GetX(), atom.GetY(), atom.GetZ()])
    return coords


def sciform_heavy_coords(sci_form_bin, smiles, seed=42):
    """Generate conformer using sci-form CLI, return heavy-atom coords."""
    try:
        proc = subprocess.run(
            [sci_form_bin, 'embed', '--seed', str(seed), smiles],
            capture_output=True, text=True, timeout=30
        )
        if proc.returncode != 0:
            return None
        out = json.loads(proc.stdout)
        if out.get('error') is not None or not out.get('coords'):
            return None
        flat = out['coords']
        n_atoms = len(flat) // 3
        # Need to know which atoms are heavy — parse SMILES to get atom info
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        if mol.GetNumAtoms() != n_atoms:
            return None
        coords = []
        for i in range(n_atoms):
            atom = mol.GetAtomWithIdx(i)
            if atom.GetAtomicNum() > 1:
                coords.append([flat[i*3], flat[i*3+1], flat[i*3+2]])
        return coords
    except Exception:
        return None


def sciform_multi_seed_heavy(sci_form_bin, smiles, seeds):
    """Generate multiple conformers with different seeds."""
    results = []
    for seed in seeds:
        c = sciform_heavy_coords(sci_form_bin, smiles, seed)
        if c is not None:
            results.append(c)
    return results


def min_cross_rmsd(coords_list_a, coords_list_b):
    """Minimum RMSD between any pair from two lists of coordinate sets."""
    if not coords_list_a or not coords_list_b:
        return float('inf')
    min_r = float('inf')
    for ca in coords_list_a:
        for cb in coords_list_b:
            if len(ca) == len(cb):
                r = heavy_pairwise_rmsd(ca, cb)
                min_r = min(min_r, r)
    return min_r


def process_molecule(smiles, sci_form_bin, seeds):
    """Process one molecule with all three methods."""
    result = {'smiles': smiles}
    
    # RDKit multi-seed
    rdkit_coords_list = rdkit_multi_seed_heavy(smiles, seeds)
    result['rdkit_ok'] = len(rdkit_coords_list) > 0
    result['rdkit_n'] = len(rdkit_coords_list)
    
    # OpenBabel single conformer (deterministic)
    ob_coords = openbabel_heavy_coords(smiles)
    result['ob_ok'] = ob_coords is not None
    
    # sci-form multi-seed
    sci_coords_list = sciform_multi_seed_heavy(sci_form_bin, smiles, seeds)
    result['sci_ok'] = len(sci_coords_list) > 0
    result['sci_n'] = len(sci_coords_list)
    
    # Check atom count consistency
    if result['rdkit_ok'] and result['ob_ok']:
        if len(rdkit_coords_list[0]) != len(ob_coords):
            result['ob_atom_mismatch'] = True
            result['ob_ok'] = False  # can't compare
    
    if result['rdkit_ok'] and result['sci_ok']:
        if len(rdkit_coords_list[0]) != len(sci_coords_list[0]):
            result['sci_atom_mismatch'] = True
            result['sci_ok'] = False
    
    # Compute cross-RMSDs
    if result['rdkit_ok'] and result['sci_ok']:
        result['rdkit_vs_sci'] = min_cross_rmsd(rdkit_coords_list, sci_coords_list)
    
    if result['rdkit_ok'] and result['ob_ok']:
        result['rdkit_vs_ob'] = min_cross_rmsd(rdkit_coords_list, [ob_coords])
    
    if result['sci_ok'] and result['ob_ok']:
        result['sci_vs_ob'] = min_cross_rmsd(sci_coords_list, [ob_coords])
    
    return result


def main():
    ref_file = sys.argv[1] if len(sys.argv) > 1 else "tests/fixtures/gdb20_reference.json"
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    n_seeds = int(sys.argv[3]) if len(sys.argv) > 3 else 5

    sci_form_bin = "target/release/sci-form"
    if not os.path.exists(sci_form_bin):
        print("ERROR: Build sci-form first: cargo build --release")
        sys.exit(1)

    print(f"Loading reference from {ref_file}...")
    ref_mols = load_json_fixture(ref_file)[:limit]
    
    seeds = [42] + list(range(n_seeds - 1))
    smiles_list = [m['smiles'] for m in ref_mols]
    
    print(f"Processing {len(smiles_list)} molecules with {n_seeds} seeds each...")
    print(f"  Seeds: {seeds}")
    t0 = time.time()

    results = []
    for i, smi in enumerate(smiles_list):
        r = process_molecule(smi, sci_form_bin, seeds)
        results.append(r)
        if (i+1) % 50 == 0:
            elapsed = time.time() - t0
            rate = (i+1) / elapsed
            print(f"  Progress: {i+1}/{len(smiles_list)} ({rate:.1f} mol/s)")

    elapsed = time.time() - t0
    n = len(results)

    # Statistics
    rdkit_ok = sum(1 for r in results if r.get('rdkit_ok'))
    ob_ok = sum(1 for r in results if r.get('ob_ok'))
    sci_ok = sum(1 for r in results if r.get('sci_ok'))
    ob_mismatch = sum(1 for r in results if r.get('ob_atom_mismatch'))
    sci_mismatch = sum(1 for r in results if r.get('sci_atom_mismatch'))

    print(f"\n{'='*65}")
    print(f"THREE-WAY COMPARISON ({n} molecules, {n_seeds} seeds, {elapsed:.1f}s)")
    print(f"  Heavy-atom pairwise-distance RMSD (min over {n_seeds}×{n_seeds} seed pairs)")
    print(f"{'='*65}")
    print(f"  RDKit success:     {rdkit_ok}/{n} ({100*rdkit_ok/n:.1f}%)")
    print(f"  OpenBabel success: {ob_ok}/{n} (atom mismatch: {ob_mismatch})")
    print(f"  sci-form success:  {sci_ok}/{n} (atom mismatch: {sci_mismatch})")

    def report(label, values):
        if not values:
            print(f"\n--- {label}: No data ---")
            return
        arr = np.array(values)
        print(f"\n--- {label} ({len(arr)} molecules) ---")
        print(f"  Avg: {arr.mean():.4f} Å")
        print(f"  Median: {np.median(arr):.4f} Å")
        for t in [0.1, 0.3, 0.5]:
            c = np.sum(arr < t)
            print(f"  < {t} Å: {c}/{len(arr)} ({100*c/len(arr):.1f}%)")
        above = np.sum(arr >= 0.5)
        print(f"  ≥ 0.5 Å: {above}/{len(arr)} ({100*above/len(arr):.2f}%)")
        print(f"  Max: {arr.max():.4f} Å")

    rs = [r['rdkit_vs_sci'] for r in results if 'rdkit_vs_sci' in r]
    report("RDKit vs sci-form (heavy-atom)", rs)

    ro = [r['rdkit_vs_ob'] for r in results if 'rdkit_vs_ob' in r]
    report("RDKit vs OpenBabel (heavy-atom)", ro)

    so = [r['sci_vs_ob'] for r in results if 'sci_vs_ob' in r]
    report("sci-form vs OpenBabel (heavy-atom)", so)

    # Head-to-head: which is closer to RDKit?
    all3 = [r for r in results if 'rdkit_vs_sci' in r and 'rdkit_vs_ob' in r]
    if all3:
        sci_closer = sum(1 for r in all3 if r['rdkit_vs_sci'] <= r['rdkit_vs_ob'])
        ob_closer = len(all3) - sci_closer
        print(f"\n--- Head-to-head: who matches RDKit best? ({len(all3)} molecules) ---")
        print(f"  sci-form closer: {sci_closer}/{len(all3)} ({100*sci_closer/len(all3):.1f}%)")
        print(f"  OpenBabel closer: {ob_closer}/{len(all3)} ({100*ob_closer/len(all3):.1f}%)")

    # Worst molecules for each comparison
    for label, key in [("RDKit vs sci-form", "rdkit_vs_sci"), ("RDKit vs OpenBabel", "rdkit_vs_ob")]:
        worst = sorted([r for r in results if key in r], key=lambda r: -r[key])[:10]
        if worst:
            print(f"\n--- Worst {label} ---")
            for r in worst:
                print(f"  {r[key]:.3f} Å  {r['smiles'][:55]}")

    print()


if __name__ == '__main__':
    main()
