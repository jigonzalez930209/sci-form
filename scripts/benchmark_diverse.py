#!/usr/bin/env python3
"""
Benchmark: sci-form vs RDKit on diverse molecule set.

Compares heavy-atom pairwise-distance RMSD between sci-form and RDKit 
conformers using the same diverse molecule set from the Rust tests.
Multi-seed ensemble comparison for fair evaluation.
"""
import json
import subprocess
import sys
import time
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)


def heavy_pairwise_rmsd(coords1, coords2):
    """Alignment-free pairwise-distance RMSD between two heavy-atom coord sets."""
    n = len(coords1)
    if n != len(coords2) or n < 2:
        return float('inf')
    c1, c2 = np.array(coords1), np.array(coords2)
    sq_sum, npairs = 0.0, 0
    for a in range(n):
        for b in range(a + 1, n):
            d1 = np.linalg.norm(c1[a] - c1[b])
            d2 = np.linalg.norm(c2[a] - c2[b])
            sq_sum += (d1 - d2) ** 2
            npairs += 1
    return np.sqrt(sq_sum / npairs) if npairs > 0 else 0.0


def rdkit_heavy_coords(smiles, seed=42):
    """Generate conformer with RDKit, return heavy-atom coords."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv2()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) < 0:
        return None
    conf = mol.GetConformer()
    return [
        [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
        for i in range(mol.GetNumAtoms())
        if mol.GetAtomWithIdx(i).GetAtomicNum() > 1
    ]


def sciform_heavy_coords(binary, smiles, seed=42):
    """Generate conformer with sci-form CLI, return heavy-atom coords."""
    try:
        proc = subprocess.run(
            [binary, 'embed', '--seed', str(seed), smiles],
            capture_output=True, text=True, timeout=30
        )
        if proc.returncode != 0:
            return None
        out = json.loads(proc.stdout)
        if out.get('error') or not out.get('coords'):
            return None
        flat = out['coords']
        elements = out['elements']
        n = len(elements)
        return [
            [flat[i * 3], flat[i * 3 + 1], flat[i * 3 + 2]]
            for i in range(n)
            if elements[i] > 1
        ]
    except Exception:
        return None


def load_diverse_molecules(path):
    """Load diverse_molecules.smi."""
    entries = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 3:
                entries.append({
                    'smiles': parts[0],
                    'name': parts[1],
                    'category': parts[2],
                })
    return entries


def main():
    mol_file = sys.argv[1] if len(sys.argv) > 1 else "tests/fixtures/diverse_molecules.smi"
    n_seeds = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    binary = "target/release/sci-form"

    entries = load_diverse_molecules(mol_file)
    seeds = [42] + list(range(n_seeds - 1))

    print(f"Benchmark: sci-form vs RDKit ({len(entries)} molecules, {n_seeds} seeds)")
    t0 = time.time()

    results = []
    categories = {}

    for i, entry in enumerate(entries):
        smi = entry['smiles']
        cat = entry['category']

        if cat not in categories:
            categories[cat] = {'total': 0, 'rdkit_ok': 0, 'sci_ok': 0, 'compared': 0, 'rmsds': []}
        categories[cat]['total'] += 1

        # Multi-seed RDKit
        rdkit_list = []
        for s in seeds:
            c = rdkit_heavy_coords(smi, s)
            if c is not None:
                rdkit_list.append(c)

        if not rdkit_list:
            results.append({'smiles': smi, 'name': entry['name'], 'cat': cat, 'rdkit_ok': False})
            continue
        categories[cat]['rdkit_ok'] += 1

        # Multi-seed sci-form
        sci_list = []
        for s in seeds:
            c = sciform_heavy_coords(binary, smi, s)
            if c is not None:
                sci_list.append(c)

        if not sci_list:
            results.append({'smiles': smi, 'name': entry['name'], 'cat': cat, 'rdkit_ok': True, 'sci_ok': False})
            continue
        categories[cat]['sci_ok'] += 1

        # Min RMSD across seed pairs
        min_rmsd = float('inf')
        for rc in rdkit_list:
            for sc in sci_list:
                if len(rc) == len(sc):
                    r = heavy_pairwise_rmsd(rc, sc)
                    min_rmsd = min(min_rmsd, r)

        if min_rmsd < float('inf'):
            categories[cat]['compared'] += 1
            categories[cat]['rmsds'].append(min_rmsd)
            results.append({
                'smiles': smi, 'name': entry['name'], 'cat': cat,
                'rdkit_ok': True, 'sci_ok': True, 'rmsd': min_rmsd,
            })
        else:
            results.append({
                'smiles': smi, 'name': entry['name'], 'cat': cat,
                'rdkit_ok': True, 'sci_ok': True, 'atom_mismatch': True,
            })

        if (i + 1) % 20 == 0:
            elapsed = time.time() - t0
            print(f"  Progress: {i + 1}/{len(entries)} ({(i + 1) / elapsed:.1f} mol/s)")

    elapsed = time.time() - t0
    all_rmsds = [r['rmsd'] for r in results if 'rmsd' in r]

    print(f"\n{'=' * 65}")
    print(f"DIVERSE MOLECULE BENCHMARK ({len(entries)} molecules, {n_seeds} seeds, {elapsed:.1f}s)")
    print(f"  Heavy-atom pairwise-distance RMSD (min over {n_seeds}×{n_seeds} seed pairs)")
    print(f"{'=' * 65}")

    total_rdkit = sum(1 for r in results if r.get('rdkit_ok'))
    total_sci = sum(1 for r in results if r.get('sci_ok'))
    print(f"  RDKit embed: {total_rdkit}/{len(entries)}")
    print(f"  sci-form embed: {total_sci}/{len(entries)}")
    print(f"  Compared: {len(all_rmsds)}/{len(entries)}")

    if all_rmsds:
        arr = np.array(all_rmsds)
        print(f"\n--- Overall RMSD ---")
        print(f"  Avg: {arr.mean():.4f} Å")
        print(f"  Median: {np.median(arr):.4f} Å")
        for t in [0.1, 0.3, 0.5, 1.0]:
            c = np.sum(arr < t)
            print(f"  < {t} Å: {c}/{len(arr)} ({100 * c / len(arr):.1f}%)")
        print(f"  Max: {arr.max():.4f} Å")

    # Per-category
    print(f"\n{'─' * 65}")
    print(f"{'Category':<16} {'N':>4} {'RDKit':>5} {'Sci':>5} {'Cmp':>4} {'Avg':>7} {'Med':>7} {'<0.5':>6}")
    print(f"{'─' * 65}")
    for cat in sorted(categories.keys()):
        s = categories[cat]
        if s['rmsds']:
            a = np.array(s['rmsds'])
            below = np.sum(a < 0.5)
            print(f"{cat:<16} {s['total']:>4} {s['rdkit_ok']:>5} {s['sci_ok']:>5} "
                  f"{s['compared']:>4} {a.mean():>7.4f} {np.median(a):>7.4f} "
                  f"{below:>3}/{s['compared']}")
        else:
            print(f"{cat:<16} {s['total']:>4} {s['rdkit_ok']:>5} {s['sci_ok']:>5} "
                  f"{s['compared']:>4}       -       -     -")

    # Worst molecules
    worst = sorted([r for r in results if 'rmsd' in r], key=lambda r: -r['rmsd'])[:15]
    if worst:
        print(f"\n--- Worst RMSD molecules ---")
        for r in worst:
            print(f"  {r['rmsd']:.3f} Å  [{r['cat']}] {r['name']} ({r['smiles'][:50]})")

    print()


if __name__ == '__main__':
    main()
