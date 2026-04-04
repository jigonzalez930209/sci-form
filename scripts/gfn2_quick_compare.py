#!/usr/bin/env python3
"""Quick comparison before/after Fermi smearing."""
import os, sys, json, subprocess
import numpy as np

sys.path.insert(0, '/home/lestad/github/sci-form/venv/lib/python3.12/site-packages')
os.environ['LD_PRELOAD'] = '/usr/lib/x86_64-linux-gnu/libgfortran.so.5'

from tblite.interface import Calculator, Result

AATOAU = 1.8897259886
CLI = '/home/lestad/github/sci-form/target/release/sci-form'
EV2HA = 1.0 / 27.21138505

def embed(smiles, seed=42):
    r = subprocess.run([CLI, 'embed', smiles, '-s', str(seed)], capture_output=True, text=True)
    data = json.loads(r.stdout)
    return data['elements'], data['coords']

def our_gfn2(elements, coords_flat):
    r = subprocess.run([CLI, 'gfn2', json.dumps(elements), json.dumps(coords_flat)],
                       capture_output=True, text=True)
    return json.loads(r.stdout)

def tblite_gfn2(elements, positions_ang):
    numbers = np.array(elements, dtype=np.int32)
    positions = np.array(positions_ang) * AATOAU
    calc = Calculator("GFN2-xTB", numbers, positions)
    res = calc.singlepoint(Result())
    return res.get("energy")

tests = [
    ("H2", [1,1], [0,0,0,0.74,0,0]),
    ("H2O", [8,1,1], [0,0,0, 0.757,0.586,0, -0.757,0.586,0]),
]

for smiles in ["C", "CC", "CCC", "CCCC", "CCCCC", "CCCCCC", "c1ccccc1", "c1ccncc1", "CCO", "CC=O"]:
    elems, coords = embed(smiles)
    tests.append((smiles, elems, coords))

print(f"{'Molecule':15s} {'N_atoms':>7s} {'diff(Ha)':>14s} {'error%':>10s}")
print("-" * 50)

for name, elems, coords in tests:
    result = our_gfn2(elems, coords)
    our_ha = result['total_energy'] * EV2HA
    positions = np.array(coords).reshape(-1, 3).tolist()
    ref_ha = tblite_gfn2(elems, positions)
    diff = our_ha - ref_ha
    pct = abs(diff / ref_ha) * 100
    n = len(elems)
    print(f"{name:15s} {n:7d} {diff:+14.10f} {pct:10.6f}%")
