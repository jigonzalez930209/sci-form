#!/usr/bin/env python3
"""Quick comparison before/after Fermi smearing.

Environment: SCIFORM_CLI, SCIFORM_VENV, GFN2_LD_PRELOAD (see
`scripts/gfn2_component_compare.py`).
"""
import os, sys, json, subprocess
from pathlib import Path
import numpy as np

# The reference python packages (tblite) come from the active interpreter or
# from the virtualenv pointed at by SCIFORM_VENV — never from a hard-coded
# path that only exists on one machine.
_VENV = os.environ.get('SCIFORM_VENV')
if _VENV:
    sys.path.insert(0, str(Path(_VENV) / 'lib' /
                            f'python{sys.version_info.major}.{sys.version_info.minor}' /
                            'site-packages'))
# The tblite wheel links against libgfortran; only preload it when the
# distribution actually ships it there.
_PRELOAD = os.environ.get('GFN2_LD_PRELOAD', '/usr/lib/x86_64-linux-gnu/libgfortran.so.5')
if Path(_PRELOAD).exists():
    os.environ['LD_PRELOAD'] = _PRELOAD

from tblite.interface import Calculator, Result

AATOAU = 1.8897259886
CLI = os.environ.get('SCIFORM_CLI',
                    str(Path(__file__).resolve().parents[1] / 'target' / 'release' / 'sci-form'))
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
