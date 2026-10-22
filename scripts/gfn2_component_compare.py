#!/usr/bin/env python3
"""Extract tblite energy components via numerical comparison.

Environment:
    SCIFORM_CLI    path to the sci-form binary
                   (default: <repo>/target/release/sci-form)
    SCIFORM_VENV   virtualenv holding the reference python packages
                   (default: the active interpreter)
    GFN2_LD_PRELOAD  libgfortran to preload for the tblite wheel
                   (default: /usr/lib/x86_64-linux-gnu/libgfortran.so.5)
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

def our_gfn2(elements, coords_flat):
    """Run our GFN2 on flat coords (Å)."""
    elems_str = json.dumps(elements)
    coords_str = json.dumps(coords_flat)
    env = os.environ.copy()
    env['GFN2_DEBUG'] = '1'
    r = subprocess.run([CLI, 'gfn2', elems_str, coords_str],
                       capture_output=True, text=True, env=env)
    if r.returncode != 0:
        raise RuntimeError(f"gfn2 failed: {r.stderr}")
    result = json.loads(r.stdout)
    
    # Parse debug for energy breakdown
    debug = {}
    for line in r.stderr.split('\n'):
        line = line.strip()
        if 'SCC iter' in line and 'e_total=' in line:
            parts = line.split()
            for p in parts:
                for key in ['e_total', 'e_band', 'e_scc', 'e_3rd', 'e_aes']:
                    if p.startswith(f'{key}='):
                        debug[key] = float(p.split('=')[1])
        if 'AES breakdown' in line:
            parts = line.split()
            for p in parts:
                for key in ['e_sd', 'e_dd', 'e_sq', 'e_dk', 'e_qk']:
                    if p.startswith(f'{key}='):
                        debug[f'aes_{key}'] = float(p.split('=')[1])
    
    return result, debug

def tblite_gfn2(elements, positions_ang):
    """Run tblite GFN2. Returns dict with energy components."""
    numbers = np.array(elements, dtype=np.int32)
    positions = np.array(positions_ang) * AATOAU
    calc = Calculator("GFN2-xTB", numbers, positions)
    res = calc.singlepoint(Result())
    
    out = {'total': res.get("energy")}
    try: out['charges'] = list(res.get("charges"))
    except: pass
    try: out['dipole'] = list(res.get("dipole"))
    except: pass
    
    return out

def compare_molecule(name, elements, coords_flat):
    """Compare all energy components for a molecule."""
    positions = np.array(coords_flat).reshape(-1, 3).tolist()
    
    result, debug = our_gfn2(elements, coords_flat)
    tbl = tblite_gfn2(elements, positions)
    
    our_ha = result['total_energy'] * EV2HA
    ref_ha = tbl['total']
    diff = our_ha - ref_ha
    pct = abs(diff / ref_ha) * 100 if ref_ha != 0 else 0
    
    print(f"\n{'='*70}")
    print(f"  {name}  (Z={elements})")
    print(f"{'='*70}")
    print(f"  Total:    ours={our_ha:+.10f}  ref={ref_ha:+.10f}  diff={diff:+.10f} Ha ({pct:.5f}%)")
    print(f"  Components (Ha):")
    print(f"    E_elec:  {result['electronic_energy']*EV2HA:+.10f}")
    print(f"    E_rep:   {result['repulsive_energy']*EV2HA:+.10f}")
    print(f"    E_disp:  {result['dispersion_energy']*EV2HA:+.10f}")
    print(f"    E_xb:    {result['halogen_bond_energy']*EV2HA:+.10f}")
    if debug:
        print(f"  SCC breakdown (Ha):")
        print(f"    e_band:  {debug.get('e_band', 0):+.10f}")
        print(f"    e_scc:   {debug.get('e_scc', 0):+.10f}")
        print(f"    e_3rd:   {debug.get('e_3rd', 0):+.10f}")
        print(f"    e_aes:   {debug.get('e_aes', 0):+.10f}")
        if 'aes_e_sd' in debug:
            print(f"  AES breakdown (Ha):")
            print(f"    e_sd:    {debug.get('aes_e_sd', 0):+.10f}")
            print(f"    e_dd:    {debug.get('aes_e_dd', 0):+.10f}")
            print(f"    e_sq:    {debug.get('aes_e_sq', 0):+.10f}")
            print(f"    e_dk:    {debug.get('aes_e_dk', 0):+.10f}")
            print(f"    e_qk:    {debug.get('aes_e_qk', 0):+.10f}")
    if 'charges' in tbl:
        print(f"  Charges: ours={[f'{q:.6f}' for q in result['mulliken_charges']]}")
        print(f"           ref ={[f'{q:.6f}' for q in tbl['charges']]}")
    
    return diff, pct

# ──────────────────────────────────────────────────────
# Test molecules with hand-specified geometries
# ──────────────────────────────────────────────────────

tests = []

# H2
tests.append(("H2", [1, 1], [0.0, 0.0, 0.0, 0.74, 0.0, 0.0]))

# Water
tests.append(("H2O", [8, 1, 1], 
    [0.0, 0.0, 0.0, 0.757, 0.586, 0.0, -0.757, 0.586, 0.0]))

# N2
tests.append(("N2", [7, 7], [0.0, 0.0, 0.0, 1.098, 0.0, 0.0]))

# CO
tests.append(("CO", [6, 8], [0.0, 0.0, 0.0, 1.128, 0.0, 0.0]))

# HF
tests.append(("HF", [9, 1], [0.0, 0.0, 0.0, 0.917, 0.0, 0.0]))

# Use sci-form CLI to embed some molecules for full comparison
def embed(smiles, seed=42):
    r = subprocess.run([CLI, 'embed', smiles, '-s', str(seed)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return None, None
    data = json.loads(r.stdout)
    return data['elements'], data['coords']

for smiles, name in [("C", "methane"), ("CC", "ethane"), ("c1ccccc1", "benzene")]:
    elems, coords = embed(smiles)
    if elems:
        tests.append((name, elems, coords))

print("GFN2-xTB Component Comparison")
print("=" * 70)

results = []
for name, elems, coords in tests:
    try:
        diff, pct = compare_molecule(name, elems, coords)
        results.append((name, diff, pct, len(elems)))
    except Exception as e:
        print(f"\nERROR for {name}: {e}")

print(f"\n\n{'='*70}")
print(f"SUMMARY")
print(f"{'='*70}")
print(f"{'Molecule':15s} {'N':>3s} {'diff(Ha)':>14s} {'error%':>10s} {'diff/N²':>14s}")
for name, diff, pct, n in results:
    n_basis_approx = n * 2  # rough
    print(f"{name:15s} {n:3d} {diff:+14.10f} {pct:10.6f}% {diff/n**2:+14.10f}")
