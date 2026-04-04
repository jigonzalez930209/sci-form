#!/usr/bin/env python3
"""Compare GFN2 energy breakdown: our Rust impl vs tblite reference."""
import subprocess, json, sys, os
import numpy as np

# Make sure tblite Python works
sys.path.insert(0, '/home/lestad/github/sci-form/venv/lib/python3.13/site-packages')
os.environ['LD_PRELOAD'] = '/usr/lib/x86_64-linux-gnu/libgfortran.so.5'

from tblite.interface import Calculator, Structure

CLI = '/home/lestad/github/sci-form/target/release/sci-form'
AATOAU = 1.8897259886

def embed(smiles, seed=42):
    """Use sci-form to embed SMILES → 3D coords."""
    r = subprocess.run([CLI, 'embed', smiles, '-s', str(seed)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"embed failed: {r.stderr}")
    data = json.loads(r.stdout)
    return data['elements'], data['coords']

def our_gfn2(elements, coords_flat):
    """Run our GFN2 CLI with debug output to get energy breakdown."""
    elems_str = json.dumps(elements)
    coords_str = json.dumps(coords_flat)
    env = os.environ.copy()
    env['GFN2_DEBUG'] = '1'
    r = subprocess.run([CLI, 'gfn2', elems_str, coords_str],
                       capture_output=True, text=True, env=env)
    if r.returncode != 0:
        raise RuntimeError(f"gfn2 failed: {r.stderr}")
    result = json.loads(r.stdout)
    
    # Parse debug output for energy breakdown
    debug = {}
    for line in r.stderr.split('\n'):
        if 'SCC iter' in line and 'e_total=' in line:
            # Last SCC iteration has the final values
            parts = line.split()
            for p in parts:
                if p.startswith('e_total='):
                    debug['e_total_ha'] = float(p.split('=')[1])
                elif p.startswith('e_band='):
                    debug['e_band_ha'] = float(p.split('=')[1])
                elif p.startswith('e_scc='):
                    debug['e_scc_ha'] = float(p.split('=')[1])
                elif p.startswith('e_3rd='):
                    debug['e_3rd_ha'] = float(p.split('=')[1])
                elif p.startswith('e_aes='):
                    debug['e_aes_ha'] = float(p.split('=')[1])
    
    return result, debug

def tblite_gfn2(elements, coords_flat):
    """Run tblite GFN2 to get reference energy."""
    numbers = np.array(elements, dtype=np.int32)
    positions = np.array(coords_flat).reshape(-1, 3) * AATOAU  # Å → bohr
    calc = Calculator("GFN2-xTB", numbers, positions)
    from tblite.interface import Result
    res = calc.singlepoint(Result())
    energy_ha = res.get("energy")
    
    # Try to get charges
    try:
        charges = res.get("charges")
    except:
        charges = None
    
    return energy_ha, charges

def compare(smiles, label=""):
    """Compare our GFN2 with tblite for a given SMILES."""
    elements, coords = embed(smiles)
    n_atoms = len(elements)
    n_heavy = sum(1 for z in elements if z > 1)
    
    result, debug = our_gfn2(elements, coords)
    ref_ha, ref_charges = tblite_gfn2(elements, coords)
    
    # Our total in Hartree
    our_total_ev = result['total_energy']
    our_elec_ev = result['electronic_energy']
    our_rep_ev = result['repulsive_energy']
    our_disp_ev = result['dispersion_energy']
    our_xb_ev = result['halogen_bond_energy']
    
    # Convert to Hartree for comparison
    EV2HA = 1.0 / 27.21138505
    our_total_ha = our_total_ev * EV2HA
    our_elec_ha = our_elec_ev * EV2HA
    our_rep_ha = our_rep_ev * EV2HA
    our_disp_ha = our_disp_ev * EV2HA
    our_xb_ha = our_xb_ev * EV2HA
    
    diff_ha = our_total_ha - ref_ha
    pct = abs(diff_ha / ref_ha) * 100 if ref_ha != 0 else 0
    
    print(f"\n{'='*60}")
    print(f"  {label or smiles}  ({n_atoms} atoms, {n_heavy} heavy)")
    print(f"{'='*60}")
    print(f"  tblite total:     {ref_ha:+.12f} Ha")
    print(f"  ours   total:     {our_total_ha:+.12f} Ha")
    print(f"  diff:             {diff_ha:+.12f} Ha  ({pct:.6f}%)")
    print(f"  diff per atom:    {diff_ha/n_atoms:+.12f} Ha")
    print(f"  diff per heavy:   {diff_ha/n_heavy:+.12f} Ha")
    print()
    print(f"  Our breakdown (Ha):")
    print(f"    E_elec:  {our_elec_ha:+.12f}")
    print(f"    E_rep:   {our_rep_ha:+.12f}")
    print(f"    E_disp:  {our_disp_ha:+.12f}")
    print(f"    E_xb:    {our_xb_ha:+.12f}")
    print(f"    Sum:     {our_elec_ha+our_rep_ha+our_disp_ha+our_xb_ha:+.12f}")
    
    if debug:
        print(f"\n  SCC breakdown (Ha, from debug):")
        for k, v in sorted(debug.items()):
            print(f"    {k:15s}: {v:+.12f}")
    
    # Compare charges if available
    if ref_charges is not None and 'mulliken_charges' in result:
        max_dq = max(abs(a - b) for a, b in zip(result['mulliken_charges'], ref_charges))
        print(f"\n  Max |Δq|: {max_dq:.8f}")
    
    return diff_ha, pct

# Test molecules: alkane series (clear size scaling)
molecules = [
    ("C", "methane"),
    ("CC", "ethane"),
    ("CCC", "propane"),
    ("CCCC", "butane"),
    ("CCCCC", "pentane"),
    ("CCCCCC", "hexane"),
    ("O", "water"),
    ("CO", "methanol"),
    ("CCO", "ethanol"),
    ("CC=O", "acetaldehyde"),
    ("c1ccccc1", "benzene"),
    ("c1ccncc1", "pyridine"),
]

print("GFN2-xTB Energy Comparison: Our Implementation vs tblite")
print("=========================================================")

results = []
for smiles, name in molecules:
    try:
        diff_ha, pct = compare(smiles, name)
        results.append((name, diff_ha, pct))
    except Exception as e:
        print(f"\nERROR for {name}: {e}")

print("\n\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"{'Molecule':20s} {'diff(Ha)':>14s} {'error%':>10s}")
for name, diff, pct in results:
    print(f"{name:20s} {diff:+14.10f} {pct:10.6f}%")
