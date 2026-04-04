#!/usr/bin/env python3
"""Compare GFN2 for H2 molecule: our impl vs tblite, with manual overlap verification."""
import subprocess, json, sys, os, math
import numpy as np

sys.path.insert(0, '/home/lestad/github/sci-form/venv/lib/python3.12/site-packages')
os.environ['LD_PRELOAD'] = '/usr/lib/x86_64-linux-gnu/libgfortran.so.5'
from tblite.interface import Calculator, Structure

CLI = '/home/lestad/github/sci-form/target/release/sci-form'
AATOAU = 1.8897259886
PI = math.pi

# ── STO-3G 1s coefficients (Stewart 1970, Table I) ──
STO3G_1S = [
    (2.227660584e+0, 1.543289673e-1),
    (4.057711562e-1, 5.353281423e-1),
    (1.098175104e-1, 4.446345422e-1),
]

def gaussian_norm_s(alpha):
    """Normalization for s-type Cartesian Gaussian."""
    return (2.0 * alpha / PI) ** 0.75

def overlap_1d_s(gamma):
    """Obara-Saika overlap for l=0: just sqrt(pi/gamma)."""
    return (PI / gamma) ** 0.5

def sto3g_overlap_ss(zeta_a, pos_a, zeta_b, pos_b):
    """Compute STO-3G 1s-1s overlap between two centers."""
    r2 = sum((a - b)**2 for a, b in zip(pos_a, pos_b))
    
    s = 0.0
    for ar, ca in STO3G_1S:
        alpha_a = ar * zeta_a**2
        na = gaussian_norm_s(alpha_a)
        for br, cb in STO3G_1S:
            alpha_b = br * zeta_b**2
            nb = gaussian_norm_s(alpha_b)
            gamma = alpha_a + alpha_b
            mu = alpha_a * alpha_b / gamma
            kab = math.exp(-mu * r2)
            sx = overlap_1d_s(gamma)  # sqrt(pi/gamma)
            # 3D overlap = sx * sy * sz = (pi/gamma)^{3/2}
            raw = (PI / gamma) ** 1.5
            s += ca * cb * na * nb * kab * raw
    return s

# ── Test H2 ──
# H₂ at 0.74 Å bond length (standard)
bond_angstrom = 0.74
h2_pos_angstrom = [[0.0, 0.0, 0.0], [bond_angstrom, 0.0, 0.0]]
h2_elements = [1, 1]

# H slater exponent in GFN2
zeta_H = 1.823  # GFN2 parameter for H 1s

# Manual STO-3G overlap (in bohr)
pos_a_bohr = [0.0, 0.0, 0.0]
pos_b_bohr = [bond_angstrom * AATOAU, 0.0, 0.0]
S_manual = sto3g_overlap_ss(zeta_H, pos_a_bohr, zeta_H, pos_b_bohr)
print(f"Manual STO-3G overlap S[0,1] for H₂ (ζ={zeta_H}, R={bond_angstrom} Å):")
print(f"  S[0,1] = {S_manual:.12f}")
print(f"  Bond distance in bohr: {bond_angstrom * AATOAU:.10f}")

# Also check with the solver.rs constant
ANGSTROM_TO_BOHR_SOLVER = 1.0 / 0.529177
pos_b_bohr_solver = [bond_angstrom * ANGSTROM_TO_BOHR_SOLVER, 0.0, 0.0]
S_manual_solver = sto3g_overlap_ss(zeta_H, pos_a_bohr, zeta_H, pos_b_bohr_solver)
print(f"\n  With solver.rs constant (1/0.529177 = {ANGSTROM_TO_BOHR_SOLVER:.12f}):")
print(f"  S[0,1] = {S_manual_solver:.12f}")
print(f"  diff = {S_manual_solver - S_manual:.3e}")

# ── Run our GFN2 CLI on H2 ──
print("\n--- Our GFN2 ---")
elems_str = json.dumps(h2_elements)
coords_flat = []
for p in h2_pos_angstrom:
    coords_flat.extend(p)
coords_str = json.dumps(coords_flat)
env = os.environ.copy()
env['GFN2_DEBUG'] = '1'
r = subprocess.run([CLI, 'gfn2', elems_str, coords_str],
                   capture_output=True, text=True, env=env)
result = json.loads(r.stdout)

# Parse overlap from debug
for line in r.stderr.split('\n'):
    if 'S[0]:' in line or 'S[1]:' in line:
        print(f"  {line.strip()}")
    if 'H0[' in line:
        print(f"  {line.strip()}")
    if 'ε_H0' in line:
        print(f"  {line.strip()}")
    if 'SCC iter' in line and 'e_total=' in line:
        print(f"  {line.strip()}")
    if 'CN:' in line:
        print(f"  {line.strip()}")

EV2HA = 1.0 / 27.21138505
our_ha = result['total_energy'] * EV2HA
print(f"\n  Total energy: {our_ha:.12f} Ha ({result['total_energy']:.8f} eV)")
print(f"  E_elec:       {result['electronic_energy'] * EV2HA:.12f} Ha")
print(f"  E_rep:        {result['repulsive_energy'] * EV2HA:.12f} Ha")
print(f"  E_disp:       {result['dispersion_energy'] * EV2HA:.12f} Ha")

# ── Run tblite on H2 ──
print("\n--- tblite GFN2 ---")
numbers = np.array(h2_elements, dtype=np.int32)
positions = np.array(h2_pos_angstrom) * AATOAU  # Å → bohr
calc = Calculator("GFN2-xTB", numbers, positions)
from tblite.interface import Result
res = calc.singlepoint(Result())
ref_ha = res.get("energy")
try:
    charges = res.get("charges")
    print(f"  Charges: {charges}")
except:
    charges = None
print(f"  Total energy: {ref_ha:.12f} Ha")

# ── Compare ──
diff = our_ha - ref_ha
pct = abs(diff / ref_ha) * 100
print(f"\n--- Comparison ---")
print(f"  Diff: {diff:+.12f} Ha ({pct:.6f}%)")
print(f"  Our charges: {result.get('mulliken_charges')}")

# Now also compare for water
print("\n\n========== Water ==========")
water_elements = [8, 1, 1]
water_pos = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]]

elems_str = json.dumps(water_elements)
coords_flat = []
for p in water_pos:
    coords_flat.extend(p)
coords_str = json.dumps(coords_flat)
r = subprocess.run([CLI, 'gfn2', elems_str, coords_str],
                   capture_output=True, text=True, env=env)
result_w = json.loads(r.stdout)
our_w_ha = result_w['total_energy'] * EV2HA

# Parse overlap
for line in r.stderr.split('\n'):
    if 'S[' in line and 'S[' in line and ']:' in line:
        print(f"  {line.strip()}")
    if line.strip().startswith('H0[') and '=' in line:
        print(f"  {line.strip()}")
    if 'SCC iter' in line and 'e_total=' in line:
        print(f"  {line.strip()}")

numbers_w = np.array(water_elements, dtype=np.int32)
positions_w = np.array(water_pos) * AATOAU
calc_w = Calculator("GFN2-xTB", numbers_w, positions_w)
res_w = calc_w.singlepoint(Result())
ref_w_ha = res_w.get("energy")

diff_w = our_w_ha - ref_w_ha
pct_w = abs(diff_w / ref_w_ha) * 100
print(f"\n  Our total:    {our_w_ha:.12f} Ha")
print(f"  tblite total: {ref_w_ha:.12f} Ha")  
print(f"  Diff: {diff_w:+.12f} Ha ({pct_w:.6f}%)")
print(f"  E_elec: {result_w['electronic_energy'] * EV2HA:.12f}")
print(f"  E_rep:  {result_w['repulsive_energy'] * EV2HA:.12f}")
print(f"  E_disp: {result_w['dispersion_energy'] * EV2HA:.12f}")
print(f"  Our charges: {result_w.get('mulliken_charges')}")
try:
    charges_w = res_w.get("charges")
    print(f"  tblite charges: {list(charges_w)}")
except:
    pass
