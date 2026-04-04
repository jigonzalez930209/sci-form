#!/usr/bin/env python3
"""
Direct matrix comparison: extract H0 and overlap from both tblite and our code for water.
waterpos in bohr: O at origin, H at known positions.
"""
import os, sys, json, subprocess
import numpy as np

sys.path.insert(0, '/home/lestad/github/sci-form/venv/lib/python3.12/site-packages')
os.environ['LD_PRELOAD'] = '/usr/lib/x86_64-linux-gnu/libgfortran.so.5'

from tblite.interface import Calculator, Result

AATOAU = 1.8897259886
CLI = '/home/lestad/github/sci-form/target/release/sci-form'

# Water geometry in Angstrom
elements = [8, 1, 1]
coords_ang = [0.0, 0.0, 0.0,  0.757, 0.586, 0.0,  -0.757, 0.586, 0.0]
positions_ang = np.array(coords_ang).reshape(-1, 3)

# ━━━ tblite ━━━
numbers = np.array(elements, dtype=np.int32)
positions_bohr = positions_ang * AATOAU
calc = Calculator("GFN2-xTB", numbers, positions_bohr)
calc.set("save-integrals", 1)
res = calc.singlepoint(Result())

tbl_energy = res.get("energy")
tbl_S = res.get("overlap-matrix")
tbl_H0 = res.get("hamiltonian-matrix")
tbl_charges = res.get("charges")
tbl_orb_e = res.get("orbital-energies")
tbl_occs = res.get("orbital-occupations")
tbl_density = res.get("density-matrix")

print(f"tblite total energy: {tbl_energy:.12f} Ha")
print(f"tblite charges: {tbl_charges}")
print(f"tblite orbital energies: {tbl_orb_e}")
print(f"tblite occupations: {tbl_occs}")
print()

# ━━━ Our code ━━━
env = os.environ.copy()
env['GFN2_DEBUG'] = '1'
r = subprocess.run([CLI, 'gfn2', json.dumps(elements), json.dumps(coords_ang)],
                   capture_output=True, text=True, env=env)
result = json.loads(r.stdout)
our_total_ev = result['total_energy']
our_total_ha = our_total_ev / 27.21138505

print(f"Our total energy: {our_total_ha:.12f} Ha")
print(f"Our charges: {result['mulliken_charges']}")
print(f"diff: {our_total_ha - tbl_energy:+.12f} Ha")
print()

# Parse our debug output for matrices
stderr = r.stderr
lines = stderr.split('\n')

# Parse overlap and H0 from debug
our_S = {}
our_H0_diag = {}
our_H0_off = {}

for line in lines:
    line = line.strip()
    # Parse S matrix: "S[i]: [v0, v1, ...]"
    if line.startswith('S['):
        idx = int(line.split('[')[1].split(']')[0])
        vals_str = line.split(': [')[1].rstrip(']')
        vals = [float(x) for x in vals_str.split(', ')]
        our_S[idx] = vals
    # Parse H0 diagonal: "H0[i,i] = val (atom=, sh=, l=, m=)"
    elif line.startswith('H0[') and ',' in line.split(']')[0]:
        parts = line.split()
        idx_str = parts[0][3:-1]  # "i,i" 
        i, j = idx_str.split(',')
        i, j = int(i), int(j)
        val = float(parts[2])
        if i == j:
            our_H0_diag[i] = val
        else:
            our_H0_off[(i, j)] = val

n_orbs = len(tbl_S)
print(f"Number of orbitals: {n_orbs}")
print()

# Compare overlap matrices
print("=== Overlap Matrix Comparison ===")
max_s_diff = 0
for i in range(n_orbs):
    for j in range(n_orbs):
        tbl_val = tbl_S[i, j]
        our_val = our_S.get(i, [0]*n_orbs)[j] if i in our_S else 0
        diff = abs(tbl_val - our_val)
        if diff > max_s_diff:
            max_s_diff = diff
        if diff > 1e-8:
            print(f"  S[{i},{j}]: tblite={tbl_val:.12f}  ours={our_val:.12f}  diff={diff:.2e}")
print(f"  Max S difference: {max_s_diff:.2e}")
print()

# Compare H0 diagonal
print("=== H0 Diagonal Comparison ===")
max_h0_diff = 0
for i in range(n_orbs):
    tbl_val = tbl_H0[i, i]
    our_val = our_H0_diag.get(i, 0)
    diff = abs(tbl_val - our_val)
    if diff > max_h0_diff:
        max_h0_diff = diff
    print(f"  H0[{i},{i}]: tblite={tbl_val:.12f}  ours={our_val:.12f}  diff={diff:.2e}")
print(f"  Max H0 diagonal difference: {max_h0_diff:.2e}")
print()

# Compare H0 off-diagonal
print("=== H0 Off-Diagonal Comparison (largest diffs) ===")
h0_diffs = []
for i in range(n_orbs):
    for j in range(i+1, n_orbs):
        tbl_val = tbl_H0[i, j]
        our_val = our_H0_off.get((i, j), 0)
        diff = abs(tbl_val - our_val)
        h0_diffs.append((i, j, tbl_val, our_val, diff))

h0_diffs.sort(key=lambda x: -x[4])
for i, j, tv, ov, d in h0_diffs[:20]:
    if d > 1e-10:
        print(f"  H0[{i},{j}]: tblite={tv:.12f}  ours={ov:.12f}  diff={d:.2e}")
if h0_diffs:
    print(f"  Max H0 off-diagonal difference: {h0_diffs[0][4]:.2e}")
print()

# Compare orbital energies
print("=== Orbital Energies Comparison (H0-only eigenvalues) ===")
# Note: tblite H0 eigenvalues would be from S^{-1/2} H0 S^{-1/2}
# Let's compute them
from scipy.linalg import eigh
evals_tbl, _ = eigh(tbl_H0, tbl_S)
print(f"  tblite H0 eigenvalues: {evals_tbl}")

# Parse our H0 eigenvalues from debug
our_h0_evals = []
for line in lines:
    if 'ε_H0[' in line:
        val = float(line.split('=')[1].split('(')[0].strip())
        our_h0_evals.append(val)

if our_h0_evals:
    our_h0_evals = sorted(our_h0_evals)
    print(f"  Our H0 eigenvalues:    {our_h0_evals}")
    print()
    print("  Per-orbital diff:")
    for i in range(min(len(evals_tbl), len(our_h0_evals))):
        d = our_h0_evals[i] - evals_tbl[i]
        print(f"    ε[{i}]: tblite={evals_tbl[i]:.12f}  ours={our_h0_evals[i]:.12f}  diff={d:+.2e}")

# Compare converged orbital energies
print()
print("=== Converged Orbital Energies ===")
our_orb_ha = [e / 27.21138505 for e in result['orbital_energies']]
for i in range(n_orbs):
    d = our_orb_ha[i] - tbl_orb_e[i]
    print(f"  ε[{i}]: tblite={tbl_orb_e[i]:.10f}  ours={our_orb_ha[i]:.10f}  diff={d:+.2e}")

# Compare density matrix
print()
print("=== Density Matrix Comparison ===")
if tbl_density is not None:
    max_p_diff = 0
    for i in range(n_orbs):
        for j in range(n_orbs):
            diff = abs(tbl_density[i, j])
            # We don't have our density exported, but we can check tblite's
    print(f"  tblite P trace: {sum(tbl_density[i,i] for i in range(n_orbs)):.10f}")
    print(f"  (Expected: {sum(tbl_occs):.1f} electrons)")
