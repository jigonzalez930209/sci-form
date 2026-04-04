#!/usr/bin/env python3
"""Verify H2 energy components to isolate 0.05 mHa discrepancy."""
import math, json, subprocess, os

CLI = '/home/lestad/github/sci-form/target/release/sci-form'
ANGSTROM_TO_BOHR = 1.0 / 0.529177  # solver.rs

# GFN2 H parameters
zeta_H = 1.230
selfenergy_H = -10.707211  # eV
kcn_H = -0.050000
shpoly_H = -0.00953618
EV_TO_HARTREE = 1.0 / 27.21138505  # gfn2_params.rs
HUBBARD_H = 0.405771
SHELL_HUBBARD_H = 1.0
DKERNEL_H = 0.05563889
QKERNEL_H = 0.00027431
REP_ALPHA_H = 2.213717
REP_ZEFF_H = 1.105388
REP_KEXP_LIGHT = 1.0
REP_REXP = 1.0
ATOMIC_RAD_H = 0.32 * 1.8897259886  # bohr
PAULING_EN_H = 2.20

# H2 geometry
bond_ang = 0.74
R_bohr = bond_ang * ANGSTROM_TO_BOHR
print(f"H2 at R = {bond_ang} Å = {R_bohr:.12f} bohr")

# 1. Repulsion energy
alpha_pair = math.sqrt(REP_ALPHA_H * REP_ALPHA_H)  # = REP_ALPHA_H
zeff_pair = REP_ZEFF_H * REP_ZEFF_H
kexp = REP_KEXP_LIGHT  # both atoms Z=1 ≤ 2
rexp = REP_REXP
e_rep = zeff_pair * math.exp(-alpha_pair * R_bohr**kexp) / R_bohr**rexp
print(f"\nRepulsion:")
print(f"  alpha_pair = {alpha_pair:.10f}")
print(f"  zeff_pair  = {zeff_pair:.10f}")
print(f"  R^kexp     = {R_bohr**kexp:.12f}")
print(f"  exp term   = {math.exp(-alpha_pair * R_bohr**kexp):.12f}")
print(f"  E_rep      = {e_rep:.12f} Ha")
print(f"  Rust E_rep = 0.039534125669 Ha")  # from our debug output

# 2. Band energy (H0 eigenvalues)
# For H2: n_occ = 1, E_band = 2 * epsilon_1
# epsilon_1 from Rust debug: -0.5189927396 Ha
e_band_rust = 2 * (-0.5189927396)
print(f"\nBand energy:")
print(f"  E_band = 2 * (-0.5189927396) = {e_band_rust:.12f}")
print(f"  Rust debug e_band = -1.037985479133")

# 3. AES kernel energy
# For H2, charges are ~0, so only kernel terms contribute
# E_AES = E_kernel_d + E_kernel_q
# E_kernel_d = Σ_A dkernel * |d_A|²
# E_kernel_q = Σ_A qkernel * Σ_c (q_A[c]² * scale[c])
# where d_A and q_A are Mulliken atomic multipoles
# For H2 this requires knowing the density matrix and multipole integrals

# From the Rust debug, E_AES = 0.016520 Ha for H2
# Since charges are ~0, this is entirely kernel energy
# E_AES = 2 * dkernel * d_x² + 2 * qkernel * Σ q²*scale
# (two H atoms)

# Let me check what d_x should be for H2
# d_A[k] = -Σ_μ_on_A Σ_ν P(ν,μ) * dpint(k,μ,ν)
# For H2 with P[0,0]=P[1,1]=P_diag, P[0,1]=P[1,0]=P_off
# d_A[x] = -P[0,0]*dpint_x(0,0) - P[1,0]*dpint_x(0,1)
# dpint_x(0,0) = 0 (same center, odd integrand)
# So d_A[x] = -P_off * dpint_x(0,1)

# S[0,1] = 0.6638
# For occupied orbital: |ψ⟩ = (|0⟩ + |1⟩) / sqrt(2*(1+S))
# P[0,0] = P[1,1] = 2/(2*(1+S)) = 1/(1+S)
# P[0,1] = P[1,0] = 1/(1+S)
S01 = 0.6637791614
P_diag = 1.0 / (1.0 + S01)
P_off = 1.0 / (1.0 + S01)
print(f"\nDensity matrix:")
print(f"  P[0,0] = P[1,1] = {P_diag:.12f}")
print(f"  P[0,1] = P[1,0] = {P_off:.12f}")

# 4. D4 dispersion
# Our code: -0.000002 Ha for H2
# tblite's D4 may differ slightly
# Let's check typical C6 for H-H
print(f"\nD4 dispersion: our value = -0.000002 Ha (very small)")

# 5. Summary of our energy components
e_elec = -1.021465936232  # from JSON
e_rep_json = 0.039534125669
e_disp = -0.000002315860
total_ours = e_elec + e_rep_json + e_disp
print(f"\nOur energy breakdown:")
print(f"  E_elec = {e_elec:.12f} Ha (E_band + E_scc + E_3rd + E_aes)")
print(f"  E_rep  = {e_rep_json:.12f} Ha")
print(f"  E_disp = {e_disp:.12f} Ha")
print(f"  Total  = {total_ours:.12f} Ha")
print(f"  tblite = -0.981983692538 Ha")
print(f"  diff   = {total_ours - (-0.981983692538):+.12f} Ha")

# Check: is our D4 energy reasonable?
# For H2 at ~1.4 bohr, with C6(H-H) ~ 4 Ha*bohr^6 (typical):
# E_D4 ~ -s6 * C6 / R^6 * fdamp ~ -1.0 * 4 / 7.5 * small_damp ≈ very small
# So -0.000002 Ha seems too small but plausible for short-range damped dispersion

# Now let me try WITH D4 = 0 (no dispersion) to see how much it matters:
e_no_disp = total_ours - e_disp
tblite_no_disp = -0.981983692538 + 0.000002  # rough estimate
print(f"\nIf we ignore dispersion:")
print(f"  Our (no disp) = {e_no_disp:.12f}")
print(f"  Error still ≈ {total_ours - (-0.981983692538):+.6f} Ha")
print(f"  D4 contribution to error << error")

# 6. Let me check if the issue is in how tblite handles rounding
# The key question: what's the PRECISE tblite value for each component?
# Without access to tblite internals, let me check D4 parameters
print("\n\n=== D4 Parameters Check ===")

# Our D4 parameters for GFN2-xTB:
# s6 = 1.0, s8 = 2.7, a1 = 0.52, a2 = 5.0
# Check against the official GFN2-xTB paper / tblite
env = os.environ.copy()
env['GFN2_DEBUG'] = '1'
r = subprocess.run([CLI, 'gfn2', '[1,1]', '[0,0,0,0.74,0,0]'],
                   capture_output=True, text=True, env=env)
# Find D4 debug info
for line in r.stderr.split('\n'):
    if 'D4' in line or 'disp' in line.lower() or 'd4' in line:
        print(f"  {line.strip()}")
    if 'e_total' in line or 'E_rep' in line or 'E_elec' in line:
        print(f"  {line.strip()}")

# Parse JSON result for precise values
result = json.loads(r.stdout)
EV2HA = 1.0 / 27.21138505
print(f"\nJSON output (converted to Ha):")
print(f"  electronic_energy: {result['electronic_energy'] * EV2HA:.12f}")
print(f"  repulsive_energy:  {result['repulsive_energy'] * EV2HA:.12f}")
print(f"  dispersion_energy: {result['dispersion_energy'] * EV2HA:.12f}")
print(f"  halogen_bond:      {result['halogen_bond_energy'] * EV2HA:.12f}")
print(f"  total_energy:      {result['total_energy'] * EV2HA:.12f}")
print(f"  homo_energy:       {result['homo_energy']:.8f} eV")
print(f"  lumo_energy:       {result['lumo_energy']:.8f} eV")
print(f"  gap:               {result['gap']:.8f} eV")
print(f"  charges:           {result['mulliken_charges']}")
print(f"  scc_iterations:    {result['scc_iterations']}")
