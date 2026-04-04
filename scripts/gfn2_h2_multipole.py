#!/usr/bin/env python3
"""Compute H2 multipole integrals and AES kernel energy, compare with Rust."""
import math

PI = math.pi

STO3G_1S = [
    (2.227660584e+0, 1.543289673e-1),
    (4.057711562e-1, 5.353281423e-1),
    (1.098175104e-1, 4.446345422e-1),
]

ANGSTROM_TO_BOHR = 1.0 / 0.529177
zeta_H = 1.230

bond_ang = 0.74
R = bond_ang * ANGSTROM_TO_BOHR  # bohr

def gaussian_norm_s(alpha):
    return (2.0 * alpha / PI) ** 0.75

# Compute overlap, dipole, and quadrupole integrals for H₁(0)-H₂(R) along x
def compute_multipole_h2():
    """Returns S, dpint[3], qpint[6] for H₁ at origin and H₂ at (R,0,0)."""
    S = 0.0
    dp = [0.0, 0.0, 0.0]  # x, y, z about A (atom of first basis)
    qp = [0.0]*6  # xx, xy, yy, xz, yz, zz about A
    
    for ar, ca in STO3G_1S:
        alpha_a = ar * zeta_H**2
        na = gaussian_norm_s(alpha_a)
        for br, cb in STO3G_1S:
            alpha_b = br * zeta_H**2
            nb = gaussian_norm_s(alpha_b)
            gamma = alpha_a + alpha_b
            mu = alpha_a * alpha_b / gamma
            kab = math.exp(-mu * R*R)
            
            # Gaussian product center
            # A at origin, B at (R,0,0)
            Px = alpha_b * R / gamma
            Py = 0.0
            Pz = 0.0
            PAx = Px  # P - A
            PAy = 0.0
            PAz = 0.0
            
            # Obara-Saika for s-s: S00 = sqrt(pi/gamma)
            s00 = (PI / gamma) ** 0.5
            
            # 3D overlap: s00_x * s00_y * s00_z = (pi/gamma)^{3/2}
            raw_s = s00 ** 3
            
            # Dipole moment about A: M_x = PA_x * S00
            # The dipole integral ⟨s_A|(r-R_A)|s_B⟩ for x-component:
            # = (PA_x * s00_x) * s00_y * s00_z
            # Because: ⟨0_A|(x-A_x)|0_B⟩_x = PA_x * s00_x
            # And y,z are just s00_y, s00_z
            dp_x_raw = PAx * raw_s  # = PAx * (pi/gamma)^{3/2}
            dp_y_raw = 0.0 * raw_s  # PAy = 0
            dp_z_raw = 0.0 * raw_s  # PAz = 0
            
            # Quadrupole about A: ⟨s_A|(r-R_A)²|s_B⟩
            # For xx component: ⟨0_A|(x-A_x)²|0_B⟩
            # From Obara-Saika: S(2,0) = PA * S(1,0) + 1/(2g) * S(0,0)
            # S(1,0)_x = PA_x * s00_x (as above)
            # S(2,0)_x = PA_x * (PA_x * s00_x) + 1/(2*gamma) * s00_x
            #          = (PA_x^2 + 1/(2*gamma)) * s00_x
            inv2g = 0.5 / gamma
            s10_x = PAx * s00
            s20_x = PAx * s10_x + inv2g * s00  # = (PA_x^2 + 1/2g) * s00
            
            # qp_xx = s20_x * s00_y * s00_z
            qp_xx_raw = s20_x * s00 * s00
            # qp_xy = s10_x * (PAy * s00) * s00 = 0 (PAy = 0)
            qp_xy_raw = 0.0
            # qp_yy = s00 * (PAy^2 + 1/2g) * s00 * s00 = (0 + 1/2g) * (pi/g)^{3/2}
            s20_y = 0.0 * 0.0 * s00 + inv2g * s00  # PAy=0, so just inv2g * s00
            qp_yy_raw = s00 * s20_y * s00
            # qp_xz = s10_x * s00 * (PAz * s00) = 0 (PAz = 0)
            qp_xz_raw = 0.0
            # qp_yz = 0
            qp_yz_raw = 0.0
            # qp_zz
            s20_z = inv2g * s00  # PAz = 0
            qp_zz_raw = s00 * s00 * s20_z
            
            coeff = ca * cb * na * nb * kab
            S += coeff * raw_s
            dp[0] += coeff * dp_x_raw
            dp[1] += coeff * dp_y_raw
            dp[2] += coeff * dp_z_raw
            qp[0] += coeff * qp_xx_raw  # xx
            qp[1] += coeff * qp_xy_raw  # xy
            qp[2] += coeff * qp_yy_raw  # yy
            qp[3] += coeff * qp_xz_raw  # xz
            qp[4] += coeff * qp_yz_raw  # yz
            qp[5] += coeff * qp_zz_raw  # zz
    
    return S, dp, qp

S01, dp01, qp01 = compute_multipole_h2()
print(f"H₂ multipole integrals (about center A):")
print(f"  S[0,1] = {S01:.15f}")
print(f"  dp_x[0,1] = {dp01[0]:.15f}")
print(f"  dp_y[0,1] = {dp01[1]:.15f}")
print(f"  dp_z[0,1] = {dp01[2]:.15f}")
print(f"  qp_xx[0,1] = {qp01[0]:.15f}")
print(f"  qp_yy[0,1] = {qp01[2]:.15f}")
print(f"  qp_zz[0,1] = {qp01[5]:.15f}")

# Diagonal multipole integrals (same center)
# ⟨s_A|(x-A)²|s_A⟩ for s-type: PA=0, so only 1/(2g) term
def compute_diagonal_qp():
    """Quadrupole diagonal: ⟨s_A|(r-R_A)²|s_A⟩."""
    qp = [0.0]*6
    for ar, ca in STO3G_1S:
        alpha_a = ar * zeta_H**2
        na = gaussian_norm_s(alpha_a)
        for br, cb in STO3G_1S:
            alpha_b = br * zeta_H**2
            nb = gaussian_norm_s(alpha_b)
            gamma = alpha_a + alpha_b
            inv2g = 0.5 / gamma
            s00 = (PI / gamma) ** 0.5
            kab = 1.0  # same center, R=0
            
            # qp_xx = (0 + 1/2g) * (pi/g)^{3/2}
            raw_xx = inv2g * s00**3
            
            coeff = ca * cb * na * nb * kab
            qp[0] += coeff * raw_xx  # xx
            qp[2] += coeff * raw_xx  # yy (same by symmetry)
            qp[5] += coeff * raw_xx  # zz (same by symmetry)
    return qp

qp_diag = compute_diagonal_qp()
print(f"\nDiagonal quadrupole (A on A):")
print(f"  qp_xx[0,0] = {qp_diag[0]:.15f}")
print(f"  qp_yy[0,0] = {qp_diag[2]:.15f}")
print(f"  qp_zz[0,0] = {qp_diag[5]:.15f}")

# Now compute Mulliken multipoles for H₂
# P matrix: P[0,0] = P[1,1] = 1/(1+S), P[0,1] = P[1,0] = 1/(1+S)
P_diag = 1.0 / (1.0 + S01)
P_off = 1.0 / (1.0 + S01)
print(f"\nDensity matrix: P_diag = P_off = {P_diag:.12f}")

# dpint about A: dpint_A(0,1) = dp_x about A
# dpint about B: dpint_B(1,0) = dp_x about A + (A_x-B_x) * S01 (shift formula)
# Since dpint[k][(j,i)] = about atom of j:
# dpint[k][(0,1)] = about atom 0 = dp01 (computed above)
# dpint[k][(1,0)] = about atom 1 = dp01 + (A-B) * S01
# For x: dpint_x[(1,0)] = dp01[0] + (0 - R) * S01 = dp01[0] - R * S01
dpint_x_10 = dp01[0] - R * S01  # about atom B for (1,0) pair

print(f"\nDipole integrals:")
print(f"  dpint_x[(0,1)] = {dp01[0]:.12f}  (about A)")
print(f"  dpint_x[(1,0)] = {dpint_x_10:.12f}  (about B)")

# Quadrupole about B: qp_B = qp_A + 2*(A-B)*dp_A + (A-B)^2 * S
# For xx: qp_xx_B(1,0) = qp_xx_A(0,1) + 2*(Ax-Bx)*dp_x_A(0,1) + (Ax-Bx)^2 * S01
qpint_xx_10 = qp01[0] + 2*(0 - R)*dp01[0] + (0-R)**2 * S01
qpint_yy_10 = qp01[2]  # y-component: (Ay-By)=0, so no shift
qpint_zz_10 = qp01[5]  # z-component: (Az-Bz)=0, so no shift

print(f"\nQuadrupole integrals:")
print(f"  qp_xx[(0,1)] = {qp01[0]:.12f}  (about A)")
print(f"  qp_xx[(1,0)] = {qpint_xx_10:.12f}  (about B)")
print(f"  qp_yy[(0,1)] = {qp01[2]:.12f}")
print(f"  qp_zz[(0,1)] = {qp01[5]:.12f}")

# Mulliken atomic dipoles:
# d_A[k] = -Σ_{μ on A} Σ_ν P(ν,μ) * dpint_k(μ,ν)  [about A]
# For atom A (basis 0): d_A[x] = -P[0,0]*dpint_x(0,0) - P[1,0]*dpint_x(0,1)
# dpint_x(0,0) = 0 (same center, odd integrand)
d_Ax = -P_off * dp01[0]
print(f"\nMulliken dipoles:")
print(f"  d_A[x] = {d_Ax:.12f}")

# For atom B (basis 1): d_B[x] = -P[0,1]*dpint_x(1,0) - P[1,1]*dpint_x(1,1)
# dpint_x(1,1) = 0 (same center)
d_Bx = -P_off * dpint_x_10
print(f"  d_B[x] = {d_Bx:.12f}")

# By symmetry: d_A = -d_B (check)
print(f"  d_A + d_B = {d_Ax + d_Bx:.12f} (should be ~0)")

# Mulliken atomic quadrupoles:
# q_A[c] = -Σ_{μ on A} Σ_ν P(ν,μ) * qpint_c(μ,ν)
# For atom A: q_A[xx] = -P[0,0]*qp_xx(0,0) - P[1,0]*qp_xx(0,1)
q_Axx = -P_diag * qp_diag[0] - P_off * qp01[0]
q_Ayy = -P_diag * qp_diag[2] - P_off * qp01[2]
q_Azz = -P_diag * qp_diag[5] - P_off * qp01[5]
print(f"\nMulliken quadrupoles:")
print(f"  q_A[xx] = {q_Axx:.12f}")
print(f"  q_A[yy] = {q_Ayy:.12f}")
print(f"  q_A[zz] = {q_Azz:.12f}")

q_Bxx = -P_diag * qp_diag[0] - P_off * qpint_xx_10
q_Byy = -P_diag * qp_diag[2] - P_off * qp01[2]  # same (no shift)
q_Bzz = -P_diag * qp_diag[5] - P_off * qp01[5]  # same
print(f"  q_B[xx] = {q_Bxx:.12f}")
print(f"  q_B[yy] = {q_Byy:.12f}")
print(f"  q_B[zz] = {q_Bzz:.12f}")

# AES kernel energy:
# E_kernel_d = Σ_A dkernel * |d_A|²
# E_kernel_q = Σ_A qkernel * Σ_c q_A[c]² * scale[c]
DKERNEL = 0.05563889
QKERNEL = 0.00027431
qp_scale = [1.0, 2.0, 1.0, 2.0, 2.0, 1.0]

e_dk = DKERNEL * (d_Ax**2 + d_Bx**2)  # only x-component nonzero
print(f"\nAES kernel energy:")
print(f"  E_dkernel = {e_dk:.12f} Ha")

# For quadrupole kernel, the qpat values include the off-diagonal xy,xz,yz = 0 for this case
q_A = [q_Axx, 0, q_Ayy, 0, 0, q_Azz]
q_B = [q_Bxx, 0, q_Byy, 0, 0, q_Bzz]
e_qk = 0.0
for c in range(6):
    e_qk += QKERNEL * (q_A[c]**2 * qp_scale[c] + q_B[c]**2 * qp_scale[c])
print(f"  E_qkernel = {e_qk:.12f} Ha")

e_aes_total = e_dk + e_qk
print(f"  E_AES total (kernel only) = {e_aes_total:.12f} Ha")
print(f"  Rust E_AES = 0.016520104100 Ha")
print(f"  diff = {e_aes_total - 0.016520104100:+.12f} Ha")

# Also check pairwise AES (should be ~0 for zero charges)
print(f"\n  Note: E_AES pairwise is ~0 because charges are ~0")

# Total energy check
e_band = -1.037985479133
e_scc = 0.0
e_3rd = 0.0
e_rep = 0.039534125669
e_disp = -0.000002315860
total = e_band + e_aes_total + e_scc + e_3rd + e_rep + e_disp
print(f"\nTotal energy (using Python AES):")
print(f"  E_band + E_AES + E_rep + E_disp = {total:.12f} Ha")
print(f"  Our Rust total = -0.981934126423 Ha") 
print(f"  tblite total   = -0.981983692538 Ha")
