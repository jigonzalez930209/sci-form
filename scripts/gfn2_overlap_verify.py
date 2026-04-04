#!/usr/bin/env python3
"""Verify H2 overlap S[0,1] from our Rust code vs analytical STO-3G computation.
H slater exponent for GFN2 is 1.23 (NOT 1.823).
"""
import math

PI = math.pi

# STO-3G 1s coefficients (Stewart 1970, Table I)
STO3G_1S = [
    (2.227660584e+0, 1.543289673e-1),
    (4.057711562e-1, 5.353281423e-1),
    (1.098175104e-1, 4.446345422e-1),
]

ANGSTROM_TO_BOHR_SOLVER = 1.0 / 0.529177  # solver.rs constant
AATOAU = 1.8897259886  # gfn2_params.rs constant

def gaussian_norm_s(alpha):
    return (2.0 * alpha / PI) ** 0.75

def sto3g_overlap_ss(zeta_a, pos_a, zeta_b, pos_b):
    """STO-3G 1s-1s overlap between two centers (bohr)."""
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
            raw = (PI / gamma) ** 1.5
            s += ca * cb * na * nb * kab * raw
    return s

# GFN2 H parameters
zeta_H = 1.230  # slater exponent from gfn2_params.rs

# H2 at 0.74 Å
bond_ang = 0.74
bond_bohr_solver = bond_ang * ANGSTROM_TO_BOHR_SOLVER  # as used by our Rust CLI
bond_bohr_exact = bond_ang * AATOAU

print(f"H2 bond: {bond_ang} Å")
print(f"  solver.rs bohr: {bond_bohr_solver:.12f}")
print(f"  AATOAU   bohr: {bond_bohr_exact:.12f}")
print(f"  difference:     {bond_bohr_solver - bond_bohr_exact:.3e}")

s_solver = sto3g_overlap_ss(zeta_H, [0,0,0], zeta_H, [bond_bohr_solver, 0, 0])
s_exact = sto3g_overlap_ss(zeta_H, [0,0,0], zeta_H, [bond_bohr_exact, 0, 0])

print(f"\nOur Python overlap (solver.rs R): {s_solver:.15f}")
print(f"Our Python overlap (AATOAU R):    {s_exact:.15f}")
print(f"Rust debug S[0,1]:                 0.663779161430429")
print(f"Ratio (Rust/Python):              {0.663779161430429 / s_solver:.15f}")

# Check self-overlap (should be ~1.0 for STO-3G approximation)
s_self = sto3g_overlap_ss(zeta_H, [0,0,0], zeta_H, [0,0,0])
print(f"\nSelf-overlap (should be ~1.0):     {s_self:.15f}")

# Also do CH4 overlap for C 2s - H 1s
# C params: slater=[1.596, 1.596], pqn=[2,2], ngauss=[4,4]
STO4G_2S = [
    (1.161525551e+01,  -3.240841700e-02),
    (2.000243111e+00,   2.346814830e-01),
    (5.530163132e-01,   8.130631817e-01),
    (1.526398440e-01,   6.319787310e-02),
]

STO4G_2P = [
    (1.798260992e+00,  4.691987220e-02),
    (5.956243758e-01,  2.918000840e-01),
    (2.190438085e-01,  5.709033824e-01),
    (8.328515172e-02,  2.273750390e-01),
]

def sto4g_overlap_sp(zeta_s, pos_s, zeta_p, pos_p, l_p, m_p):
    """STO-4G 2p - 1s overlap. l_p=1, m_p=0 means px."""
    # For p orbital: lx,ly,lz depends on m
    # m=0: px (1,0,0), m=1: py (0,1,0), m=2: pz (0,0,1)
    lx_p, ly_p, lz_p = [(1,0,0), (0,1,0), (0,0,1)][m_p]
    
    r2 = sum((a - b)**2 for a, b in zip(pos_s, pos_p))
    s = 0.0
    for ar_s, c_s in STO3G_1S:
        alpha_s = ar_s * zeta_s**2
        n_s = gaussian_norm_s(alpha_s)
        for ar_p, c_p in STO4G_2P:
            alpha_p = ar_p * zeta_p**2
            # For p-type norm: N = (2α/π)^{3/4} * sqrt(4α)
            n_p = (2*alpha_p/PI)**0.75 * (4*alpha_p)**0.5
            gamma = alpha_s + alpha_p
            mu = alpha_s * alpha_p / gamma
            kab = math.exp(-mu * r2)
            
            # Gaussian product center
            px = (alpha_s * pos_s[0] + alpha_p * pos_p[0]) / gamma
            py = (alpha_s * pos_s[1] + alpha_p * pos_p[1]) / gamma
            pz = (alpha_s * pos_s[2] + alpha_p * pos_p[2]) / gamma
            
            # Obara-Saika for s-p overlap:
            # S(0,1) for the p-direction, S(0,0) for others
            inv2g = 0.5 / gamma
            s00 = (PI/gamma)**0.5
            
            dirs = [(pos_s[0], pos_p[0], px, lx_p),
                    (pos_s[1], pos_p[1], py, ly_p),
                    (pos_s[2], pos_p[2], pz, lz_p)]
            
            prod = 1.0
            for a_coord, b_coord, p_coord, l_b in dirs:
                pa = p_coord - a_coord
                pb = p_coord - b_coord
                if l_b == 0:
                    prod *= s00
                elif l_b == 1:
                    # s[0][1] = pb * s[0][0]
                    prod *= pb * s00
            
            s += c_s * c_p * n_s * n_p * kab * prod
    return s

zeta_C = 1.596  # C slater exponent in GFN2

# Test CH with C at origin, H along x at various distances
# In methane, C-H ~ 1.09 Å
r_CH_ang = 1.09
r_CH_bohr = r_CH_ang * ANGSTROM_TO_BOHR_SOLVER

print(f"\n\n=== C 2s - H 1s overlap at {r_CH_ang} Å ({r_CH_bohr:.6f} bohr) ===")
# For s-s overlap, we need the 2s-1s overlap
def sto_overlap_ns_ms(n_a, zeta_a, pos_a, ng_a, n_b, zeta_b, pos_b, ng_b):
    """Compute ns-ms STO overlap."""
    prims_a = STO3G_1S if ng_a == 3 else STO4G_2S
    prims_b = STO3G_1S if ng_b == 3 else STO4G_2S
    
    r2 = sum((a - b)**2 for a, b in zip(pos_a, pos_b))
    s = 0.0
    for ar, ca in prims_a:
        alpha_a = ar * zeta_a**2
        na = gaussian_norm_s(alpha_a)
        for br, cb in prims_b:
            alpha_b = br * zeta_b**2
            nb = gaussian_norm_s(alpha_b)
            gamma = alpha_a + alpha_b
            mu = alpha_a * alpha_b / gamma
            kab = math.exp(-mu * r2)
            raw = (PI / gamma) ** 1.5
            s += ca * cb * na * nb * kab * raw
    return s

# C 2s (STO-4G, zeta=1.596) - H 1s (STO-3G, zeta=1.23)
s_cs_h1s = sto_overlap_ns_ms(2, zeta_C, [0,0,0], 4, 1, zeta_H, [r_CH_bohr,0,0], 3)
print(f"C(2s)-H(1s) overlap: {s_cs_h1s:.15f}")

# C 2px - H 1s (bond along x)
s_cpx_h1s = sto4g_overlap_sp(zeta_H, [r_CH_bohr,0,0], zeta_C, [0,0,0], 0, 0)  # px
print(f"C(2px)-H(1s) overlap: {s_cpx_h1s:.15f}")

# From our Rust debug for methane:
# S[0]: [1.0, 0.0, 0.0, 0.0, 0.4449504428435358, 0.4449504428435358]
# S[4]: [0.4449504428435358, 0.3287809767660397, ...]
# Basis order: C(2s), C(2px), C(2py), C(2pz), H1(1s), H2(1s), H3(1s), H4(1s)
# Our S[0,4] = C(2s)-H1(1s) overlap
print(f"\nCompare with Rust methane S[0,4]: 0.4449504428435358")
print(f"Our Python C(2s)-H(1s):          {s_cs_h1s:.15f}")
print(f"Match: {abs(s_cs_h1s - 0.4449504428435358) < 1e-10}")
