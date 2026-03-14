#!/usr/bin/env python3
"""Simulate the exact RNG and distance picking for ethane to generate the distance matrix.
Compare with RDKit's actual embedding by checking if initial coords match."""

# RDKit bounds for ethane (8 atoms, CC + 6H)
# From the comparison: bounds are identical between RDKit and our code
# So let's compute the random distances
A = 48271
M = 2147483647

def minstd_next(state):
    state = (state * A) % M
    return state

def next_double(state):
    # Both candidate divisors
    state = minstd_next(state)
    val_ours = (state - 1) / (M - 2)  # our code: M-2 = 2147483645
    val_boost = (state - 1) / (M - 1)  # boost: M-1 = 2147483646
    return state, val_ours, val_boost

# Ethane bounds (from RDKit, confirmed identical to ours)
# upper = bm[j,i] for j<i, lower = bm[i,j] for i>j  
bounds_upper = {}
bounds_lower = {}
def set_b(j, i, ub, lb):
    bounds_upper[(j,i)] = ub
    bounds_lower[(i,j)] = lb

set_b(0,1, 1.524000000000, 1.504000000000)
set_b(0,2, 1.119400794878, 1.099400794878)
set_b(1,2, 2.195066594155, 2.115066594155)
set_b(0,3, 1.119400794878, 1.099400794878)
set_b(1,3, 2.195066594155, 2.115066594155)
set_b(2,3, 1.851965580853, 1.771965580853)
set_b(0,4, 1.119400794878, 1.099400794878)
set_b(1,4, 2.195066594155, 2.115066594155)
set_b(2,4, 1.851965580853, 1.771965580853)
set_b(3,4, 1.851965580853, 1.771965580853)
set_b(0,5, 2.195066594155, 2.115066594155)
set_b(1,5, 1.119400794878, 1.099400794878)
set_b(2,5, 3.075381000435, 2.254651189939)
set_b(3,5, 3.075381000435, 2.254651189939)
set_b(4,5, 3.075381000435, 2.254651189939)
set_b(0,6, 2.195066594155, 2.115066594155)
set_b(1,6, 1.119400794878, 1.099400794878)
set_b(2,6, 3.075381000435, 2.254651189939)
set_b(3,6, 3.075381000435, 2.254651189939)
set_b(4,6, 3.075381000435, 2.254651189939)
set_b(5,6, 1.851965580853, 1.771965580853)
set_b(0,7, 2.195066594155, 2.115066594155)
set_b(1,7, 1.119400794878, 1.099400794878)
set_b(2,7, 3.075381000435, 2.254651189939)
set_b(3,7, 3.075381000435, 2.254651189939)
set_b(4,7, 3.075381000435, 2.254651189939)
set_b(5,7, 1.851965580853, 1.771965580853)
set_b(6,7, 1.851965580853, 1.771965580853)

# Distance picking: for i=1..n, for j=0..i
state = 42
print(f"{'i':>3} {'j':>3} {'rval_ours':>18} {'rval_boost':>18} {'d_ours':>14} {'d_boost':>14} {'diff':>12}")
for i in range(1, 8):
    for j in range(i):
        state, rval_ours, rval_boost = next_double(state)
        ub = bounds_upper[(j, i)]
        lb = bounds_lower[(i, j)]
        d_ours = lb + rval_ours * (ub - lb)
        d_boost = lb + rval_boost * (ub - lb)
        diff = abs(d_ours - d_boost)
        print(f"{i:3d} {j:3d} {rval_ours:18.15f} {rval_boost:18.15f} {d_ours:14.12f} {d_boost:14.12f} {diff:12.2e}")
