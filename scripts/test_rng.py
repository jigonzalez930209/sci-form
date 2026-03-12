#!/usr/bin/env python3
"""Test the exact random number sequence from boost::minstd_rand + uniform_real.
We'll simulate the RNG and compare with RDKit's actual output."""

# boost::minstd_rand parameters
A = 48271
M = 2147483647  # 2^31 - 1

def minstd_rand_next(state):
    """Returns (next_state, output_value)"""
    state = (state * A) % M
    return state, state

def uniform_real_boost(val):
    """boost::uniform_real<>(0,1) mapping for minstd_rand.
    minstd_rand::min() = 1, minstd_rand::max() = M-1 = 2147483646
    
    The boost formula for integer engines:
    result = (val - min) / (max - min + 1) * (b - a) + a  
    For [0,1): = (val - 1) / (2147483646 - 1 + 1) = (val - 1) / 2147483646
    """
    return (val - 1) / (M - 1)  # (val - 1) / 2147483646

def uniform_real_our(val):
    """Our current Rust implementation: (val - 1) / (M - 2)"""
    return (val - 1) / (M - 2)  # (val - 1) / 2147483645

# Generate first 20 values with seed=42
state = 42
print(f"{'Step':>4} {'Raw':>12} {'Boost':>18} {'Ours':>18} {'Diff':>14}")
print("-" * 70)
for i in range(20):
    state, val = minstd_rand_next(state)
    boost_val = uniform_real_boost(val)
    our_val = uniform_real_our(val)
    diff = boost_val - our_val
    print(f"{i+1:4d} {val:12d} {boost_val:18.15f} {our_val:18.15f} {diff:14.2e}")

# The difference per value is small, but accumulated over many picks it can
# shift which torsion basin a conformer falls into.
print(f"\nMax possible diff per value: {1.0/(M-1) - 1.0/(M-2):.2e}")
print(f"This is about {abs(1.0/(M-1) - 1.0/(M-2)) / (1.0/(M-1)) * 100:.4f}% relative error per draw")
