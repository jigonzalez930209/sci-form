# E9: Mobile Block Hessian (MBH)

**Module:** `sci_form::beta::mbh` (gated; F92 validated, promotion pending)
**Feature flag:** `beta-mbh` (`experimental-mbh` alias)

---

## Overview

Reduces vibrational analysis cost from $O(3N)$ energy evaluations to $O(n_\text{flex})$ by treating rigid groups (aromatic rings, fused ring systems) as rigid bodies with 6 degrees of freedom (3 translation + 3 rotation). Enables real-time IR spectra for large molecules where full Hessian computation is impractical.

<SvgDiagram src="/svg/experimental-mbh.svg" alt="MBH Pipeline" />

---

## Theory

### Rigid Block Detection

Rigid blocks are identified from the SSSR (Smallest Set of Smallest Rings):

- Aromatic rings → rigid
- Saturated rings ≤ 8 atoms → rigid
- Fused ring systems → merged into single blocks
- All other atoms → flexible (full 3 Cartesian DOF)

### DOF Reduction

Each rigid block contributes 6 collective DOF instead of $3N_\text{block}$:

$$n_\text{DOF} = 6 \cdot n_\text{blocks} + 3 \cdot n_\text{flex}$$

For a molecule with 60% of atoms in rings, the Hessian size is reduced to ~40%.

### Projection Matrix

The projection matrix $L$ maps reduced coordinates to Cartesian displacements:

- **Translation columns:** Uniform displacement of all block atoms along $x$, $y$, $z$
- **Rotation columns:** Infinitesimal rotation of block atoms about COM along $x$, $y$, $z$ axes
- **Flexible columns:** Standard Cartesian unit vectors for non-block atoms

### Generalized Eigenvalue Problem

$$H_\text{MBH} \mathbf{c} = \lambda \, M_\text{MBH} \mathbf{c}$$

where $H_\text{MBH} = L^T H L$ and $M_\text{MBH} = L^T M L$. Frequencies in cm$^{-1}$:

$$\nu_i = \frac{1}{2\pi c} \sqrt{\lambda_i} \times 108.59$$

---

## API

```rust
use sci_form::beta::mbh::*;

// Detect rigid blocks
let decomposition = detect_rigid_blocks(&elements, &positions, &ring_info);
// decomposition.blocks: Vec<RigidBlock>
// decomposition.flexible_atoms: Vec<usize>
// decomposition.n_dof_reduced: usize
// decomposition.n_dof_full: usize

// Build projection matrix
let l_matrix = build_projection_matrix(&decomposition, &positions);

// Compute MBH frequencies
let config = MbhConfig { fd_step: 0.005 };
let energy_fn = |coords: &[f64]| -> f64 { /* energy evaluation */ };
let result = compute_mbh_frequencies(&elements, &positions, &ring_info, &energy_fn, &config);
// result.frequencies: Vec<f64>    (cm⁻¹)
// result.n_dof_full: usize
// result.n_dof_reduced: usize
// result.n_blocks: usize
```

---

## Speedup

| Molecule | Full DOF | MBH DOF | Speedup |
|----------|----------|---------|---------|
| Benzene (6 ring + 6 H) | 36 | 24 | 1.5× |
| Naphthalene (10 ring + 8 H) | 54 | 30 | 1.8× |
| Aspirin (1 ring + rest) | 63 | ~45 | 1.4× |
| Penicillin (2 rings) | 99 | ~60 | 1.7× |

---

## Validation (F92, 2026-10-14)

External gate `tests/beta_mbh_validation.rs --features beta-mbh` vs
`scripts/validation/oracle_mbh.py` (diatomic stretches from SI constants,
validates the 108.59 factor): stretch 1.3e-5, benzene 1 block 24/36 DOF
speedup 1.5, zero modes ≤5 cm⁻¹ (linear-molecule quartic fd artifact).
Gap: experimental ±2 cm⁻¹ (Ghysels).

## Tests

```bash
cargo test --features "beta-cpm beta-mbh" --test experimental_gates -- test_mbh
```

8 integration tests (`tests/experimental/test_mbh.rs`) covering: rigid block detection, block COM computation, projection matrix orthogonality, DOF reduction verification, frequency computation, frequency ordering, and comparison between full and reduced Hessian results.

These tests live in the `experimental_gates` target — see [Overview](./overview.md#run-tests) for the full command matrix.
