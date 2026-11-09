# E2: Randomized Numerical Linear Algebra (RandNLA)

**Module:** `sci_form::rand_nla` (core since 2026-11-08; F93 validated)
**Feature flag:** none — `beta-randnla` and `experimental-randnla` are compatibility aliases; `sci_form::beta::rand_nla` re-exports the core module.

---

## Overview

Replaces the $O(N^3)$ EHT diagonalization with the **Randomized Nyström Approximation**, reducing cost to $O(N^2)$ or $O(N \log N)$ for sparse matrices. This enables EHT calculations on systems of thousands of atoms where full diagonalization is impractical.

<SvgDiagram src="/svg/experimental-randnla.svg" alt="RandNLA Pipeline" />

---

## Theory

### Nyström Approximation

Given a symmetric matrix $S \in \mathbb{R}^{N \times N}$:

1. **Sketch:** Generate $\Omega \in \mathbb{R}^{N \times k}$ with Gaussian entries.
2. **Project:** Compute $Y = S\Omega$ via sparse-dense multiplication
3. **Reconstruct:** $B = Q^T S Q$ and $S \approx Q B Q^T$, with $Y = QR$.

The rank parameter $k = \min(N, 50 + \lceil\log N\rceil)$ controls accuracy vs. speed.

### Randomized $S^{-1/2}$

Via randomized SVD: $S \approx U \Sigma V^T$, then $S^{-1/2} \approx U \Sigma^{-1/2} U^T$ with truncated singular values.

### Randomized Diagonalization

The transformed Hamiltonian $\tilde{H} = S^{-1/2} H S^{-1/2}$ is diagonalized
using a sketch + Rayleigh–Ritz projection. For `N > 500`,
`solve_eht_randnla_lowrank` replaces the full $N \times N$ diagonalization
with a reorthogonalized Lanczos iteration on a randomized start and returns
only the requested lowest states.

---

## API

The low-level interface works on the Hamiltonian and overlap matrices:

```rust
use nalgebra::DMatrix;
use sci_form::rand_nla::{
    solve_eht_randnla, solve_eht_randnla_lowrank, RandNlaConfig,
};

let h: DMatrix<f64> = /* EHT Hamiltonian, eV */;
let s: DMatrix<f64> = /* overlap matrix */;
let config = RandNlaConfig {
    sketch_size: Some(50), // k; None = 50 + ceil(log2 N)
    seed: 42,
    max_error: 1e-3,
    fallback_enabled: true,
    power_iterations: 2,
};

// Full randomized generalized eigensolver (all N eigenpairs).
let (energies, c, info) = solve_eht_randnla(&h, &s, &config);

// N > 500 path: only the lowest `n_eigen` states, no N×N diagonalization.
let low = solve_eht_randnla_lowrank(&h, &s, &config, 4);
// low.eigenvalues, low.eigenvectors, low.info
```

The convenience wrapper that assembles the EHT matrices is
`eht::solve_eht_randnla(elements, positions, k, config)`, which returns
`(EhtResult, RandNlaInfo)` — the usual EHT result plus rank/residual/fallback
metadata. `k` is the optional Wolfsberg–Helmholtz constant (`None` = default).

The Nyström building block is public as
`NystromApprox::from_matrix_with_power_iterations(&s, &sketch, q)` with
`GaussianSketch::new(&mut rng, n, k)`.

---

## Error Control

The module includes automatic error estimation:

$$\text{error} = \frac{\|HC - SCE\|_F}{\|SCE\|_F}$$

If the error exceeds 0.1%, the module falls back to exact diagonalization and reports in the `rand_nla_error` field.

---

## Validation (F93, 2026-10-14)

External gate `tests/beta_randnla_validation.rs` vs
`scripts/validation/oracle_randnla.py` (numpy/scipy, overlap spectra):
single-shot rank-N/4 ~0.13, near-identity regime 3.9e-3 with fallback OFF,
default config exact incl. exercised fallback at N=96.

**Subspace iterations (2026-11-07):** `NystromApprox::from_matrix_with_power_iterations`
and `RandNlaConfig::power_iterations` add re-orthonormalised power iterations.
With 384 iterations the same rank-N/4 sketches reach 3.8e-4 (F93a-bis), so
the <1e-3 accuracy criterion is met for the fixture.

**EHT backend (2026-11-07):** `eht::solve_eht_randnla` assembles the EHT
Hamiltonian/overlap exactly as `solve_eht` and routes the generalized
eigenproblem through the randomized solver, returning the usual `EhtResult`
plus the rank/residual/fallback metadata.

**N>500 path (2026-11-08):** `solve_eht_randnla_lowrank` uses Nyström for
`S^{-1/2}` and a reorthogonalized Lanczos iteration on a randomized start, so
it never forms or diagonalizes the full `N × N` transformed Hamiltonian.
Measured on `N = 600`, `k = 160`: the four lowest states have `9.6e-4`
relative error and the call takes ~80 ms versus ~143 ms for the exact
symmetric diagonalization. The track stays beta because the full-spectrum
path still needs the exact fallback.

## Tests

```bash
cargo test --test experimental_gates -- test_randnla
```

8 integration tests (`tests/experimental/test_randnla.rs`) covering: sketch matrix statistics, Nyström reconstruction error, randomized SVD accuracy, and EHT eigenvalue comparison against exact diagonalization for benzene and naphthalene.

These tests live in the `experimental_gates` target — see [Overview](./overview.md#run-tests) for the full command matrix.
