# E2: Randomized Numerical Linear Algebra (RandNLA)

**Module:** `sci_form::experimental::rand_nla`
**Feature flag:** `experimental-randnla`

---

## Overview

Replaces the $O(N^3)$ EHT diagonalization with the **Randomized Nyström Approximation**, reducing cost to $O(N^2)$ or $O(N \log N)$ for sparse matrices. This enables EHT calculations on systems of thousands of atoms where full diagonalization is impractical.

---

## Theory

### Nyström Approximation

Given a symmetric matrix $S \in \mathbb{R}^{N \times N}$:

1. **Sketch:** Generate $\Omega \in \mathbb{R}^{N \times k}$ with Gaussian entries $\mathcal{N}(0, 1/k)$
2. **Project:** Compute $Y = S\Omega$ via sparse-dense multiplication
3. **Reconstruct:** $S \approx Y(\Omega^T Y)^{-1} Y^T$

The rank parameter $k = \min(N, 50 + \lceil\log N\rceil)$ controls accuracy vs. speed.

### Randomized $S^{-1/2}$

Via randomized SVD: $S \approx U \Sigma V^T$, then $S^{-1/2} \approx U \Sigma^{-1/2} U^T$ with truncated singular values.

### Randomized Diagonalization

The transformed Hamiltonian $\tilde{H} = S^{-1/2} H S^{-1/2}$ is diagonalized using sketch + Rayleigh-Ritz refinement, extracting only the $n_{occ}/2$ occupied eigenvectors.

---

## API

```rust
use sci_form::experimental::rand_nla::*;

// Nyström approximation
let config = NystromConfig { rank: 50, seed: 42 };
let approx = nystrom_approximate(&matrix, &config);

// Randomized SVD
let rsvd = randomized_svd(&matrix, rank, seed);

// Randomized S^{-1/2}
let s_inv_sqrt = randomized_inverse_sqrt(&overlap, rank, seed);

// Full randomized EHT pipeline
let result = rand_eht_calculate(&elements, &coords, &config);
// result.homo_energy, result.gap, result.rand_nla_error
```

---

## Error Control

The module includes automatic error estimation:

$$\text{error} = \frac{\|HC - SCE\|_F}{\|SCE\|_F}$$

If the error exceeds 0.1%, the module falls back to exact diagonalization and reports in the `rand_nla_error` field.

---

## Tests

```bash
cargo test --features experimental-randnla --test regression -- test_randnla
```

Covers: sketch matrix statistics, Nyström reconstruction error, randomized SVD accuracy, and EHT eigenvalue comparison against exact diagonalization for benzene and naphthalene.
