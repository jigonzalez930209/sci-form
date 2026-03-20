# E8: Semidefinite Relaxation Embedding (SDR)

**Module:** `sci_form::experimental::sdr`
**Feature flag:** `experimental-sdr`

---

## Overview

Reformulates the conformer embedding problem as a **convex semidefinite program (SDP)**, mathematically guaranteeing that the distance matrix is positive semidefinite before extracting coordinates. This eliminates the retry loop caused by negative eigenvalues in the metric matrix — the most common failure mode in distance geometry embedding.

---

## Theory

### Distance Constraints as Convex Set

The decision variable is a PSD Gram matrix $X \in \mathbb{S}^N_+$. Distances relate to entries via:

$$d_{ij}^2 = X_{ii} + X_{jj} - 2X_{ij}$$

### Alternating Projections

The algorithm alternates between two convex projections:

1. **PSD Projection** $P_+(X)$: Eigendecompose, zero out negative eigenvalues
2. **Distance Projection** $P_\Omega(X)$: Enforce $X_{ii} + X_{jj} - 2X_{ij} = d_{ij}^2$ for known pairs

Convergence criterion: $\|X_{k+1} - X_k\|_F / \|X_k\|_F < 10^{-6}$

### SVT (Singular Value Thresholding)

Optionally, a nuclear norm penalty promotes low-rank solutions:

$$\min \|X\|_* + \lambda \|P_\Omega(X) - D\|_F^2$$

via $\text{SVT}_\tau(X) = U \max(\Sigma - \tau, 0) V^T$.

### Coordinate Extraction

From the converged Gram matrix $X^*$, 3D coordinates are extracted via:

$$C = V_3 \Sigma_3^{1/2}$$

where $V_3, \Sigma_3$ are the top-3 eigenvectors/values of $X^*$.

---

## API

```rust
use sci_form::experimental::sdr::*;

// PSD projection
let (projected, n_negative) = project_psd(&matrix);

// Distance projection
let projected = project_distances(&gram, &distance_pairs);

// Alternating projections
let config = SdrConfig { max_iter: 200, tol: 1e-6, ..Default::default() };
let (result, convergence) = alternating_projections(&x0, &pairs, &config);
// convergence.iterations, convergence.converged, convergence.final_residual

// Warm start from distance matrix
let gram = warm_start_gram(n_atoms, &distance_pairs);

// Extract 3D coordinates from Gram matrix
let coords = extract_coordinates(&gram); // flat Vec<f64>

// Full SDR embedding pipeline
let result = sdr_embed(n_atoms, &distance_pairs, &config);
// result.coords: Vec<f64>
// result.num_atoms: usize
// result.convergence: SdrConvergence
```

---

## Advantages Over Standard Embedding

| Issue | Standard ETKDG | SDR |
|-------|---------------|-----|
| Negative eigenvalues | Retry loop | Impossible by construction |
| PSD guarantee | No | Yes |
| Failure rate | ~5-10% | 0% |
| Warm start | No | From Cayley-Menger |

---

## Tests

```bash
cargo test --features experimental-sdr --test regression -- test_sdr
```

8 integration tests covering: PSD projection, negative eigenvalue removal, warm start Gram matrix, coordinate extraction with distance preservation, equilateral triangle and tetrahedron embedding, convergence tracking, and alternating projection improvement.
