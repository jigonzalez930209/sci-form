# E3: Riemannian Optimization for ETKDG

**Module:** `sci_form::riemannian` (core since 2026-10-14; `beta::riemannian` re-export)
**Feature flag:** none (core); `beta-riemannian` kept as no-op alias

---

## Overview

Replaces the Euclidean BFGS optimizer in the embedding step with **Riemannian L-BFGS** over the manifold of fixed-rank positive semidefinite (PSD) matrices. Negative eigenvalues are eliminated by design, reducing the retry loop failure rate to zero.

<SvgDiagram src="/svg/experimental-riemannian.svg" alt="Riemannian Optimization Pipeline" />

---

## Theory

### PSD Manifold

The metric matrix $G$ must be PSD for valid 3D coordinates. Instead of optimizing in $\mathbb{R}^{N \times N}$ and projecting, we optimize directly on:

$$\mathcal{M} = \{X \in \mathbb{R}^{N \times N} : X \succeq 0, \text{rank}(X) = r\}$$

### Retraction

The retraction maps a tangent vector back to the manifold:

$$\text{Retr}_X(\xi) = P_+\left(X + \xi + \frac{1}{2}\xi X^{-1} \xi\right)$$

where $P_+$ projects onto the PSD cone by zeroing negative eigenvalues.

### Riemannian Gradient

The Euclidean gradient is projected to the tangent space:

$$\text{grad}_R f = P_{T_X} \nabla f$$

For the Stiefel manifold: $P_{T_X}(\xi) = \xi - X \cdot \text{sym}(X^T \xi)$

### L-BFGS on the Manifold

- Store $m = 5$ curvature pairs $(s_k, y_k)$ in the tangent space
- Two-loop recursion with vector transport along geodesics
- Line search with Armijo condition along the manifold geodesic

---

## API

```rust
use sci_form::riemannian::*;

// PSD projection
let (projected, n_neg) = psd_projection(&matrix);

// Retraction (first-order): P₊(X + ξ)
let retracted = psd_retraction(&x, &xi);

// Full L-BFGS optimization under distance bounds
let config = RiemannianConfig {
    max_iter: 500,
    grad_tol: 1e-6,
    ..Default::default()
};
let opt = RiemannianLbfgs::new(config);
let result = opt.optimize(&x0, &constraints);
// result.metric_matrix, result.objective, result.converged
```

---

## Convergence

- **Primary:** $\|\text{grad}_R f\| < 10^{-6}$
- **Secondary:** Maximum distance violation $< 0.01$ Å
- **Fallback:** Reverts to Euclidean BFGS + PSD projection if not converged in 500 iterations

---

## Validation (F94, 2026-10-14)

External gate `tests/beta_riemannian_validation.rs` vs
`scripts/validation/oracle_riemannian.py` (numpy + scipy L-BFGS-B):
projection 8.2e-12, Frobenius distance 2.2e-16, objective/gradient 1.2e-16,
optimizer reaches the feasible set like scipy. Documented gap:
`psd_distance` is the Frobenius proxy, NOT the affine-invariant geodesic
distance (geomstats cross-check unmet). Note: the implemented retraction is
first-order P₊(X+ξ), not the second-order form once sketched in Theory.

## Tests

```bash
cargo test --features experimental-riemannian --test regression -- test_riemannian
```

Covers: PSD projection correctness, retraction stays on manifold, Riemannian gradient orthogonality, L-BFGS convergence on Rosenbrock and distance geometry objectives.
