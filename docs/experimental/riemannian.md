# E3: Riemannian Optimization for ETKDG

**Module:** `sci_form::experimental::riemannian`
**Feature flag:** `experimental-riemannian`

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
use sci_form::experimental::riemannian::*;

// PSD projection
let (projected, n_neg) = project_psd(&matrix);

// Retraction
let retracted = retract_psd(&x, &tangent_vector, epsilon);

// Riemannian gradient
let rgrad = riemannian_gradient(&x, &euclidean_grad);

// Full L-BFGS optimization
let config = RiemannianConfig {
    max_iter: 500,
    tol: 1e-6,
    memory: 5,
    ..Default::default()
};
let result = riemannian_lbfgs(&x0, &objective, &gradient, &config);
// result.x, result.iterations, result.converged, result.final_gradient_norm
```

---

## Convergence

- **Primary:** $\|\text{grad}_R f\| < 10^{-6}$
- **Secondary:** Maximum distance violation $< 0.01$ Å
- **Fallback:** Reverts to Euclidean BFGS + PSD projection if not converged in 500 iterations

---

## Tests

```bash
cargo test --features experimental-riemannian --test regression -- test_riemannian
```

Covers: PSD projection correctness, retraction stays on manifold, Riemannian gradient orthogonality, L-BFGS convergence on Rosenbrock and distance geometry objectives.
