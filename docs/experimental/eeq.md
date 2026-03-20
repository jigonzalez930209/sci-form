# E5: Dynamic EEQ Charge Model

**Module:** `sci_form::experimental::eeq`
**Feature flag:** `experimental-eeq`

---

## Overview

Implements **geometry-dependent Electronegativity Equalization (EEQ)** charges, capturing polarization effects that static charge models like Gasteiger-Marsili miss. Charges depend on the 3D molecular geometry through fractional coordination numbers, making them suitable for conformer-dependent property prediction and as input to dispersion and solvation models.

---

## Theory

### Electronegativity Equalization

At equilibrium, all atoms share a common electrochemical potential $\mu$:

$$\mu = \chi_i + \eta_i q_i + \sum_{j \neq i} \gamma(r_{ij}) q_j$$

This yields an $(N+1) \times (N+1)$ linear system (including charge neutrality constraint $\sum q_i = Q_{total}$).

### Fractional Coordination Number

Geometry dependence enters through the coordination number:

$$CN_i = \sum_{j \neq i} \frac{1}{1 + \exp\left(-16\left(\frac{r_{cov,ij}}{r_{ij}} - 1\right)\right)}$$

The Fermi-function damping provides smooth, differentiable CN values.

### Damped Coulomb Interaction

The pairwise interaction kernel uses error-function damping:

$$\gamma(r_{ij}) = \frac{\text{erf}\left(\sqrt{2} \cdot r_{ij} / \sigma_{ij}\right)}{r_{ij}}$$

where $\sigma_{ij}$ is derived from atomic charge radii.

---

## API

```rust
use sci_form::experimental::eeq::*;

// Compute EEQ charges
let result = compute_eeq_charges(&elements, &positions, total_charge);
// result.charges: Vec<f64>
// result.coordination_numbers: Vec<f64>

// Fractional coordination numbers
let cn = fractional_coordination(&elements, &positions);

// EEQ electrostatic energy
let energy_result = compute_eeq_energy(&elements, &positions, total_charge);
// energy_result.energy_kcal_mol: f64
// energy_result.charges: Vec<f64>

// Numerical gradient
let gradient = compute_eeq_gradient(&elements, &positions, total_charge, 0.001);
```

---

## Supported Elements

H, B, C, N, O, F, Si, P, S, Cl, Br, I, and transition metals supported by xTB (Ti, Cr, Mn, Fe, Co, Ni, Cu, Zn).

---

## Comparison with Gasteiger

| Property | Gasteiger | EEQ |
|----------|-----------|-----|
| Geometry-dependent | No | Yes |
| Coordination effects | No | Yes (CN) |
| Charge neutrality | Iterative | Exact (constraint) |
| Transition metals | No | Yes |
| Gradients | No | Yes (numerical) |

---

## Tests

```bash
cargo test --features experimental-eeq --test regression -- test_eeq
```

8 integration tests covering: charge neutrality, coordination numbers, energy computation, gradient accuracy, element coverage, and comparison with Gasteiger charges.
