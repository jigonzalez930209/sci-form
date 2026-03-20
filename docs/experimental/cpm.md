# E10: Constant Potential Method (CPM)

**Module:** `sci_form::experimental::cpm`
**Feature flag:** `experimental-cpm`

---

## Overview

Extends charge and electronic property calculations by coupling the molecule to a "virtual electrode" at electrochemical potential $\mu$. Charges fluctuate to equilibrate with the electrode, enabling simulation of **quantum capacitance**, charge response to potential, and electrochemical energy surfaces.

<SvgDiagram src="/svg/experimental-cpm.svg" alt="CPM Pipeline" />

---

## Theory

### Grand Potential

The grand potential functional couples electronic energy with the electrode potential:

$$\Omega(\{q_i\}, \mu) = E_\text{elec}(\{q_i\}) - \mu \sum_i q_i$$

Minimizing with respect to charges $\{q_i\}$ yields the equilibrium charge distribution at potential $\mu$.

### Charge Equilibration at Fixed $\mu$

The SCF iteration solves:

$$q_i^{(n+1)} = \frac{\mu - \chi_i - \sum_{j \neq i} J_{ij} q_j^{(n)}}{\eta_i}$$

with damped update (50/50 mixing) and Coulomb interaction:

$$J_{ij} = \frac{14.3996}{\varepsilon \cdot r_{ij}} \quad \text{(kcal/mol per e²)}$$

### Electrochemical Energy Surface

Scanning $\mu$ produces the charge response curve $Q(\mu)$ and capacitance:

$$C(\mu) = \frac{\partial Q}{\partial \mu}$$

computed via finite differences over the $\mu$ scan.

---

## API

```rust
use sci_form::experimental::cpm::*;

// CPM charges at a given potential
let config = CpmConfig {
    mu_ev: -4.44,           // SHE reference (eV)
    dielectric: 78.5,       // water
    max_iter: 100,
    charge_tol: 1e-6,
};
let result = compute_cpm_charges(&elements, &positions, &config);
// result.charges: Vec<f64>
// result.total_charge: f64
// result.grand_potential: f64
// result.converged: bool
// result.iterations: usize

// Electrochemical surface scan
let surface = compute_cpm_surface(&elements, &positions, mu_min, mu_max, n_points, dielectric);
// surface.mu_values: Vec<f64>
// surface.total_charge: Vec<f64>
// surface.free_energy: Vec<f64>
// surface.capacitance: Vec<f64>
```

---

## Potential Scale

| Reference | $\mu$ (eV) |
|-----------|-----------|
| SHE (Standard Hydrogen Electrode) | $-4.44$ |
| +1 V vs SHE | $-3.44$ |
| -1 V vs SHE | $-5.44$ |
| Typical scan range | $[-5.5, -3.5]$ |

---

## Applications

- **Quantum capacitance** of molecular electrodes
- **Redox potential** estimation: $E^0 = -\Delta\Omega / nF$
- **Conformer ranking** under electrochemical bias
- **Charge response** analysis for molecular electronics

---

## Tests

```bash
cargo test --features experimental-cpm --test regression -- test_cpm
```

8 integration tests covering: charge equilibration convergence, total charge response to potential, monotonic charge-potential relationship, grand potential computation, capacitance positivity, electrochemical surface scanning, solvent dielectric effects, and multi-element systems.
