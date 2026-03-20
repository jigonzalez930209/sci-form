# E6: Analytical Linearized Poisson-Boltzmann (ALPB)

**Module:** `sci_form::experimental::alpb`
**Feature flag:** `experimental-alpb`

---

## Overview

Implements the Klamt/Grimme **Analytical Linearized Poisson-Boltzmann (ALPB)** implicit solvation model. Computes solvation free energies $\Delta G_{solv}$ with analytical gradients at negligible computational cost, enabling thermodynamically realistic conformer ranking in solution.

<SvgDiagram src="/svg/experimental-alpb.svg" alt="ALPB Pipeline" />

---

## Theory

### Born Radii

Effective Born radii are computed using HCT (Hawkins-Cramer-Truhlar) integration enhanced with Onufriev-Bashford-Case (OBC) corrections:

$$\frac{1}{R_i^{eff}} = \frac{1}{\rho_i} - \frac{\tanh(\alpha\Psi_i - \beta\Psi_i^2 + \gamma\Psi_i^3)}{\rho_i}$$

with parameters $\alpha = 1.0$, $\beta = 0.8$, $\gamma = 4.85$.

### Generalized Born Kernel

$$f_{GB}(r_{ij}, R_i, R_j) = \sqrt{r_{ij}^2 + R_i R_j \exp\left(-\frac{r_{ij}^2}{4R_iR_j}\right)}$$

### ALPB Correction

The ALPB model improves standard GB with a universal correction:

$$E_{ALPB} = E_{GB} \cdot \frac{\varepsilon - 1}{\varepsilon + A \cdot x}$$

where $A = 0.571412$ (Klamt constant) and $x = E_{GB}/E_{Coulomb}$.

### Non-Polar Contribution

$$\Delta G_{np} = \gamma_{SA} \cdot \text{SASA}$$

where $\gamma_{SA}$ is the surface tension parameter (default: 0.005 kcal/mol/Å²).

---

## API

```rust
use sci_form::experimental::alpb::*;

// Compute Born radii
let born = compute_born_radii(&elements, &positions);
// born.radii: Vec<f64>

// Full ALPB solvation
let config = AlpbConfig {
    dielectric: 78.5,     // water
    probe_radius: 1.4,    // Å
    surface_tension: 0.005,
    ..Default::default()
};
let result = compute_alpb_solvation(&elements, &positions, &charges, &config);
// result.electrostatic_energy: f64  (kcal/mol)
// result.nonpolar_energy: f64       (kcal/mol)
// result.total_energy: f64          (kcal/mol)
// result.born_radii: Vec<f64>
```

---

## Solvent Presets

| Solvent | $\varepsilon$ |
|---------|--------------|
| Water | 78.5 |
| DMSO | 47.0 |
| Methanol | 32.7 |
| Acetone | 20.7 |
| Chloroform | 4.8 |
| Vacuum | 1.0 |

When $\varepsilon = 1.0$ (vacuum), the solvation energy is exactly zero — no division-by-zero issues.

---

## Comparison with Stable GB

| Feature | Stable GB-HCT | ALPB |
|---------|---------------|------|
| Born radii | HCT | HCT + OBC |
| Electrostatics | Standard GB | GB + ALPB correction |
| Universal constant | No | Yes ($A = 0.571412$) |
| Non-polar | SASA-based | SASA-based |
| Solvent-specific fitting | No | No |

---

## Tests

```bash
cargo test --features experimental-alpb --test regression -- test_alpb
```

8 integration tests covering: Born radii computation, GB kernel asymptotics, electrostatic solvation energy, vacuum limit (zero energy), solvent dielectric dependence, non-polar contributions, and total solvation energy.
