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

The CPM equations are linear:

$$q_i = \frac{\mu_i - \chi_i - \sum_{j \neq i} J_{ij} q_j}{\eta_i}
\iff \sum_j M_{ij} q_j = \mu_i - \chi_i,
\quad M_{ii} = \eta_i,\; M_{ij} = J_{ij}$$

with Coulomb interaction:

$$J_{ij} = \frac{14.3996}{\varepsilon \cdot r_{ij}} \quad \text{(kcal/mol per e²)}$$

The solver inverts $M$ directly (LU). The damped fixed-point iteration
($q \leftarrow (q + q_\text{new})/2$) was removed in 2026-10-13: it diverges
for electrode-sized conducting systems (charges blow up to 1e157 on a
graphene bilayer — spectral radius of $J/\eta$ exceeds 1). LAMMPS
`electrode/conp` likewise inverts the elastance matrix (`mat_inv`).

### Multi-Electrode Setups

`CpmConfig::mu_per_atom` (optional, `None` = uniform `mu_ev`) assigns one
applied potential per atom — e.g. plate A at $+V/2$, plate B at $-V/2$, as in
LAMMPS `fix electrode/conp ... couple`. The grand potential generalizes to
$\Omega = E_\text{elec} - \sum_i \mu_i q_i$.

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
    mu_per_atom: None,      // Some(vec) for multi-electrode setups
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

## Validation vs LAMMPS (F95, 2026-10-13)

Oracle `scripts/validation/oracle_cpm_lammps.py` + fixture
`tests/fixtures/validation/cpm_lammps.json` + gate
`test_cpm_bilayer_matches_lammps`: graphene bilayer capacitor (two 3×3-cell
flakes, AB @ 3.35 Å) at V = ±0.5/1/2 V vs LAMMPS `fix electrode/conp`
($\eta = 1.979$). Compared on voltage-induced charges
($\Delta Q = Q(V) - Q(0)$, where $\chi$ cancels exactly): both exactly linear
and antisymmetric; plate-Q scale 1.76x, per-atom pattern cosine similarity
0.44 (tolerances 2.5x / 0.3).

Model gap (documented): our point-charge hardness ($\eta_C = 5.0$ eV) vs the
Gaussian self-energy (~22 eV) — ours over-polarizes dense conductors
(onsite too soft, edge oscillations). Path forward: Gaussian smearing in $J$
(ELECTRODE-style) and/or retuned onsite hardness.

LAMMPS recipe notes: `atom_modify map yes`, `kspace_modify gewald` for
uncharged starts, `pair_coeff` after `read_data`, two plates with `couple`
(a lone electrode group stays neutral), charges solve while stepping
(`run 1`), dump per-atom `q`.
