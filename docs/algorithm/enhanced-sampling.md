# Enhanced sampling: well-tempered metadynamics

The first v2.x enhanced-sampling implementation adds a collective-variable
bias to a potential. It is designed to sit on top of the existing UFF force
field, but the engine is generic: any type implementing the `Potential` trait
can be sampled. The initial release covers **well-tempered metadynamics with a
Langevin thermostat**; OPES and replica exchange are follow-up implementations.

> **Status.** Implemented in `feat/wave-a`. The API is available in Rust; the
> Python, WASM and CLI wrappers are still pending. See the
> [v2.x roadmap](/v2-roadmap).

## Why metadynamics

A plain MD trajectory visits the basins it can reach on the simulation time
scale. A free-energy difference between two basins, a conformational
transition or a desorption event may need much longer. Metadynamics adds a
history-dependent bias along a collective variable (CV) so the system is
pushed out of visited regions and the bias converges to the negative free
energy.

## Collective variables

| CV | Definition | Unit |
|---|---|---|
| `Coordinate { index }` | One Cartesian component | Å |
| `Distance { i, j }` | Distance between atoms `i` and `j` | Å |
| `Angle { i, j, k }` | Angle `i–j–k` | rad |
| `Dihedral { i, j, k, l }` | Dihedral `i–j–k–l` | rad |

Distance uses an analytic Cartesian gradient. Angle and dihedral use a central
difference with a 1e-5 Å step; this is robust at the linear limits where the
analytic expression is easy to get wrong.

## Well-tempered bias

Gaussians are deposited along the CV:

$$
V_b(\mathbf{s},t) =
\sum_k h_k \exp\left(
-\frac{(\mathbf{s}-\mathbf{s}_k)^2}{2\sigma^2}
\right).
$$

The well-tempered height decays with the local bias:

$$
h_k = h_0 \exp\left(
-\frac{V_b(\mathbf{s}_k,t)}{R\,\Delta T}
\right),
\qquad
\Delta T = (\gamma-1)T.
$$

`γ > 1` is the bias factor, `h0` the initial height and `σ` the Gaussian width.
A larger `γ` keeps the bias lower and samples longer; a smaller `σ` gives a
finer resolution but needs more deposits.

The bias force added to the physical force is

$$
\mathbf{F}_b =
-\frac{\partial V_b}{\partial \mathbf{s}}
\frac{\partial \mathbf{s}}{\partial \mathbf{x}}.
$$

## Integrator

The implementation uses velocity Verlet with an Ornstein–Uhlenbeck Langevin
thermostat:

$$
\mathbf{v}_{n+1} =
e^{-\gamma\Delta t}\mathbf{v}_n +
\sqrt{\frac{k_BT}{m}\left(1-e^{-2\gamma\Delta t}\right)}\,\boldsymbol{\xi}_n.
$$

The friction is configured in ps⁻¹ and converted to fs⁻¹ internally. Masses
are in amu; the same conversion factor as the existing MD module is used.

## Free-energy surface

At the end of the run, the total bias is evaluated on a one-dimensional grid:

$$
F(s) \approx -V_b(s),
$$

shifted so that `min(F) = 0`. The result also contains the CV trace, every
deposited Gaussian and the final coordinates. The surface is only meaningful
after the bias has converged; a short run is flagged in `notes`.

## API

### Generic potential

```rust
use sci_form::enhanced_sampling::{
    CollectiveVariable, MetaDynamicsConfig, Potential, run_metadynamics,
};

struct Harmonic { k: f64 }
impl Potential for Harmonic {
    fn energy(&self, x: &[f64]) -> f64 { 0.5 * self.k * x[0] * x[0] }
    fn gradient(&self, x: &[f64], g: &mut [f64]) { g[0] = self.k * x[0]; }
}

let config = MetaDynamicsConfig {
    cv: CollectiveVariable::Coordinate { index: 0 },
    max_steps: 5_000,
    ..MetaDynamicsConfig::default()
};

let result = run_metadynamics(&Harmonic { k: 1.0 }, &[1.0], &[1.0], &config)?;
```

`masses_amu` has one value per Cartesian coordinate, matching the flat
coordinate array. The molecule helper repeats each atomic mass three times.

### UFF molecule

```rust
use sci_form::enhanced_sampling::{CollectiveVariable, MetaDynamicsConfig};
use sci_form::run_uff_metadynamics;

let conformer = sci_form::embed("CCO", 42);
let config = MetaDynamicsConfig {
    cv: CollectiveVariable::Dihedral { i: 0, j: 1, k: 2, l: 3 },
    max_steps: 20_000,
    sigma: 0.2,
    grid_min: -3.14,
    grid_max: 3.14,
    grid_bins: 121,
    ..MetaDynamicsConfig::default()
};

let result = run_uff_metadynamics("CCO", &conformer.coords, &config)?;
```

The returned `MetaDynamicsResult` has:

| Field | Meaning |
|---|---|
| `cv_trace` | Step, CV, unbiased potential and bias energy |
| `gaussians` | Step, centre, height and width of every deposit |
| `free_energy` | Grid bounds, bins and `F(s)` in kcal/mol |
| `final_coords` | Final flat coordinates in Å |
| `final_cv` | Final CV value |
| `notes` | Convergence and scope caveats |

## Validation

| Gate | Observable | Result |
|---|---|---|
| Distance CV gradient | Analytic vs finite difference | < 1e-12 |
| 1D double well | Visits both minima | `min < -0.5`, `max > 0.5` |
| UFF gradient | Analytic vs finite difference | < 1e-4 kcal/mol/Å |
| UFF metadynamics | 200-step ethanol run | finite trace and deposits |
| Invalid bias factor | `γ = 1` | rejected |

The double-well test checks exploration, not a converged barrier. A quantitative
barrier benchmark needs a longer run and is part of the promotion gate.

## Limitations

- Only one CV is biased per run. Multi-dimensional metadynamics is a follow-up.
- OPES and replica exchange are not implemented yet.
- The free-energy surface is the negative bias, not a reweighted histogram.
  Reweighting is needed for a quantitative comparison with unbiased MD.
- The Langevin thermostat is the only thermostat in this module; the existing
  NVT/NPT integrators are not yet coupled to the bias.
- The CV gradient is numerical for angles and dihedrals.
- A bad CV hides the slow mode. The API cannot detect that for the user.

## Next steps

1. Add multi-dimensional CVs and a grid in more than one dimension.
2. Add OPES and Hamiltonian replica exchange.
3. Add reweighting from the biased trajectory.
4. Wire the bias into the existing NVT/NPT dynamics module.
5. Add Python, WASM and CLI wrappers.
6. Benchmark a known free-energy landscape against long unbiased MD.
