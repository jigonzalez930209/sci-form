# Enhanced sampling: well-tempered metadynamics

The first v2.x enhanced-sampling implementation adds a collective-variable
bias to a potential. It is designed to sit on top of the existing UFF force
field, but the engine is generic: any type implementing the `Potential` trait
can be sampled. The module covers **well-tempered metadynamics**, an
**OPES-style kernel-density bias**, **histogram reweighting** and **temperature
replica exchange**, all with a Langevin thermostat.

> **Status.** Implemented in `feat/wave-a`: WTMetaD, importance-weighted OPES,
> a self-consistent 1D flow, 2D multi-CV OPES/reweighting, block reweighting
> and temperature replica exchange are in Rust. The exact OPES differential
> flow and the Python/WASM/CLI wrappers are still pending. See the
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
| `reweighted_free_energy` | Histogram-reweighted `F(s)` in kcal/mol |
| `final_coords` | Final flat coordinates in Å |
| `final_cv` | Final CV value |
| `notes` | Convergence and scope caveats |

## OPES-style bias

`run_opes` replaces the Gaussian history with an importance-weighted
kernel-density estimate of the CV values. Each deposited kernel carries
`w_i = exp(βV(s_i))`, so the density is the reweighted estimate
`P(s) = Σ w_i K(s,s_i) / Σ w_i`:

$$
V(s) = \left(1-\frac{1}{\gamma}\right)k_BT \ln\left(\frac{P(s)}{P_{\max}}+\epsilon\right),
$$

where `P` is the weighted Gaussian kernel density and `P_max` is its maximum
on the deposited kernels. The bias is zero at the most visited point and
negative where the density is low, which lowers the effective barrier. This is
an importance-weighted OPES approximation on a bounded grid; it does not
integrate the exact OPES differential flow equation.

```rust
use sci_form::enhanced_sampling::{OpesConfig, run_opes};

let config = OpesConfig {
    cv: CollectiveVariable::Coordinate { index: 0 },
    max_steps: 5_000,
    ..OpesConfig::default()
};
let result = run_opes(&potential, &initial_coords, &masses, &config)?;
```

`OpesResult` contains the CV trace, the deposited kernels, the bias-derived
surface and the reweighted surface.

## Self-consistent OPES flow

The kernel weights are not a single-pass estimate. At every deposit, the
engine iterates

$$
w_i \leftarrow e^{\beta V(s_i)},
\qquad
P(s) = \frac{\sum_i w_i K(s,s_i)}{\sum_i w_i},
$$

for `flow_iterations` fixed-point steps (`5` by default). This is the
importance-weight flow on the kernel set; it is not the exact time-dependent
OPES differential equation.

## Multi-CV OPES

`run_multi_opes` implements a two-dimensional OPES surface:

```rust
use sci_form::enhanced_sampling_multi::MultiOpesConfig;
use sci_form::run_multi_opes;

let config = MultiOpesConfig {
    max_steps: 10_000,
    pace_steps: 100,
    x_min: -2.0,
    x_max: 2.0,
    y_min: -2.0,
    y_max: 2.0,
    ..MultiOpesConfig::default()
};
let result = run_multi_opes(&potential, &initial_coords, &masses, &config)?;
```

`MultiOpesResult` contains:

| Field | Meaning |
|---|---|
| `trace` | Step, two CV values, potential and bias |
| `kernels` | 2D centres and importance weights |
| `free_energy` | 2D bias-derived surface |
| `reweighted_free_energy` | 2D histogram-reweighted surface |

`reweight_free_energy_2d` builds the reweighted surface directly from a trace.

## Reweighting

`reweight_free_energy` converts a biased CV trace into an unbiased histogram.
Each sample is weighted by `exp(+V_bias/kBT)`:

$$
F(s) = -k_BT\ln\left(\sum_{i\in s} e^{V_b(s_i)/k_BT}\right).
$$

The result is shifted so the minimum is zero. For correlated samples, use
`reweight_free_energy_correlated` with a `block_size`; it returns the same
surface plus a per-bin block standard error.

```rust
use sci_form::enhanced_sampling::{ReweightConfig, reweight_free_energy};

let surface = reweight_free_energy(&result.cv_trace, &ReweightConfig {
    grid_min: -1.5,
    grid_max: 1.5,
    grid_bins: 61,
    temperature_k: 300.0,
    block_size: 1,
})?;
```

## Replica exchange

The same module also provides temperature replica exchange. Each replica runs
at its own temperature, and adjacent pairs are exchanged every `swap_interval`
steps with

$$
P_{\mathrm{acc}} =
\min\left[1,\,
e^{(\beta_i-\beta_j)(U_i-U_j)}
\right],
\qquad
\beta = \frac{1}{k_BT}.
$$

The result contains one potential-energy trace per replica, every swap proposal
with its Metropolis exponent, and the final coordinates. It is a generic
`Potential` engine, so it works with the analytic double well and with the UFF
adapter.

```rust
use sci_form::enhanced_sampling::{ReplicaExchangeConfig, run_replica_exchange};

let config = ReplicaExchangeConfig {
    temperatures_k: vec![300.0, 400.0, 500.0, 600.0],
    max_steps: 10_000,
    ..ReplicaExchangeConfig::default()
};
let result = run_replica_exchange(&potential, &initial_coords, &masses, &config)?;
```

Replica exchange is not a substitute for a converged metadynamics bias; it
improves barrier crossing when the replicas overlap in energy. The swap
acceptance rate should be inspected in `result.swaps`.

## Validation

| Gate | Observable | Result |
|---|---|---|
| Distance CV gradient | Analytic vs finite difference | < 1e-12 |
| 1D double well | Visits both minima | `min < -0.5`, `max > 0.5` |
| UFF gradient | Analytic vs finite difference | < 1e-4 kcal/mol/Å |
| UFF metadynamics | 200-step ethanol run | finite trace and deposits |
| Invalid bias factor | `γ = 1` | rejected |
| Replica exchange | 4 replicas, 500 steps | swap proposals recorded, all final coordinates finite |
| OPES double well | 4,000 steps | both minima visited |
| 2D OPES harmonic well | 2,000 steps | finite trace and surfaces |
| 2D reweighting constant bias | flat trace | finite bins zero |
| Reweighting constant bias | flat trace | zero surface |
| Block reweighting | 3 samples, block 2 | finite non-negative standard errors |

The double-well test checks exploration, not a converged barrier. A quantitative
barrier benchmark needs a longer run and is part of the promotion gate.

## Limitations

- Only one CV is biased per run. Multi-dimensional metadynamics is a follow-up.
- The OPES implementation is an importance-weighted bounded-grid
  approximation, not the exact OPES differential flow equation.
- Multi-CV support is two-dimensional; higher-dimensional grids are not
  implemented.
- The block error is a simple block-average diagnostic, not a full
  autocorrelation analysis.
- The bias-derived surface and the reweighted surface can disagree on short runs.
- The Langevin thermostat is the only thermostat in this module; the existing
  NVT/NPT integrators are not yet coupled to the bias.
- The CV gradient is numerical for angles and dihedrals.
- A bad CV hides the slow mode. The API cannot detect that for the user.

## Next steps

1. Add multi-dimensional CVs and a grid in more than one dimension.
2. Implement the exact OPES differential flow equations.
3. Extend multi-CV grids above two dimensions.
4. Add autocorrelation-aware error propagation.
4. Wire the bias into the existing NVT/NPT dynamics module.
5. Add Python, WASM and CLI wrappers.
6. Benchmark a known free-energy landscape against long unbiased MD.
