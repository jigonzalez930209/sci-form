# Alchemical free energies: FEP, BAR and MBAR

The first v2.2 free-energy prototype implements the estimator layer for
alchemical transformations. A caller supplies samples with reduced potentials
at every alchemical state; the module returns dimensionless free-energy
differences and, when a temperature is set, kcal/mol values.

> **Status.** Implemented in `feat/wave-b`. This is the statistical estimator;
> the soft-core MD workflow and replicate exchange around it are still pending.
> See the [v2.x roadmap](/v2-roadmap).

## Why an estimator layer

Alchemical free energy is a workflow:

1. define the end states;
2. run MD at several `λ` windows;
3. collect reduced potentials;
4. estimate `ΔF` with FEP, BAR or MBAR;
5. diagnose overlap and uncertainty.

This module implements step 4 and the overlap diagnostics. It does **not**
create the λ windows or run the MD, so it can be tested against analytic
samples independently of the force field.

## Reduced potentials

For a sample `n` and state `k`,

$$
u_{kn} = \beta U_k(\mathbf{x}_n),
\qquad
\beta = \frac{1}{k_BT}.
$$

`u` and the resulting `Δf` are dimensionless (`kBT`). With a temperature `T`,
the result converts to

$$
\Delta G = \Delta f\,RT.
$$

## FEP

For samples drawn from state 0,

$$
\Delta f_{0\to1} =
-\ln \left\langle e^{-(u_1-u_0)} \right\rangle_0.
$$

FEP is exact in the infinite-sample limit but biased when the forward and
reverse work distributions do not overlap.

## BAR and MBAR

BAR is the two-state optimum; MBAR generalises it to `K` states and samples
from every state. The MBAR equations are

$$
f_k =
-\ln \sum_{n=1}^{N}
\frac{e^{-u_{kn}}}
{\sum_{j=1}^{K} N_j e^{f_j-u_{jn}}},
$$

solved self-consistently with `f_0 = 0`. The result reports:

- `delta_f_kbt` relative to state 0;
- `delta_g_kcal_mol` when the temperature is set;
- the symmetric overlap matrix `O_kl = Σ_n W_nk W_nl`;
- the effective sample size `1 / Σ_n W_nk²`;
- convergence and iteration count.

`Bar` and `Mbar` share the same multistate solver; `Bar` is the two-state case.

## API

```rust
use sci_form::free_energy::{
    AlchemicalSample, FreeEnergyConfig, FreeEnergyEstimator, estimate_free_energy,
};

let samples = vec![
    AlchemicalSample {
        origin_state: 0,
        reduced_potentials: vec![0.0, 1.5],
    },
    AlchemicalSample {
        origin_state: 1,
        reduced_potentials: vec![0.0, 1.5],
    },
];

let result = estimate_free_energy(
    &samples,
    FreeEnergyEstimator::Mbar,
    &FreeEnergyConfig::default(),
)?;

println!("Δf = {:.3} kBT", result.delta_f_kbt[1]);
println!("ΔG = {:.3} kcal/mol", result.delta_g_kcal_mol.unwrap()[1]);
println!("converged = {}", result.converged);
```

### Input rules

- Every sample has one reduced potential per state.
- Every state must contribute at least one sample.
- `origin_state` is inside `0..n_states`.
- FEP requires exactly two states.
- No non-finite reduced potentials.

## Validation

| Gate | Observable | Result |
|---|---|---|
| FEP constant offset | `Δf` from `u1-u0 = 1.5` | 1.5 `kBT` |
| MBAR constant offset | two states | 1.5 `kBT` |
| MBAR three states | offsets 0, 1, 2 | 1.0 and 2.0 `kBT` |
| Invalid shape | one state / wrong length | rejected |

These are analytic estimator tests. A soft-core MD benchmark and a hydration
free-energy comparison are still required.

## Limitations

- The module does not run MD or create λ windows.
- There is no soft-core potential, endpoint correction or Hamiltonian replica
  exchange yet.
- The overlap matrix is a diagnostic, not a proof that the estimate converged.
- Effective sample size is reported per state; it does not replace a block or
  correlation analysis.
- FEP requires state-0 samples and fails hard when none are supplied.
- The implementation assumes all samples are independent; correlated MD frames
  need subsampling or a correlation-aware estimator.
- No Python, WASM or CLI wrapper yet.

## Next steps

1. Add a soft-core λ scheduler on top of the existing dynamics.
2. Add Hamiltonian replica exchange and overlap-driven λ refinement.
3. Add block averaging and correlation-aware uncertainty.
4. Add endpoint corrections for charge changes.
5. Benchmark hydration free energies against experiment.
6. Add Python, WASM and CLI wrappers.
