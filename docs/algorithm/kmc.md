# Kinetic Monte Carlo for reaction networks

The second half of the v2.x reactivity workflow is a stochastic solver for an
elementary-step network. The existing `kinetics` module already evaluates HTST
rates and integrates a deterministic microkinetic ODE; the new module propagates
integer species counts with the Gillespie direct method.

> **Status.** Implemented in `feat/wave-a` for an **explicit** network. Automatic
> network discovery from SMIRKS/graph rules is still pending; see the
> [v2.x roadmap](/v2-roadmap).

## Why kMC

A deterministic ODE is the right tool for large, well-mixed populations. A
surface reaction, a small cluster or a catalytic site can have populations of a
few molecules, where fluctuations and rare events matter. kMC samples the same
master equation stochastically:

$$
P(\mathbf{n}, t) \rightarrow P(\mathbf{n}, t+\tau),
$$

with integer counts and a waiting time drawn from the total propensity.

## Rates and propensities

Each irreversible step `j` has an HTST forward rate `k_j`. Its propensity is

$$
a_j = k_j \prod_r \binom{n_r}{\nu_{rj}},
$$

where `n_r` is the count of reactant `r` and `ν` its stoichiometric
coefficient. The binomial factor is the number of ways to choose the reacting
molecules from the available population.

The total propensity is

$$
a_0 = \sum_j a_j.
$$

The waiting time is

$$
\tau = -\frac{\ln u}{a_0},
\qquad u \sim \mathrm{Uniform}(0,1),
$$

and the firing step is chosen with probability `a_j / a_0`.

## Algorithm

1. Evaluate every HTST rate at the requested temperature and pressure.
2. Compute the propensity of every step from the current integer counts.
3. Draw `τ` and the firing step.
4. Subtract reactants, add products, record the event.
5. Repeat until `max_events` or `max_time_s`.

Reversible chemistry is represented by two irreversible steps. The solver does
not guess a reverse direction from a reaction free energy; that keeps the
network explicit and auditable.

## API

```rust
use sci_form::kinetics::{ElementaryStep, ThermodynamicState};
use sci_form::kmc::{KmcConfig, KmcSpecies, KmcStep};
use sci_form::run_kinetic_monte_carlo;

let config = KmcConfig {
    species: vec![
        KmcSpecies { id: "A".into(), initial_count: 100 },
        KmcSpecies { id: "B".into(), initial_count: 0 },
    ],
    steps: vec![KmcStep {
        step: ElementaryStep {
            step_id: "A->B".into(),
            activation_free_energy_ev: 0.05,
            reaction_free_energy_ev: -0.1,
            prefactor_s_inv: Some(1.0e12),
        },
        reactants: vec![(0, 1)],
        products: vec![(1, 1)],
    }],
    state: ThermodynamicState { temperature_k: 500.0, pressure_bar: 1.0 },
    max_events: 10_000,
    max_time_s: 10.0,
    seed: 7,
    frame_stride: 100,
};

let result = run_kinetic_monte_carlo(&config)?;
println!("{} events, final A={}, B={}",
    result.events.len(), result.final_counts[0], result.final_counts[1]);
```

### Output

| Type | Meaning |
|---|---|
| `KmcEvent` | Event index, time, firing step, propensity and counts |
| `KmcFrame` | Sampled time, event index, species and counts |
| `KmcResult` | Events, frames, final counts, final time, stop reason and notes |

Time is in seconds, rates in s⁻¹ and counts are dimensionless integers.

## Validation

| Gate | Observable | Result |
|---|---|---|
| Mass conservation | `A + B` after irreversible conversion | unchanged |
| High barrier | Events in a 1 ps window | none; stops on time |
| Binomial coefficient | `C(5,0)`, `C(5,1)`, `C(5,2)`, `C(2,3)` | 1, 5, 10, 0 |
| Invalid species index | Network validation | rejected |

These are algorithmic gates. A quantitative comparison against a reference kMC
code or a catalytic experiment is part of the promotion work.

## Bindings

| Surface | Entry point |
|---|---|
| CLI | `sci-form kmc '<KmcConfig-json>'` |
| WASM | `run_kinetic_monte_carlo(config_json)` |
| Python | `run_kinetic_monte_carlo(config_json)` |

Python returns the JSON encoding of `KmcResult`.
## Limitations

- The network is supplied by the caller; automatic reaction discovery is not
  implemented yet.
- Only irreversible steps are accepted. Encode reversible steps explicitly.
- Counts are integers; fractional populations and continuous concentrations
  belong to the deterministic solver.
- The propensity uses binomial consumption, not a continuum mass-action term.
- There is no spatial resolution, diffusion or surface lattice. A site-based
  kMC is a follow-up.
- The RNG is seeded and deterministic for a single thread; parallel trajectories
  need one seed per replica.

## Next steps

1. Add rule-based network generation from SMIRKS and graph transforms.
2. Add site-based surface kMC with adsorption, desorption and diffusion.
3. Add uncertainty propagation from the barrier distribution.
4. Compare against the deterministic microkinetic solver in the fast-rate limit.
5. Add Python, WASM and CLI wrappers.
6. Add a curated catalytic network benchmark.
