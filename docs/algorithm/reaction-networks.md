# Rule-based reaction-network discovery

The v2.3 reactivity prototype adds a rule-based discovery engine on top of the
existing SMILES canonicalisation and Gillespie kMC solver. The engine applies
explicit rules, deduplicates canonical species, records directed elementary
edges and converts the result into a kMC configuration.

> **Status.** Implemented in `feat/wave-c`. Rules are explicit templates; SMARTS
> pattern generation and automatic barrier estimation remain follow-ups. See the
> [v2.x roadmap](/v2-roadmap).

## Why explicit rules

Automatic reaction discovery has two separate parts:

1. deciding which transformations are chemically plausible;
2. closing the graph over species and deduplicating.

This prototype implements part 2 and a clean rule contract. A rule owns the
chemistry; the engine owns canonicalisation, BFS closure, depth and species
limits. That makes a curated mechanism testable without hiding a heuristic.

## Rule contract

`ReactionRule` exposes:

| Method | Meaning |
|---|---|
| `name` | Stable rule identifier |
| `reactant_count` | 1 or 2 in the current engine |
| `barrier_ev` | Activation barrier |
| `prefactor_s_inv` | Optional prefactor |
| `apply` | Return canonical product sets |

`TemplateRule` is the built-in explicit implementation:

```rust
use sci_form::networks::TemplateRule;

let rule = TemplateRule::new(
    "C1",
    &["C"],
    &["CC"],
    0.05,
    Some(1.0e12),
)?;
```

Reactants and products are canonicalised at construction, so `"CCO"` and
`"OCC"` are the same species.

## Discovery

```rust
use sci_form::networks::{NetworkDiscoveryConfig, discover_reaction_network};

let config = NetworkDiscoveryConfig {
    initial_species: vec!["C".into()],
    rules: vec![Box::new(rule)],
    max_depth: 3,
    max_species: 10,
};
let network = discover_reaction_network(&config)?;
```

The engine runs BFS for `max_depth` layers. At each layer it applies
unimolecular rules to every species and bimolecular rules to every unordered
pair. New products are canonicalised, deduplicated and added until no new
species appear or `max_species` is reached.

### Result

| Type | Meaning |
|---|---|
| `NetworkSpecies` | Index and canonical SMILES |
| `NetworkEdge` | Reactants, products, rule, barrier and prefactor |
| `ReactionNetwork` | Species, edges and discovery notes |

## kMC integration

```rust
use sci_form::kinetics::ThermodynamicState;

let kmc = network.to_kmc_config(
    &[100, 0, 0],
    ThermodynamicState { temperature_k: 500.0, pressure_bar: 1.0 },
    1_000,
    10.0,
    7,
)?;
let result = sci_form::run_kinetic_monte_carlo(&kmc)?;
```

Every edge becomes an irreversible `ElementaryStep`. The reaction free energy is
set to zero because the discovery layer has no thermodynamic data; the barrier
and prefactor come from the rule.

## Bindings

| Surface | Entry point |
|---|---|
| CLI | `sci-form network '<TemplateNetworkConfig-json>'` |
| WASM | `discover_reaction_network(config_json)` |
| Python | pending |

## Validation

| Gate | Observable | Result |
|---|---|---|
| Linear network | C → CC → CCC | 3 species, 2 edges |
| Canonicalisation | `"CCO"` vs `"OCC"` | one species |
| kMC integration | 100 C, 1,000 events | population conserved |
| Invalid product | `"Q"` | rejected |

These are engine gates. A published mechanism and a barrier benchmark are still
required.

## Limitations

- Rules are explicit templates; no SMARTS pattern generation.
- Only unimolecular and bimolecular rules.
- No automatic barrier estimation or thermodynamic consistency.
- No reaction free energies, so kMC uses `ΔG_reaction = 0`.
- The species limit can truncate the network; skipped products are recorded in
  `notes`.
- Bimolecular enumeration is combinatorial but bounded by `max_species`.
- No Python wrapper and no curated benchmark yet.

## Next steps

1. Add a SMARTS/SMIRKS-to-product rule adapter.
2. Add GSM/MBH barrier refinement for each edge.
3. Add thermodynamic cycles and reverse rates.
4. Add a published mechanism benchmark.
5. Add site-based surface kMC.
6. Add Python bindings and a CLI kMC command.
