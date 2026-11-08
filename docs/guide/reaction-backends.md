# Reaction Simulation Backends

Reaction simulation in sci-form now has three separate layers that should not be conflated.

1. Molecule and reaction interpretation.
   This is the SMILES, SMARTS, and SMIRKS layer. It defines which atoms, bonds, mappings, and formal product transforms are being considered.
2. Path construction.
   This is where alpha GSM builds a path between reactant and product geometries and searches for a transition-state candidate.
3. Energy surface selection.
   This is the backend that assigns a scalar energy to each geometry along the path.

That split matters because reference libraries usually specialize in different layers.

- RDKit is the main reference for reaction SMARTS and transform semantics.
- Open Babel is a useful chemistry toolkit reference, but it is not the primary baseline for reaction SMARTS semantics in this repo.
- geomeTRIC is a useful chain-of-states and NEB reference for path-search behavior.
- xTB is a reference family for fast semiempirical/tight-binding reaction-path work.
- PySCF is a reference for SCF behavior and ab initio single-point expectations.

## Backend policy for GSM

The GSM backend layer normalizes energies to kcal/mol internally so barrier comparisons remain meaningful across heterogeneous methods.

Available backends are grouped by role:

- Exploratory force fields: UFF, MMFF94.
- Barrier-capable electronic methods: PM3, xTB/GFN0, GFN1-xTB, GFN2-xTB, HF-3c.
- Reactive force field: ReaxFF, currently only through the alpha CHON parameter set when the `alpha-reaxff` feature is enabled.
- Orbital-only analysis: EHT.

## Energy surface and gradient of each backend

<SvgDiagram src="/svg/energy-surface-map.svg" alt="Which energy surface each backend name evaluates and how its gradient is obtained" />

| Backend name | Energy | Gradient (NEB, MD, reaction dynamics) |
|---|---|---|
| `uff`, `mmff94` | classical force field | analytic |
| `pm3` | MOPAC-parity PM3 for `Z ≤ 17`, legacy NDDO otherwise | stationary-density gradient: fixed-density integral rebuilds, no SCF per displacement |
| `xtb` (legacy `gfn0`) | GFN2 (`solve_gfn2`) on every surface, identical to `gfn2` | warm-started central differences of the GFN2 energy (GSM evaluates energies only) |
| `gfn1`, `gfn2` | `solve_gfn1`, `solve_gfn2` | central differences |
| `hf3c` | HF-3c | central differences |

> **Absolute energies of different backends are not comparable.** For ethanol
> `neb-energy` returns −7146 kcal/mol with `xtb`/`gfn2` and −14379 kcal/mol with
> `pm3`. Compare energies only within one backend (differences along a path).
> Before v1.0, `xtb` in NEB/MD was a simplified model that returned +9117 kcal/mol.

The one rule that matters: for a given backend the gradient is always the
derivative of the energy it is paired with. That was not true for `pm3` and
`xtb` before the gradients were rebuilt (see
[PM3 and xTB gradients](/algorithm/semiempirical-gradients)). Energy-only calls
(`neb-energy`, the final image energies of the simplified NEB) do not compute a
gradient. `"xtb"` and `"gfn2"` are the same surface on every tool.

## Why EHT is not exposed as a barrier backend

EHT is retained in the planning API because users may still want to inspect whether the system is covered by the current orbital parameterization. It is intentionally marked as not suitable for GSM reaction-path energies.

Reason:

- EHT in sci-form currently provides orbital-level observables and support metadata.
- It does not provide a calibrated total-energy surface comparable to PM3, xTB, HF-3c, or even classical force fields for transition-state barrier work.
- Exposing it as if it were numerically interchangeable would make reaction barriers look comparable when they are not.

## Current public API

Rust alpha GSM now exposes:

- `plan_gsm_backends(elements)`
- `plan_gsm_backends_for_smiles(smiles)`
- `evaluate_gsm_backend(smiles, coords, backend)`
- `compare_gsm_backends(smiles, coords, backends)`
- `find_transition_state_with_backend(smiles, reactant, product, backend, config)`

Python alpha bindings expose:

- `gsm_backend_plan(smiles)`
- `gsm_compare_backends(smiles, coords, methods=None)`
- `gsm_find_ts(smiles, reactant_coords, product_coords, n_nodes=9, method="uff")`

WASM alpha bindings expose:

- `alpha_gsm_backend_plan(smiles)`
- `alpha_gsm_compare_backends(smiles, coordsFlatJson, methodsJson)`
- `alpha_gsm_find_ts_with_method(smiles, reactantJson, productJson, configJson, method)`

The legacy `alpha_gsm_find_ts(...)` entry point remains as a UFF wrapper for backward compatibility.

## Practical comparison guidance

For the kind of benchmarking already used elsewhere in sci-form, the useful comparison split is:

- Reaction interpretation and products: compare against RDKit reaction SMARTS behavior.
- Path shape and TS localization: compare against geomeTRIC/NEB-style expectations on curated reaction profiles.
- Single-point energetics along the path: compare PM3, xTB/GFN, HF-3c, and force fields on the same geometries in normalized kcal/mol.

This keeps the comparison honest: reaction semantics, path construction, and energy models are evaluated separately before being combined into a full simulation workflow.