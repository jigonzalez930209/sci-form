# Migrating from 0.15.x to 1.0.0

`sci-form` 1.0.0 is the first stable release. It closes two exhaustive audits
(100 findings in [`reports/auditoria-total.md`](../../reports/auditoria-total.md),
58 blockers in [`reports/v1.0-blockers.md`](../../reports/v1.0-blockers.md)).

Almost everything in this guide is **additive**: the deprecated names still
work, they just warn. The genuinely breaking changes are the ones where a
function that used to **panic** now returns a `Result`, because a library that
crashes on a malformed molecule is not a library you can ship.

---

## The 60-second version

| If you | Do this |
|---|---|
| Call a Python function | Nothing. Every old name still works, with a `DeprecationWarning`. |
| Call a Rust function that used to panic | Add `?` (or `.unwrap()`) — see [Rust](#rust) below. |
| `match` on a public enum | Add a `_ =>` arm. 57 enums are now `#[non_exhaustive]`. |
| Use the CLI | Stop passing `""` / `"[]"` / `0.0` to mean "absent". |
| Compare energies across methods | Read the [units](#units) section. This one can silently change your answer. |
| Compare `total_energy` between PM3 and SCF | **Stop.** One is eV, the other is Hartree. |

---

## Rust

### Functions that now return a `Result`

These used to panic on bad input. They now return `Result<_, String>`.

| Function | Before | After |
|---|---|---|
| `sci_form::embed_batch` | `Vec<ConformerResult>` | `Result<Vec<ConformerResult>, String>` |
| `sci_form::dos::compute_dos` | `DosResult` | `Result<DosResult, String>` |
| `sci_form::dos::compute_pdos` | `PdosResult` | `Result<PdosResult, String>` |
| `sci_form::eht::evaluate_orbital_on_grid` | `Vec<f64>` | `Result<Vec<f64>, String>` |
| `sci_form::eht::evaluate_orbital_on_grid_parallel` | `Vec<f64>` | `Result<Vec<f64>, String>` |
| `sci_form::eht::volume::compute_grid_extents` | `([f64; 3], [usize; 3])` | `Result<([f64; 3], [usize; 3]), String>` |
| `sci_form::eht::volume::evaluate_orbital_on_grid_chunked` | `Vec<f64>` | `Result<Vec<f64>, String>` |
| `sci_form::eht::volume::evaluate_density_on_grid` | `Vec<f64>` | `Result<Vec<f64>, String>` |
| `sci_form::xtb::gfn0::gfn0_repulsion` | `f64` | `Result<f64, String>` |
| `sci_form::xtb::gfn0::gfn0_eeq` | `Vec<f64>` | `Result<Vec<f64>, String>` |
| `sci_form::xtb::gfn0::gfn0_h0` | `Vec<f64>` | `Result<Vec<f64>, String>` |
| `sci_form::scf::uhf::run_uhf` | `UhfResult` | `Result<UhfResult, String>` |
| `sci_form::scf::uhf::run_uhf_parallel` | `UhfResult` | `Result<UhfResult, String>` |
| `sci_form::scf::uhf::run_rohf` | `UhfResult` | `Result<UhfResult, String>` |
| `sci_form::materials::simulate_powder_xrd` | `Vec<f64>` | `Result<Vec<f64>, String>` |
| `sci_form::ml::compute_3d_descriptors` | `Vec<f64>` | `Result<Vec<f64>, String>` |

```rust
// 0.15
let mesh = sci_form::eht::evaluate_orbital_on_grid(&basis, &coeffs, 0, &pos, 0.2, 3.0);

// 1.0
let mesh = sci_form::eht::evaluate_orbital_on_grid(&basis, &coeffs, 0, &pos, 0.2, 3.0)?;
```

### New `try_*` variants

Every remaining panic-on-bad-input entry point now has a fallible sibling.
The infallible one still exists and still panics, so nothing breaks — but you
should prefer the `try_*` form in new code.

| Infallible (panics) | Fallible (new) |
|---|---|
| `scf::types::MolecularSystem::from_angstrom` | `::try_from_angstrom` |
| `scf::types::MolecularSystem::from_flat_angstrom` | `::try_from_flat_angstrom` |
| `scf::scf_loop::run_scf` | `::try_run_scf` |
| `spectroscopy::hessian::compute_hessian` | `::try_compute_hessian` |
| `spectroscopy::hessian::compute_hessian_parallel` | `::try_compute_hessian_parallel` |
| `spectroscopy::ir_intensities::compute_ir_intensities` | `::try_compute_ir_intensities` |
| `spectroscopy::stda_uvvis::compute_stda` | `::try_compute_stda` |
| `ani::neighbor::CellList::new` | `::try_new` |
| `alignment::compute_rmsd` | `::try_compute_rmsd` |
| `alignment::align_coordinates` | `::try_align_coordinates` |
| `alignment::align_quaternion` | `::try_align_quaternion` |
| `gpu::context::GpuContext::new` | `::try_create` |

### Renamed root functions (old names are deprecated aliases)

The non-idiomatic `get_` prefix is gone.

```rust
// 0.15
let support = sci_form::get_eht_support(&elements);

// 1.0
let support = sci_form::eht_support(&elements);
```

| Deprecated | Canonical |
|---|---|
| `get_eht_support` | `eht_support` |
| `get_system_capabilities` | `system_capabilities` |
| `get_system_method_plan` | `system_method_plan` |

The old names still compile and emit a deprecation warning.

### New at the crate root

```rust
let r = sci_form::compute_gfn1(&elements, &positions)?;   // GFN1-xTB
let r = sci_form::compute_gfn2(&elements, &positions)?;   // GFN2-xTB
```

### `#[non_exhaustive]` on the public enums

57 of the 60 public enums are now `#[non_exhaustive]`. A downstream `match`
must have a wildcard arm:

```rust
// 0.15
match method { DosMethod::Tetra => .., DosMethod::Gauss => .. }

// 1.0
match method {
    DosMethod::Tetra => ..,
    DosMethod::Gauss => ..,
    _ => ..,   // required: the variant set may grow
}
```

Three enums were deliberately **not** marked, because their variant set is a
published taxonomy rather than an extensible API: `Mmff94AtomType` (it
mirrors the fixed MMFF94 atom-type list), `ConformerSearchStatus` (a
lifecycle enum) and `Mmff94BondType`.

### New result fields

All additive, so struct literal construction still works if you use `..`:

| Struct | Field | Meaning |
|---|---|---|
| `pm3::Pm3Result` | `notes: Vec<String>` | e.g. "the faithful solver failed, fell back to the legacy path" |
| `xtb::XtbResult` | `notes: Vec<String>` | e.g. "the GPU path was unavailable, the CPU was used" |
| `ir::VibrationalAnalysis` | `n_imaginary_modes: usize` | number of imaginary (negative) frequencies |
| `alignment::AlignmentResult` | `error: Option<String>` | the reason the alignment failed |
| `scf::types::ScfResult` | `error: Option<String>` | the reason the SCF did not converge |
| `hf::cis::CisResult` | `error: Option<String>` | the reason the CIS did not converge |

`AlignmentResult` also gained `is_valid()` and `AlignmentResult::INVALID_RMSD`
(the sentinel `f64::NAN` replaced by an explicit invalid value).

---

## Python

**Nothing breaks.** Every old name still works and emits a `DeprecationWarning`
that names its replacement.

```python
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)   # opt out
```

### The canonical names

The canonical name of a capability is now the **Rust root function name**, and
every action is a `compute_*` / `analyze_*` / `predict_*` verb. The
`*_calculate` suffix had no precedent in RDKit, PySCF, ASE or xtb and is gone.

| Deprecated | Canonical |
|---|---|
| `charges` | `compute_charges` |
| `sasa` | `compute_sasa` |
| `population` | `compute_population` |
| `dipole` | `compute_dipole` |
| `dos` | `compute_dos` |
| `bond_orders` | `compute_bond_orders` |
| `frontier_descriptors` | `compute_frontier_descriptors` |
| `fukui_descriptors` | `compute_fukui_descriptors` |
| `reactivity_ranking` | `compute_reactivity_ranking` |
| `uff_energy` | `compute_uff_energy` |
| `mmff94_energy` | `compute_mmff94_energy` |
| `empirical_pka` | `compute_empirical_pka` |
| `ml_descriptors` | `compute_ml_descriptors` |
| `ml_predict` | `predict_ml_properties` |
| `nmr_shifts` | `compute_nmr_shifts` |
| `nmr_couplings` | `compute_nmr_couplings` |
| `nmr_spectrum` | `compute_nmr_spectrum` |
| `giao_nmr` | `compute_giao_nmr` |
| `hose_codes` | `compute_hose_codes` |
| `ir_spectrum` | `compute_ir_spectrum` |
| `vibrational_analysis` | `compute_vibrational_analysis` |
| `gb_solvation` | `compute_gb_solvation` |
| `alpb_solvation` | `compute_alpb_solvation` |
| `nonpolar_solvation` | `compute_nonpolar_solvation` |
| `pm3_calculate` | `compute_pm3` |
| `xtb_calculate` | `compute_xtb` |
| `gfn1_calculate` | `compute_gfn1` |
| `gfn2_calculate` | `compute_gfn2` |
| `hf3c_calculate` | `compute_hf3c` |
| `ani_calculate` | `compute_ani` |
| `eht_calculate` | `compute_eht` |

The full table of 68 aliases is in `crates/python/src/naming.rs` and the
canonical list, grouped by capability, is `sci_form.__all__`.

### Class names lost the `Py` suffix

The suffix existed only to avoid a Rust name clash; it never helped the Python
caller.

```python
# 0.15
sci_form.ChargeResultPy
sci_form.EhtResultPy
sci_form.SasaResultPy

# 1.0
sci_form.ChargeResult
sci_form.EhtResult
sci_form.SasaResult
```

### `parse` returns an object, not a dict

```python
# 0.15
mol = sci_form.parse("CCO")
mol["num_atoms"]        # dict access

# 1.0
mol = sci_form.parse("CCO")
mol.num_atoms           # attribute access
```

### The error taxonomy

| Exception | When |
|---|---|
| `ValueError` | your input is malformed (a bad SMILES, a length mismatch, an unknown method name) |
| `RuntimeError` | the engine failed (the SCF did not converge, the GPU was unavailable) |
| `NotConvergedError` | specifically a non-convergence (a subclass of `RuntimeError`) |

```python
from sci_form import NotConvergedError

try:
    sci_form.compute_scf(elements, coords)
except NotConvergedError:
    ...   # tighten the SCF, or accept the `notes`
except ValueError as e:
    ...   # fix your input
```

### `cpm_charges` was being computed at the wrong potential

The Python module registered `cpm_charges` **twice**; the second registration
(β tier, μ = 0.0 eV) silently won over the first (μ = −4.44 eV, the standard
H⁺/H₂ electrode potential). In 1.0 the β copy is `beta_compute_cpm_charges`
and `cpm_charges` is the standard one.

```python
sci_form.cpm_charges(elements, coords)              # mu = -4.44 eV (correct)
sci_form.beta_compute_cpm_charges(elements, coords)  # mu = 0.0 eV (the beta tier)
```

### Zero-copy entry points

Two functions write into a caller-owned `float64` buffer instead of copying
the data across the FFI boundary:

```python
import numpy as np

buf = np.empty(3 * n_atoms)
sci_form.embed_coords_into("CCO", buf, seed=42)
```

Use them for a volumetric grid or a long trajectory; for a single molecule the
copy is irrelevant.

### Options mappings

The eight worst-over-argumented functions (`optimize_framework` had 14
parameters) accept an optional `options=` mapping. The positional form still
works.

```python
sci_form.optimize_framework(smiles, coords, options={"max_steps": 200, "forcefield": "uff"})
```

---

## TypeScript / WASM

### The canonical names

Same rule: the Rust root name is canonical. The old WASM names still work.

### `*_seed64`

The 32-bit `seed` argument silently truncated a 64-bit hash, so "the same
seed" gave a different conformer ensemble in JS and in Rust. Every seeded
function now has a `_seed64` variant that takes a `bigint`.

```ts
import { embedConformerSeed64 } from "sci-form-wasm";

const out = embedConformerSeed64("CCO", coords, 42n);
```

The 32-bit variants are kept for compatibility.

### `parse_smiles` returns the `atoms` array

It used to return only the counts. It now matches Python's `parse`.

### The error payload is always valid JSON

```ts
const raw = computeDose(elements, coords);
// '{"error":"..."}' — always parseable
const parsed = JSON.parse(raw);
```

---

## CLI

### One sentinel convention

An optional argument is now simply **absent**. The old sentinels (`""`,
`"[]"`, `0.0` meaning "engine default", `None`) are gone.

```bash
# 0.15 — "" meant "no coordinates"
sci-form stereo --smiles "C[C@H](N)C(=O)O" --coords ""

# 1.0 — just omit the flag
sci-form stereo --smiles "C[C@H](N)C(=O)O"
```

### One error-message style

```
<What failed> '<actual value>'. Supported: 'a', 'b', 'c'.
```

The message always quotes the offending value, always lists the valid values,
and never leaks a Rust internal name, a `None`, or a NaN.

### The developer binaries are gone from the install

`cargo install sci-form` used to install **28 executables** on your `$PATH`
(`benchmark`, `test_eigen`, `test_failures`, `check_grad`, `print_c_eri`, …).
Only the `sci-form` CLI is installed now. The sources remain in the repository
if you want to build them.

---

## Units

**This is the section that can silently change your answer.**

Before 1.0, `total_energy` meant three different things:

| Result struct | `total_energy` is |
|---|---|
| `DftResult`, `D4Result`, `ScfResult`, `UhfResult` | Hartree |
| `Pm3Result`, `XtbResult`, `Gfn1Result`, `Gfn2Result` | eV |
| `AlpbResult` | kcal/mol |

Comparing two of them was off by 27.2× (eV → Hartree) or 23.06×.

1.0 does **not** rename the fields — that would break every script at once.
Instead:

- Every result field now states its unit in the rustdoc, and the crate header
  carries a per-method unit table.
- The `_ev` / `_hartree` / `_kcal_mol` suffix convention is used for every
  *new* field.

If you compare across methods, convert explicitly:

```rust
const HARTREE_TO_EV: f64 = 27.211386245988;
let pm3_ev = pm3.total_energy;                    // already eV
let scf_ev = scf.total_energy * HARTREE_TO_EV;     // Hartree -> eV
```

And note the related trap the audit found: `Hf3cResult.homo_energy`,
`lumo_energy`, `gap` and `orbital_energies` are **Hartree**, even though they
were documented as eV in 0.15.

---

## Packaging

- The published crate is 339 files / 7.7 MiB (was 612 / 20.6 MiB).
- The 27 developer binaries are no longer installed by `cargo install sci-form`.
- The MSRV is now declared: **Rust 1.87**.
- `cargo install sci-form-cli` no longer compiles `nalgebra` (it is behind the
  `alpha-periodic-linear` feature).

---

## See also

- [`CHANGELOG.md`](../../CHANGELOG.md) — the full 1.0.0 entry
- [`reports/v1.0-blockers.md`](../../reports/v1.0-blockers.md) — the 58 blockers
  and their disposition
- [`reports/auditoria-resolucion.md`](../../reports/auditoria-resolucion.md) —
  the 100 findings of the first audit
- [`docs/api/complete-supported-reference.md`](../api/complete-supported-reference.md)
  — every capability, by surface
