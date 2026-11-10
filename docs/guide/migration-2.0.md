# Migrating from 1.x to 2.0

`sci-form` 2.0 is the release that completes the MMFF94 force field, adds typed
errors to the root Rust API, and moves the validated experimental modules into
core. The old paths and feature flags still compile; the migration is mostly a
set of renames plus one error-type change.

> Coming from 0.15? Read the historical
> [0.15 → 1.0 guide](./migration-1.0.md) first.

## The 60-second version

| If you | Do this |
|---|---|
| Use the root Rust API | Change `Result<_, String>` to `Result<_, SciFormError>`, or keep using `?` and the `Display`/`Deref<str>` impls |
| Pattern-match `Err(String)` | Match `SciFormError::Parse / Validation / Geometry / Numerical / Unsupported / Other` |
| Import `alpha::gsm`, `alpha::kinetics`, `alpha::periodic_linear` | Import `gsm`, `kinetics`, `periodic_linear` |
| Import `beta::kpm`, `beta::rand_nla`, `beta::cpm`, `beta::riemannian` | Import `kpm`, `rand_nla`, `cpm`, `riemannian` |
| Import `alpha::cga`, `alpha::sdr`, `alpha::edl`, `alpha::render_bridge` | Import `cga`, `sdr`, `edl`, `render_bridge` |
| Enable `alpha-gsm`, `alpha-kinetics`, `alpha-periodic-linear` | Remove the flag for the core module; keep it only if an old alias path needs it |
| Enable `beta-kpm` or `beta-randnla` | Remove the flag; the modules are unconditional |
| Use MBH, DFT, MLFF or reaction dynamics | Keep the matching `alpha-*` / `beta-*` flag; those remain gated |
| Use Python, WASM or CLI | Bump the package to 2.0.0; old names remain available |

## 1. The root error type

### Before

```rust
use sci_form::parse;

fn load(smiles: &str) -> Result<(), String> {
    let molecule = parse(smiles)?;
    let _ = molecule;
    Ok(())
}
```

### After

```rust
use sci_form::{parse, SciFormError};

fn load(smiles: &str) -> Result<(), SciFormError> {
    let molecule = parse(smiles)?;
    let _ = molecule;
    Ok(())
}
```

`SciFormError` is `#[non_exhaustive]` and has six variants:

| Variant | Meaning |
|---|---|
| `Parse(String)` | A SMILES, CIF or JSON input could not be parsed |
| `Validation(String)` | The input was well-formed but failed validation |
| `Geometry(String)` | Mismatched element/coordinate lengths or non-finite coordinates |
| `Numerical(String)` | SCF, EEM, decomposition or another numerical method failed |
| `Unsupported(String)` | An unsupported feature, element or parameter set |
| `Other(String)` | Any other failure |

The type implements `Display`, `std::error::Error`, `Deref<Target = str>`,
`AsRef<str>`, `From<String>` and `From<&str>`. That means most code that only
prints the error does not change:

```rust
match sci_form::parse("Q") {
    Ok(_) => unreachable!(),
    Err(error) => eprintln!("could not parse: {error}"),
}
```

Use the variants when the reaction differs:

```rust
use sci_form::SciFormError;

match sci_form::parse(input) {
    Ok(molecule) => molecule,
    Err(SciFormError::Parse(message)) => {
        eprintln!("bad SMILES: {message}");
        return Ok(None);
    }
    Err(SciFormError::Unsupported(message)) => {
        eprintln!("unsupported element: {message}");
        return Ok(None);
    }
    Err(other) => return Err(other),
}
```

If you have a helper that still returns `String`, convert at the boundary:

```rust
let molecule = sci_form::parse(input).map_err(|error| error.to_string())?;
```

The reverse conversion also exists for incremental migrations:

```rust
let typed: sci_form::SciFormError = "invalid geometry".to_string().into();
```

> The 65 root-level functions re-exported from `src/lib.rs` return
> `SciFormError`. A few deep module APIs still return `String` for legacy
> reasons; the `From<String>` implementation classifies those messages.

## 2. Promoted modules and old paths

These modules were validated in the 1.x experimental tier and are core in 2.0.
The old path is a compatibility alias, not a second implementation.

| Old path | Recommended path | Old feature flag | Status |
|---|---|---|---|
| `alpha::cga` | `cga` | `alpha-cga` | core, alias |
| `alpha::sdr` | `sdr` | `alpha-sdr` | core, alias |
| `alpha::edl` | `edl` | `alpha-edl` | core, alias |
| `alpha::gsm` | `gsm` | `alpha-gsm` | core, alias |
| `alpha::kinetics` | `kinetics` | `alpha-kinetics` | core, alias |
| `alpha::periodic_linear` | `periodic_linear` | `alpha-periodic-linear` | core, alias |
| `alpha::render_bridge` | `render_bridge` | `alpha-render-bridge` | core, alias |
| `beta::kpm` | `kpm` | `beta-kpm` | core, unconditional alias |
| `beta::rand_nla` | `rand_nla` | `beta-randnla` | core, unconditional alias |
| `beta::cpm` | `cpm` | `beta-cpm` | core, alias |
| `beta::riemannian` | `riemannian` | `beta-riemannian` | core, alias |

The core paths need no flag. The old `alpha::*` / `beta::*` aliases still
need their historical flag unless the table says "unconditional alias".

The modules that remain gated are:

| Module | Path | Feature flag |
|---|---|---|
| Kohn-Sham density functional theory | `dft` | `alpha-dft` |
| Machine-learned force field | `mlff` | `alpha-mlff` |
| Mobile Block Hessian | `beta::mbh` | `beta-mbh` |
| Full 3D reaction dynamics | `alpha::reaction_dynamics` | `alpha-reaction-dynamics` |

### Rust import migration

```rust
// Before
use sci_form::alpha::gsm::{find_transition_state, GsmConfig};
use sci_form::beta::rand_nla::solve_eht_randnla;

// After
use sci_form::gsm::{find_transition_state, GsmConfig};
use sci_form::rand_nla::solve_eht_randnla;
```

`dynamics_live` is also unconditional now. The historical
`alpha-dynamics-live` flag is a no-op alias, and the old `imd`, `kmesh`,
`forcefield::reaxff` and `scf::obara_saika` paths were already core before 2.0.

## 3. Cargo features

The default build contains the promoted modules. Remove obsolete flags from
your `Cargo.toml` unless you still need an old alias path.

```toml
# Before: the old flags selected the implementation
[dependencies]
sci-form = { version = "1", features = [
    "alpha-gsm", "alpha-kinetics", "alpha-periodic-linear",
    "beta-kpm", "beta-randnla"
] }

# After: the implementation is core
[dependencies]
sci-form = { version = "2.0", features = ["parallel"] }
```

The only flags that still select research code are:

```toml
[dependencies]
sci-form = { version = "2.0", features = [
    "alpha-dft", "alpha-mlff", "alpha-reaction-dynamics", "beta-mbh"
] }
```

Two dependency details changed:

- `alpha-periodic-linear` no longer pulls in KPM or RandNLA; both are core.
- `alpha-kinetics` now depends on `beta-mbh` because the optional
  GSM → MBH → transition-state adapter still uses the gated MBH module.

## 4. What changed in the science

The migration is not only renames. The 2.0 line also changes numbers and
coverage:

| Area | 2.0 change |
|---|---|
| MMFF94 | Complete energy expression: stretch-bend, out-of-plane, buffered 14-7 vdW and buffered Coulomb with PBCI/BCI charges; RDKit parity 2.8e-11 kcal/mol |
| ReaxFF | Analytical gradients for bond orders, bond/angle, over/under-coordination, vdW, Coulomb and EEM charge response |
| RandNLA | Re-orthonormalized power iterations and an N>500 low-rank path (9.6e-4 at N=600) |
| Cholesky ERI | `hf::integrals::CholeskyEris` reconstructs the packed ERI tensor from a pivoted Cholesky factor |
| PM3 / xTB | PM3 stationary-density gradient (no SCF per displacement); xTB warm-started central differences |
| NMR | `HoseDatabase` loader for external JSON/CSV/TSV exports; no NMRShiftDB2/SDBS data is bundled |
| Errors | Root API returns `SciFormError` instead of `String` |
| Validation | 129 recorded rows, 122 PASS, 7 method/provenance gaps, 0 FAIL |

The complete per-method report is in
[Scientific validation](/validation/scientific-validation). The raw rows are in
the [full matrix](/validation/full-matrix).

## 5. Python, WASM and CLI

- Package versions are `2.0.0` (`sciforma`, `sci-form-wasm`, `sci-form-cli`).
- Old Python and TypeScript names remain available with their deprecation
  aliases. The canonical names are the Rust root names.
- The `sci_form.alpha` and `sci_form.beta` submodules still exist; prefer the
  top-level function or the core module path in new code.
- CLI commands that wrap a promoted module still carry the historical feature
  flag in the build. The flag gates the wrapper, not the Rust implementation.
- MBH, DFT, MLFF and reaction dynamics remain genuinely gated.

## 6. Migration checklist

1. Change the dependency to `sci-form = "2.0"` (or `sci-form = { version = "2.0", ... }`).
2. Replace old module paths with the core paths from the table above.
3. Remove feature flags for promoted modules; keep only the four research flags.
4. Let `?` convert `String` errors through `From<String>`, or match
   `SciFormError` variants explicitly.
5. If you match a public enum, keep the wildcard arm: `SciFormError` is
   `#[non_exhaustive]`, and most public enums are too.
6. Re-run the validation gates relevant to your use case:

```bash
cargo test --release --test cross_validation
cargo test --release --test alpha_gsm_validation
cargo test --release --test beta_randnla_validation
cargo test --release --test alpha_edl_validation
```

The repository MSRV is **1.89**.

## 7. If you want the old behaviour

Every promoted path is still available. If a build enables an old flag, the
alias resolves to the core implementation; it does not re-enable a separate
code path. The only modules with a distinct gated implementation are the four
research lines listed in §2.

For a clean 2.0 build, remove the old flags. For a short transition period, keep
them; they are compatibility aliases and are not scheduled for removal in the
2.0 line.
