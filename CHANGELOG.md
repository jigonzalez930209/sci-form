# Changelog

This file captures the project-level highlights from `v0.12.0` onward.
GitHub Releases already exist for each tag, but most release bodies are auto-generated;
the summary below records the substantive work that landed in each release line.
Future releases should append new sections here so this file remains the maintained,
human-written changelog for the repository.

## Unreleased (1.1 line) - 2026-11-07

### Validation

- Re-ran the whole external campaign against pinned references (RDKit
  2026.09.1, PySCF 2.14.0, tblite 0.7.0, dftd4 4.3.0, xtb 6.7.1,
  multicharge 0.6.0, MOPAC 23.2.5, LAMMPS 29Aug2024_update2, ASE 3.29.0,
  pormake 0.2.3): 36/36 oracle scripts and 117/117 Rust validation tests
  pass; `reports/validation/summary.json` holds 127 rows, 120 PASS and 7
  provenance/method-gap rows, 0 FAIL. Python 3.14.4 is the oracle runtime.
  The PySCF AO->MO fixture now canonicalises degenerate MO subspaces.

### Dependencies and GPU

- nalgebra 0.35 (simba 0.10) removes the unmaintained `paste` crate;
  RUSTSEC-2024-0436 is closed and `deny.toml` has no advisory ignores.
  Workspace MSRV is 1.89.
- Native `src/gpu/` and the WASM runtime now share wgpu 29. The duplicated
  naga/ash/wgpu-core trees are gone. `tests/gpu_parity.rs` adds eight real
  GPU-vs-CPU gates (ERI, Fock, gamma, orbital/density grids, ESP, EEQ, CPM)
  and skips cleanly without an adapter (`SCIFORM_REQUIRE_GPU=1` makes that a
  failure). The ERI kernel now rejects p/d/f bases instead of returning
  silently wrong (ss|ss)-only integrals.

### Performance and docs

- Criterion infrastructure (`benches/core.rs`, `[[bench]] core`).
- Runnable doctests rise from 15 to 28 (+1 ignored).
- MMFF94 now implements the **complete** energy expression (stretch-bend,
  out-of-plane and buffered Coulomb electrostatics with PBCI/BCI charges) and
  matches RDKit to <3e-11 kcal/mol on the validation corpus.

### 2026-11-08 — 1.1 follow-up

- **MMFF94 complete**: stretch-bend, out-of-plane and buffered Coulomb
  electrostatics with MMFF partial charges. Seven-molecule RDKit parity:
  max energy deviation 2.8e-11 kcal/mol, partial charges exact. F100 added.
- **ReaxFF analytical gradients**: bond-order scaling, bond/angle,
  over/under-coordination, vdW, Coulomb and EEM charge response
  (`dq = -A^-1 dA q`); gate against central differences on H2, water,
  methane and CH3O.
- **RandNLA N>500**: `solve_eht_randnla_lowrank` (Nyström + reorthogonalized
  Lanczos). N=600, k=160: four lowest states at 9.6e-4 relative error,
  ~80 ms vs ~143 ms exact.
- **ERI Cholesky**: `hf::integrals::CholeskyEris` reconstructs the packed ERI
  tensor from a thresholded pivoted Cholesky factor; benzene STO-3G rank
  237/666, ratio 0.714, max error 9.9e-7 Ha.
- **Typed errors**: `SciFormError` (Parse, Validation, Geometry, Numerical,
  Unsupported, Other) with `Display`/`Error`/`Deref<str>`/`From<String>`;
  the 65 root-level functions in `src/lib.rs` now return it.
- **NMR HOSE**: `HoseDatabase` reads JSON/CSV/TSV exports from
  NMRShiftDB2/SDBS and overrides 1H/13C shifts on a longest-prefix match.
  No database data is bundled (CC-BY-SA review remains external).
- **PM3/xTB gradient cost**: PM3 uses the stationary-density (Pulay-free)
  form (no SCF per displacement); xTB warm-starts each displaced SCC from
  the reference charges/multipoles. Symbolic NDDO/GFN2 derivatives remain
  research work.
- **KPM promoted to core** (`src/kpm`, `beta-kpm` compatibility alias).
- **ReaxFF ffield parser**: all standard sections (general, atoms, bonds,
  off-diagonal, angles, torsions, hydrogen bonds).
- Gate: `scripts/release-check.sh` green.

### Experimental tracks

- RandNLA gains re-orthonormalised subspace iterations; the F93 rank-N/4
  top-k error drops from ~0.13 to 3.8e-4 (<1e-3) at 384 iterations.
- GSM now refines the saddle with 200 steps at `step_size`; five curved
  double wells match the analytic barrier and a discrete NEB within 0.8 %.
- SDR vs Butina is measured with the Adjusted Rand Index (1.0 on 30
  conformers from three shapes).
- Python and WASM integration tests cover the periodic KPM/RandNLA
  adapters. ReaxFF's full ffield parser/analytic gradients and the NMR
  HOSE/3J experimental set remain open.

## v1.0.0 - 2026-11-03

The first stable release. v1.0.0 closes **three exhaustive audits** (100 findings
in `reports/auditoria-total.md`, 88 findings in `reports/v1.0-blockers.md`) and
raises the floor on everything that matters for a chemistry library: the crate
no longer panics on user input, its results are reproducible bit-for-bit
regardless of the machine's core count, every public item is documented, and
the published artefact is 13 MB smaller.

**Upgrade guide:** [`docs/guide/migration-1.0.md`](docs/guide/migration-1.0.md)
lists every breaking change. The short version: the functions that can fail
now return a `Result` instead of panicking, and the canonical name of every
capability is the Rust root name in all four surfaces.

### Security and robustness

- **The library no longer panics on user input.** 164 panicking sites in the
  library code were classified; the 39 reachable from a malformed molecule, a
  degenerate geometry or an unsupported element were eliminated.
  - A new shared validator (`src/validation.rs`) rejects a length mismatch
    between the elements and the coordinates, a non-finite coordinate and an
    empty system with one actionable message, and it guards every public
    entry point that takes a geometry.
  - The same validator family now also checks the parameters that reach a
    kernel: a partial-charge vector of the wrong length, a non-finite charge,
    a dielectric constant outside `ε >= 1.0` and a negative probe radius.
  - All 30 `partial_cmp(...).unwrap()` sites became `unwrap_or(Ordering::Equal)`.
    A single NaN coordinate used to panic the whole library.
  - `Molecule::from_smiles("[")` returns an error instead of panicking, and
    the SMILES parse errors now carry the position and the offending fragment
    instead of a bare `Expected ')'`.
  - `compute_uhf` validates the charge and the multiplicity against the real
    electron count instead of casting a negative `i32` to `usize`.
  - The GFN0 public helpers, `solve_gfn0`, the PM3 faithful solver, the DOS,
    the orbital mesh and the ESP grid validate their parameters. A zero
    spacing or a zero point count used to return an EMPTY result in release
    builds — a debug-only panic let CI stay green while the published crate
    returned garbage.
  - `GasteigerConfig` rejects `max_iter = 0` and a non-positive damping; those
    used to yield all-zero charges while reporting `converged = true`.
  - `embed_batch` clamps the thread count (`1 << 20` used to hang forever) and
    propagates the pool-creation error.
  - The SMARTS parser returns parse errors instead of `panic!`.
- **The WASM error encoder produced invalid JSON.** serde_json embeds the
  offending value in quotes, so `{"error":"invalid type: string "abc"…"}` did
  not parse and `JSON.parse()` threw instead of showing the error. Every WASM
  error path was affected. It now builds the payload with serde_json, and a
  test parses the output back.
- **A panic inside a PyO3 function is a `PanicException`**, which inherits from
  `BaseException` and is therefore NOT caught by `except Exception:`. With the
  panics removed from the core, the Python surface now raises ordinary
  exceptions.

### Fixed

These are the fixes that landed after the 2026-10-22 audit freeze; they are
part of the tagged 1.0.0 build.

- **`pm3-gradient` / `xtb-gradient` (and the NEB/MD `pm3` / `xtb` backends) returned
  wrong forces.** The hand-written Hellmann-Feynman/Pulay expressions
  differentiated a legacy simplified model, not the energy the API reports: they
  deviated from the numerical gradient by up to 18 eV/Å (PM3, H₂O₂) and
  209–263 eV/Å (xTB). The gradients are now central differences (5e-4 Å) of the
  energy the matching `pm3` / `xtb` call returns — the MOPAC-parity PM3 energy
  and the GFN2-xTB energy — so they agree with it to ≤ 1e-4 eV/Å (gates F80).
  `xtb-gradient` now differentiates GFN2, matching `xtb`. Cost is `6·N` SCF
  evaluations per gradient; an analytic NDDO/GFN2 gradient remains future work.
  `compute_semianalytical_hessian` is therefore a double finite difference; prefer
  `compute_numerical_hessian`.
- **PM3 SCF (MOPAC-parity solver) converges with Pulay DIIS** (27 iterations
  instead of the full 200 on H₂O₂). The energy moves by ≲ 1e-5 eV (1e-7 on H₂O₂, 4e-6 on ethanol), but
  the scatter between finite-difference steps drops from ~2e-3 to ~6e-6 eV/Å.
  The damped iteration remains as the fallback, now with a density-residual
  criterion.
- **The KS-DFT gradient was correct.** Its F80 gate compared it with a derivative
  taken per Å (ratio exactly 0.5292 = Å/Bohr). The reference is now converted,
  the gate is enforced (3e-7 Ha/Bohr on H₂O₂) and the module no longer claims a
  missing XC response term.
- **The UKS (open-shell) KS-DFT gradient was wrong; it is fixed.** It reused the restricted
  density, occupation and orbital energies and was off by 1.4e-2 Ha/Bohr (5 %) on the OH
  radical. It now uses the spin densities and a spin-resolved Pulay term (1e-6 Ha/Bohr with
  `energy_threshold = 1e-13`), and a new F80 gate covers it. `DftResult` gains
  `orbital_energies_beta` and `coefficients_beta` (empty when restricted).
- **CI/CD: no tests on GitHub.** `ci.yml` only builds (with the cargo cache). The whole test
  suite runs locally through `scripts/release-check.sh`; `scripts/release.sh` tags only a commit
  that passed it. `release.yml` builds and publishes (crates.io, PyPI, npm) and attaches the CLI
  to the GitHub Release as `.tar.gz` (Linux, macOS), `.zip` (Windows), `.deb` and `.rpm`, with
  `SHA256SUMS`. The documentation is deployed afterwards.
- **A guided CLI.** `sci-form` with no arguments prints a guide; new `commands`, `examples
  [--run]` and `completions` commands; every command's `--help` ends with a ready-to-run example.
  The catalogue is tested (every command listed, every example executed).
- **Live examples in the documentation.** One runnable example for each of the 155 public WASM
  functions, executed in the visitor's browser with the compiled package
  (`/guide/live-examples`), plus examples embedded in the algorithm pages.
- **The 3D reaction path (`alpha-reaction-dynamics`) is Rust-core only.** It is too
  expensive for the bindings (`6·N` SCF per gradient) and is not a validated reaction path,
  so `sci_form.alpha.compute_reaction_dynamics_3d`, the `sci-form-wasm/alpha`
  `compute_reaction_dynamics_3d` / `compute_approach_dynamics` and the CLI command
  `reaction-dynamics3d` are removed, together with the feature in the three binding
  crates. The module, its tests and the feature in the core crate stay.
- **`"xtb"` means GFN2 in NEB, MD and reaction dynamics too.** It used to select a
  simplified model (`xtb::solve_xtb`) with unrelated absolute energies. `NebBackend::Xtb`
  is gone (`from_method("xtb")` → `Gfn2`); `MdBackend::Xtb` evaluates GFN2; `"gfn0"` is a legacy alias.
- Reaction-dynamics sites that only need an energy no longer compute a gradient
  (`neb_energy_kcal`). `tests/reaction_dynamics_study.rs` (17 categories) is now
  `#[ignore]`: with `xtb` = GFN2 it needs more than 10 minutes; run it with
  `cargo test --test reaction_dynamics_study -- --ignored`.
- **`xtb` is GFN2.** `compute_xtb`, `sci-form xtb` and the bindings run GFN2-xTB
  (identical to `gfn2`); the documentation, the CLI `--help` and the README
  called it GFN0. The crate-level unit table said PM3/xTB energies are in
  Hartree; they are in eV (heat of formation in kcal/mol).
- `neb-energy` and the final image energies of the simplified NEB no longer
  compute (and discard) a gradient, which would have cost `6·N` extra SCF runs
  for `pm3`/`xtb`.
- `tests/integration/test_cross_language.py` was updated to the 1.0 API
  (`embed` raises, `parse` returns a `ParseResult`, Node tests run against the
  WASM `pkg-node`); `crates/python/tests/test_bindings.py` and the WASM
  integration suite now run in CI.
- The `topology_data.rs` clippy `approx_constant` error and the stale
  `full_pipeline_ethanol_rearrangement` frame count (a single product has a
  single departure frame).

- **A misspelled method name silently changed the result.** The method names
  were parsed by every binding with its own `match`, and the copies had
  drifted apart: the orbital mesh accepted nine aliases but its error
  message listed four, the DOS multimethod had the same `match` three
  times, and the IR broadening fell back to the **Lorentzian** shape for
  anything that was not `"gaussian"` — so the typo `gauassian` changed the
  line shape without a word. The parsing now lives in the core, next to the
  enum it produces, and every surface reports the same error:
  `Invalid method 'gauassian'. Supported: 'gaussian', 'lorentzian'.`
  `compute_ir_spectrum_broadened` now returns a `Result` instead of
  silently choosing a shape.
- **The published wheel and npm package shipped fewer capabilities than CI
  tested.** The release workflow carried its own copy of the feature lists
  and that copy had fallen behind: the published Python wheel and npm
  package were missing the `alpha-edl`, `alpha-periodic-linear`,
  `alpha-kinetics`, `alpha-render-bridge` and `alpha-reaction-dynamics`
  modules, while CI built and tested all of them. Both workflows now use
  the same lists, and `scripts/check_feature_parity.py` fails the build if
  they ever drift apart again.
- **The TypeScript declarations of the alpha/beta subpaths promised functions
  the binary did not contain.** They were a hand-written file listing every
  research function, independent of the feature set the package was compiled
  with: with the release feature set, 18 of the 36 declared functions were
  missing from the binary. They are now generated from the build itself, so
  they cannot drift. `alpha_modules_info()` likewise reports the modules the
  binary actually contains.
- **The npm package could not be used from Node.js.** Only the browser ESM
  build was ever published: `require('sci-form-wasm')` — the path the
  TypeScript guide documents — failed with `self is not defined`, and the
  CommonJS build was compiled and smoke-tested in CI but never uploaded.
  The published package is now **dual**: the `require` condition of the
  `exports` map resolves to the CommonJS half (`pkg/cjs/`), the `import`
  condition to the ESM half, and both expose the same API. CI and the
  release pipeline now pack the tarball, install it and require it *by
  name*, so the artefact that is published is the artefact that is tested.
- **The Generalized Born solvation returned `NaN` for a zero dielectric
  constant and aborted the process on a short charge vector.**
  `compute_gb_solvation` divided by `ε` without checking it, so
  `solvent_dielectric = 0` produced a `NaN` energy that looked like a
  number; a charge vector shorter than the atom list indexed past its end
  and killed the process (in Python, a `PanicException`, which inherits from
  `BaseException` and is therefore *not* caught by `except Exception:`).
  The same hole was in the ALPB sibling and in the Born-radius helpers.
  The solvation family now validates the geometry, the charge vector, the
  dielectric constants and the probe radius, and returns an error instead:
  `ValueError` in Python, `{"error": "..."}` in WASM, a non-zero exit in the
  CLI. `ε = 1.0` (the vacuum limit) and a zero probe radius (the plain van
  der Waals surface) remain valid.
- **The force-field energy and gradient are now bit-identical regardless of
  the rayon thread count.** The tree reduction made every UFF/MD/NEB result
  depend on the machine's core count — three different bit patterns were
  measured for the same molecule. This broke the crate's "same seed =
  reproducible output" promise on every path that matters. The same fix was
  applied to the GB electrostatic sum.
- **The UFF/MMFF94 energy accumulated in `f32`** while the rest of the library
  is `f64`. This is the routine the conformer minimiser uses.
- **Python registered `cpm_charges` twice**, so the beta copy silently won and
  computed at μ = 0.0 eV instead of the standard −4.44 eV.
- The swallowed errors are now reported: the faithful PM3 solver's real
  failure reason, the GPU→CPU fallbacks (`notes` on the result).
- The EHT vibrational analysis now reports `n_imaginary_modes` and a loud
  note. The EHT model is one-electron, so it has no repulsive wall and its
  Hessian curvature is negative: the numbers are unchanged but the failure is
  no longer silent.
- **The reaction-dynamics path repeated the same departure frame.** With a
  single product there is no second fragment to separate, so the departure
  phase returned the *same* geometry `n_departure_frames` times. It now
  contributes one frame and says so in the result notes.
- **The reaction-dynamics path could diverge into `NaN`.** The product-guided
  alignment of acetaldehyde + water put a water oxygen 0.14 Å away from the
  carbonyl carbon; the 1/r repulsion terms diverged, the analytic GFN0-xTB
  gradients overflowed and every remaining frame of the path was `NaN`. The
  fix adds a minimum inter-fragment separation (1.5 Å), a displacement cap
  and a descent-direction check in the constrained relaxation, and a refusal
  to run the CI-NEB on a singular geometry. The test that caught it,
  `tests/reaction_dynamics_study.rs`, had never run in CI — the packaging
  audit found 9 test targets that CI skipped.
- `align_coordinates` can no longer return NaN.
- The `unsafe` block in `stda_gpu.rs` had a factually wrong SAFETY comment and
  a silent empty-slice fallback; it now returns a `Result`.
- The GPU context and the compute pipeline are cached instead of being
  re-acquired and recompiled on every dispatch (27 ms → 6 ns).
- The Python integration test that passed a 3-atom geometry for a 9-atom
  molecule now passes the right one.

### Changed

- **The k-path and Fermi-surface chart builders now live in the core.** The
  Python `compute_kpoint_path_chart` / `compute_fermi_surface_chart`, the CLI
  `kpoint-path-chart` / `fermi-surface-chart` commands and the new WASM
  exports all call `kmesh::kpath_standard` / `kmesh::fermi_surface_from_bands`;
  the duplicated interpolation and contour-tracing code in every binding was
  removed.
- **The capabilities that were Rust-only are now available in Python, WASM
  and the CLI**: the open-shell SCF (`compute_uhf`, `compute_rohf`), the
  analytical PM3 and xTB gradients, the UFF and MMFF94 energy breakdowns, the
  topological polar surface area, the ensemble property prediction and the
  CIF import/export. The Python `__all__` grows from 247 to 265 names in the
  default build (327 in the release wheel, which now enables every alpha and
  beta module) and the CLI from 90 to 100 subcommands.
- `compute_esp` is now a deprecated alias of `compute_esp_grid`, the name of
  the Rust root function that returns the full grid.
- **Every public enum is now `#[non_exhaustive]`** (57 of 60; three
  deliberately closed sets are documented as such). Adding a variant is no
  longer a breaking change for downstream users. Downstream `match`
  expressions need a wildcard arm.
- **The canonical name of every capability is the Rust root name** in all four
  surfaces (Rust, Python, WASM, CLI). 93 deprecated Python aliases keep the
  old names working with a `DeprecationWarning`; WASM gained the missing
  names plus `*_seed64` variants so a 64-bit seed is no longer truncated.
- The non-idiomatic `get_` prefix is gone: `eht_support`,
  `system_capabilities`, `system_method_plan` (the old names remain as
  deprecated aliases).
- `compute_gfn1` and `compute_gfn2` now exist at the crate root.
- The Python classes no longer carry the `Py` suffix, the module exposes
  `__all__` (grouped by capability) and `__version__`, `parse` returns a
  `ParseResult` object instead of a plain dict, and the error taxonomy is
  explicit: `ValueError` (malformed input), `RuntimeError` (engine failure),
  `NotConvergedError` (non-convergence).
- The CLI uses one sentinel convention: an optional argument is simply
  absent.
- One error-message style across the three surfaces:
  `<What failed> '<actual value>'. Supported: 'a', 'b', 'c'.`
- The published crate carries an allow-list: 612 → 340 files, 20.6 → 8.0 MiB.
  The 27 developer binaries (`benchmark`, `test_eigen`, `check_grad`, …) are
  no longer installed on every user's `$PATH`.

### Added

- **Multi-linker framework assembly.** `Topology::from_net` converts an RCSR
  `NetData` entry into the symmetry-expanded graph used by assembly,
  `Topology::from_name` now accepts the data-driven nets (`acs`, `afw`, `sod`,
  `pts`, `xah`, `tfy`), and `Topology::edge_types` records one linker type per
  edge. `assemble_framework_with_linkers` places a different linker on each
  inequivalent edge; the single-linker `assemble_framework` remains as a
  wrapper. The expanded node/edge data is generated from the pormake `.cgd`
  files and the F67 gate checks it against pormake.
- **Core k-path and Fermi-surface helpers.** `kmesh::kpath_standard` samples
  the default Γ→X→M→Γ path with cumulative distances and vertex indices, and
  `kmesh::fermi_surface_from_bands` traces the `kz = 0` contour from real-space
  operator blocks. The WASM surface exposes
  `alpha_periodic_kpath_standard` and
  `alpha_periodic_fermi_surface_from_bands`.
- **Quantitative benchmarks for the remaining beta/alpha methods** in
  `tests/beta_benchmarks.rs`: CPM charges vs the closed-form two-site solution
  (1.3e-15), Riemannian optimization vs the exact Gram matrix of a tetrahedron
  (6.0e-9), SDR vs classical MDS (6.7e-16) and GSM barriers vs five analytic
  double wells (relative error 0), alongside the existing KPM, Nyström and
  MBH gates.
- Two zero-copy Python entry points (`embed_coords_into`,
  `compute_esp_grid_into`) that write into a caller-owned buffer.
- `try_*` variants and `Result` returns for every entry point that can fail,
  plus `MolecularSystem::{try_from_angstrom, try_from_flat_angstrom, validate}`,
  `scf::scf_loop::try_run_scf`, `spectroscopy::hessian::try_compute_hessian`,
  `spectroscopy::ir_intensities::try_compute_ir_intensities`,
  `spectroscopy::stda_uvvis::try_compute_stda`, `ani::neighbor::CellList::try_new`
  and `charges::gasteiger::GasteigerConfig::validate`.
- `forcefield::calculate_total_energy_f64` and the `f64` siblings of the
  energy functions.
- `gpu::context::GpuContext::shared()` — a process-wide cached context.
- An options-mapping form for the eight worst-over-argumented Python
  functions.

### Documentation

- **The entire public API is documented.** `cargo clippy -- -W missing_docs`
  reported 1038 undocumented items out of 2121 (49 % of the public API had no
  rustdoc at all); there are now 0. `cargo doc` is warning-free.
- Every numeric field now states its unit. This was the most valuable class
  of documentation to add: the library returns energies in Hartree in some
  result structs and eV in others, and a caller comparing `total_energy`
  across methods was off by 27×. The crate header now carries a per-method
  unit table.
- Doc comments that contradicted the code were corrected rather than papered
  over. The most consequential: `Hf3cResult.orbital_energies`, `homo_energy`,
  `lumo_energy` and `gap` were documented as eV, but the values are the raw
  Hartree-based SCF output with no conversion (water's 1a1 is −20.24 Ha, not
  −20.24 eV).
- The crate header states what the library is AND IS NOT, the coordinate
  convention, the element convention, the sign convention for gradients, a
  "how do I compute a property?" lookup table and the feature-flag table.
- The API pages (Rust, Python, TypeScript, CLI) were regenerated against the
  current code, with a units table on each.
- The materials-assembly page now documents the symmetry-expanded RCSR nets,
  the `edge_types`-based multi-linker assembly, and the real Rust/Python/CLI
  signatures. The complete supported-reference and TypeScript feature tables
  list the new multi-linker and periodic k-path/Fermi-surface entries.

### Packaging and CI

- The 27 developer binaries are out of the published crate.
- 13 MB of fixtures and the internal directories are out of the published
  crate.
- `rand` 0.8.5 → 0.8.8 and `crossbeam-epoch` 0.9.18 → 0.9.21 clear two live
  advisories; a `deny.toml` and a `cargo deny` CI job were added.
- The release workflow now verifies that the git tag matches the version of
  all four packages before publishing.
- The documented WASM build path now produces a working npm package (the
  `alpha` / `beta` subpaths and the `snippets` directory that
  wasm-bindgen-rayon needs are present in both the ESM and the CommonJS
  builds).
- All 9 integration test targets that CI never ran now run in CI (77 tests).
- The 20 hard-coded absolute paths in the scripts became an environment
  variable and a repo-relative default.
- The 52 heavy benchmark tests (10 k molecules, ChemBL 10 k, 1 k molecules ×
  14 capabilities) are `#[ignore]`d, so `cargo test --workspace --all-targets`
  no longer takes tens of minutes. Run them explicitly with
  `cargo test --release --test benchmarks -- --ignored`.
- 6 of the 9 parameter generators gained a `--check` mode so they cannot
  drift silently; CI checks all of them.
- `rust-version` (MSRV 1.87) is declared in all four manifests.

### Known limitations

- The 1.0 API returns a `Result` from the functions that can fail, but the
  error type is still `String`. A typed error enum is planned for 2.0
  because changing the public `Result` signatures is a breaking change.
- The ERI tensor is stored as a full M⁴ `Vec<f64>` (12.8 GB at M = 200).
  There is no packed, Cholesky or density-fitting path yet.
- The conformer minimiser rebuilds the force field on every MD step.
- The ETKDG retry loop is 10 × N attempts of an O(N²) search.
- The open items of the first audit (12 of 100) are listed in
  `reports/auditoria-resolucion.md`.

---

This section covers the audit-remediation work of 2026-10-21 (commits `c96ffa3`,
`6d3d954`, `28323f4`, `161c2e1`), driven by `reports/auditoria-total.md`.

### Detailed audit remediation (2026-10-21)

#### Added
- **Twenty-two Rust-only capabilities are now exposed on all three surfaces**
  (Rust, Python via PyO3, WASM/TypeScript, plus CLI subcommands):
  `embed_diverse`, `to_canonical_smiles`, `compute_ecfp_batch`,
  `compute_charges_configured`, `compute_ensemble_j_couplings`,
  `compute_ir_spectrum_broadened`, `compute_cisd`, `compute_dos_multimethod`,
  `compute_eht_gradient`, `optimize_geometry_eht`, `compute_aevs_tm`,
  `compute_whim`, `compute_rdf`, `compute_getaway`, `train_random_forest`,
  `train_gradient_boosting`, `assign_peaks`, `space_group_by_number` /
  `space_group_by_symbol`, `all_space_groups`, `build_periodic_molecule`,
  `detect_hapticity`, `optimize_framework` and `filter_diverse`.
  The CLI grew from 56 to 100 subcommands (34 of them feature-gated); see the
  regenerated `docs/api/cli.md`.
- **New validation gates.** Two batches, each with a committed fixture and a
  deterministic oracle script, so CI stays hermetic (no RDKit / PySCF / ASE /
  tblite needed at test time):
  - *F64–F76* — space groups vs spglib, CIF vs pymatgen, minimum-image
    distances vs ASE, powder XRD vs pymatgen, bounds matrix vs RDKit, ETKDG
    embedding vs RDKit ETKDGv3, the ETKDG-3D force field vs RDKit, BFGS vs
    RDKit `ForceField::minimize`, SMARTS matching vs `GetSubstructMatches`,
    and a SHA-256-pinned copy of the upstream CSD torsion table.
  - *F77–F85* — sTDA UV-Vis vs PySCF TD-SCF, ESP grid vs the PySCF density
    ESP, IR intensities vs PySCF, analytic gradients vs central finite
    differences for every method that exposes one, CIP R/S and E/Z labels vs
    RDKit, ML descriptors vs RDKit, canonical-SMILES round-trip and
    reordering invariance, ALPB electrostatics vs tblite GB solvation, and
    variable-cell stress / relaxed lattice vs ASE.
  - New test targets: `tests/validation_gates.rs`, `tests/rdkit_parity.rs`,
    `tests/rdkit_chem_parity.rs`, `tests/pyscf_parity.rs`,
    `tests/gradient_parity.rs`, `tests/mmff94_validation.rs`, all wired into
    the `validation` CI job.
- `tests/experimental_gates.rs` now aggregates the whole `tests/experimental/`
  suite (149 tests) with the **current** `alpha-*` / `beta-*` feature names.
  Previously these files were gated on the legacy `experimental-*` flags, which
  compiled them to empty modules and silently skipped every test inside them.
- `tests/integration/test_python_integration.py` is now a real pytest suite
  (77 tests) covering the whole public surface, and it runs in CI.
- The four WGSL shaders that existed but were missing from the registry
  (HCT Born radii, ANI radial AEV, MMFF94 non-bonded, sTDA matmul) are now
  registered, so the tier and budget checks cover them.

#### Changed
- **MMFF94 now uses the published MMFF94 parameter set** instead of placeholders.
  - New generated table module `src/forcefield/mmff94_params.rs` (atom-type
    equivalence levels, atom-type properties, bond / angle / torsion / van der
    Waals parameters, default bond rule, Badger's-rule parameters), produced by
    `scripts/validation/gen_mmff94_params.py` from the parameter tables embedded
    in RDKit's `Code/ForceField/MMFF/Params.cpp` (BSD-3-Clause; parameter data
    © Merck and Co., Inc., 1994–1996).
  - Parameter resolution follows the published rules: the MMFF.I atom-type
    equivalence hierarchy for angles, the five-stage MMFF.IV hierarchy for
    torsions, and the MMFF.V empirical default rules for bonds, angles and
    torsions when a typed parameter is absent.
  - New builder entry point `Mmff94Builder::build_typed(mol, variant)`; the
    existing `build` / `build_variant` entry points are unchanged and now
    derive their atom types from the same MMFF94 typer.
  - Fixed the cubic bend constant (was `-0.014 deg⁻¹`, now the MMFF94/RDKit
    value `-0.4 rad⁻¹ = -0.006981317 deg⁻¹`), switched the angle force-constant
    factor to the exact `143.9325 * (pi/180)^2`, added MMFF94's linear-centre
    cosine angle form, and replaced the placeholder 1-4-only van der Waals
    terms with the full MMFF94 pair set (all pairs separated by three or more
    bonds) and the MMFF94 combination rules for `R*_ij` / `eps_ij`.
  - Verified against an independent reference implementation
    (`scripts/validation/oracle_mmff94.py`), the fixture
    `tests/fixtures/validation/mmff94_rdkit.json` and the parity test
    `tests/mmff94_validation.rs` (agreement with the Rust implementation at
    the 1e-11 kcal/mol level for methane, ethanol, acetic acid, dimethyl ether,
    benzene, water and methylamine).
- **ANI forces are now the exact negative gradient of the energy.** The angular
  species-pair index aliased 6 of the 28 blocks and the gradient back-propagated
  only the radial block, so the reported forces were not `-dE/dR`. The index is
  now bijective (with a unit test) and the analytic gradient covers the full
  AEV; it is checked against central finite differences for water and methane
  (`src/ani/gradients.rs`).
- **AEV layout corrected.** The generalized Behler-Parrinello descriptor module
  filled the radial block with a zip of the shift and eta lists, so only 16 of
  the 256 `(eta, shift)` symmetry functions were ever computed and the remaining
  240 stayed zero. Both the implementation (`src/ml/symmetry_functions.rs`) and
  the F79 oracle now use the full grid, and the F79 gate compares the whole
  radial block instead of asserting that the tail is zero.
- **PyO3 upgraded 0.22 → 0.26** in `crates/python`. PyO3 0.22 reads CPython's
  tuple layout directly, which is incompatible with Python 3.14: tuple
  arguments segfaulted and tuple fields came back with a NULL element
  (`sci_form.esp(...).origin` returned `(-3.757, -3.469, <NULL>)`). The
  `abi3-py38` feature keeps a single cp38-abi3 wheel working on 3.8–3.14.
- `sci_form.alpha` and `sci_form.beta` are now reachable as attributes of the
  top-level module (they were registered under their dotted `__name__`).
- The TypeScript declarations for the alpha/beta subpaths were regenerated from
  the Rust signatures (optional-vs-required config arguments, `apply_thermostat`
  / `apply_steering_force` / `set_temperature`, `num_atoms`) and 29 undeclared
  functions were added.
- `GpuContext::try_create` now wraps adapter/device acquisition in
  `catch_unwind`; wgpu-hal's GLES/EGL backend panics with `BadDisplay` on
  headless hosts, which broke the documented CPU-fallback contract.
- The MMFF94 WGSL kernel no longer accumulates float energies with `atomicAdd`
  on `bitcast<u32>` (integer addition of float bit patterns); it reduces in
  workgroup memory and performs signed fixed-point atomic adds, with a test
  that reproduces the kernel arithmetic on the host.
- `atomic_mass_amu` now carries the IUPAC table for Z = 1..118 instead of
  returning 12.0 amu for ~100 elements.

#### Removed
- **~63 Python / WASM / CLI binding functions that returned hard-coded data
  have been removed or wired to the real core implementations.** The removed
  stubs included `total_dos: vec![1.0; 100]`, a hard-coded
  `homo_energy_ev: -5.0`, a constant `"Generated: 2026-03-28"` timestamp, an
  empty V2000 molfile and always-true validators.
- The absolute symlinks `pkg` and `pkg-node` (which pointed into another user's
  home directory) are gone; the WASM output now lives under
  `crates/wasm/pkg{,node}`.

#### Fixed
The eight defects the new gates proved, all fixed in the same cycle:

1. **Stereo** — the R/S assignment was inverted. `assign_rs` mapped a positive
   signed volume to "R", but the CIP rule is clockwise-1-2-3-with-4-away = R
   for a *negative* signed volume. All 11 stereocentres in the corpus had the
   wrong label (E/Z was already correct).
2. **SCF gradients** — `nuclear_repulsion_gradient` returned the force
   (`-dV/dR`) although it is named and documented as the gradient.
3. **Spectroscopy / Hessian** — the frequency conversion carried a spurious
   extra `2π`, so the frequencies were 6.283× too small.
4. **Spectroscopy / IR intensities** — the conversion factor was applied to a
   derivative in e·a₀/Bohr instead of Debye/(Å·√amu), so the intensities were
   23.07× too small.
5. **Solvation / ALPB Born radii** — the descreening integral was garbled and
   produced negative ψ, hence unphysical Born radii (0.67 Å for H) and ALPB
   electrostatics off by up to 20 kcal/mol. Replaced with the OBC/ACE
   formulation; radii are now 1.6–2.7 Å.
6. **ML descriptors** — the atomic-weight fallback was `z * 1.5` (27.96 Da error
   on metal-containing molecules) and `balaban_j` was the placeholder
   `|ln μ| + 1` (6.66 off). Both are now exact against RDKit.
7. **ML ensemble** — the Ertl TPSA was a rough approximation (280 Å² off) and
   is now the real RDKit-compatible fragment table (5.7e-14 Å²).
8. **Canonical SMILES** — ring digits were emitted in `HashMap` order, so the
   output was non-deterministic; the numbering was not idempotent; formal
   charges and explicit H on bracket atoms were dropped; aromatic-aromatic
   single bonds were not written out. Rewritten with
   individualisation-refinement.

The RDKit-parity batch (F71–F76) proved three more:

- **Distance geometry / bounds** — the "squish" slack was applied when *either*
  end of a bond was flagged, but RDKit's `set12Bounds` requires *both* ends.
  Every bond touching a squished heteroatom was therefore looser than RDKit's
  (thiophene C–C: ±0.21 Å instead of ±0.01 Å).
- **SMARTS torsion table** — twelve macrocycle patterns carried `r{9-}` where
  upstream has a bare `r`, changing the semantics from "in any ring" to "in a
  ring of size ≥ 9".
- **SMARTS torsion table** — the eight sci-form-specific patterns were mixed
  into the upstream table, so it was not a faithful transcription. They now
  live in `TORSION_PREFS_EXTENSIONS`; the matching order is unchanged.

And the F64–F69 batch proved two more:

- The X-ray scattering-factor table was read one element off, so every element
  used its neighbour's coefficients (NaCl intensities were off by up to 14 %).
- The minimum-image search was only valid for orthogonal cells.

#### Known limitations
- The F80 gradient gate documents three analytic gradients that are still wrong
  and are marked `#[ignore]` with their measured numbers: PM3 (18.4 eV/Å),
  GFN0-xTB (209–263 eV/Å) and KS-DFT (0.202 Ha/Bohr, missing the adiabatic
  XC response term).
- The official ANI-2x `.netFormat` weight files still cannot be loaded — the
  loader in `src/ani/weights.rs` only reads sci-form's own binary container.
- The browser WebGPU harness (`crates/wasm/examples/vite-webgpu-hybrid/`) is
  not in the repository, so `npm run webgpu:test` cannot produce a report; the
  WebGPU path is covered by the Rust-side GPU-vs-CPU gates instead. See
  `docs/algorithm/webgpu-validation.md`.

#### Documentation
- `docs/api/cli.md` regenerated from the real `Commands` enum (90 subcommands,
  grouped by area, feature-gated ones marked); `docs/guide/cli.md` rewritten as
  a task-oriented tour.
- `docs/api/python.md`, `docs/api/typescript.md`, `docs/guide/python.md`,
  `docs/guide/rust.md`, `docs/guide/typescript.md`, `docs/guide/cli.md` and
  `docs/index.md` corrected (removed non-existent functions, fixed every code
  block, synchronised the version string to 0.15.2).
- `docs/algorithm/machine-learning-potentials.md` corrected: the non-existent
  `load_ani2x_weights` is replaced by the real `ani::weights::load_weights`, the
  two distinct AEV layouts are documented separately, and the inability to load
  the official ANI-2x weights is stated explicitly.
- `docs/algorithm/webgpu-validation.md` rewritten to state that the browser
  harness does not exist and what is actually runnable.
- `docs/guide/reaction-viewer.md` corrected: the example uses hard-coded
  analytic 1-D potentials and never imports `sci-form-wasm`; the live
  sci-form backends are not yet wired in.
- `docs/benchmarks.md` restructured into three explicit evidence tiers and the
  stale figures (embed success, GDB-20 corpus size and seed count) corrected
  against the current test output.
- `docs/experimental/*.md` test commands fixed: the suite lives in the
  `experimental_gates` target, not `regression` (where the filters matched no
  test and passed vacuously).
- The 12 orphan documentation pages are now reachable from the VitePress
  sidebar.

## v0.15.0 - 2026-04-06

### Added
- New `alpha-reaction-dynamics` pipeline for full 3D reaction-path work.
- Climbing-image NEB, IDPP path initialization, IRC tracing, constrained optimization,
  electrostatic steering, orbital-guided approach, per-frame properties, and angular sampling.
- Dedicated alpha reaction-dynamics regression tests, including checks for:
  - no hardcoded X-axis forcing,
  - SMIRKS-driven reactive-site identification,
  - reactive atom-pair extraction from bond changes,
  - end-to-end SMIRKS integration with the reaction pipeline.

### Documentation
- Added `ROADMAP_REACTION_DYNAMICS_3D.md` to track the 3D reaction-dynamics plan and gaps.

## v0.14.x - 2026-04-04 to 2026-04-05

### Added in v0.14.0
- SMIRKS parsing and reaction-transform support in Rust.
- Python bindings for SMIRKS (`crates/python/src/smirks.rs`).
- WASM reaction bindings (`crates/wasm/src/reaction.rs`).
- Reaction documentation and validation assets:
  - `SMIRKS_EXAMPLES.md`
  - `SMIRKS_REVIEW.md`
  - `SMIRKS_TESTING.md`
  - `docs/guide/reaction-backends.md`
  - `tests/test_smirks_reactions.rs`
  - `tests/integration/test_smirks_reactions.py`
  - `scripts/compare_smirks_reactions.py`
  - `scripts/compare_reaction_layers.py`
- CIF import/export support for crystallographic workflows.
- AO→MO integral transform and open-shell HF (`UHF`) infrastructure.
- xTB upgrades including Broyden SCC support, D4 dispersion data, and expanded GFN1/GFN2 internals.
- Additional spectroscopy and GPU groundwork, including Hessian and sTDA support modules.

### Follow-up work in v0.14.1-v0.14.4
- Stabilization and polish across SMIRKS, reaction tests, documentation, and CLI/WASM/Python surfaces.
- Maintenance updates in xTB, PM3, CIF, NMR, population analysis, UHF, and reaction benchmark coverage.

## v0.13.0 - 2026-03-28

### Added
- Architecture and export-matrix documentation for alpha/beta surfaces.
- New alpha modules for electrochemical double layers, kinetics, periodic linear algebra,
  and render-bridge/GPU integration.
- Expanded CLI, Python, and WASM experimental exports for those alpha modules.

### Updated
- GPU plumbing for MMFF94, PM3, and xTB-related experimental paths.
- Experimental documentation and regression tests for alpha modules.

## v0.12.x - 2026-03-28

### Added and improved across v0.12.0-v0.12.4
- Broader alpha/beta binding coverage in Python and WASM.
- Continued work on:
  - CGA conformer and motor utilities,
  - GSM string/saddle routines,
  - SDR embedding/projection paths,
  - CPM, MBH, RandNLA, and Riemannian beta modules,
  - DOS and PM3 numerical improvements.

### Release-line focus
- Consolidation and stabilization of experimental alpha/beta APIs before the larger
  documentation and reaction-system expansion that followed in `v0.13.0+`.
