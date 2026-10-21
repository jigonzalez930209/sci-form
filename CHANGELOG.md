# Changelog

This file captures the project-level highlights from `v0.12.0` onward.
GitHub Releases already exist for each tag, but most release bodies are auto-generated;
the summary below records the substantive work that landed in each release line.
Future releases should append new sections here so this file remains the maintained,
human-written changelog for the repository.

## Unreleased

This section covers the audit-remediation work of 2026-10-21 (commits `c96ffa3`,
`6d3d954`, `28323f4`, `161c2e1`), driven by `reports/auditoria-total.md`.

### Added
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
  The CLI grew from 56 to 90 subcommands (34 of them feature-gated); see the
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

### Changed
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

### Removed
- **~63 Python / WASM / CLI binding functions that returned hard-coded data
  have been removed or wired to the real core implementations.** The removed
  stubs included `total_dos: vec![1.0; 100]`, a hard-coded
  `homo_energy_ev: -5.0`, a constant `"Generated: 2026-03-28"` timestamp, an
  empty V2000 molfile and always-true validators.
- The absolute symlinks `pkg` and `pkg-node` (which pointed into another user's
  home directory) are gone; the WASM output now lives under
  `crates/wasm/pkg{,node}`.

### Fixed
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

### Known limitations
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

### Documentation
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
