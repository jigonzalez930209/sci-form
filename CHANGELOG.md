# Changelog

This file captures the project-level highlights from `v0.12.0` onward.
GitHub Releases already exist for each tag, but most release bodies are auto-generated;
the summary below records the substantive work that landed in each release line.
Future releases should append new sections here so this file remains the maintained,
human-written changelog for the repository.

## Unreleased

### Changed
- MMFF94 now uses the **published MMFF94 parameter set** instead of placeholders.
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
- Added an independent MMFF94 reference implementation
  (`scripts/validation/oracle_mmff94.py`), the fixture
  `tests/fixtures/validation/mmff94_rdkit.json` and the parity test
  `tests/mmff94_validation.rs` (agreement with the Rust implementation at the
  1e-11 kcal/mol level for methane, ethanol, acetic acid, dimethyl ether,
  benzene, water and methylamine).

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
