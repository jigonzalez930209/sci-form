# sci-form API Surface Map

This is a compact map of the major public surfaces. For exact signatures and current field names, read the source before editing.

## Geometry and topology

- `embed`
- `embed_batch`
- `parse`
- `compute_rmsd`

## Charges, surface, population, and electrostatics

- `compute_charges`
- `compute_charges_configured`
- `compute_sasa`
- `compute_population`
- `compute_dipole`
- `compute_esp`
- `compute_dos`

## Force fields and alignment

- `compute_uff_energy`
- `compute_mmff94_energy`

## Semi-empirical and tight-binding

- `compute_pm3`
- `compute_xtb`
- `solve_gfn1`
- `solve_gfn2`

## EHT and related electronic structure

- `eht_calculate`
- `eht_orbital_mesh`
- `eht_orbital_grid_typed`
- `compute_esp_grid_typed`
- `compute_esp_grid_info`
- `compute_dos_multimethod`
- `compute_eht_gradient`
- `optimize_geometry_eht`
- `compute_band_structure`

## Ab initio and excited states

- `solve_hf3c`
- `compute_cisd`

## ANI potentials

- `compute_ani`
- `compute_aevs_tm`

## ML descriptors and models

- `compute_ml_descriptors`
- `predict_ml_properties`
- `compute_whim`
- `compute_rdf`
- `compute_getaway`
- `train_random_forest`
- `train_gradient_boosting`

## Stereochemistry and reactivity

- `analyze_stereo`
- `compute_frontier_descriptors`
- `compute_fukui_descriptors`
- `compute_reactivity_ranking`
- `compute_empirical_pka`

## Spectroscopy

- `compute_stda_uvvis`
- `compute_vibrational_analysis`
- `compute_ir_spectrum`
- `assign_peaks`
- `compute_ensemble_j_couplings`
- `predict_nmr_shifts`
- `predict_nmr_couplings`
- `compute_nmr_spectrum`

## Solvation

- `compute_nonpolar_solvation`
- `compute_gb_solvation`

## Rings, fingerprints, and clustering

- `compute_sssr`
- `compute_ecfp`
- `compute_tanimoto`
- `butina_cluster`
- `compute_rmsd_matrix`
- `filter_diverse`

## Materials and periodic systems

- `create_unit_cell`
- `assemble_framework`
- `space_group_by_number`
- `space_group_by_symbol`
- `build_periodic_molecule`
- `detect_hapticity`
- `optimize_framework`

## Reactions and transport

- `parse_smirks`
- `apply_smirks`
- `pack_conformers`
- `split_worker_tasks`
- `estimate_workers`

## Surface-specific notes

- Rust returns native structs and result types.
- Python returns PyO3 wrapper objects with attribute access.
- WASM uses JSON strings for structured inputs and outputs, with typed-array fast paths where exposed.
- CLI mirrors the same capabilities through subcommands and format flags.
