# sci-form API Surface Map

This is a compact map of the major public surfaces. For exact signatures and current field names, read the source before editing.

---

## Production APIs (always available)

### Geometry and topology

- `embed` → ConformerResult
- `embed_batch` → Vec<ConformerResult>
- `parse` → {num_atoms, num_bonds, ...}
- `compute_rmsd` → f64

### Charges, surface, population, and electrostatics

- `compute_charges` → ChargeResult
- `compute_charges_configured` → ChargeResult
- `compute_sasa` → SasaResult
- `compute_population` → PopulationResult
- `compute_dipole` → DipoleResult
- `compute_esp` → EspGrid (JSON or typed grid)
- `compute_dos` → DosResult

### Force fields and alignment

- `compute_uff_energy` → f64 (kcal/mol)
- `compute_mmff94_energy` → f64 (kcal/mol)

### Semi-empirical and tight-binding

- `compute_pm3` → Pm3Result (energies, charges, orbital_energies)
- `compute_xtb` → XtbResult (GFN0)
- `solve_gfn1` → Gfn1Result (GFN1 + D3)
- `solve_gfn2` → Gfn2Result (GFN2 + D4 + XB)

### EHT and related electronic structure

- `eht_calculate` → EhtResult (Hamiltonian, overlap, orbital energies)
- `eht_orbital_mesh` → OrbitalMesh (vertices, normals, isosurface)
- `eht_orbital_grid_typed` → Float32Array (zero-copy from WASM)
- `compute_esp_grid_typed` → Float64Array (zero-copy from WASM)
- `compute_esp_grid_info` → {origin, spacing, dims}
- `compute_dos_multimethod` → DosResult (method-agnostic)
- `compute_eht_gradient` → EhtGradient (forces)
- `optimize_geometry_eht` → EhtOptResult (optimized positions, energy history)
- `compute_band_structure` → BandStructure (k-points, bands, Fermi level)

### Ab initio and excited states

- `solve_hf3c` → Hf3cResult (energies, D3/gCP/SRB corrections, optional CISD)
- `compute_cisd` → CisdResult (excited states, oscillator strengths)

### ANI potentials

- `compute_ani` → AniResult (energy, forces, atomic energies from ML net)
- `compute_aevs_tm` → Vec<Vec<f64>> (ANI-TM atomic environment vectors, 24 elements)

### ML descriptors and models

- `compute_ml_descriptors` → MolecularDescriptors (MW, Wiener, Balaban, FSP3, ...)
- `predict_ml_properties` → MlPropertyResult (LogP, MR, solubility, Lipinski, druglikeness)
- `compute_whim` → WhimDescriptors (eigenvalues, spread, anisotropy with mass/volume/EN weighting)
- `compute_rdf` → RdfDescriptors (radial distribution functions)
- `compute_getaway` → GetawayDescriptors (GETAWAY leverage, autocorrelation, information content)
- `train_random_forest` → RandomForest (regressor with variance predictions)
- `train_gradient_boosting` → GradientBoosting (regressor)

### Stereochemistry and reactivity

- `analyze_stereo` → StereoAnalysis (R/S priorities, E/Z, helical, atropisomeric)
- `compute_frontier_descriptors` → FrontierDescriptors (HOMO/LUMO atom contributions, dual descriptor)
- `compute_fukui_descriptors` → FukuiDescriptors (condensed Fukui indices for reactivity)
- `compute_reactivity_ranking` → ReactivityRanking (nucleophilic, electrophilic, radical attack sites)
- `compute_empirical_pka` → EmpiricalPkaResult (predicted acidic/basic sites)

### Spectroscopy

- `compute_stda_uvvis` → StdaUvVisSpectrum (broadened UV-Vis, peaks)
- `compute_vibrational_analysis` → VibrationalAnalysis (frequencies, modes, ZPVE, thermochemistry)
- `compute_ir_spectrum` → IrSpectrum (broadened IR, peaks with functional group assignment)
- `assign_peaks` → PeakAssignment (C=O, O-H, N-H, etc.)
- `compute_ensemble_j_couplings` → Vec<JCoupling> (Karplus ensemble-averaged)
- `predict_nmr_shifts` → NmrShiftResult (predicted ¹H, ¹³C, ¹⁹F, etc.)
- `predict_nmr_couplings` → Vec<JCoupling> (predicted ²J, ³J, etc.)
- `compute_nmr_spectrum` → NmrSpectrum (broadened NMR axis, peaks)

### Solvation

- `compute_nonpolar_solvation` → NonPolarSolvation (SASA-based energy, kcal/mol)
- `compute_gb_solvation` → GbSolvation (electrostatic + non-polar, kcal/mol)

### Rings, fingerprints, and clustering

- `compute_sssr` → SssrResult (smallest set of smallest rings, aromatic flags)
- `compute_ecfp` → ECFPFingerprint (Morgan fingerprint, sparse or dense)
- `compute_tanimoto` → f64 (0.0 to 1.0 similarity)
- `butina_cluster` → ClusterResult (centroid indices, cluster assignments)
- `compute_rmsd_matrix` → Vec<Vec<f64>> (all-pairs RMSD)
- `filter_diverse` → Vec<usize> (Butina diversity filtering)

### Materials and periodic systems

- `create_unit_cell` → UnitCell (a, b, c, α, β, γ)
- `assemble_framework` → CrystalStructure (MOF assembly from SBUs)
- `space_group_by_number` → Option<SpaceGroup> (1–230 ITC groups)
- `space_group_by_symbol` → Option<SpaceGroup> ("Fm-3m", etc.)
- `build_periodic_molecule` → PeriodicMolecule (fractional coords with PBC)
- `detect_hapticity` → HapticityAnalysis (metal-ligand interactions, metallocene flag)
- `optimize_framework` → FrameworkOptResult (BFGS/steepest descent with PBC)

### Reactions and transport

- `parse_smirks` → SmirksTransform (atom-mapped pattern)
- `apply_smirks` → SmirksResult (products from reaction template)
- `pack_conformers` → Arrow RecordBatch (columnar format for batch transport)
- `split_worker_tasks` → Vec<{smiles, seed}> (work split across N workers)
- `estimate_workers` → usize (recommended thread count)

---

## Alpha Module APIs (`alpha-*` feature gates)

### A1 — DFT (`alpha-dft`)

- Rust: `dft::ks_fock::solve_ks_dft(elements, positions, config)` → `DftResult`
- Python: `sci_form.alpha.dft_calculate(elements, coords, method="pbe")`
- WASM: `alpha_compute_dft(elements_json, coords_json, method)`
  - Methods: SVWN (LDA), PBE (GGA)
  - Returns: `{ energy, xc_energy, nuclear_repulsion, homo_energy, lumo_energy, gap, converged, scf_iterations, mulliken_charges, orbital_energies }`

### A2 — ReaxFF (`alpha-reaxff`)

- Rust: `forcefield::reaxff::compute_reaxff_gradient(positions_flat, elements, params)` → `(energy_kcal_mol, gradient_flat)`
- Python: `sci_form.alpha.reaxff_gradient(elements, coords)` and `reaxff_energy(elements, coords)`
- WASM: `alpha_compute_reaxff_gradient(...)`, `alpha_compute_reaxff_energy(...)`, `alpha_compute_eem_charges(...)`
  - **Parallelized** (with `features = ["parallel"]`): each coordinate displacement independent
  - Binding helpers expose energy breakdown `{ bonded, coulomb, van_der_waals, total }` and EEM charges

### A3 — MLFF (`alpha-mlff`)

- Rust: `mlff::compute_aevs(elements, positions, params)` → `Vec<Aev>`
- Rust: `mlff::compute_mlff(elements, positions, config)` → `MlffResult`
- Python: `alpha_compute_aevs(elements, coords, radial_cutoff=6.0, angular_cutoff=3.5)` and `alpha_compute_mlff(...)`
- WASM: `alpha_compute_aevs(elements_json, coords_json, params_json)`, `alpha_compute_mlff(...)`
  - **Parallelized** (with `features = ["parallel"]`): outer atom loop via rayon
  - Rust `Aev` stores `radial` and `angular`; Python/WASM flatten to arrays plus metadata

### A4 — Obara-Saika (`alpha-obara-saika`)

- Rust: `scf::obara_saika::boys_function`, `eri_ssss`, `schwarz_bound`
- Python/WASM: `alpha_boys_function`, `alpha_eri_ssss`, `alpha_schwarz_bound`

### A5 — CGA (`alpha-cga`)

- Rust: `alpha::cga::dihedral_motor`, `alpha::cga::apply_motor_to_subtree`, `alpha::cga::refine_torsion_cga`
- Python: `rotate_dihedral_cga(...)`, `refine_torsion_cga(...)`
- WASM: `alpha_rotate_dihedral_cga(...)`, `alpha_refine_torsion_cga(...)`

### A6 — GSM (`alpha-gsm`)

- Rust: `alpha::gsm::interpolate_node(reactant, product, t)` and `alpha::gsm::find_transition_state(...)`
- Python: `gsm_interpolate(...)`, `gsm_find_ts(smiles, reactant_coords, product_coords, n_nodes=9)`
- WASM: `alpha_gsm_interpolate(...)`, `alpha_gsm_find_ts(...)`

### A7 — SDR (`alpha-sdr`)

- `alpha::sdr::sdr_embed(distance_pairs, n_atoms, config)` → SdrResult (coordinates from SDP, convergence info)
  - Python binding name: `sdr_embed_py(...)`
  - WASM binding name: `alpha_sdr_embed(...)`
  - Input: `Vec<(i, j, d_ij)>` distance constraint triples

### A8 & A9 — Dynamics-Live & IMD (`alpha-dynamics-live`, `alpha-imd`)

- WASM exposes `LiveSimulation` from `crates/wasm/src/dynamics_live.rs`
- Core Rust implementation lives under `dynamics_live/` with IMD steering support
- Use the alpha docs for the full interactive API because this surface is class-based rather than a single free function

---

## Beta Module APIs (`beta-*` feature gates)

### B1 — KPM (`beta-kpm`)

- Rust: `beta::kpm::compute_kpm_dos(h, config, e_min, e_max, n_points)` → `KpmDosResult`
- Python: `kpm_dos(elements, coords, order=100, temperature=0.0)`
- WASM: `beta_compute_kpm_dos(elements_json, coords_json, config_json)`
  - **Parallelized** (with `features = ["parallel"]`): random probe vectors via rayon
  - O(N) scaling for sparse matrices, stochastic trace estimation
- Rust: `beta::kpm::compute_kpm_mulliken(h, s, n_electrons, nuclear_charges, config)`
- Python/WASM build EHT internally before running KPM helper functions

### B2 — MBH (`beta-mbh`)

- `beta::mbh::compute_mbh_frequencies(elements, positions, rings, energy_fn, config)` → MbhResult
  - (frequencies_cm⁻¹, n_blocks, n_flexible, n_dof_reduced, speedup factor)
  - Reduced-dimension Hessian for large molecules

### B3 — RandNLA (`beta-randnla`)

- Rust: `beta::rand_nla::solve_eht_randnla(h, overlap, config)` → `(eigenvalues, eigenvectors, info)`
- Python: `eht_randnla(elements, coords, sketch_size=None, max_error=0.001)`
- WASM: `beta_solve_eht_randnla(elements_json, coords_json, config_json)`
  - O(N k²) cost where k = sketch dimension ≈ √N
  - Includes residual error and fallback flag

### B4 — Riemannian (`beta-riemannian`)

- `beta::riemannian::psd_distance(A, B)` → f64 (affine-invariant geodesic metric)
- `beta::riemannian::psd_projection(M)` → Matrix (project symmetric M onto PSD cone)
- `beta::riemannian::riemannian_gradient(G, X)` → Matrix (map Euclidean gradient to Riemannian)

### B5 — CPM (`beta-cpm`)

- Rust: `beta::cpm::compute_cpm_charges(elements, positions, config)` → `CpmResult`
- Python: `cpm_charges(elements, coords, potential=0.0)`
- WASM: `beta_compute_cpm_charges(elements_json, coords_json, potential)`
  - (charges, total_charge, energy at fixed electrode potential)
  - Constant-potential method for electrochemistry

---

## Surface-specific notes

### Rust

- Returns native structs and result types: `Result<T, String>` for errors.
- Coordinates always flat: `Vec<f64>` or slices `&[f64]`.
- Parallel execution activated by `features = ["parallel"]` in `Cargo.toml` used at build time.
- Alpha algorithms are split between feature-gated top-level modules (`dft`, `mlff`, `dynamics_live`, `forcefield::reaxff`) and `alpha::*` research namespaces (`cga`, `gsm`, `sdr`).

### Python (PyO3 0.22)

- Returns PyO3 wrapper objects with attribute access (`.property`).
- Submodule imports: `from sci_form.alpha import ...`, `from sci_form.beta import ...`.
- Errors raised as `PyRuntimeError` for `Result<_, String>` mappings.

### WASM (wasm-bindgen)

- Structured inputs/outputs use JSON strings: `JSON.stringify(...)`, `JSON.parse(...)`.
- Performance-critical paths use typed arrays: `Float64Array`, `Float32Array` (zero-copy from WASM memory).
- Subpath imports: `import {...} from 'sci-form-wasm/alpha'`, `import {...} from 'sci-form-wasm/beta'`.
- All functions have generated TypeScript declarations (`.d.ts` files).

### CLI

- Subcommands mirror Rust APIs: `sci-form embed CCO --seed 42 --format xyz`.
- Alpha/beta commands available when library compiled with features.
- JSON output by default; format flags (`--format xyz|sdf|json`) where applicable.

---

## Feature gate reference

| Feature | What it enables | Modules |
|---------|-----------------|---------|
| `parallel` | rayon multithreading in compute-intensive loops | AEV, ReaxFF gradients, KPM probes |
| `experimental-gpu` | wgpu-based GPU acceleration (CPU fallback) | MMFF94 nonbonded, experimental AEV offload |
| `alpha-dft` | Kohn-Sham DFT solver (requires Obara-Saika) | A1 |
| `alpha-reaxff` | ReaxFF force field + EEM charges | A2 |
| `alpha-mlff` | Neural FF with AEVs | A3 |
| `alpha-obara-saika` | ERI engine, Boys function, Schwarz bounds | A4 |
| `alpha-cga` | Conformal GA motor rotations | A5 |
| `alpha-gsm` | Growing String Method | A6 |
| `alpha-sdr` | SDP-based embedding | A7 |
| `alpha-dynamics-live` | Interactive MD (required by IMD) | A8 |
| `alpha-imd` | IMD wire protocol (includes A8) | A9 |
| `beta-kpm` | Chebyshev polynomial DOS | B1 |
| `beta-mbh` | Mobile Block Hessian | B2 |
| `beta-randnla` | Randomized sparse EHT | B3 |
| `beta-riemannian` | Riemannian PSD geometry | B4 |
| `beta-cpm` | Constant Potential Method | B5 |
