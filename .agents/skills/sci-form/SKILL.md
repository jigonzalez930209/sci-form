---
name: sci-form
description: "High-performance library for 3D molecular conformer generation and quantum chemistry calculations (ETKDG, EHT, UFF, PM3, xTB, spectroscopy). Full API reference with examples in Rust, WASM/JS, Python, and CLI."
---

# `sci-form` — Complete API Reference for AI Agents

**Version:** 0.7.0 | **Repo:** [github.com/jigonzalez930209/sci-form](https://github.com/jigonzalez930209/sci-form)

`sci-form` provides high-performance molecular informatics across four surfaces: **Rust crate**, **WebAssembly (JS/TS)**, **Python (PyO3)**, and **CLI**. All heavy computation is parallelized with **Rayon**.

## Data Conventions

| Parameter | Format |
|---|---|
| `elements` | `Vec<u8>` / JSON `[8, 1, 1]` — atomic numbers |
| `coords` | Flat `Vec<f64>` / JSON `[x0,y0,z0, x1,y1,z1,…]` in Å |
| `smiles` | Standard SMILES string |
| WASM output | JSON string or zero-copy `Float32Array`/`Float64Array` |
| Python output | Typed dataclass-like `#[pyclass]` objects |

## 🔀 Parallelism

| Platform | Mechanism | How to enable |
|---|---|---|
| Rust | Rayon (default feature `parallel`) | `features = ["parallel"]` in Cargo.toml |
| CLI | Rayon auto | `--threads 0` (all cores) |
| Python | Rayon (GIL released) | `num_threads=0` in batch calls |
| WASM | `wasm-bindgen-rayon` + Web Workers | `await initThreadPool(navigator.hardwareConcurrency)` before first call |

---

## API Groups

| # | Group | Key functions |
|---|---|---|
| 1 | Geometry / Embedding | `embed`, `embed_batch`, `parse`, `rmsd` |
| 2 | Force Fields | `uff_energy`, `mmff94_energy` |
| 3 | EHT / Quantum | `eht_calculate`, `eht_orbital_mesh`, `dos`, `compare_methods` |
| 4 | Charges & Surface | `charges`, `sasa`, `dipole`, `esp` |
| 5 | Population / Bond Orders | `population`, `bond_orders` |
| 6 | Reactivity | `fukui_descriptors`, `frontier_descriptors`, `reactivity_ranking` |
| 7 | Spectroscopy UV-Vis | `uvvis_spectrum`, `stda_uvvis` |
| 8 | Spectroscopy IR | `vibrational_analysis`, `ir_spectrum` |
| 9 | Spectroscopy NMR | `nmr_shifts`, `nmr_couplings`, `nmr_spectrum`, `hose_codes` |
| 10 | Semi-empirical QM | `pm3_calculate`, `xtb_calculate` |
| 11 | ML Properties | `ml_descriptors`, `ml_predict` |
| 12 | Topology & Graph | `graph_features`, `topology_analysis` |
| 13 | Materials / Crystal | `unit_cell`, `assemble_framework` |
| 14 | Transport / Batch | `pack_conformers`, `split_worker_tasks`, `estimate_workers` |
| 15 | System Planning | `system_capabilities`, `system_method_plan`, `compare_methods` |
| 16 | Stereochemistry | `analyze_stereo`, `stereo_analysis`, `stereo` CLI |
| 17 | Solvation | `compute_nonpolar_solvation`, `compute_gb_solvation` |
| 18 | Rings & Fingerprints | `compute_sssr`, `compute_ecfp`, `compute_tanimoto` |
| 19 | Clustering | `butina_cluster`, `compute_rmsd_matrix`, `filter_diverse` |

---

## 1. Geometry / Embedding

### `embed(smiles, seed)` → ConformerResult
Generates a single 3D conformer using ETKDG.

- `smiles: &str` — SMILES string
- `seed: u64/u32` — RNG seed (default 42)

**Returns:** coords (flat Vec<f64>), elements, bonds, num_atoms, time_ms, error

### `embed_batch(smiles_list, config)` → Vec\<ConformerResult\>
Parallel batch embedding. Uses Rayon internally.

- `smiles_list: &[&str]` — list of SMILES
- `config: ConformerConfig { seed: u64, num_threads: usize }` — `num_threads=0` = auto

### `parse(smiles)` → Molecule
Topology-only parse. No 3D coordinates generated.

### `rmsd(coords, reference)` → AlignmentResult
Kabsch optimal alignment + RMSD.

- `coords/reference: Vec<f64>` — flat coordinate arrays (same length)

**Returns:** `rmsd: f64`, `aligned_coords: Vec<f64>`

---

## 2. Force Fields

### `compute_uff_energy(smiles, coords)` → f64
UFF force field energy in kcal/mol.

### `compute_uff_energy_with_aromatic_heuristics(smiles, coords)` → UffHeuristicEnergy
UFF energy with aromaticity correction. Returns `raw_energy_kcal_mol`, `aromatic_stabilization_kcal_mol`, `corrected_energy_kcal_mol`, `aromatic_bond_count`, `notes`.

### `compute_mmff94_energy(smiles, coords)` → f64
MMFF94 force field energy in kcal/mol.

---

## 3. EHT / Quantum

### `eht_calculate(elements, coords, k)` → EhtResult
Extended Hückel Theory secular equation solver.

- `k: f64` — Wolfsberg-Helmholtz constant (default 1.75, 0 = default)

**Returns:** `energies`, `n_electrons`, `homo_index`, `lumo_index`, `homo_energy`, `lumo_energy`, `gap`, `support_level`, `warnings`, `coefficients`

### `eht_orbital_mesh(elements, coords, mo_index, spacing, isovalue)` → IsoMesh
3D orbital isosurface via Marching Cubes.

- `mo_index: usize` — MO index (0-based)
- `spacing: f64` — grid spacing in Å (e.g. 0.2)
- `isovalue: f32` — isosurface cutoff (e.g. 0.02)

**Returns:** `vertices`, `normals`, `indices`, `num_triangles`, `isovalue`

### WASM-only Zero-Copy Grids
- `eht_orbital_grid_typed(elements, coords, mo_index, spacing)` → `Float32Array`
- `eht_orbital_grid_from_coefficients_typed(elements, coords, coefficients_json, mo_index, spacing)` → `Float32Array`
- `compute_esp_grid_typed(elements, coords, spacing, padding)` → `Float64Array`
- `compute_esp_grid_info(elements, coords, spacing, padding)` → JSON with origin, spacing, dims

### `compute_dos(elements, coords, sigma, e_min, e_max, n_points)` → DosResult
Density of States with Gaussian broadening.

- `sigma: f64` — broadening width in eV (default 0.3)
- `e_min/e_max: f64` — energy window in eV
- `n_points: usize` — grid points

**Returns:** `energies`, `total_dos`, `pdos`, `sigma`

### `eht_or_uff_fallback(smiles, elements, coords, allow_experimental_eht)` → ElectronicWorkflowResult
Smart router: runs EHT or falls back to UFF based on element support.

### `compare_methods(smiles, elements, coords, allow_experimental_eht)` → MethodComparisonResult
Benchmarks all available methods on the same geometry.

---

## 4. Charges & Surface

### `compute_charges(smiles)` → ChargeResult
Gasteiger-Marsili iterative partial charges.

**Returns:** `charges: Vec<f64>`, `iterations: usize`, `total_charge: f64`

### `compute_charges_configured(smiles, config)` → ChargeResult
Configurable Gasteiger-Marsili. `GasteigerConfig { max_iter: usize, initial_damping: f64, convergence_threshold: f64 }`.

### `compute_sasa(elements, coords, probe_radius)` → SasaResult
Shrake-Rupley solvent-accessible surface area.

- `probe_radius: f64` — probe radius in Å (default 1.4)

**Returns:** `total_sasa: f64`, `atom_sasa: Vec<f64>`, `probe_radius`, `num_points`

### `compute_dipole(elements, coords)` → DipoleResult
Molecular dipole moment from EHT density.

**Returns:** `vector: [f64; 3]`, `magnitude: f64`, `unit: "Debye"`

### `compute_esp(elements, coords, spacing, padding)` → VolumetricGrid
Electrostatic potential on 3D grid (Coulomb approximation).

- `spacing: f64` — grid spacing in Å
- `padding: f64` — padding beyond molecule in Å

---

## 5. Population & Bond Orders

### `compute_population(elements, coords)` → PopulationResult
Mulliken and Löwdin population analysis from EHT density matrix.

**Returns:** `mulliken_charges`, `lowdin_charges`, `total_charge_mulliken`, `total_charge_lowdin`

### `compute_bond_orders(elements, coords)` → BondOrderResult
Wiberg and Mayer bond orders from EHT density matrix.

**Returns:** `atom_pairs`, `distances`, `wiberg`, `mayer`, `wiberg_valence`, `mayer_valence`

---

## 6. Reactivity

### `compute_frontier_descriptors(elements, coords)` → FrontierDescriptors
HOMO/LUMO atom contributions and dual descriptor.

**Returns:** `homo_atom_contributions`, `lumo_atom_contributions`, `dual_descriptor`, `homo_energy`, `lumo_energy`, `gap`

### `compute_fukui_descriptors(elements, coords)` → FukuiDescriptors
Full and condensed Fukui functions f⁺, f⁻, f⁰.

**Returns:** `f_plus`, `f_minus`, `f_radical`, `dual_descriptor`, `condensed_*`, `gap`, `validity_notes`

### `compute_reactivity_ranking(elements, coords)` → ReactivityRanking
Ranked atomic sites for nucleophilic, electrophilic, and radical attack.

**Returns:** `nucleophilic_attack_sites`, `electrophilic_attack_sites`, `radical_attack_sites` — each with `atom_index`, `score`

### `compute_empirical_pka(smiles)` → EmpiricalPkaResult
Graph+heuristic estimation of acidic and basic pKa sites.

**Returns:** `acidic_sites`, `basic_sites` — each with `atom_index`, `site_type`, `environment`, `estimated_pka`, `confidence`

---

## 7. Spectroscopy — UV-Vis

### `compute_uv_vis_spectrum(elements, coords, sigma, e_min, e_max, n_points)` → UvVisSpectrum
EHT-based exploratory UV-Vis spectrum (fast, low accuracy).

**Returns:** `energies_ev`, `intensities`, `peaks[]` with `energy_ev/wavelength_nm/from_mo/to_mo/intensity`, `sigma`, `notes`

### `compute_stda_uvvis(elements, coords, sigma, e_min, e_max, n_points, broadening)` → StdaUvVisSpectrum
Simplified Tamm-Dancoff Approximation UV-Vis with proper oscillator strengths.

- `broadening: str` — `"gaussian"` (default) or `"lorentzian"`

**Returns:** `energies_ev`, `wavelengths_nm`, `absorptivity`, `excitations[]` with `oscillator_strength/transition_dipole`, `broadening`

---

## 8. Spectroscopy — IR

### `compute_vibrational_analysis(elements, coords, method, step_size)` → VibrationalAnalysis
Numerical Hessian vibrational frequencies + IR intensities.

- `method: str` — `"eht"` (default), `"pm3"`, or `"xtb"`
- `step_size: f64` — finite-difference step in Å (default 0.005)

**Returns:** `modes[]` with `frequency_cm1/ir_intensity/displacement/is_real`, `n_real_modes`, `zpve_ev`, `method`

### `compute_ir_spectrum(analysis, gamma, wn_min, wn_max, n_points)` → IrSpectrum
Lorentzian-broadened IR spectrum from vibrational analysis.

- `gamma: f64` — line width in cm⁻¹ (default 15)
- `wn_min/wn_max: f64` — wavenumber range (default 400–4000 cm⁻¹)

**Returns:** `wavenumbers`, `intensities`, `peaks[]`, `gamma`

### `compute_ir_spectrum_broadened(analysis, gamma, wn_min, wn_max, n_points, broadening)` → IrSpectrum
Like `compute_ir_spectrum` with explicit broadening type selection.

- `broadening: &str` — `"lorentzian"` (default) or `"gaussian"`

---

## 9. Spectroscopy — NMR

### `predict_nmr_shifts(smiles)` → NmrShiftResult
ML/HOSE-graph prediction of ¹H and ¹³C chemical shifts (topology-only).

**Returns:** `h_shifts[]` and `c_shifts[]` with `atom_index/element/shift_ppm/environment/confidence`

### `predict_nmr_couplings(smiles, coords)` → Vec\<JCoupling\>
Predicts vicinal and geminal J-coupling constants (Hz).

- `coords` — optional 3D coords for dihedral-aware coupling; pass `[]` for topology estimate

**Returns:** `h1_index`, `h2_index`, `j_hz`, `n_bonds`, `coupling_type`

### `compute_nmr_spectrum(smiles, nucleus, gamma, ppm_min, ppm_max, n_points)` → NmrSpectrum
Complete broadened NMR spectrum.

- `nucleus: str` — `"1H"` (default, 0–12 ppm) or `"13C"` (0–220 ppm)
- `gamma: f64` — Lorentzian line width in ppm (default 0.02)

**Returns:** `ppm_axis`, `intensities`, `peaks[]` with `multiplicity/environment`

### `compute_hose_codes(smiles, max_radius)` → Vec\<HoseCode\>
HOSE sphere codes for atom environment encoding (used in NMR prediction).

- `max_radius: usize` — sphere depth, typically 2–4

### `compute_ensemble_j_couplings(smiles, conformer_coords, energies_kcal, temperature_k)` → Vec\<JCoupling\>
Boltzmann-weighted ensemble J-couplings across multiple conformers.

- `conformer_coords: &[Vec<f64>]` — flat coord array per conformer
- `energies_kcal: &[f64]` — relative energies for Boltzmann weighting
- `temperature_k: f64` — temperature in K (default 298.15)

---

## 10. Semi-Empirical QM (PM3, xTB)

### `compute_pm3(elements, coords)` → Pm3Result
PM3 semi-empirical Hamiltonian with SCF.

**Returns:** `orbital_energies`, `electronic_energy`, `nuclear_repulsion`, `total_energy`, `heat_of_formation`, `homo_energy`, `lumo_energy`, `gap`, `mulliken_charges`, `scf_iterations`, `converged`

### `compute_xtb(elements, coords)` → XtbResult
GFN0-xTB tight-binding with SCC.

**Returns:** `orbital_energies`, `total_energy`, `repulsion_energy`, `electronic_energy`, `homo_energy`, `lumo_energy`, `gap`, `mulliken_charges`, `converged`, `scf_iterations`

---

## 11. ML Properties

### `ml_descriptors(smiles)` → MolecularDescriptors
Topology-derived descriptors for ML pipelines.

**Returns:** `molecular_weight`, `n_heavy_atoms`, `n_hydrogens`, `n_bonds`, `n_rotatable_bonds`, `n_hbd`, `n_hba`, `fsp3`, `wiener_index`, `n_rings`, `n_aromatic`, `sum_electronegativity`, `sum_polarizability`

### `ml_predict(smiles)` → MlPropertyResult
Lipinski/ADME property prediction using ML proxy models.

**Returns:** `logp`, `molar_refractivity`, `log_solubility`, `lipinski_violations`, `lipinski_passes`, `druglikeness`

---

## 12. Topology & Graph

### `analyze_graph_features(smiles)` → GraphFeatureAnalysis
Aromaticity perception and stereocenter detection.

**Returns:** `aromatic_atoms: Vec<bool>`, `aromatic_bonds`, `tagged_tetrahedral_centers`, `inferred_tetrahedral_centers`

### `compute_topology(elements, coords)` → TopologyAnalysisResult
Metal coordination geometry classification.

**Returns:** `metal_centers[]` with `atom_index/element/ligand_indices/coordination_number/geometry/geometry_score`, `warnings`

---

## 13. Materials / Crystallography

### `create_unit_cell(a, b, c, alpha, beta, gamma)` → UnitCell
Creates a crystallographic unit cell. Angles in degrees (default 90°).

**Returns:** `a/b/c/alpha/beta/gamma`, `volume`, `lattice` matrix

### `assemble_framework(topology, metal, geometry, lattice_a, supercell)` → CrystalStructure
Assembles a MOF/framework crystal structure.

- `topology: str` — `"pcu"`, `"dia"`, or `"sql"`
- `metal: u8` — metal atomic number (e.g. 30 = Zn)
- `geometry: str` — `"linear"`, `"trigonal"`, `"tetrahedral"`, `"square_planar"`, `"octahedral"`
- `lattice_a: f64` — cubic lattice parameter in Å
- `supercell: usize` — replication factor (1 = no replication)

**Returns:** `num_atoms`, `elements`, `frac_coords`, `cart_coords`, `labels`, `lattice`

---

## 14. Transport & Batch Streaming

### `pack_conformers(results)` → RecordBatch
Packs results into Arrow-compatible columnar format for zero-copy transfer.

### `split_worker_tasks(smiles, n_workers, seed)` → Vec\<WorkerTask\>
Splits SMILES batch into balanced worker sub-tasks.

### `estimate_workers(n_items, max_workers)` → usize
Returns optimal worker count for a batch size.

---

## 15. System Planning

### `system_capabilities(elements)` → SystemCapabilities
Checks EHT, UFF, embed support for a specific element set.

**Returns:** `embed/uff/eht/population/orbital_grid` — each with `available`, `confidence`, `unsupported_elements`, `warnings`

### `system_method_plan(elements)` → SystemMethodPlan
Generates recommended method + fallback plan for a molecule.

**Returns:** structured `geometry/force_field_energy/orbitals/population/orbital_grid` plans with `recommended/fallback/rationale/methods[]`

---

## 16. Stereochemistry

### `analyze_stereo(smiles, coords)` → StereoAnalysis
CIP priority assignment, R/S chirality at tetrahedral centers, E/Z configuration at double bonds.

- `smiles: &str` — SMILES input
- `coords: &[f64]` — flat conformer coords (pass `&[]` for topology-only assignment)

**Returns:** `stereocenters[]` with `atom_index/element/substituent_indices/priorities/configuration`, `double_bonds[]` with `atom1/atom2/configuration/high_priority_sub1/high_priority_sub2`, `n_stereocenters`, `n_double_bonds`

Python: `stereo_analysis(smiles, coords=[])` · WASM: `analyze_stereo(smiles, coords_flat)` · CLI: `sci-form stereo "SMILES"`

---

## 17. Solvation

### `compute_nonpolar_solvation(elements, coords, probe_radius)` → NonPolarSolvation
Non-polar solvation free energy via SASA × atomic surface tension (ASP model).

- `probe_radius: Option<f64>` — Å (default 1.4)

**Returns:** `energy_kcal_mol`, `atom_contributions`, `atom_sasa`, `total_sasa`

### `compute_gb_solvation(elements, coords, charges, solvent_dielectric, solute_dielectric, probe_radius)` → GbSolvation
Generalized Born electrostatic solvation energy (HCT model).

- `solvent_dielectric: Option<f64>` — default 78.5 (water)
- `solute_dielectric: Option<f64>` — default 1.0 (vacuum core)
- `probe_radius: Option<f64>` — Å (default 1.4)

**Returns:** `electrostatic_energy_kcal_mol`, `nonpolar_energy_kcal_mol`, `total_energy_kcal_mol`, `born_radii`, `charges`, `solvent_dielectric`, `solute_dielectric`

Python: `nonpolar_solvation(elements, coords, probe_radius=1.4)` · `gb_solvation(elements, coords, charges, solvent_dielectric=78.5, ...)` · WASM: `compute_nonpolar_solvation(...)` · CLI: `sci-form solvation [ELEMENTS] [COORDS]`

---

## 18. Rings & Fingerprints

### `compute_sssr(smiles)` → SssrResult
Smallest Set of Smallest Rings via Horton's algorithm.

**Returns:** `rings[]` with `atoms/size/is_aromatic`, `atom_ring_count`, `atom_ring_sizes`, `ring_size_histogram`

### `compute_ecfp(smiles, radius, n_bits)` → ECFPFingerprint
Extended-Connectivity Fingerprints (ECFP / Morgan algorithm).

- `radius: usize` — default 2 (ECFP4)
- `n_bits: usize` — typically 1024 or 2048

**Returns:** `n_bits`, `on_bits: BTreeSet<usize>`, `radius`, `raw_features`

### `compute_tanimoto(fp1, fp2)` → f64
Tanimoto (Jaccard) coefficient between two ECFPs. Returns 0.0–1.0.

Python: `sssr(smiles)` · `ecfp(smiles, radius=2, n_bits=2048)` · `tanimoto(fp1, fp2)` · WASM: JSON-in/JSON-out · CLI: `sci-form sssr "SMILES"`, `sci-form ecfp "SMILES"`, `sci-form tanimoto "SM1" "SM2"`

---

## 19. Clustering

### `butina_cluster(conformers, rmsd_cutoff)` → ClusterResult
Taylor-Butina leader clustering on flat conformer coordinate arrays.

- `conformers: &[Vec<f64>]` — list of flat coordinate arrays (same atom count)
- `rmsd_cutoff: f64` — cluster radius in Å (typically 1.0–2.0)

**Returns:** `n_clusters`, `assignments`, `centroid_indices`, `cluster_sizes`, `rmsd_cutoff`

### `compute_rmsd_matrix(conformers)` → Vec\<Vec\<f64\>\>
Pair-wise RMSD matrix for all conformers.

### `filter_diverse(conformers, rmsd_cutoff)` → Vec\<usize\>  *(Python-only)*
Greedy diversity filter — returns indices of maximally diverse conformers.

Python: `butina_cluster(conformers, rmsd_cutoff=1.0)` · `rmsd_matrix(conformers)` · `filter_diverse(conformers, rmsd_cutoff=1.0)` · WASM: `butina_cluster(conformers_json, rmsd_cutoff)` · `compute_rmsd_matrix(conformers_json)`

---

## Implementation Examples by Language

See the `examples/` directory in this skill for complete runnable examples:

| File | Language | Content |
|---|---|---|
| `examples/rust_examples.rs` | Rust | All major API groups |
| `examples/js_examples.js` | JavaScript/WASM | Full browser + Node usage with parallel init |
| `examples/python_examples.py` | Python | All functions with type annotations |
| `examples/cli_examples.sh` | Bash/CLI | All CLI subcommands with options |
