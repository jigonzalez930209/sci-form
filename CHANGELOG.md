# Changelog

All notable changes to sci-form are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.10.6] - 2026-03-26

### Added

#### Quantum Chemistry
- **GFN1-xTB** — shell-resolved self-consistent charges with D3-BJ dispersion correction
- **GFN2-xTB** — multipole electrostatics, D4 dispersion, and halogen-bond correction
- **PM3(tm)** — PM3 transition-metal parameter set (Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Pd, Pt, Au and more)
- **HF-3c** — minimal-basis Hartree-Fock with D3, gCP, and SRB corrections
- **CISD** — Configuration Interaction Singles+Doubles for excited states (via HF-3c orbital basis)
- **EEQ charges** — geometry-dependent electronegativity-equilibration charges with D4 coupling
- **D4 dispersion** — Caldeweyher D4 energy and analytic gradients

#### Neural Network Potentials
- **ANI-TM** — ANI neural network potential extended to 24 elements: H, C, N, O, F, Si, P, S, Cl, Ti, Cr, Mn, Fe, Co, Ni, Cu, Zn, Br, Ru, Pd, Ag, I, Pt, Au

#### Spectroscopy
- **Multi-method DOS** — unified `compute_dos_multimethod` supporting EHT, PM3, xTB, GFN1, GFN2, and HF-3c methods
- **IR peak assignment** — functional-group recognition from normal mode displacements (`assign_peaks`)
- **NMR ensemble J-coupling averaging** — Boltzmann-weighted coupling constants across conformer ensembles

#### Population & Bond Analysis
- **NPA** — Natural Population Analysis (`compute_npa`)
- **NBO** — Natural Bond Orbital analysis (`compute_nbo`)
- **Bond orders** — Wiberg and Mayer bond order matrices

#### Solvation
- **Generalized Born (HCT)** — electrostatic solvation with configurable solvent and solute dielectrics
- **ALPB** — Analytical Linearised Poisson-Boltzmann solvation model

#### Machine Learning
- **3D descriptors** — WHIM (5 weightings: unit, mass, volume, electronegativity, polarizability)
- **3D descriptors** — RDF (radial distribution function descriptors with charge/mass/electronegativity weighting)
- **3D descriptors** — GETAWAY (leverage matrix, hat autocorrelation, R autocorrelation)
- **Random Forest** — decision-tree ensemble with `predict_with_variance` uncertainty estimation
- **Gradient Boosting** — boosted tree regression with configurable learning rate and depth
- **Cross-validation** — k-fold cross-validation for ML model evaluation

#### Cheminformatics
- **Stereochemistry analysis** — CIP priorities, R/S tetrahedral centers, E/Z double bonds, atropisomeric M/P axes, helical chirality
- **SSSR** — Smallest Set of Smallest Rings via Horton's algorithm
- **ECFP fingerprints** — Extended-Connectivity Fingerprints (Morgan algorithm) with deterministic FNV-1a hashing
- **Tanimoto similarity** — Jaccard coefficient on ECFP bit vectors
- **Butina clustering** — Taylor-Butina RMSD-based molecular clustering
- **Diversity filtering** — select maximally diverse conformers by RMSD cutoff
- **SMIRKS engine** — reaction transforms via atom-mapped reactant→product SMIRKS patterns
- **Canonical SMILES** — deterministic canonical output from molecular graph

#### Materials & Crystallography
- **230 ITC space groups** — complete symmetry operations and equivalent-position generation for all 230 space groups
- **Framework geometry optimization** — BFGS and steepest-descent optimizer with periodic boundary conditions
- **Periodic systems** — PBC molecular graph construction with configurable bond tolerance for TM-TM pairs
- **Hapticity detection** — metallocene and η-coordination detection in periodic structures

#### EHT Extensions
- **Analytical gradients** — EHT energy gradients with respect to atomic coordinates
- **Geometry optimization** — EHT-based geometry optimizer (`optimize_geometry_eht`)
- **k-path band structure** — EHT band structure along high-symmetry k-paths for periodic systems

#### Molecular Dynamics
- Trajectory computation module for time-evolved molecular configurations

#### Reactivity
- **Fukui descriptors** — condensed Fukui f+, f−, and f° from finite-difference frontier-orbital densities
- **Frontier descriptors** — HOMO/LUMO atom contributions and dual descriptor
- **Empirical pKa** — rule-based acidic and basic site pKa estimation

#### Experimental / Research (feature-gated)
- **GPU** (`experimental-gpu`) — MMFF94 non-bonded GPU kernels (WGSL), EEQ GPU kernels, `GpuContext`
- **Alpha** (`alpha`) — CGA (conformal geometric algebra), GSM (growing string method for reaction paths), SDR (spectral dimensionality reduction)
- **Beta** (`beta`) — KPM (kernel polynomial method), MBH (mobile block Hessian), CPM (correlation participation matrix), RandNLA (randomized numerical linear algebra), Riemannian optimization

### Fixed
- **EHT/PM3/xTB field naming** — consistent snake_case field names across Rust, Python, WASM, and CLI bindings
- **ECFP determinism** — replaced HashMap with FNV-1a hasher for reproducible fingerprint bit assignments across platforms
- **Periodic bond tolerance** — improved default cutoffs for TM-TM and TM-ligand bond detection in PBC graphs
- **Helical chirality with heteroatoms** — corrected backbone-path detection when heteroatoms appear in the helix axis

### Changed
- **xTB field names** — `repulsive_energy` (was `repulsion_energy`), `scc_iterations` (was `scf_iterations`)
- **All result types** — implement `serde::Serialize` and `serde::Deserialize`
- **`parallel` feature** — enabled by default in library builds; opt-out via `default-features = false`
- **WASM API** — `elements` and `coords_flat` parameters are always JSON-encoded strings; typed-array paths (`embed_coords_typed`, `eht_orbital_grid_typed`, `compute_esp_grid_typed`) remain available for performance-critical paths

---

## [0.9.1] - 2025-11-14

### Added
- Initial public release of GFN0-xTB (25 elements including transition metals)
- ANI-2x neural network potential (H, C, N, O, F, S, Cl)
- sTDA UV-Vis with oscillator strengths and Gaussian/Lorentzian broadening
- IR spectroscopy via numerical Hessian (3N normal modes, ZPVE, thermochemistry via RRHO)
- NMR ¹H/¹³C/¹⁹F/³¹P chemical shifts (HOSE codes) and J-coupling (Karplus equation)
- Non-polar SASA solvation energy
- UFF and MMFF94 force field energy evaluation
- Kabsch SVD and quaternion (Coutsias 2004) molecular alignment
- Apache Arrow columnar transport for batch results
- Web Worker task-splitting utilities for WASM parallel pipelines
- 2D ML molecular descriptors (17+ topology descriptors, LogP, MR, solubility, Lipinski Ro5)
- PM3 NDDO semi-empirical SCF (heat of formation, HOMO/LUMO, Mulliken charges)
- EHT orbital grids and marching-cubes isosurfaces
- EHT Density of States (total + per-atom PDOS, Gaussian smearing)
- Electrostatic potential (Coulomb grid, `.cube` export, parallel evaluation)
- Gasteiger-Marsili partial charges
- SASA (Shrake-Rupley, Bondi radii)
- Materials: unit cells, SBU topology, MOF-type framework assembly (pcu, dia, sql)
- Python bindings (`sciforma` on PyPI) via PyO3
- WASM bindings (`sci-form-wasm` on npm) with typed-array fast paths
- CLI (`sci-form-cli` on crates.io) with JSON/XYZ/SDF output

### Benchmarks
- 60+ conformers/second on a single core (ETKDGv2)
- 0.064 Å average RMSD vs RDKit on diverse molecule set
- 0.00% heavy-atom RMSD > 0.5 Å on GDB-20 ensemble benchmark (2000 molecules)

---

*Older history predates public release and is not included.*
