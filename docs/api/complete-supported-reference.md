# Complete Supported API and Algorithm Reference

This document is the consolidated source for the **currently supported API surface** of sci-form across Rust, Python, WASM/TypeScript, and CLI, plus the **algorithmic foundations** and the **existing SVG assets** that explain those algorithms.

It is intentionally aligned to the current repository state instead of older aspirational docs. When a function is listed here, it is because it is exposed today in source. When a capability is experimental or feature-gated, that status is called out explicitly.

---

## Scope and Source of Truth

Primary source files used to build this reference:

- `src/lib.rs` for the top-level Rust API
- `crates/python/src/*.rs` for Python bindings
- `crates/wasm/src/*.rs` for WASM/TypeScript exports
- `crates/cli/src/args.rs` for CLI commands
- `docs/algorithm/*.md` for theory pages
- `docs/public/svg/*.svg` and `docs/public/*.svg` for diagrams

If any future API change diverges from this document, the source files above win.

---

## Shared Conventions

### Coordinates

- `coords` are flat arrays: `[x0, y0, z0, x1, y1, z1, ...]`
- When Rust APIs require `&[[f64; 3]]`, convert from the flat array with `chunks(3)`
- Units for coordinates are **Å**

### Elements

- `elements` are atomic numbers (`u8` in Rust, integer lists in Python and WASM JSON)

### Surface rules

- Rust returns native structs and `Result<_, String>` where appropriate
- Python returns PyO3 wrapper classes with attribute access
- WASM uses JSON strings for structured inputs/outputs, with typed-array fast paths for selected operations
- CLI exposes a subset of stable workflows as commands

### Stability model

- **Stable/default**: part of the ordinary build and intended for production usage
- **Feature-gated experimental**: only available when the corresponding feature is enabled
- **Surface-limited**: exposed only on one or more surfaces, not necessarily on all four

---

## Support Summary by Surface

| Surface | Status | Notes |
|--------|--------|-------|
| Rust core | Fullest public surface | Includes advanced functions not mirrored everywhere else |
| Python | Broad high-level surface | Closely mirrors stable Rust workflows with Pythonic wrappers |
| WASM / TypeScript | Broad browser-friendly surface | JSON-based API plus typed arrays for performance-critical paths |
| CLI | Curated operational subset | Focused on common workflows, not the full scientific surface |

---

## Function Catalog by Domain

The tables below are organized function-by-function by domain. A dash means the function is not directly exposed on that surface.

### 1. Geometry, Parsing, and Canonicalization

| Rust | Python | WASM | CLI | Purpose | Algorithmic basis |
|------|--------|------|-----|---------|-------------------|
| `version` | `version` | `version` | `info` | Library/build version and platform info | Metadata only |
| `embed` | `embed` | `embed` | `embed` | Single conformer generation from SMILES | ETKDGv2 distance geometry + refinement |
| `embed_batch` | `embed_batch` | `embed_batch` | `batch` | Batch conformer generation | ETKDGv2 + rayon parallelization |
| `embed_diverse` | - | - | - | Multiple embeddings filtered by diversity | ETKDGv2 + Butina RMSD clustering |
| `parse` | `parse` | `parse_smiles` | `parse` | Parse SMILES into graph/topology without 3D | SMILES parsing + graph construction |
| `to_canonical_smiles` | - | - | - | Canonicalize a SMILES representation | Canonical graph traversal / ordering |
| `compute_rmsd` | `rmsd` | `compute_rmsd` | `rmsd` | RMSD after optimal superposition | Kabsch alignment |

Practical note:

- `embed` always returns a result object in Rust/Python/WASM. Failure is expressed in the `error` field rather than by panicking.

### 2. Support Planning and Method Routing

| Rust | Python | WASM | CLI | Purpose | Algorithmic basis |
|------|--------|------|-----|---------|-------------------|
| `get_eht_support` | `eht_support` | `eht_support` | - | Check EHT applicability for a given element set | Element support table + confidence model |
| `get_system_capabilities` | `system_capabilities` | `system_capabilities` | - | Report what core workflows are usable for an element set | Capability aggregation |
| `get_system_method_plan` | `system_method_plan` | `system_method_plan` | - | Recommend method/fallback by property domain | Heuristic planner |
| `compute_eht_or_uff_fallback` | `eht_or_uff_fallback` | `eht_or_uff_fallback` | - | Use EHT when supported, else UFF fallback | Support-aware workflow routing |
| `compare_methods` | `compare_methods` | `compare_methods` | - | Compare method availability/results | Structured multi-method comparison |

### 3. Charges, Surface, Population, Electrostatics, and DOS

| Rust | Python | WASM | CLI | Purpose | Algorithmic basis |
|------|--------|------|-----|---------|-------------------|
| `compute_charges` | `charges` | `compute_charges` | `charges` | Gasteiger-Marsili partial charges | Electronegativity equalization |
| `compute_charges_configured` | - | - | - | Tunable Gasteiger-Marsili calculation | Same model, configurable damping/iterations |
| `compute_sasa` | `sasa` | `compute_sasa` | `sasa` | Solvent-accessible surface area | Shrake-Rupley + Fibonacci sphere |
| `compute_population` | `population` | `compute_population` | `population` | Mulliken/Löwdin population analysis | EHT + overlap-based charge partitioning |
| `compute_dipole` | `dipole` | `compute_dipole` | `dipole` | Molecular dipole moment | EHT-derived electronic distribution |
| `compute_esp` | `esp` | `compute_esp` | `esp` | Electrostatic potential on a 3D grid | Coulomb grid from atomic charges |
| `compute_esp_grid` | - | - | - | Convenience alias returning full ESP grid struct | Same as `compute_esp` |
| `compute_dos` | `dos` | `compute_dos` | `dos` | Density of states curve | Gaussian broadening over orbital energies |
| `compute_bond_orders` | `bond_orders` | `compute_bond_orders` | - | Wiberg-like and Mayer-like bond orders | EHT density/overlap analysis |
| `compute_frontier_descriptors` | `frontier_descriptors` | `compute_frontier_descriptors` | - | Atom-resolved HOMO/LUMO descriptors | EHT MO contributions |
| `compute_fukui_descriptors` | `fukui_descriptors` | `compute_fukui_descriptors` | - | Condensed Fukui descriptors | Frontier orbital-based reactivity workflow |
| `compute_reactivity_ranking` | `reactivity_ranking` | `compute_reactivity_ranking` | - | Rank reactive sites | Fukui + Mulliken hybrid heuristic |
| `compute_uv_vis_spectrum` | `uvvis_spectrum` | `compute_uv_vis_spectrum` | - | Low-cost exploratory UV-Vis-like spectrum | EHT transition approximation |

### 4. Force Fields and Heuristic Energy Models

| Rust | Python | WASM | CLI | Purpose | Algorithmic basis |
|------|--------|------|-----|---------|-------------------|
| `compute_uff_energy` | `uff_energy` | `compute_uff_energy` | `uff` | UFF energy evaluation | Universal Force Field |
| `compute_uff_energy_with_aromatic_heuristics` | `uff_energy_with_aromatic_heuristics` | `compute_uff_energy_with_aromatic_heuristics` | - | UFF with aromaticity-informed correction | UFF + aromatic post-correction |
| `compute_mmff94_energy` | `mmff94_energy` | `compute_mmff94_energy` | - | MMFF94 energy evaluation | Halgren MMFF94 |
| `compute_empirical_pka` | `empirical_pka` | `compute_empirical_pka` | - | Estimate acidic/basic pKa sites | Graph environment + charge heuristics |

### 5. Semiempirical, Tight-Binding, Ab Initio, and Neural Potentials

| Rust | Python | WASM | CLI | Purpose | Algorithmic basis |
|------|--------|------|-----|---------|-------------------|
| `compute_pm3` | `pm3_calculate` | `compute_pm3` | - | PM3 semiempirical SCF | NDDO PM3 |
| `compute_pm3_gradient` | - | - | - | Analytical PM3 gradients | Hellmann-Feynman + Pulay corrections |
| `compute_xtb` | `xtb_calculate` | `compute_xtb` | - | GFN0-xTB style tight-binding workflow | SCC tight-binding |
| `compute_xtb_gradient` | - | - | - | Analytical xTB gradients | SCC gradient formalism |
| `compute_hf3c` | `hf3c_calculate` | `compute_hf3c` / `compute_hf3c_custom` | `hf3c` | HF-3c quantum workflow | HF + D3 + gCP + SRB |
| `compute_ani` | `ani_calculate` | `compute_ani` | `ani` | ANI neural potential using available test/default models | Atomic environment vectors + feed-forward NN |
| `compute_ani_with_models` | - | - | - | ANI with explicit models/config | Same ANI formalism with custom weights |

### 6. Spectroscopy and Orbital Visualization

| Rust | Python | WASM | CLI | Purpose | Algorithmic basis |
|------|--------|------|-----|---------|-------------------|
| `compute_stda_uvvis` | `stda_uvvis` | `compute_stda_uvvis` | - | UV-Vis spectrum with oscillator strengths | sTDA-style transition workflow |
| `compute_vibrational_analysis` | `vibrational_analysis` | `compute_vibrational_analysis` | - | Vibrational frequencies, modes, IR intensities | Numerical Hessian |
| `compute_vibrational_analysis_uff` | - | - | - | UFF-based vibrational workflow | Gradient-difference Hessian |
| `compute_ir_spectrum` | `ir_spectrum` | `compute_ir_spectrum` | - | Broadened IR spectrum from vibrational analysis | Lorentzian broadening |
| `compute_ir_spectrum_broadened` | - | - | - | IR with selectable broadening type | Gaussian or Lorentzian broadening |
| `predict_nmr_shifts` | `nmr_shifts` | `predict_nmr_shifts` | - | Chemical shift prediction | HOSE/environment heuristics |
| `predict_nmr_couplings` | `nmr_couplings` | `predict_nmr_couplings` | - | J-coupling prediction | Karplus + topological heuristics |
| `compute_nmr_spectrum` | `nmr_spectrum` | `compute_nmr_spectrum` | - | Full synthetic NMR spectrum | Shift prediction + coupling model + line broadening |
| `compute_hose_codes` | `hose_codes` | `compute_hose_codes` | - | Generate HOSE codes per atom | Sphere-of-environment graph encoding |
| `compute_ensemble_j_couplings` | - | - | - | Boltzmann-averaged J-couplings across conformers | Ensemble averaging |
| `compute_orbital_mesh` | `orbital_mesh` | `compute_orbital_mesh` | - | Isosurface mesh for orbitals | Grid evaluation + marching cubes |
| - | `eht_orbital_mesh` | `eht_orbital_mesh` | - | EHT-specific orbital isosurface | Same marching-cubes family |
| - | - | `eht_orbital_grid_typed` | - | EHT orbital grid as typed array | 3D volumetric evaluation |
| - | - | `eht_orbital_grid_from_coefficients_typed` | - | Orbital grid from supplied MO coefficients | 3D volumetric evaluation |
| - | - | `eht_volumetric_grid_info` | - | Metadata for EHT volumetric grids | Grid bookkeeping |

### 7. Topology, Graph Analysis, and ML Properties

| Rust | Python | WASM | CLI | Purpose | Algorithmic basis |
|------|--------|------|-----|---------|-------------------|
| `compute_topology` | `topology_analysis` | `compute_topology` | - | Coordination/topology analysis for structures | Geometric coordination heuristics |
| `analyze_graph_features` | `graph_features` | `analyze_graph_features` | - | Aromaticity and stereotag graph analysis | Graph traversal and aromatic bond tagging |
| `compute_ml_descriptors` | `ml_descriptors` | `compute_molecular_descriptors` | - | Molecular descriptor generation | Topological descriptor extraction |
| `predict_ml_properties` | `ml_predict` | `compute_ml_properties` | - | Property prediction from descriptors | Proxy ML models |
| `predict_ensemble` | - | - | - | Ensemble property prediction with confidence | Multi-model aggregation |
| `compute_tpsa` | - | - | - | Topological polar surface area | Fragment-based TPSA heuristic |

### 8. Dynamics, Conformer Search, and Path Sampling

| Rust | Python | WASM | CLI | Purpose | Algorithmic basis |
|------|--------|------|-----|---------|-------------------|
| `compute_md_trajectory` | `md_trajectory` | `compute_md_trajectory` | - | Short exploratory MD | Velocity Verlet |
| `compute_md_trajectory_nvt` | `md_trajectory_nvt` | `compute_md_trajectory_nvt` | - | MD with thermostat | Velocity Verlet + Berendsen NVT |
| `compute_simplified_neb_path` | `simplified_neb_path` | `compute_simplified_neb_path` | - | Approximate reaction path between geometries | Simplified NEB |
| `search_conformers_with_uff` | `search_conformers` | `search_conformers_with_uff` | - | UFF-ranked conformer exploration | Multi-seed embedding + UFF scoring + clustering |

### 9. Materials, Periodic Systems, and Framework Construction

| Rust | Python | WASM | CLI | Purpose | Algorithmic basis |
|------|--------|------|-----|---------|-------------------|
| `create_unit_cell` | `unit_cell` | `create_unit_cell` | `cell` | Build a unit cell from lattice parameters | Crystallographic cell tensor construction |
| `assemble_framework` | `assemble_framework` | `assemble_framework` | `assemble` | Assemble framework structures from topology/SBUs | Topology-guided framework assembly |

Additional materials/periodic APIs currently documented in the Rust layer and repository docs include:

- `materials::space_groups::space_group_by_number`
- `materials::space_groups::space_group_by_symbol`
- `periodic::build_periodic_molecule`
- `periodic::detect_hapticity`
- `materials::geometry_opt::optimize_framework`

These are part of the broader Rust crate surface even when not mirrored one-to-one in Python, WASM, or CLI.

### 10. Stereochemistry, Solvation, Rings, Fingerprints, and Clustering

| Rust | Python | WASM | CLI | Purpose | Algorithmic basis |
|------|--------|------|-----|---------|-------------------|
| `analyze_stereo` | `stereo_analysis` | `analyze_stereo` | `stereo` | R/S and E/Z assignment | CIP priority workflow + geometry checks |
| `compute_nonpolar_solvation` | `nonpolar_solvation` | `compute_nonpolar_solvation` | `solvation` | Non-polar solvation energy | SASA + atomic solvation parameters |
| `compute_gb_solvation` | `gb_solvation` | `compute_gb_solvation` | `solvation --charges ...` | Electrostatic + nonpolar solvation | HCT-style Generalized Born |
| `compute_sssr` | `sssr` | `compute_sssr` | `sssr` | Smallest set of smallest rings | Horton-style ring perception |
| `compute_ecfp` | `ecfp` | `compute_ecfp` | `ecfp` | Morgan/ECFP fingerprint | Circular neighborhood hashing |
| `compute_ecfp_batch` | - | - | - | Batch ECFP generation | Same Morgan/ECFP workflow |
| `compute_tanimoto` | `tanimoto` | `compute_tanimoto` | `tanimoto` | Similarity between fingerprints | Tanimoto coefficient |
| `butina_cluster` | `butina_cluster` | `butina_cluster` | - | Conformer clustering | Taylor-Butina clustering |
| `compute_rmsd_matrix` | `rmsd_matrix` | `compute_rmsd_matrix` | - | Pairwise RMSD matrix | Repeated alignment metric |
| - | `filter_diverse` | - | - | Keep only diverse conformers | Diversity filtering over clusters |

### 11. Transport and Batch-Oriented Helpers

| Rust | Python | WASM | CLI | Purpose | Algorithmic basis |
|------|--------|------|-----|---------|-------------------|
| `transport` module helpers | `pack_conformers` | `pack_batch_arrow` | - | Columnar transport of conformer results | Arrow-like packing |
| `transport` module helpers | `split_worker_tasks` | `split_worker_tasks` | - | Split work for worker pools | Chunked task partitioning |
| `transport` module helpers | `estimate_workers` | `estimate_workers` | - | Estimate worker count for batch size | Heuristic load balancing |

---

## CLI Coverage

The CLI is a curated subset intended for common user workflows rather than total feature parity.

### Stable CLI commands

- `embed`
- `batch`
- `parse`
- `info`
- `eht`
- `charges`
- `sasa`
- `population`
- `dipole`
- `esp`
- `dos`
- `rmsd`
- `uff`
- `cell`
- `assemble`
- `ani`
- `hf3c`
- `stereo`
- `solvation`
- `sssr`
- `ecfp`
- `tanimoto`

### Feature-gated CLI commands

- `eeq` behind `experimental-eeq`
- `alpb` behind `experimental-alpb`
- `d4` behind `experimental-d4`
- `cpm` behind `experimental-cpm`

---

## Algorithm-by-Algorithm Guide

This section explains the supported algorithms at a practical and theoretical level and connects each one to the exposed API.

### SMILES Parsing and Molecular Graph Construction

Core idea:

- Parse SMILES into a graph with atoms, bonds, formal charges, aromaticity, and stereotags.

Main API entries:

- Rust: `parse`, `to_canonical_smiles`, `analyze_graph_features`
- Python: `parse`, `graph_features`
- WASM: `parse_smiles`, `analyze_graph_features`
- CLI: `parse`

Related theory docs:

- `/algorithm/smiles-parsing`
- `/algorithm/smarts-matching`

Related SVGs:

- `/svg/smiles-pipeline.svg`
- `/svg/smiles-phenol.svg`
- `/svg/smarts-pattern.svg`

### ETKDGv2 Conformer Generation

Core idea:

- Build a bounds matrix from graph topology, smooth it, sample distances, embed coordinates, then refine with ETKDG/UFF-like corrections and validation.

Main API entries:

- Rust: `embed`, `embed_batch`, `embed_diverse`
- Python: `embed`, `embed_batch`
- WASM: `embed`, `embed_coords`, `embed_coords_typed`, `embed_batch`
- CLI: `embed`, `batch`

Key sub-algorithms:

1. SMILES-to-graph parsing
2. Bounds-matrix construction for 1-2, 1-3, 1-4, and nonbonded pairs
3. Floyd-Warshall-style triangle smoothing
4. Distance geometry via metric matrix and eigendecomposition
5. BFGS-style bounds violation minimization
6. ETKDG torsion refinement and stereochemical validation
7. Retry loop with fallback placement when needed

Related theory docs:

- `/algorithm/overview`
- `/algorithm/bounds-matrix`
- `/algorithm/distance-geometry`
- `/algorithm/embedding`
- `/algorithm/etkdg-refinement`
- `/algorithm/validation`

Related SVGs:

- `/svg/pipeline-overview.svg`
- `/svg/etkdg-pipeline.svg`
- `/svg/bounds-matrix-overview.svg`
- `/svg/bounds-matrix-1-4.svg`
- `/svg/bounds-matrix-triangle.svg`
- `/svg/distance-geometry-pipeline.svg`
- `/svg/distance-geometry-matrix.svg`
- `/svg/embedding-eigendecomp.svg`
- `/svg/etkdg-torsion.svg`
- `/svg/validation-pipeline.svg`
- `/svg/validation-tetrahedral.svg`
- `/svg/validation-double-bond.svg`

### Gasteiger-Marsili Charges

Core idea:

- Iteratively equalize electronegativities across the graph using topology and formal charges, without a quantum calculation.

Main API entries:

- Rust: `compute_charges`, `compute_charges_configured`
- Python: `charges`
- WASM: `compute_charges`
- CLI: `charges`

Related theory docs:

- `/algorithm/population-analysis`

Related SVGs:

- `/svg/population-mulliken.svg`
- `/svg/population-lowdin.svg`

### Extended Hückel Theory (EHT)

Core idea:

- Build a minimal STO-3G-like basis, compute overlap, assemble a Wolfsberg-Helmholtz Hamiltonian, apply Löwdin orthogonalization, diagonalize, and populate orbitals.

Main API entries:

- Rust: `get_eht_support`, `compute_population`, `compute_dipole`, `compute_frontier_descriptors`, `compute_fukui_descriptors`, `compute_reactivity_ranking`, `compute_uv_vis_spectrum`, `compute_bond_orders`, `compute_esp`, `compute_dos`
- Python: `eht_calculate`, `eht_support`, `eht_or_uff_fallback`, `population`, `dipole`, `bond_orders`, `frontier_descriptors`
- WASM: `eht_calculate`, `eht_support`, `eht_or_uff_fallback`, `compute_population`, `compute_dipole`, `compute_bond_orders`, `compute_dos`
- CLI: `eht`, `population`, `dipole`, `dos`, `esp`

Why it matters in sci-form:

- EHT is the gateway calculation behind multiple “derived” properties.

Related theory docs:

- `/algorithm/extended-huckel-theory`
- `/algorithm/population-analysis`
- `/algorithm/dipole-moments`
- `/algorithm/density-of-states`
- `/algorithm/electrostatic-potential`

Related SVGs:

- `/svg/eht-pipeline.svg`
- `/svg/population-mulliken.svg`
- `/svg/population-lowdin.svg`
- `/svg/dipole-vector.svg`
- `/svg/esp-grid.svg`
- `/svg/dos-smearing.svg`

### Orbital Grids and Isosurfaces

Core idea:

- Evaluate orbital amplitudes on a 3D grid and convert scalar fields to a mesh with marching cubes.

Main API entries:

- Rust: `compute_orbital_mesh`
- Python: `eht_orbital_mesh`, `orbital_mesh`
- WASM: `eht_orbital_mesh`, `eht_orbital_grid_typed`, `eht_orbital_grid_from_coefficients_typed`, `compute_orbital_mesh`

Related theory docs:

- `/algorithm/extended-huckel-theory`

Related SVGs:

- `/svg/eht-pipeline.svg`
- `/svg/embedding-eigendecomp.svg`

### Force Fields: UFF and MMFF94

Core idea:

- Model molecular strain energy using bonded and nonbonded terms; UFF emphasizes broad element coverage, MMFF94 emphasizes medicinal-chemistry fidelity.

Main API entries:

- Rust: `compute_uff_energy`, `compute_uff_energy_with_aromatic_heuristics`, `compute_mmff94_energy`
- Python: `uff_energy`, `uff_energy_with_aromatic_heuristics`, `mmff94_energy`
- WASM: same names with `compute_...`
- CLI: `uff`

Related theory docs:

- `/algorithm/force-fields`
- `/algorithm/strain-energy`

Related SVGs:

- `/svg/force-fields-overview.svg`
- `/svg/force-fields-torsion.svg`
- `/svg/uff-energy-terms.svg`
- `/svg/mmff94-energy-terms.svg`

### PM3 and xTB

Core idea:

- PM3 uses NDDO semiempirical SCF; xTB uses self-consistent charge tight binding for faster approximate electronic structure.

Main API entries:

- Rust: `compute_pm3`, `compute_pm3_gradient`, `compute_xtb`, `compute_xtb_gradient`
- Python: `pm3_calculate`, `xtb_calculate`
- WASM: `compute_pm3`, `compute_xtb`

Related theory docs:

- currently summarized from source and repo-level docs; no dedicated standalone algorithm page yet in `docs/algorithm/`

Related SVGs:

- `/svg/spectroscopy-overview.svg`

### HF-3c

Core idea:

- Minimal-basis Hartree-Fock corrected by D3, gCP, and short-range basis corrections for an efficient ab initio-like workflow.

Main API entries:

- Rust: `compute_hf3c`
- Python: `hf3c_calculate`
- WASM: `compute_hf3c`, `compute_hf3c_custom`
- CLI: `hf3c`

Related theory docs:

- `/algorithm/hf3c-quantum-engine`

Related SVGs:

- `/hf3c-architecture.svg`

### ANI Neural Potentials

Core idea:

- Convert local atomic environments into AEV features and predict energies and forces with neural networks.

Main API entries:

- Rust: `compute_ani`, `compute_ani_with_models`
- Python: `ani_calculate`
- WASM: `compute_ani`
- CLI: `ani`

Related theory docs:

- `/algorithm/machine-learning-potentials`

Related SVGs:

- `/ani-architecture.svg`

### UV-Vis, IR, and NMR Spectroscopy

Core idea:

- UV-Vis: electronic transitions with broadening
- IR: vibrational frequencies and intensities from Hessian workflows
- NMR: chemically informed shift predictions, J-coupling estimation, and synthetic line broadening

Main API entries:

- UV-Vis: `compute_uv_vis_spectrum`, `compute_stda_uvvis`, `uvvis_spectrum`, `stda_uvvis`, `compute_uv_vis_spectrum`, `compute_stda_uvvis`
- IR: `compute_vibrational_analysis`, `compute_vibrational_analysis_uff`, `compute_ir_spectrum`, `compute_ir_spectrum_broadened`, `vibrational_analysis`, `ir_spectrum`, `compute_vibrational_analysis`, `compute_ir_spectrum`
- NMR: `predict_nmr_shifts`, `predict_nmr_couplings`, `compute_nmr_spectrum`, `compute_hose_codes`, `compute_ensemble_j_couplings`, `nmr_shifts`, `nmr_couplings`, `nmr_spectrum`, `hose_codes`, `predict_nmr_shifts`, `predict_nmr_couplings`, `compute_nmr_spectrum`, `compute_hose_codes`

Related theory docs:

- `/algorithm/uvvis-spectroscopy`
- `/algorithm/ir-spectroscopy`
- `/algorithm/nmr-spectroscopy`

Related SVGs:

- `/svg/uvvis-stda-pipeline.svg`
- `/svg/ir-hessian-pipeline.svg`
- `/svg/nmr-pipeline.svg`
- `/svg/nmr-karplus.svg`
- `/svg/broadening-functions.svg`
- `/svg/spectroscopy-overview.svg`

### Reactivity, pKa, and Site Ranking

Core idea:

- Build practical reactivity heuristics by combining frontier orbital content, Fukui-style condensed descriptors, and charge-based empirical rules.

Main API entries:

- Rust: `compute_frontier_descriptors`, `compute_fukui_descriptors`, `compute_reactivity_ranking`, `compute_empirical_pka`
- Python: `frontier_descriptors`, `fukui_descriptors`, `reactivity_ranking`, `empirical_pka`
- WASM: same `compute_...` names

Related theory docs:

- `/algorithm/population-analysis`
- `/algorithm/uvvis-spectroscopy` for transition-oriented electronic descriptors where relevant

Related SVGs:

- `/svg/population-mulliken.svg`
- `/svg/population-lowdin.svg`

### Alignment, RMSD, and Clustering

Core idea:

- Align structures optimally in Cartesian space, quantify RMSD, then cluster ensembles by geometric similarity.

Main API entries:

- Rust: `compute_rmsd`, `butina_cluster`, `compute_rmsd_matrix`, `embed_diverse`
- Python: `rmsd`, `butina_cluster`, `rmsd_matrix`, `filter_diverse`
- WASM: `compute_rmsd`, `butina_cluster`, `compute_rmsd_matrix`
- CLI: `rmsd`

Related theory docs:

- `/algorithm/molecular-alignment`

Related SVGs:

- `/svg/kabsch-alignment.svg`
- `/svg/benchmarks-rmsd.svg`

### Solvation

Core idea:

- Non-polar solvation uses SASA and atomic solvation parameters; GB adds an electrostatic implicit-solvent approximation.

Main API entries:

- Rust: `compute_nonpolar_solvation`, `compute_gb_solvation`
- Python: `nonpolar_solvation`, `gb_solvation`
- WASM: `compute_nonpolar_solvation`, `compute_gb_solvation`
- CLI: `solvation`

Related theory docs:

- `/algorithm/electrostatic-potential`

Related SVGs:

- `/svg/experimental-alpb.svg`

### Rings, Fingerprints, and Similarity

Core idea:

- Detect rings from graph structure, build circular fingerprints from local environments, and compare them by Tanimoto overlap.

Main API entries:

- Rust: `compute_sssr`, `compute_ecfp`, `compute_ecfp_batch`, `compute_tanimoto`
- Python: `sssr`, `ecfp`, `tanimoto`
- WASM: `compute_sssr`, `compute_ecfp`, `compute_tanimoto`
- CLI: `sssr`, `ecfp`, `tanimoto`

Related theory docs:

- ring and fingerprint behavior is currently summarized through API docs and examples more than a dedicated algorithm page

Related SVGs:

- `/svg/smarts-pattern.svg`

### Materials, Crystallography, and Framework Assembly

Core idea:

- Map cell parameters to lattice tensors and place node/linker building units on topology templates to assemble periodic frameworks.

Main API entries:

- Rust: `create_unit_cell`, `assemble_framework`
- Python: `unit_cell`, `assemble_framework`
- WASM: `create_unit_cell`, `assemble_framework`
- CLI: `cell`, `assemble`

Related theory docs:

- `/algorithm/materials-assembly`

Related SVGs:

- `/svg/materials-unit-cell.svg`
- `/svg/materials-assembly.svg`

### Dynamics, NEB, and Worker Transport

Core idea:

- Provide exploratory simulation and scalable transport helpers for long or distributed workflows.

Main API entries:

- Dynamics: `compute_md_trajectory`, `compute_md_trajectory_nvt`, `compute_simplified_neb_path`, `search_conformers_with_uff`
- Transport: `pack_conformers`, `pack_batch_arrow`, `split_worker_tasks`, `estimate_workers`

Related theory docs:

- `/algorithm/web-transport`

Related SVGs:

- `/svg/transport-arrow.svg`
- `/svg/transport-chunked.svg`

---

## Feature-Gated Experimental APIs

These are supported only when the corresponding feature flags are enabled.

### Experimental tracks and surfaces

| Capability | Rust feature | Python | WASM | CLI | SVG |
|-----------|--------------|--------|------|-----|-----|
| EEQ charges | `experimental-eeq` | `eeq_charges`, `eeq_energy` | `compute_eeq_charges`, `compute_eeq_energy` | `eeq` | `/svg/experimental-eeq.svg` |
| ALPB solvation | `experimental-alpb` | `alpb_solvation`, `alpb_born_radii` | `compute_alpb_solvation`, `compute_alpb_born_radii` | `alpb` | `/svg/experimental-alpb.svg` |
| D4 dispersion | `experimental-d4` | `d4_energy` | `compute_d4_energy` | `d4` | `/svg/experimental-d4.svg` |
| CPM electrochemistry | `experimental-cpm` | `cpm_charges`, `cpm_surface` | `compute_cpm_charges`, `compute_cpm_surface` | `cpm` | `/svg/experimental-cpm.svg` |
| CGA | `experimental-cga` | - | - | - | `/svg/experimental-cga.svg` |
| RandNLA | `experimental-randnla` | - | - | - | `/svg/experimental-randnla.svg` |
| Riemannian | `experimental-riemannian` | - | - | - | `/svg/experimental-riemannian.svg` |
| KPM | `experimental-kpm` | - | - | - | `/svg/experimental-kpm.svg` |
| SDR | `experimental-sdr` | - | - | - | `/svg/experimental-sdr.svg` |
| MBH | `experimental-mbh` | - | - | - | `/svg/experimental-mbh.svg` |
| GSM | `experimental-gsm` | - | - | - | `/svg/experimental-gsm.svg` |

Related theory docs:

- `/experimental/overview`
- `/experimental/eeq`
- `/experimental/alpb`
- `/experimental/d4`
- `/experimental/cpm`
- `/experimental/cga`
- `/experimental/randnla`
- `/experimental/riemannian`
- `/experimental/kpm`
- `/experimental/sdr`
- `/experimental/mbh`
- `/experimental/gsm`

---

## SVG Inventory, Diagram by Diagram

This inventory is intended to support the “SVG by SVG” documentation pass. Every listed asset exists today in `docs/public/svg/` unless otherwise noted.

### Conformer generation and validation

| SVG | Meaning |
|-----|---------|
| `pipeline-overview.svg` | Global end-to-end conformer pipeline |
| `etkdg-pipeline.svg` | ETKDG-specific stage breakdown |
| `bounds-matrix-overview.svg` | Bounds matrix concept and inputs |
| `bounds-matrix-1-4.svg` | 1-4 torsional distance constraints |
| `bounds-matrix-triangle.svg` | Triangle inequality smoothing intuition |
| `distance-geometry-pipeline.svg` | Distance geometry workflow |
| `distance-geometry-matrix.svg` | Metric matrix / embedding matrix intuition |
| `embedding-eigendecomp.svg` | Coordinate recovery from eigendecomposition |
| `etkdg-torsion.svg` | Torsion preference refinement |
| `validation-pipeline.svg` | Post-embedding validation stages |
| `validation-tetrahedral.svg` | Tetrahedral stereochemistry validation |
| `validation-double-bond.svg` | E/Z and double-bond geometry validation |

### Graphs, SMARTS, and rings

| SVG | Meaning |
|-----|---------|
| `smiles-pipeline.svg` | SMILES parsing workflow |
| `smiles-phenol.svg` | Example parsed structure |
| `smarts-pattern.svg` | SMARTS substructure pattern logic |

### EHT, population, dipole, ESP, DOS

| SVG | Meaning |
|-----|---------|
| `eht-pipeline.svg` | EHT computational pipeline |
| `population-mulliken.svg` | Mulliken partitioning concept |
| `population-lowdin.svg` | Löwdin orthogonalization / partitioning concept |
| `dipole-vector.svg` | Molecular dipole vector interpretation |
| `esp-grid.svg` | ESP grid construction and interpretation |
| `dos-smearing.svg` | Gaussian smearing for DOS |

### Force fields and alignment

| SVG | Meaning |
|-----|---------|
| `force-fields-overview.svg` | UFF/MMFF94 energy family overview |
| `force-fields-torsion.svg` | Torsional term behavior |
| `uff-energy-terms.svg` | UFF term decomposition |
| `mmff94-energy-terms.svg` | MMFF94 term decomposition |
| `kabsch-alignment.svg` | Optimal rigid alignment geometry |

### Spectroscopy

| SVG | Meaning |
|-----|---------|
| `spectroscopy-overview.svg` | High-level spectroscopy stack overview |
| `uvvis-stda-pipeline.svg` | UV-Vis / sTDA pipeline |
| `ir-hessian-pipeline.svg` | IR and vibrational Hessian workflow |
| `nmr-pipeline.svg` | NMR shift/spectrum workflow |
| `nmr-karplus.svg` | Karplus relation for J-couplings |
| `broadening-functions.svg` | Gaussian vs Lorentzian broadening |

### Materials and transport

| SVG | Meaning |
|-----|---------|
| `materials-unit-cell.svg` | Unit-cell geometry |
| `materials-assembly.svg` | Framework assembly logic |
| `transport-arrow.svg` | Arrow-style columnar transport |
| `transport-chunked.svg` | Worker/task chunking strategy |

### Benchmarks

| SVG | Meaning |
|-----|---------|
| `benchmarks-rmsd.svg` | RMSD benchmark comparison |
| `benchmarks-success.svg` | Success-rate benchmark comparison |
| `benchmarks-timing.svg` | Timing benchmark comparison |

### Experimental SVGs

| SVG | Meaning |
|-----|---------|
| `experimental-eeq.svg` | EEQ model overview |
| `experimental-alpb.svg` | ALPB implicit solvation overview |
| `experimental-d4.svg` | D4 dispersion workflow |
| `experimental-cpm.svg` | Constant-potential method workflow |
| `experimental-cga.svg` | CGA geometric algebra track |
| `experimental-kpm.svg` | Kernel Polynomial Method track |
| `experimental-mbh.svg` | Mobile Block Hessian track |
| `experimental-randnla.svg` | Randomized linear algebra track |
| `experimental-riemannian.svg` | Riemannian optimization track |
| `experimental-sdr.svg` | Semidefinite relaxation track |
| `experimental-gsm.svg` | Growing String Method track |

### Root-level public SVGs

| SVG | Meaning |
|-----|---------|
| `ani-architecture.svg` | ANI stack architecture |
| `hf3c-architecture.svg` | HF-3c stack architecture |
| `logo.svg` | Project/site branding |

---

## Practical Usage Patterns

### Recommended baseline workflow for users

1. `embed` a structure from SMILES
2. Reuse `elements` and `coords` for all geometry-dependent workflows
3. Choose the cheapest method that matches the property needed:
   - topology only: descriptors, ML properties, fingerprints, empirical pKa
   - fast electronic proxy: EHT-derived population/dipole/DOS/ESP
   - semiempirical SCF: PM3 or xTB
   - ab initio-like compact workflow: HF-3c
4. Move to spectroscopy, solvation, or materials only after geometry is validated

### Common mapping from user goal to API

| Goal | Recommended functions |
|------|-----------------------|
| Generate 3D coordinates | `embed`, `embed_batch`, `embed_diverse` |
| Get partial charges fast | `compute_charges` |
| Obtain quantum-like orbital descriptors cheaply | `compute_population`, `compute_dipole`, `compute_dos`, `compute_esp` |
| Estimate higher-fidelity electronic energies | `compute_pm3`, `compute_xtb`, `compute_hf3c` |
| Visualize orbitals | `compute_orbital_mesh`, `eht_orbital_mesh`, `eht_orbital_grid_typed` |
| Score conformers or geometry strain | `compute_uff_energy`, `compute_mmff94_energy` |
| Analyze stereochemistry | `analyze_stereo` / `stereo_analysis` / `stereo` |
| Compare molecular similarity | `compute_ecfp`, `compute_tanimoto` |
| Cluster conformers | `butina_cluster`, `compute_rmsd_matrix`, `filter_diverse` |

---

## Documentation Gaps Still Visible After This Consolidation

This file is comprehensive for the exposed surface, but the repo still has some areas that would benefit from dedicated deep-dive pages:

- PM3 theory and implementation notes
- xTB theory and implementation notes
- dedicated ring perception and fingerprint theory page
- dedicated dynamics / NEB theory page
- dedicated topology-analysis page for coordination environments
- dedicated orbital-mesh page for cross-method visualization

Those gaps are now easier to fill because the public surface and the diagram inventory are consolidated here.