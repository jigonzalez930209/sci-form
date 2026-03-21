//! E4 – Experimental Module Validation & Retention Criteria
//!
//! This module documents which experimental methods remain gated behind
//! feature flags and the criteria required for eventual promotion to core.
//!
//! ## Promoted to Core (E3)
//!
//! | Module | Core path | Justification |
//! |--------|-----------|---------------|
//! | EEQ charges | `charges_eeq` | Topology-free partial charges, O(N²), validated vs Gasteiger |
//! | D4 dispersion | `dispersion` | BJ-damped DFT-D4 energy, complements UFF/MMFF94 |
//! | ALPB solvation | `solvation_alpb` | Klamt-corrected solvation, improves over plain GB |
//! | sTDA UV-Vis | `spectroscopy::stda_uvvis` | CIS-like excitations from monopole charges, ~1000× faster than full CI |
//! | GIAO NMR | `spectroscopy::giao_nmr` | Multi-nucleus NMR shieldings, complements empirical HOSE |
//!
//! ---
//!
//! ## Retained as Experimental (E4)
//!
//! ### EXP-201: Kernel Polynomial Method (KPM)
//!
//! **Feature flag:** `experimental-kpm`
//! **Location:** `experimental/kpm/` (chebyshev.rs, density.rs)
//!
//! KPM approximates the density of states via Chebyshev expansion of the
//! spectral function. Useful for very large systems (10k+ atoms) but the
//! current implementation lacks:
//!
//! - Adaptive polynomial order selection
//! - Sparse-matrix backend for O(N) scaling
//! - Benchmarks against direct diagonalization for <500 atoms
//!
//! **Promotion criteria:**
//! 1. Demonstrate wall-clock advantage over direct DOS for >2000 atoms
//! 2. Validate spectral accuracy within 0.1 eV of full diagonalization
//! 3. Add sparse matrix support (nalgebra-sparse or sprs)
//!
//! ### EXP-202: Multi-Basin Hopping (MBH)
//!
//! **Feature flag:** `experimental-mbh`
//! **Location:** `experimental/mbh/` (blocks.rs, hessian.rs)
//!
//! Block-diagonal Hessian approach for exploring conformational basins.
//! Intended for global optimization of flexible molecules.
//!
//! **Current limitations:**
//! - No energy evaluator integration (needs UFF/MMFF94 coupling)
//! - Block decomposition heuristic is molecule-size dependent
//! - Missing Metropolis acceptance criterion
//!
//! **Promotion criteria:**
//! 1. Integrate with `forcefield` module for energy evaluation
//! 2. Reproduce known global minima for test molecules (alanine dipeptide, etc.)
//! 3. Compare basin-hopping yield vs multi-seed ETKDG
//!
//! ### EXP-203: Randomized Numerical Linear Algebra (RandNLA)
//!
//! **Feature flag:** `experimental-randnla`
//! **Location:** `experimental/rand_nla/` (nystrom.rs, solver.rs)
//!
//! Nyström approximation and randomized eigensolvers for large overlap/Fock
//! matrices. Currently a backend alternative, not a user-facing feature.
//!
//! **Current limitations:**
//! - Rank selection requires manual tuning
//! - No automatic accuracy estimation
//! - Not yet wired into SCF or EHT solvers
//!
//! **Promotion criteria:**
//! 1. Wire into EHT/SCF as optional backend with feature flag
//! 2. Show <1% error in eigenvalues for rank = N/4
//! 3. Demonstrate speedup for N>500 basis functions
//!
//! ### EXP-204: Spectral Dimensionality Reduction (SDR)
//!
//! **Feature flag:** `experimental-sdr`
//! **Location:** `experimental/sdr/` (embedding.rs, projections.rs)
//!
//! Diffusion-map and spectral embedding for chemical space visualization.
//! Alternative to fingerprint-based Tanimoto similarity.
//!
//! **Current limitations:**
//! - Kernel bandwidth selection is fixed
//! - No out-of-sample extension (Nyström embedding)
//! - Limited to datasets fitting in memory
//!
//! **Promotion criteria:**
//! 1. Validate embeddings preserve known SAR relationships
//! 2. Add adaptive bandwidth selection
//! 3. Compare clustering quality vs Butina (Tanimoto)
//!
//! ### EXP-205a: Conformal Geometric Algebra (CGA)
//!
//! **Feature flag:** `experimental-cga`
//! **Location:** `experimental/cga/` (multivector.rs, motor.rs, conformer.rs, materials.rs)
//!
//! CGA representation of molecular geometry: motors for rigid transforms,
//! conformal points for distance geometry. Research-stage alternative to
//! Cartesian coordinate manipulation.
//!
//! **Current limitations:**
//! - 5D multivector operations have overhead vs 3D rotations
//! - Conformer pipeline not competitive with ETKDG
//! - Materials assembly untested with real MOF structures
//!
//! **Promotion criteria:**
//! 1. Demonstrate advantage in distance geometry constraint handling
//! 2. Benchmark CGA conformer vs ETKDG quality and speed
//! 3. Show utility for symmetry-constrained crystal assembly
//!
//! ### EXP-205b: Growing String Method (GSM)
//!
//! **Feature flag:** `experimental-gsm`
//! **Location:** `experimental/gsm/` (string.rs, saddle.rs)
//!
//! Reaction path finding via growing string interpolation between
//! reactant and product geometries. Identifies transition states.
//!
//! **Current limitations:**
//! - No gradient source (needs coupled PES evaluation)
//! - String convergence criteria incomplete
//! - Saddle point refinement not validated
//!
//! **Promotion criteria:**
//! 1. Couple with PM3 or xTB gradient evaluation
//! 2. Reproduce known barrier heights (SN2, Diels-Alder)
//! 3. Compare path energy profiles vs NEB reference
//!
//! ### EXP-205c: Continuous Phase Methods (CPM)
//!
//! **Feature flag:** `experimental-cpm`
//! **Location:** `experimental/cpm/` (surface.rs, grand_potential.rs)
//!
//! Phase-field and grand potential methods for crystal surface
//! thermodynamics and phase boundary determination.
//!
//! **Current limitations:**
//! - Grand potential formulation needs thermodynamic data
//! - Surface energy calculations untested for ternary systems
//! - No integration with `materials` module
//!
//! **Promotion criteria:**
//! 1. Validate surface energies against DFT reference for binary alloys
//! 2. Integrate with `materials::UnitCell` for consistent crystal handling
//! 3. Add Wulff construction for equilibrium crystal shapes

// This module is documentation-only; no runtime code.
// The experimental modules remain behind their respective feature flags.
