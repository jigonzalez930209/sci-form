//! E4 – Experimental Module Validation & Retention Criteria
//!
//! This module documents which experimental methods remain gated behind
//! feature flags and the criteria required for eventual promotion to core.
//!
//! ## Promoted to Core (E3)
//!
//! | Module | Core path | Justification |
//! |--------|-----------|---------------|
//! | EEQ charges | `charges_eeq` | Topology-free partial charges, O(N²); EEQ2019 ported for D4 and validated vs dftd4 (F46) |
//! | D4 dispersion | `dispersion` | BJ-damped DFT-D4 energy, complements UFF/MMFF94; validated vs dftd4 (1.2e-10 Ha, F46) |
//! | ALPB solvation | `solvation_alpb` | Klamt-corrected solvation, improves over plain GB; shares the SASA implementation with the core (test) |
//! | sTDA UV-Vis | `spectroscopy::stda_uvvis` | CIS-like excitations from monopole charges, ~1000× faster than full CI; Cauchy-Schwarz prescreening tested |
//! | GIAO NMR | `spectroscopy::giao_nmr` | Multi-nucleus NMR shieldings, complements empirical HOSE |
//! | Obara-Saika ERIs | `scf::{obara_saika,screening}` | Boys/ERI/Schwarz vs scipy + libcint, 5.6e-14 (F80) |
//! | CGA motors | `cga` | Rotors vs Rodrigues + clifford 3.4e-15 (F81; refine/materials gaps documented) |
//! | SDR embedding | `sdr` | Projections + CVXPY embeds 7e-15 (F83) |
//! | IMD steering | `imd` | Harmonic E/F analytic 1.1e-13 (F85) |
//! | EDL | `edl` | GC/GCS/Helmholtz analytic 2e-14 (F86) |
//! | k-mesh/Bloch | `kmesh` | ASE meshes + chain bands 4.4e-16 (F88; KPM/RandNLA adapters stay alpha) |
//! | Chart payloads | `render_bridge` | JSON Schema + pack round-trip (F90; domain builders stay alpha) |
//! | Riemannian PSD | `riemannian` | Projection/optimizer vs scipy (F94; Frobenius proxy documented) |
//!
//! The promotion criteria below follow the validation culture of the repo:
//! a module is promoted when a cross-library ficha (or a documented external
//! benchmark) gates its observable, plus integration tests through the public
//! API. See `docs/validation/cross-library-validation.es.md` §13–§14.
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
//! - Benchmarks against direct diagonalization for <500 atoms (todo: gate
//!   spectral accuracy vs exact diagonalization for N > 100)
//!
//! **Promotion criteria:**
//! 1. Demonstrate wall-clock advantage over direct DOS for >2000 atoms
//! 2. Validate spectral accuracy within 0.1 eV of full diagonalization
//! 3. Add sparse matrix support (nalgebra-sparse or sprs)
//!
//! **Validación externa (2026-10-14, F91)**: momentos exactos 3.5e-15,
//! DOS Jackson L2 1.7e-14 y Mulliken (gap 3.8e-5 / finite-T 2.6e-2) vs
//! numpy eigh (`tests/beta_kpm_validation.rs`). Fixes aplicados: factor 2
//! en coeficientes FD y degeneración de espín 2 en Mulliken.
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
//! **Validación externa (2026-10-14, F92)**: stretches diatómicos vs forma
//! cerrada SI 1.3e-5 (valida el factor 108.59) y estructura de bloques en
//! benceno (`tests/beta_mbh_validation.rs`). Gap: ±2 cm⁻¹ experimental.
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
//! - Not yet wired into SCF or EHT solvers (todo: scale benchmarks)
//!
//! **Promotion criteria:**
//! 1. Wire into EHT/SCF as optional backend with feature flag
//! 2. Show <1% error in eigenvalues for rank = N/4
//! 3. Demonstrate speedup for N>500 basis functions
//!
//! **Validación externa (2026-10-14, F93)**: caracterización rank-N/4 ~0.13
//! (gap cuantificado: 1e-3 exige power iterations), régimen válido 3.9e-3,
//! producción exacta con fallback (`tests/beta_randnla_validation.rs`).
//! El criterio 2 sigue sin cumplirse.
//!
//! ### EXP-204: Spectral Dimensionality Reduction (SDR) — PROMOVIDO 2026-10-14 (`src/sdr`)
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
//! - Limited to datasets fitting in memory (todo: 10K-molecule benchmark)
//!
//! **Promotion criteria:**
//! 1. Validate embeddings preserve known SAR relationships
//! 2. Add adaptive bandwidth selection
//! 3. Compare clustering quality vs Butina (Tanimoto)
//!
//! **Validación externa (2026-10-14, F83)**: proyecciones exactas vs numpy,
//! embeds exactos 7e-15 (mejor que CVXPY/SCS 1e-8), caso ruidoso dentro de
//! 3× del óptimo (`tests/alpha_sdr_validation.rs`).
//!
//! ### EXP-205a: Conformal Geometric Algebra (CGA) — PROMOVIDO 2026-10-14 (`src/cga`)
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
//! **Validación externa (2026-10-14, F81)**: motores 3.4e-15 vs Rodrigues y
//! `clifford` Cl(4,1), diedros 6.6e-15 (`tests/alpha_cga_validation.rs`).
//! Convención de ángulo espejo documentada.
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
//! **Status update:** gradient coupling is implemented — `alpha::gsm::backend`
//! drives the string with UFF, MMFF94, PM3, xTB, GFN1, GFN2 and HF-3c.
//! **Promotion criteria:**
//! 1. ~~Couple with PM3 or xTB gradient evaluation~~ (done: multi-backend)
//! 2. Reproduce known barrier heights (SN2, Diels-Alder)
//! 3. Compare path energy profiles vs NEB reference in 5+ reactions
//!
//! **Validación externa (2026-10-14, F82)**: TS en Müller-Brown vs ASE
//! CI-NEB (barrera 1.6e-3, TS 0.028; `tests/alpha_gsm_validation.rs`).
//! Fix: reparametrización por arco cada iteración.
//!
//! ### EXP-205c: Continuous Phase Methods (CPM) — PROMOVIDO 2026-10-13 (`src/cpm`)
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
//!
//! **Nota 2026-10-13**: CPM promovido a core (`src/cpm`, F95 vs LAMMPS
//! `electrode/conp`); `beta::cpm` sobrevive como re-export compatible.

// This module is documentation-only; no runtime code.
// The experimental modules remain behind their respective feature flags.
