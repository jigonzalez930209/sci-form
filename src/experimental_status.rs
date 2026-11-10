//! E4 – Experimental Module Validation & Retention Criteria
//!
//! This module documents which experimental methods remain gated behind
//! feature flags and the criteria required for eventual promotion to core.
//!
//! # This is the authoritative "is this feature ready?" table
//!
//! When you ask whether a subsystem of this crate is production-ready, the
//! answer is the table in this file — not the module's own doc comment and
//! not the feature-flag list in `Cargo.toml`. Every experimental subsystem
//! appears exactly once, with its **feature flag**, its **location**, the
//! **validation gate** that covers it and its **promotion criteria**. A
//! subsystem moves out of the "Retained as Experimental" section only when
//! its listed criteria are met and the corresponding validation gate passes;
//! the "Promoted to Core" table above records the subsystems that already
//! cleared that bar.
//!
//! # Validation gates
//!
//! The gates are the cross-library *fichas* and integration tests under
//! `tests/` — one per subsystem, named in the tables below. The culture they
//! follow is documented in `reports/validation/cross-library-validation.es.md`
//! §13–§14: a module is promoted when a cross-library comparison (or a
//! documented external benchmark) gates its observable, plus integration
//! tests through the public API.
//!
//! ## Promoted to Core
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
//! | k-mesh/Bloch | `kmesh` | ASE meshes + chain bands 4.4e-16 (F88); periodic KPM/RandNLA adapters live in core `periodic_linear` |
//! | Chart payloads/builders | `render_bridge` | JSON Schema + pack round-trip (F90); EDL/kinetics/periodic/GSM domain builders are core |
//! | Riemannian PSD | `riemannian` | Projection/optimizer vs scipy (F94; Frobenius proxy documented); tetrahedron Gram recovery 6.0e-9 (`tests/beta_benchmarks.rs`) |
//! | CPM | `cpm` | Constant-potential electrochemistry vs LAMMPS `electrode/conp` (F95); two-site closed form 1.3e-15 (`tests/beta_benchmarks.rs`) |
//! | ReaxFF | `forcefield::reaxff` | Bond orders vs LAMMPS (F78) + analytical gradients (2026-11-08) |
//! | MMFF94 parameters | `forcefield::mmff94_params` | Real published parameter tables (F100) |
//! | KPM | `kpm` | Chebyshev DOS/moments; exact moments vs numpy, sparse CSR 5000 sites in 0.22 s, E_F <0.1 eV vs diagonalization (F91) |
//! | RandNLA | `rand_nla` | F93 3.8e-4 at rank-N/4; EHT backend; N=600 low-rank path at 9.6e-4 in ~80 ms vs ~143 ms |
//! | GSM | `gsm` | F82 vs ASE CI-NEB (1.6e-3 barrier, 0.028 TS) + five curved wells ≤0.8 % |
//! | Kinetics / HTST | `kinetics` | F87 Eyring 1.9e-16, traces vs scipy `solve_ivp`, equilibrium 8e-16 |
//! | Periodic-linear adapters | `periodic_linear` | F88 k-mesh vs ASE, phases 3.6e-16, bands 4.4e-16; KPM/RandNLA adapters tested |
//! | Live MD | `dynamics_live` | F84 forces vs RDKit UFF 3.6e-13, NVE drift 4.4e-4, reversibility 5.7e-17 |
//!
//! The promotion criteria below follow the validation culture of the repo:
//! a module is promoted when a cross-library ficha (or a documented external
//! benchmark) gates its observable, plus integration tests through the public
//! API. See `reports/validation/cross-library-validation.es.md` §13–§14.
//!
//! **This module is documentation-only; it contains no runtime code.** The
//! experimental modules themselves remain behind their respective feature
//! flags, and the table above is the single source of truth for whether a
//! given subsystem has been promoted.
//!
//! ---
//!
//! ## Module ledger — promoted entries and retained research
//!
//! ### EXP-201: Kernel Polynomial Method (KPM) — PROMOTED (`src/kpm`)
//!
//! **Feature flag:** `beta-kpm` is retained as a compatibility alias; the
//! module is now unconditional.
//! **Location:** `src/kpm/` (chebyshev.rs, density.rs, sparse.rs)
//! **Validation test:** `tests/beta_kpm_validation.rs` (F91)
//!
//! **Promotion evidence (2026-11-07):** sparse CSR backend; exact moments
//! 3.5e-15 and DOS Jackson L2 1.7e-14 vs numpy eigh; N=200/500 E_F <0.1 eV;
//! N=5000 chain in 0.22 s with 32 stochastic vectors and 400 terms. The old
//! `sci_form::beta::kpm` path is an alias.
//!
//! ### EXP-202: Multi-Basin Hopping (MBH)
//!
//! **Feature flag:** `beta-mbh`
//! **Location:** `src/beta/mbh/` (blocks.rs, hessian.rs)
//! **Validation test:** `tests/beta_mbh_validation.rs` (F92)
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
//! **External validation (2026-10-14, F92):** diatomic stretches against the
//! SI closed form 1.3e-5 (validates the 108.59 factor) and block structure in
//! benzene (`tests/beta_mbh_validation.rs`). Gap: ±2 cm⁻¹ experimental.
//!
//! **Benchmark against the stable path (2026-10-23):** MBH frequencies vs the
//! projected analytic Hessian in `tests/beta_benchmarks.rs`, maximum deviation
//! 2.6 % (near-null modes are dominated by finite-difference noise, ~5 cm⁻¹).
//!
//! ### EXP-203: Randomized Numerical Linear Algebra (RandNLA) — PROMOTED (`src/rand_nla`)
//!
//! **Feature flag:** `beta-randnla` is retained as a compatibility alias; the
//! module is now unconditional.
//! **Location:** `src/rand_nla/` (nystrom.rs, solver.rs)
//! **Validation test:** `tests/beta_randnla_validation.rs` (F93)
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
//! **External validation (2026-10-14, F93):** single-shot characterization
//! rank-N/4 ~0.13, valid regime 3.9e-3, exact production with fallback
//! (`tests/beta_randnla_validation.rs`).
//!
//! **Subspace iterations (2026-11-07):** `NystromApprox` gains
//! `from_matrix_with_power_iterations` and `RandNlaConfig::power_iterations`
//! (default 2; re-orthonormalizes each iteration). With 384 iterations the
//! rank N/4 of F93 drops to 3.8e-4, so **criterion 2 (<1e-3) is met** for the
//! fixture.
//!
//! **EHT backend (2026-11-07):** `eht::solve_eht_randnla` exposes the
//! randomized solver as an optional backend (feature `beta-randnla`), with the
//! same EHT geometry/basis and residual/fallback metadata. Criterion 1 is
//! covered.
//!
//! **N>500 path (2026-11-08):** `solve_eht_randnla_lowrank` avoids forming and
//! diagonalizing the transformed N×N Hamiltonian; it uses Nyström for
//! `S^{-1/2}` and Lanczos with full reorthogonalization from a random start.
//! At N=600 (chain, k=160) the four lowest states have 9.6e-4 relative error
//! and the time drops to ~80 ms versus ~143 ms for exact diagonalization.
//! **Criterion 3 (speedup N>500) is met** for the lowest state; higher
//! excitations still need the exact path or a larger block.
//!
//! **Benchmark against the stable path (2026-10-23):** Nyström vs exact SVD in
//! `tests/beta_benchmarks.rs`; the reconstruction error is ≤ 2× the best
//! rank-k and decreases as the rank grows.
//!
//! ### EXP-204: Spectral Dimensionality Reduction (SDR) — PROMOTED (`src/sdr`)
//!
//! **Feature flag:** `alpha-sdr`
//! **Location:** `src/sdr/` (embedding.rs, projections.rs)
//! **Validation test:** `tests/alpha_sdr_validation.rs` (F83)
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
//! **External validation (2026-10-14, F83):** exact projections vs numpy,
//! exact embeds 7e-15 (better than CVXPY/SCS 1e-8), noisy case within
//! 3× of the optimum (`tests/alpha_sdr_validation.rs`).
//!
//! **Benchmark against the stable path (2026-10-23):** SDR vs classical MDS on
//! exact distances in `tests/beta_benchmarks.rs`; maximum pair error
//! 6.7e-16 vs 2.2e-16 for MDS, without exceeding the oracle.
//!
//! **SDR vs Butina (2026-11-07):** criterion 3 is measured with the Adjusted
//! Rand Index between the k-means of the SDR embedding (on the Kabsch RMSD
//! matrix) and the Butina labels on 30 conformers of three different shapes:
//! ARI = 1.0, also against the synthetic truth (`tests/beta_benchmarks.rs`).
//! The adaptive kernel bandwidth (criterion 2) is still open.
//!
//! ### EXP-205a: Conformal Geometric Algebra (CGA) — PROMOTED (`src/cga`)
//!
//! **Feature flag:** `alpha-cga`
//! **Location:** `src/cga/` (multivector.rs, motor.rs, conformer.rs, materials.rs)
//! **Validation test:** `tests/alpha_cga_validation.rs` (F81)
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
//! **External validation (2026-10-14, F81):** motors 3.4e-15 vs Rodrigues and
//! `clifford` Cl(4,1), dihedrals 6.6e-15 (`tests/alpha_cga_validation.rs`).
//! Mirror-angle convention documented.
//!
//! ### EXP-205b: Growing String Method (GSM) — PROMOTED (`src/gsm`)
//!
//! **Feature flag:** `alpha-gsm` is retained as a compatibility alias; the
//! module is now unconditional.
//! **Location:** `src/gsm/` (string.rs, saddle.rs, backend.rs)
//! **Validation test:** `tests/alpha_gsm_validation.rs` (F82)
//!
//! Reaction path finding via growing string interpolation between
//! reactant and product geometries. Identifies transition states.
//!
//! **Current limitations:**
//! - No analytic Hessian; the saddle refinement is gradient-only
//! - String convergence criteria incomplete
//! - Saddle point refinement was under-stepped (fixed 20 steps at 0.1×
//!   `step_size`; now 200 steps at `step_size`, which brings the curved-well
//!   barriers from ~12 % high to ≤0.8 %)
//!
//! **Status update:** gradient coupling is implemented — `gsm::backend`
//! drives the string with UFF, MMFF94, PM3, xTB, GFN1, GFN2 and HF-3c.
//!
//! **Promotion criteria:**
//! 1. ~~Couple with PM3 or xTB gradient evaluation~~ (done: multi-backend)
//! 2. Reproduce known barrier heights (SN2, Diels-Alder)
//! 3. Compare path energy profiles vs NEB reference in 5+ reactions
//!
//! **External validation (2026-10-14, F82):** TS on Müller-Brown vs ASE
//! CI-NEB (barrier 1.6e-3, TS 0.028; `tests/alpha_gsm_validation.rs`).
//! Fix: arc-length reparametrization every iteration.
//!
//! **Benchmark against the stable path (2026-10-23):** five analytic double
//! wells `V(x) = (x² − a²)²` in `tests/beta_benchmarks.rs`; the barrier and the
//! TS position match the exact solution (barrier relative error 0).
//!
//! **GSM vs NEB on 5+ reactions (2026-11-07):** `tests/alpha_gsm_validation.rs`
//! adds five curved wells `V(x,y) = (x²−a²)² + ½(y−c x²)²`, whose MEP is the
//! parabola `y = c x²` and whose exact barrier is `a⁴`. GSM and a discrete NEB
//! (reference in the test) match the analytic barrier to ≤0.8 % and each
//! other to ≤0.8 %; criterion 3 is covered for five reactions.
//!
//! ### EXP-205c: Constant Potential Method (CPM) — PROMOTED (`src/cpm`)
//!
//! **Feature flag:** `beta-cpm` (alias; the core module is unconditional)
//! **Location:** `src/cpm/` (surface.rs, grand_potential.rs)
//! **Validation test:** `tests/experimental/test_cpm.rs` (F95)
//!
//! Constant-potential electrochemistry: the electrode potential is a
//! degree of freedom and the atomic charges follow from the grand potential.
//!
//! **Current limitations:**
//! - The Coulomb kernel has no Gaussian smearing (LAMMPS `electrode/conp`
//!   has), so the charge response is ~1.76x too weak
//! - The self-interaction hardness is 5 eV vs ~22 eV for a Gaussian
//! - No integration with a plane-wave / DFT reference
//!
//! **Promotion criteria:**
//! 1. Add the Gaussian-smeared Coulomb kernel and match LAMMPS
//!    `electrode/conp` on the charge pattern
//! 2. Validate the quantum capacitance curve against experiment or DFT
//! 3. Integrate with `materials::UnitCell` for electrode geometry
//!
//! **Note 2026-10-13:** CPM promoted to core (`src/cpm`, F95 vs LAMMPS
//! `electrode/conp`); `beta::cpm` survives as a compatible re-export.
//!
//! **Benchmark against the stable path (2026-10-23):** two-site charges vs the
//! closed-form solution of the CPM linear system in `tests/beta_benchmarks.rs`
//! (1.3e-15) and stationarity residual < 1e-10.
//!
//! ### EXP-206: Remaining gated subsystems
//!
//! These remain gated and are NOT promoted:
//!
//! | Subsystem | Feature flag | Location | Validation test | Notes |
//! |---|---|---|---|---|
//! | Kohn–Sham DFT | `alpha-dft` | `src/dft/` | `tests/alpha_dft_pyscf.rs` (F77) | E ≤ 1 mHa vs PySCF; `coulomb_energy` is still 0.0 and the ECP path is unfitted |
//! | ML force field | `alpha-mlff` | `src/ml/{mlff,symmetry_functions}.rs` | `tests/alpha_mlff_validation.rs` (F79) | AEV vs a numpy BP oracle; no trained weights, no MD17/SPICE benchmark |
//! | Reaction dynamics | `alpha-reaction-dynamics` | `src/alpha/reaction_dynamics/` | `tests/alpha_reaction_dynamics_3d.rs` | self-consistency only; no ASE NEB parity |
//!
//! ### EXP-207: GPU acceleration
//!
//! **Feature flag:** `experimental-gpu` (plus `experimental-gpu-rendering`)
//! **Location:** `src/gpu/`
//! **Validation test:** `tests/regression/test_gpu_candidates.rs`
//!
//! The GPU tree is a set of compute backends with CPU fallbacks; it is not a
//! "scientific method" and has no promotion criteria of its own. Two
//! properties are enforced instead:
//!
//! 1. The CPU fallback must always work. `GpuContext::try_create` wraps the
//!    adapter acquisition in `catch_unwind`, because wgpu-hal's GLES/EGL
//!    backend panics with `BadDisplay` on headless hosts.
//! 2. Every kernel that claims a GPU result must be compared against the
//!    CPU implementation. The MMFF94 kernel uses a workgroup reduction plus
//!    signed fixed-point atomic adds (WGSL has no atomic float add).
//!
//! Known gaps: the GPU results are only ever compared against sci-form's own
//! CPU code, never against an external reference (libcint for the ERI/S/T/V
//! kernels, dftd4 for the dispersion kernel, VTK for marching cubes).
//!
//! The GPU tree also honours two hard WebGPU limits that callers may rely
//! on: **256 invocations per workgroup**
//! (`max_compute_invocations_per_workgroup`) and **128 MB per storage-buffer
//! binding** (`max_storage_buffer_binding_size`).

// This module is documentation-only; no runtime code.
// The experimental modules remain behind their respective feature flags.
