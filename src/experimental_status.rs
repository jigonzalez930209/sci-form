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
//! follow is documented in `docs/validation/cross-library-validation.es.md`
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
//! | k-mesh/Bloch | `kmesh` | ASE meshes + chain bands 4.4e-16 (F88; KPM/RandNLA adapters stay alpha) |
//! | Chart payloads | `render_bridge` | JSON Schema + pack round-trip (F90; domain builders stay alpha) |
//! | Riemannian PSD | `riemannian` | Projection/optimizer vs scipy (F94; Frobenius proxy documented); tetrahedron Gram recovery 6.0e-9 (`tests/beta_benchmarks.rs`) |
//! | CPM | `cpm` | Constant-potential electrochemistry vs LAMMPS `electrode/conp` (F95); two-site closed form 1.3e-15 (`tests/beta_benchmarks.rs`) |
//! | MMFF94 parameters | `forcefield::mmff94_params` | Real published parameter tables (F100) |
//!
//! The promotion criteria below follow the validation culture of the repo:
//! a module is promoted when a cross-library ficha (or a documented external
//! benchmark) gates its observable, plus integration tests through the public
//! API. See `docs/validation/cross-library-validation.es.md` §13–§14.
//!
//! **This module is documentation-only; it contains no runtime code.** The
//! experimental modules themselves remain behind their respective feature
//! flags, and the table above is the single source of truth for whether a
//! given subsystem has been promoted.
//!
//! ---
//!
//! ## Retained as Experimental
//!
//! ### EXP-201: Kernel Polynomial Method (KPM)
//!
//! **Feature flag:** `beta-kpm`
//! **Location:** `src/beta/kpm/` (chebyshev.rs, density.rs)
//! **Validation test:** `tests/beta_kpm_validation.rs` (F91)
//!
//! KPM approximates the density of states via Chebyshev expansion of the
//! spectral function. Useful for very large systems (10k+ atoms).
//!
//! **Promotion criteria:**
//! 1. Demonstrate wall-clock advantage over direct DOS for >2000 atoms
//! 2. Validate spectral accuracy within 0.1 eV of full diagonalization
//! 3. Add sparse matrix support (nalgebra-sparse or sprs)
//!
//! **Estado 2026-11-07:** los tres criterios tienen evidencia medida:
//! (1) la cadena dispersa de 5000 sitios tarda 0.22 s con 32 vectores
//! estocásticos y 400 términos, mientras que la diagonalización densa no es
//! viable a esa escala; (2) F91d compara N=200/500 con `eigh` y deja la
//! energía de Fermi a <0.1 eV; (3) `SparseMatrix` (CSR propio) +
//! `compute_kpm_dos_sparse` cubren el backend disperso. Queda la decisión de
//! API para promoverlo fuera de `beta-kpm`.
//!
//! **Validación externa (2026-10-14, F91)**: momentos exactos 3.5e-15,
//! DOS Jackson L2 1.7e-14 y Mulliken (gap 3.8e-5 / finite-T 2.6e-2) vs
//! numpy eigh (`tests/beta_kpm_validation.rs`). Fixes aplicados: factor 2
//! en coeficientes FD y degeneración de espín 2 en Mulliken.
//!
//! **Benchmark vs ruta estable (2026-10-23)**: `tests/beta_benchmarks.rs`
//! compara la DOS KPM con diagonalización exacta (DOS integrada 40.17/40,
//! corrimiento del centro 1e-4, correlación 0.86 con el kernel Jackson).
//! **Escala (2026-11-07)**: `tests/beta_kpm_validation.rs` añade cadenas
//! desordenadas de N=200 y N=500; la DOS integrada coincide a <1e-6 relativo
//! y la energía de Fermi a <0.1 eV. Los criterios 1 (ventaja >2000 átomos) y
//! 3 (matriz dispersa) siguen abiertos.
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
//! **Validación externa (2026-10-14, F92)**: stretches diatómicos vs forma
//! cerrada SI 1.3e-5 (valida el factor 108.59) y estructura de bloques en
//! benceno (`tests/beta_mbh_validation.rs`). Gap: ±2 cm⁻¹ experimental.
//!
//! **Benchmark vs ruta estable (2026-10-23)**: frecuencias MBH vs Hessiano
//! analítico proyectado en `tests/beta_benchmarks.rs`, desviación máxima
//! 2.6 % (los modos casi nulos quedan dominados por el ruido de diferencias
//! finitas, ~5 cm⁻¹).
//!
//! ### EXP-203: Randomized Numerical Linear Algebra (RandNLA)
//!
//! **Feature flag:** `beta-randnla`
//! **Location:** `src/beta/rand_nla/` (nystrom.rs, solver.rs)
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
//! **Validación externa (2026-10-14, F93)**: caracterización single-shot
//! rank-N/4 ~0.13, régimen válido 3.9e-3, producción exacta con fallback
//! (`tests/beta_randnla_validation.rs`).
//!
//! **Subspace iterations (2026-11-07)**: `NystromApprox` gana
//! `from_matrix_with_power_iterations` y `RandNlaConfig::power_iterations`
//! (default 2; re-ortogonaliza cada iteración). Con 384 iteraciones el rank
//! N/4 de F93 baja a 3.8e-4, así que el **criterio 2 (<1e-3) queda cumplido**
//! para el fixture.
//!
//! **Backend EHT (2026-11-07)**: `eht::solve_eht_randnla` expone el solver
//! aleatorizado como backend opcional (feature `beta-randnla`), con la misma
//! geometría/basis EHT y metadatos de residuo/fallback. El criterio 1 queda
//! cubierto; el **criterio 3 (speedup N>500)** sigue abierto.
//!
//! **Benchmark vs ruta estable (2026-10-23)**: Nyström vs SVD exacta en
//! `tests/beta_benchmarks.rs`; el error de reconstrucción queda ≤ 2× el mejor
//! rank-k y decrece al aumentar el rango.
//!
//! ### EXP-204: Spectral Dimensionality Reduction (SDR) — PROMOVIDO (`src/sdr`)
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
//! **Validación externa (2026-10-14, F83)**: proyecciones exactas vs numpy,
//! embeds exactos 7e-15 (mejor que CVXPY/SCS 1e-8), caso ruidoso dentro de
//! 3× del óptimo (`tests/alpha_sdr_validation.rs`).
//!
//! **Benchmark vs ruta estable (2026-10-23)**: SDR vs MDS clásico sobre
//! distancias exactas en `tests/beta_benchmarks.rs`; error máximo de par
//! 6.7e-16 frente a 2.2e-16 del MDS, sin quedar por encima del oráculo.
//!
//! **SDR vs Butina (2026-11-07)**: el criterio 3 se mide con el Adjusted Rand
//! Index entre el k-means del embedding SDR (sobre la matriz de RMSD Kabsch)
//! y las etiquetas de Butina en 30 conformeros de tres formas distintas:
//! ARI = 1.0, también contra la verdad sintética (`tests/beta_benchmarks.rs`).
//! La banda de kernel adaptativa (criterio 2) sigue abierta.
//!
//! ### EXP-205a: Conformal Geometric Algebra (CGA) — PROMOVIDO (`src/cga`)
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
//! **Validación externa (2026-10-14, F81)**: motores 3.4e-15 vs Rodrigues y
//! `clifford` Cl(4,1), diedros 6.6e-15 (`tests/alpha_cga_validation.rs`).
//! Convención de ángulo espejo documentada.
//!
//! ### EXP-205b: Growing String Method (GSM)
//!
//! **Feature flag:** `alpha-gsm`
//! **Location:** `src/alpha/gsm/` (string.rs, saddle.rs, backend.rs)
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
//! **Status update:** gradient coupling is implemented — `alpha::gsm::backend`
//! drives the string with UFF, MMFF94, PM3, xTB, GFN1, GFN2 and HF-3c.
//!
//! **Promotion criteria:**
//! 1. ~~Couple with PM3 or xTB gradient evaluation~~ (done: multi-backend)
//! 2. Reproduce known barrier heights (SN2, Diels-Alder)
//! 3. Compare path energy profiles vs NEB reference in 5+ reactions
//!
//! **Validación externa (2026-10-14, F82)**: TS en Müller-Brown vs ASE
//! CI-NEB (barrera 1.6e-3, TS 0.028; `tests/alpha_gsm_validation.rs`).
//! Fix: reparametrización por arco cada iteración.
//!
//! **Benchmark vs ruta estable (2026-10-23)**: cinco pozos dobles analíticos
//! `V(x) = (x² − a²)²` en `tests/beta_benchmarks.rs`; la barrera y la posición
//! del TS coinciden con la solución exacta (error relativo de barrera 0).
//!
//! **GSM vs NEB en 5+ reacciones (2026-11-07)**: `tests/alpha_gsm_validation.rs`
//! añade cinco pozos curvos `V(x,y) = (x²−a²)² + ½(y−c x²)²`, cuyo MEP es la
//! parábola `y = c x²` y cuya barrera exacta es `a⁴`. GSM y un NEB discreto
//! (referencia en el test) coinciden con la barrera analítica a ≤0.8 % y entre
//! sí a ≤0.8 %; el criterio 3 queda cubierto para cinco reacciones.
//!
//! ### EXP-205c: Constant Potential Method (CPM) — PROMOVIDO (`src/cpm`)
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
//! **Nota 2026-10-13**: CPM promovido a core (`src/cpm`, F95 vs LAMMPS
//! `electrode/conp`); `beta::cpm` sobrevive como re-export compatible.
//!
//! **Benchmark vs ruta estable (2026-10-23)**: cargas de dos sitios vs la
//! solución cerrada del sistema lineal de CPM en `tests/beta_benchmarks.rs`
//! (1.3e-15) y residual de estacionariedad < 1e-10.
//!
//! ### EXP-206: Alpha-tier subsystems
//!
//! These remain gated and are NOT promoted:
//!
//! | Subsystem | Feature flag | Location | Validation test | Notes |
//! |---|---|---|---|---|
//! | Kohn–Sham DFT | `alpha-dft` | `src/dft/` | `tests/alpha_dft_pyscf.rs` (F77) | E ≤ 1 mHa vs PySCF; `coulomb_energy` is still 0.0 and the ECP path is unfitted |
//! | ReaxFF | `alpha-reaxff` | `src/forcefield/reaxff/` | `tests/alpha_reaxff_validation.rs` (F78) | bond orders vs LAMMPS `reaxff` |
//! | ML force field | `alpha-mlff` | `src/ml/{mlff,symmetry_functions}.rs` | `tests/alpha_mlff_validation.rs` (F79) | AEV vs a numpy BP oracle; no trained weights, no MD17/SPICE benchmark |
//! | Obara–Saika ERI engine | `alpha-obara-saika` | `src/scf/obara_saika.rs` | `tests/alpha_obara_saika.rs` (F80) | Boys/ERI/Schwarz vs scipy + libcint (s-type shells only) |
//! | Interactive MD | `alpha-imd` | `src/imd.rs` | `tests/alpha_imd_validation.rs` (F85) | steering forces vs an analytic oracle |
//! | Live MD | `alpha-dynamics-live` | `src/dynamics_live/` | `tests/alpha_dynamics_validation.rs` (F84) | forces vs RDKit UFF; NVE drift and reversibility |
//! | Reaction dynamics | `alpha-reaction-dynamics` | `src/alpha/reaction_dynamics/` | `tests/alpha_reaction_dynamics_3d.rs` | self-consistency only; no ASE NEB parity |
//! | EDL adapters | `alpha-edl` | `src/edl/` | `tests/alpha_edl_validation.rs` (F86) | analytic electrochemistry |
//! | Periodic adapters | `alpha-periodic-linear` | `src/alpha/periodic_linear/` | `tests/alpha_periodic_validation.rs` (F88) | k-mesh vs ASE |
//! | Kinetics (HTST) | `alpha-kinetics` | `src/alpha/kinetics/` | `tests/alpha_kinetics_validation.rs` (F87) | Eyring + scipy `solve_ivp` |
//! | Render-bridge builders | `alpha-render-bridge` | `src/alpha/render_bridge/` | `tests/alpha_render_validation.rs` (F90) | JSON Schema 2020-12 |
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
