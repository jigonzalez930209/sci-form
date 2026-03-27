# Roadmap: Dynamic Simulation & Chemical Reactivity Engine

**Goal:** Evolve sci-form from a static property calculator into a real-time dynamic simulation engine for molecular dynamics, chemical reactivity, and interactive molecular manipulation — all runnable in any browser via WebAssembly.

**Scope:** 7 implementation tasks, from WASM memory architecture through ab-initio quantum chemistry, DFT, classical MD, reactive force fields, ML potentials, and the interactive browser loop.

**Relation to existing code:** sci-form already ships conformer generation (ETKDGv2), static force fields (UFF, MMFF94), quantum methods (EHT, PM3, PM3(tm), GFN0/1/2-xTB, HF-3c, CISD), ANI neural potentials, a Velocity Verlet + Nosé-Hoover MD module (`src/dynamics.rs`), and WASM bindings (`crates/wasm/`). This roadmap builds on that foundation.

---

## Existing Foundation (already implemented)

| Capability | Module | Status |
|-----------|--------|--------|
| UFF / MMFF94 force fields | `src/forcefield/` | Production |
| EHT + orbital grids | `src/eht/` | Production |
| PM3 / PM3(tm) SCF | `src/pm3/` | Production |
| GFN0/1/2-xTB | `src/xtb/` | Production |
| HF-3c + CISD | `src/hf/`, `src/scf/` | Production |
| ANI-2x / ANI-TM | `src/ani/` | Production |
| Velocity Verlet MD (UFF/PM3/xTB) | `src/dynamics.rs` | Production |
| Nosé-Hoover NVT thermostat | `src/dynamics.rs` | Production |
| NEB path search | `src/dynamics.rs` | Production |
| WASM MD trajectory export | `crates/wasm/src/dynamics.rs` | Production |
| Gasteiger charges | `src/charges/` | Production |
| EEQ charges | `src/charges_eeq/` | Production |
| Materials / space groups | `src/materials/` | Production |
| SCF integral engine (STO-3G, 3-21G, 6-31G) | `src/scf/` | Production |
| DIIS convergence accelerator | `src/scf/diis.rs` | Production |

---

## Task 1: Zero-Copy WASM Memory Architecture for Dynamic State

**Priority:** High — prerequisite for real-time interactive simulation at 60 FPS.

**Current state:** The WASM bridge serializes molecular state as JSON strings on every call. This works for single-shot calculations but creates millisecond-scale overhead per frame, incompatible with 60 FPS rendering.

**Target:** Expose the simulation state (positions, velocities, forces) as a contiguous flat buffer in Rust's linear memory, accessible from JavaScript via `Float64Array` views — zero serialization cost per frame.

### 1.1 Dynamic atom structure

Extend the existing `Atom` struct to carry velocity and force vectors for time-stepping:

```rust
// src/dynamics/state.rs (new)
pub struct DynamicAtom {
    pub element: u8,
    pub mass: f64,
    pub charge: f64,
    pub position: [f64; 3],
    pub velocity: [f64; 3],
    pub force: [f64; 3],
}

pub struct MolecularSystem {
    pub atoms: Vec<DynamicAtom>,
    pub potential_energy: f64,
    pub kinetic_energy: f64,
    pub time_fs: f64,
    /// Flat buffer: [x0, y0, z0, x1, y1, z1, ...] — exposed to JS via pointer
    pub positions_flat: Vec<f64>,
}
```

### 1.2 Direct pointer FFI via wasm-bindgen

```rust
// crates/wasm/src/dynamics_live.rs (new)
#[wasm_bindgen]
pub struct LiveSimulation {
    system: MolecularSystem,
}

#[wasm_bindgen]
impl LiveSimulation {
    /// Returns a raw pointer to the flat positions buffer in WASM linear memory.
    /// JS reads this via: new Float64Array(wasm.memory.buffer, ptr, len)
    pub fn positions_ptr(&self) -> *const f64 {
        self.system.positions_flat.as_ptr()
    }

    pub fn positions_len(&self) -> usize {
        self.system.positions_flat.len()
    }

    /// Run N integration steps, then sync the flat buffer for rendering.
    pub fn step(&mut self, dt_fs: f64, substeps: usize) {
        for _ in 0..substeps {
            self.system.integrate(dt_fs);
        }
        self.system.sync_flat_buffer();
    }
}
```

### 1.3 Implementation steps

| Step | Action | Files affected | Depends on |
|------|--------|---------------|------------|
| 1.1a | Create `DynamicAtom` and `MolecularSystem` structs with flat buffer sync | `src/dynamics/state.rs` (new) | — |
| 1.1b | Add `From<ConformerResult>` conversion to initialize dynamics from embedding | `src/dynamics/state.rs` | 1.1a |
| 1.2a | Create `LiveSimulation` WASM wrapper with pointer-based export | `crates/wasm/src/dynamics_live.rs` (new) | 1.1a |
| 1.2b | Add TypeScript helper to wrap `Float64Array` view from WASM memory | `crates/wasm/ts/live_simulation.ts` (new) | 1.2a |
| 1.3 | Benchmark: measure frame latency (JSON vs pointer path) on 1000-atom system | `tests/benchmarks/` | 1.2a |
| 1.4 | Integration test: verify pointer stability across multiple `step()` calls | `tests/wasm/` | 1.2a |

**Performance target:** < 100 μs per-frame overhead for a 500-atom system (current JSON path: ~2–5 ms).

---

## Task 2: Ab-Initio Hartree-Fock Engine (Obara-Saika Integrals + DIIS SCF)

**Priority:** High — provides the exact-exchange foundation for DFT hybrids and excited-state chemistry.

**Current state:** sci-form already has:
- A complete SCF infrastructure in `src/scf/` (overlap, kinetic, nuclear matrices, Fock construction, DIIS, density matrix, Mulliken analysis)
- Gaussian integral evaluation in `src/scf/gaussian_integrals.rs`
- HF-3c with D3/gCP/SRB corrections in `src/hf/`
- Basis sets: STO-3G, 3-21G, 6-31G in `src/scf/extended_basis.rs`
- CISD excited states in `src/hf/cisd.rs`

**Remaining work:** The current integral engine is functional for small molecules but needs optimization for interactive use on larger systems (50–100 atoms).

### 2.1 Obara-Saika recurrence optimization

Replace or supplement the current ERI evaluation with a fully vectorized Obara-Saika recurrence engine:

| Step | Action | Files affected | Depends on |
|------|--------|---------------|------------|
| 2.1a | Implement OS vertical recurrence relation for primitive ERI with SIMD-friendly layout | `src/scf/obara_saika.rs` (new) | — |
| 2.1b | Implement OS horizontal recurrence for contracted shells | `src/scf/obara_saika.rs` | 2.1a |
| 2.1c | Add Schwarz pre-screening to skip negligible ERI shells ($Q_{ij}Q_{kl} < \tau$) | `src/scf/screening.rs` (new) | 2.1a |
| 2.1d | Benchmark: ERI evaluation speed vs current `gaussian_integrals.rs` on benzene/caffeine | `tests/benchmarks/` | 2.1b, 2.1c |

### 2.2 SCF convergence hardening

The existing DIIS in `src/scf/diis.rs` works. The target here is to harden it for difficult cases:

| Step | Action | Files affected | Depends on |
|------|--------|---------------|------------|
| 2.2a | Add level-shifting fallback when DIIS stagnates (> 30 iterations without < 10× error reduction) | `src/scf/scf_loop.rs`, `src/scf/diis.rs` | — |
| 2.2b | Add fractional occupation number (FON) option for near-degenerate systems | `src/scf/density_matrix.rs` | — |
| 2.2c | Validate convergence on transition-metal test set (Fe, Cu, Zn complexes) | `tests/regression/` | 2.2a |

**Performance target:** < 5 s for benzene/STO-3G SCF in WASM; < 200 ms for ethanol.

---

## Task 3: Kohn-Sham DFT (LDA and GGA Functionals)

**Priority:** Medium-high — DFT with correlation captures chemistry that HF misses (bond energies, hydrogen bonding, reaction barriers).

**Current state:** No DFT implementation exists yet. The SCF loop in `src/scf/` handles the one-electron + Coulomb/Exchange terms. Adding DFT requires:
1. A numerical integration grid
2. XC functional kernels (LDA, GGA)
3. Modification of the Fock matrix construction to include $V_{XC}$

### 3.1 Molecular integration grid

| Step | Action | Files affected | Depends on |
|------|--------|---------------|------------|
| 3.1a | Implement Euler-Maclaurin radial quadrature with Treutler-Ahlrichs radial mapping | `src/dft/grid.rs` (new) | — |
| 3.1b | Implement Lebedev angular quadrature (6–590 point grids, precomputed tables) | `src/dft/lebedev.rs` (new) | — |
| 3.1c | Implement Becke partitioning scheme for multi-center weights | `src/dft/becke.rs` (new) | 3.1a, 3.1b |
| 3.1d | Grid pruning: reduce angular points near nuclei where density is smooth | `src/dft/grid.rs` | 3.1c |

### 3.2 Exchange-correlation functional kernels

| Step | Action | Files affected | Depends on |
|------|--------|---------------|------------|
| 3.2a | Implement SVWN (Slater exchange + VWN-5 correlation) for LDA | `src/dft/functionals/svwn.rs` (new) | — |
| 3.2b | Implement PBE exchange-correlation for GGA (requires density gradient on grid) | `src/dft/functionals/pbe.rs` (new) | 3.1c |
| 3.2c | Build the $V_{XC}$ matrix assembly: evaluate basis functions + gradients on grid, contract with XC kernel | `src/dft/vxc_matrix.rs` (new) | 3.2a or 3.2b, 3.1c |

### 3.3 Kohn-Sham SCF loop

| Step | Action | Files affected | Depends on |
|------|--------|---------------|------------|
| 3.3a | Create KS-DFT Fock matrix builder: $F^{KS} = H_{core} + J + V_{XC}$ (no exact exchange for pure DFT) | `src/dft/ks_fock.rs` (new) | 3.2c |
| 3.3b | Wire KS-DFT into the existing SCF loop (reuse DIIS, density matrix, convergence logic) | `src/scf/scf_loop.rs` | 3.3a |
| 3.3c | Validate: atomization energies of H₂, H₂O, CH₄ within 5 kcal/mol of reference | `tests/regression/` | 3.3b |
| 3.3d | WASM export: `compute_dft(elements, coords, functional)` | `crates/wasm/src/electronic.rs` | 3.3b |

**Performance target:** < 30 s for benzene/STO-3G/PBE in WASM; < 2 s for water.

---

## Task 4: Classical Molecular Dynamics (Integrators & Ensembles)

**Priority:** Medium — extends the existing `src/dynamics.rs` with PBC and production-grade thermostats.

**Current state:** sci-form already implements:
- Velocity Verlet integrator (`simulate_velocity_verlet`)
- Nosé-Hoover NVT thermostat (`simulate_nose_hoover`)
- Multi-backend support (UFF, PM3, xTB)
- NEB path search (`compute_simplified_neb_path`)
- WASM trajectory export (`compute_md_trajectory`, `compute_md_trajectory_nvt`)

**Remaining work:**

### 4.1 Periodic boundary conditions

| Step | Action | Files affected | Depends on |
|------|--------|---------------|------------|
| 4.1a | Add `SimulationBox` struct with minimum-image convention | `src/dynamics/pbc.rs` (new) | — |
| 4.1b | Implement coordinate wrapping and minimum-image distance in the Verlet loop | `src/dynamics.rs` | 4.1a |
| 4.1c | Add neighbor list (Verlet list) to avoid O(N²) force evaluation with PBC | `src/dynamics/neighbor_list.rs` (new) | 4.1a |
| 4.1d | Wire PBC into UFF force evaluation for non-bonded terms | `src/forcefield/` | 4.1b, 4.1c |

### 4.2 Ensemble extensions

| Step | Action | Files affected | Depends on |
|------|--------|---------------|------------|
| 4.2a | Add Berendsen barostat for NPT ensemble (pressure coupling) | `src/dynamics/barostat.rs` (new) | 4.1a |
| 4.2b | Add Langevin thermostat as alternative to Nosé-Hoover (better for implicit solvent) | `src/dynamics.rs` | — |
| 4.2c | Add energy/temperature/pressure trajectory logging per step | `src/dynamics.rs` | — |
| 4.2d | WASM export: `compute_md_trajectory_npt` | `crates/wasm/src/dynamics.rs` | 4.2a |

### 4.3 Validation

| Step | Action | Files affected | Depends on |
|------|--------|---------------|------------|
| 4.3a | Energy conservation test: NVE drift < 0.1% over 10,000 steps for argon/LJ system | `tests/regression/` | 4.1b |
| 4.3b | Temperature stability test: NVT at 300K ± 5K over 50,000 steps | `tests/regression/` | 4.2a |

---

## Task 5: Reactive Force Fields (ReaxFF)

**Priority:** Medium — enables bond-breaking/formation dynamics without quantum cost.

**Current state:** No reactive force field exists. UFF/MMFF94 require fixed topology.

### 5.1 Continuous bond orders

| Step | Action | Files affected | Depends on |
|------|--------|---------------|------------|
| 5.1a | Implement continuous bond-order calculation from interatomic distances (σ + π + ππ exponentials) | `src/forcefield/reaxff/bond_order.rs` (new) | — |
| 5.1b | Implement over-coordination correction (valence penalties) | `src/forcefield/reaxff/bond_order.rs` | 5.1a |
| 5.1c | Implement Hermite tapering polynomial for smooth cutoff ($C^3$ continuity) | `src/forcefield/reaxff/taper.rs` (new) | — |
| 5.1d | Implement bonded energy terms (bond stretch, angle bend, torsion) as functions of continuous BO | `src/forcefield/reaxff/energy.rs` (new) | 5.1a, 5.1c |

### 5.2 Charge equilibration

| Step | Action | Files affected | Depends on |
|------|--------|---------------|------------|
| 5.2a | Implement EEM (Electronegativity Equalization Method) with Lagrange-constrained charge neutrality | `src/forcefield/reaxff/eem.rs` (new) | — |
| 5.2b | Wire EEM into the ReaxFF force loop (recompute charges each step) | `src/forcefield/reaxff/` | 5.2a, 5.1d |
| 5.2c | Implement non-bonded terms (van der Waals + Coulomb with shielding) | `src/forcefield/reaxff/nonbonded.rs` (new) | 5.2a |

### 5.3 ReaxFF parametrization & validation

| Step | Action | Files affected | Depends on |
|------|--------|---------------|------------|
| 5.3a | Parse ReaxFF parameter files (van Duin format) for C/H/O/N | `src/forcefield/reaxff/params.rs` (new) | — |
| 5.3b | Validate: H₂O dissociation barrier, C-C bond breaking energy vs DFT reference | `tests/regression/` | 5.1d, 5.2c, 5.3a |
| 5.3c | Analytical gradient implementation for all ReaxFF terms | `src/forcefield/reaxff/gradients.rs` (new) | 5.1d, 5.2c |
| 5.3d | Wire into Velocity Verlet loop: reactive MD simulation | `src/dynamics.rs` | 5.3c, Task 4 |
| 5.3e | WASM export: `compute_reaxff_energy`, `simulate_reactive_md` | `crates/wasm/src/` | 5.3d |

---

## Task 6: ML Force Field Inference (ONNX / Burn in WASM)

**Priority:** Medium — bridges DFT-level accuracy with classical MD speed.

**Current state:** sci-form already has ANI-2x and ANI-TM potentials in `src/ani/` with AEV descriptors, cutoff functions, and gradient computation. The existing implementation uses hand-coded neural network inference in Rust (no external ML runtime).

**Remaining work:** Generalize to ONNX model loading for user-trained potentials.

### 6.1 Behler-Parrinello symmetry functions (generalized)

| Step | Action | Files affected | Depends on |
|------|--------|---------------|------------|
| 6.1a | Extract and generalize the AEV computation from `src/ani/aev.rs` into a standalone descriptor module | `src/ml/symmetry_functions.rs` (new) | — |
| 6.1b | Add configurable radial/angular function parameter sets (beyond ANI defaults) | `src/ml/symmetry_functions.rs` | 6.1a |

### 6.2 ONNX inference via tract or burn

| Step | Action | Files affected | Depends on |
|------|--------|---------------|------------|
| 6.2a | Add `tract-onnx` as optional dependency; implement model loading and single-atom inference | `src/ml/onnx_inference.rs` (new), `Cargo.toml` | — |
| 6.2b | Implement batch inference: compute descriptors → run model → extract per-atom energies | `src/ml/onnx_inference.rs` | 6.1a, 6.2a |
| 6.2c | Implement analytical forces via chain rule (backprop through descriptors) | `src/ml/onnx_forces.rs` (new) | 6.2b |
| 6.2d | Verify WASM compatibility: compile tract with `wasm32-unknown-unknown` target, measure binary size | `crates/wasm/` | 6.2a |

### 6.3 Integration with dynamics loop

| Step | Action | Files affected | Depends on |
|------|--------|---------------|------------|
| 6.3a | Add `MdBackend::Mlff` variant to the dynamics engine | `src/dynamics.rs` | 6.2c |
| 6.3b | WASM export: `simulate_mlff_md` | `crates/wasm/src/dynamics.rs` | 6.3a |
| 6.3c | Benchmark: inference latency per atom per step (target: < 0.1 ms/atom in WASM) | `tests/benchmarks/` | 6.3b |

**Binary size target:** < 5 MB additional WASM payload for the ML inference runtime.

---

## Task 7: Interactive Molecular Dynamics (IMD) Browser Loop

**Priority:** High (after core MD is stable) — the user-facing goal of the project.

**Current state:** The existing WASM bindings export `compute_md_trajectory` and `compute_md_trajectory_nvt` which run a full trajectory and return results. These are batch APIs, not frame-by-frame.

**Target:** A persistent simulation object that advances one rendering frame at a time, synchronized with `requestAnimationFrame`, with IMD steering forces.

### 7.1 Frame-by-frame simulation API

| Step | Action | Files affected | Depends on |
|------|--------|---------------|------------|
| 7.1a | Create `LiveSimulation` WASM struct (from Task 1) with per-frame `step()` method | `crates/wasm/src/dynamics_live.rs` (new) | Task 1 |
| 7.1b | Add force backend selector: UFF, PM3, xTB, ReaxFF, MLFF | `crates/wasm/src/dynamics_live.rs` | Tasks 4–6 |
| 7.1c | Add configurable substeps-per-frame (decouple physics timestep from render rate) | `crates/wasm/src/dynamics_live.rs` | 7.1a |

### 7.2 Interactive steering (IMD forces)

| Step | Action | Files affected | Depends on |
|------|--------|---------------|------------|
| 7.2a | Add `apply_steering_force(atom_index, target_xyz, spring_k)` API to inject harmonic IMD potentials | `crates/wasm/src/dynamics_live.rs` | 7.1a |
| 7.2b | Add `clear_steering_force(atom_index)` to release IMD anchors | `crates/wasm/src/dynamics_live.rs` | 7.2a |
| 7.2c | Implement force accumulation: IMD forces are added linearly to the evaluated force field | `src/dynamics/` | 7.2a |

### 7.3 requestAnimationFrame integration

| Step | Action | Files affected | Depends on |
|------|--------|---------------|------------|
| 7.3a | Create TypeScript/JS helper: `startSimulationLoop(sim, renderer, options)` using `requestAnimationFrame` | `crates/wasm/ts/simulation_loop.ts` (new) | 7.1a |
| 7.3b | Implement pointer-based position read (no JSON): `new Float64Array(memory.buffer, sim.positions_ptr(), sim.positions_len())` | `crates/wasm/ts/simulation_loop.ts` | Task 1 |
| 7.3c | Add pause/resume/reset controls | `crates/wasm/ts/simulation_loop.ts` | 7.3a |

### 7.4 Integration with 3d-orbitals-simulations viewer

| Step | Action | Files affected | Depends on |
|------|--------|---------------|------------|
| 7.4a | Add `useSimulation` React hook wrapping `LiveSimulation` + rAF loop | `3d-orbitals-simulations/src/hooks/useSimulation.ts` (new) | 7.3a |
| 7.4b | Add simulation controls panel: play/pause, temperature, timestep, backend selector | `3d-orbitals-simulations/src/components/` | 7.4a |
| 7.4c | Wire atom picking (raycaster) → IMD steering force injection | `3d-orbitals-simulations/src/components/viewer/` | 7.2a, 7.4a |

**FPS target:** Sustained 60 FPS for systems up to 500 atoms with UFF backend; 30 FPS with xTB backend.

---

## Dependency Map

```
Task 1 (WASM Memory) ──────────────────────────────────────────┐
                                                                │
Task 2 (HF/Obara-Saika) ─────┐                                 │
                              ├─► Task 3 (DFT)                  │
                              │                                 │
Task 4 (Classical MD + PBC) ──┼─► Task 5 (ReaxFF) ──┐          │
                              │                      ├──► Task 7 (IMD Browser Loop)
                              └─► Task 6 (MLFF) ─────┘          │
                                                                │
                                                   ◄────────────┘
```

**Critical path:** Task 1 → Task 4 → Task 7 (for the fastest path to interactive visualization)

**Parallel tracks:**
- Task 2 + Task 3 can proceed independently (quantum accuracy)
- Task 5 + Task 6 can proceed after Task 4 (reactive/ML backends)

---

## Implementation Phases

### Phase A: Interactive MVP (Tasks 1, 4 partial, 7)
**Goal:** Real-time Velocity Verlet + UFF in the browser at 60 FPS with atom steering.
- Pointer-based WASM export
- Existing Verlet + Nosé-Hoover
- IMD force injection
- rAF loop integration

### Phase B: Quantum Backends (Tasks 2, 3)
**Goal:** DFT single-point calculations and SCF optimization accessible from the browser.
- Obara-Saika integral speedup
- LDA/GGA functionals
- KS-DFT SCF loop

### Phase C: Reactive Chemistry (Tasks 4 complete, 5, 6)
**Goal:** Bond breaking and formation in real time.
- PBC + neighbor lists
- ReaxFF continuous bond orders + EEM charges
- ONNX ML potential inference in WASM

### Phase D: Production (all tasks)
**Goal:** Full interactive molecular dynamics studio in the browser.
- All force backends wired to the IMD loop
- Simulation controls + analysis tools in the viewer
- Educational scenario presets (water boiling, bond breaking, protein folding)

---

## File Structure (new modules)

```
src/
├── dynamics/              # Refactored from dynamics.rs
│   ├── mod.rs
│   ├── state.rs           # DynamicAtom, MolecularSystem (Task 1)
│   ├── verlet.rs          # Existing Velocity Verlet (split from dynamics.rs)
│   ├── nose_hoover.rs     # Existing NVT (split from dynamics.rs)
│   ├── pbc.rs             # Periodic boundary conditions (Task 4)
│   ├── neighbor_list.rs   # Verlet neighbor list (Task 4)
│   └── barostat.rs        # NPT ensemble (Task 4)
├── dft/                   # New (Task 3)
│   ├── mod.rs
│   ├── grid.rs            # Molecular integration grid
│   ├── lebedev.rs         # Angular quadrature tables
│   ├── becke.rs           # Atomic partitioning weights
│   ├── vxc_matrix.rs      # XC potential matrix assembly
│   ├── ks_fock.rs         # Kohn-Sham Fock matrix builder
│   └── functionals/
│       ├── svwn.rs         # LDA
│       └── pbe.rs          # GGA
├── forcefield/
│   └── reaxff/            # New (Task 5)
│       ├── mod.rs
│       ├── bond_order.rs   # Continuous bond orders
│       ├── taper.rs        # Hermite tapering
│       ├── energy.rs       # Bonded energy terms
│       ├── nonbonded.rs    # vdW + Coulomb with shielding
│       ├── eem.rs          # Charge equilibration
│       ├── params.rs       # Parameter file parser
│       └── gradients.rs    # Analytical gradients
├── scf/
│   ├── obara_saika.rs     # New (Task 2) — optimized ERI recurrence
│   └── screening.rs       # New (Task 2) — Schwarz pre-screening
└── ml/
    ├── symmetry_functions.rs  # New (Task 6) — generalized Behler-Parrinello
    ├── onnx_inference.rs      # New (Task 6) — tract runtime
    └── onnx_forces.rs         # New (Task 6) — backprop forces

crates/wasm/
├── src/
│   └── dynamics_live.rs   # New (Task 7) — frame-by-frame simulation API
└── ts/
    ├── live_simulation.ts # New (Task 1) — Float64Array view helper
    └── simulation_loop.ts # New (Task 7) — rAF integration
```
