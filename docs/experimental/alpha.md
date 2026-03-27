# Alpha Modules

**Alpha** modules expose new computational chemistry methods that are functionally correct but not yet production-ready. They may have rough APIs, limited element coverage, or missing edge-case handling. Breaking changes between minor versions are possible.

> Alpha modules are gated by `alpha-*` feature flags. No alpha code is compiled without an explicit opt-in.

---

## Module Index

| ID | Module | Feature flag | Description |
|----|--------|-------------|-------------|
| A1 | `alpha::dft` | `alpha-dft` | Kohn-Sham DFT (SVWN / PBE) |
| A2 | `alpha::reaxff` | `alpha-reaxff` | ReaxFF reactive force field |
| A3 | `alpha::mlff` | `alpha-mlff` | Neural network force field with AEV descriptors |
| A4 | `alpha::obara_saika` | `alpha-obara-saika` | Obara-Saika ERIs and Boys function |
| A5 | `alpha::cga` | `alpha-cga` | CGA motor algebra for dihedral refinement |
| A6 | `alpha::gsm` | `alpha-gsm` | Growing String Method for reaction paths |
| A7 | `alpha::sdr` | `alpha-sdr` | Semidefinite relaxation embedding |
| A8 | `alpha::dynamics` | `alpha-dynamics-live` | Live interactive molecular dynamics |
| A9 | `alpha::imd` | `alpha-imd` | IMD wire-protocol bridge (NAMD/GROMACS) |

---

## A1 — Kohn-Sham DFT (`alpha-dft`)

**Module:** `sci_form::alpha::dft`

Density functional theory solver using a minimal Gaussian basis. Supports SVWN (LDA) and PBE (GGA) functionals. For systems up to ~20 atoms in the current implementation.

### Rust

```rust
use sci_form::alpha::dft::{compute_dft, DftConfig, DftMethod};

let conf = sci_form::embed("CCO", 42);
let pos: Vec<[f64;3]> = conf.coords.chunks(3).map(|c| [c[0],c[1],c[2]]).collect();

let config = DftConfig { method: DftMethod::Pbe, max_iter: 200, ..Default::default() };
let result = compute_dft(&conf.elements, &pos, &config)?;
println!("DFT total energy: {:.6} Ha", result.total_energy);
println!("HOMO: {:.3} eV,  LUMO: {:.3} eV,  gap: {:.3} eV",
    result.homo_energy, result.lumo_energy, result.gap);
```

### WASM / TypeScript

```typescript
import { alpha_compute_dft } from 'sci-form-wasm/alpha';

const result = JSON.parse(alpha_compute_dft(elements_json, coords_json, "pbe"));
console.log(`gap: ${result.gap.toFixed(3)} eV`);
```

### Python

```python
from sci_form.alpha import dft_calculate
result = dft_calculate(elements, coords, method="pbe")
print(f"gap: {result.gap:.3f} eV")
```

---

## A2 — ReaxFF Reactive Force Field (`alpha-reaxff`)

**Module:** `sci_form::alpha::reaxff`

Bond-order-based reactive force field for bond breaking and formation. Includes EEM charge equilibration and central-difference gradients (parallelized with the `parallel` feature).

### Rust

```rust
use sci_form::alpha::reaxff::{compute_reaxff_energy, compute_reaxff_gradient};

let (energy, grad) = compute_reaxff_gradient(&coords_flat, &elements, &ReaxffParams::default_chon())?;
println!("ReaxFF energy: {:.4} kcal/mol", energy);
// grad: [∂E/∂x0, ∂E/∂y0, ∂E/∂z0, ...]
```

### WASM / TypeScript

```typescript
import { alpha_compute_reaxff_gradient, alpha_compute_reaxff_energy } from 'sci-form-wasm/alpha';

const gradResult = JSON.parse(alpha_compute_reaxff_gradient(elements_json, coords_json));
// { energy, gradient, n_atoms }

const energyResult = JSON.parse(alpha_compute_reaxff_energy(elements_json, coords_json));
// { total_energy, bonded, nonbonded, electrostatic, bond_orders }
```

### Python

```python
from sci_form.alpha import reaxff_gradient, reaxff_energy
grad = reaxff_gradient(elements, coords)         # (energy, gradient)
breakdown = reaxff_energy(elements, coords)       # detailed energy components
```

### Parallelization

The central-difference gradient loop is parallelized under `features = ["parallel"]`. Each coordinate displacement pair runs on its own rayon thread — scales linearly with CPU cores for large systems.

---

## A3 — Neural Network Force Field (`alpha-mlff`)

**Module:** `sci_form::alpha::mlff`

Machine-learned force field using Atomic Environment Vectors (AEVs) / Behler-Parrinello symmetry functions as descriptors. The AEV computation (`compute_aevs`) is parallelized with the `parallel` feature.

### Rust

```rust
use sci_form::alpha::mlff::{compute_mlff, MlffConfig};
use sci_form::alpha::symmetry_functions::{compute_aevs, SymmetryFunctionParams};

let pos: Vec<[f64;3]> = conf.coords.chunks(3).map(|c| [c[0],c[1],c[2]]).collect();

// Compute AEV descriptors (parallelized with `parallel` feature)
let params = SymmetryFunctionParams::default();
let aevs = compute_aevs(&conf.elements, &pos, &params);

// MLFF energy + forces
let config = MlffConfig::default();
let result = compute_mlff(&conf.elements, &pos, &config)?;
println!("MLFF energy: {:.4} eV", result.energy);
```

### WASM / TypeScript

```typescript
import { alpha_compute_mlff, alpha_compute_aevs } from 'sci-form-wasm/alpha';

const aevResult = JSON.parse(alpha_compute_aevs(elements_json, coords_json, params_json));
// { aevs: [[...], ...], n_atoms, n_radial, n_angular }

const mlffResult = JSON.parse(alpha_compute_mlff(elements_json, coords_json, config_json));
// { energy, forces, atomic_energies }
```

### Python

```python
from sci_form.alpha import alpha_compute_aevs, alpha_compute_mlff

aev_result = alpha_compute_aevs(elements, coords)
mlff_result = alpha_compute_mlff(elements, coords)
print(f"MLFF energy: {mlff_result.energy:.4f} eV")
```

### Parallelization

`compute_aevs` outer atom loop: each atom's AEV is fully independent → `rayon::par_iter` activated by `features = ["parallel"]`. Scales O(1) per atom with thread count.

---

## A4 — Obara-Saika ERIs and Boys Function (`alpha-obara-saika`)

**Module:** `sci_form::alpha::obara_saika`

Two-electron repulsion integrals (ERIs) via the Obara-Saika recurrence, and the Boys function $F_n(x)$ needed for Gaussian-based SCF. Required for `alpha-dft` and forms the integral engine for HF/DFT calculations.

### Rust

```rust
use sci_form::alpha::obara_saika::{boys_function, schwarz_bound};
use sci_form::alpha::obara_saika::eri::{Shell, eri_ssss};

// Boys function F_n(x)
let f0 = boys_function(0, 2.5);  // F_0(2.5)

// Schwarz upper bound for ERI screening
let bound = schwarz_bound(1.2, &[0.0, 0.0, 0.0], 0.8, &[1.5, 0.0, 0.0]);

// (ss|ss) ERI
let shell_s = Shell { alpha: 1.0, center: [0.0, 0.0, 0.0] };
let eri = eri_ssss(&shell_s, &shell_s, &shell_s, &shell_s);
```

### WASM / TypeScript

```typescript
import { alpha_boys_function, alpha_eri_ssss, alpha_schwarz_bound } from 'sci-form-wasm/alpha';

const F0 = JSON.parse(alpha_boys_function(0, 2.5)).value;
const eri = JSON.parse(alpha_eri_ssss(shellA_json, shellB_json, shellC_json, shellD_json)).eri;
```

---

## A5 — CGA Dihedral Refinement (`alpha-cga`)

**Module:** `sci_form::alpha::cga`

Conformal Geometric Algebra (CGA) motor rotation for dihedral angle refinement. Provides numerically exact rigid-body rotation of atom subtrees around arbitrary bond axes — avoids the small-angle approximation errors of Cartesian torsion updates.

### Rust

```rust
use sci_form::alpha::cga::{rotate_dihedral, refine_torsion};

// Rotate atom subtree around bond axis by angle_rad
let new_coords = rotate_dihedral(&coords_flat, &subtree_indices, &axis_a, &axis_b, angle_rad);

// Iterative torsion refinement toward target dihedral
let refined = refine_torsion(&coords_flat, &subtree_indices, &axis_a, &axis_b, target_rad, 50)?;
```

### WASM / TypeScript

```typescript
import { alpha_rotate_dihedral_cga, alpha_refine_torsion_cga } from 'sci-form-wasm/alpha';

const rotated = JSON.parse(alpha_rotate_dihedral_cga(
    coords_json, subtree_json, axis_a_json, axis_b_json, angle_rad
));
// { coords: [...], n_rotated }
```

---

## A6 — Growing String Method (`alpha-gsm`)

**Module:** `sci_form::alpha::gsm`

Reaction path interpolation and transition-state search via the Growing String Method. Nodes are added from both reactant and product ends toward the transition state, using UFF/MMFF for node relaxation.

### Rust

```rust
use sci_form::alpha::gsm::{gsm_interpolate, gsm_find_ts, GsmConfig};

// Interpolate geometry at parameter t in [0,1]
let interp = gsm_interpolate(&reactant_coords, &product_coords, 0.5);

// Find transition state
let config = GsmConfig { max_nodes: 11, ..Default::default() };
let ts_result = gsm_find_ts("CCO>>CC=O", &reactant_coords, &product_coords, &config)?;
println!("TS energy: {:.4} kcal/mol", ts_result.ts_energy);
```

### WASM / TypeScript

```typescript
import { alpha_gsm_interpolate, alpha_gsm_find_ts } from 'sci-form-wasm/alpha';

const interp = JSON.parse(alpha_gsm_interpolate(reactant_json, product_json, 0.5));
// { coords, n_atoms }

const ts = JSON.parse(alpha_gsm_find_ts(smiles, reactant_json, product_json, config_json));
// { ts_coords, ts_energy, path_energies, converged }
```

---

## A7 — Semidefinite Relaxation Embedding (`alpha-sdr`)

**Module:** `sci_form::alpha::sdr`

Embeds atoms in 3D by solving a semidefinite program (SDP) relaxation of the distance geometry problem. Conceptually eliminates the retry loop of ETKDG for challenging ring systems by finding the global minimum of the distance bounds satisfaction.

### Rust

```rust
use sci_form::alpha::sdr::{sdr_embed, SdrConfig, DistanceBound};

let bounds = vec![
    DistanceBound { i: 0, j: 1, lower: 1.4, upper: 1.6 },
    DistanceBound { i: 1, j: 2, lower: 1.4, upper: 1.6 },
    // ...
];

let config = SdrConfig { max_iter: 500, tolerance: 1e-6, ..Default::default() };
let result = sdr_embed(&bounds, n_atoms, &config)?;
// result.coords — flat [x0,y0,z0,...]
```

### WASM / TypeScript

```typescript
import { alpha_sdr_embed } from 'sci-form-wasm/alpha';

const result = JSON.parse(alpha_sdr_embed(distance_pairs_json, n_atoms, config_json));
// { coords, converged, iterations, final_violation }
```

---

## A8 & A9 — Live MD and IMD Bridge (`alpha-dynamics-live`, `alpha-imd`)

**Module:** `sci_form::alpha::dynamics`

**A8 (`alpha-dynamics-live`):** Real-time molecular dynamics integration (Verlet / velocity Verlet) with per-frame callbacks — suitable for interactive visualization or frame streaming to a WebSocket.

**A9 (`alpha-imd`):** IMD (Interactive Molecular Dynamics) wire-protocol implementation for connecting to running NAMD or GROMACS simulations. Reads force/coordinate streams and sends external forces back.

Both are implicitly included in `alpha-imd` (A9 depends on A8).

### WASM / TypeScript

```typescript
import { LiveSimulation } from 'sci-form-wasm/alpha';

const sim = new LiveSimulation(elements_json, coords_json, config_json);
sim.on_frame(frameJson => {
    const frame = JSON.parse(frameJson);
    // frame.coords, frame.step, frame.kinetic_energy, frame.potential_energy
    render(frame.coords);
});
sim.start();
```

---

## Enabling Alpha Features

### Rust

```toml
[dependencies]
sci-form = { version = "0.11", features = ["alpha-dft", "alpha-reaxff", "alpha-mlff"] }

# All alpha modules:
sci-form = { version = "0.11", features = [
    "alpha-dft", "alpha-reaxff", "alpha-mlff",
    "alpha-obara-saika", "alpha-cga", "alpha-gsm",
    "alpha-sdr", "alpha-dynamics-live", "alpha-imd"
] }

# Alpha + parallelization:
sci-form = { version = "0.11", features = [
    "alpha-mlff", "alpha-reaxff", "parallel"
] }
```

### WASM / npm

```bash
npm install sci-form-wasm
```

```typescript
// Subpath import — alpha functions only
import { alpha_compute_dft, alpha_compute_reaxff_gradient } from 'sci-form-wasm/alpha';

// Or default import for stable + alpha mixed usage
import init, * as sf from 'sci-form-wasm';
await init();
sf.alpha_compute_dft(elements_json, coords_json, "pbe");
```

### Python

```python
# Alpha submodule
from sci_form.alpha import dft_calculate, reaxff_gradient, alpha_compute_mlff

# Or import the whole alpha namespace
import sci_form.alpha as alpha
result = alpha.dft_calculate(elements, coords)
```

### Build with all alpha features (WASM)

```bash
pnpm wasm:build:all-alpha
# Equivalent to:
wasm-pack build crates/wasm --features "alpha-dft,alpha-reaxff,alpha-mlff,alpha-obara-saika,alpha-cga,alpha-gsm,alpha-sdr,alpha-dynamics-live,alpha-imd"
```

---

## Stability Contract

- Alpha modules may **change API shape** between minor versions.
- Alpha modules are **not covered** by the stable-API backward-compatibility promise.
- Each alpha module has integration tests; check `tests/alpha_*` for current coverage.
- Graduating to **beta** requires: passing all tests, documented API, basic benchmarks.
- Graduating from beta to **stable** requires: 2+ real-world validation cases, no known edge-case failures.
