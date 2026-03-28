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
| A10 | `alpha::edl` | `alpha-edl` | Electrochemical double-layer contracts and Helmholtz CPU reference |
| A11 | `alpha::periodic_linear` | `alpha-periodic-linear` | Periodic k-mesh, Bloch-phase, and linear-scaling spectral contracts |
| A12 | `alpha::kinetics` | `alpha-kinetics` | HTST and microkinetic contracts with CPU reference rates |
| A13 | `alpha::render_bridge` | `alpha-render-bridge` | Shared chart and trajectory payload contracts for alpha workflows |

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

## A10 — Electrochemistry / EDL (`alpha-edl`)

**Module:** `sci_form::alpha::edl`

Defines serializable electrical double-layer inputs and outputs for electrochemical interface work. Models include Helmholtz (compact-layer only), Gouy-Chapman (diffuse-layer only), and Gouy-Chapman-Stern (coupled compact/diffuse). Contains CPU reference solvers and adapters for CPM, EEQ, and ALPB coupling.

### Rust

```rust
use sci_form::alpha::edl::*;

// Compute individual model profiles
let helmholtz = compute_helmholtz_profile(0.2, &EdlConfig::default())?;
let gc = compute_gouy_chapman_profile(0.2, &EdlConfig::default())?;
let gcs = compute_gcs_profile(0.2, &EdlConfig::default())?;

// Unified dispatcher (auto-selects model from config)
let config = EdlConfig {
    model: EdlModel::GouyChapmanStern,
    ionic_strength_m: 0.1,
    ..Default::default()
};
let profile = compute_edl_profile(0.2, &config)?;
println!("Capacitance: {:.3} F/m²", profile.differential_capacitance.total_f_per_m2);

// Scan capacitance over potential range
let scan = scan_edl_capacitance(-0.5, 0.5, 50, &config)?;
for (v, c) in &scan {
    println!("  V={:6.3}, C={:10.6}", v, c);
}
```

### Python

```python
from sci_form.alpha import (
    edl_profile, edl_helmholtz, edl_gouy_chapman, edl_gcs,
    edl_capacitance_scan, edl_chart, capacitance_chart
)

# Standard dispatcher
profile = edl_profile(surface_potential_v=0.2, model="gcs", ionic_strength_m=0.1)
print(f"Total drop: {profile.total_interfacial_drop_v:.4f} V")

# Model-specific functions
helmholtz = edl_helmholtz(0.2, ionic_strength_m=0.1)
gc = edl_gouy_chapman(0.2, ionic_strength_m=0.1)
gcs = edl_gcs(0.2, ionic_strength_m=0.1, stern_thickness_a=3.0)

# Capacitance scan
scan_data = edl_capacitance_scan(v_min=-0.5, v_max=0.5, n_points=50)

# Visualization payloads
edl_plot = edl_chart(surface_potential_v=0.2)
cap_plot = capacitance_chart(v_min=-0.5, v_max=0.5, n_points=50)
print(f"EDL chart: {edl_plot.title}, {len(edl_plot.series)} series")
```

### WASM / TypeScript

```typescript
import {
    compute_edl_profile,
    compute_helmholtz_profile,
    compute_gouy_chapman_profile,
    compute_gcs_profile,
    compute_edl_capacitance_scan,
} from 'sci-form-wasm/alpha';

// Profile computation
const profile = JSON.parse(compute_edl_profile(
    0.2,           // surface_potential_v
    "gcs",         // model
    0.1,           // ionic_strength_m
    298.15,        // temperature_k
    128,           // n_points
    12.0           // extent_angstrom
));
console.log(`Capacitance: ${profile.capacitance_total_f_per_m2.toFixed(6)} F/m²`);

// Scan
const scan = JSON.parse(compute_edl_capacitance_scan(-0.5, 0.5, 50, 0.1, "gcs"));
console.log(`Scan points: ${scan.length}`);
```

### CLI

```bash
# Individual model profiles
sci-form edl-helmholtz --potential 0.2 --ionic-strength 0.1 --temperature 298.15
sci-form edl-gouy-chapman --potential 0.2 --ionic-strength 0.1
sci-form edl-gcs --potential 0.2 --ionic-strength 0.1 --stern 3.0

# Standard dispatcher
sci-form edl-profile --potential 0.2 --model gcs --ionic-strength 0.1

# Scan mode
sci-form edl-scan --v-min -0.5 --v-max 0.5 --n-points 50 --model gcs
```

---

## A11 — Periodic Linear-Scaling (`alpha-periodic-linear`)

**Module:** `sci_form::alpha::periodic_linear`

Defines reusable periodic spectral contracts and CPU reference reciprocal-space helpers. Supports Monkhorst-Pack and gamma-centered k-meshes, Bloch-phase factors, k-point integration, band structure, and density-of-states calculations. Integrates with future KPM and randomized low-rank adapters.

### Rust

```rust
use sci_form::alpha::periodic_linear::*;

// Generate k-mesh
let mesh_config = KMeshConfig {
    grid: [6, 6, 6],
    centering: KMeshCentering::MonkhorstPack,
};
let kmesh = monkhorst_pack_mesh(&mesh_config)?;
println!("n_k = {}, weights sum = {:.6}", kmesh.points.len(), 
    kmesh.points.iter().map(|p| p.weight).sum::<f64>());

// Build Bloch Hamiltonian at a k-point
let phase = bloch_phase(kmesh.points[0].fractional, [1, 0, 0]);
println!("Bloch phase: {:?}", phase);

// Band structure from assembled operators
let bs = compute_band_structure_adapter(&hamiltonians, &overlaps, kmesh, 10)?;
println!("Band edges: HOMO {} eV, LUMO {} eV, gap {} eV",
    bs.homo_edge_ev, bs.lumo_edge_ev, bs.band_edges.indirect_gap_ev);

// BZ integration
let dos_integrated = bz_integrate_scalar(&dos_values, &kmesh_weights);
println!("Integrated DOS: {}", dos_integrated);
```

### Python

```python
from sci_form.alpha import kmesh

# Generate k-mesh
mesh = kmesh(grid=[4, 4, 4], centering="monkhorst-pack")  # or "gamma"
print(f"n_k = {mesh.n_points}, grid = {mesh.grid}")
print(f"Weights sum: {sum({mesh.weights)}")
```

### WASM / TypeScript

```typescript
import { compute_kmesh } from 'sci-form-wasm/alpha';

// Generate k-mesh
const mesh = JSON.parse(compute_kmesh('[4,4,4]', 'monkhorst-pack'));
console.log(`n_k = ${mesh.n_points}, weights sum = ${mesh.weights.reduce((a,b)=>a+b,0)}`);
```

### CLI

```bash
# Generate k-mesh
sci-form kmesh --grid 4 4 4 --centering monkhorst-pack
# Output: {n_points, grid, weights_sum}
```

---

---

## A12 — Kinetics (`alpha-kinetics`)

**Module:** `sci_form::alpha::kinetics`

Defines elementary-step, HTST, and microkinetic trace contracts for catalytic and surface-chemistry workflows. Computes transition-state-theory rates, temperature sweeps, and microkinetic steady-state solutions. Integrates with GSM and MBH for barrier estimation.

### Rust

```rust
use sci_form::alpha::kinetics::*;

// Single temperature HTST rate
let step = ElementaryStep {
    step_id: "A→B".into(),
    activation_free_energy_ev: 0.45,
    reaction_free_energy_ev: -0.10,
    prefactor_s_inv: None,
};
let state = ThermodynamicState {
    temperature_k: 298.15,
    pressure_bar: 1.0,
};
let rate = evaluate_htst_rate(&step, state)?;
println!("k_f = {:.3e} s⁻¹, k_r = {:.3e} s⁻¹", 
    rate.forward_rate_s_inv, rate.reverse_rate_s_inv);

// Temperature sweep
let temperatures = vec![200.0, 300.0, 400.0, 500.0];
let sweep = evaluate_htst_temperature_sweep(&step, &temperatures, 1.0)?;
for result in &sweep {
    println!("T={:6.1}K: k_f={:10.3e}, K_eq={:8.4}", 
        result.state.temperature_k, result.forward_rate_s_inv, result.equilibrium_constant);
}

// Microkinetic network solver
let network = vec![step];
let trace = solve_microkinetic_network(&network, &state)?;
println!("Converged: {}, steps: {}", trace.converged, trace.iteration_count);
```

### Python

```python
from sci_form.alpha import htst_rate, htst_temperature_sweep

# Single point
rate = htst_rate(
    activation_free_energy_ev=0.45,
    reaction_free_energy_ev=-0.10,
    temperature_k=298.15
)
print(f"k_f = {rate.forward_rate_s_inv:.3e} s⁻¹")

# Temperature sweep
results = htst_temperature_sweep(
    activation_free_energy_ev=0.45,
    reaction_free_energy_ev=-0.10,
    temperatures_k=[200.0, 300.0, 400.0, 500.0]
)
for r in results:
    print(f"T={r['temperature_k']:.1f}K: k_f={r['forward_rate_s_inv']:.3e}")
```

### WASM / TypeScript

```typescript
import { compute_htst_rate, compute_htst_sweep } from 'sci-form-wasm/alpha';

// Single point
const rate = JSON.parse(compute_htst_rate(0.45, -0.10, 298.15));
console.log(`k_f = ${rate.forward_rate_s_inv.toExponential(3)} s⁻¹`);

// Sweep
const temps = [200, 300, 400, 500];
const sweep = JSON.parse(compute_htst_sweep(0.45, -0.10, JSON.stringify(temps)));
sweep.forEach(r => console.log(`T=${r.temperature_k}K: k_f=${r.forward_rate_s_inv.toExponential(3)}`));
```

### CLI

```bash
# Single rate
sci-form htst-rate --activation 0.45 --reaction -0.10 --temperature 298.15

# Output: {step_id, forward_rate_s_inv, reverse_rate_s_inv, equilibrium_constant, temperature_k}
```

---

## A13 — Render Bridge (`alpha-render-bridge`)

**Module:** `sci_form::alpha::render_bridge`

Defines solver-agnostic payloads for charts, traces, and trajectory visualizations. Shared chart generator for EDL profiles, band structures, microkinetic population traces, Arrhenius plots, and capacitance scans. Converts computed results into Arrow RecordBatch format for transport to web frontends or Jupyter notebooks.

### Rust

```rust
use sci_form::alpha::render_bridge::*;
use sci_form::alpha::edl::*;
use sci_form::alpha::kinetics::*;

// EDL profile chart
let profile = compute_gouy_chapman_profile(0.2, &EdlConfig::default())?;
let chart = edl_profile_chart(&profile);
println!("Chart: {}", chart.title);  // "EDL Profile (Gouy-Chapman)"
for series in &chart.series {
    println!("  {}: {} points", series.label, series.x.len());
}

// Arrhenius plot from temperature sweep
let step = ElementaryStep { /* ... */ };
let temps = vec![200.0, 300.0, 400.0, 500.0];
let sweep = evaluate_htst_temperature_sweep(&step, &temps, 1.0)?;
let arrhenius = arrhenius_chart(&sweep);
println!("Arrhenius plot ready: {} series", arrhenius.series.len());

// Pack to transportable format
let batch = pack_chart_payload(&chart);
println!("Packed: {} columns", batch.float_columns.len());

// Verify round-trip
let valid = verify_chart_round_trip(&chart, 1e-10);
println!("Round-trip valid: {}", valid);
```

### Python

```python
from sci_form.alpha import edl_chart, capacitance_chart

# EDL visualization
plot = edl_chart(surface_potential_v=0.2, ionic_strength_m=0.1)
print(f"Title: {plot.title}")
print(f"Series: {plot.series}")
for s in plot.series:
    print(f"  - {s.label}: {len(s.x)} points ({s.x_unit} vs {s.y_unit})")

# Capacitance scan visualization
cap_plot = capacitance_chart(v_min=-0.5, v_max=0.5, n_points=50)
print(f"Cap plot: {cap_plot.title}, {len(cap_plot.series)} traces")
```

### WASM / TypeScript

```typescript
import { pack_chart_payload_wasm } from 'sci-form-wasm/alpha';

const chartPayload = {
    title: "My EDL Profile",
    series: [
        {
            series_id: "potential",
            label: "Electrostatic Potential",
            x: [0, 0.1, 0.2, ...],
            y: [1.2, 1.1, 0.95, ...],
            x_unit: "Å",
            y_unit: "V"
        }
    ]
};

const packed = JSON.parse(pack_chart_payload_wasm(JSON.stringify(chartPayload)));
console.log(`Packed chart: ${packed.n_columns} columns, ${packed.column_names.join(', ')}`);
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
    "alpha-sdr", "alpha-dynamics-live", "alpha-imd",
    "alpha-edl", "alpha-periodic-linear", "alpha-kinetics",
    "alpha-render-bridge"
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
wasm-pack build crates/wasm --features "alpha-dft,alpha-reaxff,alpha-mlff,alpha-obara-saika,alpha-cga,alpha-gsm,alpha-sdr,alpha-dynamics-live,alpha-imd,alpha-edl,alpha-periodic-linear,alpha-kinetics,alpha-render-bridge"
```

---

## Stability Contract

- Alpha modules may **change API shape** between minor versions.
- Alpha modules are **not covered** by the stable-API backward-compatibility promise.
- Each alpha module has integration tests; check `tests/alpha_*` for current coverage.
- Graduating to **beta** requires: passing all tests, documented API, basic benchmarks.
- Graduating from beta to **stable** requires: 2+ real-world validation cases, no known edge-case failures.
