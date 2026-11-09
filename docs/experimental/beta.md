# Beta Modules

**Beta** now contains one remaining gated module (MBH) plus compatibility aliases for methods promoted to core. The promoted paths stay available under their historical `beta-*` flags.

> Only MBH requires an explicit opt-in. KPM, CPM, RandNLA and Riemannian are core.

---

## Module Index

| ID | Module | Feature flag | Description |
|----|--------|-------------|-------------|
| B1 | `kpm` (`beta::kpm` alias) | core | Kernel Polynomial Method — O(N) DOS (promoted 2026-11-07) |
| B2 | `beta::mbh` | `beta-mbh` | Mobile Block Hessian vibrational analysis |
| B3 | `rand_nla` (core) | `beta-randnla` alias | Randomized Nyström EHT + low-rank Lanczos path for N>500 |
| B4 | `riemannian` (core, F94) | — | Riemannian geometry on the PSD cone |
| B5 | `cpm` (core, F95) | — | Constant Potential Method charge equilibration |

---

## B1 — Kernel Polynomial Method (core `kpm`, `beta-kpm` alias)

> Promoted to core on 2026-11-07. Use `sci_form::kpm`; the old
> `sci_form::beta::kpm` path remains as an alias.

**Module:** `sci_form::kpm` (core since 2026-11-07; F91 validated: moments 3.5e-15, DOS 1.7e-14, Mulliken; sparse CSR 5000 sites)

Computes the density of states (DOS) and stochastic Mulliken charges using a **Chebyshev polynomial expansion** of the Hamiltonian. Avoids diagonalization entirely, achieving $O(N)$ cost for sparse systems.

### Theory

$$\rho(E) \approx \frac{2}{\pi a} \sum_{k=0}^{M} g_k^{(M)} \mu_k \frac{T_k\!\left(\frac{E-b}{a}\right)}{\sqrt{1 - \left(\frac{E-b}{a}\right)^2}}$$

where $\mu_k = \frac{1}{n_v}\sum_v r_v^T T_k(\tilde{H}) r_v$ are stochastic moments and $g_k^{(M)}$ is the Jackson damping kernel.

### Rust

```rust
use sci_form::kpm::{compute_kpm_dos, compute_kpm_mulliken, KpmConfig};

let config = KpmConfig {
    order: 200,
    n_random: 30,
    e_min: -25.0,
    e_max: 5.0,
    n_points: 500,
    ..Default::default()
};

// DOS from geometry (builds EHT Hamiltonian internally)
let dos = compute_kpm_dos(&elements, &positions, &config)?;
// dos.energies, dos.dos, dos.homo_lumo_gap

// Stochastic Mulliken charges
let mulliken = compute_kpm_mulliken(&elements, &positions)?;
// mulliken.charges, mulliken.total_electrons
```

### WASM / TypeScript

```typescript
import { beta_compute_kpm_dos, beta_compute_kpm_mulliken } from 'sci-form-wasm/beta';

const dosResult = JSON.parse(beta_compute_kpm_dos(elements_json, coords_json, config_json));
// { energies, dos, homo_lumo_gap, gap }

const mullikenResult = JSON.parse(beta_compute_kpm_mulliken(elements_json, coords_json));
// { charges, total_electrons }
```

### Python

```python
from sci_form.beta import kpm_dos, kpm_mulliken

dos = kpm_dos(elements, coords)
print(f"KPM gap: {dos.gap:.3f} eV")

mulliken = kpm_mulliken(elements, coords)
print(f"charges: {mulliken.charges}")
```

### Parallelization

The stochastic trace estimation (`ChebyshevExpansion::from_matrix`) generates all probe vectors upfront for determinism, then processes each vector independently via `rayon::par_iter` under `features = ["parallel"]`. Speedup is linear with core count up to `n_vectors`.

---

## B2 — Mobile Block Hessian (`beta-mbh`)

**Module:** `sci_form::beta::mbh` (gated; F92 validated 2026-10-14: stretch 1.3e-5 SI; Ghysels gap)

Computes vibrational frequencies for large molecules by partitioning atoms into **mobile blocks** and computing a reduced-dimension Hessian. Dramatically reduces the cost of normal-mode analysis versus the full $3N \times 3N$ Hessian.

### Rust

```rust
use sci_form::beta::mbh::{compute_mbh_frequencies, MbhConfig};

// Blocks: each Vec<usize> is a set of atom indices treated as a rigid unit
let blocks: Vec<Vec<usize>> = vec![
    vec![0, 1, 2],    // block 0 — e.g. CH3
    vec![3, 4, 5, 6], // block 1
];

let config = MbhConfig { max_blocks: 20, ..Default::default() };
let result = compute_mbh_frequencies("CCO", &coords_flat, Some(blocks), &config)?;
// result.frequencies_cm1, result.n_modes, result.zero_point_energy
```

### WASM / TypeScript

```typescript
import { beta_compute_mbh_frequencies } from 'sci-form-wasm/beta';

const result = JSON.parse(beta_compute_mbh_frequencies(elements_json, coords_json, smiles));
// { frequencies_cm1, n_modes, zero_point_energy_hartree, converged }
```

### Python

```python
from sci_form.beta import mbh_frequencies

result = mbh_frequencies(elements, coords, smiles="CCO")
print(f"lowest mode: {result.frequencies_cm1[0]:.1f} cm⁻¹")
print(f"ZPE: {result.zero_point_energy:.4f} Ha")
```

---

## B3 — Randomized Nyström EHT (core `rand_nla`, `beta-randnla` alias)

**Module:** `sci_form::rand_nla` (core since 2026-11-08; F93 validated; low-rank N>500 path added 2026-11-08)

Approximates the EHT eigenvalue problem using a **randomized Nyström sketch**: projects the $N \times N$ Hamiltonian onto a $k$-dimensional random subspace before computing eigenvalues. Cost drops from $O(N^3)$ (full diag) to $O(N k^2)$ — enabling EHT for systems of thousands of atoms.

$$H \approx H_{\text{sketch}} = (HΩ)(Ω^T HΩ)^{-1}(HΩ)^T$$

where $Ω \in \mathbb{R}^{N \times k}$ is a random Gaussian sketch matrix.

### Rust

```rust
use sci_form::rand_nla::{
    solve_eht_randnla, solve_eht_randnla_lowrank, RandNlaConfig,
};

let config = RandNlaConfig {
    sketch_size: Some(50), // k; None = 50 + ceil(log2 N)
    seed: 42,
    max_error: 1e-3,
    fallback_enabled: true,
    power_iterations: 2,
};

// Full spectrum (N × N generalized eigensolver).
let (energies, c, info) = solve_eht_randnla(&h, &s, &config);
// info.k, info.residual_error, info.used_fallback

// N > 500: lowest states only, no N × N diagonalization.
let low = solve_eht_randnla_lowrank(&h, &s, &config, 4);
// low.eigenvalues, low.eigenvectors, low.info
```

The convenience wrapper `eht::solve_eht_randnla(elements, positions, k, &config)`
assembles the EHT matrices and returns `(EhtResult, RandNlaInfo)`.

### WASM / TypeScript

```typescript
import { beta_solve_eht_randnla } from 'sci-form-wasm/beta';

const result = JSON.parse(beta_solve_eht_randnla(elements_json, coords_json, config_json));
// { orbital_energies, homo_index, homo_energy, lumo_energy, gap, k, residual_error, used_fallback }
```

### Python

```python
from sci_form.beta import beta_solve_eht_randnla

result = beta_solve_eht_randnla(elements, coords, sketch_size=50)
print(f"RandNLA gap: {result.gap:.3f} eV  (k={result.k})")
```

---

## B4 — Riemannian PSD Geometry (`beta-riemannian`)

**Module:** `sci_form::riemannian` (core since 2026-10-14, F94: projection 8.2e-12, optimizer like scipy; Frobenius proxy)

Operations on the **manifold of positive semidefinite (PSD) matrices**, used by the SDR embedding (`alpha-sdr`) and Riemannian L-BFGS optimizer. Provides geodesic distance, projections onto the PSD cone, and Riemannian gradient maps.

### Key operations

| Function | Description |
|----------|-------------|
| `psd_distance(A, B)` | Affine-invariant Riemannian distance $d(A,B) = \|\log(A^{-1/2}B A^{-1/2})\|_F$ |
| `psd_projection(M)` | Project a symmetric matrix onto the PSD cone (clamp negative eigenvalues) |
| `riemannian_gradient(G, X)` | Map Euclidean gradient G at point X to the Riemannian gradient |

### Rust

```rust
use sci_form::riemannian::{psd_distance, psd_projection};
use nalgebra::DMatrix;

let a = DMatrix::identity(4, 4);
let b = DMatrix::from_fn(4, 4, |i, j| if i == j { 2.0 } else { 0.1 });

let dist = psd_distance(&a, &b, 4)?;
println!("geodesic distance: {dist:.6}");

let m = DMatrix::from_fn(3, 3, |i, j| if i == j { -0.5 } else { 0.1 });
let proj = psd_projection(&m, 3);
// All eigenvalues of proj are ≥ 0
```

### WASM / TypeScript

```typescript
import { beta_psd_distance, beta_psd_projection } from 'sci-form-wasm/beta';

const dist = JSON.parse(beta_psd_distance(matrix_a_row_major_json, matrix_b_row_major_json, dim)).distance;
const proj = JSON.parse(beta_psd_projection(matrix_row_major_json, dim)).matrix;
```

### Python

```python
from sci_form.beta import beta_psd_distance, beta_psd_projection

dist = beta_psd_distance(a_flat, b_flat, dim=4)
proj = beta_psd_projection(m_flat, dim=3)
```

---

## B5 — Constant Potential Method (`beta-cpm`)

**Module:** `sci_form::cpm` (core since 2026-10-13, F95 vs LAMMPS electrode/conp)

Electrochemical charge equilibration at a **fixed external electrostatic potential** $V_0$. Models electrodes and surfaces by constraining the total charge to match the applied potential, relevant for electrolyte-electrode simulations.

$$Q_{\text{total}} = \arg\min_q \left[ E_{\text{electrostatic}}(q) - V_0 \sum_i q_i \right]$$

### Rust

```rust
use sci_form::cpm::{compute_cpm_charges, CpmConfig};

let config = CpmConfig { potential: -0.5, ..Default::default() }; // V vs vacuum
let result = compute_cpm_charges(&elements, &positions, &config)?;
// result.charges, result.total_charge, result.energy
println!("total charge: {:.4} e  (target potential: {} V)", result.total_charge, config.potential);
```

### WASM / TypeScript

```typescript
import { beta_compute_cpm_charges } from 'sci-form-wasm/beta';

const result = JSON.parse(beta_compute_cpm_charges(elements_json, coords_json, -0.5));
// { charges, total_charge, energy, converged }
```

### Python

```python
from sci_form.beta import cpm_charges

result = cpm_charges(elements, coords, potential=-0.5)
print(f"CPM charges: {result.charges}")
print(f"total charge: {result.total_charge:.4f} e")
```

---

## Enabling Beta Features

### Rust

```toml
[dependencies]
# KPM is core since 2026-11-07; beta-kpm still resolves as an alias.
sci-form = { version = "2.0", features = ["beta-mbh"] }

# All still-gated beta modules:
sci-form = { version = "2.0", features = ["beta-mbh"] }

# Beta + parallelization:
sci-form = { version = "2.0", features = ["beta-mbh", "parallel"] }
```

### WASM / npm

```bash
npm install sci-form-wasm
```

```typescript
// Subpath import — beta functions only
import { beta_compute_kpm_dos, beta_solve_eht_randnla } from 'sci-form-wasm/beta';

// Mixed stable + beta
import init, * as sf from 'sci-form-wasm';
await init();
sf.beta_compute_kpm_dos(elements_json, coords_json, config_json);
```

### Python

```python
# Beta submodule
from sci_form.beta import kpm_dos, eht_randnla, cpm_charges

# Or import the whole beta namespace
import sci_form.beta as beta
result = beta.kpm_dos(elements, coords)
```

### Build with all beta features (WASM)

```bash
pnpm wasm:build:all-beta
# Equivalent to:
bash crates/wasm/build.sh --web-only --web-features "beta-mbh"
```

> Always build through `crates/wasm/build.sh` (directly or through the
> `pnpm wasm:build*` scripts).  Calling `wasm-pack` by hand regenerates
> `pkg/package.json` from scratch and drops the `snippets` entry that
> `wasm-bindgen-rayon` needs at runtime, as well as the `./alpha` and
> `./beta` subpath exports.

---

## Stability Contract

- Beta modules have **stable API shapes** — breaking changes require a minor-version bump.
- Beta modules are covered by **regression tests** and compared against reference implementations.
- Each beta module documents its **known limitations** (e.g., element coverage, system-size constraints).
- Graduating to **stable** requires: 2+ validated real-world cases, benchmarks, no known edge-case failures, and a deprecation cycle for any API changes.

---

## Comparison: Alpha vs Beta vs Stable

| Property | Alpha | Beta | Stable |
|----------|-------|------|--------|
| API stability | May break | Minor-version bumps only | Fully stable |
| Test coverage | Basic integration tests | Regression + reference validation | Comprehensive + benchmarks |
| Documentation | Module-level | Full API docs | Full API + guides |
| Real-world validation | None required | Recommended | Required |
| Element coverage | Limited | Broad | Complete within scope |
| Performance | Not yet optimized | Benchmarked | Optimized |
