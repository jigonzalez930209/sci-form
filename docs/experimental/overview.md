# Experimental Modules

sci-form organises pre-production algorithms into three tiers based on maturity:

| Tier | Stability | API guarantee | Coverage |
|------|-----------|--------------|----------|
| **Alpha** | Proof-of-concept | May break between minor versions | Basic integration tests |
| **Beta** | Validated | Minor-version bumps only | Regression + reference validation |
| **Experimental** (legacy) | Track-by-track | See individual track page | Varies |

---

## [Alpha Modules →](./alpha)

New methods gated by `alpha-*` feature flags. Functionally correct but API shapes may still change.

| ID | Module | Feature flag | Description |
|----|--------|-------------|-------------|
| A1 | `alpha::dft` | `alpha-dft` | Kohn-Sham DFT (SVWN / PBE) |
| A2 | `alpha::reaxff` | `alpha-reaxff` | ReaxFF reactive force field |
| A3 | `alpha::mlff` | `alpha-mlff` | Neural network FF with AEV descriptors |
| A4 | `alpha::obara_saika` | `alpha-obara-saika` | Obara-Saika ERIs + Boys function |
| A5 | `alpha::cga` | `alpha-cga` | CGA motor algebra dihedral refinement |
| A6 | `alpha::gsm` | `alpha-gsm` | Growing String Method |
| A7 | `alpha::sdr` | `alpha-sdr` | Semidefinite relaxation embedding |
| A8 | `alpha::dynamics` | `alpha-dynamics-live` | Live interactive MD |
| A9 | `alpha::imd` | `alpha-imd` | IMD wire-protocol bridge |

```typescript
// WASM: subpath import
import { alpha_compute_dft, alpha_compute_reaxff_gradient } from 'sci-form-wasm/alpha';
```

```python
# Python: submodule
from sci_form.alpha import dft_calculate, reaxff_gradient
```

---

## [Beta Modules →](./beta)

Validated methods gated by `beta-*` feature flags. API is stable within minor versions.

| ID | Module | Feature flag | Description |
|----|--------|-------------|-------------|
| B1 | `beta::kpm` | `beta-kpm` | Kernel Polynomial Method — O(N) DOS |
| B2 | `beta::mbh` | `beta-mbh` | Mobile Block Hessian |
| B3 | `beta::randnla` | `beta-randnla` | Randomized Nyström EHT — O(N k²) |
| B4 | `beta::riemannian` | `beta-riemannian` | Riemannian PSD cone geometry |
| B5 | `beta::cpm` | `beta-cpm` | Constant Potential Method |

```typescript
// WASM: subpath import
import { beta_compute_kpm_dos, beta_solve_eht_randnla } from 'sci-form-wasm/beta';
```

```python
# Python: submodule
from sci_form.beta import kpm_dos, eht_randnla
```

---

## Legacy Experimental Tracks

The older `experimental-*` feature flags remain available at their existing module paths:

> All experimental modules live under `sci_form::experimental::*` and require explicit feature flags to compile.

### Track Summary

| Track | Module | Feature Flag | Description |
|-------|--------|-------------|-------------|
| E1 | `cga` | `experimental-cga` | Conformal Geometric Algebra G(4,1) |
| E2 | `rand_nla` | `experimental-randnla` | Randomized NLA for O(N log N) EHT |
| E3 | `riemannian` | `experimental-riemannian` | Riemannian L-BFGS for ETKDG |
| E4 | `kpm` | `experimental-kpm` | Kernel Polynomial Method — O(N) DOS |
| E5 | `eeq` | `experimental-eeq` | Dynamic EEQ Charge Model |
| E6 | `alpb` | `experimental-alpb` | Analytical Linearized Poisson-Boltzmann |
| E7 | `d4` | `experimental-d4` | DFT-D4 Dispersion Correction |
| E8 | `sdr` | `experimental-sdr` | Semidefinite Relaxation Embedding |
| E9 | `mbh` | `experimental-mbh` | Mobile Block Hessian |
| E10 | `cpm` | `experimental-cpm` | Constant Potential Method |
| E11 | `gsm` | `experimental-gsm` | Growing String Method |

### Groups

#### Group A — Advanced Geometry

| Track | Focus |
|-------|-------|
| [E1: CGA](./cga) | Unify rotations/translations via geometric algebra motors |
| [E2: RandNLA](./randnla) | Sub-cubic EHT diagonalization with randomized sketching |
| [E3: Riemannian](./riemannian) | PSD-guaranteed embedding via manifold optimization |

#### Group B — Electronic Structure

| Track | Focus |
|-------|-------|
| [E4: KPM](./kpm) | Linear-scaling DOS and density matrix via Chebyshev expansion |
| [E5: EEQ](./eeq) | Geometry-dependent electronegativity equalization charges |
| [E6: ALPB](./alpb) | Analytical implicit solvation with Born radii |

#### Group C — Precision Methods

| Track | Focus |
|-------|-------|
| [E7: D4](./d4) | Geometry-dependent dispersion with BJ damping |
| [E8: SDR](./sdr) | Convex SDP embedding eliminating the retry loop |
| [E9: MBH](./mbh) | Reduced-DOF vibrational analysis for large molecules |
| [E10: CPM](./cpm) | Electrochemical charge equilibration at fixed potential |
| [E11: GSM](./gsm) | Reaction path and transition state search |

---

## Enabling Features

### Rust

```toml
[dependencies]
# Alpha modules
sci-form = { version = "0.11", features = ["alpha-dft", "alpha-reaxff"] }

# Beta modules
sci-form = { version = "0.11", features = ["beta-kpm", "beta-randnla"] }

# Legacy experimental
sci-form = { version = "0.11", features = ["experimental-d4", "experimental-kpm"] }

# Everything + parallel
sci-form = { version = "0.11", features = [
    "alpha-dft", "alpha-reaxff", "alpha-mlff",
    "beta-kpm", "beta-randnla",
    "experimental-d4",
    "parallel"
] }
```

### Run tests

```bash
# Alpha tests
cargo test --features alpha-dft,alpha-reaxff,alpha-mlff --lib --test 'alpha_*'

# Beta tests
cargo test --features beta-kpm,beta-mbh,beta-randnla --lib --test 'beta_*'

# Legacy experimental tests
cargo test --features experimental-d4 --test regression -- test_d4
```

---

## Design Principles

1. **Isolation** — Experimental code never modifies stable modules (`src/conformer/`, `src/eht/`, etc.)
2. **Feature gates** — Each track compiles only when its flag is enabled
3. **Independent testing** — Each track has its own test suite under `tests/experimental/`
4. **Graduation path** — Tracks may be promoted through alpha → beta → stable after validation
