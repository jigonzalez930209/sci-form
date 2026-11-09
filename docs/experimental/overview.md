# Experimental Modules

sci-form organises pre-production algorithms into three tiers based on maturity:

| Tier | Stability | API guarantee | Coverage |
|------|-----------|--------------|----------|
| **Alpha** | Proof-of-concept | May break between minor versions | Basic integration tests |
| **Beta** | Validated | Minor-version bumps only | Regression + reference validation |
| **Experimental** (legacy) | Track-by-track | See individual track page | Varies |

---

## [Alpha Modules →](./alpha)

The alpha tier now mixes compatibility aliases for promoted methods with the remaining gated research modules (DFT, MLFF, reaction dynamics).

| ID | Module | Feature flag | Description |
|----|--------|-------------|-------------|
| A1 | `alpha::dft` | `alpha-dft` | Kohn-Sham DFT (SVWN / PBE) |
| A2 | `forcefield::reaxff` (core, F78) | — | ReaxFF reactive force field |
| A3 | `alpha::mlff` | `alpha-mlff` | Neural network FF with AEV descriptors |
| A4 | `scf::obara_saika` (core, F80) | — | Obara-Saika ERIs + Boys function |
| A5 | `cga` (core, F81) | — | CGA motor algebra dihedral refinement |
| A6 | `gsm` (core) | `alpha-gsm` alias | Growing String Method |
| A7 | `sdr` (core, F83) | — | Semidefinite relaxation embedding |
| A8 | `dynamics_live` (core) | `alpha-dynamics-live` alias | Live interactive MD |
| A9 | `imd` (core, F85) | — | IMD steering forces |

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

Most of the beta tier is now core; only MBH remains gated.

| ID | Module | Feature flag | Description |
|----|--------|-------------|-------------|
| B1 | `kpm` (`beta::kpm` alias) | core | Kernel Polynomial Method — O(N) DOS (promoted 2026-11-07) |
| B2 | `beta::mbh` | `beta-mbh` | Mobile Block Hessian |
| B3 | `rand_nla` (`beta::rand_nla` alias) | core | Randomized Nyström EHT + low-rank Lanczos for N>500 (promoted 2026-11-08) |
| B4 | `riemannian` (core, F94) | — | Riemannian PSD cone geometry |
| B5 | `cpm` (core, F95) | — | Constant Potential Method |

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

> Gated modules live behind `alpha-*`/`beta-*` features (legacy
> `experimental-*` flags alias them); promoted modules need no feature.
> There is no `sci_form::experimental` Rust path.

### Track Summary

| Track | Module | Feature Flag | Description |
|-------|--------|-------------|-------------|
| E1 | `cga` | `experimental-cga` | Conformal Geometric Algebra G(4,1) |
| E2 | `rand_nla` | core (`experimental-randnla` alias) | Randomized NLA for O(N²k) EHT, low-rank N>500 path |
| E3 | `riemannian` | core (`experimental-riemannian` alias) | Riemannian L-BFGS for ETKDG |
| E4 | `kpm` | core (`experimental-kpm` alias) | Kernel Polynomial Method — O(N) DOS |
| E5 | `eeq` | `experimental-eeq` | Dynamic EEQ Charge Model |
| E6 | `alpb` | `experimental-alpb` | Analytical Linearized Poisson-Boltzmann |
| E7 | `d4` | `experimental-d4` | DFT-D4 Dispersion Correction |
| E8 | `sdr` | `experimental-sdr` | Semidefinite Relaxation Embedding |
| E9 | `mbh` | `experimental-mbh` | Mobile Block Hessian |
| E10 | `cpm` | core (`experimental-cpm` alias) | Constant Potential Method |
| E11 | `gsm` | core (`experimental-gsm` alias) | Growing String Method |

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
sci-form = { version = "2.0", features = ["alpha-dft", "alpha-mlff"] }

# Beta modules (KPM is core since 2026-11-07; beta-kpm is a compatibility alias)
sci-form = { version = "2.0", features = ["beta-mbh"] }

# Legacy experimental
sci-form = { version = "2.0", features = ["alpha-reaction-dynamics"] }

# Everything + parallel
sci-form = { version = "2.0", features = [
    "alpha-dft", "alpha-mlff", "alpha-reaction-dynamics",
    "beta-mbh",
    "parallel"
] }
```

### Run tests

The `tests/experimental/` suite is aggregated by the **`experimental_gates`** test target
(`tests/experimental_gates.rs`), not by `regression`. A command of the form
`cargo test --test regression -- test_<module>` matches **no test** and therefore
"succeeds" without running anything — always use `experimental_gates`.

```bash
# Alpha tests
cargo test --features alpha-dft,alpha-reaxff,alpha-mlff --lib --test 'alpha_*'

# Beta tests
cargo test --features beta-mbh --lib --test 'beta_*'

# The experimental/ suite (needs beta-cpm so that the aggregated CPM
# benchmark module compiles, plus the feature of the module under test)
cargo test --features "beta-mbh alpha-cga alpha-sdr" --test experimental_gates
```

| Module | Command suffix | Tests |
|--------|----------------|-------|
| ALPB | `-- test_alpb` | 8 |
| CGA | `-- test_cga` | 20 |
| CPM | `-- test_cpm` | 9 |
| D4 | `-- test_d4` | 8 |
| EEQ | `-- test_eeq` | 8 |
| GSM | `-- test_gsm` | 8 |
| KPM | `-- test_kpm` | 9 |
| MBH | `-- test_mbh` | 8 |
| RandNLA | `-- test_randnla` | 8 |
| Riemannian | `-- test_riemannian` | 11 |
| SDR | `-- test_sdr` | 8 |

The legacy `experimental-*` feature names still work as aliases
(`experimental-cga` → `alpha-cga`, `experimental-kpm` → `beta-kpm`, and so on),
but the modern `alpha-*` / `beta-*` names are the ones to use.

---

## Design Principles

1. **Isolation** — Experimental code never modifies stable modules (`src/conformer/`, `src/eht/`, etc.)
2. **Feature gates** — Each track compiles only when its flag is enabled
3. **Independent testing** — Each track has its own test suite under `tests/experimental/`
4. **Graduation path** — Tracks may be promoted through alpha → beta → stable after validation
