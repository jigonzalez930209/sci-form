# Experimental Modules

sci-form includes a set of **experimental algorithms** behind feature flags. These modules explore next-generation computational chemistry methods without affecting the stable production API.

> All experimental modules live under `sci_form::experimental::*` and require explicit feature flags to compile.

---

## Track Summary

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

---

## Groups

### Group A — Advanced Geometry

| Track | Focus |
|-------|-------|
| [E1: CGA](./cga) | Unify rotations/translations via geometric algebra motors |
| [E2: RandNLA](./randnla) | Sub-cubic EHT diagonalization with randomized sketching |
| [E3: Riemannian](./riemannian) | PSD-guaranteed embedding via manifold optimization |

### Group B — Electronic Structure

| Track | Focus |
|-------|-------|
| [E4: KPM](./kpm) | Linear-scaling DOS and density matrix via Chebyshev expansion |
| [E5: EEQ](./eeq) | Geometry-dependent electronegativity equalization charges |
| [E6: ALPB](./alpb) | Analytical implicit solvation with Born radii |

### Group C — Precision Methods

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
sci-form = { version = "0.8", features = ["experimental-d4", "experimental-kpm"] }
```

### Build with all experimental features

```bash
cargo build --features "experimental-cga,experimental-randnla,experimental-riemannian,experimental-kpm,experimental-eeq,experimental-alpb,experimental-d4,experimental-sdr,experimental-mbh,experimental-cpm,experimental-gsm"
```

### Run tests for a specific track

```bash
cargo test --features experimental-d4 --test regression -- test_d4
```

---

## Design Principles

1. **Isolation** — Experimental code never modifies stable modules (`src/conformer/`, `src/eht/`, etc.)
2. **Feature gates** — Each track compiles only when its flag is enabled
3. **Independent testing** — Each track has its own test suite under `tests/experimental/`
4. **Graduation path** — Tracks may be promoted to production after passing all validation tests and code review
