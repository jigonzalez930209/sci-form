# Experimental Modules

sci-form keeps research code separate from the stable surface. That separation is
smaller than it used to be: most of the old alpha and beta code is now core, and
the historical names remain only as compatibility aliases. This page is the
current map.

## The three states

| State | What it means | Feature flag |
|---|---|---|
| **Core** | Compiled by default, covered by a cross-library or analytic gate, and part of the stable API | none |
| **Gated research** | Validated in a narrower regime but not ready for the stable surface | an explicit `alpha-*` or `beta-*` flag |
| **Compatibility alias** | The old path or flag still resolves, but the implementation lives in core | old flag only |

## What moved to core

The following modules cleared their promotion criteria and now live at the
top-level paths. The old paths still compile.

| Old path | Core path | Evidence |
|---|---|---|
| `alpha::cga` | `sci_form::cga` | motors 3.4e-15 vs Rodrigues + `clifford` |
| `alpha::sdr` | `sci_form::sdr` | embedding 7e-15 vs CVXPY/SCS |
| `alpha::edl` | `sci_form::edl` | worst 2e-14 vs analytic Gouy-Chapman-Stern |
| `experimental-eeq` | `sci_form::charges_eeq` | charges 6.9e-12 vs multicharge |
| `experimental-alpb` | `sci_form::solvation_alpb` | electrostatic term 3.14 kcal/mol vs tblite |
| `experimental-d4` | `sci_form::dispersion` | D4(PBE) 6.2e-10 Ha vs dftd4 |
| `alpha::gsm` | `sci_form::gsm` | Müller-Brown barrier 1.6e-3 vs ASE CI-NEB |
| `alpha::kinetics` | `sci_form::kinetics` | Eyring 1.9e-16; traces vs SciPy `solve_ivp` |
| `alpha::periodic_linear` | `sci_form::periodic_linear` | k-mesh and Bloch phases vs ASE |
| `alpha::dynamics_live` | `sci_form::dynamics_live` | UFF forces 3.6e-13 vs RDKit |
| `alpha::render_bridge` | `sci_form::render_bridge` | JSON Schema 2020-12 + pack round-trip |
| `beta::kpm` | `sci_form::kpm` | moments 3.5e-15; DOS 1.7e-14 vs exact diagonalization |
| `beta::cpm` | `sci_form::cpm` | LAMMPS `fix electrode/conp` bilayer comparison |
| `beta::rand_nla` | `sci_form::rand_nla` | rank-N/4 3.8e-4; N=600 low-rank path 9.6e-4 |
| `beta::riemannian` | `sci_form::riemannian` | projection 8.2e-12 vs SciPy L-BFGS-B |
| `forcefield::reaxff` | `sci_form::forcefield::reaxff` | bond orders vs LAMMPS; analytical gradients |
| `scf::obara_saika` | `sci_form::scf::obara_saika` | Boys/ERI/Schwarz vs SciPy + libcint |

## What is still gated

Only four research lines remain behind a flag. Each one has a concrete reason,
not a vague "not finished" label.

| Module | Flag | Why it stays gated |
|---|---|---|
| **Kohn-Sham density functional theory** | `alpha-dft` | Energy error is 7.2e-5 Ha vs PySCF (target 1e-5), `coulomb_energy` is still zero, and the ECP path is not fitted |
| **Machine-learned force field** | `alpha-mlff` | Descriptors and gradients match a NumPy oracle, but there are no trained production weights and no MD17/SPICE benchmark |
| **Mobile Block Hessian** | `beta-mbh` | Stretch test 1.3e-5 SI; the experimental reference is ±2 cm⁻¹ (Ghysels) and global-minimum integration is not covered |
| **Reaction dynamics** | `alpha-reaction-dynamics` | Self-consistency only; there is no ASE NEB parity for the full 3D path |

## Plain-language names

The documentation uses the descriptive name first and the short code only when
the API requires it.

| Short name | Plain name | Where it lives |
|---|---|---|
| CGA | Conformal geometric algebra | core `cga` |
| GSM | Growing String Method | core `gsm` |
| KPM | Kernel Polynomial Method | core `kpm` |
| RandNLA | Randomized numerical linear algebra | core `rand_nla` |
| MBH | Mobile Block Hessian | gated `beta::mbh` |
| CPM | Constant Potential Method | core `cpm` |
| EDL | Electrical double layer | core `edl` |
| ALPB | Analytical linearized Poisson-Boltzmann | core `solvation_alpb` |
| EEQ | Electronegativity equalization | core `charges_eeq` |
| SDR | Semidefinite relaxation | core `sdr` |
| D4 | Dispersion correction with dynamic \(C_6/C_8\) | core `dispersion` |
| HTST | Harmonic transition-state theory | core `kinetics` |

## Enabling features

Core modules need no flag. The aliases below are only for builds that still
name the historical feature.

```toml
[dependencies]
# Gated research
sci-form = { version = "2.0", features = ["alpha-dft", "alpha-mlff"] }
sci-form = { version = "2.0", features = ["beta-mbh"] }
sci-form = { version = "2.0", features = ["alpha-reaction-dynamics"] }

# Everything + parallel
sci-form = { version = "2.0", features = [
    "alpha-dft", "alpha-mlff", "alpha-reaction-dynamics",
    "beta-mbh", "parallel"
] }
```

### Run the tests

The `tests/experimental/` suite is aggregated by the **`experimental_gates`**
target. A command that targets `--test regression -- test_<module>` matches no
test and appears to pass without running anything.

```bash
# Gated research modules
cargo test --features alpha-dft,alpha-mlff --lib --test 'alpha_*'
cargo test --features beta-mbh --lib --test 'beta_*'

# The aggregated experimental suite
cargo test --features "beta-mbh alpha-cga alpha-sdr" --test experimental_gates
```

| Module | Filter | Tests |
|---|---|---|
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

The legacy `experimental-*` names still work (`experimental-cga` →
`alpha-cga`, `experimental-kpm` → `beta-kpm`, and so on), but the modern
`alpha-*` / `beta-*` names are the ones to use.

## Design principles

1. **Isolation.** Research code does not modify the stable conformer, EHT or
   force-field paths until it passes its gate.
2. **One flag per line.** A gated module compiles only when its flag is on.
3. **Independent tests.** Every line has its own oracle or analytic reference.
4. **A promotion is a path move, not a fork.** Promoted modules keep a
   re-export at the old path so downstream code keeps compiling.
5. **Numbers, not adjectives.** "Validated" means a named reference, a measured
   deviation, a confidence region and a recorded limitation. The complete
   report is in [Scientific validation](/validation/scientific-validation).

## Where to read the numbers

- [Scientific validation](/validation/scientific-validation) — per-algorithm
  reference, measured deviation, confidence region and limitations.
- [Full validation matrix](/validation/full-matrix) — all 129 recorded rows.
- [Cross-library plan (historical, Spanish)](/validation/cross-library-validation.es)
  — the working plan that produced these gates; kept for provenance.
