# Testing Guide

This document describes how to run every test in the project, what each test validates,
and how to reproduce the full pre-release test gate locally.

---

## Prerequisites

| Tool | Minimum version | Install |
|------|-----------------|---------|
| Rust toolchain | stable (1.77+) | `rustup update stable` |
| `maturin` | 1.6+ | `pip install maturin` |
| `wasm-pack` | 0.13+ | `cargo install wasm-pack` |
| Node.js | 18+ | [nodejs.org](https://nodejs.org) |
| Python | 3.9+ | [python.org](https://python.org) |

```bash
# Add cross-compilation targets used by the release workflow
rustup target add \
  aarch64-unknown-linux-gnu \
  x86_64-apple-darwin \
  aarch64-apple-darwin \
  x86_64-pc-windows-msvc
```

---

## 1. Quality Checks

Run these before anything else. They are the first gate in CI.

```bash
# Static analysis — must produce zero warnings
cargo clippy --all-targets -- -D warnings

# Formatting check
cargo fmt --check
```

---

## 2. Unit Tests (inside `src/`)

Unit tests live in `#[cfg(test)]` modules inside each source file:

| File | Tests |
|------|-------|
| `src/smiles.rs` | SMILES round-trip, stereo, isotopes |
| `src/smarts/parser.rs` | SMARTS atom/bond query parsing |
| `src/smarts/torsion_matcher.rs` | Torsion SMARTS matching |
| `src/etkdg.rs` | ETKDGv2 entry point, seed reproducibility |

```bash
cargo test --lib
```

---

## 3. Integration Tests (`tests/`)

Tests are grouped by purpose:

- `tests/ci.rs` is the compact CI gate that exercises every major subsystem.
- `tests/regression/` holds the functional Rust integration suites.
- `tests/analysis/` holds comparison and validation helpers.
- `tests/debug/` holds investigation and trace harnesses.
- `tests/benchmarks/` holds large-data benchmarks.
- `tests/integration/` holds Python and Node.js interop scripts.

All Rust suites pass on the release profile.

The large GDB-20 reference fixtures are stored compressed in git:

- `tests/fixtures/gdb20_reference.json.gz`
- `tests/fixtures/gdb20_ensemble.json.gz`

The Rust tests and internal bins resolve `*.json` first and then fall back to `*.json.gz` automatically.
If you need to regenerate them locally, use:

```bash
python scripts/generate_gdb20_reference.py
python scripts/generate_ensemble_reference.py
```

```bash
# Compact CI gate
cargo test --release --test ci -- --nocapture

# All Rust integration tests
cargo test --tests --release

# Individual suites (useful during development)
cargo test --release --test regression test_diverse_molecules -- --nocapture
cargo test --release --test regression test_geometry_quality -- --nocapture
cargo test --release --test regression test_gradient_check -- --nocapture
cargo test --release --test regression test_ensemble_rmsd -- --nocapture
cargo test --release --test regression test_gdb20_rmsd -- --nocapture
cargo test --release --test regression test_gdb20 -- --nocapture
cargo test --release --test regression test_tet_centers -- --nocapture
```

### Test descriptions

| Test file | What it validates |
|-----------|-------------------|
| `regression/test_diverse_molecules.rs` | Embed 50+ chemically diverse molecules; all produce valid 3-D geometries |
| `regression/test_geometry_quality.rs` | Bond lengths / angles stay within ±3 σ of CSD distributions |
| `regression/test_gradient_check.rs` | Force-field gradients match numerical finite-differences (tolerance 1 × 10⁻⁴) |
| `regression/test_ensemble_rmsd.rs` | Pose diversity: multi-conformer ensembles span expected RMSD range |
| `regression/test_gdb20_rmsd.rs` | RMSD vs. RDKit/ETKDG on GDB-20 subset (benchmark accuracy) |
| `regression/test_gdb20.rs` | Success rate ≥ 99 % on GDB-20 subset |
| `regression/test_tet_centers.rs` | Tetrahedral stereocentres are preserved after embedding |
| `regression/test_session_algorithms.rs` | 33 tests covering all recent session algorithms (see below) |

### Session Algorithm Regression Tests (`test_session_algorithms.rs`)

This suite validates every algorithm implemented in the latest development iteration.

```bash
cargo test --test regression test_session_algorithms
```

| Section | Tests | What it validates |
|---------|-------|-------------------|
| UHF/ROHF (§1) | 6 | H₂ singlet/doublet/triplet, ROHF H₂⁺, configured UHF, α/β symmetry |
| CIF Import/Export (§2) | 4 | Parse NaCl-type, roundtrip fidelity, uncertainty notation `5.640(1)`, required fields |
| AO→MO Transform (§3) | 1 | `MoIntegrals` type existence and public linkage |
| GPU sTDA/Hessian (§4) | 2 | Compile-time checks (`cfg(experimental-gpu)`) |
| PM3 Gaussians (§5) | 5 | Main-group 2-Gaussian params, TM empty gaussians, water/CH₄/H₂ convergence |
| xTB Broyden (§6) | 4 | GFN0 water/H₂ convergence, GFN1 water/ethanol convergence |
| NMR 5J Coupling (§7) | 3 | Naphthalene peri-H 5J, ethane bond counts, coupling type format |
| SMIRKS Multi-Component (§8) | 4 | Multi-component parse/apply, single-component backward compat, dot-separated |
| Population Public API (§9) | 2 | Water/ethanol via `compute_population()` wrapper |
| EEQ Damping (§10) | 4 | Charge neutrality, O negative, charged molecule, H₂O symmetry |

### Population Internal Unit Tests (`src/population/population.rs`)

Internal `#[cfg(test)]` tests for `pub(crate)` functions not accessible from integration tests:

| Test | What it validates |
|------|-------------------|
| `test_valence_electrons_period1_period2` | H, He, C, N, O, F, Ne coverage |
| `test_valence_electrons_period6_main_group` | Cs, Ba, Tl, Pb, Bi, Po, At, Rn (Z=55–86) |
| `test_valence_electrons_lanthanides` | La→Lu (Z=57–71) valence electron counts |
| `test_valence_electrons_3d/4d/5d_transition_metals` | All transition metal series coverage |
| `test_valence_electrons_unknown_returns_zero` | Out-of-range Z returns 0.0 |
| `test_population_parallel_matches_serial` | Parallel ≡ serial Mulliken/Löwdin charges (feature `parallel`) |
| `test_bond_orders_parallel_matches_serial` | Parallel ≡ serial Wiberg/Mayer bond orders (feature `parallel`) |

---

## 4. SMIRKS Reaction Testing

SMIRKS (SMILES Reaction Specification) allows pattern-based chemical reaction transforms.

```bash
# Run all SMIRKS unit tests (23 tests)
cargo test --lib smirks

# Run SMIRKS integration tests (15 tests)
cargo test --test test_smirks_reactions

# Compare with RDKit (requires RDKit installation)
python scripts/compare_smirks_reactions.py

# Layered reaction comparison harness
# - RDKit for SMIRKS semantics
# - sci-form CLI for configurable-method NEB, backend planning, and energetics
# - geomeTRIC ElasticBand using configurable sci-form CLI engine
# - (optional --external) ASE NEB with tblite GFN2-xTB for truly external path validation
# - (optional --external) tblite GFN1/GFN2 and PySCF HF/STO-3G single-point cross-validation
python scripts/compare_reaction_layers.py

# Multi-method NEB (each backend runs NEB path independently)
python scripts/compare_reaction_layers.py --methods uff,pm3,xtb

# With external reference validation (requires: pip install tblite ase pyscf)
# NOTE: tblite may need LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libgfortran.so.5
python scripts/compare_reaction_layers.py --methods uff,pm3 --external

# JSON output
python scripts/compare_reaction_layers.py --methods uff,pm3 --external --json

# Python integration tests
python tests/integration/test_smirks_reactions.py
```

### SMIRKS Test Categories

| Category | Tests | What it validates |
|----------|-------|-------------------|
| Parsing | 8 tests | SMIRKS pattern parsing, atom maps, validation |
| Acid-base | 3 tests | Deprotonation, protonation reactions |
| Redox | 3 tests | Oxidation and reduction patterns |
| Substitution | 2 tests | Aromatic halogenation, nitration |
| Hydrolysis | 2 tests | Ester and amide cleavage |
| Edge cases | 5 tests | Multi-component, stereochemistry, errors |

See [SMIRKS_TESTING.md](SMIRKS_TESTING.md) for detailed documentation.

### Reaction-layer comparison harness

`scripts/compare_reaction_layers.py` splits the reaction benchmark into independent layers instead of collapsing everything into a single pass/fail number.

| Layer | Reference / target | Purpose |
|------|---------------------|---------|
| Semantics | RDKit | Compare SMIRKS interpretation and transform agreement |
| Path / TS | sci-form NEB (any method) vs geomeTRIC ElasticBand vs ASE NEB (tblite) | Compare internal vs external path optimization on curated, aligned endpoint pairs |
| Energetics | Configurable methods on common geometries | Smoke-test comparison on the same aligned endpoint geometries |
| Cross-validation | sci-form GFN1/GFN2 vs tblite GFN1/GFN2, PySCF HF/STO-3G | Truly external single-point energy validation on identical geometries |

Available NEB backend methods: `uff`, `mmff94`, `pm3`, `xtb` (GFN0), `gfn1`, `gfn2`, `hf3c`.
Methods with analytical gradients (fast): uff, mmff94, pm3, xtb.
Methods with numerical gradients (slow, for single-point use): gfn1, gfn2, hf3c.

The script uses the experimental CLI commands:

```bash
# Single-point energy with any backend
cargo run -p sci-form-cli --features alpha-gsm,alpha-kinetics -- \
    neb-energy "<smiles>" "<coords_json>" --method pm3

# Energy + gradient with any backend
cargo run -p sci-form-cli --features alpha-gsm,alpha-kinetics -- \
    neb-gradient "<smiles>" "<coords_json>" --method uff

# NEB path with configurable backend
cargo run -p sci-form-cli --features alpha-gsm,alpha-kinetics -- \
    simplified-neb-path "<smiles>" "<start_json>" "<end_json>" --method pm3 --n-images 5 --n-iter 20

# Backend planning
cargo run -p sci-form-cli --features alpha-gsm,alpha-kinetics -- \
    gsm-backend-plan "[H][H]"

# Common-geometry backend smoke comparison used by the harness energy layer:
cargo run -p sci-form-cli --features alpha-gsm,alpha-kinetics -- \
    gsm-compare-backends "<smiles>" "<coords_json>" --methods '["uff","pm3"]'

# GSM+MBH+HTST remains available for deeper investigation, but is no longer the default smoke-test path layer:
cargo run -p sci-form-cli --features alpha-gsm,alpha-kinetics -- \
    analyze-gsm-mbh-htst-step "<smiles>" "<reactant_coords_json>" "<product_coords_json>" --method xtb
```

---

## 5. Python Bindings (`crates/python/`)

```bash
# Build a wheel for the current platform
cd crates/python
maturin build --release

# Install the built wheel
pip install ../../target/wheels/sci_form-*.whl --force-reinstall

# Smoke test
python - <<'EOF'
import sci_form
coords = sci_form.embed("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", 42)  # caffeine
assert len(coords) == 24 * 3, f"unexpected length {len(coords)}"
print("Python smoke test passed — caffeine has", len(coords) // 3, "atoms")
EOF

# Full Python integration tests (requires installed wheel)
python tests/integration/test_python_integration.py
```

---

## 6. WebAssembly / Node.js (`crates/wasm/`)

```bash
# Build for Node.js (CommonJS) — sequential, no parallelization
cd crates/wasm
wasm-pack build --target nodejs --release --out-dir pkg-node

# Smoke test via Node.js
node - <<'EOF'
const sci = require("./pkg-node/sci_form_wasm.js");
const result = JSON.parse(sci.embed("CCO", 42));   // embed returns JSON string
console.assert(!result.error, result.error);
console.assert(result.num_atoms === 9, `expected 9, got ${result.num_atoms}`);
console.log("Node smoke test passed — ethanol has", result.num_atoms, "atoms");
EOF

# Build for browsers and modern dev servers (ESM / Vite) — WITH parallelization
wasm-pack build --target web --release --out-dir pkg-web --features parallel

# Full JS integration tests
node tests/integration/test_wasm_integration.js
```

---

## 7. CLI (`crates/cli/`)

```bash
# Build the CLI binary
cargo build --release --package sci-form-cli

# Smoke tests
./target/release/sci-form --version
./target/release/sci-form embed "CCO"
./target/release/sci-form parse "c1ccccc1"
```

---

## 7. Full Pre-release Check (local replica of CI)

Run this sequence to verify the entire release pipeline locally before tagging:

```bash
#!/usr/bin/env bash
set -euo pipefail

# 1. Quality gates
cargo clippy --all-targets -- -D warnings
cargo fmt --check

# 2. Tests
cargo test --lib
cargo test --tests --release

# 3. Python
cd crates/python
maturin build --release
pip install ../../target/wheels/sci_form-*.whl --force-reinstall
python - <<'EOF'
import sci_form
coords = sci_form.embed("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", 42)
assert len(coords) == 24 * 3
print("Python OK")
EOF
cd -

# 4. WASM / Node
cd crates/wasm
wasm-pack build --target nodejs --release --out-dir pkg-node
node - <<'EOF'
const sci = require("./pkg-node/sci_form_wasm.js");
const result = JSON.parse(sci.embed("CCO", 42));
console.assert(result.num_atoms === 9);
console.log("Node OK");
EOF
cd -

# 5. CLI
cargo build --release --package sci-form-cli
./target/release/sci-form --version
./target/release/sci-form embed "CCO" > /dev/null
echo "CLI OK"

echo ""
echo "All checks passed — ready to tag."
```

---

## 9. Computational Chemistry Validation and Ground Truth Testing (EHT Pipeline)

Validating a computational chemistry engine is the most critical development step. To ensure sci-form calculates quantum-mechanical tensors correctly, you must establish a "ground truth" using established reference tools and compare mathematical matrices step-by-step, not just visual output.

### Validation Methodology

#### 8.1 Select Reference Software

Generate exact reference data using tools that allow automated Python extraction:

| Tool | Purpose | Installation |
|------|---------|---------------|
| **PySCF** | Compute overlap matrix $S$ with minimal basis (STO-3G). Verify Slater integral and basis-set correctness. | `pip install pyscf` |
| **xtb** | Optimize molecular geometry with GFN-xTB or use as secondary 3D reference. | `pip install xtb` |
| **RDKit** | Generate and standardize 3D geometries. Ensures both systems start from identical coordinates. | `pip install rdkit` |
| **NumPy/SciPy** | Implement clean reference EHT logic for direct matrix comparison. | `pip install numpy scipy` |
| **Rust: `approx`** | Compare floating-point results with configurable tolerances. | Add to `Cargo.toml`: `approx = "0.5"` |

#### 8.2 What Variables to Compare (Strict Assertion Order)

If the final calculation fails, it is difficult to pinpoint the error. You must make assertions in this strict execution order:

1. **Input Geometry (XYZ Coordinates)**
   - Assert that atomic coordinates and internuclear distances match exactly.
   - Tolerance: ≤ 0.001 Ångström (difference will drastically change integrals).
   - Use `assert_relative_eq!` with `epsilon = 1e-6` in Rust tests.

2. **Overlap Matrix ($S$)**
   - Compare term-by-term: $S_{ij}$ from sci-form vs. reference (PySCF).
   - If divergence here → problem in distance calculations or basis-function integral evaluation.
   - Tolerance: ≤ 1 × 10⁻⁵ (relative error).

3. **Hamiltonian Matrix ($H$)**
   - Verify diagonal elements match exactly your VSIP (Valence State Ionization Potential) table values.
   - Verify off-diagonal elements use the correct Wolfsberg-Helmholtz formula:

   $$H_{ij} = \frac{1}{2} K S_{ij} \left(H_{ii} + H_{jj}\right)$$

   - Tolerance: ≤ 1 × 10⁻⁵.

4. **Eigenvalues ($E$, Orbital Energies)**
   - Compare HOMO (Highest Occupied Molecular Orbital) and LUMO (Lowest Unoccupied) values.
   - These are the most chemically critical values.
   - Tolerance: ≤ 1 × 10⁻⁴ (often tighter than matrix elements due to conditioning).

5. **Eigenvectors ($C$, Molecular Orbital Coefficients)**
   - **Algorithmic warning:** Eigensolvers (both Rust nalgebra and Python SciPy) can return eigenvectors with opposite sign (phase flip). Both are mathematically correct because orbital probability density $|\Psi|^2$ is identical.
   - **Solution:** Compare absolute values: $|C_{\text{Rust}}| \approx |C_{\text{Ref}}|$.
   - Tolerance: ≤ 1 × 10⁻⁴.

#### 8.3 Testing Pipeline Structure: Rust + JSON Ground Truth

The most robust and professional approach is to create a static data bridge using unit tests.

**Step A (Reference Generation):**
```bash
# Generate reference matrices for H₂O with STO-3G basis
python scripts/generate_eht_reference.py --molecule h2o --basis sto-3g --output tests/data/h2o_ref.json
python scripts/generate_eht_reference.py --molecule ethane --basis sto-3g --output tests/data/ethane_ref.json
```

The Python script exports: XYZ coordinates, $S$, $H$, $E$, and $C$ matrices in a single JSON file.

**Step B (Rust Unit Test Example):**
```rust
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use serde_json::json;
    use std::fs;

    #[derive(serde::Deserialize)]
    struct EHTReference {
        #[serde(rename = "xyz_coords")]
        xyz: Vec<Vec<f64>>,
        #[serde(rename = "overlap_matrix")]
        s_ref: Vec<Vec<f64>>,
        #[serde(rename = "hamiltonian_matrix")]
        h_ref: Vec<Vec<f64>>,
        #[serde(rename = "eigenvalues")]
        e_ref: Vec<f64>,
        #[serde(rename = "eigenvectors")]
        c_ref: Vec<Vec<f64>>,
    }

    fn load_reference(filename: &str) -> EHTReference {
        let data = fs::read_to_string(filename)
            .expect("Could not read reference file");
        serde_json::from_str(&data)
            .expect("Invalid JSON in reference file")
    }

    #[test]
    fn test_overlap_matrix_h2o() {
        let refdata = load_reference("tests/data/h2o_ref.json");
        
        // Create molecule from reference coordinates
        let molecule = Molecule::from_xyz(&refdata.xyz);
        
        // Compute sci-form overlap matrix
        let s_rust = compute_overlap_matrix(&molecule);
        
        // Element-by-element comparison
        for i in 0..s_rust.nrows() {
            for j in 0..s_rust.ncols() {
                assert_relative_eq!(
                    s_rust[(i, j)],
                    refdata.s_ref[i][j],
                    epsilon = 1e-5,
                    max_relative = 1e-5
                );
            }
        }
    }

    #[test]
    fn test_hamiltonian_matrix_ethane() {
        let refdata = load_reference("tests/data/ethane_ref.json");
        let molecule = Molecule::from_xyz(&refdata.xyz);
        let s_matrix = compute_overlap_matrix(&molecule);
        
        // Compute Hamiltonian
        let h_rust = compute_hamiltonian_matrix(&molecule, &s_matrix, 1.75);  // K=1.75
        
        // Compare
        for i in 0..h_rust.nrows() {
            for j in 0..h_rust.ncols() {
                assert_relative_eq!(
                    h_rust[(i, j)],
                    refdata.h_ref[i][j],
                    epsilon = 1e-5,
                    max_relative = 1e-5
                );
            }
        }
    }

    #[test]
    fn test_eigenvalues_homo_lumo() {
        let refdata = load_reference("tests/data/h2o_ref.json");
        let molecule = Molecule::from_xyz(&refdata.xyz);
        let s_matrix = compute_overlap_matrix(&molecule);
        let h_matrix = compute_hamiltonian_matrix(&molecule, &s_matrix, 1.75);
        
        // Solve generalized eigenproblem: HC = SCE
        let (e_rust, _c_rust) = solve_generalized_eigenproblem(&h_matrix, &s_matrix);
        
        // Sort both vectors (eigensolvers may return in different order)
        let mut e_sorted = e_rust.clone();
        e_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut e_ref_sorted = refdata.e_ref.clone();
        e_ref_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        
        // Check HOMO and LUMO specifically (highest occupied and lowest unoccupied)
        // For H₂O with STO-3G (10 basis functions, 5 occupied orbitals)
        let homo_idx = 4;  // 5th eigenvalue (0-indexed)
        let lumo_idx = 5;  // 6th eigenvalue
        
        assert_relative_eq!(
            e_sorted[homo_idx],
            e_ref_sorted[homo_idx],
            epsilon = 1e-4,
            max_relative = 1e-4
        );
        assert_relative_eq!(
            e_sorted[lumo_idx],
            e_ref_sorted[lumo_idx],
            epsilon = 1e-4,
            max_relative = 1e-4
        );
    }

    #[test]
    fn test_eigenvectors_absolute_value() {
        let refdata = load_reference("tests/data/h2o_ref.json");
        let molecule = Molecule::from_xyz(&refdata.xyz);
        let s_matrix = compute_overlap_matrix(&molecule);
        let h_matrix = compute_hamiltonian_matrix(&molecule, &s_matrix, 1.75);
        
        let (_e_rust, c_rust) = solve_generalized_eigenproblem(&h_matrix, &s_matrix);
        
        // Compare absolute values (handles sign flip)
        for i in 0..c_rust.nrows() {
            for j in 0..c_rust.ncols() {
                assert_relative_eq!(
                    c_rust[(i, j)].abs(),
                    refdata.c_ref[i][j].abs(),
                    epsilon = 1e-4,
                    max_relative = 1e-4
                );
            }
        }
    }
}
```

**Step C (Run the validation suite):**
```bash
# All EHT ground-truth tests
cargo test --release --test regression test_eht -- --nocapture

# Individual tests
cargo test --release --test regression test_overlap_matrix_h2o -- --nocapture
cargo test --release --test regression test_hamiltonian_matrix_ethane -- --nocapture
cargo test --release --test regression test_eigenvalues_homo_lumo -- --nocapture
cargo test --release --test regression test_eigenvectors_absolute_value -- --nocapture
```

#### 8.4 Integrating with CI

Add to `.github/workflows/ci.yml`:
```yaml
- name: Generate EHT reference data
  run: python scripts/generate_eht_reference.py --all

- name: Validate EHT pipeline against ground truth
    run: cargo test --release --test regression test_eht -- --nocapture
```

This ensures that any regression in the EHT tensor pipeline is caught before merge.

---

## 10. Validation Protocol: Population Analysis (Mulliken & Löwdin Charges)

Once EHT orbitals are computed, population analysis extracts per-atom partial charges. This is critical for understanding molecular electrostatics, reactivity sites, and chemical hardness.

### Ground-Truth Generation

```bash
# Generate reference Mulliken and Löwdin charges using PySCF
python scripts/generate_population_reference.py \
  --molecules h2o ethane benzene \
  --basis sto-3g \
  --output tests/data/population_ref.json
```

The reference JSON should contain:
```json
{
  "h2o": {
    "xyz_coords": [[0,0,0], [0.757,0.586,0], [-0.757,0.586,0]],
    "mulliken_charges": [-0.34, 0.17, 0.17],
    "lowdin_charges": [-0.22, 0.11, 0.11],
    "sum_mulliken": 0.0,
    "sum_lowdin": 0.0
  }
}
```

### Rust Test Structure

```rust
#[cfg(test)]
mod population_tests {
    use approx::assert_relative_eq;

    #[test]
    fn test_mulliken_charges_h2o() {
        let refdata = load_reference("tests/data/population_ref.json");
        let (elems, pos) = h2o_molecule();
        
        // Compute orbitals → density matrix → Mulliken charges
        let result = solve_eht(&elems, &pos, None).unwrap();
        let charges_mulliken = compute_mulliken_charges(&result, &elems);
        
        for i in 0..charges_mulliken.len() {
            assert_relative_eq!(
                charges_mulliken[i],
                refdata["h2o"]["mulliken_charges"][i],
                epsilon = 1e-4,
                max_relative = 1e-4
            );
        }
        
        // Charge sum must equal molecular total charge (0 for neutral)
        let sum: f64 = charges_mulliken.iter().sum();
        assert!(sum.abs() < 1e-6, "Mulliken charges must sum to 0");
    }

    #[test]
    fn test_lowdin_charges_h2o() {
        let refdata = load_reference("tests/data/population_ref.json");
        let (elems, pos) = h2o_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        let charges_lowdin = compute_lowdin_charges(&result, &elems);
        
        for i in 0..charges_lowdin.len() {
            assert_relative_eq!(
                charges_lowdin[i],
                refdata["h2o"]["lowdin_charges"][i],
                epsilon = 1e-4,
                max_relative = 1e-4
            );
        }
        
        let sum: f64 = charges_lowdin.iter().sum();
        assert!(sum.abs() < 1e-6);
    }

    #[test]
    fn test_charge_invariants_benzene() {
        let (elems, pos) = benzene_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        
        // For C₆H₆, all carbons should have identical Mulliken charges (symmetry)
        let charges = compute_mulliken_charges(&result, &elems);
        let carbon_charges: Vec<f64> = charges[..6].to_vec();
        let mean = carbon_charges.iter().sum::<f64>() / 6.0;
        
        for &q in &carbon_charges {
            assert!(
                (q - mean).abs() < 0.001,
                "Benzene carbons should be symmetric"
            );
        }
    }
}
```

---

## 11. Validation Protocol: Dipole Moments

A molecular dipole $\vec{\mu} = (\mu_x, \mu_y, \mu_z)$ is fundamental for understanding polarity and reactivity.

### Ground-Truth Generation

```bash
python scripts/generate_dipole_reference.py \
  --molecules h2o methane benzene formaldehyde \
  --basis sto-3g \
  --output tests/data/dipole_ref.json
```

Reference format:
```json
{
  "h2o": {
    "dipole_vector": [0.0, 0.0, 1.85],
    "dipole_magnitude": 1.85,
    "unit": "Debye"
  }
}
```

### Rust Test Structure

```rust
#[cfg(test)]
mod dipole_tests {
    use approx::assert_relative_eq;
    use nalgebra::Vector3;

    #[test]
    fn test_dipole_magnitude_h2o() {
        let refdata = load_reference("tests/data/dipole_ref.json");
        let (elems, pos) = h2o_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        
        let dipole_vec = compute_dipole_moment(&result, &elems, &pos);
        let magnitude = dipole_vec.norm();
        
        // Water dipole should be ~1.85 D
        assert_relative_eq!(
            magnitude,
            refdata["h2o"]["dipole_magnitude"],
            epsilon = 0.05,
            max_relative = 0.05
        );
    }

    #[test]
    fn test_dipole_zero_for_symmetric_molecules() {
        let (elems, pos) = methane_molecule(); // CH₄ is nonpolar
        let result = solve_eht(&elems, &pos, None).unwrap();
        let dipole_vec = compute_dipole_moment(&result, &elems, &pos);
        
        assert!(
            dipole_vec.norm() < 0.1,
            "Methane should be nonpolar"
        );
    }

    #[test]
    fn test_dipole_direction_h2o() {
        let refdata = load_reference("tests/data/dipole_ref.json");
        let (elems, pos) = h2o_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        let dipole_vec = compute_dipole_moment(&result, &elems, &pos);
        
        let ref_vec = Vector3::new(
            refdata["h2o"]["dipole_vector"][0],
            refdata["h2o"]["dipole_vector"][1],
            refdata["h2o"]["dipole_vector"][2],
        );
        
        // Angle between vectors should be < 5°
        let cos_angle = (dipole_vec.normalize().dot(&ref_vec.normalize())).min(1.0);
        let angle_rad = cos_angle.acos();
        assert!(angle_rad < 0.087, "Dipole direction wrong by {:.1}°", angle_rad.to_degrees());
    }
}
```

---

## 12. Validation Protocol: Electrostatic Potential Maps (ESP)

ESP is a 3D volumetric property showing the electrostatic potential at every point, essential for visualization.

### Ground-Truth Generation

```bash
# PySCF cubegen or Gaussian cube export
python scripts/generate_esp_cube.py \
  --molecule h2o \
  --basis sto-3g \
  --grid-spacing 0.2 \
  --output tests/data/h2o_esp.cube
```

Cube file format (standard):
```
 Cube file header
   Generated by PySCF
    3    0.0    0.0    0.0
   20    0.2    0.0    0.0
   30    0.0    0.2    0.0
   25    0.0    0.0    0.2
 8    0    0.0    0.0    0.0
 1    0    0.96   0.0    0.0
...
```

### Rust Test Structure

```rust
#[cfg(test)]
mod esp_tests {
    use approx::assert_relative_eq;

    #[test]
    fn test_esp_grid_generation_h2o() {
        let (elems, pos) = h2o_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        
        // Generate ESP on a 20×20×20 grid, 0.2 Å spacing
        let esp_grid = compute_esp_grid(
            &elems,
            &pos,
            &result,
            0.2,  // spacing
            (20, 20, 20),  // grid dimensions
        );
        
        assert_eq!(esp_grid.len(), 20 * 20 * 20);
        assert!(esp_grid.iter().all(|&v| !v.is_nan()), "No NaN values allowed");
    }

    #[test]
    fn test_esp_vs_reference_cube() {
        let reference_cube = read_cube_file("tests/data/h2o_esp.cube");
        let (elems, pos) = h2o_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        
        let esp_rust = compute_esp_grid(&elems, &pos, &result, 0.2, (20, 20, 20));
        
        // Compute element-wise error
        let mut max_error = 0.0;
        let mut sum_sqerror = 0.0;
        for i in 0..esp_rust.len() {
            let err = (esp_rust[i] - reference_cube.values[i]).abs();
            max_error = max_error.max(err);
            sum_sqerror += err * err;
        }
        let rmse = (sum_sqerror / esp_rust.len() as f64).sqrt();
        
        assert!(rmse < 1e-4, "ESP RMSE too high: {}", rmse);
        assert!(max_error < 0.01, "Max ESP error: {}", max_error);
    }

    #[test]
    fn test_esp_export_cube_format() {
        let (elems, pos) = h2o_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        
        // Export to .cube file
        let esp_grid = compute_esp_grid(&elems, &pos, &result, 0.2, (20, 20, 20));
        export_esp_cube(
            "target/test_esp_export.cube",
            &elems,
            &pos,
            &esp_grid,
            0.2,
            (20, 20, 20),
        ).unwrap();
        
        // Re-read and verify format
        let reread = read_cube_file("target/test_esp_export.cube");
        assert_eq!(reread.values.len(), 20 * 20 * 20);
    }
}
```

---

## 13. Validation Protocol: Density of States (DOS) and PDOS

DOS shows the number of states at each energy level; PDOS breaks this down by atom.

### Ground-Truth Generation

```bash
# Multiwfn or custom Python script
python scripts/generate_dos_reference.py \
  --molecule benzene \
  --basis sto-3g \
  --smearing 0.2 \
  --energies '-20,-15,-10,-5,0,5,10,15,20' \
  --output tests/data/dos_ref.json
```

Reference format:
```json
{
  "benzene": {
    "energies": [-20, -15, -10, ...],
    "dos_total": [0.001, 0.002, 0.05, ...],
    "pdos_carbon": [0.0001, 0.0002, 0.025, ...],
    "pdos_hydrogen": [0.00001, 0.00002, 0.025, ...]
  }
}
```

### Rust Test Structure

```rust
#[cfg(test)]
mod dos_tests {
    use approx::assert_relative_eq;

    #[test]
    fn test_dos_curve_generation() {
        let (elems, pos) = benzene_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        
        // Generate DOS with Gaussian smearing σ=0.2 eV
        let energies = linspace(-20.0, 20.0, 100);
        let dos_total = compute_dos(&result, &energies, 0.2);
        
        // DOS must be non-negative
        assert!(dos_total.iter().all(|&v| v >= 0.0));
        
        // DOS should integrate roughly to number of orbitals
        let integral: f64 = dos_total.iter().sum::<f64>() * (40.0 / 100.0);  // dE = 0.4
        assert!(
            (integral - (result.energies.len() as f64)).abs() < 1.0,
            "DOS integral should ≈ num_orbitals"
        );
    }

    #[test]
    fn test_pdos_vs_reference() {
        let refdata = load_reference("tests/data/dos_ref.json");
        let (elems, pos) = benzene_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        
        let energies = linspace(-20.0, 20.0, 100);
        let pdos_carbon = compute_pdos(&result, &elems, 6, &energies, 0.2);  // element 6 = C
        
        let ref_pdos = &refdata["benzene"]["pdos_carbon"];
        for i in 0..pdos_carbon.len() {
            assert_relative_eq!(
                pdos_carbon[i],
                ref_pdos[i],
                epsilon = 0.01,
                max_relative = 0.1
            );
        }
    }

    #[test]
    fn test_dos_peak_positions_match_eigenvalues() {
        let (elems, pos) = h2o_molecule();
        let result = solve_eht(&elems, &pos, None).unwrap();
        
        // With very small smearing, DOS peaks should align with eigenvalues
        let energies = linspace(-30.0, 10.0, 1000);
        let dos = compute_dos(&result, &energies, 0.01);  // σ=0.01 eV
        
        // For each eigenvalue, find nearest maximum in DOS
        for &eig in &result.energies {
            let idx = ((eig + 30.0) / 40.0 * 1000.0) as usize;
            assert!(
                idx > 0 && idx < dos.len(),
                "Eigenvalue {} out of energy range",
                eig
            );
            // DOS should have a local maximum near this energy
            assert!(
                dos[idx] > 0.1 * dos.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
                "DOS peak too low at eigenvalue"
            );
        }
    }
}
```

---

## 14. Validation Protocol: Molecular Alignment and RMSD (Kabsch)

RMSD is the industry standard for comparing 3D structures. Kabsch algorithm finds optimal superposition.

### Ground-Truth Generation

```bash
# Use RDKit as reference for RMSD computation
python scripts/generate_rmsd_reference.py \
  --molecules "CC(C)C,CC(CC)C,CCCC" \
  --output tests/data/rmsd_ref.json
```

Reference format:
```json
{
  "CC(C)C_vs_CC(CC)C": {
    "rmsd_heavy_atom": 0.892,
    "rmsd_all_atoms": 1.234,
    "rmsd_hydrogens_only": 1.876
  }
}
```

### Rust Test Structure

```rust
#[cfg(test)]
mod rmsd_tests {
    use approx::assert_relative_eq;

    #[test]
    fn test_kabsch_alignment_identity() {
        // Aligning a molecule to itself should give RMSD = 0
        let (elems, pos) = ethane_molecule();
        let rmsd = compute_rmsd_kabsch(&pos, &pos, &elems);
        
        assert!(rmsd < 1e-10, "Self-alignment RMSD should be ~0");
    }

    #[test]
    fn test_rmsd_vs_rdkit_reference() {
        let refdata = load_reference("tests/data/rmsd_ref.json");
        let mol1 = embed_molecule("CC(C)C", 42);
        let mol2 = embed_molecule("CC(CC)C", 42);
        
        let rmsd_heavy = compute_rmsd_kabsch(
            &mol1.coords,
            &mol2.coords,
            &mol1.elements,
            Some(true),  // heavy atoms only
        );
        
        assert_relative_eq!(
            rmsd_heavy,
            refdata["CC(C)C_vs_CC(CC)C"]["rmsd_heavy_atom"],
            epsilon = 0.01,
            max_relative = 0.05
        );
    }

    #[test]
    fn test_rmsd_invariant_under_rotation() {
        let (elems, pos) = water_molecule();
        
        // Rotate molecule by 45° around Z axis
        let rotation_matrix = rotation_z_45deg();
        let rotated = apply_rotation(&pos, &rotation_matrix);
        
        // RMSD after alignment should still be 0
        let rmsd = compute_rmsd_kabsch(&pos, &rotated, &elems);
        assert!(rmsd < 1e-10, "Rotation should not affect aligned RMSD");
    }

    #[test]
    fn test_rmsd_hydrogen_vs_heavy_atom() {
        let (elems, pos) = ethanol_molecule();
        let (elems_rot, pos_rot) = rotate_molecule(&(elems.clone(), pos.clone()), 10.0);
        
        let rmsd_all = compute_rmsd_kabsch(&pos, &pos_rot, &elems, None);15        let rmsd_heavy = compute_rmsd_kabsch(&pos, &pos_rot, &elems, Some(true));
        
        // Heavy-atom RMSD should typically be smaller (hydrogens more flexible)
        assert!(rmsd_heavy <= rmsd_all);
    }

    #[test]
    fn test_quaternion_vs_matrix_alignment_equivalence() {
        let (elems, pos) = benzene_molecule();
        let (_, pos_rot) = rotate_molecule(&(elems.clone(), pos.clone()), 5.0);
        
        // Compute RMSD using two different alignment methods
        let rmsd_kabsch = compute_rmsd_kabsch(&pos, &pos_rot, &elems);
        let rmsd_quat = compute_rmsd_quaternion(&pos, &pos_rot, &elems);
        
        // Both should agree within numerical tolerance
        assert_relative_eq!(
            rmsd_kabsch,
            rmsd_quat,
            epsilon = 1e-8,
            max_relative = 1e-8
        );
    }
}
```

---

## 15. Validation Protocol: Force Fields (UFF / MMFF94)

Force fields are parameterized to match experimental and ab initio data. Gradients must match numerical finite differences.

### Ground-Truth Generation

```bash
# Export UFF/MMFF gradients from OpenBabel or RDKit
python scripts/generate_ff_reference.py \
  --molecules "C,CCO,c1ccccc1" \
  --force-field uff \
  --output tests/data/ff_uff_ref.json
```

Reference format:
```json
{
  "C": {
    "energy": 0.0,
    "gradients": [[0,0,0], [0,0,0], [0,0,0], [0,0,0], [0,0,0]],
    "bonds": [{"a": 0, "b": 1, "order": "single", "length": 1.09}],
    "angles": [{"a": 1, "b": 0, "c": 2, "angle": 109.47}],
    "torsions": [{"a": 0, "b": 1, "c": 2, "d": 3, "phi": 180.0}]
  }
}
```

### Rust Test Structure

```rust
#[cfg(test)]
mod force_field_tests {
    use approx::assert_relative_eq;

    #[test]
    fn test_uff_energy_methane() {
        let (elems, pos) = methane_molecule();
        let energy = compute_uff_energy(&elems, &pos);
        
        // Optimized methane should have very low strain energy
        assert!(energy < 0.1, "Methane strain energy too high: {}", energy);
    }

    #[test]
    fn test_uff_gradients_vs_numerical() {
        let (elems, pos) = water_molecule();
        let gradients_analytical = compute_uff_gradients(&elems, &pos);
        
        let delta = 1e-6;
        let gradients_numerical = compute_numerical_gradients(
            &elems,
            &pos,
            delta,
            |e, p| compute_uff_energy(e, p),
        );
        
        for i in 0..gradients_analytical.len() {
            assert_relative_eq!(
                gradients_analytical[i],
                gradients_numerical[i],
                epsilon = 1e-4,
                max_relative = 1e-4
            );
        }
    }

    #[test]
    fn test_uff_vs_reference_gradients() {
        let refdata = load_reference("tests/data/ff_uff_ref.json");
        let (elems, pos) = embed_molecule("CCO", 42);
        
        let grads = compute_uff_gradients(&elems, &pos);
        let ref_grads = &refdata["CCO"]["gradients"];
        
        for i in 0..grads.len() {
            assert_relative_eq!(
                grads[i],
                ref_grads[i],
                epsilon = 0.001,
                max_relative = 0.01
            );
        }
    }

    #[test]
    fn test_mmff94_bond_parameters() {
        let refdata = load_reference("tests/data/ff_uff_ref.json");
        let (elems, pos) = ethane_molecule();
        
        // MMFF94 should have bond parameters: length and strength
        let bond_terms = compute_mmff94_bond_terms(&elems, &pos);
        
        // For C-H bond in ethane: ~1.09 Å
        let ch_bond = bond_terms.iter().find(|b| is_ch_bond(b)).unwrap();
        assert!((ch_bond.length - 1.09).abs() < 0.01, "C-H bond length off");
    }

    #[test]
    fn test_torsion_terms_benzene() {
        let (elems, pos) = benzene_molecule();
        
        let torsion_terms = compute_uff_torsion_terms(&elems, &pos);
        
        // Benzene should have planar torsions (0° or 180°)
        for torsion in torsion_terms {
            let angle_mod = torsion.phi % 180.0;
            assert!(
                angle_mod.abs() < 5.0 || (angle_mod - 180.0).abs() < 5.0,
                "Benzene should be planar"
            );
        }
    }

    #[test]
    fn test_nonbonded_vdw_repulsion() {
        let mut pos = vec![[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]];  // Very close
        let elems = vec![6, 6];  // Two carbons
        
        let energy_close = compute_uff_energy(&elems, &pos);
        
        // Move atoms farther apart
        pos[1] = [3.5, 0.0, 0.0];  // vdW distance ~3.4 Å
        let energy_far = compute_uff_energy(&elems, &pos);
        
        // Close configuration should have much higher energy (vdW repulsion)
        assert!(
            energy_close > energy_far * 10.0,
            "vdW repulsion not working: {} vs {}",
            energy_close,
            energy_far
        );
    }
}
```

---

## Summary: All-Phase Validation Matrix

| Phase | Reference Tool | What to Compare | Tolerance | Test Location |
|-------|-----------------|-----------------|-----------|---|
| A (Conformer) | RDKit ETKDG | Heavy-atom RMSD | 0.5 Å | `regression/test_geometry_quality.rs` |
| B (EHT) | PySCF/Multiwfn | S, H, E, \|C\| | 1e-4 | `regression/test_eht.rs` |
| C1 (Charges) | PySCF Gasteiger | Partial charges | 1e-4 | `regression/test_charges.rs` |
| C1 (SASA) | Shrake-Rupley ref | SASA area | ±1 Ų | `regression/test_sasa.rs` |
| C2 (Mulliken) | PySCF Mulliken | Pop. charges | 1e-4 | (Phase C2 candidate) |
| C3 (Dipole) | ORCA/Gaussian | μ magnitude + direction | 0.05 D / 5° | (Phase C3 candidate) |
| C4 (ESP) | Gaussian cubegen | .cube grid RMSE | <1e-4 | (Phase C4 candidate) |
| C5 (DOS) | Multiwfn | DOS curve MSE | <0.01 | (Phase C5 candidate) |
| C6 (RMSD) | RDKit AlignMol | Pairwise RMSD | <0.01 Å | (Phase C6 candidate) |
| C7 (UFF) | OpenBabel/RDKit | Gradients | 1e-4 | (Phase C7 candidate) |

The GitHub Actions release pipeline (`.github/workflows/release.yml`) runs the same gates
automatically when a version tag (`v*`) is pushed:

```
push tag v* 
  └─► pre-release-tests (ubuntu / macos / windows)  ─┐
  └─► lint (clippy + fmt)                             ─┴─► build-cli (5 targets)
                                                       ─► build-python (3 OS)
                                                       ─► build-wasm (bundler + node)
                                                           │
                                                      verify-python
                                                      verify-node
                                                      verify-cli
                                                           │
                                               publish-crate / publish-python / publish-npm
                                                           │
                                                       release (GitHub Release)
```

| Job | What it does |
|-----|--------------|
| `pre-release-tests` | `cargo test --release --test ci -- --nocapture` on ubuntu, macos, windows |
| `lint` | `cargo clippy -D warnings` + `cargo fmt --check` |
| `build-cli` | Cross-compiles for 5 targets (linux x86_64/aarch64, macos x86_64/aarch64, windows x86_64) |
| `build-python` | `maturin build --release` on ubuntu, macos, windows |
| `build-wasm` | `wasm-pack build` with `bundler` and `nodejs` targets |
| `verify-python` | Installs ubuntu wheel, embeds caffeine, checks atom count |
| `verify-node` | Installs node pkg, embeds ethanol via CJS require |
| `verify-cli` | Downloads linux binary, runs `--version`, `embed`, `parse` |
| `publish-crate` | `cargo publish` for lib + cli (gated on all verifications) |
| `publish-python` | `maturin upload` to PyPI |
| `publish-npm` | `wasm-pack publish` to npm |
| `release` | Creates GitHub Release with changelog + install instructions |

---

## Experimental Engine Tests (quantum chemistry regression)

Two dedicated regression suites validate the new `src/experimental_2/` Roothaan-Hall RHF engine against NIST data and legacy methods.

### `tests/regression/test_experimental_comparison.rs` — 54 tests

| Module | Count | What it tests |
|--------|-------|---------------|
| `experimental_scf_correctness` | 8 | Convergence, negative total energy, positive gap, finite Mulliken charges (H₂, HeH⁺, H₂O, CH₄, NH₃, HF, CO, C₂H₄) |
| `legacy_vs_experimental_energy` | 8 | Side-by-side EHT/PM3/xTB vs Exp SCF HOMO-LUMO gaps |
| `legacy_vs_experimental_charges` | 5 | Mulliken charge polarity across methods |
| `experimental_vs_nist` | 5 | Qualitative comparison vs NIST STO-3G reference energies |
| `quantum_engine_integrals` | 8 | Overlap matrix: diagonal > 0, off-diagonal finite, Hermitian symmetry |
| `spectroscopy_comparison` | 5 | sTDA UV-Vis and GIAO NMR legacy vs experimental |
| `timing_benchmarks` | 5 | All 5 methods timing table |
| `comprehensive_energy_table` | 10 | 11-molecule all-methods table |

```bash
cargo test --test regression test_experimental_comparison
# → test result: ok. 54 passed; 0 failed
```

### `tests/regression/test_extended_molecules.rs` — 21 tests (14 active, 7 ignored)

Extended battery with more complex molecules, NIST experimental property validation, and parallel ERI benchmarks.

**Active tests (fast, always run):**

| Test | Validates |
|------|-----------|
| `test_dipole_moments_vs_nist` | EHT dipole for 11 molecules vs NIST (water, methanol, pyridine, benzene, etc.) |
| `test_mulliken_charges_heteroatom_polarity` | PM3/xTB charges; xTB conservation; O-atom polarity |
| `test_gaps_ordering_aromatic_series` | EHT/PM3/xTB/Exp gaps for 6 aromatics; all gaps > 0 |
| `test_legacy_nmr_benzene_reference_exact` | ¹H ≈ 7.27 ppm, ¹³C ≈ 128.5 ppm (SDBS) |
| `test_legacy_nmr_carbonyl_downfield` | Carbonyl ¹³C > 150 ppm (acetone, benzaldehyde, formaldehyde) |
| `test_legacy_nmr_shifts_aromatic_series` | Aromatic ¹H 4.5–10.5 ppm, ¹³C 60–180 ppm |
| `test_legacy_stda_aromatic_series` | sTDA non-empty excitations with positive energies |
| `test_excitation_ordering_aromatic_series` | Naphthalene first exc. ≤ benzene + 1.5 eV |
| `test_ir_spectrum_extended_molecules` | EHT vibrational: ≥ 3N-6 modes for 6 molecules |
| `test_ml_logp_extended_molecules` | Crippen logP vs experiment for 11 molecules; ≥ 60% within tolerance |
| `test_ml_lipinski_druglike_molecules` | Lipinski Ro5 for aspirin, caffeine, ibuprofen, paracetamol, cyclosporin A |
| `test_fingerprint_tanimoto_similarity_series` | ECFP4 Tanimoto; benzene↔toluene > benzene↔methanol; self-sim = 1.0 |
| `test_druglike_molecules_eht_pm3_xtb` | EHT/PM3/xTB for 5 drugs all converge |
| `test_solvation_druglike_molecules` | Non-polar SASA solvation for 5 drug-like molecules |

**Heavy tests (marked `#[ignore]`, run with `--release -- --ignored`):**

| Test | Description |
|------|-------------|
| `bench_parallel_acceleration_table` | Seq vs parallel ERI speedup table: H₂O, NH₃, CH₄, methanol, formaldehyde |
| `bench_speedup_scaling_with_molecule_size` | Speedup vs N: H₂O → CH₄ → C₂H₄ |
| `test_extended_molecules_all_methods_converge` | All 5 methods × 15 molecules convergence |
| `test_homo_lumo_gap_vs_nist_ip_koopmans` | Koopmans' theorem: −ε_HOMO vs 10 NIST IPs |
| `test_giao_nmr_vs_legacy_heteroaromatics` | GIAO NMR vs HOSE codes for pyridine, furan, imidazole |
| `test_druglike_molecules_experimental_scf` | RHF/STO-3G for aspirin, paracetamol, nicotine |
| `test_master_comparison_table_extended` | Master table: all methods × 12 molecules |

```bash
# Fast tests only
cargo test --test regression test_extended_molecules
# → test result: ok. 14 passed; 0 failed; 7 ignored

# All tests including heavy (release for reasonable speed)
cargo test --release --test regression test_extended_molecules -- --include-ignored
```

### Method comparison summary

| Property | EHT | PM3 | GFN-xTB | Exp RHF/STO-3G |
|----------|-----|-----|---------|----------------|
| Theory | Extended Hückel | NDDO semi-empirical | GFN0 tight-binding | Roothaan-Hall HF |
| SCF | No | Yes | Yes | Yes (DIIS) |
| Speed (debug, ~30 atoms) | < 1 ms | 50 ms | 50 ms | 5–30 s |
| Energy unit | eV (orbital) | kcal/mol (HOF) | eV (total) | Hartree (total) |
| GPU | None | None | None | Planned (wgpu stub) |
| CPU parallel | No | No | No | rayon ERI (yes) |

### Parallel ERI acceleration

```rust
// Sequential (default)
let r = run_scf(&system, &ScfConfig::default());

// Parallel via rayon (identical results, faster for N ≥ 20 basis fns)
let r = run_scf(&system, &ScfConfig::parallel());
```

Observed on Intel i5-10500H (12 logical cores):
- N=7 (H₂O): ~0.8–1.1× (overhead dominates at small N)
- N=36 (benzene): ~3–5× speedup
- N=42 (pyridine): ~4–6× speedup

GPU acceleration (GTX 1650, CUDA 13.0) is planned via the `phase1_gpu_infrastructure` wgpu stub.
