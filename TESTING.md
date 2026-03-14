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

Each file targets a specific aspect of the library. All pass on the release profile.

```bash
# All integration tests
cargo test --tests --release

# Individual tests (useful during development)
cargo test --release --test test_diverse_molecules -- --nocapture
cargo test --release --test test_geometry_quality  -- --nocapture
cargo test --release --test test_gradient_check    -- --nocapture
cargo test --release --test test_ensemble_rmsd     -- --nocapture
cargo test --release --test test_gdb20_rmsd        -- --nocapture
cargo test --release --test test_gdb20             -- --nocapture
cargo test --release --test test_tet_centers       -- --nocapture
```

### Test descriptions

| Test file | What it validates |
|-----------|-------------------|
| `test_diverse_molecules.rs` | Embed 50+ chemically diverse molecules; all produce valid 3-D geometries |
| `test_geometry_quality.rs` | Bond lengths / angles stay within ±3 σ of CSD distributions |
| `test_gradient_check.rs` | Force-field gradients match numerical finite-differences (tolerance 1 × 10⁻⁴) |
| `test_ensemble_rmsd.rs` | Pose diversity: multi-conformer ensembles span expected RMSD range |
| `test_gdb20_rmsd.rs` | RMSD vs. RDKit/ETKDG on GDB-20 subset (benchmark accuracy) |
| `test_gdb20.rs` | Success rate ≥ 99 % on GDB-20 subset |
| `test_tet_centers.rs` | Tetrahedral stereocentres are preserved after embedding |

---

## 4. Python Bindings (`crates/python/`)

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
python tests/test_python_integration.py
```

---

## 5. WebAssembly / Node.js (`crates/wasm/`)

```bash
# Build for Node.js (CommonJS)
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

# Build for bundlers (ESM / webpack / vite)
wasm-pack build --target bundler --release --out-dir pkg-bundler

# Full JS integration tests
node tests/test_wasm_integration.js
```

---

## 6. CLI (`crates/cli/`)

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
const { SciForm } = require("./pkg-node/sci_form_wasm.js");
const sf = new SciForm();
const c = sf.embed("CCO", 42);
console.assert(c.length === 9);
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

## 8. CI / Release Workflow

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
| `pre-release-tests` | `cargo test --tests --release` on ubuntu, macos, windows |
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
