> **Historical document** — Snapshot of the unified alpha/beta architecture design. For current feature gate patterns see `.agents/skills/sci-form/SKILL.md`.

# Unified Alpha/Beta Architecture

## Overview

**Goal**: Ensure that all alpha/beta algorithms are accessible from **ALL binding layers** (Python, WASM, CLI, Rust) with consistent naming. When code migrates from alpha → beta → core, only the **import path changes**, not the function names or signatures.

## How It Works

### 1. Core Organization (Rust `src/`)

All alpha and beta algorithms live in dedicated module paths:

```
src/
├── alpha/                    # Early research (alpha stability)
│   ├── edl/                  # Electrochemistry (double layers)
│   │   ├── mod.rs            # Public API
│   │   ├── helmholtz.rs
│   │   ├── gouy_chapman.rs
│   │   └── gcs.rs
│   ├── periodic_linear/      # Periodic systems, k-mesh, band structure
│   ├── kinetics/             # HTST rates, transition state finding (GSM, MBH)
│   └── render_bridge/        # Chart builders for visualization
│
└── beta/                     # (future) Tested, documented, APIs stable
    ├── ...
```

### 2. Binding Layer Access Pattern

Each binding (Python, WASM, CLI) imports from the **same centralized paths**:

#### Python (`crates/python/src/experimental.rs`)
```rust
#[cfg(feature = "alpha-edl")]
pub mod alpha_edl_py {
    use sci_form::alpha::edl::*;  // ← Centralized import
    
    #[pyfunction]
    pub fn edl_profile(...) -> PyResult<...> { ... }
}
```

#### WASM (`crates/wasm/src/experimental.rs`)
```rust
#[cfg(feature = "alpha-edl")]
#[wasm_bindgen]
pub fn compute_edl_profile(...) -> String {
    use sci_form::alpha::edl::*;  // ← Same centralized import
    ...
}
```

#### CLI (`crates/cli/src/experimental_cmds.rs`)
```rust
#[cfg(feature = "alpha-edl")]
pub fn cmd_edl_profile(...) {
    use sci_form::alpha::edl::*;  // ← Same centralized import
    ...
}
```

### 3. Naming Convention

| Layer | Prefix | Example | Why |
|-------|--------|---------|-----|
| **Python** | None | `edl_profile()` | Clean, Pythonic |
| **WASM** | `compute_` | `compute_edl_profile()` | Distinguishes from getters |
| **CLI** | `cmd_` | `cmd_edl_profile` | Subcommand convention |
| **Rust** | None | `edl_profile()` | Direct import from lib |

**Core function name is identical** across all platforms — only **prefix changes** per binding.

## Current State (Session 5)

### Alpha Modules Exported

| Module | Function | Python | WASM | CLI | Rust |
|--------|----------|--------|------|-----|------|
| **edl** | `edl_profile` | ✅ | ✅ | ✅ | ✅ |
| | `edl_capacitance_scan` | ✅ | ✅ | ✅ | ✅ |
| | `edl_helmholtz` | ✅ | ✅ | ✅ | ✅ |
| | `edl_gouy_chapman` | ✅ | ✅ | ✅ | ✅ |
| | `edl_gcs` | ✅ | ✅ | ✅ | ✅ |
| | `edl_chart` | ✅ | ✅ | ✅ | ✅ |
| | `capacitance_chart` | ✅ | ✅ | ✅ | ✅ |
| **periodic_linear** | `kmesh` | ✅ | ✅ | ✅ | ✅ |
| **kinetics** | `htst_rate` | ✅ | ✅ | ✅ | ✅ |
| | `htst_temperature_sweep` | ✅ | ✅ | ❌ | ✅ |
| **render_bridge** | `pack_chart_payload` | ❌ | ✅ | ❌ | ✅ |

### Coverage

- **EDL module**: 100% (7/7 functions across all bindings) ✅
- **Kinetics module**: 83% (2/3 functions, needs CLI htst_temperature_sweep)
- **Periodic module**: 50% (1/11 functions, needs full export suite)
- **Render bridge**: 14% (1/15 functions, needs comprehensive export)
- **Overall alpha**: 31% (11/42 functions) — **Major gap in periodic + render**

## Migration Path (When Code Moves Core)

### Example: EDL Module → Core

**Before (Alpha)**:
```python
from sci_form import edl_profile  # import from edl module
```

**After (Core)**:
```python
from sci_form import edl_profile  # import from core module
```

**Function name unchanged** — only import source changes.

## Feature Flags

### Cargo.toml structure

```toml
[features]
alpha-edl = ["experimental-cpm", "experimental-eeq", "experimental-alpb"]
alpha-render-bridge = ["experimental-kpm", "experimental-gsm", "experimental-gpu-rendering"]
alpha-kinetics = ["experimental-gsm", "experimental-mbh"]
alpha-periodic-linear = ["experimental-kpm", "experimental-randnla"]
```

### Activation

Build with all alpha:
```bash
cargo build --features alpha-edl,alpha-render-bridge,alpha-kinetics,alpha-periodic-linear
```

Python:
```bash
maturin develop --release --features alpha-edl,alpha-render-bridge,alpha-kinetics,alpha-periodic-linear
```

WASM:
```bash
wasm-pack build --features alpha-edl,alpha-render-bridge,alpha-kinetics,alpha-periodic-linear
```

## Adding a New Alpha Function

### 1. Implement in Rust (`src/alpha/{module}/`)

```rust
pub fn my_new_function(...) -> Result<MyOutput, String> {
    // Implementation
}
```

### 2. Export to Python (`crates/python/src/experimental.rs`)

```rust
#[pyfunction]
pub fn my_new_function(...) -> PyResult<...> {
    use sci_form::alpha::{module}::*;
    // Wrapper
}

pub fn register(m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(my_new_function, m)?)?;
    Ok(())
}
```

### 3. Export to WASM (`crates/wasm/src/experimental.rs`)

```rust
#[wasm_bindgen]
pub fn compute_my_new_function(...) -> String {
    use sci_form::alpha::{module}::*;
    // JSON wrapper
}
```

### 4. Export to CLI (`crates/cli/src/experimental_cmds.rs`)

```rust
pub fn cmd_my_new_function(...) {
    use sci_form::alpha::{module}::*;
    // CLI handler
}
```

### 5. Register CLI Command (`crates/cli/src/args.rs` + `main.rs`)

```rust
// args.rs
MyNewFunction { arg1: String, #[arg(...)] arg2: f64 },

// main.rs
Commands::MyNewFunction { arg1, arg2 } => cmd_my_new_function(&arg1, arg2),
```

### 6. Add Feature Gate (`Cargo.toml`)

```toml
#[cfg(feature = "alpha-my-module")]
```

### 7. Test from All Layers

```bash
# Rust
cargo test --lib --features alpha-my-module

# Python
pytest tests/test_my_new_function.py

# WASM
npm test

# CLI
sci-form cmd-my-new-function --help
```

## Consistency Checksum

To verify all exports are consistent across bindings, run:

```bash
# Count functions in each Python
grep -c "^pub fn " crates/python/src/experimental.rs

# Count functions in WASM (all compute_edl_*)
grep -c "^pub fn compute_" crates/wasm/src/experimental.rs

# Count commands in CLI (all cmd_*)
grep -c "^pub fn cmd_" crates/cli/src/experimental_cmds.rs
```

**Expected result**: Same count for each alpha module across all bindings.

## Known Gaps (Session 5)

1. **Render bridge**: 15 functions defined, only 1-2 exported
   - Need: `edl_profile_chart()`, `capacitance_scan_chart()`, `arrhenius_chart()`, `trajectory_chart()`, etc.

2. **Periodic**: 11 functions, only `kmesh()` exported
   - Need: `bz_integration()`, `band_structure_validator()`, `fermi_surface_generator()`, etc.

3. **Kinetics**: 6 functions, only `htst_rate()` and partial `htst_temperature_sweep()`
   - Need: CLI export of `htst_temperature_sweep()`, GSM interpolation functions

4. **Naming convention**: Python uses `edl_chart()`, WASM uses `compute_edl_chart()`, CLI uses `cmd_edl_chart`
   - Should standardize (Python "clean name" is recommended)

## Next Session Priority

1. **Fill render_bridge exports** (high ROI — 14 functions)
2. **Fill periodic exports** (enables full band structure workflow)
3. **Document all 42 alpha functions** in new `ALPHA_COMPLETE_API.md`
4. **Add end-to-end tests** for each alpha track
5. **Create migration guide** for when code moves to beta/core

## References

- Rust modules: `src/alpha/{edl,periodic_linear,kinetics,render_bridge}/mod.rs`
- Python bindings: `crates/python/src/experimental.rs`
- WASM bindings: `crates/wasm/src/experimental.rs`
- CLI bindings: `crates/cli/src/{experimental_cmds.rs,args.rs,main.rs}`
- Feature gates: `Cargo.toml` (root)
