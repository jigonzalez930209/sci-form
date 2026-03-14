# Rust API Reference

## Crate: `sci-form`

```toml
[dependencies]
sci-form = "0.1"
# For parallel batch processing:
sci-form = { version = "0.1", features = ["parallel"] }
```

## Types

### `ConformerResult`

```rust
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConformerResult {
    pub smiles: String,
    pub num_atoms: usize,
    pub coords: Vec<f64>,                    // [x₀, y₀, z₀, x₁, y₁, z₁, ...]
    pub elements: Vec<u8>,                   // atomic numbers
    pub bonds: Vec<(usize, usize, String)>,  // (atom_a, atom_b, order)
    pub error: Option<String>,
    pub time_ms: f64,
}
```

| Field | Description |
|-------|-------------|
| `smiles` | The input SMILES string |
| `num_atoms` | Total number of atoms (including implicit H) |
| `coords` | Flat array of 3D coordinates, length = `num_atoms × 3` |
| `elements` | Atomic numbers: 1=H, 6=C, 7=N, 8=O, 9=F, ... |
| `bonds` | Tuples of (start, end, order). Order ∈ {"SINGLE", "DOUBLE", "TRIPLE", "AROMATIC"} |
| `error` | `None` on success, `Some(msg)` on failure |
| `time_ms` | Generation time in milliseconds |

### `ConformerConfig`

```rust
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConformerConfig {
    pub seed: u64,        // default: 42
    pub num_threads: usize, // default: 0 (auto-detect)
}
```

## Functions

### `embed`

```rust
pub fn embed(smiles: &str, seed: u64) -> ConformerResult
```

Generate a 3D conformer for a single molecule. Always returns a `ConformerResult` — check `error` field for failures.

**Parameters:**
- `smiles` — Valid SMILES string
- `seed` — RNG seed. Same seed = same output (deterministic)

### `embed_batch`

```rust
pub fn embed_batch(smiles_list: &[&str], config: &ConformerConfig) -> Vec<ConformerResult>
```

Generate conformers for multiple molecules. When compiled with the `parallel` feature, uses rayon for multi-threaded processing.

### `parse`

```rust
pub fn parse(smiles: &str) -> Result<graph::Molecule, String>
```

Parse SMILES into molecular graph. Returns `Ok(Molecule)` with atoms, bonds, ring info, or `Err(msg)` on parse failure.

### `version`

```rust
pub fn version() -> String
```

Returns `"sci-form X.Y.Z"`.

## Internal Modules

These modules are public and can be used for advanced customization:

| Module | Description |
|--------|-------------|
| `sci_form::graph` | Molecular graph types (`Molecule`, `Atom`, `Bond`, `Hybridization`) |
| `sci_form::smiles` | SMILES parser |
| `sci_form::smarts` | SMARTS parser and substructure matcher |
| `sci_form::conformer` | ETKDG pipeline |
| `sci_form::distgeom` | Distance geometry (bounds matrix, embedding) |
| `sci_form::forcefield` | Force field implementations |
| `sci_form::etkdg` | Torsion angle utilities |
| `sci_form::optimization` | BFGS optimizer |
