# Rust API

## Installation

Add sci-form to your `Cargo.toml`:

```toml
[dependencies]
sci-form = "0.1"
```

For parallel batch processing, enable the `parallel` feature:

```toml
[dependencies]
sci-form = { version = "0.1", features = ["parallel"] }
```

## Core Types

### `ConformerResult`

The result of a 3D conformer generation.

```rust
pub struct ConformerResult {
    /// Input SMILES string
    pub smiles: String,
    /// Number of atoms
    pub num_atoms: usize,
    /// Flat coordinates: [x₀, y₀, z₀, x₁, y₁, z₁, ...]
    pub coords: Vec<f64>,
    /// Atomic numbers in the same order as coordinates
    pub elements: Vec<u8>,
    /// Bonds as (atom_a, atom_b, order_string)
    pub bonds: Vec<(usize, usize, String)>,
    /// Error message (None if successful)
    pub error: Option<String>,
    /// Generation time in milliseconds
    pub time_ms: f64,
}
```

Bond order strings: `"SINGLE"`, `"DOUBLE"`, `"TRIPLE"`, `"AROMATIC"`.

### `ConformerConfig`

Configuration for batch generation.

```rust
pub struct ConformerConfig {
    /// RNG seed for reproducibility (default: 42)
    pub seed: u64,
    /// Number of threads, 0 = auto-detect (default: 0)
    pub num_threads: usize,
}
```

## Functions

### `embed(smiles, seed) → ConformerResult`

Generate a single 3D conformer.

```rust
let result = sci_form::embed("CCO", 42);
assert!(result.error.is_none());
assert_eq!(result.num_atoms, 9); // 2C + 1O + 6H
```

The same seed always produces the same coordinates (deterministic).

### `embed_batch(smiles_list, config) → Vec<ConformerResult>`

Generate conformers for multiple molecules, optionally in parallel.

```rust
let smiles = vec!["CCO", "c1ccccc1", "CC(=O)O"];
let config = sci_form::ConformerConfig { seed: 42, num_threads: 4 };
let results = sci_form::embed_batch(&smiles, &config);

for r in &results {
    println!("{}: {} atoms, {:.1}ms", r.smiles, r.num_atoms, r.time_ms);
}
```

Requires the `parallel` feature for multi-threaded execution. Without it, molecules are processed sequentially.

### `parse(smiles) → Result<Molecule, String>`

Parse a SMILES string into a molecular graph without generating 3D coordinates.

```rust
let mol = sci_form::parse("c1ccccc1").unwrap();
println!("Atoms: {}, Bonds: {}", mol.graph.node_count(), mol.graph.edge_count());
```

### `version() → String`

```rust
println!("{}", sci_form::version()); // "sci-form 0.1.0"
```

## Working with Coordinates

```rust
let result = sci_form::embed("CC(=O)O", 42);
if result.error.is_none() {
    // Access individual atom coordinates
    for i in 0..result.num_atoms {
        let x = result.coords[i * 3];
        let y = result.coords[i * 3 + 1];
        let z = result.coords[i * 3 + 2];
        let elem = result.elements[i];
        println!("Atom {} (Z={}): ({:.3}, {:.3}, {:.3})", i, elem, x, y, z);
    }
}
```

## Serialization

`ConformerResult` implements `Serialize` and `Deserialize` (serde), so you can convert to/from JSON:

```rust
let result = sci_form::embed("CCO", 42);
let json = serde_json::to_string_pretty(&result).unwrap();
println!("{}", json);
```
