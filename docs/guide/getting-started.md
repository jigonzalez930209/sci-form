# Getting Started

Choose your platform to get started with sci-form.

## Rust

```bash
cargo add sci-form
```

```rust
fn main() {
    let result = sci_form::embed("CCO", 42);
    if let Some(err) = &result.error {
        eprintln!("Error: {}", err);
        return;
    }
    println!("Ethanol: {} atoms", result.num_atoms);
    for i in 0..result.num_atoms {
        let x = result.coords[i * 3];
        let y = result.coords[i * 3 + 1];
        let z = result.coords[i * 3 + 2];
        println!("  Atom {}: ({:.3}, {:.3}, {:.3})", result.elements[i], x, y, z);
    }
}
```

→ See the [Rust guide](/guide/rust) for full API details.

## Python

The package is published as **`sciforma`** on PyPI. The import name is `sci_form`.

```bash
pip install sciforma
```

```python
import sci_form

result = sci_form.embed("c1ccccc1")  # benzene, seed=42 default
if not result.is_ok():
    raise RuntimeError(result.error)

print(f"Atoms: {result.num_atoms}, Time: {result.time_ms:.1f}ms")

for i, (x, y, z) in enumerate(result.get_positions()):
    print(f"  {result.elements[i]:2d}: ({x:7.3f}, {y:7.3f}, {z:7.3f})")
```

→ See the [Python guide](/guide/python) for full API details.

## TypeScript / JavaScript

The npm package is **`sci-form-wasm`**.

```bash
npm install sci-form-wasm
```

**Node.js (CommonJS)**
```javascript
const sci = require('sci-form-wasm');
const result = JSON.parse(sci.embed('CC(=O)O', 42));  // acetic acid
if (result.error) throw new Error(result.error);
console.log(`Atoms: ${result.num_atoms}`);
```

**ES Module / TypeScript**
```typescript
import init, { embed } from 'sci-form-wasm';
await init();
const result = JSON.parse(embed('CC(=O)O', 42));
console.log(`Atoms: ${result.num_atoms}`);
```

→ See the [TypeScript guide](/guide/typescript) for full API details.

## CLI

Download from [GitHub Releases](https://github.com/jigonzalez930209/sci-form/releases) or build from source:

```bash
cargo install sci-form-cli
```

```bash
# Single molecule → XYZ format
sci-form embed "CCO" --format xyz

# Batch processing
sci-form batch -i molecules.smi -o output.sdf --format sdf --threads 8

# Parse SMILES (no 3D generation)
sci-form parse "c1ccccc1"

# Semi-empirical properties from known coordinates
# (here: water, with elements = atomic numbers [O,H,H] and coords = flat [x0,y0,z0,...] Å array)
sci-form pm3 "[8,1,1]" "[0,0,0,0.96,0,0,-0.24,0.93,0]"
```

→ See the [CLI guide](/guide/cli) for all commands and options.
