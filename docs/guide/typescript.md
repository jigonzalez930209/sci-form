# TypeScript / JavaScript API

## Installation

The package is published as **`sci-form-wasm`** on npm.

```bash
npm install sci-form-wasm
```

The package uses WebAssembly for high-performance computation in the browser or Node.js.

## Node.js Usage

### CommonJS
```javascript
const sci = require('sci-form-wasm');

const result = JSON.parse(sci.embed('CCO', 42));
console.log(`Atoms: ${result.num_atoms}`);
```

### ES Modules / TypeScript
```typescript
import init, { embed, embed_batch, parse_smiles, version } from 'sci-form-wasm';

await init();

const result = JSON.parse(embed('CCO', 42));
console.log(`Atoms: ${result.num_atoms}`);
```

## Functions

### `embed(smiles, seed) → string`

Generate a single 3D conformer. Returns a JSON string.

```typescript
import { embed } from 'sci-form-wasm';

const result = JSON.parse(embed("CCO", 42));
console.log(`Atoms: ${result.num_atoms}`);
console.log(`Time: ${result.time_ms.toFixed(1)}ms`);
```

Result JSON structure:
```json
{
  "smiles": "CCO",
  "num_atoms": 9,
  "coords": [0.123, -0.456, 0.789, ...],
  "elements": [6, 6, 8, 1, 1, 1, 1, 1, 1],
  "bonds": [[0, 1, "SINGLE"], [1, 2, "SINGLE"], ...],
  "error": null,
  "time_ms": 15.3
}
```

### `embed_coords(smiles, seed) → string`

Compact coordinate output (no bonds/elements). Returns JSON.

```typescript
import { embed_coords } from 'sci-form-wasm';

const { coords, num_atoms } = JSON.parse(embed_coords("c1ccccc1", 42));
// coords: [x0, y0, z0, x1, y1, z1, ...]
```

### `embed_batch(smiles_list, seed) → string`

Batch-process molecules from a newline-separated string. Returns a JSON array.

```typescript
import { embed_batch } from 'sci-form-wasm';

const smiles = "CCO\nc1ccccc1\nCC(=O)O";
const results = JSON.parse(embed_batch(smiles, 42));

for (const r of results) {
  console.log(`${r.smiles}: ${r.num_atoms} atoms`);
}
```

### `parse_smiles(smiles) → string`

Parse without 3D generation. Returns JSON.

```typescript
import { parse_smiles } from 'sci-form-wasm';

const info = JSON.parse(parse_smiles("c1ccccc1"));
console.log(`Atoms: ${info.num_atoms}, Bonds: ${info.num_bonds}`);
```

### `version() → string`

```typescript
import { version } from 'sci-form-wasm';
console.log(version()); // "sci-form 0.1.0"
```

## Browser Usage

```html
<script type="module">
  import init, { embed } from './sci_form_wasm.js';
  
  await init();
  
  const result = JSON.parse(embed("CCO", 42));
  console.log(`Generated ${result.num_atoms} atoms`);
</script>
```

## React Example

```tsx
import { useEffect, useState } from 'react';
import { embed } from 'sci-form-wasm';

function MoleculeViewer({ smiles }: { smiles: string }) {
  const [result, setResult] = useState<any>(null);

  useEffect(() => {
    const r = JSON.parse(embed(smiles, 42));
    setResult(r);
  }, [smiles]);

  if (!result) return <div>Loading...</div>;
  if (result.error) return <div>Error: {result.error}</div>;

  return (
    <div>
      <p>{result.num_atoms} atoms, {result.time_ms.toFixed(1)}ms</p>
    </div>
  );
}
```

## Three.js Integration

```typescript
import * as THREE from 'three';
import { embed } from 'sci-form-wasm';

const { coords, elements, num_atoms } = JSON.parse(embed("c1ccccc1", 42));

const scene = new THREE.Scene();

for (let i = 0; i < num_atoms; i++) {
  const geometry = new THREE.SphereGeometry(0.3, 16, 16);
  const color = elements[i] === 6 ? 0x333333 :
                elements[i] === 8 ? 0xff0000 :
                elements[i] === 7 ? 0x0000ff : 0xffffff;
  const material = new THREE.MeshPhongMaterial({ color });
  const sphere = new THREE.Mesh(geometry, material);
  sphere.position.set(coords[i * 3], coords[i * 3 + 1], coords[i * 3 + 2]);
  scene.add(sphere);
}
```

## Building from Source

```bash
cd crates/wasm
wasm-pack build --target bundler --release

# For Node.js
wasm-pack build --target nodejs --release

# For browser (no bundler)
wasm-pack build --target web --release
```
