# TypeScript / JavaScript API Reference

## Installation

```bash
npm install sci-form
```

## Functions

### `embed`

```typescript
function embed(smiles: string, seed?: number): ConformerResult
```

Generate a 3D conformer. Returns a `ConformerResult` object.

```typescript
import init, { embed } from 'sci-form';

await init();
const result = embed('CCO', 42);
console.log(result.num_atoms); // 9
```

### `embed_coords`

```typescript
function embed_coords(smiles: string, seed?: number): Float64Array
```

Generate a conformer and return only the coordinate array. More efficient when you only need positions.

```typescript
const coords = embed_coords('CCO', 42);
// coords = Float64Array [x₀, y₀, z₀, x₁, y₁, z₁, ...]
```

### `embed_batch`

```typescript
function embed_batch(smiles_list: string[], seed?: number): ConformerResult[]
```

Generate conformers for multiple molecules.

```typescript
const results = embed_batch(['CCO', 'c1ccccc1', 'CC(=O)O'], 42);
results.forEach(r => {
  console.log(`${r.smiles}: ${r.num_atoms} atoms`);
});
```

### `parse_smiles`

```typescript
function parse_smiles(smiles: string): MoleculeInfo
```

Parse SMILES string and return molecular graph information.

### `version`

```typescript
function version(): string
```

Returns `"sci-form X.Y.Z"`.

## Types

### `ConformerResult`

```typescript
interface ConformerResult {
  smiles: string;
  num_atoms: number;
  coords: Float64Array;                        // [x₀, y₀, z₀, ...]
  elements: Uint8Array;                        // atomic numbers
  bonds: Array<[number, number, string]>;      // [atom_a, atom_b, order]
  error: string | null;
  time_ms: number;
}
```

### `MoleculeInfo`

```typescript
interface MoleculeInfo {
  num_atoms: number;
  num_bonds: number;
  elements: Uint8Array;
  bonds: Array<[number, number, string]>;
}
```

## Browser Usage

### Vanilla JS

```html
<script type="module">
  import init, { embed } from './node_modules/sci-form/sci_form.js';

  await init();
  const result = embed('c1ccccc1');
  console.log(result.coords);
</script>
```

### React Component

```tsx
import { useEffect, useState } from 'react';
import init, { embed, ConformerResult } from 'sci-form';

export function MoleculeViewer({ smiles }: { smiles: string }) {
  const [result, setResult] = useState<ConformerResult | null>(null);

  useEffect(() => {
    (async () => {
      await init();
      setResult(embed(smiles));
    })();
  }, [smiles]);

  if (!result) return <div>Loading...</div>;

  return (
    <div>
      <p>{result.num_atoms} atoms</p>
      <p>Generated in {result.time_ms.toFixed(1)}ms</p>
    </div>
  );
}
```

### Three.js Integration

```typescript
import init, { embed } from 'sci-form';
import * as THREE from 'three';

await init();
const result = embed('c1ccccc1');

const ELEMENT_COLORS: Record<number, number> = {
  1: 0xffffff,  // H
  6: 0x333333,  // C
  7: 0x3050f8,  // N
  8: 0xff0d0d,  // O
};

const ELEMENT_RADII: Record<number, number> = {
  1: 0.25, 6: 0.4, 7: 0.38, 8: 0.36,
};

const scene = new THREE.Scene();

// Add atoms as spheres
for (let i = 0; i < result.num_atoms; i++) {
  const el = result.elements[i];
  const geometry = new THREE.SphereGeometry(ELEMENT_RADII[el] ?? 0.3);
  const material = new THREE.MeshPhongMaterial({
    color: ELEMENT_COLORS[el] ?? 0x999999,
  });
  const sphere = new THREE.Mesh(geometry, material);
  sphere.position.set(
    result.coords[i * 3],
    result.coords[i * 3 + 1],
    result.coords[i * 3 + 2],
  );
  scene.add(sphere);
}

// Add bonds as cylinders
for (const [a, b] of result.bonds) {
  const start = new THREE.Vector3(
    result.coords[a * 3], result.coords[a * 3 + 1], result.coords[a * 3 + 2],
  );
  const end = new THREE.Vector3(
    result.coords[b * 3], result.coords[b * 3 + 1], result.coords[b * 3 + 2],
  );
  const mid = start.clone().add(end).multiplyScalar(0.5);
  const length = start.distanceTo(end);
  const geometry = new THREE.CylinderGeometry(0.05, 0.05, length);
  const material = new THREE.MeshPhongMaterial({ color: 0xcccccc });
  const cylinder = new THREE.Mesh(geometry, material);
  cylinder.position.copy(mid);
  cylinder.lookAt(end);
  cylinder.rotateX(Math.PI / 2);
  scene.add(cylinder);
}
```
