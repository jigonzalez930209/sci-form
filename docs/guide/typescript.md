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

## Molecular Properties & Analysis

### Extended Hückel Theory (EHT)

Compute electronic structure, orbital energies, and visualization:

```typescript
import { eht_calculate, eht_orbital_mesh } from 'sci-form-wasm';

const eht = JSON.parse(eht_calculate(elements, coords, 0.0));
console.log(`HOMO energy: ${eht.homo_energy.toFixed(3)} eV`);
console.log(`LUMO energy: ${eht.lumo_energy.toFixed(3)} eV`);
console.log(`Gap: ${eht.gap.toFixed(3)} eV`);

// Generate 3D orbital isosurface mesh for visualization
const mesh = JSON.parse(eht_orbital_mesh(elements, coords, eht.homo_index, 0.2, 0.02));
// mesh.vertices, mesh.normals, mesh.indices → Three.js BufferGeometry
```

### Electrostatic Potential (ESP)

```typescript
import { compute_esp_grid_typed, compute_esp_grid_info } from 'sci-form-wasm';

const spacing = 0.3;  // Ångströms
const padding = 3.0;  // padding around molecule

const info = JSON.parse(compute_esp_grid_info(elements, coords, spacing, padding));
console.log(`Grid dimensions: ${info.dims}`);  // [nx, ny, nz]
console.log(`Origin: ${info.origin}`);        // [x, y, z]

// Compute grid values (Fast: typed array, no JSON overhead)
const espValues = compute_esp_grid_typed(elements, coords, spacing, padding);
// espValues is Float64Array of length nx * ny * nz
```

### Solvent-Accessible Surface Area (SASA)

```typescript
import { compute_sasa } from 'sci-form-wasm';

const result = JSON.parse(compute_sasa(elements, coords, 1.4));
console.log(`Total SASA: ${result.total_sasa.toFixed(2)} Ų`);
console.log(`Per-atom SASA: ${result.per_atom_sasa}`);  // [a₀, a₁, ...]
```

### Population Analysis

```typescript
import { compute_population } from 'sci-form-wasm';

const pop = JSON.parse(compute_population(elements, coords));
console.log(`Mulliken charges: ${pop.mulliken_charges}`);
console.log(`Löwdin charges: ${pop.lowdin_charges}`);
console.log(`HOMO: ${pop.homo_energy.toFixed(3)} eV`);
console.log(`LUMO: ${pop.lumo_energy.toFixed(3)} eV`);
```

### Molecular Dipole Moment

```typescript
import { compute_dipole } from 'sci-form-wasm';

const dipole = JSON.parse(compute_dipole(elements, coords));
console.log(`|μ| = ${dipole.magnitude.toFixed(3)} Debye`);
console.log(`Vector: [${dipole.vector.join(', ')}] D`);
```

### Density of States (DOS)

```typescript
import { compute_dos } from 'sci-form-wasm';

const dos = JSON.parse(compute_dos(
  elements, coords,
  0.3,           // sigma (Gaussian broadening)
  -30.0,         // e_min (eV)
  5.0,           // e_max (eV)
  500            // n_points
));

console.log(`HOMO–LUMO gap: ${dos.homo_lumo_gap.toFixed(3)} eV`);
console.log(`Orbital energies: ${dos.orbital_energies}`);
console.log(`Total DOS: ${dos.total_dos}`);     // [dos₀, dos₁, ...]
console.log(`Projected DOS: ${dos.pdos}`);     // per-orbital DOS
```

### Partial Charges (Gasteiger-Marsili)

```typescript
import { compute_charges } from 'sci-form-wasm';

const charges = JSON.parse(compute_charges(smiles));
console.log(`Per-atom charges: ${charges.charges}`);
console.log(`Total charge: ${charges.total_charge.toFixed(3)}`);
```

### RMSD Alignment

```typescript
import { compute_rmsd } from 'sci-form-wasm';

const reference = JSON.stringify([x0, y0, z0, x1, y1, z1, ...]);
const coords_flat = JSON.stringify([...]);

const result = JSON.parse(compute_rmsd(coords_flat, reference));
console.log(`RMSD: ${result.rmsd.toFixed(3)} Ų`);
console.log(`Aligned coords: ${result.aligned_coords}`);
```

### Force Field Energy (UFF)

```typescript
import { compute_uff_energy } from 'sci-form-wasm';

const energy = JSON.parse(compute_uff_energy(smiles, coords));
console.log(`UFF energy: ${energy} kcal/mol`);
```

## Parallel Computation (Browser Only)

All compute functions in the **browser** use **`rayon` work-stealing thread pool** when built with the `parallel` feature (default). This provides:

- **Automatic thread scheduling** — no manual thread pool configuration needed
- **Shared work queue** — efficient load balancing across cores
- **Zero-copy data sharing** — immutable borrows avoid synchronization overhead
- **Intra-library parallelism** — grid point evaluation, population loops, force field terms all parallelized

To use: No changes required. Functions automatically use available CPU cores:

```typescript
// Browser: These run in parallel across cores
const sasa = compute_sasa(elements, coords, 1.4);           // Shrake-Rupley per-atom parallel
const esp = compute_esp_grid_typed(elements, coords, 0.3);  // Grid point evaluation parallel
const dos = compute_dos(elements, coords, 0.3, -30, 5, 500); // Energy grid parallel
const pop = compute_population(elements, coords);            // All 14 functions parallelizable
```

**Node.js Limitation:** Web Workers are not available in standard Node.js workers, so the Node.js build is compiled without the `parallel` feature. Functions execute sequentially but maintain full API compatibility. For CPU-intensive computations in Node.js, consider:
- Using Python bindings (which support true process-level parallelization)
- Batching requests across multiple Node processes
- Running the web version in a headless browser (Puppeteer, etc.)

**Performance Note:** For molecules with **< 20 atoms** or **< 1000 grid points**, overhead of thread queue scheduling may exceed speedup. Sequential execution is automatically selected in these cases.

## Building from Source

```bash
cd crates/wasm

# Browser build (Vite, Webpack, Parcel, etc.) - WITH parallelization
wasm-pack build --target web --release

# Node.js build - sequential (no parallelization)
wasm-pack build --target nodejs --release --out-dir pkg-node --no-default-features
```

**Targets:**
- `--target web` — Browser with Web Workers support for parallelization (requires `parallel` feature)
- `--target nodejs` — Server-side Node.js without Web Workers (compiled without `parallel` feature)
- ⚠️ `--target bundler` is NOT compatible with threading; use `web` instead
