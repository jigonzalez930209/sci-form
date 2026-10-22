# TypeScript / JavaScript API

## Installation

The package is published as **`sci-form-wasm`** on npm.

```bash
npm install sci-form-wasm
```

The package uses WebAssembly for high-performance computation in the browser or Node.js.

## Conventions

- Every entry point returns a **JSON string**, except the typed-array variants
  (`embed_coords_typed`, `embed_coords_typed_seed64`, `compute_esp_grid_typed`,
  `eht_orbital_grid_typed`, `eht_orbital_grid_from_coefficients_typed`), which return a typed
  array, and `estimate_workers`, which returns a number.
- On failure the payload is always `{"error": "<message>"}` and nothing else — it always parses
  as JSON.
- The exported names match the **Rust root function names**, so the same capability has the same
  name in Rust, Python, WASM and (where exposed) the CLI.
- Seeds are 64-bit. Every seed-taking entry point has a **`*_seed64` variant** that takes a
  JavaScript `BigInt` for the full 64-bit range; the 32-bit variants are kept for compatibility.
- Units: coordinates in **Å**, charges in **e**, dipoles in **Debye**, EHT/PM3/xTB orbital
  energies in **eV**, UFF/MMFF94 energies in **kcal/mol**, HF-3c and ANI energies in **Hartree**,
  IR wavenumbers in **cm⁻¹**, NMR shifts in **ppm** and J-couplings in **Hz**.

## Node.js usage

### CommonJS

```javascript
const sci = require('sci-form-wasm');

const result = JSON.parse(sci.embed('CCO', 42));
console.log(`Atoms: ${result.num_atoms}`);
```

### ES modules / TypeScript

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
  "coords": [-0.8930637240409851, 0.06926389038562775, 0.06360100954771042, "..."],
  "elements": [6, 6, 8, 1, 1, 1, 1, 1, 1],
  "bonds": [[0, 1, "SINGLE"], [1, 2, "SINGLE"], "..."],
  "error": null,
  "time_ms": 0
}
```

The WASM build does not measure wall-clock time, so `time_ms` is always `0`; time the call from
JavaScript if you need a duration.

`embed_seed64(smiles, seed: bigint)` is the 64-bit variant.

### `embed_coords(smiles, seed) → string`

Compact coordinate output (no bonds/elements). Returns JSON
`{"coords": number[], "num_atoms": number}`.

```typescript
import { embed_coords } from 'sci-form-wasm';

const { coords, num_atoms } = JSON.parse(embed_coords("c1ccccc1", 42));
// coords: [x0, y0, z0, x1, y1, z1, ...] in Å
```

`embed_coords_seed64` takes a `BigInt` seed.

### `embed_coords_typed(smiles, seed) → Float64Array`

Coordinates as a typed array — **no JSON parsing overhead**. Preferred for high-throughput use.

```typescript
import { embed_coords_typed } from 'sci-form-wasm';

const coords: Float64Array = embed_coords_typed('CCO', 42);
```

`embed_coords_typed_seed64` takes a `BigInt` seed.

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

`embed_batch_seed64` takes a `BigInt` seed.

### `parse_smiles(smiles) → string`

Parse without 3D generation. Returns JSON
`{"num_atoms": N, "num_bonds": M, "atoms": [{"element": 6, "hybridization": "SP3", "formal_charge": 0}, ...]}`.

```typescript
import { parse_smiles } from 'sci-form-wasm';

const info = JSON.parse(parse_smiles("c1ccccc1"));
console.log(`Atoms: ${info.num_atoms}, Bonds: ${info.num_bonds}`);
```

### `to_canonical_smiles(smiles) → string`

Returns the canonical SMILES as a plain string, e.g. `"BrC(Cl)(F)I"`.

### `version() → string`

```typescript
import { version } from 'sci-form-wasm';
console.log(version()); // "sci-form 1.0.0"
```

## Browser usage

```html
<script type="module">
  import init, { embed } from './sci_form_wasm.js';

  await init();

  const result = JSON.parse(embed("CCO", 42));
  console.log(`Generated ${result.num_atoms} atoms`);
</script>
```

## React example

```tsx
import { useEffect, useState } from 'react';
import init, { embed } from 'sci-form-wasm';

function MoleculeViewer({ smiles }: { smiles: string }) {
  const [result, setResult] = useState<any>(null);

  useEffect(() => {
    (async () => {
      await init();
      setResult(JSON.parse(embed(smiles, 42)));
    })();
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

## Three.js integration

```typescript
import * as THREE from 'three';
import { embed, solve_eht, eht_orbital_mesh } from 'sci-form-wasm';

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

## Molecular properties & analysis

### Extended Hückel Theory (EHT)

```typescript
import { solve_eht, eht_orbital_mesh } from 'sci-form-wasm';

const eht = JSON.parse(solve_eht(elements, coords, 1.75));
console.log(`HOMO energy: ${eht.homo_energy.toFixed(3)} eV`);
console.log(`LUMO energy: ${eht.lumo_energy.toFixed(3)} eV`);
console.log(`Gap: ${eht.gap.toFixed(3)} eV`);

// Generate the 3D orbital isosurface mesh for visualization
const mesh = JSON.parse(eht_orbital_mesh(elements, coords, eht.homo_index, 0.2, 0.02));
// mesh.vertices, mesh.normals, mesh.indices → Three.js BufferGeometry
```

`eht_calculate` is the older alias of `solve_eht` and returns the same payload.

### Electrostatic potential (ESP)

```typescript
import { compute_esp_grid_typed, compute_esp_grid_info } from 'sci-form-wasm';

const spacing = 0.3;   // Å
const padding = 3.0;   // padding around the molecule, Å

const info = JSON.parse(compute_esp_grid_info(elements, coords, spacing, padding));
console.log(`Grid dimensions: ${info.dims}`);   // [nx, ny, nz]
console.log(`Origin: ${info.origin}`);           // [x, y, z]

// Compute grid values (fast: typed array, no JSON overhead)
const espValues = compute_esp_grid_typed(elements, coords, spacing, padding);
// espValues is a Float64Array of length nx * ny * nz, in Hartree/e
```

### WebGPU-accelerated volumetric grids

When the browser build is compiled with `parallel experimental-gpu`, the WASM layer can route
volumetric kernels through WebGPU while keeping CPU fallback and an explicit
`cpu | gpu | hybrid | auto` execution mode.

The accelerated APIs return Promises because WebGPU device creation and dispatch are
asynchronous in the browser.

```typescript
import init, {
  initThreadPool,
  init_webgpu,
  compute_esp_grid_accelerated,
  eht_orbital_grid_accelerated,
  eht_density_grid_accelerated,
} from 'sci-form-wasm';

await init();
await initThreadPool(navigator.hardwareConcurrency ?? 4);

const gpuStatus = JSON.parse(await init_webgpu());
console.log(gpuStatus);

const esp = JSON.parse(
  await compute_esp_grid_accelerated(elements, coords, 0.3, 3.0, 'hybrid')
);
console.log(esp.backend, esp.used_gpu, esp.dims);

const orbital = JSON.parse(
  await eht_orbital_grid_accelerated(elements, coords, 0, 0.25, 3.0, 'gpu')
);
console.log(orbital.backend, orbital.note);

const density = JSON.parse(
  await eht_density_grid_accelerated(elements, coords, 0.25, 3.0, 'auto')
);
console.log(density.backend, density.note);
```

Current WebGPU volumetric coverage in the browser includes:

- `compute_esp_grid_accelerated(...)`
- `eht_orbital_grid_accelerated(...)`
- `eht_density_grid_accelerated(...)`

If you want to avoid serializing large volumetric payloads through JSON, use the typed async
variants and request metadata separately:

```typescript
import {
  compute_esp_grid_accelerated_typed,
  compute_esp_grid_info,
  eht_orbital_grid_accelerated_typed,
  eht_density_grid_accelerated_typed,
  eht_volumetric_grid_info,
} from 'sci-form-wasm';

const espInfo = JSON.parse(compute_esp_grid_info(elements, coords, 0.3, 3.0));
const espValues = await compute_esp_grid_accelerated_typed(
  elements, coords, 0.3, 3.0, 'hybrid'
);

const orbitalInfo = JSON.parse(eht_volumetric_grid_info(elements, coords, 0.25, 3.0));
const orbitalValues = await eht_orbital_grid_accelerated_typed(
  elements, coords, 0, 0.25, 3.0, 'gpu'
);

const densityValues = await eht_density_grid_accelerated_typed(
  elements, coords, 0.25, 3.0, 'auto'
);
```

The metadata calls are cheap because they only compute the bounding box and grid shape; the
large value buffer stays in `Float32Array` or `Float64Array`.

### Accelerated orbital mesh

The browser layer also exposes an accelerated mesh path that reuses the WebGPU orbital-grid
stage before running marching cubes:

```typescript
import { compute_orbital_mesh_accelerated } from 'sci-form-wasm';

const mesh = JSON.parse(
  await compute_orbital_mesh_accelerated(
    elements, coords, 'eht', 0, 0.25, 3.0, 0.02, 'hybrid'
  )
);

console.log(mesh.backend, mesh.used_gpu, mesh.note);
```

This removes the previous browser bottleneck where the full volumetric orbital grid had to be
rebuilt on CPU just to extract the mesh. The current marching-cubes extraction itself still runs
on CPU after the accelerated grid stage.

`hybrid` uses WebGPU for a front slice of the grid and CPU for the remainder. If the build also
enables `parallel` and `initThreadPool()` has been called, the CPU portion can use the Rayon
worker pool.

### Solvent-accessible surface area (SASA)

```typescript
import { compute_sasa } from 'sci-form-wasm';

const result = JSON.parse(compute_sasa(elements, coords, 1.4));
console.log(`Total SASA: ${result.total_sasa.toFixed(2)} Å²`);
console.log(`Per-atom SASA: ${result.atom_sasa}`);   // [a₀, a₁, ...]
```

### Population analysis

```typescript
import { compute_population } from 'sci-form-wasm';

const pop = JSON.parse(compute_population(elements, coords));
console.log(`Mulliken charges: ${pop.mulliken_charges}`);
console.log(`Löwdin charges: ${pop.lowdin_charges}`);
console.log(`HOMO: ${pop.homo_energy.toFixed(3)} eV`);
console.log(`LUMO: ${pop.lumo_energy.toFixed(3)} eV`);
console.log(`Gap: ${pop.gap.toFixed(3)} eV`);
```

### Molecular dipole moment

```typescript
import { compute_dipole } from 'sci-form-wasm';

const dipole = JSON.parse(compute_dipole(elements, coords));
console.log(`|μ| = ${dipole.magnitude.toFixed(3)} Debye`);
console.log(`Vector: [${dipole.vector.join(', ')}] D`);
```

### Density of states (DOS)

```typescript
import { compute_dos } from 'sci-form-wasm';

const dos = JSON.parse(compute_dos(
  elements, coords,
  0.3,           // sigma (Gaussian broadening, eV)
  -30.0,         // e_min (eV)
  5.0,           // e_max (eV)
  500            // n_points
));

console.log(`Energy axis: ${dos.energies}`);   // eV
console.log(`Total DOS: ${dos.total_dos}`);     // 1/eV
console.log(`Smearing: ${dos.sigma} eV`);
```

To build the DOS from another method use
`compute_dos_multimethod(elements, coords, "xtb", 0.3, -30, 5, 500)`.

### Partial charges (Gasteiger-Marsili)

```typescript
import { compute_charges } from 'sci-form-wasm';

const charges = JSON.parse(compute_charges(smiles));
console.log(`Per-atom charges: ${charges.charges}`);
console.log(`Total charge: ${charges.total_charge.toFixed(3)}`);
```

### RMSD alignment

```typescript
import { compute_rmsd } from 'sci-form-wasm';

const reference = JSON.stringify([x0, y0, z0, x1, y1, z1, ...]);
const coordsFlat = JSON.stringify([...]);

const result = JSON.parse(compute_rmsd(coordsFlat, reference));
console.log(`RMSD: ${result.rmsd.toFixed(3)} Å`);
```

### Force field energy (UFF / MMFF94)

```typescript
import { compute_uff_energy, compute_mmff94_energy } from 'sci-form-wasm';

const uff = JSON.parse(compute_uff_energy(smiles, coords));
console.log(`UFF energy: ${uff.energy} ${uff.unit}`);   // kcal/mol

const mmff = JSON.parse(compute_mmff94_energy(smiles, coords));
console.log(`MMFF94 energy: ${mmff.energy} ${mmff.unit}`);
```

## Quantum chemistry

### PM3 semi-empirical SCF

PM3 provides NDDO-level semi-empirical quantum chemistry. PM3(tm) parameters automatically
activate for transition metals (Ti–Au). Energies are in **eV**; `heat_of_formation` is in
**kcal/mol**.

```typescript
import { embed, compute_pm3 } from 'sci-form-wasm';

const conf = JSON.parse(embed('CCO', 42));
const elements = JSON.stringify(conf.elements);
const coords = JSON.stringify(conf.coords);

const pm3 = JSON.parse(compute_pm3(elements, coords));
console.log(`Heat of formation: ${pm3.heat_of_formation.toFixed(2)} kcal/mol`);
console.log(`Total energy: ${pm3.total_energy.toFixed(3)} eV`);
console.log(`HOMO: ${pm3.homo_energy.toFixed(3)} eV`);
console.log(`LUMO: ${pm3.lumo_energy.toFixed(3)} eV`);
console.log(`Gap:  ${pm3.gap.toFixed(3)} eV`);
console.log(`Converged: ${pm3.converged} (${pm3.scf_iterations} iterations)`);
```

### GFN-xTB tight binding

GFN0/GFN1/GFN2-xTB support 25 elements including transition metals (Ti, Cr, Mn, Fe, Co, Ni,
Cu, Zn). All energies are in **eV**.

```typescript
import { compute_xtb, compute_gfn1, compute_gfn2 } from 'sci-form-wasm';

const xtb = JSON.parse(compute_xtb(elements, coords));
console.log(`Total energy: ${xtb.total_energy.toFixed(4)} eV`);
console.log(`Electronic:   ${xtb.electronic_energy.toFixed(4)} eV`);
console.log(`Repulsive:    ${xtb.repulsive_energy.toFixed(4)} eV`);
console.log(`Gap:          ${xtb.gap.toFixed(3)} eV`);
console.log(`SCC converged: ${xtb.converged} (${xtb.scc_iterations} iters)`);
```

### HF-3c

HF-3c energies are in **Hartree**.

```typescript
import { compute_hf3c } from 'sci-form-wasm';

const hf = JSON.parse(compute_hf3c(elements, coords));
console.log(`HF-3c energy: ${hf.energy.toFixed(4)} Hartree`);
```

## Spectroscopy

### UV-Vis (sTDA)

```typescript
import { compute_stda_uvvis } from 'sci-form-wasm';

const uvvis = JSON.parse(compute_stda_uvvis(
  elements, coords,
  0.3,        // sigma (broadening, eV)
  1.0,        // e_min (eV)
  8.0,        // e_max (eV)
  500,        // n_points
  'gaussian'  // broadening type
));

for (const ex of uvvis.excitations) {
  console.log(`  ${ex.energy_ev.toFixed(2)} eV (${ex.wavelength_nm.toFixed(0)} nm), f = ${ex.oscillator_strength.toFixed(3)}`);
}
```

### IR spectroscopy

`compute_vibrational_analysis` runs the numerical Hessian; `compute_ir_spectrum` broadens the
**analysis JSON** it returns.

```typescript
import { compute_vibrational_analysis, compute_ir_spectrum } from 'sci-form-wasm';

// Step 1: Vibrational analysis (numerical Hessian)
const analysis = compute_vibrational_analysis(elements, coords, 'pm3', 0.005);

// Step 2: Broadened IR spectrum from the analysis
const ir = JSON.parse(compute_ir_spectrum(
  analysis,   // JSON string from step 1
  15.0,       // gamma (broadening width, cm⁻¹)
  400.0,      // wn_min (cm⁻¹)
  4000.0,     // wn_max (cm⁻¹)
  1000,       // n_points
));

console.log(`Peaks: ${ir.peaks.length}`);
for (const p of ir.peaks) {
  console.log(`  ${p.frequency_cm1.toFixed(0)} cm⁻¹, intensity: ${p.ir_intensity.toFixed(2)}`);
}
```

`compute_ir_spectrum_broadened` does both steps in one call and lets you choose the
broadening type.

### NMR

```typescript
import { predict_nmr_shifts, predict_nmr_couplings, compute_nmr_spectrum } from 'sci-form-wasm';

// Chemical shifts (ppm)
const shifts = JSON.parse(predict_nmr_shifts('CCO'));
for (const h of shifts.h_shifts) {
  console.log(`H#${h.atom_index}: ${h.shift_ppm.toFixed(2)} ppm [${h.environment}]`);
}

// J-coupling constants (Hz); pass "[]" for the topology-only average
const couplings = JSON.parse(predict_nmr_couplings('CCO', '[]'));
for (const j of couplings) {
  console.log(`  H${j.h1_index}–H${j.h2_index}: ${j.j_hz.toFixed(1)} Hz (${j.coupling_type})`);
}

// Full 1D NMR spectrum
const spectrum = JSON.parse(compute_nmr_spectrum(
  'CCO',
  '1H',     // nucleus: '1H' | '13C' | '35Cl' | '79Br' | '195Pt'
  0.02,     // gamma (broadening, ppm)
  0.0,      // ppm_min
  12.0,     // ppm_max
  1000      // n_points
));
console.log(`Peaks: ${spectrum.peaks.length}`);
```

## Bond orders & reactivity

### Bond order analysis

```typescript
import { compute_bond_orders } from 'sci-form-wasm';

const bo = JSON.parse(compute_bond_orders(elements, coords));
for (const b of bo.bonds) {
  console.log(`Bond ${b.atom_i}–${b.atom_j}: Wiberg = ${b.wiberg.toFixed(2)}, Mayer = ${b.mayer.toFixed(2)} (${b.distance.toFixed(3)} Å)`);
}
```

### Frontier molecular orbital descriptors

```typescript
import { compute_frontier_descriptors } from 'sci-form-wasm';

const fmo = JSON.parse(compute_frontier_descriptors(elements, coords));
console.log(`HOMO: ${fmo.homo_energy.toFixed(3)} eV`);
console.log(`LUMO: ${fmo.lumo_energy.toFixed(3)} eV`);
console.log(`Gap:  ${fmo.gap.toFixed(3)} eV`);
console.log(`HOMO contributions: ${fmo.homo_atom_contributions}`);
console.log(`LUMO contributions: ${fmo.lumo_atom_contributions}`);
```

### Fukui reactivity indices

```typescript
import { compute_fukui_descriptors } from 'sci-form-wasm';

const fukui = JSON.parse(compute_fukui_descriptors(elements, coords));
console.log(`f⁺ (nucleophilic attack): ${fukui.f_plus}`);
console.log(`f⁻ (electrophilic attack): ${fukui.f_minus}`);
console.log(`f⁰ (radical attack): ${fukui.f_radical}`);
```

### Reactivity site ranking

```typescript
import { compute_reactivity_ranking } from 'sci-form-wasm';

const ranking = JSON.parse(compute_reactivity_ranking(elements, coords));
console.log('Nucleophilic attack sites:', ranking.nucleophilic_attack_sites);
console.log('Electrophilic attack sites:', ranking.electrophilic_attack_sites);
console.log('Radical attack sites:', ranking.radical_attack_sites);
```

### Empirical pKa prediction

```typescript
import { compute_empirical_pka } from 'sci-form-wasm';

const pka = JSON.parse(compute_empirical_pka('CC(=O)O'));   // acetic acid
console.log('Acidic sites:');
for (const site of pka.acidic_sites) {
  console.log(`  Atom ${site.atom_index}: pKa ≈ ${site.estimated_pka.toFixed(1)} (${site.site_type})`);
}
```

## ML properties & descriptors

ML property prediction works from SMILES alone — no 3D coordinates needed.

```typescript
import { predict_ml_properties, compute_ml_descriptors } from 'sci-form-wasm';

// Predict druglikeness, LogP, solubility, Lipinski RO5
const props = JSON.parse(predict_ml_properties('CC(=O)Nc1ccc(O)cc1'));   // acetaminophen
console.log(`LogP: ${props.logp.toFixed(2)}`);
console.log(`Molar refractivity: ${props.molar_refractivity.toFixed(2)}`);
console.log(`Log solubility: ${props.log_solubility.toFixed(2)}`);
console.log(`Druglikeness: ${props.druglikeness.toFixed(3)}`);
console.log(`Lipinski passes: ${props.lipinski_passes}`);

// Compute molecular descriptors for custom ML models
const desc = JSON.parse(compute_ml_descriptors('CC(=O)Nc1ccc(O)cc1'));
console.log(`MW: ${desc.molecular_weight.toFixed(1)}`);
console.log(`Rotatable bonds: ${desc.n_rotatable_bonds}`);
console.log(`HBD: ${desc.n_hbd}, HBA: ${desc.n_hba}`);
console.log(`FSP3: ${desc.fsp3.toFixed(2)}`);
console.log(`Wiener index: ${desc.wiener_index}`);
console.log(`Balaban J: ${desc.balaban_j.toFixed(3)}`);
```

`compute_ml_properties` and `compute_molecular_descriptors` are the deprecated aliases of
`predict_ml_properties` and `compute_ml_descriptors`.

## Stereochemistry, fingerprints & solvation

### Stereochemistry analysis

```typescript
import { analyze_stereo } from 'sci-form-wasm';

const stereo = JSON.parse(analyze_stereo('C(F)(Cl)(Br)I', '[]'));
console.log(`Stereocenters: ${stereo.n_stereocenters}`);
for (const sc of stereo.stereocenters) {
  console.log(`  Atom ${sc.atom_index}: ${sc.configuration}`);   // R or S
}
console.log(`Double bonds with E/Z: ${stereo.n_double_bonds}`);
```

### ECFP fingerprints & Tanimoto similarity

```typescript
import { compute_ecfp, compute_tanimoto } from 'sci-form-wasm';

const fp1 = compute_ecfp('c1ccccc1', 2, 2048);
const fp2 = compute_ecfp('Cc1ccccc1', 2, 2048);
const similarity = JSON.parse(compute_tanimoto(fp1, fp2));
console.log(`Tanimoto: ${similarity.tanimoto.toFixed(3)}`);
```

### Ring perception (SSSR)

```typescript
import { compute_sssr } from 'sci-form-wasm';

const sssr = JSON.parse(compute_sssr('c1ccc2ccccc2c1'));   // naphthalene
console.log(`Rings: ${sssr.rings.length}`);
console.log(`Ring sizes: ${JSON.stringify(sssr.ring_size_histogram)}`);
```

### Solvation

```typescript
import { compute_nonpolar_solvation, compute_gb_solvation } from 'sci-form-wasm';

// Non-polar (SASA-based)
const np = JSON.parse(compute_nonpolar_solvation(elements, coords, 1.4));
console.log(`Non-polar solvation: ${np.energy_kcal_mol.toFixed(2)} kcal/mol`);

// Generalized Born electrostatic solvation
const charges = JSON.stringify([-0.8, 0.4, 0.4]);
const gb = JSON.parse(compute_gb_solvation(elements, coords, charges, 78.5, 1.0, 1.4));
if (gb.error) throw new Error(gb.error);   // invalid input → {"error": "..."}
console.log(`GB electrostatic: ${gb.electrostatic_energy_kcal_mol.toFixed(2)} kcal/mol`);
console.log(`Total solvation: ${gb.total_energy_kcal_mol.toFixed(2)} kcal/mol`);
```

### Clustering (Butina RMSD)

```typescript
import { butina_cluster } from 'sci-form-wasm';

const conformers = JSON.stringify([coords1, coords2, coords3, ...]);
const clusters = JSON.parse(butina_cluster(conformers, 1.0));
console.log(`Clusters: ${clusters.n_clusters}`);
console.log(`Assignments: ${clusters.assignments}`);
```

### Materials & crystallography

```typescript
import { create_unit_cell, assemble_framework } from 'sci-form-wasm';

const cell = JSON.parse(create_unit_cell(5.43, 5.43, 5.43, 90, 90, 90));
console.log(`Volume: ${cell.volume.toFixed(2)} Å³`);

const framework = JSON.parse(assemble_framework('pcu', 30, 'octahedral', 10.0, 1));
console.log(`Framework atoms: ${framework.atoms.length}`);
```

## Parallelism and WebGPU in the browser

Browser-side CPU parallelism uses the **rayon work-stealing thread pool** when the WASM package
is built with the `parallel` feature and the page is running with `SharedArrayBuffer` available.
This provides:

- **Automatic thread scheduling** — no manual thread pool configuration needed
- **Shared work queue** — efficient load balancing across cores
- **Zero-copy data sharing** — immutable borrows avoid synchronization overhead
- **Intra-library parallelism** — grid point evaluation, population loops, force field terms all parallelized

At runtime you must initialize the worker pool explicitly:

```typescript
import init, { initThreadPool } from 'sci-form-wasm';

await init();
await initThreadPool(navigator.hardwareConcurrency ?? 4);
```

After that, CPU-heavy functions automatically use the worker pool when the library path has a
parallel implementation:

```typescript
// Browser: these run in parallel across cores
const sasa = compute_sasa(elements, coords, 1.4);               // Shrake-Rupley per atom
const esp = compute_esp_grid_typed(elements, coords, 0.5, 3.0); // grid point evaluation
const dos = compute_dos(elements, coords, 0.3, -30, 5, 500);    // energy grid
const pop = compute_population(elements, coords);              // per-atom charge partitioning
```

For browser-side hybrid CPU + GPU execution, the recommended sequence is:

```typescript
import init, { initThreadPool, init_webgpu } from 'sci-form-wasm';

await init();
await initThreadPool(navigator.hardwareConcurrency ?? 4);
await init_webgpu();
```

Required response headers for browser threads:

- `Cross-Origin-Opener-Policy: same-origin`
- `Cross-Origin-Embedder-Policy: require-corp`

Without those headers, `SharedArrayBuffer` is unavailable and the browser cannot start the
Rayon worker pool even if the package was compiled with `parallel`.

**Node.js limitation:** Web Workers are not available in standard Node.js workers, so the
Node.js build is compiled without the `parallel` feature. Functions execute sequentially but
maintain full API compatibility. For CPU-intensive computations in Node.js, consider:

- Using the Python bindings (which support true process-level parallelization)
- Batching requests across multiple Node processes
- Running the web version in a headless browser (Puppeteer, etc.)

**Performance note:** for molecules with **< 20 atoms** or **< 1000 grid points**, the overhead of
thread-queue scheduling may exceed the speedup. Sequential execution is automatically selected in
these cases.

## Building from source

```bash
cd crates/wasm

# Browser build for threaded CPU + WebGPU
./build.sh --web-only --web-features "parallel experimental-gpu"

# Browser build for threaded CPU + WebGPU + extra experimental kernels
./build.sh --web-only --web-features "parallel experimental-gpu experimental-eeq experimental-d4 experimental-cpm"

# Node.js build - sequential (no parallelization)
./build.sh --node-only
```

**Targets:**

- `./build.sh --web-only --web-features "parallel experimental-gpu"` — canonical browser path for threaded CPU + WebGPU
- `--target nodejs` (via `./build.sh --node-only`) — server-side Node.js without Web Workers (sequential execution)
- `--target web` remains the correct target for bundlers like Vite; import the generated `pkg/` output from the bundler app
- `--target bundler` is not the correct target for browser threads here

### Vite reference setup

There is a minimal working example in
[crates/wasm/examples/vite-webgpu-hybrid](https://github.com/jigonzalez930209/sci-form/tree/main/crates/wasm/examples/vite-webgpu-hybrid).

Its `vite.config.mjs` sets both development and preview headers:

```javascript
server: {
  headers: {
    'Cross-Origin-Opener-Policy': 'same-origin',
    'Cross-Origin-Embedder-Policy': 'require-corp',
  },
}
```

That example imports the `pkg/` artifacts generated by `./build.sh` and initializes both
`initThreadPool()` and `init_webgpu()` before calling the accelerated volumetric APIs.

## Full export index

Every export of the default `sci-form-wasm` package, grouped by capability. The full reference,
with payload shapes and units, is in the [TypeScript / JavaScript API reference](/api/typescript).

| Area | Exports |
|---|---|
| Metadata | `version`, `alpha_modules_info`, `beta_modules_info` |
| Conformer generation | `embed`, `embed_seed64`, `embed_coords`, `embed_coords_seed64`, `embed_coords_typed`, `embed_coords_typed_seed64`, `embed_batch`, `embed_batch_seed64`, `embed_diverse`, `parse_smiles`, `to_canonical_smiles` |
| Charges, surface, electrostatics | `compute_charges`, `compute_charges_configured`, `compute_sasa`, `compute_population`, `compute_dipole`, `compute_esp`, `compute_esp_grid_info`, `compute_esp_grid_typed` |
| DOS and alignment | `compute_dos`, `compute_dos_multimethod`, `compute_rmsd`, `compute_rmsd_matrix` |
| Extended Hückel theory | `eht_calculate`, `solve_eht`, `eht_support`, `get_eht_support`, `compute_eht_or_uff_fallback`, `eht_or_uff_fallback`, `compute_eht_gradient`, `optimize_geometry_eht` |
| Electronic structure | `compute_pm3`, `compute_xtb`, `compute_gfn1`, `compute_gfn2`, `compute_hf3c`, `compute_hf3c_custom`, `compute_ani`, `compute_cisd`, `compute_aevs_tm` |
| Force fields | `compute_uff_energy`, `compute_uff_energy_with_aromatic_heuristics`, `compute_mmff94_energy`, `compute_empirical_pka` |
| Spectroscopy | `compute_stda_uvvis`, `compute_fluorescence`, `compute_vibrational_analysis`, `compute_vibrational_analysis_uff`, `compute_ir_spectrum`, `compute_ir_spectrum_broadened`, `assign_peaks` |
| NMR | `predict_nmr_shifts`, `predict_nmr_shifts_for_nucleus`, `predict_nmr_couplings`, `compute_ensemble_j_couplings`, `compute_nmr_spectrum`, `compute_nmr_spectrum_with_coords`, `compute_giao_nmr`, `compute_hose_codes` |
| Solvation | `compute_nonpolar_solvation`, `compute_gb_solvation` |
| Rings, fingerprints, clustering | `compute_sssr`, `compute_ecfp`, `compute_ecfp_batch`, `compute_tanimoto`, `butina_cluster`, `filter_diverse` |
| Descriptors and machine learning | `compute_ml_descriptors`, `compute_molecular_descriptors`, `predict_ml_properties`, `compute_ml_properties`, `compute_whim`, `compute_rdf`, `compute_getaway`, `train_random_forest`, `train_gradient_boosting` |
| Reactivity and topology | `compute_bond_orders`, `compute_frontier_descriptors`, `compute_fukui_descriptors`, `compute_reactivity_ranking`, `compute_uv_vis_spectrum`, `compute_topology`, `analyze_graph_features`, `analyze_stereo` |
| Orbital meshes and grids | `compute_orbital_mesh`, `eht_orbital_mesh`, `eht_orbital_grid_typed`, `eht_orbital_grid_from_coefficients_typed`, `eht_volumetric_grid_info` |
| Materials and periodic systems | `create_unit_cell`, `assemble_framework`, `space_group_by_number`, `space_group_by_symbol`, `all_space_groups`, `build_periodic_molecule`, `detect_hapticity`, `optimize_framework` |
| Dynamics, NEB, reaction paths | `compute_md_trajectory`, `compute_md_trajectory_seed64`, `compute_md_trajectory_nvt`, `compute_md_trajectory_nvt_seed64`, `compute_md_trajectory_with_method`, `compute_md_trajectory_with_method_seed64`, `compute_simplified_neb_path`, `compute_neb_path_with_method`, `neb_single_point_energy`, `neb_energy_and_gradient`, `list_neb_backends`, `search_conformers_with_uff`, `search_conformers_with_uff_seed64`, `compute_reaction_dynamics` |
| SMIRKS | `parse_smirks`, `apply_smirks` |
| Transport and worker helpers | `pack_conformers`, `pack_batch_arrow`, `split_worker_tasks`, `split_worker_tasks_seed64`, `estimate_workers` |
| System planning | `system_capabilities`, `get_system_capabilities`, `system_method_plan`, `get_system_method_plan`, `compare_methods` |
| Experimental (always built) | `compute_eeq_charges`, `compute_eeq_energy` |

Total: **128** exports in the default build. Feature-gated exports (`compute_alpb_solvation`, `compute_d4_energy`, `compute_cpm_charges`, the `alpha-periodic-linear`, `alpha-kinetics` and `alpha-render-bridge` families, and the `alpha_*` / `beta_*` research functions) appear only when the matching Cargo feature is enabled at build time.

---

## See also

- [TypeScript / JavaScript API reference](/api/typescript)
- [Python reference](/api/python)
- [Rust reference](/api/rust)
