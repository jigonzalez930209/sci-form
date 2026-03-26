# TypeScript / JavaScript API Reference

## Installation

```bash
npm install sci-form-wasm
```

All functions accept and return JSON strings (or typed arrays), since WASM cannot pass complex objects directly.

---

## Initialization

All platforms require initializing the WASM module once:

```typescript
import init, { embed, eht_calculate, /* ... */ } from 'sci-form-wasm';
await init();
```

For Node.js CommonJS:

```javascript
const sci = require('sci-form-wasm');
// No async init required for Node.js pkg
```

---

## Conformer Generation

### `embed`

```typescript
function embed(smiles: string, seed?: number): string
```

Returns a JSON string with the full conformer result.

```typescript
const json = embed('CCO', 42);
const result = JSON.parse(json);
// result.num_atoms, result.coords, result.elements, result.bonds, result.error
```

### `embed_coords`

```typescript
function embed_coords(smiles: string, seed?: number): string
```

Returns compact JSON: `{"coords": [x0,y0,z0,...], "num_atoms": N}`.

### `embed_coords_typed`

```typescript
function embed_coords_typed(smiles: string, seed: number): Float64Array
```

Returns coordinates as a typed array directly — **no JSON parsing overhead**. Preferred for high-throughput use.

```typescript
const coords: Float64Array = embed_coords_typed('CCO', 42);
// coords = Float64Array [x₀, y₀, z₀, x₁, y₁, z₁, ...]
```

### `embed_batch`

```typescript
function embed_batch(smiles_list: string, seed?: number): string
```

Batch embed from newline-separated SMILES. Returns JSON array.

```typescript
const smiles = "CCO\nc1ccccc1\nCC(=O)O";
const results = JSON.parse(embed_batch(smiles, 42));
results.forEach(r => console.log(`${r.smiles}: ${r.num_atoms} atoms`));
```

### `parse_smiles`

```typescript
function parse_smiles(smiles: string): string
```

Returns `{"num_atoms": N, "num_bonds": M}` or `{"error": "..."}`.

---

## EHT — Extended Hückel Theory

### `eht_calculate`

```typescript
function eht_calculate(
  elements: string,    // JSON: [8, 6, 6, 1, 1, 1, 1, 1, 1]
  coords_flat: string, // JSON: [x0,y0,z0,x1,y1,z1,...]
  k: number            // Wolfsberg-Helmholtz constant (0.0 → default 1.75)
): string
```

Returns JSON:
```typescript
interface EhtResult {
  energies: number[];     // eV, all MO energies
  n_electrons: number;
  homo_index: number;
  lumo_index: number;
  homo_energy: number;    // eV
  lumo_energy: number;    // eV
  gap: number;            // HOMO-LUMO gap eV
}
```

```typescript
const coords = embed_coords_typed('CCO', 42);
const elements = JSON.stringify([8, 6, 6, 1, 1, 1, 1, 1, 1]);
const coordsJson = JSON.stringify(Array.from(coords));

const eht = JSON.parse(eht_calculate(elements, coordsJson, 0.0));
console.log(`HOMO: ${eht.homo_energy.toFixed(3)} eV`);
console.log(`Gap: ${eht.gap.toFixed(3)} eV`);
```

### `eht_orbital_mesh`

```typescript
function eht_orbital_mesh(
  elements: string,
  coords_flat: string,
  mo_index: number,
  spacing: number,    // Å, e.g. 0.2
  isovalue: number    // e.g. 0.02
): string
```

Returns JSON for rendering the molecular orbital as a 3D mesh:

```typescript
interface OrbitalMesh {
  vertices: number[];     // flat [x,y,z, ...]
  normals: number[];      // flat [nx,ny,nz, ...]
  indices: number[];      // flat triangle indices [i0,i1,i2, ...]
  num_triangles: number;
  isovalue: number;
}
```

```typescript
const mesh = JSON.parse(eht_orbital_mesh(elements, coordsJson, homoIdx, 0.2, 0.02));

// Three.js example:
const geometry = new THREE.BufferGeometry();
geometry.setAttribute('position',
  new THREE.Float32BufferAttribute(mesh.vertices, 3));
geometry.setAttribute('normal',
  new THREE.Float32BufferAttribute(mesh.normals, 3));
geometry.setIndex(mesh.indices);
```

### `eht_orbital_grid_typed`

```typescript
function eht_orbital_grid_typed(
  elements: string,
  coords_flat: string,
  mo_index: number,
  spacing: number
): Float32Array
```

Returns the raw volumetric grid as a `Float32Array` for GPU-based volume rendering (WebGL, WebGPU). Use `compute_esp_grid_info` to get grid dimensions.

---

## Electrostatic Potential

### `compute_esp_grid_typed`

```typescript
function compute_esp_grid_typed(
  elements: string,   // JSON array of atomic numbers
  coords_flat: string, // JSON flat xyz array
  spacing: number,    // grid interval in Å
  padding: number     // extra space around molecule in Å
): Float64Array
```

Returns ESP grid values as `Float64Array` (no JSON). Pair with `compute_esp_grid_info` for grid metadata.

```typescript
const espData = compute_esp_grid_typed(elements, coordsJson, 0.5, 3.0);
const espInfo = JSON.parse(compute_esp_grid_info(elements, coordsJson, 0.5, 3.0));
// espInfo = {origin:[x,y,z], spacing:0.5, dims:[nx,ny,nz]}
```

### `compute_esp_grid_info`

```typescript
function compute_esp_grid_info(
  elements: string,
  coords_flat: string,
  spacing: number,
  padding: number
): string
```

Returns JSON:
```typescript
interface EspInfo {
  origin: [number, number, number];
  spacing: number;
  dims: [number, number, number];  // [nx, ny, nz]
}
```

---

## Gasteiger Charges

### `compute_charges`

```typescript
function compute_charges(smiles: string): string
```

Returns JSON:
```typescript
interface ChargeResult {
  charges: number[];      // per-atom partial charges
  iterations: number;
  total_charge: number;
}
```

```typescript
const result = JSON.parse(compute_charges('CCO'));
console.log(result.charges);       // [-0.387, -0.042, -0.228, ...]
```

---

## SASA

### `compute_sasa`

```typescript
function compute_sasa(
  elements: string,
  coords_flat: string,
  probe_radius: number   // 0.0 → default 1.4 Å
): string
```

Returns JSON:
```typescript
interface SasaResult {
  total_sasa: number;    // Å²
  atom_sasa: number[];   // per-atom Å²
  probe_radius: number;
  num_points: number;
}
```

---

## Population Analysis

### `compute_population`

```typescript
function compute_population(
  elements: string,
  coords_flat: string
): string
```

Returns JSON:
```typescript
interface PopulationResult {
  mulliken_charges: number[];
  lowdin_charges: number[];
  num_atoms: number;
  total_charge_mulliken: number;
  total_charge_lowdin: number;
}
```

---

## Dipole Moment

### `compute_dipole`

```typescript
function compute_dipole(
  elements: string,
  coords_flat: string
): string
```

Returns JSON:
```typescript
interface DipoleResult {
  vector: [number, number, number];  // Debye [x, y, z]
  magnitude: number;                 // Debye
  unit: 'Debye';
}
```

---

## Density of States

### `compute_dos`

```typescript
function compute_dos(
  elements: string,
  coords_flat: string,
  sigma: number,      // Gaussian smearing width (eV)
  e_min: number,      // energy range start (eV)
  e_max: number,      // energy range end (eV)
  n_points: number    // number of energy grid points
): string
```

Returns JSON:
```typescript
interface DosResult {
  sigma: number;
  energies: number[];    // eV axis
  total_dos: number[];   // DOS curve
  n_atoms_pdos: number;  // number of atoms with PDOS
}
```

```typescript
import { Chart } from 'chart.js';

const dos = JSON.parse(compute_dos(elements, coordsJson, 0.2, -20, 5, 200));

new Chart(ctx, {
  type: 'line',
  data: {
    labels: dos.energies,
    datasets: [{ data: dos.total_dos, label: 'DOS' }],
  },
});
```

---

## Molecular Alignment

### `compute_rmsd`

```typescript
function compute_rmsd(
  coords: string,     // JSON flat array
  reference: string   // JSON flat array
): string
```

Returns `{"rmsd": 0.034281}`.

---

## Force Field Energy

### `compute_uff_energy`

```typescript
function compute_uff_energy(smiles: string, coords: string): string
```

Returns `{"energy": 12.345, "unit": "kcal/mol"}` or `{"error": "..."}`.

---

## Spectroscopy (Track D)

### `compute_stda_uvvis`

```typescript
function compute_stda_uvvis(
  elements: string,      // JSON number[]
  coords_flat: string,   // JSON number[] flat [x0,y0,z0,...]
  sigma: number,
  e_min: number,
  e_max: number,
  n_points: number,
  broadening: string     // "gaussian" | "lorentzian"
): string
```

Returns JSON `StdaUvVisSpectrum`:

```typescript
interface StdaUvVisSpectrum {
  wavelengths_nm: number[];    // wavelength axis in nm
  intensities: number[];       // ε(λ) values
  excitations: {
    energy_ev: number;
    wavelength_nm: number;
    oscillator_strength: number;
    from_mo: number;
    to_mo: number;
  }[];
  homo_energy: number;         // eV
  lumo_energy: number;         // eV
  gap: number;                 // eV
  n_transitions: number;
  broadening_type: string;
  notes: string[];
}
```

```typescript
const r = JSON.parse(embed('c1ccccc1', 42));
const el = JSON.stringify(r.elements);
const co = JSON.stringify(r.coords);
const spec = JSON.parse(compute_stda_uvvis(el, co, 0.3, 1.0, 8.0, 500, 'gaussian'));
const maxIdx = spec.intensities.indexOf(Math.max(...spec.intensities));
console.log(`λ_max ≈ ${spec.wavelengths_nm[maxIdx].toFixed(1)} nm, gap ${spec.gap.toFixed(2)} eV`);
```

---

### `compute_vibrational_analysis`

```typescript
function compute_vibrational_analysis(
  elements: string,   // JSON number[]
  coords_flat: string, // JSON number[] flat
  method: string      // "eht" | "pm3" | "xtb"
): string
```

Returns JSON `VibrationalAnalysis`:

```typescript
interface VibrationalAnalysis {
  n_atoms: number;
  modes: {
    frequency_cm1: number;   // negative = imaginary (transition state)
    ir_intensity: number;    // km/mol
    displacement: number[];  // 3N normal-coordinate displacements
    is_real: boolean;
    label: string | null;    // functional-group annotation, e.g. "C=O stretch"
  }[];
  n_real_modes: number;
  zpve_ev: number;
  method: string;
  notes: string[];
  thermochemistry: {
    zpve_kcal: number;
    thermal_energy_kcal: number;
    entropy_vib_cal: number;
    gibbs_correction_kcal: number;
  };
}
```

```typescript
const vib = JSON.parse(compute_vibrational_analysis(el, co, 'xtb'));
console.log(`ZPVE: ${vib.zpve_ev.toFixed(4)} eV, ${vib.n_real_modes} real modes`);
console.log(`Gibbs correction: ${vib.thermochemistry.gibbs_correction_kcal.toFixed(3)} kcal/mol`);
```

### `compute_vibrational_analysis_uff`

```typescript
function compute_vibrational_analysis_uff(
  smiles: string,
  elements: string,
  coords_flat: string,
  step_size: number
): string
```

Fast vibrational-analysis path backed by the UFF analytical Hessian. This is useful for interactive IR previews where the numerical PM3/xTB Hessian would be too slow.

```typescript
const vibFast = JSON.parse(compute_vibrational_analysis_uff('CCO', el, co, 0.005));
console.log(vibFast.method); // "UFF"
```

---

### `compute_ir_spectrum`

```typescript
function compute_ir_spectrum(
  elements: string,    // JSON number[]
  coords_flat: string, // JSON number[] flat
  method: string       // "eht" | "pm3" | "xtb"
): string
```

Convenience wrapper: runs vibrational analysis + default Lorentzian broadening (γ = 15 cm⁻¹, 600–4000 cm⁻¹, 1000 points).

Returns JSON `IrSpectrum { wavenumbers, intensities, transmittance, peaks, gamma, notes }`.

```typescript
const spec = JSON.parse(compute_ir_spectrum(el, co, 'xtb'));
const peak_wn = spec.wavenumbers[spec.intensities.indexOf(Math.max(...spec.intensities))];
console.log(`Dominant IR band: ${peak_wn.toFixed(1)} cm⁻¹`);
```

---

### `compute_ir_spectrum_broadened`

```typescript
function compute_ir_spectrum_broadened(
  analysis_json: string,  // JSON-serialised VibrationalAnalysis
  gamma: number,          // HWHM (Lorentzian) or σ (Gaussian) in cm⁻¹
  wn_min: number,
  wn_max: number,
  n_points: number,
  broadening: string      // "lorentzian" | "gaussian"
): string
```

Full-control IR broadening from an existing `VibrationalAnalysis` JSON. Returns the same `IrSpectrum` shape with both `.intensities` (absorbance) and `.transmittance` axes.

```typescript
const spec = JSON.parse(
  compute_ir_spectrum_broadened(JSON.stringify(vib), 10.0, 400.0, 4000.0, 2000, 'gaussian')
);
console.log(`${spec.peaks.length} labelled peaks`);
```

---

### `predict_nmr_shifts`

```typescript
function predict_nmr_shifts(smiles: string): string
```

Returns JSON `NmrShiftResult`:

```typescript
interface NmrShiftResult {
  h_shifts: { atom_index: number; element: number; shift_ppm: number; environment: string; confidence: number }[];
  c_shifts: { atom_index: number; element: number; shift_ppm: number; environment: string; confidence: number }[];
  other_shifts: { nucleus: string; shifts: { atom_index: number; element: number; shift_ppm: number; environment: string; confidence: number }[] }[];
  notes: string[];
}
```

```typescript
const shifts = JSON.parse(predict_nmr_shifts('CCO'));
shifts.h_shifts.forEach((h: any) =>
  console.log(`H#${h.atom_index}: ${h.shift_ppm.toFixed(2)} ppm [${h.environment}]`));
```

---

### `predict_nmr_shifts_for_nucleus`

```typescript
function predict_nmr_shifts_for_nucleus(smiles: string, nucleus: string): string
```

Returns JSON `ChemicalShift[]` for the requested nucleus. Use this for the expanded registry such as `2H`, `35Cl`, `79Br`, `195Pt`, or `207Pb`.

```typescript
const cl35 = JSON.parse(predict_nmr_shifts_for_nucleus('[Cl]', '35Cl'));
console.log(cl35.map((peak: any) => peak.shift_ppm.toFixed(1)).join(', '));
```

---

### `compute_giao_nmr`

```typescript
function compute_giao_nmr(
  elements_json: string,
  coords_flat_json: string,
  nucleus: string,
  charge: number,
  multiplicity: number,
  max_scf_iter: number,
  allow_basis_fallback: boolean,
  use_parallel_eri: boolean
): string
```

Runs the public SCF-backed GIAO route for one requested nucleus and returns JSON `GiaoNmrResult { chemical_shifts, shieldings, scf_converged, scf_iterations, fallback_elements, notes, ... }`.

The current public path is singlet closed-shell only. By default it rejects elements that only have the hydrogen-like fallback basis in the SCF stack.

```typescript
const waterElements = JSON.stringify([8, 1, 1]);
const waterCoords = JSON.stringify([
  0.0, 0.0, 0.1173,
  0.0, 0.7572, -0.4692,
  0.0, -0.7572, -0.4692,
]);
const giao = JSON.parse(compute_giao_nmr(waterElements, waterCoords, '1H', 0, 1, 100, false, false));
console.log(giao.chemical_shifts);
```

---

### `predict_nmr_couplings`

```typescript
function predict_nmr_couplings(
  smiles: string,
  coords_flat: string   // JSON number[] or "[]" for free-rotation average
): string
```

Returns JSON array of `JCoupling { h1_index, h2_index, j_hz, n_bonds, coupling_type }`. Parameterized pathways: H-C-C-H (Altona-Sundaralingam), H-C-N-H (Bystrov), H-C-O-H, H-C-S-H.

---

### `compute_ensemble_j_couplings`

```typescript
function compute_ensemble_j_couplings(
  smiles: string,
  conformers_json: string,  // JSON number[][] (array of flat coord arrays)
  energies_json: string,    // JSON number[] in kcal/mol
  temperature_k: number
): string
```

Boltzmann-average ³J couplings over a conformer ensemble. Returns JSON array of `JCoupling`.

```typescript
const allCoords = results.map((r: any) => r.coords);
const energies = new Array(allCoords.length).fill(0);
const couplings = JSON.parse(
  compute_ensemble_j_couplings('CCCC', JSON.stringify(allCoords), JSON.stringify(energies), 298.15)
);
couplings.forEach((jc: any) =>
  console.log(`H${jc.h1_index}–H${jc.h2_index}: ${jc.j_hz.toFixed(2)} Hz`));
```

---

### `compute_nmr_spectrum`

```typescript
function compute_nmr_spectrum(
  smiles: string,
  nucleus: string,    // e.g. "1H" | "13C" | "35Cl" | "79Br" | "195Pt"
  gamma: number,
  ppm_min: number,
  ppm_max: number,
  n_points: number
): string
```

Full NMR spectrum pipeline. Returns JSON `NmrSpectrum { ppm_axis, intensities, peaks, nucleus, gamma }`.

The ¹H path still has the richest splitting model. Other nuclei use fast relative inference with nucleus-specific linewidth defaults and are intended for quick qualitative inspection.

```typescript
const spec = JSON.parse(compute_nmr_spectrum('CCO', '1H', 0.01, -2.0, 12.0, 2000));
const maxIdx = spec.intensities.indexOf(Math.max(...spec.intensities));
console.log(`Most intense peak: δ ${spec.ppm_axis[maxIdx].toFixed(2)} ppm`);
```

### `compute_nmr_spectrum_with_coords`

```typescript
function compute_nmr_spectrum_with_coords(
  smiles: string,
  coords_flat: string,
  nucleus: string,
  gamma: number,
  ppm_min: number,
  ppm_max: number,
  n_points: number
): string
```

Coordinate-aware variant of the NMR spectrum API. When 3D coordinates are provided, vicinal ³J couplings use the Karplus equation instead of the topology-only fallback.

```typescript
const spec3d = JSON.parse(
  compute_nmr_spectrum_with_coords('CCO', co, '1H', 0.01, -2.0, 12.0, 2000)
);
console.log(`${spec3d.peaks.length} resolved peaks from the 3D coupling model`);
```

---

### `compute_hose_codes`

```typescript
function compute_hose_codes(smiles: string, max_radius: number): string
```

Returns JSON array of `HoseCode { atom_index, element, code_string, radius }`.

---

## Materials

### `create_unit_cell`

```typescript
function create_unit_cell(
  a: number, b: number, c: number,
  alpha: number, beta: number, gamma: number
): string
```

Returns JSON: `{"a":..., "b":..., "c":..., "alpha":..., "beta":..., "gamma":..., "volume":...}`.

### `assemble_framework`

```typescript
function assemble_framework(
  topology: string,   // "pcu" | "dia" | "sql"
  metal: number,      // atomic number
  geometry: string,   // "linear" | "trigonal" | "tetrahedral" | "square_planar" | "octahedral"
  lattice_a: number,  // Å
  supercell: number   // replication factor (1 = no), number
): string
```

Returns JSON `CrystalStructure` with `elements`, `cart_coords`, `frac_coords`, `labels`, `lattice`.

```typescript
const mof = JSON.parse(assemble_framework("pcu", 30, "octahedral", 26.3, 1));
console.log(`${mof.num_atoms} atoms in framework`);
```

---

## Web Workers / Batch Transport

### `pack_batch_arrow`

```typescript
function pack_batch_arrow(results_json: string): string
```

Pack a JSON array of conformer results into Arrow-compatible columnar format for efficient Worker postMessage transfer.

### `split_worker_tasks`

```typescript
function split_worker_tasks(
  smiles_json: string,  // JSON array of SMILES
  n_workers: number,
  seed: number
): string
```

Split SMILES into balanced task batches for Web Workers.

### `estimate_workers`

```typescript
function estimate_workers(n_items: number, max_workers: number): number
```

Heuristic to pick the optimal worker count.

---

## Semi-Empirical QM

### `compute_pm3`

```typescript
function compute_pm3(elements: string, coords_flat: string): string
```

PM3 NDDO semi-empirical SCF. Supports PM3(tm) transition metals.

**Returns JSON:** `{orbital_energies, electronic_energy, nuclear_repulsion, total_energy, heat_of_formation, n_basis, n_electrons, homo_energy, lumo_energy, gap, mulliken_charges, scf_iterations, converged}`

```typescript
const pm3 = JSON.parse(compute_pm3(elements, coords));
console.log(`HOF: ${pm3.heat_of_formation.toFixed(2)} kcal/mol`);
```

### `compute_xtb`

```typescript
function compute_xtb(elements: string, coords_flat: string): string
```

GFN0-xTB tight-binding with SCC. Supports 25 elements.

**Returns JSON:** `{orbital_energies, electronic_energy, repulsive_energy, total_energy, n_basis, n_electrons, homo_energy, lumo_energy, gap, mulliken_charges, scc_iterations, converged}`

### `compute_gfn1`

```typescript
function compute_gfn1(elements: string, coords_flat: string): string
```

GFN1-xTB tight-binding with shell-resolved charges and D3-style dispersion.

**Returns JSON:** `{orbital_energies, electronic_energy, repulsive_energy, dispersion_energy, total_energy, n_basis, n_electrons, homo_energy, lumo_energy, gap, mulliken_charges, shell_charges, scc_iterations, converged}`

### `compute_gfn2`

```typescript
function compute_gfn2(elements: string, coords_flat: string): string
```

GFN2-xTB tight-binding with multipole electrostatics, D4-style dispersion, and halogen-bond corrections.

**Returns JSON:** `{orbital_energies, electronic_energy, repulsive_energy, dispersion_energy, halogen_bond_energy, total_energy, n_basis, n_electrons, homo_energy, lumo_energy, gap, mulliken_charges, atomic_dipoles, atomic_quadrupoles, scc_iterations, converged}`

---

## Bond Orders & Reactivity

### `compute_bond_orders`

```typescript
function compute_bond_orders(elements: string, coords_flat: string): string
```

**Returns JSON:** `{atom_pairs, distances, wiberg, mayer, wiberg_valence, mayer_valence}`

### `compute_frontier_descriptors`

```typescript
function compute_frontier_descriptors(elements: string, coords_flat: string): string
```

**Returns JSON:** `{homo_atom_contributions, lumo_atom_contributions, dual_descriptor, homo_energy, lumo_energy, gap}`

### `compute_fukui_descriptors`

```typescript
function compute_fukui_descriptors(elements: string, coords_flat: string): string
```

**Returns JSON:** `{condensed_atom_indices, condensed_f_plus, condensed_f_minus, condensed_f_radical, dual_descriptor, gap}`

### `compute_reactivity_ranking`

```typescript
function compute_reactivity_ranking(elements: string, coords_flat: string): string
```

**Returns JSON:** `{nucleophilic_attack_sites, electrophilic_attack_sites, radical_attack_sites}`

### `compute_empirical_pka`

```typescript
function compute_empirical_pka(smiles: string): string
```

**Returns JSON:** `{acidic_sites, basic_sites}`

---

## ML Properties

### `compute_ml_properties`

```typescript
function compute_ml_properties(smiles: string): string
```

**Returns JSON:** `{logp, molar_refractivity, log_solubility, lipinski_violations, lipinski_passes, druglikeness}`

### `compute_molecular_descriptors`

```typescript
function compute_molecular_descriptors(smiles: string): string
```

**Returns JSON:** `{molecular_weight, n_heavy_atoms, n_hydrogens, n_bonds, n_rotatable_bonds, n_hbd, n_hba, fsp3, wiener_index, n_rings, n_aromatic, balaban_j, ...}`

---

## Stereochemistry

### `analyze_stereo`

```typescript
function analyze_stereo(
  smiles: string,
  coords_flat: string  // JSON number[] flat, or "[]" for topology-only
): string
```

Assign CIP priorities, R/S at stereocenters, and E/Z at double bonds. Returns JSON:

```typescript
interface StereoAnalysis {
  stereocenters: {
    atom_index: number;
    element: number;
    substituent_indices: number[];
    priorities: number[];          // CIP rank (1 = highest)
    configuration: 'R' | 'S' | null;
  }[];
  double_bonds: {
    atom1: number;
    atom2: number;
    configuration: 'E' | 'Z' | null;
    high_priority_sub1: number | null;
    high_priority_sub2: number | null;
  }[];
  n_stereocenters: number;
  n_double_bonds: number;
}
```

```typescript
const r = JSON.parse(embed('C(F)(Cl)(Br)I', 42));
const stereo = JSON.parse(analyze_stereo('C(F)(Cl)(Br)I', JSON.stringify(r.coords)));
console.log(`${stereo.n_stereocenters} stereocenters`);
stereo.stereocenters.forEach((sc: any) =>
  console.log(`  Atom ${sc.atom_index}: ${sc.configuration}`));
```

---

## Implicit Solvation

### `compute_nonpolar_solvation`

```typescript
function compute_nonpolar_solvation(
  elements: string,     // JSON number[]
  coords_flat: string,  // JSON number[] flat
  probe_radius: number  // Å, typically 1.4
): string
```

Returns JSON:

```typescript
interface NonPolarSolvation {
  energy_kcal_mol: number;
  atom_contributions: number[];
  atom_sasa: number[];
  total_sasa: number;
}
```

### `compute_gb_solvation`

```typescript
function compute_gb_solvation(
  elements: string,          // JSON number[]
  coords_flat: string,       // JSON number[] flat
  charges: string,           // JSON number[]
  solvent_dielectric: number, // default 78.5
  solute_dielectric: number,  // default 1.0
  probe_radius: number        // default 1.4
): string
```

GB/SA electrostatic solvation. Returns JSON:

```typescript
interface GbSolvation {
  electrostatic_energy_kcal_mol: number;
  nonpolar_energy_kcal_mol: number;
  total_energy_kcal_mol: number;
  born_radii: number[];
  charges: number[];
  solvent_dielectric: number;
  solute_dielectric: number;
}
```

```typescript
const r = JSON.parse(embed('CCO', 42));
const el        = JSON.stringify(r.elements);
const co        = JSON.stringify(r.coords);
const chg       = JSON.parse(compute_charges('CCO')).charges;
const chgJson   = JSON.stringify(chg);

const np = JSON.parse(compute_nonpolar_solvation(el, co, 1.4));
console.log(`Non-polar ΔG: ${np.energy_kcal_mol.toFixed(2)} kcal/mol`);

const gb = JSON.parse(compute_gb_solvation(el, co, chgJson, 78.5, 1.0, 1.4));
console.log(`Total solvation: ${gb.total_energy_kcal_mol.toFixed(2)} kcal/mol`);
```

---

## Ring Perception

### `compute_sssr`

```typescript
function compute_sssr(smiles: string): string
```

Smallest Set of Smallest Rings (Horton's algorithm). Returns JSON:

```typescript
interface SssrResult {
  rings: { atoms: number[]; size: number; is_aromatic: boolean }[];
  atom_ring_count: number[];
  atom_ring_sizes: number[][];
  ring_size_histogram: Record<string, number>;
}
```

```typescript
const r = JSON.parse(compute_sssr('c1ccc2ccccc2c1'));  // naphthalene
console.log(`${r.rings.length} rings`);
r.rings.forEach((ring: any) =>
  console.log(`  ${ring.size}-membered, aromatic=${ring.is_aromatic}`));
```

---

## Fingerprints (ECFP)

### `compute_ecfp`

```typescript
function compute_ecfp(
  smiles: string,
  radius: number,  // ECFP4 = 2, ECFP6 = 3
  n_bits: number   // 1024 or 2048
): string
```

Extended-Connectivity Fingerprint (Morgan algorithm). Returns JSON:

```typescript
interface ECFPFingerprint {
  n_bits: number;
  on_bits: number[];   // set bit indices
  radius: number;
  raw_features: number[];
}
```

### `compute_tanimoto`

```typescript
function compute_tanimoto(
  fp1_json: string,
  fp2_json: string
): string
```

Returns JSON `{ "tanimoto": number }`.

```typescript
const fp1 = compute_ecfp('c1ccccc1', 2, 2048);
const fp2 = compute_ecfp('Cc1ccccc1', 2, 2048);
const t   = JSON.parse(compute_tanimoto(fp1, fp2));
console.log(`Tanimoto: ${t.tanimoto.toFixed(3)}`);  // ~0.5–0.7
```

---

## Conformer Clustering

### `butina_cluster`

```typescript
function butina_cluster(
  conformers_json: string,  // JSON number[][] (array of flat coord arrays)
  rmsd_cutoff: number
): string
```

Taylor-Butina single-linkage RMSD clustering. Returns JSON:

```typescript
interface ClusterResult {
  n_clusters: number;
  assignments: number[];       // cluster index per conformer
  centroid_indices: number[];  // representative conformer per cluster
  cluster_sizes: number[];
  rmsd_cutoff: number;
}
```

### `compute_rmsd_matrix`

```typescript
function compute_rmsd_matrix(conformers_json: string): string
// Returns JSON number[][] — pairwise RMSD matrix
```

```typescript
const coords = results.map((r: any) => r.coords);
const cl = JSON.parse(butina_cluster(JSON.stringify(coords), 1.0));
console.log(`${cl.n_clusters} clusters from ${coords.length} conformers`);
```

---

## Usage Patterns

### High-Throughput Conformer Generation

```typescript
import init, { embed_coords_typed, estimate_workers, split_worker_tasks } from 'sci-form-wasm';

await init();

// Single molecule — typed array (fastest)
const coords = embed_coords_typed('c1ccccc1', 42);
// coords[0], coords[1], coords[2] = x₀, y₀, z₀ of atom 0

// Check auto worker estimate
const nWorkers = estimate_workers(1000, navigator.hardwareConcurrency ?? 4);
const tasks = JSON.parse(split_worker_tasks(
  JSON.stringify(smilesList), nWorkers, 42
));
```

### Three.js Molecule Viewer

```typescript
import init, { embed, eht_calculate, eht_orbital_mesh } from 'sci-form-wasm';
import * as THREE from 'three';

await init();

const smiles = 'c1ccccc1';
const result = JSON.parse(embed(smiles, 42));
const elements = JSON.stringify(Array.from(result.elements));
const coordsJson = JSON.stringify(result.coords);

// Atom spheres
const scene = new THREE.Scene();
for (let i = 0; i < result.num_atoms; i++) {
  const sphere = new THREE.Mesh(
    new THREE.SphereGeometry(0.3),
    new THREE.MeshPhongMaterial({ color: 0x336699 })
  );
  sphere.position.set(result.coords[i*3], result.coords[i*3+1], result.coords[i*3+2]);
  scene.add(sphere);
}

// HOMO orbital
const eht = JSON.parse(eht_calculate(elements, coordsJson, 0.0));
const mesh = JSON.parse(eht_orbital_mesh(elements, coordsJson, eht.homo_index, 0.2, 0.02));

const orbGeom = new THREE.BufferGeometry();
orbGeom.setAttribute('position', new THREE.Float32BufferAttribute(mesh.vertices, 3));
orbGeom.setAttribute('normal', new THREE.Float32BufferAttribute(mesh.normals, 3));
orbGeom.setIndex(mesh.indices);
scene.add(new THREE.Mesh(orbGeom, new THREE.MeshPhongMaterial({
  color: 0x0077ff, transparent: true, opacity: 0.5, side: THREE.DoubleSide
})));
```

### React Hook

```tsx
import { useEffect, useState } from 'react';
import init, { embed } from 'sci-form-wasm';

export function useConformer(smiles: string) {
  const [result, setResult] = useState<any>(null);

  useEffect(() => {
    (async () => {
      await init();
      setResult(JSON.parse(embed(smiles, 42)));
    })();
  }, [smiles]);

  return result;
}
```

