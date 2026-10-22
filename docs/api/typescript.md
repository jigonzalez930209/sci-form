# TypeScript / JavaScript API Reference

The npm package is published as **`sci-form-wasm`**.

```bash
npm install sci-form-wasm
```

It is a dual package: `import` resolves to the ESM (browser/bundler) build and
`require` resolves to the CommonJS (Node.js) build. Both halves expose the
same API.

Every entry point returns a **JSON string**, except the typed-array variants
(`embed_coords_typed`, `embed_coords_typed_seed64`, `compute_esp_grid_typed`,
`eht_orbital_grid_typed`, `eht_orbital_grid_from_coefficients_typed`), which return a typed
array, and `estimate_workers`, which returns a number.

---

## Initialisation

Browser / bundler (ESM):

```typescript
import init, { embed, solve_eht } from 'sci-form-wasm';
await init();
```

Node.js (CommonJS) needs no async initialisation:

```javascript
const sci = require('sci-form-wasm');
```

Both entry points come from the same tarball: the `require` condition of the
package points at the CommonJS half in `cjs/`, which loads its own `.wasm`
file from disk.

The `alpha` and `beta` research modules are published under the `sci-form-wasm/alpha` and
`sci-form-wasm/beta` sub-paths.

---

## Conventions

### Error payload

On failure **every** entry point returns the canonical payload

```json
{ "error": "<message>" }
```

and nothing else. The payload is always valid JSON — a NaN or a bare Rust message is never
emitted. So a caller can always write:

```typescript
const parsed = JSON.parse(sci.compute_something(...));
if (parsed.error) { /* handle */ }
```

### Seeds are 64-bit

`wasm-bindgen` maps a Rust `u64` to a JavaScript `BigInt`, so every seed-taking entry point has
a **`*_seed64` variant** that takes the full 64-bit range. The older 32-bit variants are kept
for compatibility; they no longer truncate silently — a seed above 2³² is rejected.

The `*_seed64` variants are:

| 64-bit variant | 32-bit counterpart |
|---|---|
| `embed_seed64` | `embed` |
| `embed_coords_seed64` | `embed_coords` |
| `embed_coords_typed_seed64` | `embed_coords_typed` |
| `embed_batch_seed64` | `embed_batch` |
| `compute_md_trajectory_seed64` | `compute_md_trajectory` |
| `compute_md_trajectory_nvt_seed64` | `compute_md_trajectory_nvt` |
| `compute_md_trajectory_with_method_seed64` | `compute_md_trajectory_with_method` |
| `search_conformers_with_uff_seed64` | `search_conformers_with_uff` |
| `split_worker_tasks_seed64` | `split_worker_tasks` |

`embed_diverse`, `train_random_forest` and `train_gradient_boosting` only exist in the 64-bit
form.

### Units

| Quantity | Unit |
|---|---|
| Geometry (coordinates, bond lengths, cell parameters) | Å (cell angles in degrees) |
| UFF / MMFF94 energies | kcal/mol |
| EHT, PM3, xTB orbital energies and HOMO–LUMO gaps | eV |
| PM3 / xTB total, electronic, repulsive, dispersion energies | eV |
| PM3 heat of formation | kcal/mol |
| HF-3c energy, HF energy, nuclear repulsion, D3/gCP/SRB | Hartree |
| HF-3c orbital energies, HOMO/LUMO, gap | Hartree |
| ANI energy, per-atom energies | Hartree (forces Hartree/Å) |
| CISD `energy`, `correlation_energy` | Hartree (`energy_ev` is eV) |
| EHT gradients, `optimize_geometry_eht` energy | eV, eV/Å |
| `optimize_framework` energy and forces | eV, eV/Å |
| Vibrational frequencies, IR wavenumbers | cm⁻¹ |
| NMR chemical shifts and shieldings | ppm |
| NMR J-couplings | Hz |
| Dipole moment | Debye |
| Partial charges | elementary charge (e) |
| Electrostatic potential | Hartree/e |
| DOS / PDOS | 1/eV |
| Surface area | Å² |
| MD energies | kcal/mol; times in fs; temperatures in K |

Energies are always total system energies unless the field name says otherwise.

### Naming

The exported names match the **Rust root function names**, so the same capability has the same
name in Rust, Python, WASM and (where exposed) the CLI. A handful of names are kept as aliases
for the older spelling: `eht_calculate` (alias of `solve_eht`), `eht_or_uff_fallback` (alias of
`compute_eht_or_uff_fallback`), `get_eht_support` / `get_system_capabilities` /
`get_system_method_plan` (aliases of the `get_`-less names), `pack_batch_arrow` (alias of
`pack_conformers`), `compute_ml_properties` and `compute_molecular_descriptors` (deprecated
aliases of `predict_ml_properties` and `compute_ml_descriptors`).

---

## Metadata

### `version`

```typescript
function version(): string
```

Returns `"sci-form 1.0.0"`.

### `alpha_modules_info` / `beta_modules_info`

```typescript
function alpha_modules_info(): string
function beta_modules_info(): string
```

Returns JSON `{"alpha_modules": [...], "stability": "...", "version": "..."}` and the `beta_`
equivalent.

---

## Conformer generation

### `embed`

```typescript
function embed(smiles: string, seed: number): string
```

JSON payload:

```typescript
interface ConformerResult {
  smiles: string;
  num_atoms: number;
  coords: number[];            // flat [x0,y0,z0,...], Å
  elements: number[];          // atomic numbers
  bonds: [number, number, string][];
  error: string | null;
  time_ms: number;   // always 0 in the WASM build (no wall-clock timer)
}
```

```typescript
const result = JSON.parse(embed('CCO', 42));
console.log(result.num_atoms, result.coords.length);
```

### `embed_seed64`

```typescript
function embed_seed64(smiles: string, seed: bigint): string
```

Same payload; `seed` is a `BigInt` (e.g. `42n`). Numbers are also accepted below 2⁵³.

### `embed_coords`

```typescript
function embed_coords(smiles: string, seed: number): string
```

Compact JSON: `{"coords": number[], "num_atoms": number}`.

### `embed_coords_seed64`

```typescript
function embed_coords_seed64(smiles: string, seed: bigint): string
```

### `embed_coords_typed`

```typescript
function embed_coords_typed(smiles: string, seed: number): Float64Array
```

The coordinates as a `Float64Array` — **no JSON parsing overhead**. Preferred for
high-throughput use. Empty on failure.

```typescript
const coords: Float64Array = embed_coords_typed('CCO', 42);
```

### `embed_coords_typed_seed64`

```typescript
function embed_coords_typed_seed64(smiles: string, seed: bigint): Float64Array
```

### `embed_batch`

```typescript
function embed_batch(smiles_list: string, seed: number): string
```

Batch-embed from a newline-separated SMILES string. Returns a JSON array of
`ConformerResult`. A batch-level failure returns a JSON array holding a single
`{"error": "..."}` object, so the array shape is preserved.

### `embed_batch_seed64`

```typescript
function embed_batch_seed64(smiles_list: string, seed: bigint): string
```

### `embed_diverse`

```typescript
function embed_diverse(smiles: string, n_conformers: number, rmsd_cutoff: number, base_seed: bigint): string
```

Embeds `n_conformers` times with consecutive seeds and keeps the Butina RMSD cluster centroids
(cutoff in Å). Returns a JSON array of conformer records (`seed`, `cluster_id`, `coords`,
`energy_kcal_mol`).

### `parse_smiles`

```typescript
function parse_smiles(smiles: string): string
```

Returns JSON:

```typescript
interface ParseResult {
  num_atoms: number;
  num_bonds: number;
  atoms: { element: number; hybridization: string; formal_charge: number }[];
}
```

The `atoms` array matches the payload of Python's `parse()`.

### `to_canonical_smiles`

```typescript
function to_canonical_smiles(smiles: string): string
```

Returns the canonical SMILES as a plain string, e.g. `"BrC(Cl)(F)I"`.

---

## Charges, surface and electrostatics

### `compute_charges`

```typescript
function compute_charges(smiles: string): string
```

```typescript
interface ChargeResult {
  charges: number[];        // e
  iterations: number;
  total_charge: number;     // e
  converged: boolean;
}
```

### `compute_charges_configured`

```typescript
function compute_charges_configured(
  smiles: string,
  max_iter: number,
  initial_damping: number,
  convergence_threshold: number,
): string
```

Same `ChargeResult` payload with an explicit convergence configuration.

### `compute_sasa`

```typescript
function compute_sasa(elements: string, coords_flat: string, probe_radius: number): string
```

```typescript
interface SasaResult {
  total_sasa: number;    // Å²
  atom_sasa: number[];   // Å²
  probe_radius: number;  // Å
  num_points: number;
}
```

### `compute_population`

```typescript
function compute_population(elements: string, coords_flat: string): string
```

```typescript
interface PopulationResult {
  mulliken_charges: number[];      // e
  lowdin_charges: number[];        // e
  mulliken_populations: number[];
  num_atoms: number;
  total_charge_mulliken: number;
  total_charge_lowdin: number;
  charge_conservation_error: number;
  homo_energy: number;             // eV
  lumo_energy: number;             // eV
  gap: number;                     // eV
  warnings: string[];
}
```

### `compute_dipole`

```typescript
function compute_dipole(elements: string, coords_flat: string): string
```

```typescript
interface DipoleResult {
  vector: [number, number, number];  // Debye
  magnitude: number;                 // Debye
  unit: 'Debye';
  nuclear_dipole: [number, number, number];
  electronic_dipole: [number, number, number];
}
```

### `compute_esp_grid`

```typescript
function compute_esp_grid(elements: string, coords_flat: string, spacing: number, padding: number): string
```

Electrostatic potential on a 3D grid, in **Hartree/e**. `compute_esp` is the
historical alias of this function.

### `compute_esp_grid_info`

```typescript
function compute_esp_grid_info(elements: string, coords_flat: string, spacing: number, padding: number): string
```

Returns `{"dims": [...], "origin": [...], "spacing": ...}` — the grid metadata without the
values.

### `compute_esp_grid_typed`

```typescript
function compute_esp_grid_typed(elements: string, coords_flat: string, spacing: number, padding: number): Float64Array
```

The grid values as a `Float64Array` of length `dims[0] * dims[1] * dims[2]` (Hartree/e). Pair
with `compute_esp_grid_info` for the metadata.

---

## Density of states and alignment

### `compute_dos`

```typescript
function compute_dos(elements: string, coords_flat: string, sigma: number, e_min: number, e_max: number, n_points: number): string
```

```typescript
interface DosResult {
  sigma: number;        // eV
  energies: number[];   // eV
  total_dos: number[];  // 1/eV
  n_atoms_pdos: number;
}
```

### `compute_dos_multimethod`

```typescript
function compute_dos_multimethod(elements: string, coords_flat: string, method: string, sigma: number, e_min: number, e_max: number, n_points: number): string
```

`method` is `eht`, `pm3`, `xtb` (GFN2), `gfn1`, `gfn2` or `hf3c`.

```typescript
interface MultiMethodDosResult {
  energies: number[];         // eV
  total_dos: number[];        // 1/eV
  method: string;
  homo_energy: number;        // eV
  lumo_energy: number;        // eV
  gap: number;                // eV
  orbital_energies: number[]; // eV
  sigma: number;              // eV
}
```

### `compute_rmsd`

```typescript
function compute_rmsd(coords: string, reference: string): string
```

Returns `{"rmsd": 0.034281}` (Å). A NaN RMSD is never emitted.

### `compute_rmsd_matrix`

```typescript
function compute_rmsd_matrix(conformers_json: string): string
```

Returns a JSON 2-D array of pairwise RMSD values (Å).

---

## Extended Hückel theory

### `eht_calculate`

```typescript
function eht_calculate(elements: string, coords_flat: string, k: number): string
```

`k` is the Wolfsberg-Helmholtz constant; `0.0` selects the default 1.75.

```typescript
interface EhtResult {
  energies: number[];      // eV, all MO energies
  coefficients: number[][];
  n_electrons: number;
  homo_index: number;
  lumo_index: number;
  homo_energy: number;     // eV
  lumo_energy: number;     // eV
  gap: number;             // eV
  support: EhtSupport;
}

interface EhtSupport {
  level: string;
  has_transition_metals: boolean;
  supported_elements: number[];
  provisional_elements: number[];
  unsupported_elements: number[];
  warnings: string[];
}
```

### `solve_eht`

```typescript
function solve_eht(elements: string, coords_flat: string, k: number): string
```

The canonical name of the same calculation (`sci_form::eht::solve_eht`); `eht_calculate` is
kept as an alias.

### `eht_support` / `get_eht_support`

```typescript
function eht_support(elements: string): string
function get_eht_support(elements: string): string
```

Returns the `EhtSupport` payload. `get_eht_support` is the Rust-root spelling, kept as an
alias.

### `compute_eht_or_uff_fallback` / `eht_or_uff_fallback`

```typescript
function compute_eht_or_uff_fallback(smiles: string, elements: string, coords_flat: string, allow_experimental_eht: boolean): string
function eht_or_uff_fallback(smiles: string, elements: string, coords_flat: string, allow_experimental_eht: boolean): string
```

Returns `{"mode": "eht" | "uff", "result": {...}}`.

### `compute_eht_gradient`

```typescript
function compute_eht_gradient(elements: string, coords_flat: string): string
```

```typescript
interface EhtGradient {
  gradients: [number, number, number][];  // eV/Å
  rms_gradient: number;                   // eV/Å
  max_gradient: number;                   // eV/Å
  energy: number;                         // eV
}
```

### `optimize_geometry_eht`

```typescript
function optimize_geometry_eht(elements: string, coords_flat: string, max_steps: number, grad_threshold: number, energy_threshold: number, step_size: number): string
```

```typescript
interface EhtOptResult {
  positions: [number, number, number][];  // Å
  energy: number;                         // eV
  n_steps: number;
  converged: boolean;
  rms_gradient: number;                   // eV/Å
  energies: number[];                     // eV
}
```

---

## Electronic structure

### `compute_pm3`

```typescript
function compute_pm3(elements_json: string, coords_flat_json: string): string
```

```typescript
interface Pm3Result {
  orbital_energies: number[];   // eV
  electronic_energy: number;    // eV
  nuclear_repulsion: number;    // eV
  total_energy: number;         // eV
  heat_of_formation: number;    // kcal/mol
  n_basis: number;
  n_electrons: number;
  homo_energy: number;          // eV
  lumo_energy: number;          // eV
  gap: number;                  // eV
  mulliken_charges: number[];   // e
  scf_iterations: number;
  converged: boolean;
  notes: string[];
}
```

### `compute_xtb`

```typescript
function compute_xtb(elements_json: string, coords_flat_json: string): string
```

GFN0-xTB: `{orbital_energies, electronic_energy, repulsive_energy, total_energy, n_basis,
n_electrons, homo_energy, lumo_energy, gap, mulliken_charges, scc_iterations, converged,
notes}` — all energies in **eV**.

### `compute_gfn1`

```typescript
function compute_gfn1(elements_json: string, coords_flat_json: string): string
```

GFN1-xTB: as `compute_xtb` plus `dispersion_energy` (eV), `dispersion_e6`, `dispersion_e8`,
`halogen_bond_energy` and `shell_charges` (e, per atom and shell).

### `compute_gfn2`

```typescript
function compute_gfn2(elements_json: string, coords_flat_json: string): string
```

GFN2-xTB: as `compute_xtb` plus `dispersion_energy` (eV), `dispersion_2_body_energy` (eV),
`halogen_bond_energy` (eV), `atomic_dipoles` (Debye) and `atomic_quadrupoles`.

### `compute_hf3c` / `compute_hf3c_custom`

```typescript
function compute_hf3c(elements_json: string, coords_flat_json: string): string
function compute_hf3c_custom(elements_json: string, coords_flat_json: string, max_scf_iter: number, n_cis_states: number, corrections: boolean): string
```

```typescript
interface Hf3cResult {
  energy: number;             // Hartree
  hf_energy: number;          // Hartree
  nuclear_repulsion: number;  // Hartree
  d3_energy: number;          // Hartree
  gcp_energy: number;         // Hartree
  srb_energy: number;         // Hartree
  orbital_energies: number[]; // Hartree
  scf_iterations: number;
  converged: boolean;
  cis: { excitations: unknown[]; error: string | null };
  n_basis: number;
  n_electrons: number;
  homo_energy: number;        // Hartree
  lumo_energy: number;        // Hartree
  gap: number;                // Hartree
  mulliken_charges: number[]; // e
}
```

### `compute_ani`

```typescript
function compute_ani(elements_json: string, coords_flat_json: string): string
```

```typescript
interface AniResult {
  energy: number;             // Hartree
  forces: [number, number, number][];  // Hartree/Å
  species: number[];
  atomic_energies: number[];  // Hartree
  aevs: number[][] | null;
}
```

### `compute_cisd`

```typescript
function compute_cisd(elements_json: string, coords_flat_json: string, n_states: number): string
```

```typescript
interface CisdResult {
  excitations: CisdExcitation[];
  correlation_energy: number;  // Hartree
  n_csfs: number;
  approximation: string;       // "full" or "perturbative"
}

interface CisdExcitation {
  energy: number;             // Hartree
  energy_ev: number;          // eV
  wavelength_nm: number;
  oscillator_strength: number;
  dominant_transition: [number, number][];
  character: string;
  singles_weight: number;
}
```

### `compute_aevs_tm`

```typescript
function compute_aevs_tm(elements_json: string, coords_flat_json: string): string
```

`{"aev_length": number, "aevs": number[][], "n_atoms": number}`.

---

### `compute_uhf` / `compute_rohf`

```typescript
function compute_uhf(elements: string, coords_flat: string, charge: number, multiplicity: number): string
function compute_rohf(elements: string, coords_flat: string, charge: number, multiplicity: number): string
```

Unrestricted and restricted open-shell Hartree-Fock. The energies are in
**Hartree**, the gap in **eV**. The dense matrices of the result are not part
of the payload; the scalar results, the orbital energies and the Mulliken
charges are.

### `compute_pm3_gradient` / `compute_xtb_gradient`

```typescript
function compute_pm3_gradient(elements: string, coords_flat: string): string
function compute_xtb_gradient(elements: string, coords_flat: string): string
```

Analytical energy gradients in **eV/Å** per atom.

## Force fields

### `compute_uff_energy`

```typescript
function compute_uff_energy(smiles: string, coords: string): string
```

`{"energy": 6.993048, "unit": "kcal/mol"}`.

### `compute_uff_energy_with_aromatic_heuristics`

```typescript
function compute_uff_energy_with_aromatic_heuristics(smiles: string, coords: string): string
```

Returns the raw energy, the aromatic correction and the corrected energy, all in **kcal/mol**.

### `compute_mmff94_energy`

```typescript
function compute_mmff94_energy(smiles: string, coords: string): string
```

`{"energy": 3.058318, "unit": "kcal/mol"}`.

### `compute_empirical_pka`

```typescript
function compute_empirical_pka(smiles: string): string
```

`{"acidic_sites": [...], "basic_sites": [...], "notes": [...]}` — coarse screening hints, not
publication-grade predictions.

---

### `compute_uff_energy_breakdown` / `compute_mmff94_energy_breakdown`

```typescript
function compute_uff_energy_breakdown(smiles: string, coords_flat: string): string
function compute_mmff94_energy_breakdown(smiles: string, coords_flat: string): string
```

The force-field energy split by term type, in **kcal/mol**:

```json
{ "terms": [{ "term": "bond", "energy": -12.3 }, { "term": "vdw", "energy": 4.1 }],
  "unit": "kcal/mol" }
```

## Spectroscopy

### `compute_stda_uvvis`

```typescript
function compute_stda_uvvis(elements_json: string, coords_flat_json: string, sigma: number, e_min: number, e_max: number, n_points: number, broadening: string): string
```

```typescript
interface StdaUvVisSpectrum {
  energies_ev: number[];      // eV
  wavelengths_nm: number[];   // nm
  absorptivity: number[];
  excitations: {
    energy_ev: number;
    wavelength_nm: number;
    oscillator_strength: number;
    from_mo: number;
    to_mo: number;
  }[];
  sigma: number;               // eV
  broadening: string;          // "Gaussian" | "Lorentzian"
  notes: string[];
}
```

### `compute_fluorescence`

```typescript
function compute_fluorescence(elements_json: string, coords_flat_json: string, sigma: number, e_min: number, e_max: number, n_points: number, stokes_ev: number, refractive_index: number, dark_threshold: number, broadening: string): string
```

```typescript
interface FluorescenceSpectrum {
  energies_ev: number[];          // eV
  wavelengths_nm: number[];       // nm
  intensity: number[];
  absorption_energy_ev: number;   // eV
  emission_energy_ev: number;     // eV
  emission_wavelength_nm: number; // nm
  oscillator_strength: number;
  radiative_rate_hz: number;      // Hz
  radiative_lifetime_ns: number;  // ns
  kasha: boolean;
  from_mo: number;
  to_mo: number;
  notes: string[];
}
```

### `compute_vibrational_analysis`

```typescript
function compute_vibrational_analysis(elements_json: string, coords_flat_json: string, method: string, step_size: number): string
```

```typescript
interface VibrationalAnalysis {
  n_atoms: number;
  elements: number[];
  modes: {
    frequency_cm1: number;    // cm⁻¹
    ir_intensity: number;     // km/mol
    displacement: number[];
    is_real: boolean;
  }[];
  n_real_modes: number;
  n_imaginary_modes: number;
  zpve_ev: number;            // eV
  thermochemistry: {
    temperature_k: number;
    zpve_kcal: number;        // kcal/mol
    thermal_energy_kcal: number;
    enthalpy_correction_kcal: number;
    entropy_vib_cal: number;  // cal/mol/K
    gibbs_correction_kcal: number;
  };
  method: string;
  notes: string[];
}
```

::: warning
The EHT Hessian is a one-electron model with no electron-electron repulsion, so the surface has
no minimum and imaginary (negative) frequencies are expected. Use `pm3` or `xtb`.
:::

### `compute_vibrational_analysis_uff`

```typescript
function compute_vibrational_analysis_uff(smiles: string, elements_json: string, coords_flat_json: string, step_size: number): string
```

The same payload with `method: "UFF"` — the fast analytical-Hessian path, useful for
interactive IR previews.

### `compute_ir_spectrum`

```typescript
function compute_ir_spectrum(analysis_json: string, gamma: number, wn_min: number, wn_max: number, n_points: number): string
```

Broadens an **existing** vibrational-analysis JSON string (the output of
`compute_vibrational_analysis` or `compute_vibrational_analysis_uff`).

```typescript
interface IrSpectrum {
  wavenumbers: number[];   // cm⁻¹
  intensities: number[];   // km/mol
  peaks: { frequency_cm1: number; ir_intensity: number; mode_index: number }[];
  gamma: number;           // cm⁻¹
  notes: string[];
}
```

```typescript
const analysis = compute_vibrational_analysis(elements, coords, 'xtb', 0.005);
const ir = JSON.parse(compute_ir_spectrum(analysis, 15.0, 400.0, 4000.0, 1000));
```

### `compute_ir_spectrum_broadened`

```typescript
function compute_ir_spectrum_broadened(elements_json: string, coords_flat_json: string, method: string, gamma: number, wn_min: number, wn_max: number, n_points: number, broadening: string, step_size: number, smiles: string): string
```

Runs the vibrational analysis internally and broadens the result. `broadening` is `gaussian` or
`lorentzian`; `step_size` is the finite-difference step in Å (`0.0` uses the engine default);
`smiles` is required when `method` is `uff` and ignored otherwise.

### `assign_peaks`

```typescript
function assign_peaks(frequencies_json: string, intensities_json: string, elements_json: string, mode_displacements_json: string): string
```

`frequencies_json` and `intensities_json` are JSON arrays of numbers; `elements_json` is a JSON
array of atomic numbers; `mode_displacements_json` is a JSON array of modes (each an array of
`[x,y,z]` triples) — pass `""` to skip the normal-mode refinement.

Returns `{"assignments": [...], "functional_groups": [...]}`.

---

## NMR

### `predict_nmr_shifts`

```typescript
function predict_nmr_shifts(smiles: string): string
```

```typescript
interface NmrShiftResult {
  h_shifts: ChemicalShift[];
  c_shifts: ChemicalShift[];
  f_shifts: ChemicalShift[];
  p_shifts: ChemicalShift[];
  n_shifts: ChemicalShift[];
  b_shifts: ChemicalShift[];
  si_shifts: ChemicalShift[];
  se_shifts: ChemicalShift[];
  o_shifts: ChemicalShift[];
  s_shifts: ChemicalShift[];
  other_shifts: unknown[];
  notes: string[];
}

interface ChemicalShift {
  atom_index: number;
  element: number;
  shift_ppm: number;   // ppm
  environment: string;
  confidence: number;
}
```

### `predict_nmr_shifts_for_nucleus`

```typescript
function predict_nmr_shifts_for_nucleus(smiles: string, nucleus: string): string
```

Returns a JSON array of `ChemicalShift` for the requested nucleus — `2H`, `35Cl`, `79Br`,
`195Pt`, `207Pb`, …

### `predict_nmr_couplings`

```typescript
function predict_nmr_couplings(smiles: string, coords_flat_json: string): string
```

Returns a JSON array of `{h1_index, h2_index, j_hz, n_bonds, coupling_type}` — couplings in
**Hz**. Pass `"[]"` for the topology-only free-rotation average.

### `compute_ensemble_j_couplings`

```typescript
function compute_ensemble_j_couplings(smiles: string, conformer_coords_json: string, energies_kcal_json: string, temperature_k: number): string
```

J-couplings averaged over a conformer ensemble with Boltzmann weighting.
`conformer_coords_json` is a JSON array of conformers (each a flat xyz array) and
`energies_kcal_json` the matching energies in **kcal/mol**.

### `compute_nmr_spectrum`

```typescript
function compute_nmr_spectrum(smiles: string, nucleus: string, gamma: number, ppm_min: number, ppm_max: number, n_points: number): string
```

```typescript
interface NmrSpectrum {
  ppm_axis: number[];    // ppm
  intensities: number[];
  peaks: { shift_ppm: number; intensity: number; atom_index: number; atom_indices: number[]; multiplicity: string; environment: string }[];
  nucleus: string;
  gamma: number;         // ppm
  integrations: unknown[];
  notes: string[];
}
```

### `compute_nmr_spectrum_with_coords`

```typescript
function compute_nmr_spectrum_with_coords(smiles: string, coords_flat_json: string, nucleus: string, gamma: number, ppm_min: number, ppm_max: number, n_points: number): string
```

Coordinate-aware variant: with 3D coordinates the vicinal ³J couplings use the Karplus relation
instead of the topology-only fallback.

### `compute_giao_nmr`

```typescript
function compute_giao_nmr(elements_json: string, coords_flat_json: string, nucleus: string, charge: number, multiplicity: number, max_scf_iter: number, allow_basis_fallback: boolean, use_parallel_eri: boolean): string
```

```typescript
interface GiaoNmrResult {
  nucleus: string;
  target_atomic_number: number;
  method: string;             // "RHF/GIAO"
  basis_set: string;          // "STO-3G"
  charge: number;
  multiplicity: number;
  scf_converged: boolean;
  scf_iterations: number;
  total_energy_hartree: number;   // Hartree
  n_basis: number;
  n_target_atoms: number;
  chemical_shifts: number[];      // ppm
  shieldings: { atom_index: number; element: number; isotropic: number; anisotropy: number; tensor: number[]; chemical_shift: number }[];
  fallback_elements: number[];
  notes: string[];
}
```

The public path is restricted to singlet closed-shell systems.

### `compute_hose_codes`

```typescript
function compute_hose_codes(smiles: string, max_radius: number): string
```

Returns a JSON array of `{atom_index, element, spheres: string[], full_code}`.

---

## Solvation

### `compute_nonpolar_solvation`

```typescript
function compute_nonpolar_solvation(elements: string, coords_flat: string, probe_radius: number): string
```

```typescript
interface NonPolarSolvation {
  energy_kcal_mol: number;      // kcal/mol
  atom_contributions: number[]; // kcal/mol
  atom_sasa: number[];          // Å²
  total_sasa: number;           // Å²
}
```

### `compute_gb_solvation`

```typescript
function compute_gb_solvation(elements: string, coords_flat: string, charges: string, solvent_dielectric: number, solute_dielectric: number, probe_radius: number): string
```

Returns a JSON string. On invalid input the result is
`{"error": "<message>"}`: the coordinate array must hold `3 * atoms` values,
the charge array must hold one value per atom, a dielectric constant must be
finite and `>= 1.0` (`1.0` is the vacuum limit), and the probe radius must be
finite and `>= 0`.

```typescript
interface GbSolvation {
  electrostatic_energy_kcal_mol: number;  // kcal/mol
  nonpolar_energy_kcal_mol: number;       // kcal/mol
  total_energy_kcal_mol: number;          // kcal/mol
  born_radii: number[];                   // Å
  charges: number[];                      // e
  solvent_dielectric: number;
  solute_dielectric: number;
}
```

---

## Rings, fingerprints and clustering

### `compute_sssr`

```typescript
function compute_sssr(smiles: string): string
```

```typescript
interface SssrResult {
  rings: { atoms: number[]; size: number; is_aromatic: boolean }[];
  atom_ring_count: number[];
  atom_ring_sizes: number[][];
  ring_size_histogram: number[];
}
```

### `compute_ecfp`

```typescript
function compute_ecfp(smiles: string, radius: number, n_bits: number): string
```

```typescript
interface ECFPFingerprint {
  n_bits: number;
  on_bits: number[];
  radius: number;          // 2 = ECFP4, 3 = ECFP6
  raw_features: number[];
}
```

### `compute_ecfp_batch`

```typescript
function compute_ecfp_batch(smiles_list_json: string, radius: number, n_bits: number): string
```

Returns a JSON array of `ECFPFingerprint`, in input order.

### `compute_tanimoto`

```typescript
function compute_tanimoto(fp1_json: string, fp2_json: string): string
```

Returns `{"tanimoto": 0.333333}`.

### `butina_cluster`

```typescript
function butina_cluster(conformers_json: string, rmsd_cutoff: number): string
```

```typescript
interface ClusterResult {
  n_clusters: number;
  assignments: number[];
  centroid_indices: number[];
  cluster_sizes: number[];
  rmsd_cutoff: number;   // Å
}
```

### `filter_diverse`

```typescript
function filter_diverse(conformers_json: string, rmsd_cutoff: number): string
```

Returns a JSON array of the retained (cluster centroid) conformer indices.

---

## Descriptors and machine learning

### `compute_ml_descriptors`

```typescript
function compute_ml_descriptors(smiles: string): string
```

```typescript
interface MolecularDescriptors {
  molecular_weight: number;
  n_heavy_atoms: number;
  n_hydrogens: number;
  n_bonds: number;
  n_rotatable_bonds: number;
  n_hbd: number;
  n_hba: number;
  fsp3: number;
  total_abs_charge: number;
  max_charge: number;
  min_charge: number;
  wiener_index: number;
  n_rings: number;
  n_aromatic: number;
  balaban_j: number;
  sum_electronegativity: number;
  sum_polarizability: number;
}
```

> **Deprecated alias:** `compute_molecular_descriptors` → `compute_ml_descriptors`.

### `predict_ml_properties`

```typescript
function predict_ml_properties(smiles: string): string
```

`{"logp", "molar_refractivity", "log_solubility", "lipinski_violations", "lipinski_passes",
"druglikeness"}`.

> **Deprecated alias:** `compute_ml_properties` → `predict_ml_properties`.

### `compute_tpsa`

```typescript
function compute_tpsa(smiles: string): string
```

Topological polar surface area in **Å²** (Ertl fragment contributions). It is a
topological property: no 3D geometry is needed.

### `predict_ensemble`

```typescript
function predict_ensemble(smiles: string): string
```

Ensemble property prediction: the consensus LogP with its uncertainty, the
solubility, the TPSA, the pKa estimates, the Veber rules, the BBB permeability
and an overall confidence score (0–1).

### `compute_whim`

```typescript
function compute_whim(elements_json: string, coords_flat_json: string, weighting: string): string
```

`weighting` is `unit`, `mass`, `volume`, `electronegativity` or `polarizability`. Returns
`{eigenvalues, lambda, nu, gamma, eta, total_spread, anisotropy, directional, kurtosis}`.

### `compute_rdf`

```typescript
function compute_rdf(elements_json: string, coords_flat_json: string, charges_json: string, beta: number, r_max: number, n_points: number): string
```

`charges_json` is a JSON array of partial charges; pass `""` to derive them internally with
Gasteiger. Returns `{radii (Å), rdf_unweighted, rdf_mass, rdf_electronegativity, rdf_charge,
beta (Å⁻²), n_pairs}`.

### `compute_getaway`

```typescript
function compute_getaway(elements_json: string, coords_flat_json: string, bonds_json: string, max_lag: number): string
```

`bonds_json` is a JSON array of `[i, j]` pairs; pass `""` to perceive the bonds from the
geometry. Returns `{leverages, hat_autocorrelation, r_autocorrelation, hat_total, r_total,
h_information, h_max, h_mean, max_lag}`.

### `train_random_forest`

```typescript
function train_random_forest(features_json: string, targets_json: string, n_trees: number, max_depth: number, min_samples_split: number, min_samples_leaf: number, sample_fraction: number, seed: bigint): string
```

`{"model_type": "random_forest", "n_estimators", "n_samples", "predictions", "train_rmse"}`.

### `train_gradient_boosting`

```typescript
function train_gradient_boosting(features_json: string, targets_json: string, n_estimators: number, learning_rate: number, max_depth: number, min_samples_split: number, min_samples_leaf: number, subsample: number, seed: bigint): string
```

`{"model_type": "gradient_boosting", "n_estimators", "n_samples", "predictions", "train_rmse"}`.

---

## Reactivity and topology

### `compute_bond_orders`

```typescript
function compute_bond_orders(elements: string, coords_flat: string): string
```

`{bonds: [{atom_i, atom_j, distance (Å), wiberg, mayer}], num_atoms, wiberg_valence,
mayer_valence}`.

### `compute_frontier_descriptors`

```typescript
function compute_frontier_descriptors(elements: string, coords_flat: string): string
```

`{num_atoms, homo_atom_contributions, lumo_atom_contributions, dual_descriptor, homo_energy
(eV), lumo_energy (eV), gap (eV)}`.

### `compute_fukui_descriptors`

```typescript
function compute_fukui_descriptors(elements: string, coords_flat: string): string
```

`{num_atoms, f_plus, f_minus, f_radical, dual_descriptor, condensed, homo_energy, lumo_energy,
gap, validity_notes}`.

### `compute_reactivity_ranking`

```typescript
function compute_reactivity_ranking(elements: string, coords_flat: string): string
```

`{nucleophilic_attack_sites, electrophilic_attack_sites, radical_attack_sites, notes}` — the
scores are empirical composites, not calibrated kinetic barriers.

### `compute_uv_vis_spectrum`

```typescript
function compute_uv_vis_spectrum(elements: string, coords_flat: string, sigma: number, e_min: number, e_max: number, n_points: number): string
```

`{energies_ev (eV), intensities, peaks, sigma (eV), notes}`.

### `compute_topology`

```typescript
function compute_topology(elements: string, coords_flat: string): string
```

`{metal_centers, warnings}`.

### `analyze_graph_features`

```typescript
function analyze_graph_features(smiles: string): string
```

`{aromaticity: {aromatic_atoms, aromatic_bonds}, stereocenters: {tagged_tetrahedral_centers,
inferred_tetrahedral_centers}}`.

### `analyze_stereo`

```typescript
function analyze_stereo(smiles: string, coords_flat: string): string
```

```typescript
interface StereoAnalysis {
  stereocenters: {
    atom_index: number;
    element: number;
    substituent_indices: number[];
    priorities: number[];
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
  atropisomeric_axes: unknown[];
  prochiral_centers: unknown[];
  helical_chirality: unknown[];
}
```

Pass `"[]"` for the topology-only analysis.

---

## Orbital meshes and grids

### `compute_orbital_mesh`

```typescript
function compute_orbital_mesh(elements_json: string, coords_flat_json: string, method: string, mo_index: number, spacing: number, padding: number, isovalue: number): string
```

`method` is `eht`, `pm3`, `xtb` or `hf3c`; `spacing` and `padding` are in Å.

```typescript
interface OrbitalMeshResult {
  mesh: { vertices: number[]; normals: number[]; indices: number[]; num_triangles: number; isovalue: number };
  grid: { origin: [number, number, number]; spacing: number; dims: [number, number, number]; values: number[] };
  method: string;
  mo_index: number;
  homo_index: number;
  orbital_energies: number[];  // eV
  gap: number;                 // eV
}
```

### `eht_orbital_mesh`

```typescript
function eht_orbital_mesh(elements: string, coords_flat: string, mo_index: number, spacing: number, isovalue: number): string
```

The flat EHT mesh payload: `{vertices, normals, indices, num_triangles, isovalue}`.

### `eht_orbital_grid_typed`

```typescript
function eht_orbital_grid_typed(elements: string, coords_flat: string, mo_index: number, spacing: number): Float32Array
```

The raw volumetric orbital grid as a `Float32Array`. Use `eht_volumetric_grid_info` for the
grid dimensions.

### `eht_orbital_grid_from_coefficients_typed`

```typescript
function eht_orbital_grid_from_coefficients_typed(elements: string, coords_flat: string, coefficients_json: string, mo_index: number, spacing: number): Float32Array
```

The same grid, computed from supplied MO coefficients.

### `eht_volumetric_grid_info`

```typescript
function eht_volumetric_grid_info(elements: string, coords_flat: string, spacing: number, padding: number): string
```

`{"dims": [...], "origin": [...], "spacing": ...}`.

---

## Materials and periodic systems

### `parse_cif` / `write_cif`

```typescript
function parse_cif(cif_text: string): string
function write_cif(cif_json: string): string
```

Read a CIF file into a crystal structure (as JSON) and write it back. The cell
lengths are in **Å**, the cell angles in **degrees** and the atom-site
coordinates are fractional.

### `create_unit_cell`

```typescript
function create_unit_cell(a: number, b: number, c: number, alpha: number, beta: number, gamma: number): string
```

`{"a", "b", "c" (Å), "alpha", "beta", "gamma" (degrees), "volume" (Å³)}`.

### `assemble_framework`

```typescript
function assemble_framework(topology: string, metal: number, geometry: string, lattice_a: number, supercell: number): string
```

`topology` is `pcu`, `dia` or `sql`; `geometry` is `linear`, `trigonal`, `tetrahedral`,
`octahedral`, `square_planar` or `trigonal`. Returns
`{cell: {lattice}, atoms: [{element, frac_coords}], labels}`.

### `space_group_by_number` / `space_group_by_symbol`

```typescript
function space_group_by_number(number: number): string
function space_group_by_symbol(symbol: string): string
```

```typescript
interface SpaceGroup {
  number: number;
  symbol: string;
  crystal_system: string;
  lattice_type: string;
  point_group: string;
  operations: { rotation: number[][]; translation: number[] }[];
}
```

### `all_space_groups`

```typescript
function all_space_groups(): string
```

Returns a JSON array of `[number, symbol]` pairs — all 230 groups.

### `build_periodic_molecule`

```typescript
function build_periodic_molecule(elements_json: string, frac_coords_json: string, cell_json: string, bond_tolerance: number): string
```

`frac_coords_json` is a JSON array of `[x, y, z]` triples and `cell_json` the JSON array
`[a, b, c, alpha, beta, gamma]` (Å, degrees). `bond_tolerance` is extra bonding tolerance in
Å (`0.0` = engine default). Returns
`{elements, frac_coords, cart_coords, bonds, cell, n_atoms}`.

### `detect_hapticity`

```typescript
function detect_hapticity(elements_json: string, coords_flat_json: string, bonds_json: string): string
```

`bonds_json` is a JSON array of `[i, j]` pairs; pass `""` to perceive the bonds from the
geometry. Returns `{interactions, is_metallocene, n_interactions}`.

### `optimize_framework`

```typescript
function optimize_framework(smiles: string, coords_flat_json: string, method: string, max_iter: number, force_tol: number, energy_tol: number, max_step: number, fixed_atoms_json: string, cell_json: string): string
```

`cell_json` optionally enables periodic (minimum-image) optimisation and is the JSON array
`[a, b, c, alpha, beta, gamma]`. Energies are reported in **eV** and forces in **eV/Å**.

```typescript
interface FrameworkOptResult {
  positions: [number, number, number][];  // Å
  energy: number;                        // eV
  forces: [number, number, number][];    // eV/Å
  max_force: number;                     // eV/Å
  n_iterations: number;
  converged: boolean;
  energy_history: number[];              // eV
}
```

---

## Dynamics, NEB and reaction paths

### `compute_md_trajectory`

```typescript
function compute_md_trajectory(smiles: string, coords_json: string, n_steps: number, dt_fs: number, seed: number): string
```

Short NVE molecular dynamics with the UFF force field. `dt_fs` in fs; energies in **kcal/mol**.

```typescript
interface MdTrajectory {
  frames: { step: number; time_fs: number; coords: number[]; potential_energy_kcal_mol: number; kinetic_energy_kcal_mol: number; temperature_k: number }[];
  dt_fs: number;
  notes: string[];
  energy_drift_percent: number | null;
}
```

### `compute_md_trajectory_seed64`

```typescript
function compute_md_trajectory_seed64(smiles: string, coords_json: string, n_steps: number, dt_fs: number, seed: bigint): string
```

### `compute_md_trajectory_nvt`

```typescript
function compute_md_trajectory_nvt(smiles: string, coords_json: string, n_steps: number, dt_fs: number, seed: number, target_temp_k: number, thermostat_tau_fs: number): string
```

NVT with the Berendsen thermostat; `target_temp_k` in K, `thermostat_tau_fs` in fs.

### `compute_md_trajectory_nvt_seed64`

```typescript
function compute_md_trajectory_nvt_seed64(smiles: string, coords_json: string, n_steps: number, dt_fs: number, seed: bigint, target_temp_k: number, thermostat_tau_fs: number): string
```

### `compute_md_trajectory_with_method`

```typescript
function compute_md_trajectory_with_method(smiles: string, coords_json: string, elements_json: string, n_steps: number, dt_fs: number, seed: number, target_temp_k: number, backend: string): string
```

`backend` is `uff`, `pm3` or `xtb`; `target_temp_k > 0` enables the Berendsen thermostat with a
100 fs coupling.

### `compute_md_trajectory_with_method_seed64`

```typescript
function compute_md_trajectory_with_method_seed64(smiles: string, coords_json: string, elements_json: string, n_steps: number, dt_fs: number, seed: bigint, target_temp_k: number, backend: string): string
```

### `compute_simplified_neb_path`

```typescript
function compute_simplified_neb_path(smiles: string, start_coords_json: string, end_coords_json: string, n_images: number, n_iter: number, spring_k: number, step_size: number): string
```

`{"images": [...], "notes": [...]}` — a low-cost exploratory path, not a full
climbing-image / tangent-projected NEB.

### `compute_neb_path_with_method`

```typescript
function compute_neb_path_with_method(smiles: string, start_coords_json: string, end_coords_json: string, n_images: number, n_iter: number, spring_k: number, step_size: number, method: string): string
```

`method` is `uff`, `mmff94`, `pm3`, `xtb`, `gfn1`, `gfn2` or `hf3c`.

### `neb_single_point_energy` / `neb_energy_and_gradient`

```typescript
function neb_single_point_energy(method: string, smiles: string, coords_json: string): string
function neb_energy_and_gradient(method: string, smiles: string, coords_json: string): string
```

`{"energy": <f64> (kcal/mol), "method": "<str>"}` and
`{"energy": <f64> (kcal/mol), "gradient": [<f64>...], "method": "<str>"}`.

### `list_neb_backends`

```typescript
function list_neb_backends(elements_json: string): string
```

Returns a JSON array of `{"method": string, "available": bool, "note": string}`.

### `search_conformers_with_uff`

```typescript
function search_conformers_with_uff(smiles: string, n_samples: number, seed: number, rmsd_threshold: number): string
```

`{"generated", "unique", "rotatable_bonds", "conformers", "clusters", "notes"}` — energies in
kcal/mol, RMSD threshold in Å.

### `search_conformers_with_uff_seed64`

```typescript
function search_conformers_with_uff_seed64(smiles: string, n_samples: number, seed: bigint, rmsd_threshold: number): string
```

### `compute_reaction_dynamics`

```typescript
function compute_reaction_dynamics(reactant_smiles_json: string, product_smiles_json: string, method: string, config_json: string): string
```

Embeds the reactant and product fragments, builds oriented complexes, runs NEB for the reactive
zone and generates approach/departure frames. `method` is `uff`, `mmff94`, `pm3`, `xtb`, `gfn1`,
`gfn2` or `hf3c`; `config_json` is an optional JSON configuration (pass `""` or `"{}"` for the
defaults).

---

## SMIRKS

### `parse_smirks`

```typescript
function parse_smirks(smirks: string): string
```

```typescript
interface SmirksTransform {
  reactant_smarts: string[];
  product_smarts: string[];
  atom_map: Record<string, number>;
  bond_changes: unknown[];
  smirks: string;
}
```

### `apply_smirks`

```typescript
function apply_smirks(smirks: string, smiles: string): string
```

```typescript
interface SmirksResult {
  products: string[];
  atom_mapping: Record<string, number>;
  n_transforms: number;
  success: boolean;
  messages: string[];
}
```

---

## Transport and worker helpers

### `pack_conformers` / `pack_batch_arrow`

```typescript
function pack_conformers(results_json: string): string
function pack_batch_arrow(results_json: string): string
```

Pack a JSON array of conformer results into sci-form's columnar (Arrow-like) format. This is
**not** the Apache Arrow IPC format — it is a JSON serialisation of sci-form's own columnar
layout.

```typescript
interface RecordBatch {
  num_rows: number;
  schema: { name: string; dtype: string }[];
  float_columns: Record<string, number[]>;
  int_columns: Record<string, number[]>;
  uint8_columns: Record<string, number[]>;
}
```

### `split_worker_tasks`

```typescript
function split_worker_tasks(smiles_json: string, n_workers: number, seed: number): string
```

Returns a JSON array of worker tasks: `{id: number, kind: string, smiles: string[],
params: unknown[]}`.

### `split_worker_tasks_seed64`

```typescript
function split_worker_tasks_seed64(smiles_json: string, n_workers: number, seed: bigint): string
```

### `estimate_workers`

```typescript
function estimate_workers(n_items: number, max_workers: number): number
```

Heuristic to pick the optimal worker count.

---

## System planning

### `system_capabilities` / `get_system_capabilities`

```typescript
function system_capabilities(elements: string): string
function get_system_capabilities(elements: string): string
```

```typescript
interface MethodCapability {
  available: boolean;
  confidence: string;   // "supported" | "experimental" | "unsupported"
  unsupported_elements: number[];
  warnings: string[];
}

interface SystemCapabilities {
  embed: MethodCapability;
  uff: MethodCapability;
  eht: MethodCapability;
  population: MethodCapability;
  orbital_grid: MethodCapability;
}
```

### `system_method_plan` / `get_system_method_plan`

```typescript
function system_method_plan(elements: string): string
function get_system_method_plan(elements: string): string
```

```typescript
interface PropertyMethodPlan {
  property: string;
  recommended: string | null;
  fallback: string | null;
  rationale: string[];
  methods: MethodMetadata[];
}

interface SystemMethodPlan {
  capabilities: SystemCapabilities;
  geometry: PropertyMethodPlan;
  force_field_energy: PropertyMethodPlan;
  orbitals: PropertyMethodPlan;
  population: PropertyMethodPlan;
  orbital_grid: PropertyMethodPlan;
}
```

### `compare_methods`

```typescript
function compare_methods(smiles: string, elements: string, coords_flat: string, allow_experimental_eht: boolean): string
```

`{"plan": SystemMethodPlan, "comparisons": MethodComparisonEntry[]}`.

---

## Experimental (feature-gated)

These exports only exist when the corresponding Cargo feature is enabled on the
`sci-form-wasm` package. The always-available experimental entry points are:

| Function | Feature |
|---|---|
| `compute_eeq_charges`, `compute_eeq_energy` | always built |
| `compute_alpb_solvation`, `compute_alpb_born_radii` | `experimental-alpb` |
| `compute_d4_energy` | `experimental-d4` |
| `compute_cpm_charges`, `compute_cpm_surface` | `experimental-cpm` |
| `compute_kmesh`, `compute_periodic_dos`, `solve_periodic_randnla`, `bloch_phase`, `build_bloch_hamiltonian`, `assemble_periodic_operators`, `validate_electron_count`, `validate_kmesh_weights` | `alpha-periodic-linear` |
| `compute_htst_rate`, `compute_htst_sweep`, `extract_kinetics_diagnostics`, `solve_microkinetic_network`, `solve_microkinetic_steady_state`, `analyze_gsm_mbh_htst_step` | `alpha-kinetics` |
| `chart_to_json`, `chart_from_json`, `validate_chart_schema`, `pack_chart_payload`, `compute_edl_chart`, `compute_capacitance_chart`, `compute_arrhenius_chart`, `compute_thermal_prop_chart`, `compute_phase_portrait_chart`, `compute_band_structure_chart`, `compute_dos_chart`, `compute_kpoint_path_chart`, `compute_fermi_surface_chart`, `compute_trajectory_chart`, `compute_reaction_coordinate_chart` | `alpha-render-bridge` |
| `alpha_*`, `beta_*` research functions | the matching `alpha-*` / `beta-*` features |

---

## Usage patterns

### High-throughput conformer generation

```typescript
import init, { embed_coords_typed, estimate_workers, split_worker_tasks } from 'sci-form-wasm';

await init();

const coords = embed_coords_typed('c1ccccc1', 42);
const nWorkers = estimate_workers(1000, navigator.hardwareConcurrency ?? 4);
const tasks = JSON.parse(split_worker_tasks(JSON.stringify(smilesList), nWorkers, 42));
```

### Three.js molecule viewer

```typescript
import init, { embed, solve_eht, eht_orbital_mesh } from 'sci-form-wasm';
import * as THREE from 'three';

await init();

const result = JSON.parse(embed('c1ccccc1', 42));
const elements = JSON.stringify(result.elements);
const coordsJson = JSON.stringify(result.coords);

const scene = new THREE.Scene();
for (let i = 0; i < result.num_atoms; i++) {
  const sphere = new THREE.Mesh(
    new THREE.SphereGeometry(0.3),
    new THREE.MeshPhongMaterial({ color: 0x336699 })
  );
  sphere.position.set(result.coords[i*3], result.coords[i*3+1], result.coords[i*3+2]);
  scene.add(sphere);
}

const eht = JSON.parse(solve_eht(elements, coordsJson, 1.75));
const mesh = JSON.parse(eht_orbital_mesh(elements, coordsJson, eht.homo_index, 0.2, 0.02));

const orbGeom = new THREE.BufferGeometry();
orbGeom.setAttribute('position', new THREE.Float32BufferAttribute(mesh.vertices, 3));
orbGeom.setAttribute('normal', new THREE.Float32BufferAttribute(mesh.normals, 3));
orbGeom.setIndex(mesh.indices);
scene.add(new THREE.Mesh(orbGeom, new THREE.MeshPhongMaterial({
  color: 0x0077ff, transparent: true, opacity: 0.5, side: THREE.DoubleSide
})));
```

### React hook

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

---

## See also

- [CLI reference](./cli.md)
- [Python reference](./python.md)
- [Rust reference](./rust.md)
- [Complete supported reference](./complete-supported-reference.md)
- [TypeScript guide](../guide/typescript.md)
