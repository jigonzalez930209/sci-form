# NMR Spectroscopy — HOSE Codes & Karplus Equation

sci-form predicts ¹H and ¹³C NMR chemical shifts from the molecular graph using HOSE codes (Hierarchical Organization of Spherical Environments), and estimates vicinal J-coupling constants via the Karplus equation.

## Pipeline

<SvgDiagram src="/svg/nmr-pipeline.svg" alt="NMR pipeline: SMILES → Mol.Graph → HOSE codes → Shift lookup → J-Coupling → Peak splitting → Lorentzian → NMR Spectrum" />

## Theory

### HOSE Codes — Chemical Environment Fingerprint

A **HOSE code** encodes the spherical chemical environment around each atom by performing a breadth-first traversal of the molecular graph out to radius $R$ (default 4 bonds):

1. Start at the central atom (sphere 0).
2. At each sphere $r = 1, 2, \ldots, R$, record the neighbor atoms sorted by atomic number, bond order, and connectivity count.
3. Concatenate these strings separated by sphere delimiters `(` and `)`.

Example for the **carbonyl carbon** in acetone (CH₃COCH₃):

```
C-6;C;=O,C/C,C//H,H,H,H,H,H/
```

Reading: central C → sphere 1 has C=O and C–C → sphere 2 has C and C → sphere 3 has H×3 and H×3.

The resulting string is matched against a reference database of >10,000 shifts to find similar chemical environments and predict the shift by analogy.

### Chemical Shift Prediction

The predicted shift for atom $i$ is a combination of:

1. **Hybridization base value** — sp³, sp², sp (alkyne), and aromatic carbons/hydrogens have distinct chemical shift ranges.
2. **Electronegativity correction** — adjacent electronegative atoms (O, N, F, halogens) deshield via inductive effects:

$$
\delta_\text{ind} = \sum_{j \in \text{neighbors}} k_\alpha (\chi_j - \chi_\text{ref})
$$

3. **Ring current effect** — aromatic rings cause through-space shielding/deshielding:

$$
\delta_\text{ring} = \begin{cases} +1.5 \; \text{ppm} & \text{H} \text{ in aromatic ring} \\ -0.5 \; \text{ppm} & \text{H} \text{ above aromatic ring (anisotropy)} \end{cases}
$$

4. **Functional group corrections** — carbonyl (C=O), carboxylic acid, aldehyde, ether, amine corrections applied by pattern matching.

#### Reference Ranges

| Nucleus | Environment | δ range (ppm) |
|---------|-------------|----------------|
| ¹H | sp³ C–H (methyl) | 0.8–1.0 |
| ¹H | sp³ C–H (CH₂, CH) | 1.2–2.0 |
| ¹H | sp² C–H (alkene) | 4.8–6.5 |
| ¹H | aromatic C–H | 6.5–8.5 |
| ¹H | aldehyde C–H | 9.0–10.5 |
| ¹H | O–H (alcohol) | 1.0–5.0 |
| ¹H | O–H (carboxylic acid) | 10–13 |
| ¹³C | sp³ (alkyl) | 10–40 |
| ¹³C | sp³ with O/N | 40–90 |
| ¹³C | sp² alkene/aromatic | 100–160 |
| ¹³C | carbonyl C=O (ketone) | 195–215 |
| ¹³C | carboxylic COOH | 160–185 |

### Karplus Equation — ³J(H,H) Vicinal Coupling

Vicinal (three-bond) ¹H–¹H coupling constants are predicted using the Karplus equation:

$$
{}^3J(\phi) = A \cos^2\phi - B \cos\phi + C
$$

with parameters $A = 7.76$, $B = 1.10$, $C = 1.40$ Hz (Altona–Sundaralingam, JACS 1972).

<SvgDiagram src="/svg/nmr-karplus.svg" alt="Karplus curve: J(φ) vs dihedral angle, maximum near 0° and 180°, minimum near 90°" />

When 3D coordinates are available, the actual H–C–C–H dihedral angle $\phi$ is computed from the positions. When only topology is available, a **free-rotation average** is used:

$$
\langle {}^3J \rangle = \int_0^{2\pi} (A\cos^2\phi - B\cos\phi + C) \frac{d\phi}{2\pi} = \frac{A}{2} + C \approx 5.3 \; \text{Hz}
$$

#### Geminal ²J and Long-Range ⁴J Couplings

- **²J(H,H)**: $-10$ to $+15$ Hz depending on hybridization and substituents; sci-form uses topology-based estimates (sp³ ~$-12$ Hz, sp² ~$+3$ Hz).
- **⁴J(H,H)**: typically < 3 Hz; included for zig-zag (W-shaped) paths.

### Spectrum Generation — Pascal's Triangle Splitting

For ¹H nuclei with adjacent vicinal couplings $J_1, J_2, \ldots$, the multiplet pattern is built by iteratively convolving the peak position with first-order splittings:

$$
\text{doublet}(\delta, J) : \left\{\delta - \frac{J}{2}, \delta + \frac{J}{2}\right\} \quad \text{relative heights} \; 1:1
$$

$$
\text{triplet}(\delta, J) : \left\{\delta - J, \delta, \delta + J\right\} \quad \text{relative heights} \; 1:2:1
$$

$$
\text{quartet}(\delta, J) : \left\{\delta - \frac{3J}{2}, \delta - \frac{J}{2}, \delta + \frac{J}{2}, \delta + \frac{3J}{2}\right\} \quad \text{relative heights} \; 1:3:3:1
$$

### Lorentzian Broadening

Each multiplet line is broadened by a Lorentzian:

$$
S(\delta) = \sum_k w_k \cdot \frac{\gamma/\pi}{(\delta - \delta_k)^2 + \gamma^2}
$$

where $w_k$ are Pascal's triangle weights and $\gamma$ is the half-width (default 0.01 ppm for ¹H, 0.5 ppm for ¹³C).

## Parameters

### `predict_nmr_shifts`

| Parameter | Type | Description |
|-----------|------|-------------|
| `smiles` | `&str` | SMILES string |

### `predict_nmr_couplings`

| Parameter | Type | Description |
|-----------|------|-------------|
| `smiles` | `&str` | SMILES string |
| `positions` | `&[[f64;3]]` | 3D coordinates (Å); pass `&[]` for topological estimate |

### `compute_nmr_spectrum`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `smiles` | `&str` | — | SMILES string |
| `nucleus` | `&str` | — | `"1H"` or `"13C"` |
| `gamma` | `f64` | `0.01` | Lorentzian HWHM in ppm |
| `ppm_min` | `f64` | `-2.0` (`"1H"`) | Spectral window start |
| `ppm_max` | `f64` | `14.0` (`"1H"`) | Spectral window end |
| `n_points` | `usize` | `2000` | Grid resolution |

## Output Types

### `NmrShiftResult`

| Field | Type | Description |
|-------|------|-------------|
| `h_shifts` | `Vec<ChemicalShift>` | Predicted ¹H shifts |
| `c_shifts` | `Vec<ChemicalShift>` | Predicted ¹³C shifts |
| `notes` | `Vec<String>` | Caveats |

### `ChemicalShift`

| Field | Type | Description |
|-------|------|-------------|
| `atom_index` | `usize` | Atom index in molecule |
| `element` | `u8` | Atomic number (1 or 6) |
| `shift_ppm` | `f64` | Predicted chemical shift (ppm) |
| `environment` | `String` | Environment classification |
| `confidence` | `f64` | Confidence 0.0–1.0 |

### `JCoupling`

| Field | Type | Description |
|-------|------|-------------|
| `h1_index` | `usize` | First H atom index |
| `h2_index` | `usize` | Second H atom index |
| `j_hz` | `f64` | Coupling constant (Hz) |
| `n_bonds` | `usize` | Number of bonds ($^2$J, $^3$J, $^4$J) |
| `coupling_type` | `String` | Type name (e.g., `"vicinal_3J"`) |

### `NmrSpectrum`

| Field | Type | Description |
|-------|------|-------------|
| `ppm_axis` | `Vec<f64>` | Chemical shift grid (ppm) |
| `intensities` | `Vec<f64>` | Broadened spectrum intensities |
| `peaks` | `Vec<NmrPeak>` | Discrete multiplet lines |
| `nucleus` | `String` | `"1H"` or `"13C"` |
| `gamma` | `f64` | Broadening HWHM (ppm) |

### `NmrPeak`

| Field | Type | Description |
|-------|------|-------------|
| `shift_ppm` | `f64` | Peak position (ppm) |
| `intensity` | `f64` | Relative intensity |
| `atom_index` | `usize` | Source atom index |
| `multiplicity` | `String` | e.g., `"s"`, `"d"`, `"t"`, `"q"`, `"m"` |

## API

### Rust

```rust
use sci_form::{embed, predict_nmr_shifts, predict_nmr_couplings, compute_nmr_spectrum};

// Shifts only (no 3D needed)
let shifts = predict_nmr_shifts("CCO").unwrap();
for h in &shifts.h_shifts {
    println!("H#{}: {:.2} ppm  ({})", h.atom_index, h.shift_ppm, h.environment);
}

// With 3D coordinates for accurate Karplus coupling
let conf = embed("CC", 42);
let pos: Vec<[f64;3]> = conf.coords.chunks(3).map(|c| [c[0],c[1],c[2]]).collect();
let couplings = predict_nmr_couplings("CC", &pos).unwrap();
for j in &couplings {
    println!("{}J(H{},H{}): {:.2} Hz  ({})",
        j.n_bonds, j.h1_index, j.h2_index, j.j_hz, j.coupling_type);
}

// Full ¹H spectrum
let spec = compute_nmr_spectrum("CCO", "1H", 0.01, -2.0, 12.0, 2000).unwrap();
println!("{} ¹H peaks", spec.peaks.len());
```

### Python

```python
from sci_form import nmr_shifts, nmr_couplings, nmr_spectrum, embed

# Chemical shifts
shifts = nmr_shifts("CCO")
for h in shifts.h_shifts:
    print(f"H#{h.atom_index}: {h.shift_ppm:.2f} ppm  [{h.environment}]")
for c in shifts.c_shifts:
    print(f"C#{c.atom_index}: {c.shift_ppm:.2f} ppm  [{c.environment}]")

# Optional: accurate J with 3D
conf = embed("CC", seed=42)
couplings = nmr_couplings("CC", conf.coords)
for j in couplings:
    print(f"{j.n_bonds}J(H{j.h1_index},H{j.h2_index}) = {j.j_hz:.2f} Hz")

# Full ¹H spectrum
spec = nmr_spectrum("CCO", nucleus="1H", gamma=0.01, ppm_min=-2.0, ppm_max=12.0)

import matplotlib.pyplot as plt
plt.plot(spec.ppm_axis, spec.intensities)
plt.xlabel("δ (ppm)")
plt.gca().invert_xaxis()
plt.show()
```

### TypeScript/WASM

```typescript
import init, { predict_nmr_shifts, compute_nmr_spectrum } from 'sci-form-wasm';
await init();

const shifts = JSON.parse(predict_nmr_shifts('CCO'));
console.log('1H shifts:', shifts.h_shifts.map((h: any) =>
  `H#${h.atom_index}: ${h.shift_ppm.toFixed(2)} ppm`).join(', '));

const spec = JSON.parse(compute_nmr_spectrum('CCO', '1H', 0.01, -2.0, 12.0, 2000));
const maxI = Math.max(...spec.intensities);
const ppm  = spec.ppm_axis[spec.intensities.indexOf(maxI)];
console.log(`Most intense peak at δ ${ppm.toFixed(2)} ppm`);
```

## Limitations and Caveats

- Predictions are empirical; typical accuracy is ±0.3 ppm for ¹H and ±5 ppm for ¹³C.
- Stereochemical effects (axial/equatorial, diastereotopic protons) are partially handled through 3D Karplus.
- Solvent and temperature effects are not modeled.
- ¹⁵N, ³¹P, ¹⁹F, and other heteronuclear shifts are not supported.
- For HOSE code matching, rare environments with no database match fall back to hybridization-based defaults.
