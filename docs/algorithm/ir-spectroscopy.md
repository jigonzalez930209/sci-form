# IR Spectroscopy — Numerical Hessian

sci-form computes infrared spectra from first principles via a numerical second-derivative (Hessian) approach. The resulting normal-mode frequencies and dipole derivatives directly yield IR intensities.

## Pipeline

<SvgDiagram src="/svg/ir-hessian-pipeline.svg" alt="IR Hessian pipeline: Geometry → 6N Energy Evaluations → Hessian → Mass-weighted → Diagonalize → Modes + Frequencies → Dipole Derivatives → IR spectrum" />

## Theory

### Numerical Hessian via Central Finite Differences

For a molecule with $N$ atoms, the Cartesian Hessian $H_{ij}$ is computed using central finite differences with step size $\delta$:

$$
H_{ij} = \frac{\partial^2 E}{\partial q_i \partial q_j} \approx \frac{E(\mathbf{q} + \delta\hat{e}_i + \delta\hat{e}_j) - E(\mathbf{q} + \delta\hat{e}_i - \delta\hat{e}_j) - E(\mathbf{q} - \delta\hat{e}_i + \delta\hat{e}_j) + E(\mathbf{q} - \delta\hat{e}_i - \delta\hat{e}_j)}{4\delta^2}
$$

where $q_i$ runs over the $3N$ Cartesian degrees of freedom and $\hat{e}_i$ are unit displacement vectors. This requires $2 \times (3N)^2$ energy evaluations for the full Hessian, reduced to $6N$ evaluations using the mixed-partial shortcut:

$$
H_{ij} = \frac{E(\mathbf{q} + \delta_i + \delta_j) - E(\mathbf{q} + \delta_i - \delta_j) - E(\mathbf{q} - \delta_i + \delta_j) + E(\mathbf{q} - \delta_i - \delta_j)}{4\delta^2}
$$

Default step size: $\delta = 0.01$ Å.

### Mass-Weighted Hessian

To account for atomic masses in Newton's equations of motion, the Hessian is transformed to mass-weighted coordinates:

$$
\tilde{H}_{ij} = \frac{H_{ij}}{\sqrt{m_i \, m_j}}
$$

where $m_i$ is the mass of the atom associated with degree of freedom $i$.

### Normal Mode Diagonalization

The mass-weighted Hessian is diagonalized to obtain eigenvalues $\lambda_k$ (in units of eV/(amu·Å²)) and eigenvectors $\mathbf{L}$ (normal coordinates):

$$
\tilde{H} \mathbf{L} = \boldsymbol{\Lambda} \mathbf{L}
$$

Vibrational frequencies in cm⁻¹ are obtained from the eigenvalues:

$$
\tilde{\nu}_k = \frac{1}{2\pi c} \sqrt{\frac{F \lambda_k}{1 \, \text{amu} \cdot \text{Å}^2}} = \frac{\sqrt{9.6485 \times 10^{27} \, \lambda_k}}{2\pi c}
$$

where $c = 2.998 \times 10^{10}$ cm/s. Negative eigenvalues (imaginary frequencies) indicate structural instabilities; they are reported with a negative sign by convention.

### Translation/Rotation Separation

The $3N$ modes include 3 translations and 3 (or 2 for linear molecules) rotations. These appear as near-zero frequency modes (|ν| < 10 cm⁻¹) and are excluded from the IR spectrum. The number of real vibrational modes is:

- $3N - 6$ for non-linear molecules
- $3N - 5$ for linear molecules

### IR Intensities — Dipole Derivatives

The IR intensity of mode $k$ is proportional to the squared gradient of the dipole moment along the normal coordinate:

$$
I_k \propto \left|\frac{\partial \boldsymbol{\mu}}{\partial Q_k}\right|^2
$$

Numerically, this is estimated via finite differences of the electronic dipole along each displacement, then projected onto the normal mode:

$$
\frac{\partial \mu_\alpha}{\partial Q_k} \approx \sum_i L_{ik} \frac{\mu_\alpha(\mathbf{q} + \delta_i) - \mu_\alpha(\mathbf{q} - \delta_i)}{2\delta}
$$

Intensities are reported in km/mol (standard IUPAC unit: 1 km/mol = 10 L mol⁻¹ cm⁻²).

### Zero-Point Vibrational Energy (ZPVE)

The ZPVE is summed over real modes:

$$
E_\text{ZPVE} = \frac{1}{2} \sum_k h c \tilde{\nu}_k = \frac{1}{2} \sum_k \hbar \omega_k
$$

In eV, using 1 cm⁻¹ = 1.2399 × 10⁻⁴ eV:

$$
E_\text{ZPVE}\,[\text{eV}] = \frac{1}{2} \sum_k \tilde{\nu}_k \times 1.2399 \times 10^{-4}
$$

### Lorentzian Broadened Spectrum

The discrete IR peaks are broadened using Lorentzian line shapes:

$$
A(\tilde{\nu}) = \sum_k I_k \cdot \frac{\gamma/\pi}{(\tilde{\nu} - \tilde{\nu}_k)^2 + \gamma^2}
$$

where $\gamma$ is the half-width at half-maximum (HWHM) in cm⁻¹.

## Parameters

### `compute_vibrational_analysis`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `elements` | `&[u8]` | — | Atomic numbers |
| `positions` | `&[[f64;3]]` | — | Cartesian coordinates (Å) |
| `method` | `&str` | `"eht"` | Semiempirical method: `"eht"`, `"pm3"`, or `"xtb"` |
| `step_size` | `Option<f64>` | `None` (→ 0.01 Å) | Finite-difference step Δ in Å |

### `compute_ir_spectrum`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `analysis` | `&VibrationalAnalysis` | — | Result from `compute_vibrational_analysis` |
| `gamma` | `f64` | — | Lorentzian HWHM in cm⁻¹ (typical: 10–30) |
| `wn_min` | `f64` | — | Low-wavenumber limit (cm⁻¹) |
| `wn_max` | `f64` | — | High-wavenumber limit (cm⁻¹) |
| `n_points` | `usize` | — | Grid resolution |

## Output Types

### `VibrationalAnalysis`

| Field | Type | Description |
|-------|------|-------------|
| `n_atoms` | `usize` | Number of atoms |
| `modes` | `Vec<VibrationalMode>` | All $3N$ modes sorted by frequency |
| `n_real_modes` | `usize` | Real modes (excluding translation/rotation) |
| `zpve_ev` | `f64` | Zero-point vibrational energy (eV) |
| `method` | `String` | Method used for Hessian |
| `notes` | `Vec<String>` | Caveats |

### `VibrationalMode`

| Field | Type | Description |
|-------|------|-------------|
| `frequency_cm1` | `f64` | Frequency (cm⁻¹); negative = imaginary |
| `ir_intensity` | `f64` | IR intensity (km/mol) |
| `displacement` | `Vec<f64>` | $3N$ displacement vector |
| `is_real` | `bool` | True if real vibrational mode |

### `IrSpectrum`

| Field | Type | Description |
|-------|------|-------------|
| `wavenumbers` | `Vec<f64>` | Wavenumber grid (cm⁻¹) |
| `intensities` | `Vec<f64>` | Broadened absorption (km/mol) |
| `peaks` | `Vec<IrPeak>` | Discrete peaks |
| `gamma` | `f64` | Broadening HWHM (cm⁻¹) |
| `notes` | `Vec<String>` | Caveats |

### `IrPeak`

| Field | Type | Description |
|-------|------|-------------|
| `frequency_cm1` | `f64` | Peak centre (cm⁻¹) |
| `ir_intensity` | `f64` | Raw intensity (km/mol) |
| `mode_index` | `usize` | Index into `modes` |

## API

### Rust

```rust
use sci_form::{embed, compute_vibrational_analysis, compute_ir_spectrum};

let conf = embed("CCO", 42);
let pos: Vec<[f64;3]> = conf.coords.chunks(3).map(|c| [c[0],c[1],c[2]]).collect();

// ~6N energy evaluations (may take a few seconds for large molecules)
let vib = compute_vibrational_analysis(&conf.elements, &pos, "xtb", None).unwrap();
println!("ZPVE: {:.4} eV, {} real modes", vib.zpve_ev, vib.n_real_modes);

// Dominant IR mode
if let Some(peak) = vib.modes.iter()
    .filter(|m| m.is_real)
    .max_by(|a,b| a.ir_intensity.partial_cmp(&b.ir_intensity).unwrap())
{
    println!("Strongest band: {:.1} cm⁻¹  ({:.1} km/mol)", peak.frequency_cm1, peak.ir_intensity);
}

// Broaden into a 600–4000 cm⁻¹ spectrum
let spec = compute_ir_spectrum(&vib, 15.0, 600.0, 4000.0, 1000);
println!("{} wavenumber points, γ = {} cm⁻¹", spec.wavenumbers.len(), spec.gamma);
```

### Python

```python
from sci_form import embed, vibrational_analysis, ir_spectrum

conf = embed("CCO", seed=42)
vib  = vibrational_analysis(conf.elements, conf.coords, method="xtb")
spec = ir_spectrum(vib, gamma=15.0, wn_min=600.0, wn_max=4000.0, n_points=1000)

print(f"ZPVE: {vib.zpve_ev:.4f} eV, {vib.n_real_modes} real modes")

import matplotlib.pyplot as plt
plt.plot(spec.wavenumbers, spec.intensities)
plt.xlabel("Wavenumber (cm⁻¹)")
plt.ylabel("Absorbance (km/mol)")
plt.title("IR Spectrum of Ethanol (GFN-xTB)")
plt.gca().invert_xaxis()
plt.show()
```

### TypeScript/WASM

```typescript
import init, { embed, compute_vibrational_analysis, compute_ir_spectrum } from 'sci-form-wasm';
await init();

const r = JSON.parse(embed('CCO', 42));
const elements = JSON.stringify(r.elements);
const coords   = JSON.stringify(r.coords);

const vib = JSON.parse(compute_vibrational_analysis(elements, coords, 'xtb', null));
const spec = JSON.parse(compute_ir_spectrum(JSON.stringify(vib), 15.0, 600.0, 4000.0, 1000));

console.log(`ZPVE: ${vib.zpve_ev.toFixed(4)} eV`);
const maxI = Math.max(...spec.intensities);
const wn   = spec.wavenumbers[spec.intensities.indexOf(maxI)];
console.log(`Dominant band at ${wn.toFixed(1)} cm⁻¹`);
```

## Performance Notes

- The Hessian requires $6N$ energy evaluations for central differences, where $N$ is the number of atoms.
- For EHT: ~0.3 ms/eval → ~1 s for 50 atoms.
- For PM3: ~5 ms/eval → ~15 s for 50 atoms.
- For GFN-xTB: ~20 ms/eval → ~60 s for 50 atoms.
- Use `method = "eht"` for rapid screening; `"xtb"` for higher fidelity.
