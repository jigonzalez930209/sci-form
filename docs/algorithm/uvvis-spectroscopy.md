# UV-Vis Spectroscopy (sTDA)

sci-form implements a simplified Tamm-Dancoff approximation (sTDA) for UV-Vis absorption spectra using EHT molecular orbital transitions with oscillator strength weighting and configurable broadening.

## Overview

<SvgDiagram src="/svg/spectroscopy-overview.svg" alt="Spectroscopy overview — three tracks D1 UV-Vis, D2 IR, D3 NMR" />

## sTDA Pipeline

<SvgDiagram src="/svg/uvvis-stda-pipeline.svg" alt="sTDA UV-Vis pipeline: Geometry → EHT/xTB → MO Transitions → Oscillator Strength → Broadening → ε(λ)" />

## Theory

### Molecular Orbital Transitions

Starting from a converged EHT calculation with $N_\text{occ}$ occupied and $N_\text{virt}$ virtual orbitals, all single-excitation transitions $i \to a$ are enumerated:

$$
\Delta E_{ia} = \varepsilon_a - \varepsilon_i \quad (\text{HOMO-LUMO and all sub-gap excitations})
$$

The energy is converted to wavelength:

$$
\lambda_{ia} = \frac{hc}{\Delta E_{ia}} = \frac{1239.84}{\Delta E_{ia}[\text{eV}]} \; \text{nm}
$$

### Oscillator Strength

The oscillator strength $f_{ia}$ accounts for the transition probability via the electric dipole matrix element:

$$
f_{ia} = \frac{2 m_e \Delta E_{ia}}{3 \hbar^2} |\langle \psi_i | \mathbf{r} | \psi_a \rangle|^2
$$

In the sTDA approximation, the transition dipole moment is evaluated as a sum over orbital coefficients:

$$
|\mu_{ia}|^2 = \left| \sum_{\mu\nu} C_{\mu i} \, C_{\nu a} \, \langle \phi_\mu | \mathbf{r} | \phi_\nu \rangle \right|^2
$$

### Broadened Absorption Spectrum

The molar extinction coefficient $\varepsilon(\lambda)$ is built by summing Gaussian or Lorentzian line shapes centred at each transition:

**Gaussian broadening** (default):

$$
\varepsilon(\lambda) = \sum_{ia} f_{ia} \cdot G(\lambda - \lambda_{ia}, \sigma)
$$

$$
G(x, \sigma) = \frac{1}{\sigma\sqrt{2\pi}} \exp\!\left(-\frac{x^2}{2\sigma^2}\right)
$$

**Lorentzian broadening**:

$$
\varepsilon(\lambda) = \sum_{ia} f_{ia} \cdot L(\lambda - \lambda_{ia}, \gamma)
$$

$$
L(x, \gamma) = \frac{\gamma/\pi}{x^2 + \gamma^2}
$$

### Broadening Functions Compared

<SvgDiagram src="/svg/broadening-functions.svg" alt="Gaussian vs Lorentzian broadening: Gaussian has faster wing decay; Lorentzian has heavier tails" />

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sigma` | `f64` | `0.3` | Broadening width in eV (Gaussian σ or Lorentzian γ) |
| `e_min` | `f64` | `0.5` | Minimum transition energy in eV |
| `e_max` | `f64` | `12.0` | Maximum transition energy in eV |
| `n_points` | `usize` | `500` | Number of grid points in the energy axis |
| `broadening` | `BroadeningType` | `Gaussian` | `BroadeningType::Gaussian` or `BroadeningType::Lorentzian` |

## Output Types

### `StdaUvVisSpectrum`

| Field | Type | Description |
|-------|------|-------------|
| `wavelengths_nm` | `Vec<f64>` | Wavelength axis (nm) |
| `intensities` | `Vec<f64>` | ε(λ) — broadened extinction (arb. units) |
| `excitations` | `Vec<StdaExcitation>` | Individual sTDA excitations |
| `homo_energy` | `f64` | HOMO energy (eV) |
| `lumo_energy` | `f64` | LUMO energy (eV) |
| `gap` | `f64` | HOMO-LUMO gap (eV) |
| `n_transitions` | `usize` | Number of MO→MO transitions used |
| `broadening_type` | `BroadeningType` | Applied broadening type |
| `notes` | `Vec<String>` | Method caveats |

### `StdaExcitation`

| Field | Type | Description |
|-------|------|-------------|
| `energy_ev` | `f64` | Vertical excitation energy (eV) |
| `wavelength_nm` | `f64` | Corresponding wavelength (nm) |
| `oscillator_strength` | `f64` | Dimensionless oscillator strength $f$ |
| `from_mo` | `usize` | Occupied MO index |
| `to_mo` | `usize` | Virtual MO index |

## API

### Rust

```rust
use sci_form::{embed, compute_stda_uvvis};
use sci_form::reactivity::BroadeningType;

let conf = embed("c1ccc(cc1)C=O", 42);
let pos: Vec<[f64;3]> = conf.coords.chunks(3).map(|c| [c[0],c[1],c[2]]).collect();

let spectrum = compute_stda_uvvis(
    &conf.elements,
    &pos,
    0.3,             // sigma in eV
    1.0,             // e_min eV
    8.0,             // e_max eV
    500,             // n_points
    BroadeningType::Gaussian,
).unwrap();

println!("HOMO-LUMO gap: {:.2} eV", spectrum.gap);
println!("Strongest transition: {:.1} nm, f={:.4}",
    spectrum.excitations.iter()
        .max_by(|a,b| a.oscillator_strength.partial_cmp(&b.oscillator_strength).unwrap())
        .map(|e| e.wavelength_nm).unwrap_or(0.0),
    spectrum.excitations.iter()
        .map(|e| e.oscillator_strength).fold(0.0_f64, f64::max),
);
```

### Python

```python
from sci_form import stda_uvvis, embed

conf = embed("c1ccc(cc1)C=O", seed=42)
spec = stda_uvvis(conf.elements, conf.coords, sigma=0.3, e_min=1.0, e_max=8.0)
print(f"Gap: {spec.gap:.2f} eV, {spec.n_transitions} transitions")
print(f"λ_max ≈ {spec.wavelengths_nm[spec.intensities.index(max(spec.intensities))]:.1f} nm")
```

### TypeScript/WASM

```typescript
import init, { embed, compute_stda_uvvis } from 'sci-form-wasm';
await init();

const r = JSON.parse(embed('c1ccc(cc1)C=O', 42));
const elements = JSON.stringify(r.elements);
const coords   = JSON.stringify(r.coords);

const spec = JSON.parse(compute_stda_uvvis(elements, coords, 0.3, 1.0, 8.0, 500, 'gaussian'));
console.log(`Gap: ${spec.gap.toFixed(2)} eV`);
const maxIdx = spec.intensities.indexOf(Math.max(...spec.intensities));
console.log(`λ_max ≈ ${spec.wavelengths_nm[maxIdx].toFixed(1)} nm`);
```

## Limitations

- Excitation energies from EHT are systematically blue-shifted compared to TDDFT or experimental values; use for **trends and fingerprinting**, not absolute accuracy.
- Transition dipole moments are computed in the minimal STO basis; heavy elements may have reduced accuracy.
- No solvent effects are included.
- For transition-metal chromophores, GFN-xTB MOs are used automatically when EHT unsupported elements are detected.
