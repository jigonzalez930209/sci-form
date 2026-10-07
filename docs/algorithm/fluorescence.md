# Fluorescence

Fluorescence is built from the same EHT sTDA singlets as the UV-Vis absorption spectrum. It is a vertical emission model with a spontaneous-emission lineshape and an Einstein A radiative rate. It does not optimize the S1 geometry and it does not report a quantum yield.

## Selection of the emitting state

Excitations are the occupied→virtual gaps of the EHT orbitals, including dark states. They are ordered by energy.

1. If the lowest singlet has oscillator strength \(f\) at or above `dark_threshold` (default \(10^{-3}\)), it is the emitter. That is Kasha's rule, and `kasha` is true.
2. If that state is dark, the lowest brighter singlet is used and `kasha` is false. The note records the departure.
3. If every singlet in the window is dark, the call fails. No band is invented.

## Lineshape

Absorption weights each band by \(f\). Emission does not reuse that envelope. The relative intensity is

\[
I(E) \propto n^2 E^3 |\mu|^2 g(E - E_\text{em})
\]

with \(|\mu|^2 \propto f / \Delta E_\text{abs}\) from \(f = \tfrac{2}{3}\Delta E|\mu|^2\). \(g\) is the same Gaussian or Lorentzian used for UV-Vis. The curve is scaled so its maximum is 1.

The emission energy is the vertical gap of the emitter minus `stokes_ev`. The default is 0: vertical emission, no excited-state relaxation. A positive `stokes_ev` is a shift you supply. It is not computed from an S1 optimization.

## Radiative rate

The Einstein A coefficient, with \(\tilde\nu\) in cm⁻¹, is

\[
A = 0.667025\, n^2 \tilde\nu^2 f \quad (\mathrm{s}^{-1})
\]

`radiative_lifetime_ns` is \(10^9 / A\). For \(f = 1\), \(E = 3.1\,\mathrm{eV}\) and \(n = 1\), the lifetime is about 2.4 ns. The fluorescence quantum yield is not reported, because non-radiative decay is not computed.

## API

Rust:

```rust
use sci_form::{compute_fluorescence, spectroscopy::FluorescenceConfig};

let mut config = FluorescenceConfig::default();
config.stokes_ev = 0.2;
let spec = compute_fluorescence(&elements, &positions, config)?;
```

Python: `fluorescence_spectrum(elements, coords, stokes_ev=0.2)`.

WASM: `compute_fluorescence(elements_json, coords_json, sigma, e_min, e_max, n_points, stokes_ev, refractive_index, broadening)`.

## What an organometallic spectrum is

Passing a 3D organometallic structure returns a spectrum when EHT can build the orbitals. That spectrum is a screening estimate.

- UV-Vis and fluorescence use minimal-basis EHT. Literature parameters cover the d-block already in the table. Every other metal through Z = 103 uses a hydrogenic fallback: H_ii is the first ionization potential and ζ = sqrt(IP / 13.6 eV). Lanthanides and actinides keep the f shell in the core (valence ns/np/(n−1)d), because the overlap is s/p/d. Shells with n = 6 and n = 7 reuse the STO-3G contraction of n = 5, scaled by (5/n)². Odd-electron molecules keep the SOMO occupied. None of this is TD-DFT.
- NMR shifts exist for the NMR-active isotopes in the nucleus registry, including the metal catalog. Cerium is omitted: every natural Ce isotope has \(I = 0\). Shifts other than ¹H/¹³C are relative screening values.
- IR frequencies use atomic masses through \(Z = 103\). The Hessian still needs a method that can evaluate the energy of those elements (`eht`, `pm3`, `xtb`, or `uff`).
