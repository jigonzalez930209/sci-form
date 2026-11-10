# VCD and ROA screening models

The v2.1 spectroscopy prototype adds vibrational circular dichroism (VCD) and
Raman optical activity (ROA) outputs for the existing UFF normal modes. The
magnetic response is a **nuclear charge-current approximation**: the electronic
magnetic-dipole response is not solved self-consistently. Relative signs and
trends are the useful output.

> **Status.** Implemented in `feat/wave-a`. This is a qualitative screening
> model, not a coupled-perturbed magnetic response. See the
> [v2.x roadmap](/v2-roadmap).

## VCD

For a normal mode `k`, the electric and magnetic dipole derivatives are

$$
\boldsymbol{\mu}_k' = \sum_i q_i \mathbf{e}_{i,k},
\qquad
\mathbf{m}_k = \frac{1}{2}\sum_i q_i\,
\mathbf{R}_i \times \mathbf{e}_{i,k},
$$

where `e_{i,k}` is the normalized Cartesian displacement. The rotational
strength is

$$
R_k = \boldsymbol{\mu}_k'\cdot \mathbf{m}_k.
$$

This is the **nuclear charge-current** contribution. Real VCD requires the
electronic magnetic-dipole derivative from response theory; that is the main
missing physics.

## ROA

The ROA screening intensity is

$$
I_k^{\mathrm{ROA}} = R_k \cdot A_k^{\mathrm{Raman}},
$$

where `A_k` is the Raman activity from the existing prototype. This is a
qualitative product, not a full electric-dipole/magnetic-dipole invariant
calculation.

## API

```rust
use sci_form::vcd_roa::VcdRoaConfig;
use sci_form::{compute_vcd_spectrum, compute_roa_spectrum};

let vcd = compute_vcd_spectrum(smiles, &elements, &positions, &VcdRoaConfig::default())?;
for mode in vcd.modes {
    println!("{:.1} cm^-1  R={:.3e}", mode.frequency_cm1, mode.rotational_strength);
}

let roa = compute_roa_spectrum(smiles, &elements, &positions, &VcdRoaConfig::default())?;
```

### Output

| Type | Meaning |
|---|---|
| `VcdMode` | Frequency, electric/magnetic derivatives, rotational strength |
| `VcdSpectrum` | Modes and notes |
| `RoaMode` | Raman activity and screening circular intensity |
| `RoaSpectrum` | Modes, Raman spectrum and notes |

## Validation

| Gate | Observable | Result |
|---|---|---|
| Chiral molecule | Nonzero rotational strengths | 1-chloroethanol passes |
| ROA smoke | Finite circular intensities | passes (ignored by default: heavy) |
| Invalid geometry | Length mismatch | rejected |

The chiral test checks that the model produces a nonzero pseudoscalar. It does
not establish quantitative VCD accuracy.

## Limitations

- Nuclear charge-current model only; no electronic magnetic response.
- Gasteiger charges, not a QM charge density.
- UFF normal modes, not a QM Hessian.
- ROA is a qualitative product of Raman activity and the VCD term.
- No gauge correction, no magnetic-field perturbation, no coupled-perturbed HF.
- Absolute intensities are not calibrated; only relative signs/trends.
- The ROA test is `#[ignore]` because it runs two polarizability finite
  differences per mode.
- No Python, WASM or CLI wrapper.

## Next steps

1. Add coupled-perturbed magnetic-dipole response.
2. Use the same electronic method for Hessian and response.
3. Derive the full ROA invariants.
4. Benchmark against ORCA/Gaussian for chiral molecules.
5. Add VCD/ROA spectrum broadening.
6. Add Python, WASM and CLI wrappers.
