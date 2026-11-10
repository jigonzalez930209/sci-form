# Raman spectroscopy from response polarizability

The first v2.x response-property prototype adds Raman activities for the
existing UFF normal modes. It is deliberately scoped: the polarizability is an
uncoupled HF/STO-3G sum-over-states tensor and the Raman intensity is a finite
difference along each normal mode. VCD and ROA need magnetic-dipole response and
remain follow-ups.

> **Status.** Implemented in `feat/wave-a`. This is a screening prototype, not a
> calibrated Raman cross-section. See the [v2.x roadmap](/v2-roadmap).

## Why a response polarizability

IR intensities need the derivative of the dipole moment along a normal mode.
Raman needs the derivative of the **polarizability** tensor. The two spectra
therefore share the normal modes but not the response property.

The first implementation avoids a new coupled-perturbed SCF by using the
uncoupled sum-over-states expression

$$
\alpha_{ij} =
\sum_{\mathrm{spin}}
\sum_{a}^{\mathrm{occ}}
\sum_{b}^{\mathrm{vir}}
\frac{
\langle a|\mu_i|b\rangle
\langle b|\mu_j|a\rangle
}{
\epsilon_b-\epsilon_a
}.
$$

The dipole integrals are evaluated in the existing STO-3G basis and transformed
to the UHF molecular orbitals. The result is in Å³ after the Bohr³ conversion.
This is the uncoupled HF approximation; electron correlation and coupled
response are not included.

## Raman activity

For each real normal mode, displace the geometry by `±δ` along the normalized
mass-weighted mode, compute the polarizability at both points, and form

$$
\alpha_k' = \frac{\alpha(\delta)-\alpha(-\delta)}{2\delta},
\qquad
\beta_k' = \frac{\beta(\delta)-\beta(-\delta)}{2\delta}.
$$

The Placzek activity is

$$
A_k = 45(\alpha_k')^2 + 7(\beta_k')^2,
$$

and the depolarization ratio is

$$
\rho_k =
\frac{3(\beta_k')^2}
{45(\alpha_k')^2 + 4(\beta_k')^2}.
$$

Relative activities are meaningful. Absolute cross-sections are not calibrated.

## API

```rust
use sci_form::raman::RamanConfig;
use sci_form::compute_raman_spectrum;

let water = [
    [0.0, 0.0, 0.0],
    [0.9572, 0.0, 0.0],
    [-0.239987, 0.927297, 0.0],
];

let spectrum = compute_raman_spectrum(
    "O",
    &[8, 1, 1],
    &water,
    &RamanConfig::default(),
)?;

for mode in spectrum.modes {
    println!("{:.1} cm^-1  activity={:.3}", mode.frequency_cm1, mode.activity);
}
```

### Output

| Type | Meaning |
|---|---|
| `Polarizability` | 3×3 tensor in Å³, isotropic and anisotropic invariants |
| `RamanMode` | Mode index, frequency cm⁻¹, activity and depolarization |
| `RamanSpectrum` | Modes, input polarizability, step and notes |

## Validation

| Gate | Observable | Result |
|---|---|---|
| Tensor symmetry | `α_ij = α_ji` | < 1e-10 |
| Water polarizability | Isotropic value | positive and finite |
| Water Raman smoke | Real modes and finite activities | passes |
| Invalid step | `δ = 0` | rejected |

The water test checks that the response path runs and produces a symmetric
tensor. A quantitative benchmark against ORCA/Gaussian Raman intensities is
still required before promotion.

## Limitations

- The polarizability is uncoupled HF/STO-3G, not a coupled-perturbed or
  correlated response.
- The normal modes come from the UFF Hessian, not the same electronic method as
  the polarizability.
- VCD and ROA are not implemented; they need magnetic-dipole derivatives.
- The finite-difference step is a fixed 0.01 Å. A step convergence study is
  required for production.
- Absolute Raman cross-sections are not calibrated; only relative activities
  are meaningful.
- The tensor origin is the coordinate origin. For neutral molecules the
  polarizability is origin-independent; charged systems need a documented
  origin convention.

## Next steps

1. Replace the uncoupled sum-over-states tensor with a coupled-perturbed HF
   response.
2. Use the same electronic method for the Hessian and the polarizability.
3. Add VCD and ROA magnetic-dipole response.
4. Add a step-size convergence study.
5. Benchmark against ORCA/Gaussian for a small molecule set.
6. Add Python, WASM and CLI wrappers.
