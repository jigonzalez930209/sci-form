# Dipole Moments

The molecular dipole moment $\vec{\mu}$ measures the overall polarity of a molecule. sci-form computes it from partial charges and atomic positions.

## Dipole Calculation

<img src="/svg/dipole-vector.svg" alt="Dipole moment vector" class="svg-diagram" />

### Formula

The molecular dipole moment is computed as the charge-weighted sum of atomic positions:

$$
\vec{\mu} = \sum_A q_A \, \vec{R}_A
$$

where $q_A$ is the partial charge on atom $A$ (from Mulliken or Löwdin analysis) and $\vec{R}_A$ is the Cartesian position of atom $A$.

### Unit Conversion

The raw dipole is in atomic units (e·Å). Conversion to **Debye** (the standard unit for molecular dipoles):

$$
|\vec{\mu}|_{\text{Debye}} = |\vec{\mu}|_{\text{a.u.}} \times 4.80321
$$

### Component Form

The dipole has three Cartesian components:

$$
\mu_x = \sum_A q_A \, x_A, \quad
\mu_y = \sum_A q_A \, y_A, \quad
\mu_z = \sum_A q_A \, z_A
$$

The magnitude is:

$$
|\vec{\mu}| = \sqrt{\mu_x^2 + \mu_y^2 + \mu_z^2}
$$

## Physical Interpretation

- **Nonpolar molecules**: $|\vec{\mu}| \approx 0$ D — symmetric charge distribution (CO₂, CH₄, H₂)
- **Polar molecules**: $|\vec{\mu}| > 0$ — asymmetric charge distribution
  - H₂O ≈ 1.85 D
  - HF ≈ 1.83 D
  - NH₃ ≈ 1.47 D
- **Direction**: Points from the center of negative charge toward the center of positive charge

## Conformer Dependence

The dipole moment depends on the 3D geometry. Different conformers of the same molecule can have different dipole magnitudes and orientations. sci-form computes the dipole for whatever conformer geometry is provided.

## API

### Rust

```rust
use sci_form::compute_dipole;

let result = compute_dipole("O", None);
// result.components: [f64; 3] — (μx, μy, μz)
// result.magnitude: f64 — |μ| in Debye
```

### CLI

```bash
sci-form dipole "O"
# Output: JSON with components and magnitude
```

### Python

```python
import sci_form
result = sci_form.dipole("O")
print(result.magnitude)    # ~1.85 D for water
print(result.components)   # [μx, μy, μz]
```

## Validation

- **H₂O**: magnitude in range 1.5–2.5 D
- **H₂**: magnitude ≈ 0 D (symmetric)
- **CH₄**: magnitude ≈ 0 D (tetrahedral symmetry)
- **HF**: dipole points from H toward F (more electronegative)
- **Consistency**: $|\vec{\mu}| = \sqrt{\mu_x^2 + \mu_y^2 + \mu_z^2}$ always holds
