# Force Fields and Strain Energy (UFF)

sci-form implements the **Universal Force Field (UFF)** for molecular energy evaluation and geometry optimization. The force field describes the potential energy surface as a sum of bonded and non-bonded interaction terms.

## UFF Energy Terms

<img src="/svg/uff-energy-terms.svg" alt="UFF energy terms" class="svg-diagram" />

### Total Energy

$$
E_{\text{total}} = E_{\text{stretch}} + E_{\text{bend}} + E_{\text{torsion}} + E_{\text{vdW}} + E_{\text{oop}} + E_{\text{elec}}
$$

### 1. Bond Stretch (Harmonic)

$$
E_{\text{stretch}} = \frac{1}{2} k_{ij} (r - r_0)^2
$$

**Force constant** from UFF:

$$
k_{ij} = 664.12 \cdot \frac{Z_i^* Z_j^*}{r_0^3} \quad \text{(kcal/mol/Å}^2\text{)}
$$

**Natural bond length**:

$$
r_0 = r_i + r_j - r_{ij}(\text{EN}) - r_{ij}(\text{BO})
$$

where $r_{ij}(\text{EN})$ is the electronegativity correction and $r_{ij}(\text{BO})$ is the bond-order correction.

### 2. Angle Bend (Cosine Fourier)

$$
E_{\text{bend}} = k_\theta \left(C_0 + C_1 \cos\theta + C_2 \cos 2\theta\right)
$$

The Fourier coefficients $C_0, C_1, C_2$ are derived from the natural angle $\theta_0$:

| Hybridization | $\theta_0$ | Notes |
|---------------|------------|-------|
| sp³ | 109.47° | Tetrahedral |
| sp² | 120° | Trigonal planar |
| sp | 180° | Linear |

### 3. Torsion (Cosine)

$$
E_{\text{torsion}} = \frac{1}{2} V_\varphi \left(1 - \cos(n\varphi_0) \cos(n\varphi)\right)
$$

| Central bond | $n$ | $V_\varphi$ |
|-------------|-----|-------------|
| sp³–sp³ | 3 | $\sqrt{V_j V_k}$ |
| sp²–sp² | 2 | 5.0 kcal/mol |
| sp²–sp³ | 6 | 1.0 kcal/mol |

### 4. Van der Waals (Lennard-Jones 12-6)

$$
E_{\text{vdW}} = \varepsilon_{ij} \left[\left(\frac{x_{ij}}{r}\right)^{12} - 2\left(\frac{x_{ij}}{r}\right)^6\right]
$$

Geometric combining rules:

$$
x_{ij} = \sqrt{x_i x_j}, \quad \varepsilon_{ij} = \sqrt{\varepsilon_i \varepsilon_j}
$$

### 5. Inversion (Out-of-Plane)

For sp² centers with three neighbors:

$$
E_{\text{oop}} = k_\omega \left(C_0 + C_1 \cos\omega + C_2 \cos 2\omega\right)
$$

where $\omega$ is the angle the central atom makes with the plane of its three neighbors. For ideal sp² centers, $\omega_0 = 0°$ (planar).

### 6. Electrostatic (Coulomb)

$$
E_{\text{elec}} = 332.0637 \cdot \frac{q_i q_j}{\varepsilon \cdot r_{ij}}
$$

where $332.0637$ converts from elementary charges and Ångströms to kcal/mol.

## Geometry Optimization

sci-form uses **L-BFGS** for energy minimization:

1. Compute energy $E(\vec{x})$ and gradient $\nabla E(\vec{x})$
2. Approximate inverse Hessian from history of positions and gradients
3. Line search along descent direction
4. Repeat until $|\nabla E| < \text{tolerance}$

The gradient for each energy term is computed analytically, not numerically.

## Strain Energy

The **strain energy** is the total potential energy of the molecule in its current geometry. Higher strain energy indicates a geometry further from the force field's ideal.

Distorting bonds, angles, or torsions from their equilibrium values increases strain energy. This is useful for:

- Conformer quality assessment (lower energy → better conformer)
- Ring strain analysis
- Post-embedding refinement

## API

### Rust

```rust
use sci_form::compute_uff_energy;

let energy = compute_uff_energy("CCO", None);
// energy: f64 in kcal/mol
```

### CLI

```bash
sci-form energy "CCO"
# Output: total UFF energy in kcal/mol
```

### Python

```python
import sci_form
energy = sci_form.uff_energy("CCO")
print(energy)  # kcal/mol
```

## Validation

- **Finite energies**: All molecules produce finite (non-NaN, non-infinite) energies
- **Distortion increases energy**: Stretching bonds increases $E_{\text{stretch}}$
- **Finite gradients**: All gradient components are finite
- **Optimization convergence**: L-BFGS converges within max iterations for well-formed molecules
