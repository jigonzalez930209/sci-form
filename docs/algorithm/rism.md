# 3D-RISM solver core

The v2.3 solvation prototype implements the numerical core of the
three-dimensional reference interaction site model (3D-RISM): a periodic-grid
Ornstein–Zernike solver with the Kovalenko–Hirata closure and a 3D FFT. The
solute–solvent potential and the solvent susceptibility are caller-supplied.

> **Status.** Solver core, two simple susceptibility models and the mean
> interaction energy are implemented in `feat/wave-c`. A real solvent
> parameterisation, the excess chemical potential and a hydration-free-energy
> benchmark are follow-ups. See the [v2.x roadmap](/v2-roadmap).

## Equations

In Fourier space, for a one-component solvent,

$$
h(\mathbf{k}) =
\frac{c(\mathbf{k})}{1-\chi(\mathbf{k})c(\mathbf{k})},
$$

where `h` is the total correlation function, `c` the direct correlation
function and `χ` the solvent susceptibility.

The Kovalenko–Hirata closure is

$$
g(\mathbf{r}) =
\begin{cases}
\exp\left[-\beta u(\mathbf{r})+h(\mathbf{r})-c(\mathbf{r})\right],
& h-c \le 0,\\[4pt]
1-\beta u(\mathbf{r})+h(\mathbf{r})-c(\mathbf{r}),
& h-c > 0,
\end{cases}
$$

and

$$
c(\mathbf{r}) = g(\mathbf{r}) - 1 - h(\mathbf{r}).
$$

The solver iterates the closure and the Fourier-space OZ equation with Picard
mixing. The inverse FFT is normalized by `1/(nx·ny·nz)`.

## API

```rust
use sci_form::rism::{RismConfig, solve_3d_rism};

let config = RismConfig::default();
let n = config.nx * config.ny * config.nz;
let beta_u = vec![0.0; n];              // solute-solvent potential
let chi_k = vec![config.bulk_density_angstrom3; n]; // ideal-gas susceptibility

let result = solve_3d_rism(&beta_u, &chi_k, &config)?;
println!("converged={} mean_g={:.4} max_g={:.4}",
    result.converged, result.mean_g, result.max_g);
```

### Inputs

| Input | Meaning |
|---|---|
| `solute_beta_potential` | `βu(r)` on an `nx × ny × nz` grid |
| `solvent_susceptibility_k` | `χ(k)` on the same grid |
| `RismConfig` | Grid, box, temperature, density, mixing, tolerance |

For an ideal-gas solvent, use a constant `χ = ρ_bulk` on every k point. A real
solvent model must supply the Fourier-space susceptibility from its solvent
correlation functions.

### Solvent models

`solvent_susceptibility_k` provides two test models:

| Model | `χ(k)` |
|---|---|
| `IdealGas` | constant `ρ` |
| `Gaussian { correlation_length_angstrom }` | `ρ exp(-k²ξ²/2)` |

The Gaussian model is a numerical response for tests, not a parameterisation of
water. A real SPC/E/TIP3P model must compute `χ(k)` from its solvent correlation
functions.

```rust
use sci_form::rism::{RismSolventModel, solvent_susceptibility_k};

let chi = solvent_susceptibility_k(&config, RismSolventModel::Gaussian {
    correlation_length_angstrom: 3.0,
})?;
```

### Output

| Field | Meaning |
|---|---|
| `g`, `h`, `c` | Pair distribution, total and direct correlation arrays |
| `residual` | Maximum `|h_new − h|` at the last iteration |
| `iterations` | Picard iterations |
| `converged` | Tolerance reached |
| `mean_g`, `max_g` | Mean and contact-peak density |
| `mean_interaction_energy_kcal_mol` | Mean solute–solvent interaction energy |
| `shape` | Grid dimensions |

## Validation

| Gate | Observable | Result |
|---|---|---|
| Ideal gas | Mean `g` | 1.0 ± 1e-8 |
| Attractive solute | Contact peak | `max_g > 1` |
| Non-cubic grid | 4 × 6 × 8 ideal gas | runs |
| Ideal-gas χ | Constant density | constant |
| Gaussian χ | Decay with k | decreasing |
| Interaction energy | Finite for a test potential | finite |
| Wrong array length | Grid size mismatch | rejected |

These are solver-core tests. A solvent-model test and a hydration free-energy
comparison are still required.

## Bindings

| Surface | Entry point |
|---|---|
| CLI | `sci-form rism '<beta-u-json>' '<susceptibility-json>' '<RismConfig-json>'` |
| WASM | `solve_3d_rism(beta_u, susceptibility, config_json)` |
| Python | `solve_3d_rism(beta_u, susceptibility, config_json)` |

The grid arrays are row-major `nx × ny × nz` JSON arrays.
## Limitations

- No real solvent force field is bundled; the ideal-gas and Gaussian models are
  test responses, not water.
- `mean_interaction_energy_kcal_mol` is an interaction energy, not the excess
  chemical potential or a solvation free energy.
- No excess chemical potential or hydration free energy yet.
- No long-range charge treatment or Ewald split; the solver assumes a supplied
  periodic susceptibility.
- The FFT is a simple pure-Rust 3D transform; grid sizes are not optimized for
  large boxes.
- The KH closure only; HNC and other closures are follow-ups.
- No solute-site decomposition or spatial distribution output beyond `g`.
- No Python, WASM or CLI wrapper.
- No experimental or reference-code benchmark.

## Next steps

1. Add a real solvent model (SPC/E, TIP3P) that generates `χ(k)`.
2. Add the excess chemical potential and hydration free energy.
3. Add Ewald splitting for charged solutes.
4. Add HNC and other closures.
5. Add grid convergence and box-size tests.
6. Benchmark against Amber/Placevent or a reference 3D-RISM code.
7. Add Python, WASM and CLI wrappers.
