# Density of States (DOS) and Projected DOS

The density of states (DOS) describes how many electronic states exist at each energy level. sci-form computes DOS from EHT orbital energies using Gaussian smearing.

## DOS via Gaussian Smearing

<img src="/svg/dos-smearing.svg" alt="DOS Gaussian smearing" class="svg-diagram" />

### Total DOS Formula

Given a set of orbital energies $\{\varepsilon_i\}$ from the EHT solver, the density of states at energy $E$ is:

$$
D(E) = \sum_{i=1}^{N_{\text{orb}}} G(E - \varepsilon_i, \sigma)
$$

where $G$ is a Gaussian kernel:

$$
G(x, \sigma) = \frac{1}{\sigma \sqrt{2\pi}} \exp\left(-\frac{x^2}{2\sigma^2}\right)
$$

### Parameters

| Parameter | Symbol | Typical Range | Default |
|-----------|--------|---------------|---------|
| Smearing width | $\sigma$ | 0.1–0.5 eV | 0.3 eV |
| Energy grid points | $N_{\text{grid}}$ | 500–2000 | 1000 |
| Energy range | $[E_{\min}, E_{\max}]$ | Auto from eigenvalues ± 5σ | Auto |

### Properties

1. **Integral**: $\int D(E) \, dE \approx N_{\text{orbitals}}$
2. **Non-negative**: $D(E) \geq 0$ for all $E$
3. **σ dependence**: Smaller σ → sharper peaks near eigenvalue positions

## Projected DOS (PDOS)

PDOS decomposes the total DOS into contributions from individual atoms, revealing which atoms contribute to states at each energy.

### PDOS Formula

The projected density of states for atom $A$ at energy $E$:

$$
D_A(E) = \sum_{i=1}^{N_{\text{orb}}} w_{Ai} \cdot G(E - \varepsilon_i, \sigma)
$$

where the **Mulliken weight** for atom $A$ in orbital $i$ is:

$$
w_{Ai} = \sum_{\mu \in A} |C_{\mu i}|^2
$$

### PDOS Sum Rule

The atom-projected contributions always sum to the total DOS:

$$
\sum_A D_A(E) = D(E) \quad \forall E
$$

This is guaranteed because $\sum_A w_{Ai} = 1$ for each orbital $i$ (the coefficients squared sum to 1 in the orthonormal basis).

## Energy Grid

The energy grid is constructed automatically:

$$
E_{\min} = \min_i \varepsilon_i - 5\sigma, \quad E_{\max} = \max_i \varepsilon_i + 5\sigma
$$

Grid points are evenly spaced:

$$
E_k = E_{\min} + k \cdot \frac{E_{\max} - E_{\min}}{N_{\text{grid}} - 1}, \quad k = 0, \ldots, N_{\text{grid}}-1
$$

## API

### Rust

```rust
use sci_form::compute_dos;

let result = compute_dos("O", 0.3, 1000, None);
// result.energies: Vec<f64> — energy grid points
// result.total_dos: Vec<f64> — D(E) values
// result.pdos: Vec<Vec<f64>> — per-atom PDOS
// result.eigenvalues: Vec<f64> — orbital energies
```

### CLI

```bash
sci-form dos "O" --sigma 0.3 --grid-points 1000
# Output: JSON with energies, total_dos, pdos arrays
```

### Python

```python
import sci_form
result = sci_form.dos("O", sigma=0.3, grid_points=1000)
print(len(result.energies))   # 1000
print(len(result.total_dos))  # 1000
```

## Validation

- **Grid match**: `len(energies)` = `len(total_dos)` = `grid_points`
- **Non-negative**: All DOS values ≥ 0
- **PDOS sum**: Sum of PDOS arrays equals total DOS at each energy point
- **σ effect**: Narrower σ produces sharper, taller peaks
- **Energy range**: Grid spans all eigenvalues with sufficient margin

## DOS Comparison (MSE)

Two DOS curves can be quantitatively compared using the **Mean Squared Error**:

$$
\text{MSE}(D_A, D_B) = \frac{1}{N} \sum_{k=1}^{N} \left(D_A(E_k) - D_B(E_k)\right)^2
$$

This is used for benchmarking and to compare computed DOS against reference spectra. Both curves must be on the same energy grid (same length).

```rust
use sci_form::dos::{compute_dos, dos_mse};

let res_a = compute_dos("O", 0.3, 1000, None)?;
let res_b = compute_dos("O", 0.5, 1000, None)?;  // different sigma
let mse = dos_mse(&res_a.total_dos, &res_b.total_dos);
println!("MSE = {mse:.6}");
```

## JSON Export for Web Visualization

The `export_dos_json` function serializes a `DosResult` into a compact JSON string suitable for sending to a browser or plotting library:

```json
{
  "energies": [-35.2, -34.8, ..., 5.0],
  "total_dos": [0.0, 0.12, ..., 0.0],
  "sigma": 0.3,
  "pdos": {
    "0": [0.0, 0.08, ...],
    "1": [0.0, 0.04, ...]
  }
}
```

```rust
use sci_form::dos::{compute_dos, export_dos_json};

let result = compute_dos("c1ccccc1", 0.3, 1000, None)?;
let json = export_dos_json(&result);
// Send to frontend over WebSocket or WASM call
```

The PDOS dictionary keys are **atom indices** (0-based). This format is consumed directly by the TypeScript WASM API:

```typescript
const result = JSON.parse(compute_dos("c1ccccc1", 0.3, 1000));
// result.energies: number[]
// result.total_dos: number[]
// result.pdos: { [atomIndex: string]: number[] }
```
