# Electrostatic Potential (ESP) Maps

The electrostatic potential (ESP) describes the electric potential energy a unit positive test charge would experience at each point in space around a molecule. sci-form generates ESP maps on 3D voxel grids.

## ESP Grid Generation

<SvgDiagram src="/svg/esp-grid.svg" alt="ESP volumetric grid" />

### Coulomb Potential Formula

At each grid point $\vec{r}$, the ESP is the sum of Coulomb contributions from all atoms:

$$
\Phi(\vec{r}) = \sum_A \frac{q_A}{|\vec{r} - \vec{R}_A|}
$$

where:
- $q_A$ = partial charge on atom $A$ (from population analysis)
- $\vec{R}_A$ = position of atom $A$
- $|\vec{r} - \vec{R}_A|$ = distance from the grid point to atom $A$

### Grid Construction

The 3D grid is built from the molecular bounding box:

1. **Bounding box**: Find min/max atomic coordinates along each axis
2. **Padding**: Extend by a padding distance (default 3.0 Å) in all directions
3. **Spacing**: Place grid points at regular intervals (default 0.5 Å)
4. **Dimensions**: $N_x \times N_y \times N_z$ grid points

$$
N_i = \left\lfloor \frac{(\text{max}_i - \text{min}_i) + 2 \cdot \text{padding}}{\text{spacing}} \right\rfloor + 1
$$

### ESP Data Structure

```rust
pub struct EspGrid {
    pub origin: [f64; 3],     // grid corner
    pub spacing: f64,          // grid spacing in Å
    pub dims: [usize; 3],     // (Nx, Ny, Nz)
    pub values: Vec<f64>,     // Nx × Ny × Nz ESP values
    pub elements: Vec<u8>,    // atomic numbers
    pub positions: Vec<[f64; 3]>, // atom positions
}
```

## Parallel Grid Evaluation

The full 3D grid can involve millions of points. sci-form evaluates grid points in parallel using **Rayon** when the `parallel` feature is enabled:

```rust
// Serial (default)
let grid = compute_esp_grid(&mol, 0.5, 3.0);

// Parallel (requires feature = "parallel")
let grid = compute_esp_grid_parallel(&mol, 0.5, 3.0);
```

The parallel version partitions grid points into chunks across all available threads, giving near-linear speedup on multi-core systems.

## Color Mapping

ESP values are mapped to an RGB color scale using `esp_color_map(value, range)`:

$$
t = \text{clamp}\!\left(\frac{\Phi}{R}, -1, +1\right)
$$

$$
\text{RGB} = \begin{cases}
(255,\ (1+t)\cdot 255,\ (1+t)\cdot 255) & \text{if } t < 0 \quad (\text{red, nucleophilic}) \\
((1-t)\cdot 255,\ (1-t)\cdot 255,\ 255) & \text{if } t \geq 0 \quad (\text{blue, electrophilic}) \\
\end{cases}
$$

where $R$ is the saturation range (values beyond ±$R$ map to full red or full blue).

| ESP Value | Color | Chemical Meaning |
|-----------|-------|------------------|
| $\Phi \ll 0$ | Deep red | Strongly electron-rich, nucleophilic site |
| $\Phi < 0$ | Light red → white | Mildly electron-rich |
| $\Phi \approx 0$ | White | Neutral region |
| $\Phi > 0$ | Light blue → white | Mildly electron-poor |
| $\Phi \gg 0$ | Deep blue | Strongly electron-poor, electrophilic site |

The full grid can be converted to an RGB byte array:

```rust
use sci_form::esp::{esp_grid_to_colors};

let colors = esp_grid_to_colors(&grid, 0.05); // range = 0.05 a.u.
// returns Vec<u8> of length 3 × Nx × Ny × Nz  (r,g,b per point)
```

## Gaussian Cube Format

sci-form supports import/export of the `.cube` file format for interoperability with Gaussian, PySCF, and visualization tools.

### Cube File Structure

```
Title line 1
Title line 2
N_atoms  origin_x  origin_y  origin_z
N_x  dx_x  dx_y  dx_z
N_y  dy_x  dy_y  dy_z
N_z  dz_x  dz_y  dz_z
Z_1  charge_1  x_1  y_1  z_1
...
value_1  value_2  value_3  ...  (Nx × Ny × Nz float values)
```

Units in cube files are **Bohr** (1 Bohr = 0.529177 Å).

## API

### Rust

```rust
use sci_form::compute_esp;

let grid = compute_esp("O", 0.5, 3.0, None);
// grid.values: Vec<f64> — ESP at each grid point
// grid.dims: [usize; 3] — grid dimensions
// grid.to_cube() -> String — export as .cube format
```

### CLI

```bash
sci-form esp "O" --spacing 0.5 --padding 3.0
# Output: JSON with grid metadata and values
```

### Python

```python
import sci_form
grid = sci_form.esp("O", spacing=0.5, padding=3.0)
print(grid.dims)       # (Nx, Ny, Nz)
print(len(grid.values)) # Nx * Ny * Nz
cube_str = grid.to_cube()
```

## Validation

- **Consistency**: $N_x \times N_y \times N_z$ = `len(values)`
- **Finiteness**: All ESP values are finite (no infinities or NaN)
- **Polarity**: Polar molecules have both positive and negative ESP regions
- **Cube roundtrip**: Export to `.cube` and re-import produces identical grid
- **Resolution**: Finer spacing produces more grid points
