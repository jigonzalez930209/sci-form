# Extended Hückel Theory (EHT)

sci-form implements a full **Extended Hückel Theory** (EHT) pipeline for semiempirical electronic-structure calculations. EHT calculates molecular orbital (MO) energies, coefficients, and volumetric orbital densities from a 3D geometry without self-consistent field iteration.

## Overview

<img src="/svg/eht-pipeline.svg" alt="EHT pipeline phases B1–B5" class="svg-diagram" />

EHT was introduced by Hoffmann (1963) as an extension of Hückel theory to all valence electrons (not just π electrons). It requires only:

1. The molecular geometry (Cartesian coordinates)
2. A set of Valence State Ionization Potentials (VSIP) per element
3. Slater-type orbital exponents per element

The pipeline has five phases:

| Phase | Description |
|-------|-------------|
| B1 | Build atomic orbital basis from element parameters |
| B2 | Compute overlap matrix S and Hamiltonian H |
| B3 | Solve generalized eigenproblem HC = SCE |
| B4 | Evaluate MO densities on 3D volumetric grid |
| B5 | Extract isosurfaces via Marching Cubes |

---

## Phase B1: Atomic Basis Set

Each atom contributes one or two valence orbitals (s and p) described as **Slater-type orbitals (STOs)** expanded in a minimal STO-3G Gaussian basis.

### Slater-Type Orbitals

A Slater orbital of principal quantum number $n$ and exponent $\zeta$ (bohr⁻¹):

$$
\chi_{nlm}(r,\theta,\phi) = N \cdot r^{n-1} e^{-\zeta r} \cdot Y_l^m(\theta,\phi)
$$

where $N$ is a normalization constant and $Y_l^m$ are spherical harmonics.

### STO-3G Expansion

Each Slater orbital is approximated by 3 Gaussian primitives fitted to reproduce the STO shape:

$$
\chi(\mathbf{r}) = \sum_{k=1}^{3} d_k \cdot G_k(\alpha_k, \mathbf{r} - \mathbf{R}_A)
$$

where $G_k(\alpha_k, \mathbf{r}) = N_k \exp(-\alpha_k |\mathbf{r}|^2)$ is a normalized Gaussian primitive.

### VSIP Parameter Table

The following Hoffmann parameters are used as diagonal Hamiltonian elements:

| Element | Orbital | VSIP (eV) | Slater $\zeta$ (bohr⁻¹) |
|---------|---------|-----------|--------------------------|
| H  | 1s | −13.6  | 1.300 |
| C  | 2s | −21.4  | 1.625 |
| C  | 2p | −11.4  | 1.625 |
| N  | 2s | −26.0  | 1.950 |
| N  | 2p | −13.4  | 1.950 |
| O  | 2s | −32.3  | 2.275 |
| O  | 2p | −14.8  | 2.275 |
| F  | 2s | −40.0  | 2.425 |
| F  | 2p | −18.1  | 2.425 |
| Cl | 3s | −26.3  | 2.183 |
| Cl | 3p | −14.2  | 1.733 |
| Br | 4s | −22.07 | 2.588 |
| Br | 4p | −13.1  | 2.131 |
| I  | 5s | −18.0  | 2.679 |
| I  | 5p | −12.7  | 2.322 |

---

## Phase B2: Overlap Matrix and Hamiltonian

### Overlap Matrix $S$

$$
S_{ij} = \int \chi_i(\mathbf{r})\, \chi_j(\mathbf{r})\, d\mathbf{r}
$$

Computed analytically using the **Gaussian product theorem**: the product of two Gaussians centered at $\mathbf{A}$ and $\mathbf{B}$ is a single Gaussian centered at the weighted midpoint $\mathbf{P}$:

$$
G_a(\alpha, \mathbf{r} - \mathbf{A}) \cdot G_b(\beta, \mathbf{r} - \mathbf{B}) = K_{AB} \cdot G_{a+b}(\gamma, \mathbf{r} - \mathbf{P})
$$

where $\gamma = \alpha + \beta$, $\mathbf{P} = \frac{\alpha \mathbf{A} + \beta \mathbf{B}}{\gamma}$, and $K_{AB} = \exp\!\left(-\frac{\alpha\beta}{\gamma}|\mathbf{A}-\mathbf{B}|^2\right)$.

This allows exact analytical evaluation of all $\langle s|s \rangle$, $\langle s|p \rangle$, and $\langle p|p \rangle$ integrals without numerical quadrature.

### Hamiltonian Matrix $H$

**Diagonal elements** (Hückel approximation):

$$
H_{ii} = \text{VSIP}_i
$$

**Off-diagonal elements** (Wolfsberg-Helmholtz approximation):

$$
H_{ij} = \frac{1}{2} K \cdot S_{ij} \cdot (H_{ii} + H_{jj})
$$

where $K = 1.75$ is the **Wolfsberg-Helmholtz constant**. This is the key approximation of EHT: the resonance integral between two orbitals is proportional to their overlap and the average of their ionization potentials.

---

## Phase B3: Generalized Eigenproblem

The EHT secular equation is the generalized eigenproblem:

$$
\mathbf{H}\mathbf{C} = \mathbf{S}\mathbf{C}\mathbf{E}
$$

where $\mathbf{E}$ is a diagonal matrix of orbital energies and the columns of $\mathbf{C}$ are MO coefficient vectors. sci-form solves this using **Löwdin symmetric orthogonalization**:

### Step 1: Diagonalize S

$$
\mathbf{S} = \mathbf{U}\,\mathbf{\Lambda}\,\mathbf{U}^T
$$

### Step 2: Compute $S^{-1/2}$

$$
\mathbf{S}^{-1/2} = \mathbf{U}\,\text{diag}\!\left(\lambda_i^{-1/2}\right)\,\mathbf{U}^T
$$

Eigenvalues smaller than $10^{-10}$ are discarded to avoid numerical instability.

### Step 3: Transform to Standard Form

$$
\mathbf{H}' = \mathbf{S}^{-1/2}\,\mathbf{H}\,\mathbf{S}^{-1/2}
$$

### Step 4: Diagonalize $H'$

$$
\mathbf{H}'\mathbf{C}' = \mathbf{C}'\mathbf{E}
$$

### Step 5: Back-Transform Coefficients

$$
\mathbf{C} = \mathbf{S}^{-1/2}\,\mathbf{C}'
$$

Eigenvalues are sorted in ascending order. The HOMO and LUMO indices are identified from the number of valence electrons (filled two electrons per orbital).

### Output: EhtResult

```rust
pub struct EhtResult {
    pub energies: Vec<f64>,      // MO energies (eV), ascending
    pub coefficients: Vec<Vec<f64>>, // C matrix: coefficients[ao][mo]
    pub n_electrons: usize,       // total valence electrons
    pub homo_index: usize,        // HOMO (0-based)
    pub lumo_index: usize,        // LUMO (0-based)
    pub homo_energy: f64,         // eV
    pub lumo_energy: f64,         // eV
    pub gap: f64,                 // HOMO-LUMO gap (eV)
}
```

---

## Phase B4: Orbital Grid Evaluation

Once the MO coefficients are known, any orbital $\psi_m$ can be evaluated on a 3D volumetric grid:

$$
\psi_m(\mathbf{r}) = \sum_{\mu=1}^{N_\text{AO}} C_{\mu m}\, \chi_\mu(\mathbf{r})
$$

The **orbital density** (for visualization) is:

$$
\rho_m(\mathbf{r}) = |\psi_m(\mathbf{r})|^2
$$

The grid is parallelized over chunks using Rayon, evaluating all AO basis functions at each grid point.

### Grid Parameters

```rust
pub struct VolumetricGrid {
    pub origin: [f64; 3],    // lower-left corner (Å)
    pub spacing: f64,         // grid spacing (Å), default 0.15
    pub dims: [usize; 3],    // (Nx, Ny, Nz)
    pub values: Vec<f32>,    // scalar field (f32 for compact transfer)
}
```

---

## Phase B5: Marching Cubes Isosurface

Isosurfaces of MO density are extracted using the **Marching Cubes** algorithm (Lorensen & Cline 1987):

1. For each voxel cube in the grid, determine which vertices are inside/outside the isovalue threshold
2. Use the 256-case lookup table to identify the edge intersections
3. Interpolate vertex positions on crossed edges
4. Return a triangle mesh

The isosurface mesh is returned as:

```rust
pub struct IsosurfaceMesh {
    pub vertices: Vec<[f32; 3]>,    // xyz positions
    pub normals: Vec<[f32; 3]>,     // per-vertex normals
    pub indices: Vec<u32>,           // triangle indices
}
```

---

## API

### Rust

```rust
use sci_form::{eht_calculate, eht_orbital_mesh};

// Run EHT calculation
let result = eht_calculate("C6H6", None)?;  // benzene
println!("HOMO: {:.2} eV", result.homo_energy);
println!("LUMO: {:.2} eV", result.lumo_energy);
println!("Gap:  {:.2} eV", result.gap);

// Get orbital isosurface as JSON mesh
let mesh_json = eht_orbital_mesh("C6H6", result.homo_index, 0.05, 0.15, None)?;
```

### Python

```python
import sci_form

result = sci_form.eht_calculate("c1ccccc1")  # benzene (aromatic SMILES)
print(f"HOMO = {result.homo_energy:.2f} eV")
print(f"LUMO = {result.lumo_energy:.2f} eV")
print(f"Gap  = {result.gap:.2f} eV")
print(f"Coefficients shape: {len(result.coefficients)} x {len(result.coefficients[0])}")

# Get orbital grid for visualization
grid = sci_form.eht_orbital_mesh("c1ccccc1", result.homo_index, 0.05, 0.15)
```

### TypeScript / WASM

```typescript
import init, { eht_calculate, eht_orbital_grid_typed } from 'sci-form-wasm';

await init();

const result = JSON.parse(eht_calculate("c1ccccc1"));
console.log(`HOMO = ${result.homo_energy.toFixed(2)} eV`);
console.log(`Gap  = ${result.gap.toFixed(2)} eV`);

// Get orbital volume as a flat Float32Array (Nx × Ny × Nz)
const gridInfo = JSON.parse(compute_esp_grid_info("c1ccccc1", ...));
const orbital = eht_orbital_grid_typed("c1ccccc1", result.homo_index, 0.05, 0.15);
// orbital: Float32Array — ready for Three.js Data3DTexture
```

---

## Capabilities and Limitations

| Property | Value |
|----------|-------|
| Supported elements | H, C, N, O, F, Si, P, S, Cl, Br, I |
| Basis set | STO-3G (minimal, 1–2 orbitals per atom) |
| Correlation | None (single-determinant) |
| Relativistic effects | No |
| Self-consistency | No (one-shot diagonalization) |
| Cost | $O(N^3)$ diagonalization + $O(N^2)$ overlap |

EHT is appropriate for qualitative visualization of frontier orbitals (HOMO/LUMO), rough HOMO-LUMO gap estimates, and as a fast electronic descriptor. For quantitative accuracy, DFT or post-Hartree-Fock methods are required.

---

## References

1. Hoffmann, R. *An Extended Hückel Theory. I. Hydrocarbons.* J. Chem. Phys. **1963**, 39, 1397–1412.
2. Wolfsberg, M.; Helmholtz, L. *The Spectra and Electronic Structure of the Tetrahedral Ions MnO₄⁻, CrO₄²⁻, and ClO₄⁻.* J. Chem. Phys. **1952**, 20, 837–843.
3. Löwdin, P.-O. *On the Non-Orthogonality Problem Connected with the Use of Atomic Wave Functions in the Theory of Molecules and Crystals.* J. Chem. Phys. **1950**, 18, 365–375.
4. Lorensen, W. E.; Cline, H. E. *Marching Cubes: A High Resolution 3D Surface Construction Algorithm.* SIGGRAPH **1987**, 21, 163–169.
