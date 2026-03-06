# Definitive RDKit Geometry & Force Field Implementation Guide

This document provides every technical detail extracted from the RDKit C++ core to enable a faithful replication of Distance Geometry and Force Field (UFF/MMFF94) algorithms in Rust.

---

---
# RDKit Force Field Research: Definitive Guide (Volume 1, 2 & 3)

## 1. Constants & Unit Conversions

| Constant | Value | Description |
| :--- | :--- | :--- |
| `MDYNE_A_TO_KCAL_MOL` | 143.9325 | Conversion from millidyne-Å to kcal/mol |
| `ELE_DIEL` | 332.0716 | Dielectric constant for electrostatic terms |
| `DEG2RAD` | $\pi / 180$ | Degrees to Radians |
| `RAD2DEG` | $180 / \pi$ | Radians to Degrees |
| `VDW_1` | 1.07 | Buffered 14-7 VdW constant |
| `VDW_2` | 1.12 | Buffered 14-7 VdW constant |

---

## 2. Distance Geometry (DG) Algorithms

### 2.1 Triangle Inequality Smoothing (Floyd-Warshall)
RDKit implementation in `TriangleSmooth.cpp`:
- **Upper Bound Check**: $U_{ij} \le U_{ik} + U_{kj}$
- **Lower Bound Check**: $L_{ij} \ge L_{ik} - U_{kj}$
- **Lower Bound Check 2**: $L_{ij} \ge L_{kj} - U_{ik}$
- Returns `false` if $L_{ij} > U_{ij}$ (impossible bounds).

### 2.2 Metric Matrix & Embedding
RDKit implementation in `DistGeomUtils.cpp`:
1. **Squared Distances**: $sqD_{ij} = d_{ij}^2$
2. **Double Centering**:
   - $sumSqD2 = \frac{1}{N^2} \sum_{i,j} sqD_{ij}$
   - $sqD0i = \frac{1}{N} \sum_j sqD_{ij} - sumSqD2$
   - $T_{ij} = \frac{1}{2} (sqD0i + sqD0j - sqD_{ij})$
3. **Eigen Decomposition**: $T = V \Lambda V^T$.
4. **Coords**: $X = V \Lambda^{1/2}$ (take top 3 eigenvalues).

### 2.3 Refinement Potentials
During refinement, specialized terms are used:
- **Distance Violation**: $E = w \cdot (d^2 / u^2 - 1)^2$ if $d > u$; $E = w \cdot (2l^2 / (l^2 + d^2) - 1)^2$ if $d < l$.
- **Chiral Violation**: $E = w \cdot (V - V_{target})^2$.
  - Volume $V = (\vec{r}_1 - \vec{r}_4) \cdot [(\vec{r}_2 - \vec{r}_4) \times (\vec{r}_3 - \vec{r}_4)]$.

---

## 3. Force Field & Optimizer (BFGS)

### 3.1 Optimizer (BFGS)
RDKit uses a robust BFGS in `BFGSOpt.h`:
- **Direction**: $\vec{p}_k = -H_k \nabla f_k$ (where $H_k$ is the inverse Hessian).
- **Line Search**: Cubic backtracking that ensures Wolf/Armijo conditions are satisfied.
- **Gradient Scaling Hack**: Scale gradient by $1.0 / \sqrt{|\nabla f|}$ if norm $> 1.0$ (highly important for numerical stability).

---

## 4. Rust Implementation Strategy

### 4.1 Data Structures
```rust
// Use faer for linear algebra (highly performant)
use faer::Mat;

pub trait ForceField {
    fn energy(&self, coords: &[f64]) -> f64;
    fn gradient(&self, coords: &[f64], grad: &mut [f64]) -> f64;
}

pub struct BoundsMatrix {
    pub matrix: Mat<f64>, // N x N
}

impl BoundsMatrix {
    pub fn smooth(&mut self) -> Result<(), &str> {
        // Floyd-Warshall implementation...
    }
}
```

### 4.2 Potential Terms in Rust
Each term (Bond, Angle, etc.) should be a `struct` implementing a `Contribution` trait:
```rust
pub trait Contribution {
    fn add_energy(&self, coords: &[f64]) -> f64;
    fn add_gradient(&self, coords: &[f64], grad: &mut [f64]);
}
```

---

## 5. Verification Roadmap
1. Generate RDKit dump of internal `T` (Metric Matrix) and compare against Rust.
2. Cross-verify `g_scale` behavior in `BFGS` steps.
3. Compare `check_grad.rs` outputs for both UFF and MMFF analytically vs numerically.
3. **Chiral Constraints**: Volume of the tetrahedra $V = \frac{1}{6} |(r_1 \cdot (r_2 \times r_3))|$.
   - Penalty: $E_{chiral} = w \cdot (V - V_{target})^2$ if $V$ is on the wrong side.

### B. Coordinate Embedding
- **Metric Matrix $G$**: Centered at origin. $G_{ij} = 0.5 \cdot (d_{i0}^2 + d_{j0}^2 - d_{ij}^2)$.
- **Power Eigen Solver**: RDKit uses an iterative power solver for the top 3 eigenvalues. Useful for large $N$.

---

## 5. Optimization & BFGS Hooks

### Gradient Scaling (Critical for Convergence)
```cpp
// Logic from ForceField.cpp
double gradScale = 0.1;
double maxG = find_max_abs_grad(grad);
if (maxG > 10.0) {
    while (maxG * gradScale > 10.0) {
        gradScale *= 0.5;
    }
}
// Final gradient used index-wise: grad[i] *= gradScale;
```

### BFGS Termination
- `FUNCTOL = 1e-4` (Energy change relative to value).
- `MOVETOL = 1e-7` (Coordinate change).
- `EPS = 3e-8` (Gradient magnitude).

---

## 6. Implementation Checklist for Rust

1. [ ] **SMARTS Typer**: Port `AtomTyper.cpp` logic. Essential for aromaticity and N-oxide detection.
2. [ ] **Matrix Ops**: Use `faer` or `nalgebra` for EVD.
3. [ ] **Triangle Smoothing**: Implement $O(N^3)$ triangle inequality or $O(N^2 \log N)$ sparse version.
4. [ ] **Force Field**: Implementation must be modular (contrib-based) to allow UFF and MMFF swapping.

---

## 7. Rust Integration Roadmap

### Recommended Data Structures:

```rust
trait ForceFieldContrib {
    fn get_energy(&self, pos: &[f64]) -> f64;
    fn get_grad(&self, pos: &[f64], grad: &mut [f64]);
}

struct ForceField {
    dimension: usize,
    positions: Vec<f64>, // Flat array [x1, y1, z1, x2, y2, z2, ...]
    contribs: Vec<Box<dyn ForceFieldContrib>>,
}

impl ForceField {
    fn calc_energy(&self) -> f64 {
        self.contribs.iter().map(|c| c.get_energy(&self.positions)).sum()
    }
    
    fn minimize(&mut self, max_its: usize) -> i32 {
        // Implement L-BFGS or BFGS here
    }
}
```

### Critical Implementation Details:
1.  **Atom Typing**: Port the `AtomTyper` logic. UFF uses SMARTS-like patterns for hybridization (e.g., `C_3` is tetrahedral carbon).
2.  **Distance Matrix**: RDKit uses a lazy-evaluated distance matrix for speed. In Rust, use `ndarray` or a simple symmetric flat vector.
3.  **Neighbor Matrix**: Pre-calculate 1-2, 1-3, and 1-4 relationships to exclude them from vdW calculations (usually 1-2 and 1-3 are excluded, 1-4 is scaled).

---

---

## 8. Volume 4: Current Gaps & Implementation Roadmap

To achieve full parity with RDKit, the following components are prioritized for the next development phases:

### 8.1 System de Tipado (Atom Typer) — **CRITICAL**
Currently, the mathematical potentials are implemented, but the logic to assign parameters from a molecular graph is missing:
- **UFF Typer**: Assigning types (e.g., `C_3`, `C_R`, `N_2`) based on hybridization, aromaticity, and connectivity.
- **MMFF94 Typer**: Complex engine requiring:
    - Aromatic ring detection (5 and 6 members).
    - Special types (e.g., `NPYL`, `C5A`).
    - Partial Bond Charge Increments (PBCI) and formal charge adjustments.

### 8.2 Missing Potentials & Parameters
- **UFF Van der Waals**: Implement Lennard-Jones with Waldman-Hagler combination rules.
- **MMFF94 Suite**:
    - **Bond Stretching**: Quartic potential.
    - **Angle Bending**: Cubic potential.
    - **Torsions**: Specific MMFF Fourier terms.
    - **Out-of-Plane**: Wilson inversion terms.
    - **Electrostatics**: Distance-dependent dielectric Coulomb term.
- **Parameter Tables**: Implementing the massive look-up tables for both UFF and MMFF94.

### 8.3 Distance Geometry (DG) Enhancements
- **ETKDG Knowledge Terms**: Knowledge-based 1-4, 1-5 distances and torsion constraints.
- **Gradual 4D to 3D Projection**: Implement the "squeeze" mechanism for the 4th dimension during refinement to bypass local minima.

### 8.4 Infrastructure "Glue"
- **ForceFieldBuilder**: A central high-level API that automates the flow: `Molecule` → `Typer` → `Parameter Search` → `MolecularForceField` → `Minimize`.

---

## 9. References
- **UFF**: Rappe, A. K., et al. "UFF, a full periodic table force field for molecular mechanics and molecular dynamics simulations." *J. Am. Chem. Soc.* 114.25 (1992): 10024-10035.
- **MMFF94**: Halgren, T. A. "Merck molecular force field. I. Basis, form, scope, parameterization, and performance of MMFF94." *J. Comput. Chem.* 17.5-6 (1996): 490-519.
- **RDKit Source**: `Code/ForceField`, `Code/DistGeom`, `Code/GraphMol/ForceFieldHelpers`.
