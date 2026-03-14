# Force Fields

sci-form uses two force fields in sequence: a **Bounds Violation Force Field** for the initial embedding stage, and an **ETKDG 3D Force Field** for the refinement stage. Both are minimized using the BFGS quasi-Newton optimizer.

## 1. Bounds Violation Force Field

This force field operates in 4D (or 3D) space and penalizes violations of the distance bounds. Its goal is to move atoms from the approximate eigendecomposition coordinates to positions that satisfy all pairwise distance constraints.

### Energy Function

$$E_{\text{total}} = E_{\text{dist}} + w_{\text{chiral}} \cdot E_{\text{chiral}} + w_{4D} \cdot E_{4D}$$

#### Distance Violation Term

For each pair $(i, j)$ with bounds $[l_{ij}, u_{ij}]$ and current distance $d_{ij}$, only considering pairs where $u_{ij} - l_{ij} \leq 5.0$ (the "basin" filter):

$$E_{\text{dist}} = \sum_{i<j} \begin{cases}
\left(\frac{d_{ij}^2}{u_{ij}^2} - 1\right)^2 & \text{if } d_{ij} > u_{ij} \\[8pt]
\left(\frac{2\,l_{ij}^2}{l_{ij}^2 + d_{ij}^2} - 1\right)^2 & \text{if } d_{ij} < l_{ij} \\[8pt]
0 & \text{otherwise}
\end{cases}$$

<img src="/svg/force-fields-overview.svg" alt="force-fields-overview" class="svg-diagram" />

The **basin filter** ($u - l \leq 5.0$) skips widely constrained pairs (typically distant non-bonded atoms) to focus optimization on well-determined distances.

#### Chiral Volume Term

For each chiral center with neighbors $(a, b, c, d)$:

$$E_{\text{chiral}} = \sum \begin{cases}
(V - V_{\text{upper}})^2 & \text{if } V > V_{\text{upper}} \\
(V - V_{\text{lower}})^2 & \text{if } V < V_{\text{lower}} \\
0 & \text{otherwise}
\end{cases}$$

where $V = \vec{v}_1 \cdot (\vec{v}_2 \times \vec{v}_3)$ is the signed volume of the tetrahedron formed by the four neighbor atoms, with $\vec{v}_k = \mathbf{x}_{n_k} - \mathbf{x}_{n_0}$.

#### 4D Penalty Term

$$E_{4D} = \sum_{i=1}^{N} x_{i4}^2$$

This penalizes the 4th dimension, driving coordinates toward a 3D embedding.

### Two-Phase Minimization

| Phase | $w_{\text{chiral}}$ | $w_{4D}$ | Purpose |
|-------|---------------------|----------|---------|
| 1 | 1.0 | 0.1 | Establish correct chirality |
| 2 | 0.2 | 1.0 | Collapse 4th dimension |

Phase 2 only runs when 4D embedding is used (molecule has chiral centers).

### BFGS Optimizer (Bounds FF)

The bounds force field uses a **full BFGS** optimizer with an $N_{dim} \times N_{dim}$ inverse Hessian approximation:

- **$N_{dim}$** = $N \times d$ (number of atoms × number of dimensions)
- **Maximum iterations**: 400 per pass
- **Force tolerance**: $10^{-3}$
- **Line search**: Cubic/quadratic interpolation backtracking (matching RDKit)
- **Gradient scaling**: $\times 0.1$, capped at 10.0 per component
- **Restarts**: up to 50 passes of 400 iterations each

The optimizer converges when the maximum gradient component falls below the force tolerance.

## 2. ETKDG 3D Force Field

After projecting to 3D, this force field refines the geometry using experimentally-derived torsion preferences.

### Energy Function

$$E_{\text{3D}} = E_{\text{torsion}} + E_{\text{inversion}} + E_{\text{dist12}} + E_{\text{angle}} + E_{\text{dist13}} + E_{\text{longrange}}$$

Terms are added in this exact order, matching RDKit's `construct3DForceField`.

#### M6 Torsion Terms

For each rotatable bond with matched CSD pattern:

$$E_{\text{torsion}} = \sum_{k=1}^{6} V_k\left(1 + s_k \cos(k\phi)\right)$$

where $V_k$ and $s_k$ are the Fourier coefficients from the CSD pattern. The cosines are computed using **Chebyshev polynomial recurrence** (no trigonometric functions beyond the initial $\cos\phi$):

$$
\begin{aligned}
\cos(2\phi) &= 2\cos^2\phi - 1 \\
\cos(3\phi) &= 4\cos^3\phi - 3\cos\phi \\
\cos(4\phi) &= 8\cos^4\phi - 8\cos^2\phi + 1 \\
\cos(5\phi) &= 16\cos^5\phi - 20\cos^3\phi + 5\cos\phi \\
\cos(6\phi) &= 32\cos^6\phi - 48\cos^4\phi + 18\cos^2\phi - 1
\end{aligned}
$$

<img src="/svg/force-fields-torsion.svg" alt="force-fields-torsion" class="svg-diagram" />

#### UFF Inversion Terms

For SP2 atoms (C, N, O) with 3 neighbors, out-of-plane bending is penalized:

$$E_{\text{inv}} = \frac{K_{\text{base}} \cdot 10}{3} \sum_{\text{perm}} (1 - \sin Y)$$

where $Y$ is the Wilson out-of-plane angle, and:
- $K_{\text{base}} = 50.0$ for C bonded to SP2 oxygen
- $K_{\text{base}} = 6.0$ for all other cases

Three permutations are evaluated for each improper center (cycling through the neighbor triples).

#### Distance Constraints

| Type | Force Constant $k$ | Application |
|------|-------------------|-------------|
| 1-2 (bond lengths) | 100.0 | All bonded pairs, flat-bottom ±0.01 Å |
| 1-3 (angle distances) | 100.0 | Pairs involving improper centers |
| Long-range | 10.0 | All remaining pairs from bounds matrix |

**Flat-bottom form:**
$$E = \frac{k}{2}(d - d_{\text{bound}})^2 \quad \text{when } d \text{ outside } [d_0 - \epsilon, d_0 + \epsilon]$$

#### Angle Constraints

For SP/linear atoms, constrain the angle to 180°:

$$E_{\text{angle}} = k(\theta - 180°)^2$$

### BFGS Optimizer (ETKDG 3D)

The 3D force field uses a similar BFGS optimizer:

- **Maximum iterations**: 300 (single pass)
- **Early termination**: skip if initial energy < $10^{-5}$
- **No restarts** — single optimization pass

## Gradient Computation

Both force fields use **analytical gradients** for efficient optimization. Key gradient formulas:

### Distance Term Gradient

For the distance between atoms $i$ and $j$:

$$\frac{\partial d_{ij}}{\partial \mathbf{x}_i} = \frac{\mathbf{x}_i - \mathbf{x}_j}{d_{ij}}$$

### Torsion Angle Gradient

The torsion angle $\phi$ defined by atoms $(a, b, c, d)$:

$$\frac{\partial \phi}{\partial \mathbf{x}_a} = \frac{|\mathbf{b}_1|}{|\mathbf{n}_1|^2} \mathbf{n}_1$$

$$\frac{\partial \phi}{\partial \mathbf{x}_d} = -\frac{|\mathbf{b}_2|}{|\mathbf{n}_2|^2} \mathbf{n}_2$$

where $\mathbf{b}_1 = \mathbf{x}_b - \mathbf{x}_a$, $\mathbf{b}_2 = \mathbf{x}_c - \mathbf{x}_d$, $\mathbf{n}_1 = \mathbf{b}_1 \times \mathbf{b}_{\text{mid}}$, $\mathbf{n}_2 = \mathbf{b}_2 \times \mathbf{b}_{\text{mid}}$.

### Chiral Volume Gradient

For the chiral volume $V = \vec{v}_1 \cdot (\vec{v}_2 \times \vec{v}_3)$:

$$\frac{\partial V}{\partial \mathbf{x}_{n_1}} = \vec{v}_2 \times \vec{v}_3$$

$$\frac{\partial V}{\partial \mathbf{x}_{n_2}} = \vec{v}_3 \times \vec{v}_1$$

$$\frac{\partial V}{\partial \mathbf{x}_{n_3}} = \vec{v}_1 \times \vec{v}_2$$

$$\frac{\partial V}{\partial \mathbf{x}_{center}} = -\left(\frac{\partial V}{\partial \mathbf{x}_{n_1}} + \frac{\partial V}{\partial \mathbf{x}_{n_2}} + \frac{\partial V}{\partial \mathbf{x}_{n_3}}\right)$$

## BFGS Algorithm

Both force fields use the **Broyden-Fletcher-Goldfarb-Shanno** quasi-Newton method:

$$\mathbf{H}_{k+1}^{-1} = \left(I - \frac{\mathbf{s}_k \mathbf{y}_k^T}{\mathbf{y}_k^T \mathbf{s}_k}\right) \mathbf{H}_k^{-1} \left(I - \frac{\mathbf{y}_k \mathbf{s}_k^T}{\mathbf{y}_k^T \mathbf{s}_k}\right) + \frac{\mathbf{s}_k \mathbf{s}_k^T}{\mathbf{y}_k^T \mathbf{s}_k}$$

where:
- $\mathbf{s}_k = \mathbf{x}_{k+1} - \mathbf{x}_k$ (step)
- $\mathbf{y}_k = \nabla f_{k+1} - \nabla f_k$ (gradient change)
- $\mathbf{H}_k^{-1}$ is the inverse Hessian approximation

### Line Search

The line search uses the following strategy:
1. **Initial step**: from the BFGS direction
2. **Backtracking**: cubic/quadratic interpolation
3. **Sufficient decrease**: Armijo condition with $c = 10^{-4}$
4. **Maximum step**: scaled to not exceed 100.0
