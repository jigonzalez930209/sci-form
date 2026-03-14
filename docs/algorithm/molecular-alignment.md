# Molecular Alignment and RMSD (Kabsch Algorithm)

The Kabsch algorithm finds the optimal rotation to align two molecular structures, minimizing the root-mean-square deviation (RMSD). sci-form uses SVD-based Kabsch alignment for conformer comparison and benchmarking.

## Kabsch SVD Pipeline

<img src="/svg/kabsch-alignment.svg" alt="Kabsch alignment algorithm" class="svg-diagram" />

### Step 1: Center Both Structures

Given reference structure $P$ and target structure $Q$ (both $N \times 3$ matrices):

$$
\bar{P} = P - \text{centroid}(P), \quad \bar{Q} = Q - \text{centroid}(Q)
$$

where $\text{centroid}(P) = \frac{1}{N}\sum_{i=1}^N P_i$.

### Step 2: Cross-Covariance Matrix

Compute the $3 \times 3$ cross-covariance matrix:

$$
H = \bar{P}^T \bar{Q}
$$

### Step 3: Singular Value Decomposition

Decompose $H$ into:

$$
H = U \Sigma V^T
$$

### Step 4: Optimal Rotation

The rotation matrix that minimizes RMSD:

$$
d = \text{sign}(\det(V U^T))
$$
$$
R = V \cdot \text{diag}(1, 1, d) \cdot U^T
$$

The sign correction $d$ ensures a proper rotation (not a reflection).

### Step 5: RMSD Calculation

After applying rotation $R$:

$$
\text{RMSD} = \sqrt{\frac{1}{N} \sum_{i=1}^N |R \bar{P}_i - \bar{Q}_i|^2}
$$

## RMSD Properties

The RMSD metric satisfies:

| Property | Formula | Meaning |
|----------|---------|---------|
| Non-negative | $\text{RMSD} \geq 0$ | Always positive or zero |
| Identity | $\text{RMSD}(P, P) = 0$ | Identical structures |
| Symmetry | $\text{RMSD}(P, Q) = \text{RMSD}(Q, P)$ | Order doesn't matter |
| Triangle inequality | $\text{RMSD}(P, R) \leq \text{RMSD}(P, Q) + \text{RMSD}(Q, R)$ | Metric property |
| Rotation invariance | $\text{RMSD}(RP, Q) = \text{RMSD}(P, Q)$ | Optimal alignment removes rotation |

## RMSD Variants

sci-form supports:

- **All-atom RMSD**: Includes all atoms including hydrogens
- **Heavy-atom RMSD**: Excludes hydrogen atoms (more chemically meaningful for conformer comparison)
- **Type-specific RMSD**: Only certain element types (e.g., backbone atoms)

## Alignment Result

The `AlignmentResult` contains both the RMSD and the transformation:

```rust
pub struct AlignmentResult {
    pub rmsd: f64,
    pub rotation: [[f64; 3]; 3],
    pub translation_p: [f64; 3],
    pub translation_q: [f64; 3],
    pub aligned_coords: Vec<[f64; 3]>,
}
```

## API

### Rust

```rust
use sci_form::compute_rmsd;

let rmsd = compute_rmsd(
    &elements,
    &coords_reference,
    &coords_target,
    false, // heavy_atom_only
);
```

### CLI

```bash
sci-form rmsd --reference ref.xyz --target target.xyz
# Output: RMSD value and alignment metadata
```

### Python

```python
import sci_form
result = sci_form.rmsd(elements, coords_ref, coords_target)
print(result.rmsd)             # 0.12 Å
print(result.aligned_coords)   # aligned reference coordinates
```

## Benchmarking Use

RMSD is the primary metric for evaluating conformer quality:

- **< 0.5 Å**: Excellent agreement with reference
- **0.5–1.0 Å**: Good conformer, minor deviations
- **> 1.0 Å**: Significant geometric differences
- **> 2.0 Å**: Likely incorrect conformer

sci-form targets RMSD < 0.5 Å for the majority of molecules when compared to RDKit reference conformers.

---

## Quaternion Alignment (Coutsias Method)

sci-form also implements an independent **quaternion eigenvector method** for optimal rotation, based on Coutsias et al. (2004). This method is numerically more stable than SVD in near-degenerate cases and avoids explicit matrix decomposition.

### The Davenport Key Matrix

From the same cross-covariance matrix $R = \bar{P}^T \bar{Q}$ (where $S_{ij}$ denotes element $R_{ij}$), build the **4×4 symmetric key matrix** $F$:

$$
F = \begin{pmatrix}
S_{xx}+S_{yy}+S_{zz} & S_{yz}-S_{zy} & S_{zx}-S_{xz} & S_{xy}-S_{yx} \\
S_{yz}-S_{zy} & S_{xx}-S_{yy}-S_{zz} & S_{xy}+S_{yx} & S_{zx}+S_{xz} \\
S_{zx}-S_{xz} & S_{xy}+S_{yx} & -S_{xx}+S_{yy}-S_{zz} & S_{yz}+S_{zy} \\
S_{xy}-S_{yx} & S_{zx}+S_{xz} & S_{yz}+S_{zy} & -S_{xx}-S_{yy}+S_{zz}
\end{pmatrix}
$$

The **optimal rotation quaternion** $\mathbf{q} = (q_0, q_1, q_2, q_3)$ is the eigenvector of $F$ corresponding to its **largest eigenvalue** $\lambda_{\max}$.

This works because the RMSD objective can be rewritten as:

$$
\text{RMSD}^2 = \frac{1}{N}\left(\|\bar{P}\|_F^2 + \|\bar{Q}\|_F^2 - 2\lambda_{\max}\right)
$$

### Quaternion to Rotation Matrix

Given the unit quaternion $(q_0, q_1, q_2, q_3)$:

$$
R = \begin{pmatrix}
q_0^2+q_1^2-q_2^2-q_3^2 & 2(q_1q_2-q_0q_3) & 2(q_1q_3+q_0q_2) \\
2(q_1q_2+q_0q_3) & q_0^2-q_1^2+q_2^2-q_3^2 & 2(q_2q_3-q_0q_1) \\
2(q_1q_3-q_0q_2) & 2(q_2q_3+q_0q_1) & q_0^2-q_1^2-q_2^2+q_3^2
\end{pmatrix}
$$

This is a unit rotation matrix (det = +1, no reflection correction needed) because it is derived directly from a unit quaternion.

### Numerical Stability

The quaternion method avoids the sign-correction step required in Kabsch SVD and handles symmetric molecules (where multiple equivalent alignments exist) without numerical issues. Both methods are verified to produce identical RMSD values:

```rust
use sci_form::alignment::{align_coordinates, align_quaternion};

let kabsch = align_coordinates(&coords, &reference);
let quat   = align_quaternion(&coords, &reference);

assert!((kabsch.rmsd - quat.rmsd).abs() < 1e-10);
```

### API

```rust
use sci_form::alignment::align_quaternion;

let result = align_quaternion(&mobile_coords, &reference_coords);
println!("RMSD = {:.4} Å", result.rmsd);
println!("Rotation matrix: {:?}", result.rotation);
```

### Reference

- Coutsias, E. A.; Seok, C.; Dill, K. A. *Using quaternions to calculate RMSD.* J. Comput. Chem. **2004**, 25, 1849–1857.
