# Validation

After minimization, each candidate conformer must pass a series of geometric validation checks. If any check fails, the conformer is rejected and the pipeline retries with a new random embedding.

## Validation Pipeline

<SvgDiagram src="/svg/validation-pipeline.svg" alt="validation-pipeline" />

## Check 1: Energy Per Atom

After bounds force field minimization, the energy per atom must be below a threshold:

$$\frac{E_{\text{total}}}{N} < 0.05$$

High energy indicates the optimizer couldn't satisfy the distance constraints — typically from a bad initial embedding where atoms are placed too close together or too far apart.

## Check 2: Tetrahedral Centers

For each carbon or nitrogen atom that:
- Has degree 4 (4 neighbors)
- Is in 2 or more SSSR rings
- Is NOT in any 3-membered ring

The volume of the tetrahedron formed by its 4 neighbors must be substantial:

<SvgDiagram src="/svg/validation-tetrahedral.svg" alt="validation-tetrahedral" />

### Volume Computation

For center atom $C$ with neighbors $A, B, D, E$:

1. Compute normalized direction vectors from each neighbor to the center:

$$\hat{v}_{AC} = \frac{\mathbf{x}_C - \mathbf{x}_A}{|\mathbf{x}_C - \mathbf{x}_A|}$$

2. Compute 4 triple products (one for each triple of neighbors):

$$V_1 = \hat{v}_{AC} \cdot (\hat{v}_{BC} \times \hat{v}_{DC})$$

3. All 4 triple products must exceed the threshold:

$$|V_k| > \begin{cases} 0.50 & \text{normal atoms} \\ 0.25 & \text{atoms in small rings} \end{cases}$$

4. Additionally, the center must be **inside** the tetrahedron — all 4 face normals must point outward relative to the center (tolerance: 0.30).

## Check 3: Chiral Volume Signs

For atoms with `@` or `@@` stereo specification:

$$\text{sign}(V_{\text{actual}}) = \text{sign}(V_{\text{expected}})$$

where:
- `@` (counterclockwise) → volume should be positive
- `@@` (clockwise) → volume should be negative

The volume is $V = \vec{v}_1 \cdot (\vec{v}_2 \times \vec{v}_3)$ with vectors computed from the chiral center to its neighbors, in the order specified by the SMILES.

A **20% tolerance** is applied — the absolute volume must be at least 20% of the target:

$$|V| > 0.2 \cdot |V_{\text{target}}|$$

## Check 4: Planarity

For SP2 centers (C=C, C=O, aromatic C, amide N), the out-of-plane energy is computed using UFF inversions:

$$E_{\text{oop}} = K \sum_{\text{impropers}} (1 - \sin Y)$$

where $Y$ is the Wilson out-of-plane angle.

**Rejection criterion:**

$$E_{\text{oop}} > N_{\text{improper}} \times 0.7$$

This accepts conformers where most SP2 centers are planar, rejecting only those with severely distorted geometry.

Additionally, SP-hybridized atoms (triple bonds, allenes, linear) are checked for linearity:

$$\theta_{\text{actual}} > 175°$$

## Check 5: Double-Bond Geometry

For each double bond, the substituents on each side must not be linear with the bond axis:

$$\cos\theta + 1 > 10^{-3}$$

where $\theta$ is the angle between a substituent and the double bond vector. A value near 180° (cos θ ≈ −1) would mean the substituent is collinear with the double bond — a physically impossible configuration.

<SvgDiagram src="/svg/validation-double-bond.svg" alt="validation-double-bond" />

## Summary of Thresholds

| Check | Threshold | What It Catches |
|-------|-----------|-----------------|
| Energy/atom | < 0.05 | Bad embedding, stuck optimization |
| Tetrahedral volume | > 0.50 (0.25 small ring) | Collapsed ring junctions |
| Chiral sign | ±20% tolerance | Wrong stereochemistry |
| Planarity | < 0.7 × $N_{\text{improper}}$ | Puckered SP2 centers |
| Double bond | $\cos\theta + 1 > 10^{-3}$ | Linear substituents |

## Retry Budget

The pipeline allows up to $10N$ total retry iterations, where $N$ is the number of atoms. For a typical drug-like molecule with 30 atoms, this gives 300 attempts — which is almost always sufficient.

In practice, most molecules succeed within 1–5 attempts. Molecules with many chiral centers or strained ring systems may need 10–50 attempts.
