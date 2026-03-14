# Distance Geometry

Distance Geometry (DG) is the mathematical framework for determining point positions from interpoint distances. In molecular conformer generation, we know approximate distance ranges between all atom pairs from the molecular topology, and we need to find 3D coordinates that satisfy these constraints.

## The Core Idea

Given $N$ atoms, we want to find coordinates $\mathbf{x}_1, \dots, \mathbf{x}_N \in \mathbb{R}^3$ such that:

$$l_{ij} \leq \|\mathbf{x}_i - \mathbf{x}_j\| \leq u_{ij} \quad \forall\, i < j$$

where $l_{ij}$ and $u_{ij}$ are the lower and upper distance bounds.

<img src="/svg/distance-geometry-pipeline.svg" alt="distance-geometry-pipeline" class="svg-diagram" />

## Distance Bounds Sources

The bounds matrix is populated from several sources, each corresponding to the topological distance between atoms:

| Topological Distance | Source | Precision |
|---------------------|--------|-----------|
| 1 (bonded) | UFF bond lengths | ±0.01 Å |
| 2 (1-3 path) | Law of cosines from bond angles | ±0.04 Å |
| 3 (1-4 path) | Torsion angle cis/trans extremes | computed |
| 4 (1-5 path) | Chained 1-4 distances | ±0.08 Å |
| ≥5 (non-bonded) | Van der Waals radii | 0.7–1.0× |

Details on each: [Bounds Matrix](/algorithm/bounds-matrix)

## The Triangle Inequality

For any three points in metric space:

$$d_{ij} \leq d_{ik} + d_{kj}$$

Applied to bounds, this means for all triples $(i, j, k)$:

$$
\begin{aligned}
u_{ij} &\leftarrow \min(u_{ij}, \, u_{ik} + u_{kj}) \\
l_{ij} &\leftarrow \max(l_{ij}, \, l_{ik} - u_{kj}, \, l_{kj} - u_{ik})
\end{aligned}
$$

These updates are applied iteratively via **Floyd-Warshall** until convergence. This tightens the bounds and ensures feasibility — if any $l_{ij} > u_{ij}$ after smoothing, the bounds are inconsistent (the atom arrangement is geometrically impossible).

## Distance Picking

Given smoothed bounds $[l_{ij}, u_{ij}]$, we pick a random distance for each pair:

$$d_{ij} = l_{ij} + r_{ij} \cdot (u_{ij} - l_{ij})$$

where $r_{ij} \sim \text{Uniform}(0, 1)$ from the MinstdRand LCG:

$$s_{n+1} = 48271 \cdot s_n \mod (2^{31} - 1)$$

This is the same RNG as RDKit's `boost::minstd_rand`, ensuring reproducible, bit-identical outputs for the same seed.

## From Distances to Coordinates

### The Metric Matrix (Cayley-Menger Transform)

Given a distance matrix $D$ where $D_{ij} = d_{ij}^2$, we construct the **metric matrix** $T$:

$$T_{ij} = \frac{1}{2}\left(D_{0i} + D_{0j} - D_{ij}\right)$$

where:

$$D_{0i} = \frac{1}{N}\sum_{k=1}^{N} D_{ik} - \frac{1}{N^2}\sum_{k<l} D_{kl}$$

Intuitively, $T_{ij} = \mathbf{x}_i \cdot \mathbf{x}_j$ when coordinates are centered at the centroid. The metric matrix is the **Gram matrix** of the coordinate vectors.

### Eigendecomposition

If the distances correspond to an exact Euclidean embedding in $d$ dimensions, then $T$ has exactly $d$ positive eigenvalues and the rest are zero.

$$T = V \Lambda V^T$$

The coordinates are recovered as:

$$x_{ik} = \sqrt{\lambda_k} \cdot v_{ik}$$

where $\lambda_k$ are the top $d$ eigenvalues and $v_{ik}$ are the corresponding eigenvector components.

<img src="/svg/distance-geometry-matrix.svg" alt="distance-geometry-matrix" class="svg-diagram" />

### Why 4D?

When the molecule has **chiral centers** (`@`/`@@` in SMILES), we embed in **4 dimensions** instead of 3. This gives the optimizer additional freedom to satisfy chiral volume constraints. After bounds minimization, the 4th dimension is collapsed:

1. Phase 1 bounds FF: $w_{\text{chiral}} = 1.0$, $w_{4D} = 0.1$ — establish chirality
2. Phase 2 bounds FF: $w_{\text{chiral}} = 0.2$, $w_{4D} = 1.0$ — collapse 4th dim
3. Take first 3 columns of the coordinate matrix

## Power Iteration

Instead of a full eigendecomposition ($O(N^3)$), sci-form uses **power iteration** to extract only the top $d$ eigenpairs:

1. Start with a random vector $\mathbf{v}^{(0)}$
2. Iterate: $\mathbf{v}^{(k+1)} = \frac{T \mathbf{v}^{(k)}}{\|T \mathbf{v}^{(k)}\|}$
3. Eigenvalue: $\lambda = \mathbf{v}^T T \mathbf{v}$
4. Deflate: $T \leftarrow T - \lambda \mathbf{v}\mathbf{v}^T$
5. Repeat for next eigenpair

This is more efficient for large molecules since we only need 3–4 eigenpairs, not all $N$.

## Rejection Criteria

Not every random distance sample yields a valid embedding. Rejections happen when:

| Condition | Meaning |
|-----------|---------|
| $D_{0i} < 10^{-3}$ for any $i$ | Degenerate metric matrix |
| $\lambda_k \leq 0$ for $k \leq d$ | Distances not embeddable in $d$D |
| Many consecutive failures | Switch to random-box fallback |

After $N/4$ consecutive embedding failures, the algorithm falls back to randomly placing atoms in a $[-5, 5]^3$ box and relying entirely on the force field minimization.
