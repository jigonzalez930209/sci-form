# Embedding

This page covers the embedding step in detail: converting a smoothed bounds matrix into initial 3D (or 4D) atomic coordinates.

## Overview

The embedding process has four stages:

1. **Distance matrix construction** — pick random distances from bounds
2. **Metric matrix computation** — Cayley-Menger transform
3. **Eigendecomposition** — power iteration for top eigenpairs
4. **Coordinate extraction** — scale eigenvectors by eigenvalue square roots

## Step 1: Distance Picking

Given the smoothed bounds $[l_{ij}, u_{ij}]$, we construct a complete distance matrix by sampling:

$$d_{ij} = l_{ij} + r \cdot (u_{ij} - l_{ij})$$

where $r \sim \text{Uniform}(0, 1)$ from the **MinstdRand** linear congruential generator:

$$s_{n+1} = 48271 \cdot s_n \mod (2^{31} - 1)$$

$$r = \frac{s_{n+1}}{2^{31} - 1}$$

This is the same RNG as RDKit's `boost::minstd_rand`, ensuring bit-identical outputs for the same seed. The choice of LCG over Mersenne Twister for distance picking is deliberate — it's simple, fast, and matches the RDKit reference implementation.

::: info
The distances are sampled for all $\binom{N}{2}$ pairs, in the order $d_{01}, d_{02}, \dots, d_{0,N-1}, d_{12}, d_{13}, \dots$, matching RDKit's iteration order.
:::

## Step 2: Metric Matrix (Cayley-Menger Transform)

The **metric matrix** $T$ converts squared interpoint distances into inner products:

$$T_{ij} = \frac{1}{2}\left(D_{0i} + D_{0j} - d_{ij}^2\right)$$

where:

$$D_{0i} = \frac{1}{N}\sum_{k=1}^{N} d_{ik}^2 - \bar{d^2}$$

$$\bar{d^2} = \frac{1}{N^2}\sum_{k<l} d_{kl}^2$$

### Geometric Interpretation

If we place the centroid of all atoms at the origin, then:

$$T_{ij} = \mathbf{x}_i \cdot \mathbf{x}_j$$

This is the **Gram matrix**. Its eigenvalues correspond to the variance of coordinates along each principal axis, and its eigenvectors give the principal directions.

### Positive Definite Check

For exact Euclidean distances in $d$ dimensions, the metric matrix has:
- Exactly $d$ positive eigenvalues
- All remaining eigenvalues are zero

Since our distances are randomly sampled (not exact), $T$ may have small negative eigenvalues. We take the top 3 (or 4 for chiral molecules) positive eigenvalues and ignore the rest.

**Rejection criterion:** If $D_{0i} < 10^{-3}$ for any atom $i$, the metric matrix is degenerate and we reject this sample.

## Step 3: Power Iteration

Instead of computing the full eigendecomposition, we use **power iteration with deflation** to extract only the needed eigenpairs:

### Algorithm

For each eigenpair $k = 1, \dots, d$:

::: details Power Iteration Pseudocode
```
function power_iteration(T, max_iter=200, tol=1e-6):
    v = random_unit_vector(N)  // seeded from Mt19937
    
    for iter in 1..max_iter:
        w = T × v                    // matrix-vector product
        λ = v · w                     // Rayleigh quotient
        v_new = w / ||w||             // normalize
        if ||v_new - v|| < tol:
            break
        v = v_new
    
    return (λ, v)
```
:::

After extracting eigenpair $(\lambda_k, \mathbf{v}_k)$, we **deflate** the matrix:

$$T \leftarrow T - \lambda_k \, \mathbf{v}_k \, \mathbf{v}_k^T$$

This removes the contribution of the found eigenvector, so the next power iteration converges to the next-largest eigenvalue.

### RNG for Initial Vectors

The initial random vectors for power iteration use **Mt19937** (Mersenne Twister), seeded from the same base seed. This differs from the MinstdRand used for distance picking — Mt19937 provides better uniformity for high-dimensional random vectors.

## Step 4: Coordinate Extraction

The coordinates are recovered as:

$$x_{ik} = \sqrt{\lambda_k} \cdot v_{ik}$$

where:
- $i$ indexes atoms ($1 \leq i \leq N$)
- $k$ indexes spatial dimensions ($1 \leq k \leq d$)
- $\lambda_k$ is the $k$-th eigenvalue
- $v_{ik}$ is the $i$-th component of the $k$-th eigenvector

### Validation

After extraction, we check:

1. **All eigenvalues positive:** $\lambda_k > 10^{-3}$ for $k = 1, \dots, d$
2. **No NaN/Inf coordinates**
3. **Reasonable scale:** coordinates should not be excessively large

If any check fails, we reject this embedding and retry with the next RNG state.

## Complete Embedding Example

For a 4-atom molecule with bounds:

<img src="/svg/embedding-eigendecomp.svg" alt="embedding-eigendecomp" class="svg-diagram" />

## Random Box Fallback

After $N/4$ consecutive embedding failures, the algorithm switches to **random box placement**:

$$x_{id} \sim \text{Uniform}(-5, 5) \quad \text{for } d \in \{1, 2, 3\}$$

This always succeeds but produces a poor initial geometry. The bounds force field minimization then moves atoms to satisfy distance constraints from scratch. While slower, this ensures the pipeline never gets stuck on difficult molecules.
