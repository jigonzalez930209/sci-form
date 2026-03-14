# Bounds Matrix

The bounds matrix encodes the geometric constraints between all atom pairs. It is the central data structure of distance geometry — all subsequent steps derive from it.

## Matrix Layout

For $N$ atoms, the bounds matrix $B$ is an $N \times N$ matrix:

- **Upper triangle** ($i < j$): upper bounds $u_{ij}$
- **Lower triangle** ($j < i$): lower bounds $l_{ij}$ (stored as $B[j][i]$)
- **Diagonal**: zero

<img src="/svg/bounds-matrix-overview.svg" alt="bounds-matrix-overview" class="svg-diagram" />

## 1-2 Bounds (Bonded Pairs)

For directly bonded atoms, the distance is the **bond length** computed from UFF (Universal Force Field) parameters:

$$r_{ij} = r_i + r_j + r_{\text{BO}} - r_{\text{EN}}$$

where:

**Bond order correction:**
$$r_{\text{BO}} = -0.1332(r_i + r_j) \ln(\text{BO})$$

**Electronegativity correction:**
$$r_{\text{EN}} = \frac{r_i \, r_j \left(\sqrt{\chi_i} - \sqrt{\chi_j}\right)^2}{\chi_i \, r_i + \chi_j \, r_j}$$

These terms come from the UFF parameterization. $r_i$ is the covalent radius, $\chi_i$ is the GMP electronegativity, and BO is the bond order.

The bounds are very tight: $[r_{ij} - 0.01,\; r_{ij} + 0.01]$ Å.

### Special Case: Conjugated Heterocycles

For conjugated bonds in 5-membered rings involving heteroatoms (like pyrrole, furan), an additional **squish factor** of 0.2 Å is applied to the lower bound to accommodate aromatic bond length variation:

$$l_{ij} = r_{ij} - 0.01 - 0.2$$

### UFF Parameters

Selected covalent radii and electronegativities:

| Element | $r_i$ (Å) | $\chi_i$ |
|---------|-----------|----------|
| H | 0.354 | 2.20 |
| C (SP3) | 0.757 | 2.55 |
| C (SP2) | 0.732 | 2.55 |
| C (SP) | 0.706 | 2.55 |
| N (SP3) | 0.700 | 3.04 |
| O (SP3) | 0.658 | 3.44 |
| F | 0.668 | 3.98 |
| S | 1.047 | 2.58 |
| Cl | 1.044 | 2.87 |
| Br | 1.214 | 2.74 |

## 1-3 Bounds (Angle Paths)

For atoms separated by two bonds (sharing a common neighbor), the distance is computed via the **law of cosines**:

<img src="/svg/bounds-matrix-triangle.svg" alt="bounds-matrix-triangle" class="svg-diagram" />

$$d_{13} = \sqrt{d_1^2 + d_2^2 - 2 d_1 d_2 \cos\theta}$$

The equilibrium angle $\theta$ depends on hybridization:

| Hybridization | Angle $\theta$ |
|---------------|------|
| SP | 180° |
| SP2 | 120° |
| SP3 | 109.47° |

The tolerance for 1-3 bounds is ±0.04 Å, doubled for each "large SP2 atom" in the path (atomic number > 13 that is SP2 and in a ring).

When multiple paths connect the same atom pair, **union semantics** apply: the lower bound is the minimum lower bound and the upper bound is the maximum upper bound across all paths.

## 1-4 Bounds (Torsion Paths)

For atoms separated by three bonds, the distance depends on the **torsion angle** between the two outer atoms. We compute the extreme distances at cis ($\phi = 0°$) and trans ($\phi = 180°$) configurations.

<img src="/svg/bounds-matrix-1-4.svg" alt="bounds-matrix-1-4" class="svg-diagram" />

The cis and trans distances are computed from planar geometry:

$$d_{\text{cis}} = \sqrt{(x_3 - x_0)^2 + (y_3 - y_0)^2}$$

where:
- $x_0 = r_{12}\cos\alpha$, $\;y_0 = r_{12}\sin\alpha$
- $x_3 = r_{23} - r_{34}\cos\beta$, $\;y_3 = r_{34}\sin\beta$

For the trans case, $y_3 = -r_{34}\sin\beta$.

### Stereochemistry Constraints

| Bond Type | Lower Bound | Upper Bound |
|-----------|-------------|-------------|
| E/Z stereo | forced cis or trans | same |
| SP2–SP2 (not in ring) | $d_{\text{cis}}$ | $d_{\text{trans}}$ |
| SP2–SP2 (in same ring) | $d_{\text{cis}}$ | $d_{\text{cis}}$ (planar) |
| Amide / Ester | $d_{\text{cis}}$ | $d_{\text{trans}}$ with bias |
| Default | min($d_{\text{cis}}, d_{\text{trans}}$) | max($d_{\text{cis}}, d_{\text{trans}}$) |

## VDW Bounds (Non-bonded)

For atom pairs with topological distance ≥ 4, the lower bound comes from **van der Waals radii**:

$$l_{ij} = f \cdot (r_{\text{vdw}}^i + r_{\text{vdw}}^j)$$

The scaling factor $f$ depends on topological distance:

| Topological Distance | Factor $f$ |
|---------------------|-----------|
| 4 | 0.70 |
| 5 | 0.85 |
| ≥6 | 1.00 |

The upper bound for non-bonded pairs defaults to $10^6$ Å (effectively unconstrained) until tightened by triangle smoothing.

## Triangle Inequality Smoothing

After populating bounds from all sources, **Floyd-Warshall** tightens them:

```
for k in 0..N:
    for i in 0..N:
        for j in (i+1)..N:
            u[i][j] = min(u[i][j], u[i][k] + u[k][j])
            l[i][j] = max(l[i][j], l[i][k] - u[k][j], l[k][j] - u[i][k])
            if l[i][j] > u[i][j]:
                INFEASIBLE
```

This is $O(N^3)$ and typically tightens many VDW upper bounds from $10^6$ to reasonable values.

### Fallback Strategy

If triangle smoothing detects infeasibility ($l_{ij} > u_{ij}$), sci-form tries a sequence of fallbacks:

1. **Strict bounds + set15** (1-5 bounds enabled) — most constrained
2. **No set15** — relax 1-5 constraints
3. **Soft smoothing** — allow $l_{ij} > u_{ij}$ with tolerance 0.05 Å

This ensures even strained molecules can produce initial coordinates, with the force field correcting violations later.
