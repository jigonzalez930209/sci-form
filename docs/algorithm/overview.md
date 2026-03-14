# Algorithm Overview

sci-form is a computational chemistry library with two main areas: **3D conformer generation** and **quantum-chemistry-inspired property computation**. This page describes both pipelines at a high level.

---

## Part 1: ETKDGv2 Conformer Pipeline

sci-form implements **ETKDGv2** (Experimental Torsion Knowledge Distance Geometry v2) to generate 3D molecular conformers from SMILES strings.

### The 9-Step Pipeline

<img src="/svg/pipeline-overview.svg" alt="pipeline-overview" class="svg-diagram" />

### Phase 1: Topology → Bounds (Steps 1–3)

Build a **molecular graph** and derive **distance constraints** between all atom pairs. Constraints form a bounds matrix $B$ where $l_{ij} \leq d_{ij} \leq u_{ij}$.

- **1-2 bounds**: Bond lengths from UFF parameters
- **1-3 bounds**: From bond lengths + equilibrium angles (law of cosines)
- **1-4 bounds**: From torsion angle cis/trans extremes
- **VDW bounds**: Lower bounds from van der Waals radii
- **Smoothing**: Floyd-Warshall to enforce the triangle inequality

→ Details: [Bounds Matrix](/algorithm/bounds-matrix), [SMILES Parsing](/algorithm/smiles-parsing)

### Phase 2: Embedding → Optimization (Steps 4–6)

Generate 3D (or 4D) coordinates from the distance constraints using **distance geometry**:

1. Pick random distances from the smoothed bounds (MinstdRand RNG)
2. Convert distances to a **metric matrix** via the Cayley-Menger transform
3. Extract coordinates via **eigendecomposition** (power iteration solver)
4. Minimize distance violations using **BFGS** with a bounds violation force field

$$T_{ij} = \frac{1}{2}\left(D_{0i} + D_{0j} - d_{ij}^2\right)$$

where $D_{0i} = \frac{1}{N}\sum_{k} d_{ik}^2 - \frac{1}{N^2}\sum_{k<l} d_{kl}^2$

→ Details: [Distance Geometry](/algorithm/distance-geometry), [Embedding](/algorithm/embedding)

### Phase 3: Refinement → Output (Steps 7–9)

After obtaining valid 3D coordinates, refine using the **ETKDG force field**:

- **CSD torsion preferences**: 846 SMARTS patterns with Fourier coefficients
- **UFF inversions**: Out-of-plane energy for SP2 centers
- **Distance constraints**: Maintain bond lengths and angles
- **Validation**: Reject conformers failing tetrahedral/planarity/stereo checks

→ Details: [ETKDG Refinement](/algorithm/etkdg-refinement), [Force Fields](/algorithm/force-fields), [Validation](/algorithm/validation)

### The Retry Loop

The pipeline uses a **retry loop** with up to $10N$ iterations:

1. Metric matrix has zero or negative eigenvalues → **retry**
2. Energy/atom after bounds minimization exceeds 0.05 → **retry**
3. Tetrahedral centers fail volume test → **retry**
4. Chiral volume signs don't match → **retry**
5. Planarity check fails → **retry**
6. Double-bond geometry is wrong → **retry**

After $N/4$ consecutive failures, fall back to **random box placement** (uniform from $[-5, 5]^3$).

---

## Part 2: Property Computation Pipeline

Once a conformer is generated, sci-form can compute a range of molecular properties.

### Gasteiger Charges (no QM)

Fast empirical electronegativity equalization — requires only topology + atomic numbers, runs in microseconds.

→ Details: [Population Analysis](/algorithm/population-analysis)

### Extended Hückel Theory (EHT)

EHT is the gateway to most quantum properties:

1. **Build STO-3G basis** — contracted Gaussian orbitals for each atom
2. **Overlap matrix S** — $S_{\mu\nu} = \langle \phi_\mu | \phi_\nu \rangle$
3. **Wolfsberg-Helmholtz H** — $H_{\mu\nu} = \frac{K}{2} S_{\mu\nu}(H_{\mu\mu} + H_{\nu\nu})$
4. **Löwdin ortho** — $\tilde{H} = S^{-1/2} H S^{-1/2}$, diagonalize → orbital energies $\varepsilon_i$ and MO coefficients
5. Fill $N_e$ electrons from lowest orbital up → HOMO/LUMO

From EHT:

| Property | Method | Key Output |
|----------|--------|------------|
| Population analysis | Mulliken/Löwdin | Per-atom charges |
| Dipole moment | Bond + lone-pair | Vector + Debye magnitude |
| Orbital grids | STO-3G on 3D grid | Float32 volumetric array |
| Isosurface mesh | Marching cubes | Vertices, normals, triangles |
| DOS/PDOS | Gaussian smearing | Total DOS, per-atom DOS |

→ Details: [Density of States](/algorithm/density-of-states), [Population Analysis](/algorithm/population-analysis), [Dipole Moments](/algorithm/dipole-moments)

### Electrostatic Potential

Coulomb ESP on a 3D grid from Mulliken charges:

$$V(\mathbf{r}) = \sum_i \frac{q_i^{\text{Mulliken}}}{|\mathbf{r} - \mathbf{r}_i|}$$

Color-mapped (red = negative, white = zero, blue = positive). Gaussian Cube export.

→ Details: [Electrostatic Potential](/algorithm/electrostatic-potential)

### Force Fields

- **UFF** — Universal Force Field, 50+ element types (including transition metals)
- **MMFF94** — Merck Molecular Force Field: quartic stretch, cubic bend, 3-term Fourier torsion, Halgren 14-7 vdW

→ Details: [Force Fields](/algorithm/force-fields), [Strain Energy](/algorithm/strain-energy)

### Molecular Alignment

Two algorithms:
- **Kabsch SVD** — optimal rotation, $O(N \cdot \min(N,3)^2)$
- **Quaternion alignment** — Coutsias 2004 4×4 eigenproblem, faster for large $N$

→ Details: [Molecular Alignment](/algorithm/molecular-alignment)

### Materials Assembly

Node/linker SBUs + topology → periodic crystal structure (MOF-type).

→ Details: [Materials Assembly](/algorithm/materials-assembly)

---

## References

- Riniker & Landrum, *J. Chem. Inf. Model.* **2015**, 55, 2562 (ETKDGv2)
- Wang et al., *J. Chem. Inf. Model.* **2020**, 60, 2044 (ETKDGv3)
- Blaney & Dixon, *Rev. Comput. Chem.* **1994**, 5, 299 (Distance geometry)
- Wolfsberg & Helmholtz, *J. Chem. Phys.* **1952**, 20, 837 (EHT)
- Mulliken, *J. Chem. Phys.* **1955**, 23, 1833 (Population analysis)
- Gasteiger & Marsili, *Tetrahedron* **1980**, 36, 3219 (Charge equalization)
- Coutsias et al., *J. Comput. Chem.* **2004**, 25, 1849 (Quaternion alignment)
- Halgren, *J. Comput. Chem.* **1996**, 17, 490 (MMFF94)
