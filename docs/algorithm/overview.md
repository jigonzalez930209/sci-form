# Pipeline Overview

sci-form implements the **ETKDGv2** (Experimental Torsion Knowledge Distance Geometry version 2) algorithm to generate 3D molecular conformers from SMILES strings. This page describes the complete pipeline at a high level. Each stage is covered in detail in its own section.

## The 9-Step Pipeline

<img src="/svg/pipeline-overview.svg" alt="pipeline-overview" class="svg-diagram" />

## Phase Breakdown

The pipeline can be understood in three conceptual phases:

### Phase 1: Topology → Bounds (Steps 1–3)

Starting from a SMILES string, we build a **molecular graph** and derive **distance constraints** between all atom pairs. These constraints form a bounds matrix $B$ where $l_{ij} \leq d_{ij} \leq u_{ij}$.

- **1-2 bounds**: Bond lengths from UFF parameters
- **1-3 bounds**: From bond lengths + equilibrium angles (law of cosines)
- **1-4 bounds**: From torsion angle cis/trans extremes
- **VDW bounds**: Lower bounds from van der Waals radii
- **Smoothing**: Floyd-Warshall to enforce the triangle inequality

→ Details: [Bounds Matrix](/algorithm/bounds-matrix)

### Phase 2: Embedding → Optimization (Steps 4–6)

We generate 3D (or 4D) coordinates from the distance constraints using **distance geometry**:

1. Pick random distances from the smoothed bounds
2. Convert distances to a **metric matrix** via the Cayley-Menger transform
3. Extract coordinates via **eigendecomposition**
4. Minimize distance violations using a **BFGS optimizer** with a bounds violation force field

$$T_{ij} = \frac{1}{2}\left(D_{0i} + D_{0j} - d_{ij}^2\right)$$

where $D_{0i} = \frac{1}{N}\sum_{k} d_{ik}^2 - \frac{1}{N^2}\sum_{k<l} d_{kl}^2$

→ Details: [Distance Geometry](/algorithm/distance-geometry), [Embedding](/algorithm/embedding)

### Phase 3: Refinement → Output (Steps 7–9)

After obtaining valid 3D coordinates, we refine them using the **ETKDG force field** that incorporates:

- **CSD torsion preferences**: 837 SMARTS patterns with Fourier coefficients derived from the Cambridge Structural Database
- **UFF inversions**: Out-of-plane energy for SP2 centers
- **Distance constraints**: Maintain bond lengths and angles
- **Validation**: Reject conformers with bad tetrahedral geometry, non-planar SP2 centers, or wrong double-bond configuration

→ Details: [ETKDG Refinement](/algorithm/etkdg-refinement), [Force Fields](/algorithm/force-fields)

## The Retry Loop

The embedding process is stochastic — not every random distance sample produces a valid conformer. The pipeline uses a **retry loop** with up to $10N$ iterations (where $N$ is the number of atoms):

1. If the metric matrix has zero or negative eigenvalues → **retry** with next RNG state
2. If energy per atom after bounds minimization exceeds 0.05 → **retry**
3. If tetrahedral centers fail the volume test → **retry**
4. If chiral volume signs don't match → **retry**
5. If planarity check fails → **retry**
6. If double-bond geometry is wrong → **retry**

Each retry advances the RNG state, producing different random distances. After $N/4$ consecutive failures, the algorithm falls back to **random box placement** (coordinates drawn uniformly from $[-5, 5]^3$).

## Algorithm Origin

The ETKDG algorithm was developed by Riniker and Landrum (2015) as an improvement to the DG (Distance Geometry) approach of Blaney and Dixon (1994). The key innovation is the use of **experimental torsion angle preferences** from the Cambridge Structural Database to guide conformer generation toward chemically realistic geometries.

**References:**
- Riniker, S.; Landrum, G. A. *J. Chem. Inf. Model.* **2015**, *55*, 2562–2574
- Blaney, J. M.; Dixon, J. S. *Rev. Comput. Chem.* **1994**, *5*, 299–335
- Wang, S.; Witek, J.; Landrum, G. A.; Riniker, S. *J. Chem. Inf. Model.* **2020**, *60*, 2044–2058
