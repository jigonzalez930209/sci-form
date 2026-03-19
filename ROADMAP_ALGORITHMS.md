# Roadmap: Algorithm Implementation Plan for sci-form

**Goal:** Implement 10 strategic algorithm tracks to advance sci-form into a competitive computational chemistry platform. Each track adds production-grade algorithms with validation tests enforcing <0.5% deviation against reference implementations.

**Validation standard:** All numerical algorithms must pass automated tests with max deviation < 0.5% against established reference values (PySCF, NIST, RDKit, OpenBabel, or peer-reviewed literature).

---

## Implementation Status Overview

| Track | Description | Status |
|-------|-------------|--------|
| H1 | NMR Spectroscopy | ✅ Complete |
| H2 | Vibrational Analysis & IR Spectroscopy | 🟡 Mostly complete (H2.1c pending) |
| H3 | Charge Modeling & Force-Field Parametrization | 🟡 H3.1 complete, H3.2 partial |
| H4 | Conformational Sampling Refinement | 🟡 H4.1 complete, H4.2 pending |
| H5 | Topology & Molecular Descriptors | ✅ Complete |
| H6 | Molecular Dynamics & Thermodynamics | 🔲 Not started (enhancements) |
| H7 | Implicit Solvation Modeling | ✅ Complete |
| H8 | Quantum Mechanics: Iterative SCF Evolution | 🟡 H8.1b/H8.2a partial |
| H9 | Stereochemistry & Chirality Perception | 🟡 H9.1/H9.2a-b complete, H9.2c pending |
| H10 | Surface Generation & Scientific Visualization | 🟡 H10.1a/H10.1d complete, rest pending |

---

## Track H1: NMR Spectroscopy

Enhance the existing NMR module (`src/nmr/`) with higher-fidelity algorithms for chemical shift prediction, spin-spin coupling, and spectral simulation.

### Phase H1.1: HOSE Code Enhancement (Hierarchical Organisation of Spherical Environments)

**Algorithm:** Graph-traversal encoding of the concentric chemical environment around each atom up to a configurable radius. Each sphere layer records element, bond order, and branching topology, producing a canonical string descriptor.

**Current state:** Basic HOSE code generation exists in `src/nmr/hose.rs`.

- [x] **H1.1a** Extend HOSE code generation to support sphere radii up to 5 (currently limited)
- [x] **H1.1b** Add canonical ordering within each sphere (Morgan-algorithm-based tie-breaking)
- [x] **H1.1c** Embed a compile-time HOSE → chemical-shift lookup database (using `phf` or `include!` macro) covering ¹H and ¹³C shifts from NMRShiftDB2 reference data
- [x] **H1.1d** Implement database-backed shift prediction: look up the most specific HOSE match (radius 4 → 3 → 2 → 1 fallback) to predict δ(ppm)
- [x] **H1.1e** Validate: ¹H shift predictions within ±0.3 ppm MAE on a 50-molecule organic test set; ¹³C within ±3.0 ppm MAE

### Phase H1.2: 3D Karplus Equation for J-Couplings

**Algorithm:** Trigonometric evaluation of dihedral angles from 3D conformer coordinates to compute vicinal spin-spin coupling constants using the Karplus relation:

$$
{}^3J = A\cos^2\phi + B\cos\phi + C
$$

where $A$, $B$, $C$ are empirical parameters depending on the coupling pathway (H-C-C-H, H-C-O-H, etc.).

**Current state:** Basic Karplus implementation exists in `src/nmr/coupling.rs`.

- [x] **H1.2a** Expand parameterization tables for all common H-C-X-H pathways (X = C, N, O)
- [x] **H1.2b** Add Altona-modified Karplus parameters for sugar/nucleotide systems
- [x] **H1.2c** Add ensemble averaging: compute ³J over multiple conformer geometries and Boltzmann-weight the result
- [x] **H1.2d** Validate: ³J(H,H) predictions within ±0.5 Hz MAE on ethane, cyclohexane, and butane reference values

### Phase H1.3: Lorentzian Broadening for NMR Spectra

**Algorithm:** Convolve discrete chemical shifts with Lorentzian line shapes:

$$
L(\nu) = \frac{1}{\pi} \frac{\gamma/2}{(\nu - \nu_0)^2 + (\gamma/2)^2}
$$

where $\gamma$ is the full-width at half-maximum.

**Current state:** Lorentzian broadening implemented in `src/nmr/spectrum.rs`.

- [x] **H1.3a** Add configurable FWHM per nucleus type (¹H: 1 Hz, ¹³C: 5 Hz defaults)
- [x] **H1.3b** Add first-order multiplet splitting from J-couplings (doublet, triplet, quartet patterns)
- [x] **H1.3c** Add integration (relative area) calculation for each peak group
- [x] **H1.3d** Validate: spectral peak positions match input shifts to <0.01 ppm; multiplet patterns match analytical Pascal's triangle splitting

---

## Track H2: Vibrational Analysis & IR Spectroscopy

Harden the existing IR module (`src/ir/`) with improved accuracy and analytical gradient support.

### Phase H2.1: Numerical Hessian (Central Finite Differences)

**Algorithm:** For each of the $3N$ Cartesian coordinates, displace by $\pm h$ and evaluate the energy to approximate second derivatives:

$$
H_{ij} \approx \frac{E(x_i + h, x_j + h) - E(x_i + h, x_j - h) - E(x_i - h, x_j + h) + E(x_i - h, x_j - h)}{4h^2}
$$

**Current state:** Numerical Hessian implemented in `src/ir/hessian.rs` using central differences.

- [x] **H2.1a** Add automatic step-size selection based on atomic mass and force-field stiffness
- [x] **H2.1b** Add symmetry enforcement: $H_{ij} = H_{ji}$ averaging to reduce numerical noise
- [ ] **H2.1c** Add analytical Hessian path for UFF bond-stretch and angle-bend terms (avoiding 6N energy evaluations for these)
- [x] **H2.1d** Validate: Hessian symmetry error $\|H - H^T\|_F / \|H\|_F < 10^{-6}$

### Phase H2.2: Mass-Weighted Diagonalization and Normal Modes

**Algorithm:** Construct the mass-weighted Hessian $\tilde{H}_{ij} = H_{ij} / \sqrt{m_i m_j}$, diagonalize to get eigenvalues $\lambda_k$ and eigenvectors (normal modes), then convert to vibrational frequencies:

$$
\tilde{\nu}_k = \frac{1}{2\pi c} \sqrt{\lambda_k}
$$

**Current state:** Mass-weighted normal mode analysis implemented in `src/ir/vibrations.rs`.

- [x] **H2.2a** Add translation/rotation mode removal (project out 6 zero-frequency modes for nonlinear, 5 for linear molecules)
- [x] **H2.2b** Add zero-point vibrational energy (ZPVE) with anharmonic correction factor (0.9 × harmonic)
- [x] **H2.2c** Add thermochemical properties: thermal energy, entropy, and Gibbs free energy at 298.15 K using the rigid-rotor harmonic-oscillator (RRHO) approximation
- [x] **H2.2d** Validate: H₂O frequencies within 5% of NIST values (3657, 1595, 3756 cm⁻¹); CO₂ within 5% of NIST (667, 1388, 2349 cm⁻¹)

### Phase H2.3: IR Spectrum Generation with Intensity Calculation

**Algorithm:** IR intensities from the derivative of the dipole moment along each normal mode:

$$
I_k \propto \left| \frac{\partial \vec{\mu}}{\partial Q_k} \right|^2
$$

**Current state:** Numerical dipole derivative intensities implemented.

- [x] **H2.3a** Add spectral assignment: label peaks with functional-group annotations (C=O stretch, O-H stretch, etc.) using frequency ranges
- [x] **H2.3b** Add Gaussian broadening option alongside existing Lorentzian
- [x] **H2.3c** Add transmittance (%T) output in addition to absorbance
- [x] **H2.3d** Validate: key functional group frequencies (C=O ~1700 cm⁻¹, O-H ~3400 cm⁻¹, C-H ~2900 cm⁻¹) within 10% of experimental values

---

## Track H3: Charge Modeling & Force-Field Parametrization

### Phase H3.1: Gasteiger-Marsili Enhancement (PEOE)

**Algorithm:** Iterative partial equalization of orbital electronegativity. On each iteration, charge flows along bonds proportional to the electronegativity difference, damped by a decay factor:

$$
\Delta q_{A \to B}^{(k)} = \frac{\chi_B(q_B) - \chi_A(q_A)}{\chi_B^+(q_B)} \cdot d^k
$$

where $\chi(q) = a + bq + cq^2$ and $d$ is the damping factor (typically 0.5).

**Current state:** Gasteiger-Marsili fully implemented in `src/charges/gasteiger.rs` for 12 elements.

- [x] **H3.1a** Extend element coverage to all main-group elements up to period 4 (Li, Be, Na, Mg, Al, etc.)
- [x] **H3.1b** Add configurable damping factor and convergence threshold
- [x] **H3.1c** Add formal-charge-aware initialization for charged species (protonated amines, carboxylates)
- [x] **H3.1d** Validate: charges on ethanol, acetic acid, benzene within 0.5% of OpenBabel Gasteiger reference

### Phase H3.2: SMARTS-Based Atom Typing for Force Fields

**Algorithm:** Hierarchical substructure matching using SMARTS patterns to assign atom types for MMFF94 and UFF force fields. Each atom is classified by iterating through a priority-ordered list of SMARTS patterns until the first match.

**Current state:** Basic SMARTS parser and matcher exist in `src/smarts/`. MMFF94 atom typing (28 types) exists in `src/forcefield/`.

- [ ] **H3.2a** Add complete MMFF94 atom-type SMARTS rules with priority ordering (75 rules)
- [x] **H3.2b** Add UFF atom-type assignment covering coordination number and hybridization
- [ ] **H3.2c** Add validation: atom-type assignment for 20 diverse molecules matches RDKit MMFF94 typing exactly
- [x] **H3.2d** Add automatic parameter lookup after typing: bond/angle/torsion/vdW parameters from typed atoms

---

## Track H4: Conformational Sampling Refinement

### Phase H4.1: Butina Clustering by RMSD

**Algorithm:** Greedy distance-based clustering (Taylor-Butina):
1. Compute the all-pairs RMSD matrix after Kabsch alignment
2. For each conformer, count neighbors within a distance threshold $d_{\text{cut}}$
3. Select the conformer with the most neighbors as a cluster centroid
4. Remove all its neighbors from the pool; repeat until empty

- [x] **H4.1a** Implement Butina clustering with configurable RMSD cutoff (default 1.0 Å)
- [x] **H4.1b** Add O(N²) RMSD matrix computation with Kabsch alignment from `src/alignment/`
- [x] **H4.1c** Return cluster centroids, cluster sizes, and representative conformer indices
- [ ] **H4.1d** Integrate into `embed_batch` pipeline as a post-processing filter
- [x] **H4.1e** Validate: clustering of 50 ethane conformers yields ≤3 clusters (gauche, anti, eclipsed)

### Phase H4.2: Torsional Sampling (Systematic + Simulated Annealing)

**Algorithm:** Identify rotatable bonds, enumerate torsion angles at fixed intervals (e.g., 30°), and evaluate force-field energy at each combination. For molecules with many rotatable bonds (>4), switch to simulated annealing with Metropolis acceptance:

$$
P(\text{accept}) = \min\left(1, e^{-\Delta E / k_B T}\right)
$$

- [ ] **H4.2a** Add rotatable-bond detection (exclude ring bonds, double bonds, terminal groups)
- [ ] **H4.2b** Implement systematic rotor search for ≤4 rotatable bonds (12 angles per rotor)
- [ ] **H4.2c** Implement simulated annealing with exponential cooling schedule for >4 rotors
- [ ] **H4.2d** Integrate with Butina clustering to filter the generated ensemble
- [ ] **H4.2e** Validate: n-butane sampling finds gauche (±60°) and anti (180°) minima; energy difference within 10% of MMFF94 reference

---

## Track H5: Topology & Molecular Descriptors

### Phase H5.1: Smallest Set of Smallest Rings (SSSR)

**Algorithm:** Identify the fundamental cycle basis of the molecular graph. Uses Horton's or Vismara's algorithm:
1. Compute shortest paths between all pairs of atoms
2. For each edge not in the spanning tree, extract the associated fundamental cycle
3. Reduce to the smallest linearly independent set of rings

**Current state:** Ring detection exists via the `petgraph` cycle basis, but SSSR-specific logic is not fully canonical.

- [x] **H5.1a** Implement canonical SSSR using Horton's algorithm (shortest-path-based cycle extraction)
- [x] **H5.1b** Add ring-size histogram and ring-membership annotation per atom/bond
- [x] **H5.1c** Add aromaticity perception: Hückel's rule (4n+2 electrons) applied to SSSR rings with alternating single/double bonds or full conjugation
- [x] **H5.1d** Validate: benzene → 1 ring (size 6, aromatic); naphthalene → 2 rings (both size 6, aromatic); cubane → 5 rings (all size 4, non-aromatic)

### Phase H5.2: Morgan Algorithm & ECFP Fingerprints

**Algorithm:** Extended-Connectivity Fingerprints (ECFP, Rogers & Hahn 2010):
1. Initialize each atom with an invariant (atomic number, degree, H-count, charge, ring membership)
2. For each iteration (radius), hash the atom's current identifier with its neighbors' identifiers
3. Collect all unique identifiers across all radii as the fingerprint bit set
4. Fold to a fixed-length bit vector (1024 or 2048 bits) via modular hashing

- [x] **H5.2a** Implement atom-invariant initialization (element, degree, valence, ring, aromaticity, formal charge)
- [x] **H5.2b** Implement iterative neighborhood aggregation with deterministic hashing (FNV-1a or xxHash)
- [x] **H5.2c** Add configurable radius (ECFP4 = radius 2, ECFP6 = radius 3)
- [x] **H5.2d** Add bit-vector folding to 1024/2048 bits with Tanimoto similarity function
- [x] **H5.2e** Validate: Tanimoto(benzene, toluene) > 0.5; Tanimoto(benzene, hexane) < 0.2; self-similarity = 1.0

---

## Track H6: Molecular Dynamics & Thermodynamics

### Phase H6.1: Velocity Verlet Integration (Enhancement)

**Algorithm:** Symplectic integration of Newton's equations:

$$
\begin{aligned}
v(t + \tfrac{1}{2}\Delta t) &= v(t) + \tfrac{1}{2}a(t)\Delta t \\
x(t + \Delta t) &= x(t) + v(t + \tfrac{1}{2}\Delta t)\Delta t \\
a(t + \Delta t) &= F(x(t+\Delta t)) / m \\
v(t + \Delta t) &= v(t + \tfrac{1}{2}\Delta t) + \tfrac{1}{2}a(t+\Delta t)\Delta t
\end{aligned}
$$

**Current state:** Velocity Verlet with UFF already implemented in `src/dynamics.rs`.

- [ ] **H6.1a** Add support for PM3 and xTB force backends (not just UFF)
- [ ] **H6.1b** Add energy conservation monitoring (drift < 0.1% over 1000 steps for NVE)
- [ ] **H6.1c** Add trajectory output in XYZ format for external visualization
- [ ] **H6.1d** Validate: harmonic oscillator test — period matches analytical $T = 2\pi\sqrt{m/k}$ within 0.5%

### Phase H6.2: Nosé-Hoover Thermostat

**Algorithm:** Extended system dynamics that couples the physical system to a thermal reservoir:

$$
\dot{\xi} = \frac{1}{Q}\left(\sum_i \frac{p_i^2}{m_i} - (3N-6)k_BT\right)
$$

where $\xi$ is the thermostat variable and $Q$ is the thermostat mass.

**Current state:** Berendsen thermostat implemented; Nosé-Hoover not yet.

- [ ] **H6.2a** Implement Nosé-Hoover chain thermostat (chain length 3) for rigorous NVT sampling
- [ ] **H6.2b** Add temperature monitoring and fluctuation analysis
- [ ] **H6.2c** Add Nosé-Hoover coupling to velocity Verlet integration
- [ ] **H6.2d** Validate: average temperature within 1% of target; temperature fluctuations $\propto 1/\sqrt{N_{DOF}}$

---

## Track H7: Implicit Solvation Modeling

### Phase H7.1: Shrake-Rupley SASA (Enhancement)

**Current state:** Full Shrake-Rupley implementation in `src/surface/sasa.rs`.

- [x] **H7.1a** Add per-element solvation parameters (atomic solvation parameter, ASP)
- [x] **H7.1b** Add non-polar solvation energy: $\Delta G_{\text{np}} = \sum_i \sigma_i \cdot A_i$ where $\sigma_i$ is the ASP
- [x] **H7.1c** Add effective Born radii computation from SASA burial fraction
- [x] **H7.1d** Validate: SASA of methane within 0.5% of analytical sphere area; water SASA matches literature

### Phase H7.2: Generalized Born (GB) Model

**Algorithm:** Compute electrostatic solvation energy using the Still equation:

$$
\Delta G_{\text{GB}} = -\frac{1}{2}\left(\frac{1}{\epsilon_{\text{in}}} - \frac{1}{\epsilon_{\text{out}}}\right) \sum_{i,j} \frac{q_i q_j}{f_{\text{GB}}(r_{ij}, R_i, R_j)}
$$

where $f_{\text{GB}} = \sqrt{r_{ij}^2 + R_iR_j \exp(-r_{ij}^2/4R_iR_j)}$ and $R_i$ are effective Born radii.

- [x] **H7.2a** Implement effective Born radii from pairwise descreening (HCT model)
- [x] **H7.2b** Implement Still GB equation for electrostatic solvation free energy
- [x] **H7.2c** Add SA+GB combined solvation energy
- [x] **H7.2d** Validate: solvation energies for small molecules (water, methanol, acetic acid) within 1 kcal/mol of SMD/PCM reference

---

## Track H8: Quantum Mechanics — Iterative SCF Evolution

### Phase H8.1: Self-Consistent Field (SCF) Cycle

**Current state:** Roothaan-Hall SCF implemented in `src/hf/scf.rs`; PM3 SCF in `src/pm3/solver.rs`; xTB SCC in `src/xtb/`.

- [ ] **H8.1a** Unify SCF infrastructure: extract a common `ScfSolver` trait for HF, PM3, and xTB backends
- [x] **H8.1b** Add configurable convergence criteria (energy: 1e-8 Eh, density: 1e-6)
- [ ] **H8.1c** Add level shifting for difficult-to-converge metallic systems ($\Delta = 0.3$ Eh shift on virtuals)
- [ ] **H8.1d** Validate: H₂ restricted HF energy within 0.1% of PySCF reference; H₂O within 0.5%

### Phase H8.2: DIIS Convergence Acceleration

**Algorithm:** Direct Inversion in the Iterative Subspace — maintain a history of Fock matrices and error vectors, solve the Pulay equation to extrapolate the optimal Fock matrix:

$$
\min_{\{c_i\}} \left\| \sum_i c_i e_i \right\|^2 \quad \text{s.t.} \quad \sum_i c_i = 1
$$

**Current state:** DIIS implemented in `src/hf/scf.rs`.

- [x] **H8.2a** Add DIIS history size control (default 8, max 15 matrices)
- [ ] **H8.2b** Add ADIIS (energy-based DIIS) for initial iterations before switching to DIIS
- [ ] **H8.2c** Add automatic DIIS reset on numerical instability (condition number monitoring)
- [ ] **H8.2d** Validate: metalloporphyrin SCF converges in <50 iterations with DIIS vs >200 without

---

## Track H9: Stereochemistry & Chirality Perception

### Phase H9.1: Cahn-Ingold-Prelog (CIP) Priority Assignment

**Algorithm:** Depth-first graph traversal with priority rules:
1. Higher atomic number → higher priority
2. At tie, explore outward and compare next-shell atoms
3. Phantom atoms (lone pairs) have priority 0
4. Double bonds → duplicate the bonded atom in the priority tree

- [x] **H9.1a** Implement CIP priority tree construction via DFS with duplicate-atom expansion for multiple bonds
- [x] **H9.1b** Handle tie-breaking with recursive exploration (up to depth 10)
- [x] **H9.1c** Validate: alanine → S-configuration; D-glucose → R at C2

### Phase H9.2: R/S and E/Z Assignment from 3D Coordinates

**Algorithm:** Using CIP priorities and 3D coordinates:
- **R/S:** Compute the signed volume of the tetrahedron formed by the 4 substituents. Positive volume with CIP ordering 1→2→3 (4th behind) = R.
- **E/Z:** Compute the dihedral angle across the double bond between highest-priority substituents. Same side = Z, opposite = E.

- [x] **H9.2a** Implement R/S assignment using the triple-product sign test
- [x] **H9.2b** Implement E/Z assignment for double bonds using the dihedral angle of CIP-ordered substituents
- [ ] **H9.2c** Add stereo descriptor output in SMILES string form (@, @@, /, \)
- [x] **H9.2d** Validate: (R)-2-bromobutane → R; trans-2-butene → E; cis-2-butene → Z

---

## Track H10: Surface Generation & Scientific Visualization

### Phase H10.1: Enhanced Marching Cubes

**Current state:** Marching cubes implemented in `src/eht/marching_cubes.rs` for orbital isosurfaces.

- [x] **H10.1a** Generalize marching cubes to accept any `VolumetricGrid` (not just EHT orbitals)
- [ ] **H10.1b** Add support for dual-phase isosurfaces (positive/negative lobes of orbitals, ESP coloring)
- [ ] **H10.1c** Add mesh simplification: vertex welding and degenerate-triangle removal
- [x] **H10.1d** Validate: sphere isosurface area matches analytical $4\pi r^2$ within 2%

### Phase H10.2: Vertex Normal Computation

**Algorithm:** For each vertex, compute the normal as the weighted average of face normals of all adjacent triangles:

$$
\hat{n}_v = \text{normalize}\left(\sum_{f \ni v} \alpha_f \hat{n}_f\right)
$$

where $\alpha_f$ is the angle at vertex $v$ in face $f$.

- [ ] **H10.2a** Implement angle-weighted vertex normals for smooth shading
- [ ] **H10.2b** Add normal flipping for consistent outward orientation
- [ ] **H10.2c** Export mesh data in WASM-ready format (interleaved position+normal Float32Arrays)
- [ ] **H10.2d** Validate: normals on sphere surface match analytical $\hat{r}$ to cosine similarity > 0.99

---

## Priority Order (Implementation Sequence)

1. **Track H1** — NMR Spectroscopy (builds on existing `src/nmr/`)
2. **Track H2** — Vibrational Analysis & IR (builds on existing `src/ir/`)
3. **Track H3** — Charge Modeling & Atom Typing (builds on existing `src/charges/`, `src/smarts/`)
4. **Track H4** — Conformational Sampling (builds on existing `src/alignment/`, `src/conformer.rs`)
5. **Track H5** — Topology & Descriptors (builds on existing `src/ml/`, `src/topology.rs`)
6. **Track H9** — Stereochemistry (critical for correctness, smaller scope)
7. **Track H10** — Visualization (builds on existing `src/eht/marching_cubes.rs`)
8. **Track H6** — Molecular Dynamics (builds on existing `src/dynamics.rs`)
9. **Track H7** — Implicit Solvation (builds on existing `src/surface/`)
10. **Track H8** — SCF Evolution (builds on existing `src/hf/`, `src/pm3/`, `src/xtb/`)

---

## Validation Framework

All tests live in `tests/test_algorithm_validation.rs` and enforce:

```
Maximum allowed deviation: 0.5%
Metric: |computed - reference| / |reference| × 100 < 0.5%
```

Reference sources:
- **NMR shifts:** NMRShiftDB2, Bruker spectral database
- **IR frequencies:** NIST WebBook vibrational data
- **Charges:** OpenBabel Gasteiger output
- **RMSD/Clustering:** RDKit conformer ensemble reference
- **Ring perception:** IUPAC SSSR definition
- **Fingerprints:** RDKit Morgan fingerprint reference
- **MD energy conservation:** analytical harmonic oscillator
- **Solvation:** SMD/PCM literature values
- **SCF energies:** PySCF minimal-basis reference
- **Stereochemistry:** IUPAC CIP rules, textbook assignments
- **Mesh quality:** analytical sphere geometry
