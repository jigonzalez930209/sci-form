# Experimental Roadmap — Next-Generation Algorithms

> **Important notice:** All tasks in this document are **experimental** and are developed **in parallel without interfering** with the currently implemented and working production algorithms. Each experimental module will live in its own namespace (`sci_form::experimental::*`) or behind a feature flag during the research phase, ensuring that stable code is never affected.

---

## Track Summary

| Track | Description | Group |
|-------|-------------|-------|
| E1 | Conformal Geometric Algebra (CGA) | Advanced Geometry |
| E2 | RandNLA for EHT O(N log N) | Advanced Geometry |
| E3 | Riemannian Optimization for ETKDG | Advanced Geometry |
| E4 | Kernel Polynomial Method (KPM) O(N) | Electronic Structure |
| E5 | Dynamic EEQ Force Field | Electronic Structure |
| E6 | Analytical Implicit Solvation (ALPB) | Electronic Structure |
| E7 | D4 Dispersion Correction | Precision Methods |
| E8 | Semidefinite Relaxation Embedding (SDR) | Precision Methods |
| E9 | Mobile Block Hessian (MBH) | Precision Methods |
| E10 | Constant Potential Method (CPM) | Precision Methods |
| E11 | Growing String Method (GSM) | Precision Methods |

---

## Group A — Advanced Geometry and Algebra

---

## Track E1: Conformal Geometric Algebra (CGA)

**Goal:** Replace the internal geometry engine (separate quaternions + vectors) with Conformal Geometric Algebra G(4,1), unifying rotations and translations into a single multiplicative operator called a **Motor** ($M X \tilde{M}$). Eliminates gimbal lock and enables massively parallelizable operations for conformer generation and MOF assembly.

**Experimental module:** `sci_form::experimental::cga`
**Feature flag:** `experimental-cga`
**Does not modify:** `src/conformer/`, `src/materials/`, or any stable public API.

---

### E1.1 — CGA Algebraic Layer

#### E1.1a — G(4,1) Multivector Type
- Internal representation: 32-component array `[f64; 32]` (algebra basis 2^5)
- Implement `Multivector` with operators `*` (geometric product), `^` (outer), `|` (inner)
- Conformal basis sign table: $e_1^2 = e_2^2 = e_3^2 = e_+^2 = 1$, $e_-^2 = -1$
- Tests: verify $e_i e_j = -e_j e_i$ for $i \neq j$; $e_+ e_- + e_- e_+ = 0$

#### E1.1b — Fundamental Operations
- Full geometric product (32×32 multiplication table precomputed as `const`)
- Reverse $\tilde{M}$: reverse factor order in each blade
- Norm and normalization: $\|M\| = \sqrt{M \tilde{M}}$
- Grade involution: $\hat{M}$ (sign by blade grade)
- Tests: verify identity $M \tilde{M} = 1$ for unit motors

#### E1.1c — Motor (Rotor + Translator)
- Axis-angle rotor: $R = \cos(\theta/2) - \sin(\theta/2) \hat{L}$ where $\hat{L}$ is the line in CGA
- Translator: $T = 1 - \frac{1}{2} t e_\infty$ with translation vector $t$
- Composite motor: $M = T R$
- Sandwich application: $X' = M X \tilde{M}$ on a conformal point
- Tests: 90° rotation about Z axis, translation (1,0,0), composition of both

---

### E1.2 — CGA-based Conformer Generation

#### E1.2a — Point Embedding in CGA
- Lift atomic coordinates (`[f64;3]`) → conformal points: $P = p + \frac{1}{2}|p|^2 e_\infty + e_o$
- Extract coordinates from multivector after transformation
- Test suite: round-trip embed/extract with error < 1e-10 Å

#### E1.2b — Dihedral Torsion as CGA Motor
- Given torsional axis (bond A→B), build dihedral rotation Motor
- Identify subtree of atoms to rotate (same logic as current `find_rotatable_bonds`)
- Apply $X'_i = M X_i \tilde{M}$ to all subtree atoms in parallel (SIMD/rayon)
- Tests: resulting dihedral angle matches requested angle ±0.001°

#### E1.2c — Full CGA Refinement Branch
- Replace torsional rotation step in experimental `torsion_scan.rs` with CGA Motor
- Benchmark: wall-clock time vs current implementation on n-butane and aspirin
- Geometric validation: RMSD of CGA vs current coordinates < 0.001 Å
- Verify absence of gimbal lock at torsion ω ≈ ±90°

---

### E1.3 — CGA-based Materials Assembly (MOFs)

#### E1.3a — Unit Cell as CGA Frame
- Encode the three lattice vectors $(a, b, c)$ as translation motors
- Crystallographic frame: $M_{cell} = T_a T_b T_c$
- Fractional ↔ Cartesian conversion via inverse of the composite motor

#### E1.3b — SBU Placement via Motors
- Each Secondary Building Unit (SBU) is positioned with a Motor $M_{SBU} = T_i R_i$
- Motor composition for supercells: $M_{super} = M_{cell}^{n_a \cdot n_b \cdot n_c}$
- Framework construction: loop over topology, apply motor to node coordinates

#### E1.3c — Validation Against Current Assembler
- Compare stable `assemble_framework` vs CGA version for topology `pcu`
- Verify: node-linker bond distances identical ±0.001 Å
- Performance benchmark for 2×2×2 and 3×3×3 supercells

---

## Track E2: RandNLA for EHT O(N log N)

**Goal:** Replace the $O(N^3)$ EHT diagonalization with the **Randomized Nyström Approximation**, reducing cost to $O(N^2)$ or $O(N \log N)$ for sparse matrices, enabling systems of thousands of atoms.

**Experimental module:** `sci_form::experimental::rand_nla`
**Feature flag:** `experimental-randnla`
**Does not modify:** stable `src/eht/`.

---

### E2.1 — Nyström Approximation Infrastructure

#### E2.1a — Gaussian Sketch Matrix Generation
- Generate $\Omega \in \mathbb{R}^{N \times k}$ with entries $\mathcal{N}(0, 1/k)$
- Configurable $k$ (default: $k = \min(N, 50 + \lceil\log N\rceil)$)
- Alternative: SRHT (Subsampled Randomized Hadamard Transform) for $N > 5000$
- Tests: verify $E[\Omega^T \Omega] \approx I_k$ with tolerance $O(1/\sqrt{N})$

#### E2.1b — Product Y = SΩ and Factorization
- Compute $Y = S \Omega$ via sparse-dense multiplication if $S$ is sparse
- Orthogonalize $Y$: thin QR decomposition → $Y = Q R$
- Compute $B = Q^T S Q$ (small $k \times k$ matrix)
- Tests: verify $\|Y - Q(Q^TY)\|_F < 10^{-12}$

#### E2.1c — Reconstruction $S \approx Y(\Omega^T Y)^{-1} Y^T$
- Factor $\Omega^T Y$ with Cholesky (must be SPD), regularize if $\kappa > 10^{12}$
- Stabilized pseudoinverse over the $k$ dominant singular values
- Error bound: $\|S - \hat{S}\|_F \leq \|S - S_k\|_F (1 + \epsilon)$ with high probability
- Tests: Frobenius error < 1% for benzene overlap matrix $S$

---

### E2.2 — Randomized S^{-1/2} and Diagonalization

#### E2.2a — Low-Rank S^{-1/2} via rSVD
- Randomized SVD: $S \approx U \Sigma V^T$ with $O(Nk)$ flops
- $S^{-1/2} \approx U \Sigma^{-1/2} U^T$ (rank $k$, truncate singular values < $\epsilon_{mach}$)
- Stability comparison vs Cholesky + standard triangular solve

#### E2.2b — Randomized Diagonalization of $\tilde{H}$
- $\tilde{H} = S^{-1/2} H S^{-1/2}$: use low-rank versions of both factors
- Randomized eigensolver: sketch with $\Omega$, Rayleigh-Ritz refinement
- Extract only the $n_{occ}/2$ occupied eigenvectors (full spectrum not needed)

#### E2.2c — Error Estimation and Confidence Bound
- A posteriori error monitor: $\|H C - S C E\|_F / \|S C E\|_F$
- Automatic fallback: if error > 0.1%, escalate to exact diagonalization
- Confidence report in result: field `rand_nla_error: f64` in `EhtResult`

---

### E2.3 — Validation and Benchmarks

#### E2.3a — Accuracy: Benzene and Naphthalene
- HOMO/LUMO energies from RandNLA vs exact diagonalization: deviation < 0.1%
- HOMO-LUMO gap: error < 0.01 eV
- Mulliken charges: MAE < 0.001 e

#### E2.3b — Scaling: 100 and 1000 Atom Systems
- Nanotube (C100H20) and fullerene C60: measure wall time for $k = 20, 50, 100$
- Plot time vs $N$ to confirm sub-cubic scaling
- Crossover threshold: identify $N$ where RandNLA outperforms exact diagonalization

#### E2.3c — Memory Footprint and WASM Usage
- Measure peak RSS memory in standard vs RandNLA releases
- Verify WASM binary does not exceed 4 GB linear memory limit
- Node.js benchmark: `eht_calculate` latency on a 200-atom system

---

## Track E3: Riemannian Optimization for ETKDG

**Goal:** Replace the Euclidean BFGS optimizer in the embedding step with **Riemannian L-BFGS** over the manifold of fixed-rank PSD matrices, eliminating negative eigenvalues by design and reducing the retry loop failure rate to zero.

**Experimental module:** `sci_form::experimental::riemannian`
**Feature flag:** `experimental-riemannian`
**Does not modify:** stable `src/conformer/`.

---

### E3.1 — Riemannian Geometry Infrastructure

#### E3.1a — Retraction onto the PSD Cone
- Implement projection $P_+(X) = V \max(\Lambda, 0) V^T$ (zero out negative eigenvalues)
- First-order retraction: $\text{Retr}_X(\xi) = X + \xi + \frac{1}{2}\xi X^{-1} \xi$ (Burer-Monteiro)
- Manifold metrics: geodesic distance between two PSD matrices

#### E3.1b — Riemannian Gradient from Euclidean Gradient
- Project gradient to tangent space: $\text{grad}_R f = P_{T_X} \nabla f$
- For Stiefel manifold: $P_{T_X}(\xi) = \xi - X \text{sym}(X^T \xi)$
- Tests: verify $\langle \text{grad}_R f, \xi \rangle = D f(X)[\xi]$ for random perturbations

#### E3.1c — Exponential Map on Fixed-Rank Manifold
- $\text{Exp}_X(\xi) = (X + \xi + \frac{1}{2} \xi X^{-1} \xi)(I + X^{-1}\xi)^{-1}$
- Implement numerically stable version with regularization $X + \epsilon I$
- Parallel transport: move L-BFGS history along geodesic

---

### E3.2 — Riemannian L-BFGS for the Metric Matrix

#### E3.2a — L-BFGS History on the Manifold
- Store $m = 5$ pairs $(s_k, y_k)$ in tangent space
- Adapted two-loop recursion: transport $s_k, y_k$ to current point before recursion
- Riemannian curvature condition: $\langle y_k, s_k \rangle_R > 0$

#### E3.2b — Line Search on the Manifold
- Armijo along geodesic: $f(\text{Exp}_X(-\alpha d)) \leq f(X) - c_1 \alpha \|\text{grad} f\|^2$
- Strong Wolfe: also check curvature of transported gradient
- Fallback to steepest descent if L-BFGS fails Wolfe in 20 iterations

#### E3.2c — Convergence Criteria and Fallback
- Primary criterion: $\|\text{grad}_R f\| < \epsilon_{tol} = 10^{-6}$
- Secondary criterion: distance violation < 0.01 Å (per-pair maximum)
- Fallback to Euclidean BFGS + PSD projection if Riemannian does not converge in 500 iter

---

### E3.3 — Validation

#### E3.3a — Retry Loop Failure Rate
- Run 10k ChEMBL molecules with Riemannian embedding
- Metric: percentage of molecules requiring more than 1 attempt (target: < 0.1%)
- Compare against stable ETKDG baseline rate

#### E3.3b — Geometric Quality
- RMS deviation of bond lengths vs MMFF94 ideal values
- Bond angles: MAE vs ideal angles
- Compare distributions with violin plot on a 1000-molecule benchmark

#### E3.3c — Computational Performance
- Wall-clock time per molecule: Riemannian vs current BFGS
- Compare number of required energy evaluations
- Parallel transport overhead: measure as fraction of total time

---

## Group B — Electronic Structure and Force Field

---

## Track E4: Kernel Polynomial Method (KPM) — O(N)

**Goal:** Compute the EHT density matrix, DOS, and electronic populations via **Chebyshev polynomial expansion**, avoiding diagonalization entirely and achieving $O(N)$ scaling for sparse Hamiltonians.

**Experimental module:** `sci_form::experimental::kpm`
**Feature flag:** `experimental-kpm`

---

### E4.1 — Chebyshev Expansion Infrastructure

#### E4.1a — Spectral Radius Estimation
- Fast estimator of $[E_{min}, E_{max}]$ via Gershgorin or 10-step Lanczos
- Map to $[-1, 1]$: $\tilde{H} = (H - b I) / a$ with $a = (E_{max}-E_{min})/2$
- Tests: verify $\|\tilde{H}\|_2 \leq 1$ with high probability for 20 molecules

#### E4.1b — Chebyshev Recursion $T_k(\tilde{H})$
- Recurrence: $T_0 = I$, $T_1 = \tilde{H}$, $T_{k+1} = 2\tilde{H}T_k - T_{k-1}$
- Implement as sparse matvec: cost $O(N \cdot nnz/N) = O(nnz)$ per step
- Numerical stability: monitor $\max |T_k|_{ij}$ to detect divergence

#### E4.1c — Coefficients with Jackson Kernel
- Gibbs damping: $g_k^{(M)} = \frac{(M-k+1)\cos(\pi k/(M+1)) + \sin(\pi k/(M+1))\cot(\pi/(M+1))}{M+1}$
- Precompute $c_k = \int_{-1}^{1} f_{FD}(\epsilon) T_k(\epsilon) w(\epsilon) d\epsilon$ by quadrature
- Tests: 1D Hückel chain DOS vs analytical cosine solution

---

### E4.2 — Density, DOS and Populations via KPM

#### E4.2a — Fermi-Dirac Coefficients for Occupation Number
- Compute $\{c_k\}$ for Fermi-Dirac at electronic temperature $T_{el}$
- Electron count: $N_e = \text{Tr}[\rho] = \sum_k c_k \text{Tr}[T_k(H)]$ via stochastic traces
- Fermi level $\mu$ adjustment by bisection until $|N_e(\mu) - N_{target}| < 10^{-6}$

#### E4.2b — Stochastic Trace Estimation for DOS
- Stochastic trace: $\text{Tr}[A] \approx \frac{1}{n_v} \sum_{i=1}^{n_v} r_i^T A r_i$ with $n_v = 20$ random vectors
- Energy-resolved DOS: $g(\epsilon) \approx \frac{2}{N} \sum_k g_k c_k(\epsilon) \text{Tr}[T_k]$
- pDOS by element: project onto subset of atomic indices

#### E4.2c — Mulliken Charges from KPM Density Matrix
- $\rho_{ij} \approx \sum_k c_k (T_k)_{ij}$: evaluate only diagonal and needed entries
- Mulliken population: $q_i = Z_i - (PS)_{ii}$ (same formula as diagonalization, different $P$)
- Tests: KPM Mulliken charges vs exact in benzene, MAE < 0.01 e

---

### E4.3 — Validation and Scaling

#### E4.3a — DOS vs Exact Diagonalization
- Benzene (N=12), naphthalene (N=18), coronene (N=36): compare broadened DOS
- Metric: L2 distance between DOS curves in [-20, 5] eV, target < 1%
- Required expansion order: find $M^*$ achieving acceptable quality

#### E4.3b — 1000-Atom System: Graphene
- Graphene sheet (N=1000 atoms, H-terminated edges): wall time KPM vs O(N³)
- Measure speedup and memory footprint
- Verify DOS and semiconductor gap are correct

#### E4.3c — Population Analysis Accuracy
- 20 diverse organic molecules: KPM Mulliken vs exact Mulliken
- HOMO/LUMO energies: deviation < 0.05 eV
- Integration with public `compute_dos` via flag without breaking API

---

## Track E5: Dynamic EEQ Force Field

**Goal:** Implement a molecular mechanics engine with dynamic charges based on **Geometry-Dependent Electronegativity Equalization (EEQ)**, capturing geometry-dependent polarization and many-body dispersion. More accurate than MMFF94 for organometallic systems.

**Experimental module:** `sci_form::experimental::eeq`
**Feature flag:** `experimental-eeq`

---

### E5.1 — EEQ Charge Model

#### E5.1a — Per-Element Parameters
- For each element: electronegativity $\chi_i$, chemical hardness $\eta_i$, charge radius $r_i^{EEQ}$
- Source: Grimme GFN-FF parameters (permissive license) or independent rederivation
- Table: H, C, N, O, F, P, S, Cl, Br, I (10 key organic elements)

#### E5.1b — Fractional Coordination Number
- $CN_i = \sum_{j \neq i} \frac{1}{1 + e^{-16(r_{cov,ij}/r_{ij} - 1)}}$
- Analytical derivative $\partial CN_i / \partial r_{ij}$ for gradients
- Periodicity correction for crystal unit cells (lattice image summation)

#### E5.1c — Linear Charge System Solve
- $(N+1)\times(N+1)$ system $A q = b$: $A_{ii} = \eta_i$, $A_{ij} = \gamma(r_{ij})$, constraint $\sum q_i = Q_{total}$
- Solve with Cholesky (A is SPD for valid geometries) or regularized LDLT
- Tests: EEQ ethanol charges vs Gasteiger MAE < 0.05 e; charge neutrality guaranteed

---

### E5.2 — EEQ Energy and Gradients

#### E5.2a — EEQ Electrostatic Energy
- $E_{elec} = \sum_i \chi_i q_i + \frac{1}{2} \sum_{i,j} q_i \gamma(r_{ij}) q_j$
- Damping function $\gamma(r_{ij}) = \text{erf}(\sqrt{2}/\sigma_{ij} \cdot r_{ij})/r_{ij}$
- Analytical gradient $\partial E_{elec} / \partial R_k$ including charge response $\partial q_i/\partial R_k$

#### E5.2b — Many-Body D3/D4 Dispersion
- Integrate D3-BJ dispersion correction (optional 3-body term)
- CN-dependent $C_6^{ij}$ coefficients: bilinear interpolation from atomic references
- Full EEQ force field: $E_{total} = E_{rep} + E_{elec} + E_{disp} + E_{cov}$

#### E5.2c — Analytical Gradients and Numerical Hessian
- Total gradients $\nabla_R E$ for geometry minimization and MD
- Hessian: 2D finite differences over analytical gradients (flexible bonds only)
- Optimization convergence criterion: $\|\nabla E\| < 10^{-4}$ Hartree/Bohr

---

### E5.3 — Validation

#### E5.3a — Conformer Energies vs MMFF94
- 20 drug-like molecules: compare EEQ vs MMFF94 rotational energy profiles
- Correlation $R^2 > 0.90$ for rotational barriers
- Identify cases where EEQ improves on MMFF94 (C-O bond in alcohols, C-S in thioethers)

#### E5.3b — Charges vs Gasteiger on 20 Molecules
- Systematic comparison using EHT Mulliken charges as reference
- MAE EEQ vs Mulliken: target < 0.05 e; RMSE < 0.08 e
- Behavior on transition metals: Zn, Cu (not available in Gasteiger)

#### E5.3c — Geometry Optimization Convergence
- BFGS steps to minimize geometry: EEQ vs UFF vs MMFF94
- Final geometry: RMSD vs crystallography for 5 molecules with experimental structure

---

## Track E6: Analytical Implicit Solvation (ALPB)

**Goal:** Implement the Klamt/Grimme **Analytical Linearized Poisson-Boltzmann (ALPB)** model to compute $\Delta G_{solv}$ and its analytical gradients at negligible cost, enabling thermodynamically realistic conformers in solution.

**Experimental module:** `sci_form::experimental::alpb`
**Feature flag:** `experimental-alpb`
**Does not modify:** existing `src/solvation/` (stable GB).

---

### E6.1 — Born Radii and GB Kernel

#### E6.1a — Analytical Born Radii
- Hawkins-Cramer-Truhlar (HCT) integration enhanced with Onufriev-Bashford-Case
- $1/R_i^{eff} = 1/\rho_i - \tanh(\alpha \Psi_i - \beta \Psi_i^2 + \gamma \Psi_i^3)/\rho_i$
- Sensitivity to parameters $\alpha, \beta, \gamma$ tabulated per element

#### E6.1b — Function $f_{GB}(r_{ij}, R_i, R_j)$
- $f_{GB} = \sqrt{r_{ij}^2 + R_i R_j \exp(-r_{ij}^2 / (4 R_i R_j))}$
- Vectorized evaluation over all $O(N^2)$ pairs
- Correct asymptotic form: $f_{GB} \to r_{ij}$ as $r_{ij} \to \infty$

#### E6.1c — ALPB Correction on Top of GB
- ALPB correction term: $E_{ALPB} = E_{GB} \cdot \frac{\epsilon - 1}{\epsilon + A\cdot x}$ where $x = E_{GB}/E_{Coulomb}$
- Parameter $A = 0.571412$ (universal Klamt constant, no per-solvent fitting)
- Non-polar correction: $\Delta G_{np} = \gamma_{SA} \cdot SASA + \beta$ (reuse existing SASA)

---

### E6.2 — ALPB Energy and Gradients

#### E6.2a — Electrostatic Solvation Energy
- $\Delta G_{el} = -\frac{1}{2}(1 - 1/\epsilon) \sum_{ij} q_i q_j / f_{GB}$
- Configurable $\epsilon$: water (78.5), DMSO (47), chloroform (4.8), vacuum (1.0)
- Tests: $\Delta G_{solv}$ of Na⁺ ion in water vs experimental ≈ -98 kcal/mol

#### E6.2b — Analytical Gradient $\partial \Delta G_{el} / \partial R_k$
- Three contributions: direct derivative, derivative via $f_{GB}$, derivative via $R_i^{eff}(R_k)$
- Finite-difference verification of analytical gradients: error < 1e-6 kcal/(mol·Å)

#### E6.2c — Integration with Optimization Engine
- Wrapper `alpb_energy_and_gradient(elements, coords, charges, epsilon) -> (f64, Vec<f64>)`
- Compatible with minimize_geometry and MD backends
- Supplies gradients so Nosé-Hoover and Velocity Verlet can accept ALPB

---

### E6.3 — Validation

#### E6.3a — Comparison with Stable GB Model
- 10 neutral amino acids: $\Delta G_{solv}$ ALPB vs stable GB-HCT, difference < 5%
- Compute time: ALPB must be < 2x the cost of the existing GB

#### E6.3b — Correlation with Experimental Solubility
- Experimental $\log S$ from FreeSolv DB (20 small molecules) vs $\Delta G_{solv}$ ALPB
- Correlation $R^2 > 0.80$ target
- Compare slope with expected value (-1/RT ln 10 ≈ -0.73 units)

#### E6.3c — Conformer Re-ranking in Water vs Vacuum
- Aspirin (3 conformers), phenol (2 conformers): compare vacuum vs aqueous energy ranking
- Dominant conformer in water must match hydrated crystal structure

---

## Group C — Advanced Precision Methods

---

## Track E7: Grimme D4 Dispersion Correction

**Goal:** Implement the **DFT-D4** correction with geometry-dependent $C_6$/$C_8$ coefficients to augment EHT, UFF, and MMFF94 with quantum-precision van der Waals interactions.

**Experimental module:** `sci_form::experimental::d4`
**Feature flag:** `experimental-d4`

---

### E7.1 — D4 Parameter Infrastructure

#### E7.1a — Atomic Reference Polarizabilities
- Per element: table of polarizabilities $\alpha_0^{A,ref}$ for different reference states
- $C_6^{AA,ref}$ and $C_8^{AA,ref}$ coefficients from Grimme's D4 library (public data)
- Support: H, B, C, N, O, F, Si, P, S, Cl, Br, I + xTB transition metals

#### E7.1b — D4 Fractional Coordination Number
- $CN_i^{D4} = \sum_{j \neq i} f(r_{ij} / r_{cov,ij})$ with Fermi function
- Difference vs EEQ CN: D4 uses electronegativity-weighted covalent radii
- Analytical derivative for gradients: $\partial CN_i / \partial \mathbf{r}_j$

#### E7.1c — Dynamic $C_6$/$C_8$ Interpolation
- Mulliken partial charges as oxidation state proxy → interpolation weights
- $C_6^{AB} = \sum_{a \in A, b \in B} W_a W_b C_6^{AB,ref_{ab}}$ (generalized harmonic mean)
- Tests: C-C pair $C_6$ in methane, benzene, graphite: difference < 5% from CCSD(T) reference

---

### E7.2 — D4 Energy and Gradients

#### E7.2a — Two-Body Dispersion Sum
- $E_{disp}^{(2)} = -\sum_{A<B} \left[ s_6 \frac{C_6^{AB}}{r^6} f_{damp}^{(6)} + s_8 \frac{C_8^{AB}}{r^8} f_{damp}^{(8)} \right]$
- Becke-Johnson damping: $f_{damp}^{(n)} = r^n / (r^n + (a_1 \sqrt{C_8/C_6} + a_2)^n)$
- Method-specific $s_6, s_8, a_1, a_2$ (tables for EHT-D4, UFF-D4)

#### E7.2b — Three-Body Term (Axilrod-Teller-Muto)
- $E_{disp}^{(3)} = s_9 \sum_{A<B<C} C_9^{ABC} \frac{(3\cos\theta_A \cos\theta_B \cos\theta_C + 1)}{(r_{AB}r_{BC}r_{CA})^3}$
- Enable only for systems > 50 atoms (overhead justified)
- Tests: benzene trimer $E^{(3)}$ within 10% of CCSD(T)/aug-cc-pVDZ

#### E7.2c — Analytical D4 Gradients
- $\partial E_{disp} / \partial \mathbf{r}_A$: include $C_6, C_8$ response through $CN$
- Compact chain-rule form: $\partial C_6 / \partial \mathbf{r} = (\partial C_6/\partial CN)(\partial CN/\partial \mathbf{r})$
- Validation: $\|\nabla_{analytical} - \nabla_{fd}\| < 10^{-6}$ Hartree/Bohr

---

### E7.3 — Validation

#### E7.3a — S22 Non-Covalent Interaction Benchmark
- 22 reference dimers (Hobza S22): MAE of $\Delta E_{int}$ vs CCSD(T)/CBS
- Target: MAE < 0.3 kcal/mol for EHT-D4; < 0.5 kcal/mol for UFF-D4
- Breakdown: aromatic stacking, H-bonding, pure dispersion

#### E7.3b — Crystal Packing
- 5 simple organic crystals: D4 lattice energy vs experimental calorimetry
- Error < 5 kJ/mol for organic crystal lattice energies
- Compare π-π stacking distances: D4 vs experimental geometry

#### E7.3c — UFF-D4 and MMFF94-D4 Integration
- Wrapper `compute_uff_d4_energy(smiles, coords) -> f64` and MMFF94 analogue
- Benchmark: UFF-D4 vs plain UFF wall-clock time on 100 molecules
- Overhead target: < 15% additional time for N < 100 atoms

---

## Track E8: Semidefinite Relaxation Embedding (SDR)

**Goal:** Reformulate the conformer generation problem as a **convex optimization problem (SDP)**, mathematically guaranteeing that the distance matrix is positive semidefinite before extracting coordinates, eliminating the retry loop caused by negative eigenvalues.

**Experimental module:** `sci_form::experimental::sdr`
**Feature flag:** `experimental-sdr`

---

### E8.1 — SDP Formulation

#### E8.1a — Distance Constraints as Convex Set
- Decision variable: $X \in \mathbb{S}^N_+$ (PSD Gram matrix)
- Transformation: $d_{ij}^2 = X_{ii} + X_{jj} - 2X_{ij}$ (distances from Gram matrix)
- Projection $P_\Omega(X)$: enforce known bound distances onto entries of $X$

#### E8.1b — Alternating Projections Algorithm
- Projection 1: $P_+(X)$ — zero out negative eigenvalues (PSD subspace)
- Projection 2: $P_\Omega(X)$ — assign known distances (constraint subspace)
- Convergence: $\|X_{k+1} - X_k\|_F / \|X_k\|_F < 10^{-6}$, max 200 iterations

#### E8.1c — SVT (Singular Value Thresholding)
- SVT step: $\text{SVT}_\tau(X) = U \max(\Sigma - \tau, 0) V^T$
- Balance distance fidelity and low rank: $\min \|X\|_* + \lambda \|P_\Omega(X) - D\|_F^2$
- Adaptive $\tau$: reduce by factor 0.9 each improving iteration

---

### E8.2 — Integration with the Embedding Pipeline

#### E8.2a — Retry Loop Replacement
- Pre-process bounds distance matrix with SDR before eigendecomposition
- If SDR converges ($\|res\| < tol$): directly use $X^*$ to extract coordinates $C = V_3 \Sigma_3^{1/2}$
- If SDR does not converge in 200 iter: transparent fallback to current method with retry

#### E8.2b — Warm-Start from Cayley-Menger
- Initialize $X_0$ from the Cayley-Menger distance matrix (initial approximation)
- Compare convergence rate of warm-start vs cold start (zeros)
- Target: convergence in < 50 iterations with warm-start vs 150 without

#### E8.2c — Positive Semidefiniteness Guarantee
- Post-SDR verification: all eigenvalues of $X^*$ ≥ -1e-8 (numerical rounding)
- Geometric error bound: $\max_{ij} |d_{ij}^{SDR} - d_{ij}^{target}| < 0.1$ Å
- Log of eliminated retries: `sdr_attempts_saved: usize` in result metadata

---

### E8.3 — Validation

#### E8.3a — Zero Failure Rate on 10k Benchmark
- Run SDR on `data/chembl_10k_largest.smi`
- Metric: 0 molecules with negative eigenvalues in $X^*$ after SDR
- Compare against current stable pipeline retry rate (baseline)

#### E8.3b — Geometric Quality vs Stable ETKDG
- 1000 molecules: RMSD of SDR coordinates vs ETKDG after UFF minimization
- Torsion angles: compare Ramachandran-style clustering distributions
- Bond lengths and angles: equal or better than current method

#### E8.3c — Performance Overhead
- SDR must add < 2x time per molecule vs direct embedding (no retry)
- Measure net gain: (time_with_current_retry) vs (time_SDR_always_first_try)
- Parallelization: SDR is trivially parallelizable across molecules (no dependencies)

---

## Track E9: Mobile Block Hessian (MBH)

**Goal:** Reduce vibrational analysis cost from $O(3N)$ energy evaluations to $O(n_{flex})$ by treating rigid groups (aromatic rings, clusters) as bodies with 6 DOF, enabling real-time IR spectra for large molecules.

**Experimental module:** `sci_form::experimental::mbh`
**Feature flag:** `experimental-mbh`

---

### E9.1 — Rigid Block Detection

#### E9.1a — Ring Identification via SSSR
- Use stable `compute_sssr` to find rings ≤ 8 atoms
- Rigidity criterion: aromatic=true, or saturated ring without heteroatoms with substituents
- Classify: rigid, semi-rigid (1 flex bond), flexible

#### E9.1b — Molecular Graph Decomposition into Fragments
- Build fragment graph: each rigid block is a node
- Bridging atoms (between fragments): belong to no block, full DOF
- Mapping atom_index → (block_id or ∅) for the whole molecule

#### E9.1c — DOF Assignment per Block
- Each rigid block: 6 collective DOF (3 translation + 3 rotation of the block)
- Block DOF basis: 6 collective displacement vectors for the entire block
- Reduced system size: $n_{DOF} = 6 \cdot n_{blocks} + 3 \cdot n_{atoms_{flex}}$

---

### E9.2 — MBH Construction

#### E9.2a — Low-Level Hessian: Flexible DOF Only
- 2D finite differences over reduced DOF (not $3N$ Cartesians)
- Block displacement: move the entire block rigidly by ±h Å in each coordinate
- Savings: if 60% of atoms are in rigid blocks, Hessian is 40% the full size

#### E9.2b — MBH Projection and Assembly
- $H_{MBH} = L^T H_{full} L$ where $L$ is the projection matrix reduced DOF → Cartesian
- Columns of $L$: normalized collective displacement vectors per block (translations and rotations)
- Ensure $L^T L = I$ (orthogonality of the MBH basis)

#### E9.2c — Frequency Recovery via Generalized Eigensolution
- Solve $H_{MBH} C = M_{MBH} C \Lambda$ where $M_{MBH} = L^T M_{mass} L$
- Convert eigenvalues to $cm^{-1}$: $\nu_i = \frac{1}{2\pi c} \sqrt{\lambda_i}$
- Project MBH modes back to Cartesian displacements: $q^{cart} = L \cdot q^{MBH}$

---

### E9.3 — Validation

#### E9.3a — Benzene Frequencies: MBH vs Full Hessian
- Use UFF/MMFF94 as energy backend
- MBH vs full frequencies: target MAE < 5 cm⁻¹ for modes > 200 cm⁻¹
- Low-frequency modes (ring deformation < 200 cm⁻¹): extended tolerance of 20 cm⁻¹

#### E9.3b — Speedup on Penicillin (N=33 heavy atoms, 2 rings)
- Count energy evaluations: MBH vs full
- Expected speedup: ~3-5x for a molecule with 60% atoms in rings
- Validate broadened IR spectrum: intensity correlation

#### E9.3c — Normal Mode Shapes
- Compare MBH vs full eigenvectors: cosine similarity > 0.98 for the 30 highest-frequency modes
- Mode visualization: output as multi-frame XYZ trajectory (animation)

---

## Track E10: Constant Potential Method (CPM)

**Goal:** Extend charge and electronic property calculations by coupling the molecule to a "virtual electrode" with electrochemical potential $\mu$, allowing charges to fluctuate and enabling simulation of quantum capacitance and impedance.

**Experimental module:** `sci_form::experimental::cpm`
**Feature flag:** `experimental-cpm`

---

### E10.1 — Grand Potential Framework

#### E10.1a — Fermi Level $\mu$ as External Parameter
- Interface: `compute_cpm_charges(elements, coords, mu_ev: f64, epsilon: f64) -> CpmResult`
- $\mu$ in eV vs vacuum; relation to SHE: $\mu_{SHE} \approx -4.44$ eV
- Typical range: [-5.5, -3.5] eV (corresponds to -1 V to +1 V vs SHE)

#### E10.1b — Grand Potential $\Omega$ Functional
- $\Omega(\{q_i\}, \mu) = E_{elec}(\{q_i\}) - \mu \sum_i q_i$
- Minimize with respect to $\{q_i\}$ under constraint $\sum q_i = Q$ (free, not fixed)
- $Q^* = -\partial \Omega / \partial \mu |_{min}$ → total charge at electrochemical equilibrium

#### E10.1c — Charge Fluctuation During Optimization
- At each geometry optimization step: solve EEQ with fixed $\mu$, obtain $\{q_i(\mathbf{R})\}$
- Update gradient: include $\partial E / \partial q_i \cdot \partial q_i / \partial \mathbf{R}$ (charge response)
- Simultaneous geometric (force) and electronic (charge) convergence: unified criterion

---

### E10.2 — CPM Integration

#### E10.2a — Coupling with EEQ Model
- Replace charge neutrality constraint in EEQ with $\mu$ constraint
- Modify EEQ linear system: add column/row for electrode potential
- Tests: total charge vs $\mu$ for naphthalene in water should be monotonically decreasing

#### E10.2b — Electrochemical Energy Surface
- Scan $\mu$ in [-2, 2] V range: compute $\Omega(\mu)$, $Q(\mu)$, $\partial Q/\partial \mu$ (capacitance)
- Output: `CpmSurface { mu_values, total_charge, free_energy, capacitance }`
- $Q(\mu)$ curve: S-sigmoid shape for simple redox molecules

#### E10.2c — Conformer Ranking with Variable Potential
- For 3+ conformers of a redox molecule: compare $\Omega(\mu)$ vs $\mu$
- Dominant conformer crossover: identify $\mu^*$ where the minimum changes
- Application: anthraquinone (2 conformers): verify conformer switch at -0.6 V vs SHE

---

### E10.3 — Validation

#### E10.3a — Capacitance Curve for a Simple Redox Molecule
- Ferrocene / viologen: $C(\mu) = \partial Q / \partial \mu$ vs voltage curve
- Compare with experimental EIS data: capacitance peak position ±0.2 V
- Gaussian peak width: related to $k_BT/e$ at room temperature

#### E10.3b — Charge Response vs Potential
- Verify $\partial Q / \partial \mu > 0$ across the full range (capacitance always positive)
- Weak-field limit: $Q(\mu) \approx C_{geom} \cdot \mu$ recovers geometric capacitance
- Thermodynamic consistency test: $\int d\mu Q(\mu) = \Delta \Omega$ (Maxwell identity)

#### E10.3c — Comparison with Electrochemical DFT
- 3 molecules with published DFT-U data: $\Delta G_{redox}$ bias CPM vs DFT
- CPM vs DFT error: target < 0.3 eV for ionization energy in solution
- Computational cost: CPM must be > 100× faster than equivalent DFT/PCM

---

## Track E11: Growing String Method (GSM)

**Goal:** Implement a reaction path and transition state finder driven purely by force-field gradients, without requiring an initial guess of the full path as NEB does, enabling activation barrier prediction within the library.

**Experimental module:** `sci_form::experimental::gsm`
**Feature flag:** `experimental-gsm`

---

### E11.1 — Core GSM Algorithm

#### E11.1a — Reactant/Product String Initialization
- Input: reactant `R` and product `P` coordinates (minimized 3D geometries)
- Internal representation: redundant internal coordinates (generalized Z-matrix)
- Initialize each string: 3 nodes (R, interpolated midpoint, P)
- Data structure: `GsmPath { nodes: Vec<Vec<f64>>, energies: Vec<f64>, gradients: Vec<Vec<f64>> }`

#### E11.1b — Growth Step: Tangent Direction and Projection
- String tangent at interior node: normalized central difference
- Perpendicular gradient: $g^\perp = g - (g \cdot \hat{\tau})\hat{\tau}$
- Add new node at fixed distance $\Delta s = 0.1$ Å along the downhill gradient
- Growth stopping criterion: $|g^\perp| < 0.1$ kcal/(mol·Å) at all nodes

#### E11.1c — Junction Criterion and Saddle Point Detection
- Junction: when distance between string tips < $\Delta s / 2$
- Saddle point: node with $g^\perp \approx 0$ and zero tangential component (1D profile maximum)
- Refinement: 20 Dimer Method or P-RFO steps around the saddle candidate

---

### E11.2 — Energy Backend Integration

#### E11.2a — UFF/MMFF94 Backend for Fast GSM
- Use `compute_uff_energy` and numerical gradients (finite differences) as oracle
- Force field parameters fixed throughout the trajectory (no re-typing)
- Target time: transition state of a 20-atom reaction in < 5 seconds

#### E11.2b — EHT Backend for Electronic Effects
- Use `eht_calculate` as energy oracle (includes conjugation effects)
- Relevant for reactions in π-conjugated systems (cyclization, ring opening)
- EHT gradient: numerical with step h=0.001 Å (consistent with EHT Hessian)

#### E11.2c — Transition State Output
- `GsmResult { ts_coords: Vec<f64>, ts_energy: f64, activation_energy: f64, path_energies: Vec<f64>, path_coords: Vec<Vec<f64>>, irc_forward: Vec<Vec<f64>>, irc_backward: Vec<Vec<f64>> }`
- Topological TS validation: must have exactly 1 imaginary frequency (MBH Hessian)
- Export as multi-frame XYZ for reaction path visualization

---

### E11.3 — Validation

#### E11.3a — Hydrogen Transfer: Barrier vs NEB
- H-shift reaction in malonaldehyde (standard benchmark system)
- GSM-UFF vs NEB-UFF barrier: difference < 5%
- Energy evaluation count: GSM must require < 50% of a converged NEB

#### E11.3b — Ring Opening: Cyclobutane → 2 Ethylenes
- Concerted retro-[2+2] path: GSM must find the C2h saddle point
- Verify TS geometry has the correct imaginary normal mode (square collapse)
- Compare $\Delta E^\ddagger$ vs experimental kinetic data (62 kcal/mol thermolytic UFF)

#### E11.3c — Performance vs NEB
- 5 diverse reactions (H-transfer, ring opening, model SN2 substitution, cis-trans isomerization, tautomerization)
- Compare: convergence steps, energy evaluations, TS quality
- GSM target: < 300 energy evaluations per reaction with N < 30 atoms

---

## Implementation Guide and Priorities

### Recommended Implementation Order

| Priority | Track | Reason |
|----------|-------|--------|
| High | **E7** (D4) | Minimal new infrastructure; immediate impact on energy quality |
| High | **E4** (KPM) | Reuses existing EHT; O(N) is a key differentiator |
| High | **E8** (SDR) | Eliminates the retry loop; tangible pipeline improvement |
| Medium | **E6** (ALPB) | Extends existing GB; high user demand |
| Medium | **E5** (EEQ) | Prerequisite for CPM; improves charges on metals |
| Medium | **E9** (MBH) | Extends existing Hessian; real-time IR |
| Low | **E1** (CGA) | Fundamental research; high long-term impact |
| Low | **E2** (RandNLA) | Requires molecules > 500 atoms to justify the overhead |
| Low | **E3** (Riemannian) | Complex; E8 already solves the same problem |
| Research | **E10** (CPM) | Electrochemical niche; requires complete E5 |
| Research | **E11** (GSM) | High algorithmic complexity; high impact in catalysis |

### Development Conventions

1. **Never in the same module** as stable code. Path: `src/experimental/<track>/`
2. **Feature flags** enable modules: `cargo build --features experimental-d4`
3. **No breaking changes**: current public APIs must not be touched under any circumstances
4. **Isolated tests**: each track has its own file `tests/experimental/test_<track>.rs`
5. **Benchmarks**: comparative results in `benches/experimental/<track>_bench.rs`
6. **Graduation to production**: a track may leave experimental once it passes all its E_.3 validation tests and a full code review

---

*This document is the living research roadmap of sci-form. The tasks described here represent the state of the art in computational chemistry and must not be merged into production code until each track has completed its full validation phase.*
