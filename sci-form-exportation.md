> **Historical document** — This is a session exploration log, not current documentation. For the current API surface see `.github/copilot-instructions.md` and `README.md`.

# sci-form Exploration Summary

## Key Files Read Successfully

### 1. **XTB (GFN0-xTB Tight-Binding) Module**
- **File**: `src/xtb/solver.rs` and `src/xtb/mod.rs`
- **Result struct**: `XtbResult`
  - `orbital_energies: Vec<f64>` (eV, ascending)
  - `electronic_energy: f64` (eV)
  - `repulsive_energy: f64` (eV)
  - `total_energy: f64` (eV)
  - `n_basis: usize`, `n_electrons: usize`
  - `homo_energy: f64`, `lumo_energy: f64`, `gap: f64`
  - `mulliken_charges: Vec<f64>` (per-atom)
  - `scc_iterations: usize`, `converged: bool`

- **Function**: `solve_xtb(elements: &[u8], positions: &[[f64; 3]]) -> Result<XtbResult, String>`

- **Implementation Details**:
  - SCC (self-consistent charges) loop with damped charge mixing
  - Löwdin S^{-1/2} orthogonalization
  - Supports 25 elements (H–Rn, including transition metals)
  - Uses charge-shifted Hamiltonian within SCC iterations

### 2. **Extended Hückel Theory (EHT) Module**
- **File**: `src/eht/solver.rs` and `src/eht/mod.rs`
- **Result struct**: `EhtResult`
  - `energies: Vec<f64>` (sorted ascending, eV)
  - `coefficients: Vec<Vec<f64>>` (MO coefficient matrix)
  - `n_electrons: usize`
  - `homo_index: usize`, `lumo_index: usize`
  - `homo_energy: f64`, `lumo_energy: f64`, `gap: f64`
  - `support: EhtSupport` (metadata on capabilities)

- **Function**: `solve_eht(elements: &[u8], positions: &[[f64; 3]], k: Option<f64>) -> Result<EhtResult, String>`

- **Solver Algorithm**:
  1. Build basis (STO Gaussian expansions)
  2. Compute S (overlap) and H (Hamiltonian) matrices
  3. Solve generalized eigenproblem HC = SCE via Löwdin orthogonalization:
     - Diagonalize S → S^{-1/2}
     - Transform H' = S^{-1/2} H S^{-1/2}
     - Diagonalize H' → eigenvalues E, eigenvectors C'
     - Back-transform C = S^{-1/2} C'
  4. Sort by energy (ascending)

### 3. **PM3 (NDDO Semi-empirical) Module**
- **File**: `src/pm3/solver.rs` and `src/pm3/mod.rs`
- **Result struct**: `Pm3Result`
  - `orbital_energies: Vec<f64>` (eV)
  - `electronic_energy: f64` (eV)
  - `nuclear_repulsion: f64` (eV)
  - `total_energy: f64` (eV)
  - `heat_of_formation: f64` (**kcal/mol**)
  - `n_basis: usize`, `n_electrons: usize`
  - `homo_energy: f64`, `lumo_energy: f64`, `gap: f64`
  - `mulliken_charges: Vec<f64>`
  - `scf_iterations: usize`, `converged: bool`

- **Function**: `solve_pm3(elements: &[u8], positions: &[[f64; 3]]) -> Result<Pm3Result, String>`

- **Implementation**:
  - NDDO (neglect of diatomic differential overlap) approximation
  - Stewart's PM3 parameterization
  - Hartree-Fock SCF with damped mixing
  - One-center and two-center Coulomb integrals
  - PM3 core-core repulsion with exponential damping

### 4. **Density of States (DOS) Module**
- **Files**: `src/dos/dos.rs` and `src/dos/mod.rs`
- **Result struct**: `DosResult`
  - `energies: Vec<f64>` (eV grid)
  - `total_dos: Vec<f64>` (states/eV)
  - `pdos: Vec<Vec<f64>>` (per-atom projected DOS)
  - `sigma: f64` (Gaussian smearing width)

- **Functions**:
  - `compute_dos(energies, sigma, e_min, e_max, n_points) -> DosResult`
  - `compute_pdos(elements, positions, orbital_energies, coefficients, n_electrons, sigma, e_min, e_max, n_points) -> DosResult`
  - `compute_dos_parallel(...)` (rayon-based)
  - `compute_pdos_parallel(...)` (rayon-based)

- **Implementation**:
  - Gaussian smearing: `norm * exp(-(e - ei)² / (2σ²))`
  - PDOS computed via Mulliken-weighted orbital contributions
  - Parallel version uses rayon

### 5. **Dipole Moment Module**
- **Files**: `src/dipole/dipole.rs` and `src/dipole/mod.rs`
- **Result struct**: `DipoleResult`
  - `vector: [f64; 3]` (Debye)
  - `magnitude: f64` (Debye)
  - `unit: String` (always "Debye")

- **Functions**:
  - `compute_dipole(mulliken_charges: &[f64], positions: &[[f64; 3]]) -> DipoleResult`
  - `compute_dipole_from_eht(elements, positions, coefficients, n_electrons) -> DipoleResult`

- **Formula**: 
  - `μ = Σ_A q_A * R_A` (charge × position)
  - Conversion: 1 e·Å = 4.80321 Debye
  - Uses Mulliken partial charges

### 6. **Population Analysis Module**
- **Files**: `src/population/mod.rs` and `src/population/population.rs`
- **Result struct**: `PopulationResult`
  - `mulliken_charges: Vec<f64>` (per-atom)
  - `lowdin_charges: Vec<f64>` (per-atom)
  - `mulliken_populations: Vec<f64>` (per-AO gross populations)
  - `total_charge_mulliken: f64`, `total_charge_lowdin: f64`

- **Bond Order Result struct**: `BondOrderResult`
  - `bonds: Vec<BondOrderEntry>` (Wiberg & Mayer indices per pair)
  - `wiberg_valence: Vec<f64>` (sum for each atom)
  - `mayer_valence: Vec<f64>` (sum for each atom)

- **Functions**:
  - `mulliken_charges(elements, basis, overlap, coefficients, n_electrons) -> Vec<f64>`
  - `lowdin_charges(elements, basis, overlap, coefficients, n_electrons) -> Vec<f64>`
  - `compute_population(elements, positions, coefficients, n_electrons) -> PopulationResult`
  - `compute_bond_orders(elements, positions, coefficients, n_electrons) -> BondOrderResult`

- **Implementation**:
  - Mulliken: q_A = Z_A - Σ_{μ∈A} (PS)_{μμ}
  - Löwdin: q_A = Z_A - Σ_{μ∈A} (S^{1/2} P S^{1/2})_{μμ}
  - Uses density matrix P = Σ_i^{occ} n_i c_i c_i^T

### 7. **Gradient Computation (Force Field Module)**
- **File**: `src/forcefield/gradients.rs`
- **Functions** (all analytical):
  - `analytical_grad_bond(p1, p2, kb, r_eq, grad, idx1, idx2)`
  - `analytical_grad_angle(p1, pc, p2, k_theta, theta_eq, g1, gc, g2)`
  - `analytical_grad_bounds(p1, p2, lower, upper, k_bounds, grad, idx1, idx2)`
  - `analytical_grad_torsion(p1, p2, p3, p4, v, n_fold, gamma, grad, idx1, idx2, idx3, idx4)` 
  - `analytical_grad_oop(pc, p1, p2, p3, k_oop, grad, idx_c, idx1, idx2, idx3)` (out-of-plane)
  - `analytical_grad_distance_constraint(p1, p2, min_len, max_len, k, grad, idx1, idx2)`
  - `compute_analytical_gradient(coords, mol, params, bounds_matrix) -> DMatrix<f32>`

- **Trait**: `ForceFieldContribution`
  - `evaluate_energy_and_inject_gradient(coords: &[f64], grad: &mut [f64]) -> f64`

- **Force Field Wrapper**: `MolecularForceField`
  - `compute_system_energy_and_gradients(coords: &[f64], grad: &mut [f64]) -> f64`
  - Parallel support via rayon

### 8. **API Functions in lib.rs (> line 500)**
Key public functions:
- `embed(smiles: &str, seed: u64) -> ConformerResult`
- `embed_batch(smiles_list: &[&str], config: &ConformerConfig) -> Vec<ConformerResult>`
- `parse(smiles: &str) -> Result<Molecule, String>`
- `compute_charges(smiles: &str) -> Result<ChargeResult, String>`
- `compute_sasa(elements, coords_flat, probe_radius) -> Result<SasaResult, String>`
- `compute_population(elements, positions) -> Result<PopulationResult, String>`
- `compute_dipole(elements, positions) -> Result<DipoleResult, String>`
- `compute_esp(elements, positions, spacing, padding) -> Result<EspGrid, String>`
- `compute_dos(elements, positions, sigma, e_min, e_max, n_points) -> Result<DosResult, String>`
- `compute_rmsd(coords, reference) -> f64`
- `compute_uff_energy(smiles, coords) -> Result<f64, String>`
- `compute_mmff94_energy(smiles, coords) -> Result<f64, String>`
- `compute_pm3(elements, positions) -> Result<Pm3Result, String>`
- `compute_xtb(elements, positions) -> Result<XtbResult, String>`
- `compute_ml_descriptors(elements, bonds, charges, aromatic_atoms) -> MolecularDescriptors`
- `predict_ml_properties(desc) -> MlPropertyResult`
- `compute_md_trajectory(...) -> Result<MdTrajectory, String>`
- `compute_md_trajectory_nvt(...) -> Result<MdTrajectory, String>`
- `search_conformers_with_uff(smiles, n_samples, seed, rmsd_threshold) -> Result<ConformerSearchResult, String>`

## Tests Directory
Located at `/home/lestad/github/sci-form/tests/`
Key test files:
- `test_eht.rs` - EHT solver tests
- `test_charges.rs` - Gasteiger charge tests
- `test_sasa.rs` - SASA computation tests
- `test_10k.rs`, `test_chembl_10k.rs`, `test_gdb20.rs` - Large molecule sets
- `test_gradient_*` - Gradient verification tests
- `test_cross_language.py` - Python integration tests
- `test_wasm_integration.js` - WebAssembly tests
- `debug_*.rs` - Various debugging utilities

## Unit Conversion Summary
| Quantity | Unit |
|----------|------|
| Coordinates | Å |
| Energies (EHT, PM3, xTB) | eV |
| PM3 heat of formation | **kcal/mol** |
| Force field energy (UFF/MMFF94) | kcal/mol |
| Dipole | Debye |
| SASA | Ų |
| DOS grid spacing | Å |

## Key Observations for Numerical Hessian Research
1. **Gradient functions exist** for UFF force field but NOT for EHT/PM3/xTB directly
2. **Finite differences** would need to use energy functions (compute_eht, compute_pm3, compute_xtb)
3. **Force field gradients** are fully analytical in `forcefield/gradients.rs`
4. **Coordinate format**: Always flat `[x0,y0,z0, x1,y1,z1, ...]` as `Vec<f64>`
5. **No explicit frequency analysis** - would need to implement via numerical Hessian
6. **SCC/SCF convergence** varies: xTB (max 50 iter), PM3 (max 100 iter)
