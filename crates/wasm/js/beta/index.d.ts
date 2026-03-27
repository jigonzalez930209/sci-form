/**
 * sci-form-wasm/beta — TypeScript declarations.
 *
 * All functions accept/return JSON strings.
 * Call `await init()` from the main `sci-form-wasm` package once before use.
 */

// ── B1: KPM ──────────────────────────────────────────────────────────────

/**
 * Compute EHT Density of States using Kernel Polynomial Method (O(N)).
 *
 * KPM scales linearly with system size vs. O(N³) exact diagonalization.
 *
 * @param elements_json  JSON array of atomic numbers
 * @param coords_flat_json  JSON flat coords [x0,y0,z0,...]
 * @param config_json  Optional JSON config:
 *   `{ order?: number, n_points?: number, e_min?: number, e_max?: number, temperature?: number }`
 * @returns JSON:
 *   `{ energies: number[], total_dos: number[], e_min, e_max, order }`
 */
export declare function beta_compute_kpm_dos(
  elements_json: string,
  coords_flat_json: string,
  config_json?: string
): string;

/**
 * Compute Mulliken populations using KPM stochastic trace (O(N)).
 * @returns JSON `{ mulliken_charges: number[], orbital_populations: number[] }`
 */
export declare function beta_compute_kpm_mulliken(
  elements_json: string,
  coords_flat_json: string
): string;

// ── B2: MBH ──────────────────────────────────────────────────────────────

/**
 * Compute vibrational frequencies using Mobile Block Hessian.
 *
 * Reduces degrees of freedom by treating rigid groups as blocks,
 * giving speedup = n_dof_full / n_dof_reduced.
 *
 * @param smiles  SMILES string for ring detection and UFF topology
 * @returns JSON:
 *   `{ frequencies: number[], n_blocks, n_flexible,
 *      n_dof_reduced, n_dof_full, speedup }`
 */
export declare function beta_compute_mbh_frequencies(
  elements_json: string,
  coords_flat_json: string,
  smiles: string
): string;

// ── B3: RandNLA ───────────────────────────────────────────────────────────

/**
 * Solve the EHT generalized eigenvalue problem using randomized Nyström.
 *
 * Scales as O(N k²) instead of O(N³), where k ≈ √N by default.
 *
 * @param config_json  Optional JSON:
 *   `{ sketch_size?: number, seed?: number, max_error?: number, fallback_enabled?: boolean }`
 * @returns JSON:
 *   `{ orbital_energies, homo_index, homo_energy, lumo_energy, gap,
 *      k, residual_error, used_fallback }`
 */
export declare function beta_solve_eht_randnla(
  elements_json: string,
  coords_flat_json: string,
  config_json?: string
): string;

// ── B4: Riemannian ────────────────────────────────────────────────────────

/**
 * Compute Riemannian distance between two PSD matrices on the manifold.
 * @param matrix_a_json / matrix_b_json  JSON flat row-major arrays
 * @param dim  Matrix dimension (matrices are dim×dim)
 * @returns JSON `{ distance: number }`
 */
export declare function beta_psd_distance(
  matrix_a_json: string,
  matrix_b_json: string,
  dim: number
): string;

/**
 * Project a symmetric matrix onto the PSD cone (zero negative eigenvalues).
 * @param matrix_json  JSON flat row-major array
 * @param dim  Matrix dimension
 * @returns JSON `{ matrix: number[], was_psd: boolean }`
 */
export declare function beta_psd_projection(
  matrix_json: string,
  dim: number
): string;

// ── B5: CPM ───────────────────────────────────────────────────────────────

/**
 * Compute constant-potential electrochemical charges.
 * @param potential  Target electrode potential in eV
 * @returns JSON `{ charges: number[], total_charge: number, energy: number }`
 */
export declare function beta_compute_cpm_charges(
  elements_json: string,
  coords_flat_json: string,
  potential: number
): string;

// ── Metadata ─────────────────────────────────────────────────────────────

/**
 * Return available beta module list and stability info.
 * @returns JSON `{ beta_modules: string[], stability: string, version: string }`
 */
export declare function beta_modules_info(): string;
