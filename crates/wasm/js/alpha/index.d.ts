/**
 * sci-form-wasm/alpha — TypeScript declarations.
 *
 * All functions accept/return JSON strings except where noted.
 * Call `await init()` from the main `sci-form-wasm` package once before use.
 */

// ── A1: Kohn-Sham DFT ─────────────────────────────────────────────────────

/**
 * Run a Kohn-Sham DFT single-point calculation.
 * @param elements_json  JSON array of atomic numbers e.g. `"[6,1,1,1,1]"`
 * @param coords_flat_json  JSON flat coords `"[x0,y0,z0,...]"`
 * @param method  `"svwn"` (LDA) | `"pbe"` (GGA, default)
 * @returns JSON DftResult:
 *   `{ energy, homo_energy, lumo_energy, gap, converged, n_basis,
 *      scf_iterations, mulliken_charges, orbital_energies,
 *      nuclear_repulsion, xc_energy }`
 */
export declare function alpha_compute_dft(
  elements_json: string,
  coords_flat_json: string,
  method?: string
): string;

// ── A2: ReaxFF ────────────────────────────────────────────────────────────

/**
 * Compute ReaxFF energy and analytical gradient.
 * @returns JSON `{ energy_kcal_mol: number, gradient: number[] }`
 */
export declare function alpha_compute_reaxff_gradient(
  elements_json: string,
  coords_flat_json: string
): string;

/**
 * Compute ReaxFF energy components (no gradient).
 * @returns JSON `{ bonded, coulomb, van_der_waals, total }` in kcal/mol
 */
export declare function alpha_compute_reaxff_energy(
  elements_json: string,
  coords_flat_json: string
): string;

/**
 * Compute EEM equilibrium point charges.
 * @returns JSON `{ charges: number[], total_charge: number }`
 */
export declare function alpha_compute_eem_charges(
  elements_json: string,
  coords_flat_json: string
): string;

// ── A3: MLFF ─────────────────────────────────────────────────────────────

/**
 * Compute MLFF energy and forces using element-specific neural networks.
 * @param config_json JSON `MlffConfig` object
 * @returns JSON `{ energy: number, atomic_energies: number[], forces: number[] }`
 */
export declare function alpha_compute_mlff(
  elements_json: string,
  coords_flat_json: string,
  config_json: string
): string;

/**
 * Compute Atomic Environment Vectors (AEVs) for each atom.
 * @param params_json JSON `SymmetryFunctionParams` or `"{}"`  for defaults
 * @returns JSON `{ n_atoms, aev_length, aevs: number[][] }`
 */
export declare function alpha_compute_aevs(
  elements_json: string,
  coords_flat_json: string,
  params_json?: string
): string;

// ── A4: Obara-Saika ERIs ─────────────────────────────────────────────────

/**
 * Compute Boys function F_n(x).
 * @returns JSON `{ value: number }`
 */
export declare function alpha_boys_function(n: number, x: number): string;

/**
 * Compute a (ss|ss) two-electron repulsion integral.
 * Each shell JSON: `{ "alpha": number, "center": [x, y, z] }`
 * @returns JSON `{ eri: number }`
 */
export declare function alpha_eri_ssss(
  shell_a_json: string,
  shell_b_json: string,
  shell_c_json: string,
  shell_d_json: string
): string;

/**
 * Compute Schwarz pre-screening upper bound Q_ab = √⟨ab|ab⟩.
 * @returns JSON `{ schwarz_bound: number }`
 */
export declare function alpha_schwarz_bound(
  alpha: number,
  center_a_json: string,
  beta: number,
  center_b_json: string
): string;

// ── A5: CGA Motor Algebra ────────────────────────────────────────────────

/**
 * Rotate a sub-tree of atoms around a bond using a CGA motor.
 * @param subtree_indices_json JSON array of atom indices to rotate
 * @param axis_a_json / axis_b_json `[x, y, z]` of axis endpoints
 * @param angle_rad Rotation angle in radians
 * @returns JSON `{ coords: number[] }`
 */
export declare function alpha_rotate_dihedral_cga(
  coords_flat_json: string,
  subtree_indices_json: string,
  axis_a_json: string,
  axis_b_json: string,
  angle_rad: number
): string;

/**
 * Refine a dihedral to a target angle using CGA.
 * @param torsion_indices_json JSON `[i, j, k, l]`
 * @returns JSON `{ coords: number[] }`
 */
export declare function alpha_refine_torsion_cga(
  elements_json: string,
  coords_flat_json: string,
  smiles: string,
  torsion_indices_json: string,
  target_angle_rad: number
): string;

// ── A6: Growing String Method ─────────────────────────────────────────────

/**
 * Interpolate a reaction path node at t ∈ [0, 1].
 * @returns JSON `{ coords: number[] }`
 */
export declare function alpha_gsm_interpolate(
  reactant_json: string,
  product_json: string,
  t: number
): string;

/**
 * Find transition state using Growing String Method with UFF energy.
 * @param config_json JSON `GsmConfig` or `"{}"`
 * @returns JSON `GsmResult`:
 *   `{ ts_coords, ts_energy, path_energies, n_nodes, converged }`
 */
export declare function alpha_gsm_find_ts(
  smiles: string,
  reactant_coords_json: string,
  product_coords_json: string,
  config_json?: string
): string;

// ── A7: SDR Embedding ────────────────────────────────────────────────────

/**
 * Embed a molecule from pairwise distances using Semidefinite Relaxation.
 * @param distance_pairs_json JSON `[[i, j, d_ij], ...]`
 * @param n_atoms Total atom count
 * @param config_json JSON `SdrConfig` or `"{}"`
 * @returns JSON `{ coords, residual, converged, iterations }`
 */
export declare function alpha_sdr_embed(
  distance_pairs_json: string,
  n_atoms: number,
  config_json?: string
): string;

// ── A8: Live MD ───────────────────────────────────────────────────────────

/** Live interactive molecular dynamics simulation handle. */
export declare class LiveSimulation {
  constructor(smiles: string, coords_flat: string, backend: string);
  step(n_steps: number, dt_fs: number): void;
  apply_thermostat(target_temp: number, coupling: number): void;
  apply_steering_force(atom_index: number, target: [number, number, number], k: number): void;
  clear_steering(): void;
  positions_json(): string;
  positions_ptr(): number;
  temperature(): number;
  potential_energy(): number;
  kinetic_energy(): number;
  step_count(): number;
  set_temperature(temp: number): void;
}

// ── Metadata ─────────────────────────────────────────────────────────────

/**
 * Return available alpha module list and stability info.
 * @returns JSON `{ alpha_modules: string[], stability: string, version: string }`
 */
export declare function alpha_modules_info(): string;
