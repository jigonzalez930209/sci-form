/**
 * sci-form-wasm/alpha — TypeScript declarations.
 *
 * All functions accept/return JSON strings except where noted.
 * Call `await init()` from the main `sci-form-wasm` package once before use.
 *
 * Signatures mirror `crates/wasm/src/alpha.rs` and
 * `crates/wasm/src/dynamics_live.rs` exactly: every parameter declared by the
 * Rust `#[wasm_bindgen]` function is required here, including the JSON
 * payload arguments (pass `"{}"` or `""` for the engine defaults).
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
  method: string
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
  params_json: string
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
 * Plan which GSM backends are available for a molecule.
 * @returns JSON GsmBackendPlan
 */
export declare function alpha_gsm_backend_plan(smiles: string): string;

/**
 * Compare one or more GSM backends on the same geometry.
 * @param methods_json JSON array of backend names, or `"[]"` for all
 * @returns JSON GsmBackendComparison
 */
export declare function alpha_gsm_compare_backends(
  smiles: string,
  coords_flat_json: string,
  methods_json: string
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
  config_json: string
): string;

/**
 * Find transition state using a specific GSM energy backend.
 * @param method Backend name, e.g. `"uff"`, `"pm3"`, `"gfn2"`, `"hf3c"`.
 * @returns JSON `GsmResult`
 */
export declare function alpha_gsm_find_ts_with_method(
  smiles: string,
  reactant_coords_json: string,
  product_coords_json: string,
  config_json: string,
  method: string
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
  config_json: string
): string;

// ── A8: Live MD ───────────────────────────────────────────────────────────

/** Live interactive molecular dynamics simulation handle. */
export declare class LiveSimulation {
  /** Create a simulation from a SMILES, flat coords and a backend name. */
  constructor(smiles: string, coords_flat: string, backend: string);
  /** Advance the simulation by `n_steps` with a `dt_fs` femtosecond step. */
  step(n_steps: number, dt_fs: number): void;
  /** Apply a Berendsen thermostat coupling towards `target_temp` K. */
  apply_thermostat(target_temp: number, tau_fs: number, dt_fs: number): void;
  /**
   * Apply a harmonic steering force pulling `atom_index` towards
   * `(tx, ty, tz)` with a force constant of `spring_k`.
   */
  apply_steering_force(
    atom_index: number,
    tx: number,
    ty: number,
    tz: number,
    spring_k: number
  ): void;
  /** Clear the accumulated steering forces (call before the next IMD frame). */
  clear_steering(): void;
  /** Number of atoms in the simulation. */
  num_atoms(): number;
  /** Current flat positions as a JSON array `[x0, y0, z0, ...]`. */
  positions_json(): string;
  /**
   * Pointer to the flat position buffer inside WASM memory; build a
   * `Float64Array` view of length `num_atoms() * 3` over it.
   */
  positions_ptr(): number;
  /** Current kinetic temperature in Kelvin. */
  temperature(): number;
  /** Current potential energy (from the last force computation). */
  potential_energy(): number;
  /** Current kinetic energy. */
  kinetic_energy(): number;
  /** Number of integration steps taken so far. */
  step_count(): number;
  /**
   * Reset the velocities from a Maxwell-Boltzmann distribution at
   * `temperature` K.
   *
   * `seed` is a **64-bit** seed, so it is a JavaScript `BigInt` (e.g. `42n`).
   */
  set_temperature(temperature: number, seed: bigint): void;
}

// ── A9: Reaction dynamics (3D) ────────────────────────────────────────────

/**
 * Compute a full 3D reaction-dynamics path (approach + NEB + departure).
 * @param reactant_smiles_json JSON array of reactant SMILES
 * @param product_smiles_json JSON array of product SMILES
 * @param config_json JSON `ReactionDynamics3DConfig`, `""` or `"{}"` for defaults
 * @returns JSON `ReactionDynamics3DResult`
 */
export declare function compute_reaction_dynamics_3d(
  reactant_smiles_json: string,
  product_smiles_json: string,
  config_json: string
): string;

/**
 * Compute the approach phase only (no products required).
 * @param reactant_smiles_json JSON array of reactant SMILES
 * @param config_json JSON `ReactionDynamics3DConfig`, `""` or `"{}"` for defaults
 * @returns JSON `ReactionDynamics3DResult` with `phase = "approach"`
 */
export declare function compute_approach_dynamics(
  reactant_smiles_json: string,
  config_json: string
): string;

// ── Render bridge: chart utilities ────────────────────────────────────────

/** Build an EDL profile chart from an `EdlProfileResult` JSON payload. */
export declare function alpha_render_edl_profile_chart(profile_json: string): string;

/** Build a capacitance-scan chart from a `CpmEdlScanResult` JSON payload. */
export declare function alpha_render_capacitance_scan_chart(scan_json: string): string;

/** Build a periodic-DOS chart from a `PeriodicKpmDosResult` JSON payload. */
export declare function alpha_render_periodic_dos_chart(dos_json: string): string;

/** Build a band-structure chart from a `BandStructureAdapterResult` JSON payload. */
export declare function alpha_render_band_structure_chart(bs_json: string): string;

/** Build a microkinetic population chart from a `MicrokineticTrace` JSON payload. */
export declare function alpha_render_microkinetic_population_chart(
  trace_json: string
): string;

/** Build an Arrhenius chart from a JSON array of `ElementaryRateResult`. */
export declare function alpha_render_arrhenius_chart(results_json: string): string;

/** Build GSM path trajectory frames from coords + energies JSON payloads. */
export declare function alpha_render_gsm_path_trajectory(
  coords_json: string,
  energies_json: string
): string;

/** Pack a chart payload into the columnar transport format. */
export declare function alpha_render_pack_chart_payload(chart_json: string): string;

/** Verify a chart payload survives a pack/unpack round trip. */
export declare function alpha_render_verify_chart_round_trip(
  chart_json: string,
  tol: number
): boolean;

// ── Periodic linear (k-meshes, BZ integration, EHT band structure) ─────────

/** Build a Monkhorst-Pack (or gamma-centred) k-mesh. */
export declare function alpha_periodic_monkhorst_pack(config_json: string): string;

/** Bloch phase for a fractional k-point and a lattice translation. */
export declare function alpha_periodic_bloch_phase(
  k_fractional_json: string,
  translation_json: string
): string;

/** Brillouin-zone weighted integral of a scalar field. */
export declare function alpha_periodic_bz_integrate_scalar(
  values_json: string,
  weights_json: string
): number;

/** Brillouin-zone weighted integral of a vector field. */
export declare function alpha_periodic_bz_integrate_vector(
  vectors_json: string,
  weights_json: string
): string;

/** Validate that the k-mesh weights sum to one. */
export declare function alpha_periodic_validate_kmesh_weights(
  kmesh_json: string,
  tol: number
): boolean;

/** EHT band structure along a k-path. */
export declare function alpha_periodic_band_structure_eht(
  elements_json: string,
  coords_flat_json: string,
  lattice_json: string,
  config_json: string,
  n_electrons: number
): string;

/** Validate the electron count implied by a periodic DOS result. */
export declare function alpha_periodic_validate_electron_count(
  dos_json: string,
  fermi_energy_ev: number,
  expected_electrons: number
): number;

// ── Metadata ─────────────────────────────────────────────────────────────

/**
 * Return available alpha module list and stability info.
 * @returns JSON `{ alpha_modules: string[], stability: string, version: string }`
 */
export declare function alpha_modules_info(): string;
