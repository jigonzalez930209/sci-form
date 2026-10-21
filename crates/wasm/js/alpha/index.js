/**
 * sci-form-wasm/alpha — Alpha-tier experimental modules.
 *
 * Early-stage research implementations. APIs subject to breaking changes.
 *
 * Usage:
 * ```js
 * import init from 'sci-form-wasm';
 * import {
 *   alpha_compute_dft,
 *   alpha_compute_reaxff_gradient,
 *   alpha_compute_reaxff_energy,
 *   alpha_compute_eem_charges,
 *   alpha_compute_mlff,
 *   alpha_compute_aevs,
 *   alpha_boys_function,
 *   alpha_eri_ssss,
 *   alpha_schwarz_bound,
 *   alpha_rotate_dihedral_cga,
 *   alpha_refine_torsion_cga,
 *   alpha_gsm_interpolate,
 *   alpha_gsm_backend_plan,
 *   alpha_gsm_compare_backends,
 *   alpha_gsm_find_ts,
 *   alpha_gsm_find_ts_with_method,
 *   alpha_sdr_embed,
 *   alpha_modules_info,
 *   LiveSimulation,
 * } from 'sci-form-wasm/alpha';
 *
 * await init();
 * const result = JSON.parse(alpha_compute_dft(elements, coords, 'pbe'));
 * ```
 *
 * Note: `init()` must be called once from the main `sci-form-wasm` entry point
 * before using any alpha functions.
 */

export {
  // ── A1: Kohn-Sham DFT ─────────────────────────────────────────────────────
  alpha_compute_dft,

  // ── A2: ReaxFF ────────────────────────────────────────────────────────────
  alpha_compute_reaxff_gradient,
  alpha_compute_reaxff_energy,
  alpha_compute_eem_charges,

  // ── A3: MLFF Neural Network Force Field ───────────────────────────────────
  alpha_compute_mlff,
  alpha_compute_aevs,

  // ── A4: Obara-Saika ERIs ──────────────────────────────────────────────────
  alpha_boys_function,
  alpha_eri_ssss,
  alpha_schwarz_bound,

  // ── A5: CGA Motor Algebra ─────────────────────────────────────────────────
  alpha_rotate_dihedral_cga,
  alpha_refine_torsion_cga,

  // ── A6: Growing String Method ─────────────────────────────────────────────
  alpha_gsm_interpolate,
  alpha_gsm_backend_plan,
  alpha_gsm_compare_backends,
  alpha_gsm_find_ts,
  alpha_gsm_find_ts_with_method,

  // ── A7: SDR Embedding ─────────────────────────────────────────────────────
  alpha_sdr_embed,

  // ── A8: Live MD (also available via main entry) ───────────────────────────
  LiveSimulation,

  // ── A9: 3D reaction dynamics ─────────────────────────────────────────────
  compute_reaction_dynamics_3d,
  compute_approach_dynamics,

  // ── Render bridge: chart utilities ────────────────────────────────────────
  alpha_render_edl_profile_chart,
  alpha_render_capacitance_scan_chart,
  alpha_render_periodic_dos_chart,
  alpha_render_band_structure_chart,
  alpha_render_microkinetic_population_chart,
  alpha_render_arrhenius_chart,
  alpha_render_gsm_path_trajectory,
  alpha_render_pack_chart_payload,
  alpha_render_verify_chart_round_trip,

  // ── Periodic linear ──────────────────────────────────────────────────────
  alpha_periodic_monkhorst_pack,
  alpha_periodic_bloch_phase,
  alpha_periodic_bz_integrate_scalar,
  alpha_periodic_bz_integrate_vector,
  alpha_periodic_validate_kmesh_weights,
  alpha_periodic_band_structure_eht,
  alpha_periodic_validate_electron_count,

  // ── Metadata ──────────────────────────────────────────────────────────────
  alpha_modules_info,
} from '../sci_form_wasm.js';
