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
 *   alpha_gsm_find_ts,
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
  alpha_gsm_find_ts,

  // ── A7: SDR Embedding ─────────────────────────────────────────────────────
  alpha_sdr_embed,

  // ── A8: Live MD (also available via main entry) ───────────────────────────
  LiveSimulation,

  // ── Metadata ──────────────────────────────────────────────────────────────
  alpha_modules_info,
} from '../sci_form_wasm.js';
