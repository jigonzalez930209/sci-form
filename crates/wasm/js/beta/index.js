/**
 * sci-form-wasm/beta — Beta-tier experimental modules.
 *
 * Validated methods approaching stable release. APIs are relatively stable
 * but may have minor changes before promotion to core.
 *
 * Usage:
 * ```js
 * import init from 'sci-form-wasm';
 * import {
 *   beta_compute_kpm_dos,
 *   beta_compute_kpm_mulliken,
 *   beta_compute_mbh_frequencies,
 *   beta_solve_eht_randnla,
 *   beta_psd_distance,
 *   beta_psd_projection,
 *   beta_compute_cpm_charges,
 *   beta_modules_info,
 * } from 'sci-form-wasm/beta';
 *
 * await init();
 * const dos = JSON.parse(beta_compute_kpm_dos(elements, coords, '{}'));
 * ```
 */

export {
  // ── B1: KPM — Kernel Polynomial Method ────────────────────────────────────
  beta_compute_kpm_dos,
  beta_compute_kpm_mulliken,

  // ── B2: MBH — Mobile Block Hessian ─────────────────────────────────────────
  beta_compute_mbh_frequencies,

  // ── B3: RandNLA — Randomized EHT Eigensolver ──────────────────────────────
  beta_solve_eht_randnla,

  // ── B4: Riemannian Geometry ────────────────────────────────────────────────
  beta_psd_distance,
  beta_psd_projection,

  // ── B5: CPM — Constant Potential Method ───────────────────────────────────
  beta_compute_cpm_charges,

  // ── Metadata ──────────────────────────────────────────────────────────────
  beta_modules_info,
} from '../sci_form_wasm.js';
