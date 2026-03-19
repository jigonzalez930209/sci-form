use crate::distgeom::{
    calculate_bounds_matrix_opts, check_chiral_centers, check_double_bond_geometry,
    check_tetrahedral_centers, compute_initial_coords_rdkit, identify_chiral_sets,
    identify_tetrahedral_centers, pick_rdkit_distances, triangle_smooth_tol, MinstdRand,
    MAX_MINIMIZED_E_PER_ATOM,
};
use crate::forcefield::bounds_ff::minimize_bfgs_rdkit;
use crate::forcefield::etkdg_3d::{build_etkdg_3d_ff_with_torsions, minimize_etkdg_3d_bfgs};
use crate::graph::Molecule;
/// End-to-end 3D conformer generation pipeline.
///
/// Algorithm matching RDKit's ETKDG embedPoints() with retry-on-failure:
///   1. Build distance-bounds matrix → Floyd-Warshall smoothing
///   2. Identify chiral sets + tetrahedral centers (for validation)
///   3. Retry loop (up to 10×N attempts):
///      a. Sample random distances → metric matrix → 3D/4D embedding (4D only if chiral)
///      b. First minimization: bounds FF (chiral_w=1.0, 4d_w=0.1, basin=5.0), loop until converged
///      c. Energy/atom check: reject if energy/N ≥ 0.05
///      d. Tetrahedral center volume check
///      e. Chiral center sign check
///      f. Second minimization: bounds FF (chiral_w=0.2, 4d_w=1.0), loop until converged
///      g. Drop to 3D → ETKDG 3D FF minimization (300 iters, single pass)
///      h. Planarity check (OOP energy)
///      i. Double bond geometry check
use nalgebra::DMatrix;

const BASIN_THRESH: f32 = 5.0;
const FORCE_TOL: f32 = 1e-3;
const PLANARITY_TOLERANCE: f32 = 1.0;
const ERROR_TOL: f64 = 1e-5; // RDKit's ERROR_TOL for energy pre-check

/// Generate a 3D conformer from a SMILES string.
pub fn generate_3d_conformer_from_smiles(smiles: &str, seed: u64) -> Result<DMatrix<f32>, String> {
    let mol = Molecule::from_smiles(smiles)?;
    generate_3d_conformer(&mol, seed)
}

/// Generate a 3D conformer for an already-parsed `Molecule`.
///
/// Implements RDKit's embedPoints() retry-on-failure loop.
/// Returns the first valid 3D conformer, or an error if all attempts fail.
pub fn generate_3d_conformer(mol: &Molecule, seed: u64) -> Result<DMatrix<f32>, String> {
    let csd_torsions = crate::smarts::match_experimental_torsions(mol);
    generate_3d_conformer_with_torsions(mol, seed, &csd_torsions)
}

/// Generate multiple conformers with different seeds and return the one
/// with the lowest ETKDG 3D force field energy (best geometry).
pub fn generate_3d_conformer_best_of_k(
    mol: &Molecule,
    seed: u64,
    csd_torsions: &[crate::forcefield::etkdg_3d::M6TorsionContrib],
    num_seeds: usize,
) -> Result<DMatrix<f32>, String> {
    if num_seeds <= 1 {
        return generate_3d_conformer_with_torsions(mol, seed, csd_torsions);
    }

    let mut best: Option<(DMatrix<f32>, f64)> = None;
    let mut last_err = String::new();

    // Pre-compute bounds + FF scaffold once (topology-dependent, not seed-dependent)
    let bounds = {
        let raw = calculate_bounds_matrix_opts(mol, true);
        let mut b = raw;
        if triangle_smooth_tol(&mut b, 0.0) {
            b
        } else {
            let raw2 = calculate_bounds_matrix_opts(mol, false);
            let mut b2 = raw2.clone();
            if triangle_smooth_tol(&mut b2, 0.0) {
                b2
            } else {
                let mut b3 = raw2;
                triangle_smooth_tol(&mut b3, 0.05);
                b3
            }
        }
    };

    for k in 0..num_seeds {
        let s = seed.wrapping_add(k as u64 * 1000);
        match generate_3d_conformer_with_torsions(mol, s, csd_torsions) {
            Ok(coords) => {
                // Score using ETKDG 3D energy (topology-dependent, comparable across seeds)
                let n = mol.graph.node_count();
                let coords_f64 = coords.map(|v| v as f64);
                let ff = build_etkdg_3d_ff_with_torsions(mol, &coords_f64, &bounds, csd_torsions);
                let mut flat = vec![0.0f64; n * 3];
                for a in 0..n {
                    flat[a * 3] = coords[(a, 0)] as f64;
                    flat[a * 3 + 1] = coords[(a, 1)] as f64;
                    flat[a * 3 + 2] = coords[(a, 2)] as f64;
                }
                let energy = crate::forcefield::etkdg_3d::etkdg_3d_energy_f64(&flat, n, mol, &ff);
                match &best {
                    Some((_, best_e)) if energy >= *best_e => {}
                    _ => {
                        best = Some((coords, energy));
                    }
                }
            }
            Err(e) => {
                last_err = e;
            }
        }
    }

    match best {
        Some((coords, _)) => Ok(coords),
        None => Err(last_err),
    }
}

/// Generate a 3D conformer with optional CSD torsion overrides.
///
/// If `csd_torsions` is non-empty, they replace the default torsion terms
/// in the 3D force field, providing much better torsion angle quality.
pub fn generate_3d_conformer_with_torsions(
    mol: &Molecule,
    seed: u64,
    csd_torsions: &[crate::forcefield::etkdg_3d::M6TorsionContrib],
) -> Result<DMatrix<f32>, String> {
    let n = mol.graph.node_count();
    if n == 0 {
        return Err("Empty molecule".to_string());
    }

    // RDKit-style two-pass bounds: try with set15 + strict smoothing,
    // fall back to no set15 if triangle smoothing fails.
    let bounds = {
        let raw = calculate_bounds_matrix_opts(mol, true);
        let mut b = raw;
        if triangle_smooth_tol(&mut b, 0.0) {
            b
        } else {
            #[cfg(test)]
            eprintln!("  [FALLBACK] strict smoothing failed, retrying without set15");
            let raw2 = calculate_bounds_matrix_opts(mol, false);
            let mut b2 = raw2.clone();
            if triangle_smooth_tol(&mut b2, 0.0) {
                b2
            } else {
                #[cfg(test)]
                eprintln!("  [FALLBACK] second smoothing also failed, using soft smooth");
                let mut b3 = raw2;
                triangle_smooth_tol(&mut b3, 0.05);
                b3
            }
        }
    };
    let chiral_sets = identify_chiral_sets(mol);
    let tetrahedral_centers = identify_tetrahedral_centers(mol);

    let max_iterations = 10 * n;
    let mut rng = MinstdRand::new(seed as u32);

    // RDKit: 4D embedding only when chiral centers (CW/CCW) are present
    // Otherwise 3D embedding (no 4th dimension overhead)
    let use_4d = !chiral_sets.is_empty();
    let embed_dim = if use_4d { 4 } else { 3 };

    // Track consecutive embedding failures for random coord fallback
    // For large molecules (100+ atoms), eigendecomposition is O(N³) and fails more often,
    // so we switch to random coords sooner. For small molecules, we allow more eigen attempts.
    let mut consecutive_embed_fails = 0u32;
    let embed_fail_threshold = if n > 100 {
        (n as u32 / 8).max(10)
    } else {
        (n as u32 / 4).max(20)
    };
    let mut random_coord_attempts = 0u32;
    let max_random_coord_attempts = if n > 100 { 80u32 } else { 150u32 };
    // Scale BFGS restart limit: large molecules converge with fewer restarts
    let bfgs_restart_limit = if n > 100 { 20 } else { 50 };

    // Track energy-check failures for progressive relaxation
    let mut energy_check_failures = 0u32;
    // After 30% of iterations fail at the energy check, relax from 0.05 to 0.125
    let energy_relax_threshold = (max_iterations as f64 * 0.3) as u32;

    for _iter in 0..max_iterations {
        // Log attempt number if requested (works in both lib and integration tests)
        let _log_attempts = std::env::var("LOG_ATTEMPTS").is_ok();

        // Step 1: Generate initial coords
        let use_random_coords = consecutive_embed_fails >= embed_fail_threshold
            && random_coord_attempts < max_random_coord_attempts;
        let (mut coords, basin_thresh) = if use_random_coords {
            random_coord_attempts += 1;
            // Random coordinate fallback: RDKit uses boxSizeMult * cube_root(N)
            let box_size = 2.0 * (n as f64).cbrt().max(2.5);
            let mut c = DMatrix::from_element(n, embed_dim, 0.0f64);
            for i in 0..n {
                for d in 0..embed_dim {
                    c[(i, d)] = box_size * (rng.next_double() - 0.5);
                }
            }
            // RDKit disables basin threshold for random coords
            (c, 1e8f64)
        } else {
            let dists = pick_rdkit_distances(&mut rng, &bounds);
            let coords_opt = compute_initial_coords_rdkit(&mut rng, &dists, embed_dim);
            match coords_opt {
                Some(c) => {
                    consecutive_embed_fails = 0;
                    (c, BASIN_THRESH as f64)
                }
                None => {
                    consecutive_embed_fails += 1;
                    // If we've exhausted random coord attempts and still failing, give up
                    if random_coord_attempts >= max_random_coord_attempts {
                        break;
                    }
                    if _log_attempts && consecutive_embed_fails == embed_fail_threshold {
                        eprintln!(
                            "  attempt {} → switching to random coords after {} failures",
                            _iter, embed_fail_threshold
                        );
                    } else if _log_attempts {
                        eprintln!("  attempt {} → embedding failed", _iter);
                    }
                    continue;
                }
            }
        };

        // Step 2: First minimization (bounds FF with chiral_w=1.0, 4d_w=0.1)
        // RDKit: while(needMore) { needMore = field->minimize(400, forceTol); }
        // Safety limit: max 50 restarts to prevent infinite loops (RDKit typically needs < 10)
        {
            let bt = basin_thresh as f32;
            let initial_energy =
                compute_total_bounds_energy_f64(&coords, &bounds, &chiral_sets, bt, 0.1, 1.0);
            if initial_energy > ERROR_TOL {
                let mut need_more = 1;
                let mut restarts = 0;
                while need_more != 0 && restarts < bfgs_restart_limit {
                    need_more = minimize_bfgs_rdkit(
                        &mut coords,
                        &bounds,
                        &chiral_sets,
                        400,
                        FORCE_TOL as f64,
                        bt,
                        0.1,
                        1.0,
                    );
                    restarts += 1;
                }
            }
        }

        // Step 3: Energy per atom check with progressive relaxation
        // For difficult molecules (complex polycyclic), relax energy threshold after many failures
        let bt = basin_thresh as f32;
        let energy = compute_total_bounds_energy_f64(&coords, &bounds, &chiral_sets, bt, 0.1, 1.0);
        let effective_e_thresh = if energy_check_failures >= energy_relax_threshold {
            MAX_MINIMIZED_E_PER_ATOM as f64 * 2.5 // 0.125 — still rejects truly bad conformers
        } else {
            MAX_MINIMIZED_E_PER_ATOM as f64
        };
        if energy / n as f64 >= effective_e_thresh {
            energy_check_failures += 1;
            if _log_attempts {
                eprintln!(
                    "  attempt {} → energy check failed: {:.6}/atom",
                    _iter,
                    energy / n as f64
                );
            }
            continue;
        }

        // Step 4: Check tetrahedral centers (f64 coords matching RDKit's Point3D)
        if !check_tetrahedral_centers(&coords, &tetrahedral_centers) {
            if _log_attempts {
                eprintln!("  attempt {} → tetrahedral check failed", _iter);
            }
            continue;
        }

        // Step 5: Check chiral center volumes (f64 coords)
        if !chiral_sets.is_empty() && !check_chiral_centers(&coords, &chiral_sets) {
            if _log_attempts {
                eprintln!("  attempt {} → chiral check failed", _iter);
            }
            continue;
        }

        // Step 6: Second minimization (chiral_w=0.2, 4d_w=1.0) — only if 4D embedding
        // RDKit: while(needMore) { needMore = field2->minimize(200, forceTol); }
        if use_4d {
            let energy2 =
                compute_total_bounds_energy_f64(&coords, &bounds, &chiral_sets, bt, 1.0, 0.2);
            if energy2 > ERROR_TOL {
                let mut need_more = 1;
                let mut restarts = 0;
                while need_more != 0 && restarts < bfgs_restart_limit {
                    need_more = minimize_bfgs_rdkit(
                        &mut coords,
                        &bounds,
                        &chiral_sets,
                        200,
                        FORCE_TOL as f64,
                        bt,
                        1.0,
                        0.2,
                    );
                    restarts += 1;
                }
            }
        }

        // Step 7: Drop to 3D (no-op for 3D embedding)
        let coords3d = coords.columns(0, 3).into_owned();

        // Step 8: ETKDG 3D FF minimization — single pass of 300 iterations (matching RDKit)
        // Build FF using f64 coords for distance computation (matching RDKit's Point3D)
        let ff = build_etkdg_3d_ff_with_torsions(mol, &coords3d, &bounds, csd_torsions);
        // RDKit: only minimize if energy > ERROR_TOL
        let e3d = crate::forcefield::etkdg_3d::etkdg_3d_energy_f64(
            &{
                let n = mol.graph.node_count();
                let mut flat = vec![0.0f64; n * 3];
                for a in 0..n {
                    flat[a * 3] = coords3d[(a, 0)];
                    flat[a * 3 + 1] = coords3d[(a, 1)];
                    flat[a * 3 + 2] = coords3d[(a, 2)];
                }
                flat
            },
            mol.graph.node_count(),
            mol,
            &ff,
        );
        let refined = if e3d > ERROR_TOL {
            minimize_etkdg_3d_bfgs(mol, &coords3d, &ff, 300, FORCE_TOL)
        } else {
            coords3d
        };

        // Step 9: Planarity check — matching RDKit's construct3DImproperForceField
        // Uses UFF inversion energy + SP angle constraint energy (k=10.0) in f64
        {
            let n_improper_atoms = ff.inversion_contribs.len() / 3;
            let flat_f64: Vec<f64> = {
                let nr = refined.nrows();
                let mut flat = vec![0.0f64; nr * 3];
                for a in 0..nr {
                    flat[a * 3] = refined[(a, 0)];
                    flat[a * 3 + 1] = refined[(a, 1)];
                    flat[a * 3 + 2] = refined[(a, 2)];
                }
                flat
            };
            let planarity_energy =
                crate::forcefield::etkdg_3d::planarity_check_energy_f64(&flat_f64, n, &ff);
            if planarity_energy > n_improper_atoms as f64 * PLANARITY_TOLERANCE as f64 {
                if _log_attempts {
                    eprintln!(
                        "  attempt {} → planarity check failed (energy={:.4} > threshold={:.4})",
                        _iter,
                        planarity_energy,
                        n_improper_atoms as f64 * PLANARITY_TOLERANCE as f64
                    );
                }
                continue;
            }
        }

        // Step 10: Double bond geometry check (f64 coords)
        if !check_double_bond_geometry(mol, &refined) {
            if _log_attempts {
                eprintln!("  attempt {} → double bond check failed", _iter);
            }
            continue;
        }

        if _log_attempts {
            eprintln!("  attempt {} → SUCCESS", _iter);
        }

        let refined_f32 = refined.map(|v| v as f32);
        return Ok(refined_f32);
    }

    Err(format!(
        "Failed to generate valid conformer after {} iterations",
        max_iterations
    ))
}

/// Compute bounds FF energy in f64, matching RDKit's field->calcEnergy() exactly.
/// Used for the energy pre-check before minimization.
pub fn compute_total_bounds_energy_f64(
    coords: &DMatrix<f64>,
    bounds: &DMatrix<f64>,
    chiral_sets: &[crate::forcefield::bounds_ff::ChiralSet],
    basin_thresh: f32,
    weight_4d: f32,
    weight_chiral: f32,
) -> f64 {
    let n = coords.nrows();
    let dim_coords = coords.ncols();
    let basin_thresh_f64 = basin_thresh as f64;
    let weight_4d_f64 = weight_4d as f64;
    let weight_chiral_f64 = weight_chiral as f64;

    let mut energy = 0.0f64;
    for i in 1..n {
        for j in 0..i {
            let ub = bounds[(j, i)];
            let lb = bounds[(i, j)];
            if ub - lb > basin_thresh_f64 {
                continue;
            }
            let mut d2 = 0.0f64;
            for d in 0..dim_coords {
                let diff = coords[(i, d)] - coords[(j, d)];
                d2 += diff * diff;
            }
            let ub2 = ub * ub;
            let lb2 = lb * lb;
            let val = if d2 > ub2 {
                d2 / ub2 - 1.0
            } else if d2 < lb2 {
                2.0 * lb2 / (lb2 + d2) - 1.0
            } else {
                0.0
            };
            if val > 0.0 {
                energy += val * val;
            }
        }
    }
    if !chiral_sets.is_empty() {
        // Flatten coords to flat array for f64 chiral energy
        let mut flat = vec![0.0f64; n * dim_coords];
        for i in 0..n {
            for d in 0..dim_coords {
                flat[i * dim_coords + d] = coords[(i, d)];
            }
        }
        energy += weight_chiral_f64
            * crate::forcefield::bounds_ff::chiral_violation_energy_f64(
                &flat,
                dim_coords,
                chiral_sets,
            );
    }
    if dim_coords == 4 {
        for i in 0..n {
            let x4 = coords[(i, 3)];
            energy += weight_4d_f64 * x4 * x4;
        }
    }
    energy
}
