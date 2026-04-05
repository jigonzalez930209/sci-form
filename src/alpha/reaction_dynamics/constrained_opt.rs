//! Constrained geometry optimization for reaction frames.
//!
//! At each point along the reaction coordinate, the internal geometry of
//! the molecular system is relaxed while keeping the reaction coordinate
//! (typically the distance between reactive atoms) frozen.
//!
//! This produces physically meaningful geometries where bonds, angles,
//! and dihedrals adjust naturally to the changing inter-fragment distance,
//! rather than naive rigid-body translations.
//!
//! Two modes:
//! 1. **Distance-constrained relaxation**: Fix distance(s) between specific
//!    atom pairs, relax everything else. Used for approach/departure frames.
//! 2. **NEB-image relaxation**: Fix the reaction coordinate (projection along
//!    the NEB tangent), relax perpendicular DOFs. Used to clean up NEB images.

/// A distance constraint between two atoms.
#[derive(Debug, Clone)]
pub struct DistanceConstraint {
    /// First atom index
    pub atom_a: usize,
    /// Second atom index
    pub atom_b: usize,
    /// Target distance (Å)
    pub target_distance: f64,
}

/// Configuration for constrained optimization.
#[derive(Debug, Clone)]
pub struct ConstrainedOptConfig {
    /// Maximum optimization steps per frame.
    pub max_steps: usize,
    /// Gradient norm convergence threshold (kcal/mol/Å).
    pub grad_threshold: f64,
    /// Step size for L-BFGS-like updates (Å).
    pub step_size: f64,
    /// NEB method backend for energy/gradient evaluation.
    pub method: String,
    /// SHAKE-like constraint tolerance for distance enforcement.
    pub constraint_tolerance: f64,
}

impl Default for ConstrainedOptConfig {
    fn default() -> Self {
        Self {
            max_steps: 50,
            grad_threshold: 0.05,
            step_size: 0.01,
            method: "gfn2".into(),
            constraint_tolerance: 1e-4,
        }
    }
}

/// Result of constrained optimization.
#[derive(Debug, Clone)]
pub struct ConstrainedOptResult {
    /// Optimized coordinates (flat).
    pub coords: Vec<f64>,
    /// Final energy (kcal/mol).
    pub energy: f64,
    /// Number of steps taken.
    pub n_steps: usize,
    /// Whether it converged.
    pub converged: bool,
    /// Final gradient norm (kcal/mol/Å).
    pub grad_norm: f64,
}

/// Optimize a molecular geometry with distance constraints.
///
/// Uses projected gradient descent: compute the full gradient, then
/// remove the components that would violate the constraints (SHAKE-like
/// projection), and take a step along the projected gradient.
///
/// This is physically equivalent to a "relaxed scan" — at each fixed
/// value of the reaction coordinate, the system sits at a local minimum
/// on the constrained PES.
pub fn optimize_with_constraints(
    smiles: &str,
    elements: &[u8],
    initial_coords: &[f64],
    constraints: &[DistanceConstraint],
    config: &ConstrainedOptConfig,
) -> Result<ConstrainedOptResult, String> {
    let backend = crate::dynamics::NebBackend::from_method(&config.method)?;
    let mol = crate::graph::Molecule::from_smiles(smiles)?;
    let n_xyz = initial_coords.len();

    let mut x = initial_coords.to_vec();

    // Enforce constraints on initial geometry
    enforce_distance_constraints(&mut x, constraints, config.constraint_tolerance);

    let mut converged = false;
    let mut n_steps = 0;
    let mut energy = 0.0;
    let mut grad_norm = f64::INFINITY;

    // L-BFGS memory (limited to last m steps)
    let m = 5; // memory depth
    let mut s_history: Vec<Vec<f64>> = Vec::new();
    let mut y_history: Vec<Vec<f64>> = Vec::new();
    let mut rho_history: Vec<f64> = Vec::new();
    let mut prev_x: Option<Vec<f64>> = None;
    let mut prev_grad: Option<Vec<f64>> = None;

    for step in 0..config.max_steps {
        n_steps = step + 1;

        let mut grad = vec![0.0; n_xyz];
        energy = crate::dynamics::neb_energy_and_gradient(
            backend, smiles, elements, &mol, &x, &mut grad,
        )?;

        // Project gradient to remove constrained components
        project_gradient_onto_constraints(&mut grad, &x, constraints);

        grad_norm = grad.iter().map(|g| g * g).sum::<f64>().sqrt() / (n_xyz as f64).sqrt();

        if grad_norm < config.grad_threshold {
            converged = true;
            break;
        }

        // L-BFGS two-loop recursion to compute search direction
        let direction = if let (Some(ref px), Some(ref pg)) = (&prev_x, &prev_grad) {
            // Update L-BFGS memory
            let s: Vec<f64> = x.iter().zip(px.iter()).map(|(xi, pi)| xi - pi).collect();
            let y: Vec<f64> = grad.iter().zip(pg.iter()).map(|(gi, pi)| gi - pi).collect();
            let sy: f64 = s.iter().zip(y.iter()).map(|(a, b)| a * b).sum();

            if sy > 1e-10 {
                let rho = 1.0 / sy;
                s_history.push(s);
                y_history.push(y);
                rho_history.push(rho);

                // Keep only last m entries
                if s_history.len() > m {
                    s_history.remove(0);
                    y_history.remove(0);
                    rho_history.remove(0);
                }
            }

            lbfgs_direction(&grad, &s_history, &y_history, &rho_history)
        } else {
            // First step: use negative gradient
            grad.iter().map(|g| -g).collect()
        };

        prev_x = Some(x.clone());
        prev_grad = Some(grad.clone());

        // Line search with backtracking
        let mut alpha = config.step_size;
        let g_dot_d: f64 = grad.iter().zip(direction.iter()).map(|(g, d)| g * d).sum();

        for _ in 0..10 {
            let mut x_trial: Vec<f64> = x
                .iter()
                .zip(direction.iter())
                .map(|(&xi, &di)| xi + alpha * di)
                .collect();

            // Enforce constraints after each step
            enforce_distance_constraints(&mut x_trial, constraints, config.constraint_tolerance);

            let mut grad_trial = vec![0.0; n_xyz];
            let e_trial = crate::dynamics::neb_energy_and_gradient(
                backend,
                smiles,
                elements,
                &mol,
                &x_trial,
                &mut grad_trial,
            )?;

            // Sufficient decrease (Armijo)
            if e_trial <= energy + 1e-4 * alpha * g_dot_d || alpha < 1e-6 {
                x = x_trial;
                energy = e_trial;
                break;
            }
            alpha *= 0.5;
        }

        // Final constraint enforcement
        enforce_distance_constraints(&mut x, constraints, config.constraint_tolerance);
    }

    Ok(ConstrainedOptResult {
        coords: x,
        energy,
        n_steps,
        converged,
        grad_norm,
    })
}

/// L-BFGS two-loop recursion.
fn lbfgs_direction(
    grad: &[f64],
    s_hist: &[Vec<f64>],
    y_hist: &[Vec<f64>],
    rho_hist: &[f64],
) -> Vec<f64> {
    let n = grad.len();
    let k = s_hist.len();
    if k == 0 {
        return grad.iter().map(|g| -g).collect();
    }

    let mut q = grad.to_vec();
    let mut alphas = vec![0.0; k];

    // Forward pass
    for i in (0..k).rev() {
        alphas[i] = rho_hist[i]
            * s_hist[i]
                .iter()
                .zip(q.iter())
                .map(|(s, qi)| s * qi)
                .sum::<f64>();
        for j in 0..n {
            q[j] -= alphas[i] * y_hist[i][j];
        }
    }

    // Initial Hessian scaling
    let sy: f64 = s_hist[k - 1]
        .iter()
        .zip(y_hist[k - 1].iter())
        .map(|(s, y)| s * y)
        .sum();
    let yy: f64 = y_hist[k - 1].iter().map(|y| y * y).sum();
    let gamma = if yy > 1e-12 { sy / yy } else { 1.0 };

    let mut r: Vec<f64> = q.iter().map(|qi| gamma * qi).collect();

    // Backward pass
    for i in 0..k {
        let beta = rho_hist[i]
            * y_hist[i]
                .iter()
                .zip(r.iter())
                .map(|(y, ri)| y * ri)
                .sum::<f64>();
        for j in 0..n {
            r[j] += s_hist[i][j] * (alphas[i] - beta);
        }
    }

    // Return negative direction (descent)
    r.iter().map(|ri| -ri).collect()
}

/// Enforce distance constraints using SHAKE-like iterative projection.
///
/// Adjusts atom positions to satisfy distance constraints while minimally
/// perturbing the geometry. Uses iterative projection (RATTLE/SHAKE algorithm).
fn enforce_distance_constraints(
    coords: &mut [f64],
    constraints: &[DistanceConstraint],
    tolerance: f64,
) {
    let max_shake_iter = 100;

    for _ in 0..max_shake_iter {
        let mut max_violation = 0.0f64;

        for c in constraints {
            let a = c.atom_a;
            let b = c.atom_b;
            let dx = coords[b * 3] - coords[a * 3];
            let dy = coords[b * 3 + 1] - coords[a * 3 + 1];
            let dz = coords[b * 3 + 2] - coords[a * 3 + 2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();

            let violation = (dist - c.target_distance).abs();
            if violation > max_violation {
                max_violation = violation;
            }

            if dist < 1e-10 {
                continue;
            }

            // Correction: move atoms symmetrically to satisfy constraint
            let correction = 0.5 * (dist - c.target_distance) / dist;
            coords[a * 3] += correction * dx;
            coords[a * 3 + 1] += correction * dy;
            coords[a * 3 + 2] += correction * dz;
            coords[b * 3] -= correction * dx;
            coords[b * 3 + 1] -= correction * dy;
            coords[b * 3 + 2] -= correction * dz;
        }

        if max_violation < tolerance {
            break;
        }
    }
}

/// Project gradient to remove components along constrained directions.
///
/// For each distance constraint (a,b), removes the gradient component
/// along the a→b direction from both atoms. This prevents the optimizer
/// from trying to change the constrained distance.
fn project_gradient_onto_constraints(
    grad: &mut [f64],
    coords: &[f64],
    constraints: &[DistanceConstraint],
) {
    for c in constraints {
        let a = c.atom_a;
        let b = c.atom_b;
        let dx = coords[b * 3] - coords[a * 3];
        let dy = coords[b * 3 + 1] - coords[a * 3 + 1];
        let dz = coords[b * 3 + 2] - coords[a * 3 + 2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

        if dist < 1e-10 {
            continue;
        }

        // Unit vector along constraint
        let ux = dx / dist;
        let uy = dy / dist;
        let uz = dz / dist;

        // Remove gradient component along constraint for atom a
        let dot_a = grad[a * 3] * ux + grad[a * 3 + 1] * uy + grad[a * 3 + 2] * uz;
        grad[a * 3] -= dot_a * ux;
        grad[a * 3 + 1] -= dot_a * uy;
        grad[a * 3 + 2] -= dot_a * uz;

        // Remove gradient component along constraint for atom b
        let dot_b = grad[b * 3] * ux + grad[b * 3 + 1] * uy + grad[b * 3 + 2] * uz;
        grad[b * 3] -= dot_b * ux;
        grad[b * 3 + 1] -= dot_b * uy;
        grad[b * 3 + 2] -= dot_b * uz;
    }
}

/// Optimize a sequence of frames along a reaction path.
///
/// For each frame, identifies the reactive atom pair (closest inter-fragment
/// pair) and performs constrained optimization holding that distance fixed
/// while relaxing all internal DOFs.
///
/// This is the core "relaxed scan" that makes approach/departure frames
/// chemically realistic.
pub fn optimize_frame_sequence(
    smiles: &str,
    elements: &[u8],
    frames: &[Vec<f64>],
    mol_ranges: &[(usize, usize)],
    config: &ConstrainedOptConfig,
) -> Result<Vec<ConstrainedOptResult>, String> {
    if mol_ranges.len() < 2 {
        // Single molecule: no inter-fragment constraints, just optimize each
        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            let results: Vec<Result<ConstrainedOptResult, String>> = frames
                .par_iter()
                .map(|frame| optimize_with_constraints(smiles, elements, frame, &[], config))
                .collect();
            return results.into_iter().collect();
        }
        #[cfg(not(feature = "parallel"))]
        {
            return frames
                .iter()
                .map(|frame| optimize_with_constraints(smiles, elements, frame, &[], config))
                .collect();
        }
    }

    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        let results: Vec<Result<ConstrainedOptResult, String>> = frames
            .par_iter()
            .map(|frame| {
                let constraints = find_reactive_constraints(frame, mol_ranges);
                optimize_with_constraints(smiles, elements, frame, &constraints, config)
            })
            .collect();
        results.into_iter().collect()
    }
    #[cfg(not(feature = "parallel"))]
    {
        frames
            .iter()
            .map(|frame| {
                let constraints = find_reactive_constraints(frame, mol_ranges);
                optimize_with_constraints(smiles, elements, frame, &constraints, config)
            })
            .collect()
    }
}

/// Find the reactive atom pair (closest inter-fragment pair) and create constraints.
fn find_reactive_constraints(
    coords: &[f64],
    mol_ranges: &[(usize, usize)],
) -> Vec<DistanceConstraint> {
    if mol_ranges.len() < 2 {
        return vec![];
    }

    let mut best_a = 0usize;
    let mut best_b = 0usize;
    let mut best_d2 = f64::INFINITY;

    // Check all inter-fragment atom pairs between first two fragments
    let (start_a, end_a) = mol_ranges[0];
    let (start_b, end_b) = mol_ranges[1];

    for i in start_a..end_a {
        for j in start_b..end_b {
            let dx = coords[j * 3] - coords[i * 3];
            let dy = coords[j * 3 + 1] - coords[i * 3 + 1];
            let dz = coords[j * 3 + 2] - coords[i * 3 + 2];
            let d2 = dx * dx + dy * dy + dz * dz;
            if d2 < best_d2 {
                best_d2 = d2;
                best_a = i;
                best_b = j;
            }
        }
    }

    vec![DistanceConstraint {
        atom_a: best_a,
        atom_b: best_b,
        target_distance: best_d2.sqrt(),
    }]
}

/// Optimize NEB images with constrained relaxation.
///
/// For each interior NEB image:
///   1. Project the image position onto the NEB path tangent (reaction coordinate)
///   2. Optimize perpendicular DOFs with L-BFGS
///   3. Enforce that the image stays on the NEB path ribbon
///
/// This ensures each image represents a true constrained minimum perpendicular
/// to the reaction path, not just an artifact of the elastic band forces.
pub fn relax_neb_images(
    smiles: &str,
    elements: &[u8],
    images: &mut [Vec<f64>],
    config: &ConstrainedOptConfig,
) -> Result<Vec<f64>, String> {
    let backend = crate::dynamics::NebBackend::from_method(&config.method)?;
    let mol = crate::graph::Molecule::from_smiles(smiles)?;
    let n_images = images.len();
    let n_xyz = images[0].len();

    // Compute initial energies (parallelized)
    let mut energies = vec![0.0f64; n_images];
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        let par_e: Vec<Result<f64, String>> = images
            .par_iter()
            .map(|img| {
                let mut grad = vec![0.0; n_xyz];
                crate::dynamics::neb_energy_and_gradient(
                    backend, smiles, elements, &mol, img, &mut grad,
                )
            })
            .collect();
        for (i, r) in par_e.into_iter().enumerate() {
            energies[i] = r?;
        }
    }
    #[cfg(not(feature = "parallel"))]
    {
        for i in 0..n_images {
            let mut grad = vec![0.0; n_xyz];
            energies[i] = crate::dynamics::neb_energy_and_gradient(
                backend, smiles, elements, &mol, &images[i], &mut grad,
            )?;
        }
    }

    // Compute tangents before parallel relaxation (need neighbors)
    let mut tangents = Vec::with_capacity(n_images);
    tangents.push(vec![0.0f64; n_xyz]); // endpoint placeholder
    for i in 1..(n_images - 1) {
        let mut tau = vec![0.0f64; n_xyz];
        for k in 0..n_xyz {
            tau[k] = images[i + 1][k] - images[i - 1][k];
        }
        let tau_norm = tau.iter().map(|t| t * t).sum::<f64>().sqrt();
        if tau_norm > 1e-12 {
            for t in &mut tau {
                *t /= tau_norm;
            }
        }
        tangents.push(tau);
    }
    tangents.push(vec![0.0f64; n_xyz]); // endpoint placeholder

    // Relax interior images in parallel
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        let interior_results: Vec<Result<(Vec<f64>, f64), String>> = (1..(n_images - 1))
            .into_par_iter()
            .map(|i| {
                let tau = &tangents[i];
                let mut img = images[i].clone();

                for _step in 0..config.max_steps {
                    let mut grad = vec![0.0; n_xyz];
                    let _e = crate::dynamics::neb_energy_and_gradient(
                        backend, smiles, elements, &mol, &img, &mut grad,
                    )?;

                    let dot: f64 = grad.iter().zip(tau.iter()).map(|(g, t)| g * t).sum();
                    for k in 0..n_xyz {
                        grad[k] -= dot * tau[k];
                    }

                    let gnorm = (grad.iter().map(|g| g * g).sum::<f64>() / n_xyz as f64).sqrt();
                    if gnorm < config.grad_threshold {
                        break;
                    }

                    for k in 0..n_xyz {
                        img[k] -= config.step_size * grad[k];
                    }
                }

                let mut grad = vec![0.0; n_xyz];
                let e = crate::dynamics::neb_energy_and_gradient(
                    backend, smiles, elements, &mol, &img, &mut grad,
                )?;
                Ok((img, e))
            })
            .collect();

        for (idx, r) in interior_results.into_iter().enumerate() {
            let i = idx + 1;
            let (new_img, e) = r?;
            images[i] = new_img;
            energies[i] = e;
        }
    }
    #[cfg(not(feature = "parallel"))]
    {
        for i in 1..(n_images - 1) {
            let tau = &tangents[i];

            for _step in 0..config.max_steps {
                let mut grad = vec![0.0; n_xyz];
                let _e = crate::dynamics::neb_energy_and_gradient(
                    backend, smiles, elements, &mol, &images[i], &mut grad,
                )?;

                let dot: f64 = grad.iter().zip(tau.iter()).map(|(g, t)| g * t).sum();
                for k in 0..n_xyz {
                    grad[k] -= dot * tau[k];
                }

                let gnorm = (grad.iter().map(|g| g * g).sum::<f64>() / n_xyz as f64).sqrt();
                if gnorm < config.grad_threshold {
                    break;
                }

                for k in 0..n_xyz {
                    images[i][k] -= config.step_size * grad[k];
                }
            }

            let mut grad = vec![0.0; n_xyz];
            energies[i] = crate::dynamics::neb_energy_and_gradient(
                backend, smiles, elements, &mol, &images[i], &mut grad,
            )?;
        }
    }

    Ok(energies)
}
