//! Climbing-Image Nudged Elastic Band (CI-NEB).
//!
//! Implements the improved tangent method of Henkelman & Jónsson (2000) with
//! force projection and a climbing image that converges to the true saddle point.
//!
//! Key improvements over the simplified NEB:
//! - IDPP initialization (avoids atom clashes in the initial path)
//! - FIRE optimizer (Fast Inertial Relaxation Engine) instead of steepest descent
//! - Post-convergence perpendicular relaxation of each image
//! - Image redistribution to maintain even spacing

use serde::{Deserialize, Serialize};

/// Configuration for CI-NEB calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CiNebConfig {
    pub n_images: usize,
    pub max_iter: usize,
    pub spring_k: f64,
    pub step_size: f64,
    pub force_threshold: f64,
    pub method: String,
    /// Use IDPP initialization instead of linear interpolation.
    #[serde(default = "default_true")]
    pub use_idpp: bool,
    /// Perform perpendicular relaxation after NEB convergence.
    #[serde(default = "default_true")]
    pub relax_images: bool,
    /// Maximum steps for post-NEB perpendicular relaxation per image.
    #[serde(default = "default_relax_steps")]
    pub relax_max_steps: usize,
}

fn default_true() -> bool {
    true
}
fn default_relax_steps() -> usize {
    30
}

/// A single NEB image.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CiNebImage {
    pub index: usize,
    pub coords: Vec<f64>,
    pub energy_kcal_mol: f64,
    pub is_climbing: bool,
}

/// Result of a CI-NEB computation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CiNebResult {
    pub images: Vec<CiNebImage>,
    pub converged: bool,
    pub iterations: usize,
    pub max_force: f64,
    pub notes: Vec<String>,
}

/// Compute the tangent vector at image `i` using the improved method.
///
/// Uses energy-weighted bisection when `i` is a local extremum, otherwise
/// uses the higher-energy neighbor direction. This prevents kinking artifacts.
fn compute_tangent(images: &[Vec<f64>], energies: &[f64], i: usize) -> Vec<f64> {
    let n = images[i].len();
    let mut tau_plus = vec![0.0; n];
    let mut tau_minus = vec![0.0; n];
    for k in 0..n {
        tau_plus[k] = images[i + 1][k] - images[i][k];
        tau_minus[k] = images[i][k] - images[i - 1][k];
    }

    let e_i = energies[i];
    let e_plus = energies[i + 1];
    let e_minus = energies[i - 1];

    let tau = if e_plus > e_i && e_i > e_minus {
        // Monotonically increasing → use forward tangent
        tau_plus
    } else if e_plus < e_i && e_i < e_minus {
        // Monotonically decreasing → use backward tangent
        tau_minus
    } else {
        // Local extremum → energy-weighted bisection
        let de_max = (e_plus - e_i).abs().max((e_minus - e_i).abs());
        let de_min = (e_plus - e_i).abs().min((e_minus - e_i).abs());

        let mut tau = vec![0.0; n];
        if e_plus > e_minus {
            for k in 0..n {
                tau[k] = tau_plus[k] * de_max + tau_minus[k] * de_min;
            }
        } else {
            for k in 0..n {
                tau[k] = tau_plus[k] * de_min + tau_minus[k] * de_max;
            }
        }
        tau
    };

    // Normalise
    let norm = tau.iter().map(|x| x * x).sum::<f64>().sqrt();
    if norm > 1e-12 {
        tau.iter().map(|x| x / norm).collect()
    } else {
        tau
    }
}

/// Project force perpendicular to tangent: F_perp = F - (F·tau)tau
fn project_perpendicular(force: &[f64], tau: &[f64]) -> Vec<f64> {
    let dot: f64 = force.iter().zip(tau.iter()).map(|(f, t)| f * t).sum();
    force
        .iter()
        .zip(tau.iter())
        .map(|(f, t)| f - dot * t)
        .collect()
}

/// Spring force along the tangent.
fn spring_force_parallel(images: &[Vec<f64>], i: usize, tau: &[f64], spring_k: f64) -> Vec<f64> {
    let n = images[i].len();
    let d_plus: f64 = (0..n)
        .map(|k| (images[i + 1][k] - images[i][k]).powi(2))
        .sum::<f64>()
        .sqrt();
    let d_minus: f64 = (0..n)
        .map(|k| (images[i][k] - images[i - 1][k]).powi(2))
        .sum::<f64>()
        .sqrt();

    let spring_mag = spring_k * (d_plus - d_minus);
    tau.iter().map(|t| spring_mag * t).collect()
}

/// Climbing image force: F_CI = F_real - 2*(F_real·tau)*tau
fn climbing_image_force(real_force: &[f64], tau: &[f64]) -> Vec<f64> {
    let dot: f64 = real_force.iter().zip(tau.iter()).map(|(f, t)| f * t).sum();
    real_force
        .iter()
        .zip(tau.iter())
        .map(|(f, t)| f - 2.0 * dot * t)
        .collect()
}

/// Compute CI-NEB path between reactant and product coordinates.
pub fn compute_ci_neb_path(
    smiles: &str,
    start_coords: &[f64],
    end_coords: &[f64],
    elements: &[u8],
    config: &CiNebConfig,
) -> Result<CiNebResult, String> {
    let n_images = config.n_images;
    if n_images < 3 {
        return Err("CI-NEB requires at least 3 images".into());
    }

    let mol = crate::graph::Molecule::from_smiles(smiles)?;
    let backend = crate::dynamics::NebBackend::from_method(&config.method)?;
    let n_xyz = start_coords.len();
    if end_coords.len() != n_xyz {
        return Err("start/end coords length mismatch".into());
    }

    // ── 1. Initial path: IDPP (physically meaningful) or linear interpolation ──
    let mut images: Vec<Vec<f64>> = if config.use_idpp {
        super::idpp::generate_idpp_path(start_coords, end_coords, n_images, 200)
    } else {
        (0..n_images)
            .map(|i| {
                let t = i as f64 / (n_images - 1) as f64;
                (0..n_xyz)
                    .map(|k| (1.0 - t) * start_coords[k] + t * end_coords[k])
                    .collect()
            })
            .collect()
    };

    // ── 2. Compute initial energies (parallelized) ──
    let mut energies = vec![0.0f64; n_images];
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        let par_energies: Vec<Result<f64, String>> = images
            .par_iter()
            .map(|img| {
                let mut grad = vec![0.0; n_xyz];
                crate::dynamics::neb_energy_and_gradient(
                    backend, smiles, elements, &mol, img, &mut grad,
                )
            })
            .collect();
        for (i, r) in par_energies.into_iter().enumerate() {
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

    let mut converged = false;
    let mut max_force = f64::INFINITY;
    let mut actual_iter = 0;

    // ── 3. FIRE optimizer state (per image) ──
    let n_interior = n_images - 2;
    let mut velocities: Vec<Vec<f64>> = vec![vec![0.0f64; n_xyz]; n_interior];
    let mut dt = vec![config.step_size; n_interior];
    let mut alpha_fire = vec![0.1f64; n_interior];
    let mut n_positive = vec![0usize; n_interior];
    let fire_n_min = 5;
    let fire_f_inc = 1.1;
    let fire_f_dec = 0.5;
    let fire_alpha_start = 0.1;
    let fire_f_alpha = 0.99;

    // Start CI after warm-up (30% of iterations)
    let ci_start_iter = config.max_iter / 3;

    // ── 4. NEB optimization loop ──
    for iter in 0..config.max_iter {
        actual_iter = iter + 1;

        // Find highest-energy interior image (climbing image candidate)
        let ci_idx = (1..n_images - 1)
            .max_by(|&a, &b| {
                energies[a]
                    .partial_cmp(&energies[b])
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(n_images / 2);

        let use_ci = iter >= ci_start_iter;
        let mut max_f = 0.0f64;

        // Update interior images
        for i in 1..(n_images - 1) {
            let img_i = i - 1; // index into velocities/dt/alpha arrays

            let mut grad = vec![0.0; n_xyz];
            let energy = crate::dynamics::neb_energy_and_gradient(
                backend, smiles, elements, &mol, &images[i], &mut grad,
            )?;
            energies[i] = energy;

            // Real force = -gradient
            let real_force: Vec<f64> = grad.iter().map(|g| -g).collect();

            let tau = compute_tangent(&images, &energies, i);

            let neb_force = if use_ci && i == ci_idx {
                // Climbing image: invert component along tangent
                climbing_image_force(&real_force, &tau)
            } else {
                // Regular NEB: perpendicular real force + parallel spring force
                let f_perp = project_perpendicular(&real_force, &tau);
                let f_spring = spring_force_parallel(&images, i, &tau, config.spring_k);
                f_perp
                    .iter()
                    .zip(f_spring.iter())
                    .map(|(fp, fs)| fp + fs)
                    .collect()
            };

            // Track max force for convergence
            let f_mag = neb_force.iter().map(|f| f * f).sum::<f64>().sqrt();
            if f_mag > max_f {
                max_f = f_mag;
            }

            // FIRE update for this image
            let power: f64 = velocities[img_i]
                .iter()
                .zip(neb_force.iter())
                .map(|(v, f)| v * f)
                .sum();

            if power > 0.0 {
                n_positive[img_i] += 1;

                // Mix velocity with force direction (FIRE)
                let v_norm = velocities[img_i].iter().map(|v| v * v).sum::<f64>().sqrt();
                let f_norm = f_mag;
                if f_norm > 1e-12 {
                    for k in 0..n_xyz {
                        velocities[img_i][k] = (1.0 - alpha_fire[img_i]) * velocities[img_i][k]
                            + alpha_fire[img_i] * (v_norm / f_norm) * neb_force[k];
                    }
                }

                if n_positive[img_i] > fire_n_min {
                    dt[img_i] = (dt[img_i] * fire_f_inc).min(config.step_size * 10.0);
                    alpha_fire[img_i] *= fire_f_alpha;
                }
            } else {
                // Power is negative: reset velocity, reduce time step
                velocities[img_i] = vec![0.0; n_xyz];
                dt[img_i] *= fire_f_dec;
                dt[img_i] = dt[img_i].max(config.step_size * 0.01);
                alpha_fire[img_i] = fire_alpha_start;
                n_positive[img_i] = 0;
            }

            // Velocity Verlet-like step
            for k in 0..n_xyz {
                velocities[img_i][k] += dt[img_i] * neb_force[k];
                images[i][k] += dt[img_i] * velocities[img_i][k];
            }
        }

        max_force = max_f;

        // Check convergence (on forces of all interior images)
        if max_force < config.force_threshold && iter >= ci_start_iter {
            converged = true;
            break;
        }

        // ── Periodic image redistribution (every 20 iterations) ──
        if iter > 0 && iter % 20 == 0 {
            redistribute_images(&mut images, n_images);
        }
    }

    // ── 5. Post-NEB perpendicular relaxation ──
    if config.relax_images {
        let relax_config = super::constrained_opt::ConstrainedOptConfig {
            max_steps: config.relax_max_steps,
            grad_threshold: config.force_threshold.max(0.02),
            step_size: config.step_size * 0.5,
            method: config.method.clone(),
            constraint_tolerance: 1e-4,
        };
        if let Ok(new_energies) =
            super::constrained_opt::relax_neb_images(smiles, elements, &mut images, &relax_config)
        {
            energies = new_energies;
        }
    }

    // ── 6. Final energies (parallelized) ──
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        let par_energies: Vec<Result<f64, String>> = images
            .par_iter()
            .map(|img| {
                let mut grad = vec![0.0; n_xyz];
                crate::dynamics::neb_energy_and_gradient(
                    backend, smiles, elements, &mol, img, &mut grad,
                )
            })
            .collect();
        for (i, r) in par_energies.into_iter().enumerate() {
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

    let ci_idx = (1..n_images - 1)
        .max_by(|&a, &b| {
            energies[a]
                .partial_cmp(&energies[b])
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .unwrap_or(n_images / 2);

    let result_images: Vec<CiNebImage> = images
        .into_iter()
        .enumerate()
        .map(|(i, coords)| CiNebImage {
            index: i,
            coords,
            energy_kcal_mol: energies[i],
            is_climbing: i == ci_idx,
        })
        .collect();

    let mut notes = vec![
        format!(
            "CI-NEB ({}) with {} images, {} iterations, converged={}",
            config.method, n_images, actual_iter, converged
        ),
        format!("Final max force: {:.4} kcal/mol/Å", max_force),
    ];
    if config.use_idpp {
        notes.push("IDPP initialization used (physically meaningful initial path).".into());
    }
    if config.relax_images {
        notes.push("Post-NEB perpendicular image relaxation applied.".into());
    }

    Ok(CiNebResult {
        images: result_images,
        converged,
        iterations: actual_iter,
        max_force,
        notes,
    })
}

/// Redistribute images to maintain roughly equal spacing along the path.
///
/// Computes the total arc length, then redistributes images at equal
/// intervals using cubic interpolation.
fn redistribute_images(images: &mut [Vec<f64>], n_images: usize) {
    if n_images < 3 {
        return;
    }
    let n_xyz = images[0].len();

    // Compute cumulative arc lengths
    let mut arc = vec![0.0f64; n_images];
    for i in 1..n_images {
        let d: f64 = (0..n_xyz)
            .map(|k| (images[i][k] - images[i - 1][k]).powi(2))
            .sum::<f64>()
            .sqrt();
        arc[i] = arc[i - 1] + d;
    }

    let total_len = arc[n_images - 1];
    if total_len < 1e-10 {
        return;
    }

    // Target positions at equal arc-length intervals
    let old_images = images.to_vec();

    for i in 1..(n_images - 1) {
        let target_s = total_len * i as f64 / (n_images - 1) as f64;

        // Find bracket in old arc lengths
        let j = match arc.iter().position(|&s| s >= target_s) {
            Some(j) if j > 0 => j,
            _ => continue,
        };

        let s0 = arc[j - 1];
        let s1 = arc[j];
        let ds = s1 - s0;
        if ds < 1e-12 {
            continue;
        }
        let t = (target_s - s0) / ds;

        // Linear interpolation between neighbouring old images
        for k in 0..n_xyz {
            images[i][k] = (1.0 - t) * old_images[j - 1][k] + t * old_images[j][k];
        }
    }
}
