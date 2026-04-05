//! Multi-angular approach sampling on SO(3).
//!
//! Generates candidate approach orientations uniformly sampled on the
//! rotation group, then scores each with a short NEB probe to find the
//! lowest-energy approach direction.

/// Result of angular sampling.
pub struct SamplingResult {
    /// Best approach direction found.
    pub best_direction: [f64; 3],
    /// Energy barrier for the best direction (eV).
    pub best_barrier: f64,
    /// All sampled directions with their barriers, sorted ascending.
    pub candidates: Vec<SamplingCandidate>,
}

/// A single candidate approach.
#[derive(Debug, Clone)]
pub struct SamplingCandidate {
    /// Direction unit vector.
    pub direction: [f64; 3],
    /// Estimated barrier from short NEB probe (eV).
    pub barrier: f64,
    /// Whether the probe converged.
    pub converged: bool,
}

/// Sample approach orientations using quasi-uniform quaternion sampling.
///
/// Generates `n_samples` orientations on the unit sphere, positions fragment B
/// around fragment A at each orientation, then runs a short NEB probe (few images)
/// to estimate the barrier. Returns the candidate with the lowest barrier.
pub fn sample_approach_orientations(
    smiles: &str,
    elements_a: &[u8],
    coords_a: &[f64],
    elements_b: &[u8],
    coords_b: &[f64],
    reactive_dist: f64,
    n_samples: usize,
    method: &str,
) -> Result<SamplingResult, String> {
    let directions = fibonacci_sphere(n_samples);

    // Parallelise NEB probes — each direction is independent
    #[cfg(feature = "parallel")]
    let mut candidates: Vec<SamplingCandidate> = {
        use rayon::prelude::*;
        directions
            .par_iter()
            .map(|dir| {
                let barrier = probe_neb_barrier(
                    smiles,
                    elements_a,
                    coords_a,
                    elements_b,
                    coords_b,
                    *dir,
                    reactive_dist,
                    method,
                );
                let (b, conv) = match barrier {
                    Ok(b) => (b, true),
                    Err(_) => (f64::INFINITY, false),
                };
                SamplingCandidate {
                    direction: *dir,
                    barrier: b,
                    converged: conv,
                }
            })
            .collect()
    };
    #[cfg(not(feature = "parallel"))]
    let mut candidates: Vec<SamplingCandidate> = {
        directions
            .iter()
            .map(|dir| {
                let barrier = probe_neb_barrier(
                    smiles,
                    elements_a,
                    coords_a,
                    elements_b,
                    coords_b,
                    *dir,
                    reactive_dist,
                    method,
                );
                let (b, conv) = match barrier {
                    Ok(b) => (b, true),
                    Err(_) => (f64::INFINITY, false),
                };
                SamplingCandidate {
                    direction: *dir,
                    barrier: b,
                    converged: conv,
                }
            })
            .collect()
    };

    // Sort by barrier (ascending)
    candidates.sort_by(|a, b| {
        a.barrier
            .partial_cmp(&b.barrier)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let best = candidates.first().ok_or("No candidates generated")?;
    Ok(SamplingResult {
        best_direction: best.direction,
        best_barrier: best.barrier,
        candidates,
    })
}

/// Generate approximately uniform points on a sphere using the Fibonacci spiral.
fn fibonacci_sphere(n: usize) -> Vec<[f64; 3]> {
    let golden_ratio = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let mut points = Vec::with_capacity(n);

    for i in 0..n {
        let theta = 2.0 * std::f64::consts::PI * (i as f64) / golden_ratio;
        let phi = (1.0 - 2.0 * (i as f64 + 0.5) / (n as f64)).acos();

        let x = phi.sin() * theta.cos();
        let y = phi.sin() * theta.sin();
        let z = phi.cos();
        points.push([x, y, z]);
    }

    points
}

/// Run a short NEB probe (5 images, few iterations) to estimate the barrier.
fn probe_neb_barrier(
    smiles: &str,
    elements_a: &[u8],
    coords_a: &[f64],
    elements_b: &[u8],
    coords_b: &[f64],
    direction: [f64; 3],
    reactive_dist: f64,
    method: &str,
) -> Result<f64, String> {
    let n_a = elements_a.len();
    let n_b = elements_b.len();
    let n_total = n_a + n_b;

    // Build start configuration: fragments far apart along direction
    let far_dist = reactive_dist + 4.0;
    let mut start_coords = Vec::with_capacity(n_total * 3);

    // Fragment A at origin (centred)
    let mut ca = coords_a.to_vec();
    centre_at_origin(&mut ca);
    start_coords.extend_from_slice(&ca);

    // Fragment B offset far along direction
    let mut cb = coords_b.to_vec();
    centre_at_origin(&mut cb);
    for k in 0..n_b {
        start_coords.push(cb[k * 3] + far_dist * direction[0]);
        start_coords.push(cb[k * 3 + 1] + far_dist * direction[1]);
        start_coords.push(cb[k * 3 + 2] + far_dist * direction[2]);
    }

    // Build end configuration: fragments at reactive distance apart
    let mut end_coords = Vec::with_capacity(n_total * 3);
    end_coords.extend_from_slice(&ca);
    for k in 0..n_b {
        end_coords.push(cb[k * 3] + reactive_dist * direction[0]);
        end_coords.push(cb[k * 3 + 1] + reactive_dist * direction[1]);
        end_coords.push(cb[k * 3 + 2] + reactive_dist * direction[2]);
    }

    // Combined elements
    let mut elements = Vec::with_capacity(n_total);
    elements.extend_from_slice(elements_a);
    elements.extend_from_slice(elements_b);

    // Build linearly interpolated images (5 images)
    let n_images = 5;
    let n_xyz = n_total * 3;
    let mut images = Vec::with_capacity(n_images);
    for img in 0..n_images {
        let t = img as f64 / (n_images - 1) as f64;
        let mut frame = Vec::with_capacity(n_xyz);
        for k in 0..n_xyz {
            frame.push(start_coords[k] * (1.0 - t) + end_coords[k] * t);
        }
        images.push(frame);
    }

    // Evaluate energy at each image (parallelized)
    let mol = crate::graph::Molecule::from_smiles(smiles)?;
    let backend = crate::dynamics::NebBackend::from_method(method)?;

    #[cfg(feature = "parallel")]
    let energies: Vec<f64> = {
        use rayon::prelude::*;
        let par_results: Vec<Result<f64, String>> = images
            .par_iter()
            .map(|image| {
                let mut grad_buf = vec![0.0f64; n_xyz];
                crate::dynamics::neb_energy_and_gradient(
                    backend,
                    smiles,
                    &elements,
                    &mol,
                    image,
                    &mut grad_buf,
                )
            })
            .collect();
        let mut e = Vec::with_capacity(n_images);
        for r in par_results {
            e.push(r?);
        }
        e
    };
    #[cfg(not(feature = "parallel"))]
    let energies: Vec<f64> = {
        let mut e = Vec::with_capacity(n_images);
        let mut grad_buf = vec![0.0f64; n_xyz];
        for image in &images {
            let energy = crate::dynamics::neb_energy_and_gradient(
                backend,
                smiles,
                &elements,
                &mol,
                image,
                &mut grad_buf,
            )?;
            e.push(energy);
        }
        e
    };

    // Barrier = max energy - start energy
    let e_start = energies[0];
    let e_max = energies.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    Ok(e_max - e_start)
}

fn centre_at_origin(coords: &mut [f64]) {
    let n = coords.len() / 3;
    if n == 0 {
        return;
    }
    let mut c = [0.0; 3];
    for i in 0..n {
        c[0] += coords[i * 3];
        c[1] += coords[i * 3 + 1];
        c[2] += coords[i * 3 + 2];
    }
    let nf = n as f64;
    for i in 0..n {
        coords[i * 3] -= c[0] / nf;
        coords[i * 3 + 1] -= c[1] / nf;
        coords[i * 3 + 2] -= c[2] / nf;
    }
}
