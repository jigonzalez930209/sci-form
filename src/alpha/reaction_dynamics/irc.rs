//! Intrinsic Reaction Coordinate (IRC).
//!
//! Follows the minimum energy path from a transition state toward
//! reactants and products using mass-weighted steepest descent.

/// IRC result for one direction (forward or reverse).
#[derive(Debug, Clone)]
pub struct IrcPathSegment {
    /// Sequence of coordinate frames along the IRC.
    pub frames: Vec<Vec<f64>>,
    /// Energies at each frame (eV or kcal/mol depending on backend).
    pub energies: Vec<f64>,
    /// Arc-length along the path (mass-weighted).
    pub arc_lengths: Vec<f64>,
    /// Whether the path converged to a minimum.
    pub converged: bool,
}

/// Full IRC result (forward + reverse).
#[derive(Debug, Clone)]
pub struct IrcResult {
    /// Forward path (TS → products).
    pub forward: IrcPathSegment,
    /// Reverse path (TS → reactants).
    pub reverse: IrcPathSegment,
}

/// IRC configuration.
pub struct IrcConfig {
    /// Step size in mass-weighted coordinates (amu^(1/2) Å).
    pub step_size: f64,
    /// Maximum number of steps per direction.
    pub max_steps: usize,
    /// Convergence threshold on gradient norm.
    pub grad_threshold: f64,
    /// NEB backend for energy/gradient evaluation.
    pub method: String,
}

impl Default for IrcConfig {
    fn default() -> Self {
        Self {
            step_size: 0.05,
            max_steps: 100,
            grad_threshold: 0.005,
            method: "uff".into(),
        }
    }
}

/// Compute IRC from a transition state in both directions.
///
/// Uses the numerical Hessian at the TS to determine the imaginary-frequency
/// mode, then follows steepest descent in mass-weighted coordinates.
pub fn compute_irc(
    smiles: &str,
    elements: &[u8],
    ts_coords: &[f64],
    config: &IrcConfig,
) -> Result<IrcResult, String> {
    let _n_xyz = ts_coords.len();
    let mol = crate::graph::Molecule::from_smiles(smiles)?;
    let backend = crate::dynamics::NebBackend::from_method(&config.method)?;

    // 1. Compute numerical Hessian at the TS
    let ts_mode = compute_ts_mode(backend, smiles, elements, &mol, ts_coords)?;

    // 2. Follow IRC in the forward direction (+ mode displacement)
    let forward = follow_irc_direction(
        backend, smiles, elements, &mol, ts_coords, &ts_mode, 1.0, config,
    )?;

    // 3. Follow IRC in the reverse direction (– mode displacement)
    let reverse = follow_irc_direction(
        backend, smiles, elements, &mol, ts_coords, &ts_mode, -1.0, config,
    )?;

    Ok(IrcResult { forward, reverse })
}

/// Compute the transition-state mode (imaginary frequency eigenvector).
///
/// Uses a numerical Hessian (central finite differences) and finds the
/// eigenvector corresponding to the lowest (most negative) eigenvalue.
fn compute_ts_mode(
    backend: crate::dynamics::NebBackend,
    smiles: &str,
    elements: &[u8],
    mol: &crate::graph::Molecule,
    coords: &[f64],
) -> Result<Vec<f64>, String> {
    let n = coords.len();
    let h = 0.005; // finite difference step (Å)

    // Compute reference gradient
    let mut g0 = vec![0.0f64; n];
    crate::dynamics::neb_energy_and_gradient(backend, smiles, elements, mol, coords, &mut g0)?;

    // Compute Hessian rows in parallel — each displacement is independent
    #[cfg(feature = "parallel")]
    let hessian_rows: Vec<Result<Vec<f64>, String>> = {
        use rayon::prelude::*;
        (0..n)
            .into_par_iter()
            .map(|i| {
                let mut x_plus = coords.to_vec();
                x_plus[i] += h;
                let mut g_plus = vec![0.0f64; n];
                crate::dynamics::neb_energy_and_gradient(
                    backend,
                    smiles,
                    elements,
                    mol,
                    &x_plus,
                    &mut g_plus,
                )?;
                // Central difference: H_ij ≈ (g_j(x+h_i) - g_j(x)) / h
                let row: Vec<f64> = (0..n).map(|j| (g_plus[j] - g0[j]) / h).collect();
                Ok(row)
            })
            .collect()
    };
    #[cfg(not(feature = "parallel"))]
    let hessian_rows: Vec<Result<Vec<f64>, String>> = {
        (0..n)
            .map(|i| {
                let mut x_plus = coords.to_vec();
                x_plus[i] += h;
                let mut g_plus = vec![0.0f64; n];
                crate::dynamics::neb_energy_and_gradient(
                    backend,
                    smiles,
                    elements,
                    mol,
                    &x_plus,
                    &mut g_plus,
                )?;
                let row: Vec<f64> = (0..n).map(|j| (g_plus[j] - g0[j]) / h).collect();
                Ok(row)
            })
            .collect()
    };

    let mut hessian = vec![0.0f64; n * n];
    for (i, row_result) in hessian_rows.into_iter().enumerate() {
        let row = row_result?;
        for j in 0..n {
            hessian[i * n + j] = row[j];
        }
    }

    // Symmetrise
    for i in 0..n {
        for j in (i + 1)..n {
            let avg = 0.5 * (hessian[i * n + j] + hessian[j * n + i]);
            hessian[i * n + j] = avg;
            hessian[j * n + i] = avg;
        }
    }

    // Mass-weight the Hessian
    let masses = atomic_masses(elements);
    for i in 0..n {
        let mi = masses[i / 3].sqrt();
        for j in 0..n {
            let mj = masses[j / 3].sqrt();
            hessian[i * n + j] /= mi * mj;
        }
    }

    // Find the eigenvector with the lowest eigenvalue
    // Use power iteration on -H to find the largest eigenvalue of -H = most negative of H
    let mut v = vec![0.0f64; n];
    // Initialise with a random-ish direction
    for i in 0..n {
        v[i] = ((i as f64 + 1.0) * 0.618).sin();
    }
    normalise(&mut v);

    let neg_hessian: Vec<f64> = hessian.iter().map(|&x| -x).collect();

    for _ in 0..200 {
        let mut new_v = vec![0.0f64; n];
        for i in 0..n {
            for j in 0..n {
                new_v[i] += neg_hessian[i * n + j] * v[j];
            }
        }
        normalise(&mut new_v);
        v = new_v;
    }

    // Un-mass-weight the eigenvector
    for i in 0..n {
        v[i] /= masses[i / 3].sqrt();
    }
    normalise(&mut v);

    Ok(v)
}

/// Follow IRC in one direction using mass-weighted steepest descent.
fn follow_irc_direction(
    backend: crate::dynamics::NebBackend,
    smiles: &str,
    elements: &[u8],
    mol: &crate::graph::Molecule,
    ts_coords: &[f64],
    mode: &[f64],
    sign: f64,
    config: &IrcConfig,
) -> Result<IrcPathSegment, String> {
    let n = ts_coords.len();
    let masses = atomic_masses(elements);

    // Initial displacement along the TS mode
    let mut x: Vec<f64> = ts_coords
        .iter()
        .enumerate()
        .map(|(i, &c)| c + sign * config.step_size * mode[i])
        .collect();

    let mut frames = vec![ts_coords.to_vec(), x.clone()];
    let mut energies = Vec::new();
    let mut arc_lengths = vec![0.0];
    let mut total_arc = 0.0f64;

    // Get TS energy
    let mut g0 = vec![0.0f64; n];
    let e_ts = crate::dynamics::neb_energy_and_gradient(
        backend, smiles, elements, mol, ts_coords, &mut g0,
    )?;
    energies.push(e_ts);

    // Get first displaced frame energy
    let mut grad = vec![0.0f64; n];
    let e0 =
        crate::dynamics::neb_energy_and_gradient(backend, smiles, elements, mol, &x, &mut grad)?;
    energies.push(e0);

    // Compute initial arc length
    let ds: f64 = (0..n)
        .map(|i| {
            let dx = x[i] - ts_coords[i];
            masses[i / 3] * dx * dx
        })
        .sum::<f64>()
        .sqrt();
    total_arc += ds;
    arc_lengths.push(total_arc);

    let mut converged = false;

    for _ in 0..config.max_steps {
        let mut grad = vec![0.0f64; n];
        let energy = crate::dynamics::neb_energy_and_gradient(
            backend, smiles, elements, mol, &x, &mut grad,
        )?;

        // Mass-weighted gradient norm
        let mw_gnorm: f64 = grad
            .iter()
            .enumerate()
            .map(|(i, &g)| g * g / masses[i / 3])
            .sum::<f64>()
            .sqrt();

        if mw_gnorm < config.grad_threshold {
            converged = true;
            break;
        }

        // Steepest descent step in mass-weighted coordinates
        let mut x_new = vec![0.0f64; n];
        for i in 0..n {
            x_new[i] = x[i] - config.step_size * grad[i] / (masses[i / 3] * mw_gnorm);
        }

        // Check for NaN
        if !x_new.iter().all(|v| v.is_finite()) {
            break;
        }

        let ds: f64 = (0..n)
            .map(|i| {
                let dx = x_new[i] - x[i];
                masses[i / 3] * dx * dx
            })
            .sum::<f64>()
            .sqrt();
        total_arc += ds;

        x = x_new;
        frames.push(x.clone());
        energies.push(energy);
        arc_lengths.push(total_arc);

        // Check if energy is rising (we passed through the minimum)
        if energies.len() >= 3 {
            let last = energies.len() - 1;
            if energies[last] > energies[last - 1] && energies[last - 1] > energies[last - 2] {
                converged = true;
                break;
            }
        }
    }

    Ok(IrcPathSegment {
        frames,
        energies,
        arc_lengths,
        converged,
    })
}

fn normalise(v: &mut [f64]) {
    let len: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    if len > 1e-14 {
        for x in v.iter_mut() {
            *x /= len;
        }
    }
}

fn atomic_masses(elements: &[u8]) -> Vec<f64> {
    elements
        .iter()
        .map(|&z| match z {
            1 => 1.008,
            6 => 12.011,
            7 => 14.007,
            8 => 15.999,
            9 => 18.998,
            15 => 30.974,
            16 => 32.06,
            17 => 35.45,
            35 => 79.904,
            53 => 126.904,
            _ => z as f64 * 2.0, // rough estimate
        })
        .collect()
}
