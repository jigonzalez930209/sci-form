//! Three-dimensional reference interaction site model (3D-RISM) solver core.
//!
//! This module solves the 3D Ornstein–Zernike equation with the
//! Kovalenko–Hirata closure on a periodic grid. The solute–solvent potential and
//! the solvent susceptibility are supplied by the caller, so the numerical core
//! can be tested independently of a particular solvent force field.
//!
//! # Equations
//!
//! In Fourier space, for a one-component solvent,
//!
//! ```text
//! h(k) = c(k) / (1 - χ(k) c(k)).
//! ```
//!
//! The KH closure is
//!
//! ```text
//! g(r) = exp(-βu(r) + h(r) - c(r))   when h - c <= 0
//! g(r) = 1 - βu(r) + h(r) - c(r)     otherwise
//! c(r) = g(r) - 1 - h(r).
//! ```
//!
//! # Units
//!
//! `βu` and `h`, `c`, `g` are dimensionless. The grid spacing follows from
//! `box_length_angstrom / n`. The bulk density is in Å⁻³ and the susceptibility
//! is the Fourier-space density response in Å³.
//!
//! # Scope
//!
//! The solver core is implemented. A solvent model that produces `χ(k)`, the
//! excess chemical potential and a hydration-free-energy benchmark are
//! follow-ups; the caller must supply a physically meaningful `χ(k)`.

use rustfft::num_complex::Complex;
use rustfft::FftPlanner;
use serde::{Deserialize, Serialize};

use crate::error::SciFormError;

/// Closure used by the solver.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RismClosure {
    /// Kovalenko–Hirata closure.
    Kh,
}

/// Simple solvent susceptibility models for tests and demonstrations.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum RismSolventModel {
    /// Constant `χ(k) = ρ`; the ideal-gas limit.
    IdealGas,
    /// Gaussian response `χ(k) = ρ exp(-k²ξ²/2)`.
    Gaussian {
        /// Correlation length in Å.
        correlation_length_angstrom: f64,
    },
}

/// Solver configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RismConfig {
    /// Grid points along x.
    pub nx: usize,
    /// Grid points along y.
    pub ny: usize,
    /// Grid points along z.
    pub nz: usize,
    /// Periodic box length in Å (cubic box).
    pub box_length_angstrom: f64,
    /// Temperature in K.
    pub temperature_k: f64,
    /// Bulk solvent density in Å⁻³.
    pub bulk_density_angstrom3: f64,
    /// Maximum Picard iterations.
    pub max_iterations: usize,
    /// Picard mixing fraction in `(0, 1]`.
    pub mixing: f64,
    /// Convergence tolerance on the maximum `|h_new - h|`.
    pub tolerance: f64,
    /// Closure.
    pub closure: RismClosure,
}

impl Default for RismConfig {
    fn default() -> Self {
        Self {
            nx: 32,
            ny: 32,
            nz: 32,
            box_length_angstrom: 16.0,
            temperature_k: 298.15,
            bulk_density_angstrom3: 0.0334,
            max_iterations: 500,
            mixing: 0.2,
            tolerance: 1e-6,
            closure: RismClosure::Kh,
        }
    }
}

impl RismConfig {
    fn n_grid(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    fn validate(&self) -> Result<(), SciFormError> {
        if self.nx < 2 || self.ny < 2 || self.nz < 2 {
            return Err(SciFormError::validation(
                "each RISM grid dimension must be >= 2",
            ));
        }
        if !self.box_length_angstrom.is_finite() || self.box_length_angstrom <= 0.0 {
            return Err(SciFormError::validation(
                "box_length_angstrom must be finite and > 0",
            ));
        }
        if !self.temperature_k.is_finite() || self.temperature_k <= 0.0 {
            return Err(SciFormError::validation(
                "temperature_k must be finite and > 0",
            ));
        }
        if !self.bulk_density_angstrom3.is_finite() || self.bulk_density_angstrom3 < 0.0 {
            return Err(SciFormError::validation(
                "bulk_density_angstrom3 must be finite and >= 0",
            ));
        }
        if self.max_iterations == 0 {
            return Err(SciFormError::validation("max_iterations must be >= 1"));
        }
        if !self.mixing.is_finite() || self.mixing <= 0.0 || self.mixing > 1.0 {
            return Err(SciFormError::validation(
                "mixing must be finite and in (0, 1]",
            ));
        }
        if !self.tolerance.is_finite() || self.tolerance <= 0.0 {
            return Err(SciFormError::validation("tolerance must be finite and > 0"));
        }
        Ok(())
    }
}

/// 3D-RISM solver result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RismResult {
    /// Solute–solvent pair distribution function `g(r)`.
    pub g: Vec<f64>,
    /// Total correlation function `h(r)`.
    pub h: Vec<f64>,
    /// Direct correlation function `c(r)`.
    pub c: Vec<f64>,
    /// Maximum `|h_new - h|` at the last iteration.
    pub residual: f64,
    /// Number of Picard iterations.
    pub iterations: usize,
    /// Whether the solver converged.
    pub converged: bool,
    /// Mean `g` over the grid.
    pub mean_g: f64,
    /// Maximum `g` over the grid (contact peak proxy).
    pub max_g: f64,
    /// Mean solute-solvent interaction energy in kcal/mol.
    pub mean_interaction_energy_kcal_mol: f64,
    /// Grid dimensions.
    pub shape: [usize; 3],
    /// Notes and caveats.
    pub notes: Vec<String>,
}

/// Solve the 3D-RISM equations on a periodic grid.
///
/// `solute_beta_potential` is `βu(r)` on the `nx × ny × nz` grid, in row-major
/// order. `solvent_susceptibility_k` is the Fourier-space solvent
/// susceptibility `χ(k)` on the same grid. For an ideal-gas solvent, use a
/// constant `bulk_density_angstrom3` for every k point.
pub fn solve_3d_rism(
    solute_beta_potential: &[f64],
    solvent_susceptibility_k: &[f64],
    config: &RismConfig,
) -> Result<RismResult, SciFormError> {
    config.validate()?;
    let n = config.n_grid();
    if solute_beta_potential.len() != n || solvent_susceptibility_k.len() != n {
        return Err(SciFormError::validation(format!(
            "grid arrays must contain {n} values (nx*ny*nz)"
        )));
    }
    if solute_beta_potential.iter().any(|value| !value.is_finite())
        || solvent_susceptibility_k
            .iter()
            .any(|value| !value.is_finite())
    {
        return Err(SciFormError::validation(
            "grid arrays must contain finite values",
        ));
    }

    let mut h = vec![0.0f64; n];
    let mut c = vec![0.0f64; n];
    let mut h_new = vec![0.0f64; n];
    let mut c_complex = vec![Complex::new(0.0, 0.0); n];
    let mut h_complex = vec![Complex::new(0.0, 0.0); n];
    let mut residual = f64::INFINITY;
    let mut converged = false;
    let mut iterations = 0;

    for iteration in 1..=config.max_iterations {
        iterations = iteration;
        for index in 0..n {
            let hc = h[index] - c[index];
            let g = match config.closure {
                RismClosure::Kh => {
                    if hc <= 0.0 {
                        (-solute_beta_potential[index] + hc).exp()
                    } else {
                        1.0 - solute_beta_potential[index] + hc
                    }
                }
            };
            c[index] = g - 1.0 - h[index];
            c_complex[index] = Complex::new(c[index], 0.0);
        }
        fft3d(&mut c_complex, config, false);
        for index in 0..n {
            let chi = solvent_susceptibility_k[index];
            let denominator = 1.0 - chi * c_complex[index].re;
            if denominator.abs() < 1e-12 {
                return Err(SciFormError::numerical(
                    "3D-RISM OZ denominator vanished; reduce the coupling or density",
                ));
            }
            // For a one-component susceptibility the OZ solution is scalar in
            // Fourier space; the imaginary part of c is retained.
            h_complex[index] = c_complex[index] / Complex::new(denominator, 0.0);
        }
        fft3d(&mut h_complex, config, true);
        residual = 0.0;
        for index in 0..n {
            let new_value = h_complex[index].re;
            residual = residual.max((new_value - h[index]).abs());
            h_new[index] = new_value;
        }
        for index in 0..n {
            h[index] = (1.0 - config.mixing) * h[index] + config.mixing * h_new[index];
        }
        if residual < config.tolerance {
            converged = true;
            break;
        }
    }

    let mut g = vec![0.0f64; n];
    for index in 0..n {
        g[index] = 1.0 + h[index];
    }
    if g.iter().any(|value| !value.is_finite()) {
        return Err(SciFormError::numerical(
            "3D-RISM produced a non-finite density",
        ));
    }
    let mean_g = g.iter().sum::<f64>() / n as f64;
    let max_g = g.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let cell_volume = config.box_length_angstrom.powi(3);
    let cell_volume_per_point = cell_volume / n as f64;
    let kbt_kcal_mol = 0.001_987_204_258_640_83 * config.temperature_k;
    let mean_interaction_energy_kcal_mol = config.bulk_density_angstrom3
        * kbt_kcal_mol
        * cell_volume_per_point
        * g.iter()
            .zip(solute_beta_potential)
            .map(|(density, beta_u)| density * beta_u)
            .sum::<f64>();

    Ok(RismResult {
        g,
        h,
        c,
        residual,
        iterations,
        converged,
        mean_g,
        max_g,
        mean_interaction_energy_kcal_mol,
        shape: [config.nx, config.ny, config.nz],
        notes: vec![
            "3D-RISM core with the Kovalenko–Hirata closure.".to_string(),
            "The solvent susceptibility is caller-supplied; no solvent force field is bundled."
                .to_string(),
            "The excess chemical potential and a hydration-free-energy benchmark are follow-ups."
                .to_string(),
        ],
    })
}

/// Build a solvent susceptibility grid for a simple model.
///
/// The grid uses the same periodic k-vector convention as the solver. The
/// Gaussian model is a numerical test response, not a parameterisation of a
/// real solvent.
pub fn solvent_susceptibility_k(
    config: &RismConfig,
    model: RismSolventModel,
) -> Result<Vec<f64>, SciFormError> {
    config.validate()?;
    if let RismSolventModel::Gaussian {
        correlation_length_angstrom,
    } = model
    {
        if !correlation_length_angstrom.is_finite() || correlation_length_angstrom <= 0.0 {
            return Err(SciFormError::validation(
                "Gaussian solvent correlation length must be finite and > 0",
            ));
        }
    }
    let n = config.n_grid();
    let mut values = Vec::with_capacity(n);
    let two_pi_over_l = 2.0 * std::f64::consts::PI / config.box_length_angstrom;
    for z in 0..config.nz {
        for y in 0..config.ny {
            for x in 0..config.nx {
                let kx = two_pi_over_l
                    * if x <= config.nx / 2 {
                        x as f64
                    } else {
                        x as f64 - config.nx as f64
                    };
                let ky = two_pi_over_l
                    * if y <= config.ny / 2 {
                        y as f64
                    } else {
                        y as f64 - config.ny as f64
                    };
                let kz = two_pi_over_l
                    * if z <= config.nz / 2 {
                        z as f64
                    } else {
                        z as f64 - config.nz as f64
                    };
                let k2 = kx * kx + ky * ky + kz * kz;
                let value = match model {
                    RismSolventModel::IdealGas => config.bulk_density_angstrom3,
                    RismSolventModel::Gaussian {
                        correlation_length_angstrom,
                    } => {
                        config.bulk_density_angstrom3
                            * (-0.5
                                * k2
                                * correlation_length_angstrom
                                * correlation_length_angstrom)
                                .exp()
                    }
                };
                values.push(value);
            }
        }
    }
    Ok(values)
}

fn fft3d(data: &mut [Complex<f64>], config: &RismConfig, inverse: bool) {
    let [nx, ny, nz] = [config.nx, config.ny, config.nz];
    let mut planner = FftPlanner::new();
    let fft_x = if inverse {
        planner.plan_fft_inverse(nx)
    } else {
        planner.plan_fft_forward(nx)
    };
    let fft_y = if inverse {
        planner.plan_fft_inverse(ny)
    } else {
        planner.plan_fft_forward(ny)
    };
    let fft_z = if inverse {
        planner.plan_fft_inverse(nz)
    } else {
        planner.plan_fft_forward(nz)
    };

    // x lines
    let mut line_x = vec![Complex::new(0.0, 0.0); nx];
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                line_x[x] = data[(z * ny + y) * nx + x];
            }
            fft_x.process(&mut line_x);
            for x in 0..nx {
                data[(z * ny + y) * nx + x] = line_x[x];
            }
        }
    }
    // y lines
    let mut line_y = vec![Complex::new(0.0, 0.0); ny];
    for z in 0..nz {
        for x in 0..nx {
            for y in 0..ny {
                line_y[y] = data[(z * ny + y) * nx + x];
            }
            fft_y.process(&mut line_y);
            for y in 0..ny {
                data[(z * ny + y) * nx + x] = line_y[y];
            }
        }
    }
    // z lines
    let mut line_z = vec![Complex::new(0.0, 0.0); nz];
    for y in 0..ny {
        for x in 0..nx {
            for z in 0..nz {
                line_z[z] = data[(z * ny + y) * nx + x];
            }
            fft_z.process(&mut line_z);
            for z in 0..nz {
                data[(z * ny + y) * nx + x] = line_z[z];
            }
        }
    }
    if inverse {
        let norm = 1.0 / (nx * ny * nz) as f64;
        for value in data.iter_mut() {
            *value *= norm;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ideal_gas_has_unit_density() {
        let config = RismConfig {
            nx: 8,
            ny: 8,
            nz: 8,
            max_iterations: 20,
            mixing: 0.5,
            tolerance: 1e-8,
            ..RismConfig::default()
        };
        let n = config.n_grid();
        let potential = vec![0.0; n];
        let chi = vec![config.bulk_density_angstrom3; n];
        let result = solve_3d_rism(&potential, &chi, &config).unwrap();
        assert!(result.converged);
        assert!(
            (result.mean_g - 1.0).abs() < 1e-8,
            "mean g = {}",
            result.mean_g
        );
    }

    #[test]
    fn attractive_potential_increases_the_contact_peak() {
        let config = RismConfig {
            nx: 8,
            ny: 8,
            nz: 8,
            max_iterations: 100,
            mixing: 0.3,
            ..RismConfig::default()
        };
        let n = config.n_grid();
        let potential: Vec<f64> = (0..n)
            .map(|index| {
                let x = (index % 8) as f64 / 8.0;
                let y = ((index / 8) % 8) as f64 / 8.0;
                let z = (index / 64) as f64 / 8.0;
                let r2 = (x - 0.5).powi(2) + (y - 0.5).powi(2) + (z - 0.5).powi(2);
                -2.0 * (-r2 / 0.02).exp()
            })
            .collect();
        let chi = vec![0.0; n];
        let result = solve_3d_rism(&potential, &chi, &config).unwrap();
        assert!(result.max_g > 1.0, "max g = {}", result.max_g);
    }

    #[test]
    fn non_cubic_grid_runs() {
        let config = RismConfig {
            nx: 4,
            ny: 6,
            nz: 8,
            max_iterations: 20,
            ..RismConfig::default()
        };
        let n = config.n_grid();
        let result = solve_3d_rism(
            &vec![0.0; n],
            &vec![config.bulk_density_angstrom3; n],
            &config,
        )
        .unwrap();
        assert_eq!(result.shape, [4, 6, 8]);
        assert!((result.mean_g - 1.0).abs() < 1e-8);
    }

    #[test]
    fn ideal_gas_susceptibility_is_constant() {
        let config = RismConfig::default();
        let values = solvent_susceptibility_k(&config, RismSolventModel::IdealGas).unwrap();
        assert!(values
            .iter()
            .all(|value| (*value - config.bulk_density_angstrom3).abs() < 1e-15));
    }

    #[test]
    fn gaussian_susceptibility_decays() {
        let config = RismConfig::default();
        let values = solvent_susceptibility_k(
            &config,
            RismSolventModel::Gaussian {
                correlation_length_angstrom: 3.0,
            },
        )
        .unwrap();
        assert!(values[0] > values[values.len() - 1]);
    }

    #[test]
    fn interaction_energy_is_finite() {
        let config = RismConfig {
            nx: 8,
            ny: 8,
            nz: 8,
            max_iterations: 20,
            ..RismConfig::default()
        };
        let n = config.n_grid();
        let potential: Vec<f64> = (0..n).map(|index| -0.1 * (index % 5) as f64).collect();
        let chi = solvent_susceptibility_k(&config, RismSolventModel::IdealGas).unwrap();
        let result = solve_3d_rism(&potential, &chi, &config).unwrap();
        assert!(result.mean_interaction_energy_kcal_mol.is_finite());
    }

    #[test]
    fn wrong_grid_length_is_rejected() {
        let config = RismConfig::default();
        assert!(solve_3d_rism(&[0.0], &[0.0], &config).is_err());
    }
}
