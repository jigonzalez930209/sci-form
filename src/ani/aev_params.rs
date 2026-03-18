//! ANI-2x symmetry function parameter tables.
//!
//! Defines the radial (η, Rs) and angular (η, Rs, ζ, θs) parameters
//! that characterize the Behler-Parrinello atomic environment vectors.

use serde::{Deserialize, Serialize};

/// Supported chemical elements for ANI potentials.
pub const ANI_ELEMENTS: [u8; 7] = [1, 6, 7, 8, 9, 16, 17]; // H,C,N,O,F,S,Cl

/// Number of supported element types.
pub const N_SPECIES: usize = ANI_ELEMENTS.len();

/// AEV hyperparameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AevParams {
    /// Cutoff radius for radial terms (Å).
    pub radial_cutoff: f64,
    /// Cutoff radius for angular terms (Å).
    pub angular_cutoff: f64,
    /// η values for radial symmetry functions.
    pub radial_eta: Vec<f64>,
    /// Shifted center Rs values for radial symmetry functions.
    pub radial_rs: Vec<f64>,
    /// η values for angular symmetry functions.
    pub angular_eta: Vec<f64>,
    /// Shifted center Rs values for angular symmetry functions.
    pub angular_rs: Vec<f64>,
    /// ζ values for angular symmetry functions.
    pub angular_zeta: Vec<f64>,
    /// θs angular shift values (radians).
    pub angular_theta_s: Vec<f64>,
}

impl AevParams {
    /// Dimension of the radial sub-vector per species pair.
    pub fn radial_length(&self) -> usize {
        self.radial_eta.len() * self.radial_rs.len()
    }

    /// Dimension of the angular sub-vector per species-pair combination.
    pub fn angular_length(&self) -> usize {
        self.angular_eta.len() * self.angular_rs.len()
            * self.angular_zeta.len() * self.angular_theta_s.len()
    }

    /// Total AEV length for one atom.
    pub fn total_aev_length(&self) -> usize {
        // Radial: N_SPECIES components × radial_length
        // Angular: N_SPECIES*(N_SPECIES+1)/2 components × angular_length
        let n_rad = N_SPECIES * self.radial_length();
        let n_ang = N_SPECIES * (N_SPECIES + 1) / 2 * self.angular_length();
        n_rad + n_ang
    }
}

/// Map atomic number to species index for AEV computation.
pub fn species_index(z: u8) -> Option<usize> {
    ANI_ELEMENTS.iter().position(|&e| e == z)
}

/// Default ANI-2x-like parameters.
///
/// These are representative values inspired by the ANI-2x publication:
/// Devereux et al., JCTC 16 (2020) 4192.
pub fn default_ani2x_params() -> AevParams {
    use std::f64::consts::PI;

    // 8 radial functions
    let radial_eta = vec![19.7; 8];
    let radial_rs: Vec<f64> = (0..8).map(|i| 0.8 + 0.5625 * i as f64).collect();

    // 4 angular shells
    let angular_eta = vec![12.5; 4];
    let angular_rs: Vec<f64> = (0..4).map(|i| 0.8 + 0.95 * i as f64).collect();
    let angular_zeta = vec![14.1; 1];
    let angular_theta_s: Vec<f64> = (0..8)
        .map(|i| PI * i as f64 / 8.0)
        .collect();

    AevParams {
        radial_cutoff: 5.2,
        angular_cutoff: 3.5,
        radial_eta,
        radial_rs,
        angular_eta,
        angular_rs,
        angular_zeta,
        angular_theta_s,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_species_index() {
        assert_eq!(species_index(1), Some(0));  // H
        assert_eq!(species_index(6), Some(1));  // C
        assert_eq!(species_index(8), Some(3));  // O
        assert_eq!(species_index(26), None);    // Fe not supported
    }

    #[test]
    fn test_aev_dimensions() {
        let params = default_ani2x_params();
        assert!(params.radial_length() > 0);
        assert!(params.angular_length() > 0);
        assert!(params.total_aev_length() > 0);
    }
}
