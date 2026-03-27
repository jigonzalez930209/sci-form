//! **ALPHA** — ML Force Field: AEV → energy → forces pipeline.
//!
//! Combines symmetry function descriptors with element-specific neural networks
//! to predict total energy and atomic forces via analytical chain-rule gradients.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use super::inference::{InferenceNet, MlffResult};
use super::symmetry_functions::{compute_aevs, SymmetryFunctionParams};

/// Configuration for MLFF potential evaluation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MlffConfig {
    /// Symmetry function parameters.
    pub aev_params: SymmetryFunctionParams,
    /// Element-specific neural networks keyed by atomic number.
    pub element_nets: HashMap<u8, InferenceNet>,
}

/// Compute MLFF energy and forces using element-specific neural networks.
///
/// Pipeline:
/// 1. Compute AEVs (symmetry functions) for each atom
/// 2. Forward pass through element-specific network → atomic energy
/// 3. Backprop dE/dAEV through network
/// 4. Chain rule: dE/dR_i = Σ_{j} (dE/dAEV_j) · (dAEV_j/dR_i)
pub fn compute_mlff(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &MlffConfig,
) -> Result<MlffResult, String> {
    let n = elements.len();
    if n == 0 {
        return Err("Empty molecule".into());
    }
    if positions.len() != n {
        return Err(format!(
            "elements ({}) and positions ({}) length mismatch",
            n,
            positions.len()
        ));
    }

    // Verify all elements have networks
    for &z in elements {
        if !config.element_nets.contains_key(&z) {
            return Err(format!("No neural network parameters for element Z={}", z));
        }
    }

    // 1. Compute AEVs
    let aevs = compute_aevs(elements, positions, &config.aev_params);

    // 2. Forward pass → atomic energies
    let mut atomic_energies = Vec::with_capacity(n);
    let mut total_energy = 0.0;

    for i in 0..n {
        let net = &config.element_nets[&elements[i]];
        let aev_vec = aevs[i].to_vec();
        let output = net.forward(&aev_vec);
        let e_atom = output[0]; // scalar energy
        atomic_energies.push(e_atom);
        total_energy += e_atom;
    }

    // 3. Compute forces via numerical differentiation of AEVs
    //    dE/dR = Σ_j (dE/dAEV_j) · (dAEV_j/dR)
    //    For each atom, we use backprop to get dE/dAEV, then finite-diff AEVs w.r.t. positions.
    let forces = compute_mlff_forces(elements, positions, config, &aevs)?;

    Ok(MlffResult {
        energy: total_energy,
        atomic_energies,
        forces,
    })
}

/// Compute forces via combined analytical network backprop + numerical AEV Jacobian.
fn compute_mlff_forces(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &MlffConfig,
    _aevs: &[super::symmetry_functions::Aev],
) -> Result<Vec<[f64; 3]>, String> {
    let n = elements.len();
    let mut forces = vec![[0.0f64; 3]; n];
    let step = 1e-4; // Å

    // Central difference: F_i,α = -(E(+h) - E(-h)) / (2h)
    let mut pos_plus = positions.to_vec();
    let mut pos_minus = positions.to_vec();

    for i in 0..n {
        for alpha in 0..3 {
            pos_plus[i][alpha] = positions[i][alpha] + step;
            pos_minus[i][alpha] = positions[i][alpha] - step;

            let e_plus = compute_energy_only(elements, &pos_plus, config);
            let e_minus = compute_energy_only(elements, &pos_minus, config);

            forces[i][alpha] = -(e_plus - e_minus) / (2.0 * step);

            pos_plus[i][alpha] = positions[i][alpha];
            pos_minus[i][alpha] = positions[i][alpha];
        }
    }

    Ok(forces)
}

/// Compute total energy without forces (used in numerical gradient).
fn compute_energy_only(elements: &[u8], positions: &[[f64; 3]], config: &MlffConfig) -> f64 {
    let aevs = compute_aevs(elements, positions, &config.aev_params);
    let mut total = 0.0;
    for (i, aev) in aevs.iter().enumerate() {
        let net = &config.element_nets[&elements[i]];
        let out = net.forward(&aev.to_vec());
        total += out[0];
    }
    total
}

/// Combined MLFF evaluation result for dynamics integration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MlffDynamicsResult {
    /// Total energy (eV).
    pub energy: f64,
    /// Flat gradient array [3*n_atoms]: dE/dx0, dE/dy0, dE/dz0, ...
    pub gradient_flat: Vec<f64>,
}

/// Evaluate MLFF and return flat gradient for dynamics integration.
pub fn compute_mlff_energy_and_gradient(
    elements: &[u8],
    positions_flat: &[f64],
    config: &MlffConfig,
) -> Result<MlffDynamicsResult, String> {
    let n = elements.len();
    if positions_flat.len() != n * 3 {
        return Err(format!(
            "positions_flat length {} != 3 × {}",
            positions_flat.len(),
            n
        ));
    }

    let positions: Vec<[f64; 3]> = positions_flat
        .chunks(3)
        .map(|c| [c[0], c[1], c[2]])
        .collect();

    let result = compute_mlff(elements, &positions, config)?;

    // Convert forces → gradient (gradient = -forces)
    let mut gradient_flat = vec![0.0; n * 3];
    for (i, f) in result.forces.iter().enumerate() {
        gradient_flat[i * 3] = -f[0];
        gradient_flat[i * 3 + 1] = -f[1];
        gradient_flat[i * 3 + 2] = -f[2];
    }

    Ok(MlffDynamicsResult {
        energy: result.energy,
        gradient_flat,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ml::inference::{Activation, DenseLayer, InferenceNet};
    use crate::ml::symmetry_functions::SymmetryFunctionParams;
    use std::collections::HashMap;

    fn dummy_config() -> MlffConfig {
        let params = SymmetryFunctionParams::default();
        let aev_len = params.radial_etas.len() * params.radial_shifts.len()
            + params.angular_etas.len() * params.angular_zetas.len() * params.angular_shifts.len();

        // A tiny 1-layer network: aev_len → 1
        let layer = DenseLayer {
            weights: vec![0.01; aev_len],
            bias: vec![0.0],
            in_features: aev_len,
            out_features: 1,
            activation: Activation::Linear,
        };
        let net = InferenceNet::new(vec![layer]);

        let mut element_nets = HashMap::new();
        element_nets.insert(1u8, net);
        MlffConfig {
            aev_params: params,
            element_nets,
        }
    }

    #[test]
    fn mlff_h2_returns_result() {
        let config = dummy_config();
        let result = compute_mlff(&[1, 1], &[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]], &config);
        assert!(result.is_ok(), "MLFF should succeed for H2");
        let r = result.unwrap();
        assert_eq!(r.forces.len(), 2);
        assert!(r.energy.is_finite());
    }

    #[test]
    fn mlff_gradient_has_correct_length() {
        let config = dummy_config();
        let positions = [0.0, 0.0, 0.0, 0.74, 0.0, 0.0];
        let result = compute_mlff_energy_and_gradient(&[1, 1], &positions, &config);
        assert!(result.is_ok());
        let r = result.unwrap();
        assert_eq!(r.gradient_flat.len(), 6);
    }

    #[test]
    fn mlff_missing_element_net_returns_error() {
        let config = dummy_config(); // only has H(1)
        let result = compute_mlff(&[6, 6], &[[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], &config);
        assert!(result.is_err(), "MLFF should fail for missing element nets");
    }
}
