//! Public API for ANI machine-learning potentials.
//!
//! Provides the top-level functions `compute_ani()` for energy evaluation
//! and force computation from molecular geometries.

use super::aev::compute_aevs;
use super::aev_params::{default_ani2x_params, species_index};
use super::gradients::compute_forces;
use super::neighbor::CellList;
use super::nn::FeedForwardNet;
use nalgebra::DVector;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Configuration for ANI computation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AniConfig {
    /// Cutoff radius (Å). Default: 5.2.
    pub cutoff: f64,
    /// Whether to compute forces.
    pub compute_forces: bool,
    /// Whether to return the computed AEVs for validation.
    pub output_aevs: bool,
}

impl Default for AniConfig {
    fn default() -> Self {
        AniConfig {
            cutoff: 5.2,
            compute_forces: true,
            output_aevs: false,
        }
    }
}

/// Result of an ANI energy/force calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AniResult {
    /// Total potential energy (Hartree).
    pub energy: f64,
    /// Atomic forces [x,y,z] per atom (Hartree/Å). Empty if not requested.
    pub forces: Vec<[f64; 3]>,
    /// Atomic species (atomic numbers).
    pub species: Vec<u8>,
    /// Per-atom energy contributions.
    pub atomic_energies: Vec<f64>,
    /// Optional AEVs for each atom.
    pub aevs: Option<Vec<Vec<f64>>>,
}

/// Compute ANI energy (and optionally forces) for a molecular geometry.
///
/// `elements`: atomic numbers.
/// `positions`: [x,y,z] coordinates in Ångström.
/// `config`: computation parameters.
/// `models`: pre-loaded element→network map (from `weights::load_weights`).
pub fn compute_ani(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &AniConfig,
    models: &HashMap<u8, FeedForwardNet>,
) -> Result<AniResult, String> {
    if elements.len() != positions.len() {
        return Err(format!(
            "elements ({}) and positions ({}) length mismatch",
            elements.len(),
            positions.len()
        ));
    }

    // Validate species support
    for &z in elements {
        if species_index(z).is_none() {
            return Err(format!("Unsupported element Z={z} for ANI potential"));
        }
        if !models.contains_key(&z) {
            return Err(format!("No model weights for element Z={z}"));
        }
    }

    let params = default_ani2x_params();

    // Build neighbor list
    let cell_list = CellList::new(positions, config.cutoff);
    let neighbors = cell_list.find_neighbors(positions);

    // Compute AEVs
    let aevs = compute_aevs(elements, positions, &neighbors, &params);

    // Neural network inference: per-atom energies
    let mut atomic_energies = Vec::with_capacity(elements.len());
    for (i, &z) in elements.iter().enumerate() {
        let net = &models[&z];
        let input = DVector::from_vec(aevs[i].clone());
        let e_atom = net.forward(&input);
        atomic_energies.push(e_atom);
    }

    let energy: f64 = atomic_energies.iter().sum();

    // Forces
    let forces = if config.compute_forces {
        compute_forces(elements, positions, &neighbors, &params, models)
    } else {
        Vec::new()
    };

    Ok(AniResult {
        energy,
        forces,
        species: elements.to_vec(),
        atomic_energies,
        aevs: if config.output_aevs { Some(aevs) } else { None },
    })
}

/// Compute ANI energy using internally-generated test weights.
///
/// This is for testing and demonstration only — not physically meaningful.
pub fn compute_ani_test(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<AniResult, String> {
    let params = default_ani2x_params();
    let aev_len = params.total_aev_length();

    let mut models = HashMap::new();
    for &z in elements {
        if !models.contains_key(&z) {
            models.insert(z, super::weights::make_test_model(aev_len));
        }
    }

    compute_ani(elements, positions, &AniConfig::default(), &models)
}

/// Batch-compute ANI energies for multiple molecules in parallel.
#[cfg(feature = "parallel")]
pub fn compute_ani_batch(
    molecules: &[(&[u8], &[[f64; 3]])],
    config: &AniConfig,
    models: &HashMap<u8, FeedForwardNet>,
) -> Vec<Result<AniResult, String>> {
    use rayon::prelude::*;
    molecules
        .par_iter()
        .map(|(els, pos)| compute_ani(els, pos, config, models))
        .collect()
}

/// Batch-compute ANI energies for multiple molecules sequentially.
#[cfg(not(feature = "parallel"))]
pub fn compute_ani_batch(
    molecules: &[(&[u8], &[[f64; 3]])],
    config: &AniConfig,
    models: &HashMap<u8, FeedForwardNet>,
) -> Vec<Result<AniResult, String>> {
    molecules
        .iter()
        .map(|(els, pos)| compute_ani(els, pos, config, models))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_api_water() {
        let elements = [8u8, 1, 1];
        let positions = [
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469],
        ];
        let result = compute_ani_test(&elements, &positions).unwrap();
        assert_eq!(result.species, vec![8, 1, 1]);
        assert_eq!(result.atomic_energies.len(), 3);
        assert_eq!(result.forces.len(), 3);
        assert!(result.aevs.is_none());
        assert!(result.energy.is_finite());
    }

    #[test]
    fn test_unsupported_element() {
        let elements = [26u8]; // Fe not in ANI
        let positions = [[0.0, 0.0, 0.0]];
        let result = compute_ani_test(&elements, &positions);
        assert!(result.is_err());
    }
}
