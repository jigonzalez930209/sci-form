//! GSM Saddle point detection — E11.1c / E11.2c
//!
//! Identifies the transition state from the grown string and
//! refines it using gradient information.

use super::string::{grow_string, GsmConfig, GsmPath};
use serde::{Deserialize, Serialize};

/// Result of a GSM transition state search.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GsmResult {
    /// Transition state coordinates (flat).
    pub ts_coords: Vec<f64>,
    /// Transition state energy (kcal/mol).
    pub ts_energy: f64,
    /// Activation energy: E(TS) - E(reactant) (kcal/mol).
    pub activation_energy: f64,
    /// Reverse barrier: E(TS) - E(product) (kcal/mol).
    pub reverse_barrier: f64,
    /// Energy along the path at each node.
    pub path_energies: Vec<f64>,
    /// Coordinates along the path.
    pub path_coords: Vec<Vec<f64>>,
    /// Index of the TS node in the path.
    pub ts_node_index: usize,
    /// Number of nodes in the path.
    pub n_nodes: usize,
    /// Total energy evaluations.
    pub energy_evaluations: usize,
}

/// Find the highest-energy interior node as the TS candidate.
fn find_ts_candidate(path: &GsmPath) -> (usize, f64) {
    let mut max_idx = 1;
    let mut max_e = f64::NEG_INFINITY;

    for i in 1..(path.nodes.len() - 1) {
        if path.energies[i] > max_e {
            max_e = path.energies[i];
            max_idx = i;
        }
    }

    (max_idx, max_e)
}

/// Refine the saddle point by optimizing the TS candidate.
///
/// Uses a few steps of steepest descent along the perpendicular gradient
/// to minimize the force at the TS.
pub fn refine_saddle(
    ts_coords: &[f64],
    energy_fn: &dyn Fn(&[f64]) -> f64,
    neighbors: (&[f64], &[f64]), // (prev_node, next_node)
    max_steps: usize,
    step_size: f64,
    fd_step: f64,
) -> (Vec<f64>, f64) {
    let n = ts_coords.len();
    let mut coords = ts_coords.to_vec();

    for _ in 0..max_steps {
        // Compute tangent from neighbors
        let mut tangent: Vec<f64> = neighbors
            .1
            .iter()
            .zip(neighbors.0)
            .map(|(a, b)| a - b)
            .collect();
        let norm: f64 = tangent.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm > 1e-15 {
            for t in tangent.iter_mut() {
                *t /= norm;
            }
        }

        // Numerical gradient
        let mut grad = vec![0.0; n];
        for i in 0..n {
            let mut cp = coords.clone();
            let mut cm = coords.clone();
            cp[i] += fd_step;
            cm[i] -= fd_step;
            grad[i] = (energy_fn(&cp) - energy_fn(&cm)) / (2.0 * fd_step);
        }

        // Perpendicular component
        let dot: f64 = grad.iter().zip(&tangent).map(|(g, t)| g * t).sum();
        let perp: Vec<f64> = grad
            .iter()
            .zip(&tangent)
            .map(|(g, t)| g - dot * t)
            .collect();

        let perp_norm: f64 = perp.iter().map(|x| x * x).sum::<f64>().sqrt();
        if perp_norm < 0.01 {
            break;
        }

        // Move along negative perpendicular gradient
        for i in 0..n {
            coords[i] -= step_size * perp[i];
        }
    }

    let energy = energy_fn(&coords);
    (coords, energy)
}

/// Full GSM transition state search.
///
/// 1. Grow string between reactant and product
/// 2. Find highest-energy node (TS candidate)
/// 3. Refine TS using perpendicular gradient minimization
pub fn find_transition_state(
    reactant: &[f64],
    product: &[f64],
    energy_fn: &dyn Fn(&[f64]) -> f64,
    config: &GsmConfig,
) -> GsmResult {
    // Step 1: Grow string
    let path = grow_string(reactant, product, energy_fn, config);

    // Step 2: Find TS candidate
    let (ts_idx, _ts_energy) = find_ts_candidate(&path);

    // Step 3: Refine
    let (prev, next) = if path.nodes.len() >= 3 {
        let prev_idx = if ts_idx > 0 { ts_idx - 1 } else { 0 };
        let next_idx = if ts_idx < path.nodes.len() - 1 {
            ts_idx + 1
        } else {
            path.nodes.len() - 1
        };
        (&path.nodes[prev_idx][..], &path.nodes[next_idx][..])
    } else {
        (reactant, product)
    };

    let (ts_coords, ts_energy) = refine_saddle(
        &path.nodes[ts_idx],
        energy_fn,
        (prev, next),
        20,
        config.step_size * 0.1,
        config.fd_step,
    );

    let e_r = energy_fn(reactant);
    let e_p = energy_fn(product);

    // Estimate energy evaluations:
    // String growth: ~iterations * nodes * 2*n_coords (for gradients)
    // + refinement: ~20 * 2 * n_coords
    let n_coords = reactant.len();
    let energy_evals = path.iterations * path.nodes.len() * 2 * n_coords + 20 * 2 * n_coords;

    GsmResult {
        ts_coords,
        ts_energy,
        activation_energy: ts_energy - e_r,
        reverse_barrier: ts_energy - e_p,
        path_energies: path.energies.clone(),
        path_coords: path.nodes.clone(),
        ts_node_index: ts_idx,
        n_nodes: path.nodes.len(),
        energy_evaluations: energy_evals,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn double_well(coords: &[f64]) -> f64 {
        // 1D double well along x: V = (x² - 1)²
        let x = coords[0];
        (x * x - 1.0).powi(2)
    }

    #[test]
    fn test_find_ts_double_well() {
        let reactant = vec![-1.0, 0.0, 0.0];
        let product = vec![1.0, 0.0, 0.0];
        let config = GsmConfig {
            max_nodes: 9,
            max_iter: 30,
            step_size: 0.01,
            ..Default::default()
        };

        let result = find_transition_state(&reactant, &product, &double_well, &config);

        // TS should be near x=0 where V=1
        assert!(
            result.ts_energy >= 0.0,
            "TS energy should be positive: {}",
            result.ts_energy
        );
        // Activation barrier should be positive
        assert!(
            result.activation_energy > 0.0,
            "Barrier should be positive: {}",
            result.activation_energy
        );
    }

    #[test]
    fn test_gsm_path_endpoints() {
        let reactant = vec![-1.0, 0.0, 0.0];
        let product = vec![1.0, 0.0, 0.0];
        let config = GsmConfig {
            max_nodes: 7,
            max_iter: 10,
            ..Default::default()
        };

        let result = find_transition_state(&reactant, &product, &double_well, &config);

        // First and last nodes should be near reactant/product
        let d_r = (result.path_coords[0][0] - reactant[0]).abs();
        let d_p = (result.path_coords.last().unwrap()[0] - product[0]).abs();
        assert!(d_r < 0.5, "First node should be near reactant: d={}", d_r);
        assert!(d_p < 0.5, "Last node should be near product: d={}", d_p);
    }

    #[test]
    fn test_gsm_result_consistency() {
        let reactant = vec![0.0, 0.0, 0.0];
        let product = vec![2.0, 0.0, 0.0];
        let energy = |c: &[f64]| -> f64 {
            let x = c[0];
            (x - 1.0).powi(2) // single well at x=1
        };

        let config = GsmConfig {
            max_nodes: 5,
            max_iter: 10,
            ..Default::default()
        };

        let result = find_transition_state(&reactant, &product, &energy, &config);
        assert_eq!(result.path_energies.len(), result.n_nodes);
        assert_eq!(result.path_coords.len(), result.n_nodes);
        assert!(result.ts_node_index < result.n_nodes);
    }
}
