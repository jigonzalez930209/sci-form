//! GSM String growth — E11.1
//!
//! Core algorithm for growing the string of images between reactant and product.
//!
//! # Stability tier
//!
//! **Alpha** — early research tier. Gated behind the `alpha-gsm` feature
//! flag; the validation gate is `tests/alpha_gsm_validation.rs` (F82).
//!
//! # Units
//!
//! Coordinates are flat `3N` arrays in **Å**, energies in **kcal/mol** and
//! gradients in **kcal/(mol·Å)**.

use serde::{Deserialize, Serialize};

/// GSM configuration.
///
/// Controls the string growth: how many images the string may hold, how far
/// a new image is inserted, and when the growth is considered converged.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct GsmConfig {
    /// Maximum number of nodes in the string (including endpoints).
    pub max_nodes: usize,
    /// Step size for growing new nodes (Å).
    pub step_size: f64,
    /// Convergence tolerance for perpendicular gradient (kcal/(mol·Å)).
    pub grad_tol: f64,
    /// Maximum growth iterations.
    pub max_iter: usize,
    /// Finite difference step for numerical gradient (Å).
    pub fd_step: f64,
}

impl Default for GsmConfig {
    fn default() -> Self {
        Self {
            max_nodes: 11,
            step_size: 0.1,
            grad_tol: 0.5,
            max_iter: 100,
            fd_step: 0.005,
        }
    }
}

/// A string of images along the reaction path.
#[derive(Debug, Clone)]
pub struct GsmPath {
    /// Coordinates of each node (each node is a flat Vec of 3N coordinates in
    /// **Å**).
    pub nodes: Vec<Vec<f64>>,
    /// Energy at each node (kcal/mol).
    pub energies: Vec<f64>,
    /// Perpendicular gradient magnitude at each node, in **kcal/(mol·Å)**.
    pub perp_grad_norms: Vec<f64>,
    /// Whether the string has joined (tips are close enough).
    pub joined: bool,
    /// Number of growth iterations performed.
    pub iterations: usize,
}

/// Linearly interpolate between two coordinate sets.
///
/// Returns: reactant + t * (product - reactant)
pub fn interpolate_node(reactant: &[f64], product: &[f64], t: f64) -> Vec<f64> {
    reactant
        .iter()
        .zip(product.iter())
        .map(|(r, p)| r + t * (p - r))
        .collect()
}

/// Compute numerical gradient of energy function.
fn numerical_gradient(coords: &[f64], energy_fn: &dyn Fn(&[f64]) -> f64, h: f64) -> Vec<f64> {
    let n = coords.len();
    let mut grad = vec![0.0; n];

    for i in 0..n {
        let mut cp = coords.to_vec();
        let mut cm = coords.to_vec();
        cp[i] += h;
        cm[i] -= h;
        grad[i] = (energy_fn(&cp) - energy_fn(&cm)) / (2.0 * h);
    }

    grad
}

/// Arc-length reparametrization (core GSM step).
///
/// Keeps the endpoints fixed and redistributes interior nodes uniformly
/// along the polygonal path. Without this, interior nodes slide downhill
/// along the tangent into the minima and the string maximum collapses
/// to an endpoint instead of converging to the MEP/TS.
fn reparametrize_nodes(nodes: &mut Vec<Vec<f64>>) {
    let m = nodes.len();
    if m <= 2 {
        return;
    }
    let dim = nodes[0].len();
    let seg_len = |a: &[f64], b: &[f64]| {
        a.iter()
            .zip(b)
            .map(|(x, y)| (x - y) * (x - y))
            .sum::<f64>()
            .sqrt()
    };
    let mut cum = vec![0.0; m];
    for i in 1..m {
        cum[i] = cum[i - 1] + seg_len(&nodes[i - 1], &nodes[i]);
    }
    let total = cum[m - 1];
    if total < 1e-15 {
        return;
    }
    let mut out = Vec::with_capacity(m);
    out.push(nodes[0].clone());
    let mut seg = 0;
    for i in 1..(m - 1) {
        let target = total * i as f64 / (m - 1) as f64;
        while seg + 1 < m - 1 && cum[seg + 1] < target {
            seg += 1;
        }
        let lo = cum[seg];
        let hi = cum[seg + 1];
        let t = if hi > lo {
            ((target - lo) / (hi - lo)).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let mut p = vec![0.0; dim];
        for k in 0..dim {
            p[k] = nodes[seg][k] + t * (nodes[seg + 1][k] - nodes[seg][k]);
        }
        out.push(p);
    }
    out.push(nodes[m - 1].clone());
    *nodes = out;
}
/// Compute string tangent at a given node (central difference).
fn compute_tangent(nodes: &[Vec<f64>], idx: usize) -> Vec<f64> {
    if idx == 0 {
        // Forward difference
        let mut t: Vec<f64> = nodes[1].iter().zip(&nodes[0]).map(|(a, b)| a - b).collect();
        normalize(&mut t);
        t
    } else if idx == nodes.len() - 1 {
        // Backward difference
        let mut t: Vec<f64> = nodes[idx]
            .iter()
            .zip(&nodes[idx - 1])
            .map(|(a, b)| a - b)
            .collect();
        normalize(&mut t);
        t
    } else {
        // Central difference
        let mut t: Vec<f64> = nodes[idx + 1]
            .iter()
            .zip(&nodes[idx - 1])
            .map(|(a, b)| a - b)
            .collect();
        normalize(&mut t);
        t
    }
}

/// Compute perpendicular gradient: g_perp = g - (g·τ)τ
fn perpendicular_gradient(gradient: &[f64], tangent: &[f64]) -> Vec<f64> {
    let dot: f64 = gradient.iter().zip(tangent).map(|(g, t)| g * t).sum();
    gradient
        .iter()
        .zip(tangent)
        .map(|(g, t)| g - dot * t)
        .collect()
}

fn normalize(v: &mut [f64]) {
    let norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    if norm > 1e-15 {
        for x in v.iter_mut() {
            *x /= norm;
        }
    }
}

fn vec_norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn distance(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// Grow a string of images between reactant and product.
///
/// The string grows from both ends toward the middle, adding new nodes
/// along the negative perpendicular gradient direction.
pub fn grow_string(
    reactant: &[f64],
    product: &[f64],
    energy_fn: &dyn Fn(&[f64]) -> f64,
    config: &GsmConfig,
) -> GsmPath {
    let n_coords = reactant.len();

    // Initialize with 3 nodes: R, midpoint, P
    let mid = interpolate_node(reactant, product, 0.5);
    let mut nodes = vec![reactant.to_vec(), mid, product.to_vec()];

    let mut joined = false;
    let mut iterations = 0;

    for it in 0..config.max_iter {
        iterations = it + 1;

        if nodes.len() >= config.max_nodes {
            break;
        }

        // Try to add nodes from the reactant side
        if nodes.len() < config.max_nodes {
            // Interpolate a new node between first interior node and its neighbor
            let idx = 1; // after reactant
            let new = interpolate_node(&nodes[0], &nodes[idx], 0.5);
            nodes.insert(1, new);
        }

        // Try to add from the product side
        if nodes.len() < config.max_nodes {
            let idx = nodes.len() - 2; // before product
            let new = interpolate_node(&nodes[idx], &nodes[nodes.len() - 1], 0.5);
            nodes.insert(nodes.len() - 1, new);
        }

        // Optimize interior nodes: move along perpendicular gradient
        for i in 1..(nodes.len() - 1) {
            let tangent = compute_tangent(&nodes, i);
            let grad = numerical_gradient(&nodes[i], energy_fn, config.fd_step);
            let perp = perpendicular_gradient(&grad, &tangent);

            // Step along negative perpendicular gradient
            let step = config.step_size;
            for k in 0..n_coords {
                nodes[i][k] -= step * perp[k];
            }
        }

        // Arc-length reparametrization keeps interior nodes on the MEP.
        reparametrize_nodes(&mut nodes);

        // Check junction: distance between adjacent interior nodes
        let mut all_close = true;
        for i in 0..(nodes.len() - 1) {
            let d = distance(&nodes[i], &nodes[i + 1]);
            if d > config.step_size * 2.0 {
                all_close = false;
                break;
            }
        }

        if all_close && nodes.len() >= 5 {
            joined = true;
            break;
        }

        // Check convergence: all perpendicular gradients small
        let mut max_perp = 0.0;
        for i in 1..(nodes.len() - 1) {
            let tangent = compute_tangent(&nodes, i);
            let grad = numerical_gradient(&nodes[i], energy_fn, config.fd_step);
            let perp = perpendicular_gradient(&grad, &tangent);
            let norm = vec_norm(&perp);
            if norm > max_perp {
                max_perp = norm;
            }
        }

        if max_perp < config.grad_tol {
            joined = true;
            break;
        }
    }

    // Compute final energies and gradients
    let energies: Vec<f64> = nodes.iter().map(|n| energy_fn(n)).collect();

    let mut perp_norms = Vec::new();
    for i in 0..nodes.len() {
        if i == 0 || i == nodes.len() - 1 {
            perp_norms.push(0.0);
        } else {
            let tangent = compute_tangent(&nodes, i);
            let grad = numerical_gradient(&nodes[i], energy_fn, config.fd_step);
            let perp = perpendicular_gradient(&grad, &tangent);
            perp_norms.push(vec_norm(&perp));
        }
    }

    GsmPath {
        nodes,
        energies,
        perp_grad_norms: perp_norms,
        joined,
        iterations,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interpolate_midpoint() {
        let r = vec![0.0, 0.0, 0.0];
        let p = vec![2.0, 2.0, 2.0];
        let mid = interpolate_node(&r, &p, 0.5);
        assert_eq!(mid, vec![1.0, 1.0, 1.0]);
    }

    #[test]
    fn test_interpolate_endpoints() {
        let r = vec![1.0, 2.0, 3.0];
        let p = vec![4.0, 5.0, 6.0];
        let r0 = interpolate_node(&r, &p, 0.0);
        let p1 = interpolate_node(&r, &p, 1.0);
        assert_eq!(r0, r);
        assert_eq!(p1, p);
    }

    #[test]
    fn test_perpendicular_gradient() {
        let g = vec![1.0, 1.0, 0.0];
        let t = vec![1.0, 0.0, 0.0]; // tangent along x
        let perp = perpendicular_gradient(&g, &t);
        // Should remove x-component
        assert!((perp[0]).abs() < 1e-10);
        assert!((perp[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_grow_string_basic() {
        // Simple 1D double-well potential
        let energy_fn = |coords: &[f64]| -> f64 {
            let x = coords[0];
            (x * x - 1.0).powi(2)
        };

        let reactant = vec![-1.0, 0.0, 0.0];
        let product = vec![1.0, 0.0, 0.0];
        let config = GsmConfig {
            max_nodes: 7,
            max_iter: 20,
            step_size: 0.01,
            ..Default::default()
        };

        let path = grow_string(&reactant, &product, &energy_fn, &config);
        assert!(path.nodes.len() >= 3);
        assert_eq!(path.energies.len(), path.nodes.len());
    }
}
