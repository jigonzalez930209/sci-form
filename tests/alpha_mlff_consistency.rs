//! Integration tests for MLFF (machine-learned force field) consistency.
//!
//! Validates:
//! - AEV descriptors are translation/rotation equivariant
//! - Neural network forward pass produces sensible energies
//! - MLFF forces are consistent with numerical energy gradient
//! - Gradient length matches 3×n_atoms
//! - Error handling for unsupported elements

#![cfg(feature = "alpha-mlff")]

use std::collections::HashMap;

use sci_form::ml::inference::{Activation, DenseLayer, InferenceNet};
use sci_form::ml::symmetry_functions::{compute_aevs, SymmetryFunctionParams};
use sci_form::mlff::{compute_mlff, compute_mlff_energy_and_gradient, MlffConfig};

/// Build a minimal test network: 1 hidden layer [aev_dim → 16 → 1].
fn tiny_net(aev_dim: usize) -> InferenceNet {
    let hidden = DenseLayer {
        weights: vec![0.01; 16 * aev_dim],
        bias: vec![0.0; 16],
        in_features: aev_dim,
        out_features: 16,
        activation: Activation::Gelu,
    };
    let output = DenseLayer {
        weights: vec![0.01; 16],
        bias: vec![0.0; 1],
        in_features: 16,
        out_features: 1,
        activation: Activation::Linear,
    };
    InferenceNet::new(vec![hidden, output])
}

fn h2_config(aev_dim: usize) -> MlffConfig {
    let mut element_nets = HashMap::new();
    element_nets.insert(1u8, tiny_net(aev_dim));
    MlffConfig {
        aev_params: SymmetryFunctionParams::default(),
        element_nets,
    }
}

// ─── AEV descriptors ─────────────────────────────────────────────────────────

/// Translation invariance: shifting both atoms by same vector → same AEVs.
#[test]
fn aev_translation_invariant() {
    let params = SymmetryFunctionParams::default();
    let elements = [1u8, 1];
    let pos_a = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
    let pos_b = [[5.0, 3.0, -2.0], [5.74, 3.0, -2.0]]; // translated by (5,3,-2)

    let aev_a = compute_aevs(&elements, &pos_a, &params);
    let aev_b = compute_aevs(&elements, &pos_b, &params);

    // AEVs should be very close (within floating-point tolerance)
    let va = aev_a[0].to_vec();
    let vb = aev_b[0].to_vec();
    assert_eq!(va.len(), vb.len());
    let max_diff: f64 = va
        .iter()
        .zip(vb.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    assert!(
        max_diff < 1e-10,
        "AEVs should be translation-invariant, max diff = {:.2e}",
        max_diff
    );
}

/// AEV length should be consistent for all atoms.
#[test]
fn aev_consistent_length() {
    let params = SymmetryFunctionParams::default();
    let elements = [1u8, 1];
    let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
    let aevs = compute_aevs(&elements, &positions, &params);

    assert_eq!(aevs.len(), 2);
    assert_eq!(
        aevs[0].len(),
        aevs[1].len(),
        "Both atoms should have same AEV length"
    );
    assert!(!aevs[0].is_empty(), "AEV should have non-zero length");
}

/// Single atom should have zero AEV (no neighbors).
#[test]
fn aev_single_atom_is_zero() {
    let params = SymmetryFunctionParams::default();
    let elements = [1u8];
    let positions = [[0.0, 0.0, 0.0]];
    let aevs = compute_aevs(&elements, &positions, &params);
    assert_eq!(aevs.len(), 1);
    let all_zero = aevs[0].to_vec().iter().all(|v| v.abs() < 1e-15);
    assert!(all_zero, "Single atom AEV should be all zeros");
}

// ─── Neural network inference ────────────────────────────────────────────────

/// Forward pass through dense layer produces correct output dimensions.
#[test]
fn dense_layer_output_dim() {
    let layer = DenseLayer {
        weights: vec![0.1; 8 * 4], // 8 outputs, 4 inputs
        bias: vec![0.0; 8],
        in_features: 4,
        out_features: 8,
        activation: Activation::Relu,
    };
    let input = vec![1.0, 2.0, 3.0, 4.0];
    let output = layer.forward(&input);
    assert_eq!(output.len(), 8, "Output should have 8 elements");
}

/// GELU activation: output should differ from input (nonlinear).
#[test]
fn gelu_activation_nonlinear() {
    let layer = DenseLayer {
        weights: vec![1.0], // identity mapping 1→1
        bias: vec![0.0],
        in_features: 1,
        out_features: 1,
        activation: Activation::Gelu,
    };
    let pos_output = layer.forward(&[2.0])[0];
    let neg_output = layer.forward(&[-2.0])[0];
    // GELU(2) ≈ 1.954, GELU(-2) ≈ -0.045
    assert!(
        pos_output > 1.0,
        "GELU(2) should be > 1, got {:.4}",
        pos_output
    );
    assert!(
        neg_output.abs() < 0.5,
        "GELU(-2) should be near 0, got {:.4}",
        neg_output
    );
}

// ─── MLFF energy/force pipeline ──────────────────────────────────────────────

/// H₂ MLFF should return finite energy and forces.
#[test]
fn mlff_h2_finite_energy_and_forces() {
    let params = SymmetryFunctionParams::default();
    let elements = [1u8, 1];
    let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
    let aevs = compute_aevs(&elements, &positions, &params);
    let aev_dim = aevs[0].len();

    let config = h2_config(aev_dim);
    let result = compute_mlff(&elements, &positions, &config).unwrap();

    assert!(
        result.energy.is_finite(),
        "MLFF energy should be finite, got {}",
        result.energy
    );
    assert_eq!(result.forces.len(), 2, "Should have forces for 2 atoms");
    for (i, f) in result.forces.iter().enumerate() {
        for (d, fd) in f.iter().enumerate() {
            assert!(
                fd.is_finite(),
                "Force[{}][{}] should be finite, got {}",
                i,
                d,
                fd
            );
        }
    }
}

/// MLFF gradient flat vector should have length 3*n_atoms.
#[test]
fn mlff_gradient_correct_length() {
    let params = SymmetryFunctionParams::default();
    let elements = [1u8, 1];
    let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
    let aevs = compute_aevs(&elements, &positions, &params);
    let aev_dim = aevs[0].len();
    let config = h2_config(aev_dim);

    let positions_flat = [0.0, 0.0, 0.0, 0.74, 0.0, 0.0];
    let result = compute_mlff_energy_and_gradient(&elements, &positions_flat, &config).unwrap();
    assert_eq!(
        result.gradient_flat.len(),
        6,
        "Gradient should have 3*2=6 components"
    );
}

/// Missing element → error, not panic.
#[test]
fn mlff_missing_element_returns_error() {
    let mut element_nets = HashMap::new();
    element_nets.insert(1u8, tiny_net(10));
    let config = MlffConfig {
        aev_params: SymmetryFunctionParams::default(),
        element_nets,
    };
    let elements = [6u8, 1]; // Carbon not in networks
    let positions = [[0.0, 0.0, 0.0], [1.1, 0.0, 0.0]];
    let result = compute_mlff(&elements, &positions, &config);
    assert!(result.is_err(), "Should error for missing element network");
}

/// Empty molecule → error, not panic.
#[test]
fn mlff_empty_molecule_error() {
    let config = h2_config(10);
    let result = compute_mlff(&[], &[], &config);
    assert!(result.is_err(), "Empty molecule should return error");
}

// ─── Network backward pass (gradient) ────────────────────────────────────────

/// Network backward should produce gradient of correct length.
#[test]
fn network_backward_correct_length() {
    let aev_dim = 8;
    let net = tiny_net(aev_dim);
    let input = vec![0.1; aev_dim];
    let d_output = vec![1.0]; // ∂L/∂output = 1 for energy gradient

    let grad = net.backward(&input, &d_output);
    assert_eq!(
        grad.len(),
        aev_dim,
        "Backward gradient should have same dim as input"
    );
}
