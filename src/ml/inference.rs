//! **ALPHA** — Lightweight feed-forward neural network inference.
//!
//! Provides a minimal `FeedForwardNet` that performs forward passes using
//! pre-loaded weights. Designed for ML force field potentials that store
//! element-specific networks (e.g., ANI-style).

use serde::{Deserialize, Serialize};

/// Activation function applied after each hidden layer.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum Activation {
    /// Gaussian Error Linear Unit  (GELU).
    Gelu,
    /// Rectified Linear Unit.
    Relu,
    /// Continuously Differentiable Exponential Linear Unit.
    Celu,
    /// No activation (identity).
    Linear,
}

impl Activation {
    fn apply(&self, x: f64) -> f64 {
        match self {
            Activation::Gelu => x * 0.5 * (1.0 + erf_approx(x / std::f64::consts::SQRT_2)),
            Activation::Relu => x.max(0.0),
            Activation::Celu => {
                let alpha = 1.0;
                if x >= 0.0 {
                    x
                } else {
                    alpha * ((x / alpha).exp() - 1.0)
                }
            }
            Activation::Linear => x,
        }
    }

    fn derivative(&self, x: f64) -> f64 {
        match self {
            Activation::Gelu => {
                let cdf = 0.5 * (1.0 + erf_approx(x / std::f64::consts::SQRT_2));
                let pdf = (-0.5 * x * x).exp() / (2.0 * std::f64::consts::PI).sqrt();
                cdf + x * pdf
            }
            Activation::Relu => {
                if x > 0.0 {
                    1.0
                } else {
                    0.0
                }
            }
            Activation::Celu => {
                let alpha = 1.0;
                if x >= 0.0 {
                    1.0
                } else {
                    (x / alpha).exp()
                }
            }
            Activation::Linear => 1.0,
        }
    }
}

/// Dense layer: output = activation(W · input + bias).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DenseLayer {
    /// Weight matrix, row-major [out_features × in_features].
    pub weights: Vec<f64>,
    /// Bias vector [out_features].
    pub bias: Vec<f64>,
    /// Input dimension.
    pub in_features: usize,
    /// Output dimension.
    pub out_features: usize,
    /// Activation function.
    pub activation: Activation,
}

impl DenseLayer {
    /// Forward pass for a single input vector.
    pub fn forward(&self, input: &[f64]) -> Vec<f64> {
        debug_assert_eq!(input.len(), self.in_features);
        let mut output = vec![0.0; self.out_features];
        for o in 0..self.out_features {
            let mut sum = self.bias[o];
            let row_start = o * self.in_features;
            for i in 0..self.in_features {
                sum += self.weights[row_start + i] * input[i];
            }
            output[o] = self.activation.apply(sum);
        }
        output
    }

    /// Forward pass storing pre-activation values (for gradient computation).
    pub fn forward_with_cache(&self, input: &[f64]) -> (Vec<f64>, Vec<f64>) {
        debug_assert_eq!(input.len(), self.in_features);
        let mut pre_act = vec![0.0; self.out_features];
        let mut output = vec![0.0; self.out_features];
        for o in 0..self.out_features {
            let mut sum = self.bias[o];
            let row_start = o * self.in_features;
            for i in 0..self.in_features {
                sum += self.weights[row_start + i] * input[i];
            }
            pre_act[o] = sum;
            output[o] = self.activation.apply(sum);
        }
        (output, pre_act)
    }
}

/// A simple feed-forward network composed of sequential dense layers.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InferenceNet {
    pub layers: Vec<DenseLayer>,
}

impl InferenceNet {
    /// Create from a list of layers.
    pub fn new(layers: Vec<DenseLayer>) -> Self {
        Self { layers }
    }

    /// Forward pass: input → final output.
    pub fn forward(&self, input: &[f64]) -> Vec<f64> {
        let mut x = input.to_vec();
        for layer in &self.layers {
            x = layer.forward(&x);
        }
        x
    }

    /// Forward pass that stores all intermediate activations (for backprop).
    pub fn forward_with_intermediates(
        &self,
        input: &[f64],
    ) -> (Vec<f64>, Vec<Vec<f64>>, Vec<Vec<f64>>) {
        let mut activations = Vec::with_capacity(self.layers.len() + 1);
        let mut pre_activations = Vec::with_capacity(self.layers.len());
        activations.push(input.to_vec());

        let mut x = input.to_vec();
        for layer in &self.layers {
            let (out, pre) = layer.forward_with_cache(&x);
            pre_activations.push(pre);
            activations.push(out.clone());
            x = out;
        }
        (x, activations, pre_activations)
    }

    /// Compute dE/d(input) via backpropagation given dE/d(output).
    /// Returns gradient w.r.t. the network input.
    pub fn backward(&self, input: &[f64], d_output: &[f64]) -> Vec<f64> {
        let (_output, activations, pre_activations) = self.forward_with_intermediates(input);

        let mut d_next = d_output.to_vec();

        for l in (0..self.layers.len()).rev() {
            let layer = &self.layers[l];
            let pre_act = &pre_activations[l];

            // d_pre = d_next ⊙ activation'(pre_act)
            let d_pre: Vec<f64> = d_next
                .iter()
                .zip(pre_act.iter())
                .map(|(&dn, &pa)| dn * layer.activation.derivative(pa))
                .collect();

            // d_input = Wᵀ · d_pre
            let input_act = &activations[l];
            let mut d_input = vec![0.0; input_act.len()];
            for o in 0..layer.out_features {
                let row_start = o * layer.in_features;
                for i in 0..layer.in_features {
                    d_input[i] += layer.weights[row_start + i] * d_pre[o];
                }
            }
            d_next = d_input;
        }

        d_next
    }
}

/// Approximate erf using Abramowitz & Stegun formula 7.1.26.
fn erf_approx(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    sign * (1.0 - poly * (-x * x).exp())
}

/// Result from MLFF inference.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MlffResult {
    /// Total predicted energy.
    pub energy: f64,
    /// Per-atom energies.
    pub atomic_energies: Vec<f64>,
    /// Forces on each atom [n_atoms × 3]. Computed via chain rule through AEVs.
    pub forces: Vec<[f64; 3]>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dense_layer_forward() {
        let layer = DenseLayer {
            weights: vec![1.0, 0.0, 0.0, 1.0],
            bias: vec![0.1, -0.1],
            in_features: 2,
            out_features: 2,
            activation: Activation::Linear,
        };
        let out = layer.forward(&[2.0, 3.0]);
        assert!((out[0] - 2.1).abs() < 1e-10);
        assert!((out[1] - 2.9).abs() < 1e-10);
    }

    #[test]
    fn test_inference_net_roundtrip() {
        let net = InferenceNet::new(vec![
            DenseLayer {
                weights: vec![1.0, 2.0, 3.0, 4.0],
                bias: vec![0.0, 0.0],
                in_features: 2,
                out_features: 2,
                activation: Activation::Relu,
            },
            DenseLayer {
                weights: vec![1.0, 1.0],
                bias: vec![0.0],
                in_features: 2,
                out_features: 1,
                activation: Activation::Linear,
            },
        ]);
        let out = net.forward(&[1.0, 1.0]);
        // Layer 1: [1*1+2*1, 3*1+4*1] = [3, 7], relu → [3, 7]
        // Layer 2: [1*3+1*7] = [10]
        assert!((out[0] - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_erf_approx() {
        assert!((erf_approx(0.0)).abs() < 1e-6);
        assert!((erf_approx(1.0) - 0.8427).abs() < 0.001);
        assert!((erf_approx(-1.0) + 0.8427).abs() < 0.001);
    }
}
