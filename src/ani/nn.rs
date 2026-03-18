//! Pure-Rust feed-forward neural network for ANI atomic energy prediction.
//!
//! Minimal inference engine: matrix multiply + bias + activation.
//! Supports GELU and CELU activation functions used in ANI models.

use nalgebra::{DMatrix, DVector};

/// A single dense (fully-connected) layer.
#[derive(Debug, Clone)]
pub struct DenseLayer {
    /// Weight matrix (output_dim × input_dim).
    pub weights: DMatrix<f64>,
    /// Bias vector (output_dim).
    pub bias: DVector<f64>,
    /// Activation function.
    pub activation: Activation,
}

/// Activation function type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Activation {
    /// Gaussian Error Linear Unit.
    Gelu,
    /// Continuously Differentiable Exponential Linear Unit.
    Celu,
    /// No activation (identity, used on the output layer).
    None,
}

/// Feed-forward neural network with multiple dense layers.
#[derive(Debug, Clone)]
pub struct FeedForwardNet {
    pub layers: Vec<DenseLayer>,
}

impl FeedForwardNet {
    /// Create a new network from a list of layers.
    pub fn new(layers: Vec<DenseLayer>) -> Self {
        FeedForwardNet { layers }
    }

    /// Forward pass: compute scalar output from input vector.
    pub fn forward(&self, input: &DVector<f64>) -> f64 {
        let mut x = input.clone();
        for layer in &self.layers {
            x = &layer.weights * &x + &layer.bias;
            apply_activation(&mut x, layer.activation);
        }
        assert_eq!(x.len(), 1, "Output layer must produce a scalar");
        x[0]
    }

    /// Forward pass returning all intermediate activations (for backprop).
    pub fn forward_with_intermediates(
        &self,
        input: &DVector<f64>,
    ) -> Vec<DVector<f64>> {
        let mut activations = Vec::with_capacity(self.layers.len() + 1);
        activations.push(input.clone());

        let mut x = input.clone();
        for layer in &self.layers {
            let z = &layer.weights * &x + &layer.bias;
            let mut a = z.clone();
            apply_activation(&mut a, layer.activation);
            x = a.clone();
            activations.push(a);
        }
        activations
    }

    /// Backward pass: compute gradient of output w.r.t. input.
    pub fn backward(&self, input: &DVector<f64>) -> DVector<f64> {
        // Forward pass storing pre-activation values
        let mut pre_acts = Vec::with_capacity(self.layers.len());
        let mut acts = Vec::with_capacity(self.layers.len() + 1);
        acts.push(input.clone());

        let mut x = input.clone();
        for layer in &self.layers {
            let z = &layer.weights * &x + &layer.bias;
            pre_acts.push(z.clone());
            let mut a = z;
            apply_activation(&mut a, layer.activation);
            x = a.clone();
            acts.push(a);
        }

        // Backward: dL/dz for output layer is 1.0
        let n_layers = self.layers.len();
        let mut grad = DVector::from_element(1, 1.0);

        for l in (0..n_layers).rev() {
            // Apply activation derivative
            let act_deriv = activation_derivative(&pre_acts[l], self.layers[l].activation);
            grad = grad.component_mul(&act_deriv);
            // Propagate through weights
            grad = self.layers[l].weights.transpose() * &grad;
        }
        grad
    }

    /// Input dimension expected by the network.
    pub fn input_dim(&self) -> usize {
        if self.layers.is_empty() {
            0
        } else {
            self.layers[0].weights.ncols()
        }
    }
}

fn apply_activation(x: &mut DVector<f64>, act: Activation) {
    match act {
        Activation::Gelu => {
            for v in x.iter_mut() {
                *v = gelu(*v);
            }
        }
        Activation::Celu => {
            for v in x.iter_mut() {
                *v = celu(*v, 1.0);
            }
        }
        Activation::None => {}
    }
}

fn activation_derivative(z: &DVector<f64>, act: Activation) -> DVector<f64> {
    match act {
        Activation::Gelu => DVector::from_iterator(z.len(), z.iter().map(|&v| gelu_deriv(v))),
        Activation::Celu => DVector::from_iterator(z.len(), z.iter().map(|&v| celu_deriv(v, 1.0))),
        Activation::None => DVector::from_element(z.len(), 1.0),
    }
}

/// GELU activation: x · Φ(x) where Φ is the standard normal CDF.
#[inline]
fn gelu(x: f64) -> f64 {
    0.5 * x * (1.0 + erf(x / std::f64::consts::SQRT_2))
}

/// Approximate GELU derivative.
#[inline]
fn gelu_deriv(x: f64) -> f64 {
    let s2 = std::f64::consts::SQRT_2;
    let phi = 0.5 * (1.0 + erf(x / s2));
    let pdf = (-0.5 * x * x).exp() / (2.0 * std::f64::consts::PI).sqrt();
    phi + x * pdf
}

/// CELU activation: max(0,x) + min(0, α(e^{x/α} - 1)).
#[inline]
fn celu(x: f64, alpha: f64) -> f64 {
    if x >= 0.0 {
        x
    } else {
        alpha * ((x / alpha).exp() - 1.0)
    }
}

/// CELU derivative.
#[inline]
fn celu_deriv(x: f64, alpha: f64) -> f64 {
    if x >= 0.0 {
        1.0
    } else {
        (x / alpha).exp()
    }
}

/// Fast erf approximation (Abramowitz & Stegun 7.1.26, max error 1.5e-7).
fn erf(x: f64) -> f64 {
    let sign = if x >= 0.0 { 1.0 } else { -1.0 };
    let x = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736
                + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    sign * (1.0 - poly * (-x * x).exp())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_net() -> FeedForwardNet {
        let l1 = DenseLayer {
            weights: DMatrix::from_row_slice(3, 2, &[1.0, 0.5, -0.3, 0.8, 0.2, -0.1]),
            bias: DVector::from_vec(vec![0.1, -0.2, 0.05]),
            activation: Activation::Gelu,
        };
        let l2 = DenseLayer {
            weights: DMatrix::from_row_slice(1, 3, &[0.4, -0.6, 0.3]),
            bias: DVector::from_vec(vec![0.0]),
            activation: Activation::None,
        };
        FeedForwardNet::new(vec![l1, l2])
    }

    #[test]
    fn test_forward_deterministic() {
        let net = make_test_net();
        let input = DVector::from_vec(vec![1.0, -0.5]);
        let out1 = net.forward(&input);
        let out2 = net.forward(&input);
        assert!((out1 - out2).abs() < 1e-15);
    }

    #[test]
    fn test_backward_numerical() {
        let net = make_test_net();
        let input = DVector::from_vec(vec![1.0, -0.5]);
        let grad = net.backward(&input);
        let h = 1e-6;
        for d in 0..input.len() {
            let mut inp_p = input.clone();
            let mut inp_m = input.clone();
            inp_p[d] += h;
            inp_m[d] -= h;
            let num = (net.forward(&inp_p) - net.forward(&inp_m)) / (2.0 * h);
            assert!(
                (num - grad[d]).abs() < 1e-4,
                "dim {d}: numerical={num}, analytical={}",
                grad[d]
            );
        }
    }
}
