//! Minimal OPES example on a one-dimensional harmonic surface.
//!
//! Run with:
//!
//! ```bash
//! cargo run --example opes
//! ```

use sci_form::enhanced_sampling::{run_opes, CollectiveVariable, OpesConfig, Potential};

struct Harmonic;

impl Potential for Harmonic {
    fn energy(&self, coords: &[f64]) -> f64 {
        0.5 * coords[0] * coords[0]
    }

    fn gradient(&self, coords: &[f64], gradient: &mut [f64]) {
        gradient[0] = coords[0];
    }
}

fn main() {
    let config = OpesConfig {
        cv: CollectiveVariable::Coordinate { index: 0 },
        sigma: 0.2,
        bias_factor: 10.0,
        pace_steps: 200,
        temperature_k: 300.0,
        max_steps: 20_000,
        grid_min: -1.5,
        grid_max: 1.5,
        grid_bins: 61,
        trace_stride: 1_000,
        flow_iterations: 50,
        ..OpesConfig::default()
    };
    let result = run_opes(&Harmonic, &[0.0], &[10.0], &config).expect("OPES");
    println!(
        "kernels={} flow_residual={:.3e} converged={}",
        result.kernels.len(),
        result.flow_residual,
        result.flow_converged
    );
    for (index, value) in result.free_energy.values.iter().enumerate() {
        let x = result.free_energy.coordinate(index);
        println!("{x:6.2} {value:8.3}");
    }
}
