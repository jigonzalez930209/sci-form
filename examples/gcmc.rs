//! Minimal GCMC adsorption example.
//!
//! Run with:
//!
//! ```bash
//! cargo run --example gcmc --features alpha-v2
//! ```

use sci_form::gcmc::{compute_gcmc_isotherm, GcmcConfig, GcmcFramework, GuestModel};

fn main() {
    let framework = GcmcFramework::empty_cubic(30.0);
    let guest = GuestModel::single_site(0.0, 3.0, 0.0);
    let config = GcmcConfig {
        temperature_k: 300.0,
        pressures_bar: vec![1.0, 10.0],
        equilibration_steps: 2_000,
        production_steps: 10_000,
        seed: 42,
        ..GcmcConfig::default()
    };
    let isotherm = compute_gcmc_isotherm(&framework, &guest, &config).expect("GCMC");
    for point in &isotherm.points {
        println!(
            "{:6.2} bar -> {:8.3} +/- {:.3} molecules/cell ({} samples)",
            point.pressure_bar, point.mean_loading, point.standard_error, point.n_samples
        );
    }
}
