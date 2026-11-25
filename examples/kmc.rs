//! Minimal kMC reaction-network example.
//!
//! Run with:
//!
//! ```bash
//! cargo run --example kmc --features alpha-v2
//! ```

use sci_form::kinetics::{ElementaryStep, ThermodynamicState};
use sci_form::kmc::{run_kinetic_monte_carlo, KmcConfig, KmcSpecies, KmcStep};

fn main() {
    let config = KmcConfig {
        species: vec![
            KmcSpecies {
                id: "A".into(),
                initial_count: 1_000,
            },
            KmcSpecies {
                id: "B".into(),
                initial_count: 0,
            },
        ],
        steps: vec![KmcStep {
            step: ElementaryStep {
                step_id: "A->B".into(),
                activation_free_energy_ev: 0.3,
                reaction_free_energy_ev: 0.0,
                prefactor_s_inv: Some(1.0e12),
            },
            reactants: vec![(0, 1)],
            products: vec![(1, 1)],
        }],
        state: ThermodynamicState {
            temperature_k: 500.0,
            pressure_bar: 1.0,
        },
        max_events: 10_000,
        max_time_s: 1.0,
        seed: 7,
        frame_stride: 1_000,
    };
    let result = run_kinetic_monte_carlo(&config).expect("kMC");
    println!(
        "{} events in {:.3e} s: A={} B={}",
        result.events.len(),
        result.final_time_s,
        result.final_counts[0],
        result.final_counts[1]
    );
}
