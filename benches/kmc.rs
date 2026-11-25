//! kMC throughput benchmark for a two-species first-order network.
//!
//! Run with:
//!
//! ```bash
//! cargo bench --bench kmc --features alpha-v2
//! ```

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use sci_form::kinetics::{ElementaryStep, ThermodynamicState};
use sci_form::kmc::{run_kinetic_monte_carlo, KmcConfig, KmcSpecies, KmcStep};

fn kmc_network(c: &mut Criterion) {
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
    let mut group = c.benchmark_group("kmc");
    group.sample_size(20);
    group.bench_function("first_order_10k_events", |b| {
        b.iter(|| black_box(run_kinetic_monte_carlo(black_box(&config)).expect("kMC")))
    });
    group.finish();
}

criterion_group!(benches, kmc_network);
criterion_main!(benches);
