//! Enhanced-sampling benchmarks: OPES flow and replica exchange.
//!
//! Run with:
//!
//! ```bash
//! cargo bench --bench enhanced_sampling
//! ```

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use sci_form::enhanced_sampling::{
    run_opes, run_replica_exchange, CollectiveVariable, DoubleWell, OpesConfig, Potential,
    ReplicaExchangeConfig,
};

struct Harmonic;

impl Potential for Harmonic {
    fn energy(&self, coords: &[f64]) -> f64 {
        0.5 * coords[0] * coords[0]
    }

    fn gradient(&self, coords: &[f64], gradient: &mut [f64]) {
        gradient[0] = coords[0];
    }
}

fn opes(c: &mut Criterion) {
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
    let mut group = c.benchmark_group("enhanced_sampling");
    group.sample_size(10);
    group.bench_function("opes_1d_20k_steps", |b| {
        b.iter(|| {
            black_box(
                run_opes(black_box(&Harmonic), &[0.0], &[10.0], black_box(&config)).expect("OPES"),
            )
        })
    });
    group.finish();
}

fn replica_exchange(c: &mut Criterion) {
    let potential = DoubleWell::new(5.0, 1.0);
    let config = ReplicaExchangeConfig {
        temperatures_k: vec![300.0, 600.0, 1300.0],
        max_steps: 10_000,
        swap_interval: 100,
        trace_stride: 1_000,
        ..ReplicaExchangeConfig::default()
    };
    let mut group = c.benchmark_group("enhanced_sampling");
    group.sample_size(10);
    group.bench_function("replica_exchange_3x10k", |b| {
        b.iter(|| {
            black_box(
                run_replica_exchange(black_box(&potential), &[1.0], &[10.0], black_box(&config))
                    .expect("replica exchange"),
            )
        })
    });
    group.finish();
}

criterion_group!(benches, opes, replica_exchange);
criterion_main!(benches);
