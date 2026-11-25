//! GCMC throughput benchmark on a small ideal-gas box.
//!
//! Run with:
//!
//! ```bash
//! cargo bench --bench gcmc --features alpha-v2
//! ```

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use sci_form::gcmc::{compute_gcmc_isotherm, GcmcConfig, GcmcFramework, GuestModel};

fn gcmc_isotherm(c: &mut Criterion) {
    let framework = GcmcFramework::empty_cubic(25.0);
    let guest = GuestModel::single_site(0.0, 3.0, 0.0);
    let config = GcmcConfig {
        temperature_k: 300.0,
        pressures_bar: vec![1.0],
        equilibration_steps: 100,
        production_steps: 1_000,
        sample_interval: 50,
        seed: 7,
        ..GcmcConfig::default()
    };
    let mut group = c.benchmark_group("gcmc");
    group.sample_size(20);
    group.bench_function("ideal_gas_1000_moves", |b| {
        b.iter(|| {
            black_box(
                compute_gcmc_isotherm(black_box(&framework), black_box(&guest), black_box(&config))
                    .expect("GCMC"),
            )
        })
    });
    group.finish();
}

criterion_group!(benches, gcmc_isotherm);
criterion_main!(benches);
