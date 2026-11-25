//! Active-learning acquisition benchmark.
//!
//! Run with:
//!
//! ```bash
//! cargo bench --bench active_learning --features alpha-v2
//! ```

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use sci_form::active_learning::{propose_active_learning, Acquisition, LabeledPoint};
use sci_form::ml::advanced_models::RandomForest;

fn acquisition(c: &mut Criterion) {
    let model = RandomForest {
        trees: Vec::new(),
        n_trees: 0,
        oob_score: None,
    };
    let labeled: Vec<LabeledPoint> = (0..20)
        .map(|index| LabeledPoint {
            features: vec![index as f64 / 20.0],
            target: 0.0,
        })
        .collect();
    let candidates: Vec<Vec<f64>> = (0..1_000).map(|i| vec![i as f64 / 1_000.0]).collect();
    let mut group = c.benchmark_group("active_learning");
    group.sample_size(20);
    group.bench_function("distance_1k_candidates_10_labels", |b| {
        b.iter(|| {
            black_box(
                propose_active_learning(
                    black_box(&model),
                    black_box(&labeled),
                    black_box(&candidates),
                    black_box(&Acquisition::Distance),
                    10,
                )
                .expect("active learning"),
            )
        })
    });
    group.finish();
}

criterion_group!(benches, acquisition);
criterion_main!(benches);
