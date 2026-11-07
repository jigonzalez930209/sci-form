//! Core benchmarks: parsing, conformer generation, force fields and the
//! semi-empirical/electronic-structure entry points.
//!
//! Run with:
//!
//! ```bash
//! cargo bench --bench core
//! cargo bench --bench core -- embed          # one group
//! ```
//!
//! The benchmark inputs are the same small molecule set used by the quick
//! validation gates, so a regression here is directly comparable with the
//! numbers in `reports/validation/summary.json`.

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use std::time::Duration;

const PARSE_SMILES: &[(&str, &str)] = &[
    ("ethanol", "CCO"),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("caffeine", "Cn1c(=O)c2c(ncn2C)n(C)c1=O"),
    (
        "cholesterol",
        "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C",
    ),
];

const EMBED_SMILES: &[(&str, &str)] = &[
    ("ethanol", "CCO"),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    (
        "cholesterol",
        "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C",
    ),
];

fn parse(c: &mut Criterion) {
    let mut group = c.benchmark_group("parse");
    for (name, smiles) in PARSE_SMILES {
        group.bench_function(*name, |b| {
            b.iter(|| black_box(sci_form::parse(black_box(smiles)).expect("valid SMILES")))
        });
    }
    group.finish();
}

fn embed(c: &mut Criterion) {
    let mut group = c.benchmark_group("embed");
    for (name, smiles) in EMBED_SMILES {
        group.bench_function(*name, |b| {
            b.iter(|| black_box(sci_form::embed(black_box(smiles), 42)))
        });
    }
    group.finish();
}

fn force_fields(c: &mut Criterion) {
    let smiles = "CC(=O)Oc1ccccc1C(=O)O";
    let conformer = sci_form::embed(smiles, 42);
    assert!(
        conformer.error.is_none(),
        "embed failed: {:?}",
        conformer.error
    );
    let coords = conformer.coords;

    let mut group = c.benchmark_group("force_field");
    group.bench_function("uff_energy_aspirin", |b| {
        b.iter(|| {
            black_box(
                sci_form::compute_uff_energy(black_box(smiles), black_box(&coords))
                    .expect("UFF energy"),
            )
        })
    });
    group.bench_function("mmff94_energy_aspirin", |b| {
        b.iter(|| {
            black_box(
                sci_form::compute_mmff94_energy(black_box(smiles), black_box(&coords))
                    .expect("MMFF94 energy"),
            )
        })
    });
    group.finish();
}

fn elements_and_positions(smiles: &str) -> (Vec<u8>, Vec<[f64; 3]>) {
    let conformer = sci_form::embed(smiles, 42);
    assert!(
        conformer.error.is_none(),
        "embed failed for {smiles}: {:?}",
        conformer.error
    );
    let positions = conformer
        .coords
        .chunks_exact(3)
        .map(|p| [p[0], p[1], p[2]])
        .collect();
    (conformer.elements, positions)
}

fn quantum(c: &mut Criterion) {
    let (water_z, water_xyz) = elements_and_positions("O");
    let (methane_z, methane_xyz) = elements_and_positions("C");
    let (benzene_z, benzene_xyz) = elements_and_positions("c1ccccc1");

    let mut group = c.benchmark_group("quantum");
    group.sample_size(30);

    group.bench_function("eht_water", |b| {
        b.iter(|| {
            black_box(
                sci_form::eht::solve_eht(black_box(&water_z), black_box(&water_xyz), None)
                    .expect("EHT water"),
            )
        })
    });
    group.bench_function("eht_benzene", |b| {
        b.iter(|| {
            black_box(
                sci_form::eht::solve_eht(black_box(&benzene_z), black_box(&benzene_xyz), None)
                    .expect("EHT benzene"),
            )
        })
    });
    group.bench_function("pm3_water", |b| {
        b.iter(|| {
            black_box(
                sci_form::compute_pm3(black_box(&water_z), black_box(&water_xyz))
                    .expect("PM3 water"),
            )
        })
    });
    group.bench_function("pm3_methane", |b| {
        b.iter(|| {
            black_box(
                sci_form::compute_pm3(black_box(&methane_z), black_box(&methane_xyz))
                    .expect("PM3 methane"),
            )
        })
    });
    group.bench_function("gfn2_water", |b| {
        b.iter(|| {
            black_box(
                sci_form::compute_xtb(black_box(&water_z), black_box(&water_xyz))
                    .expect("GFN2 water"),
            )
        })
    });
    group.bench_function("gfn2_methane", |b| {
        b.iter(|| {
            black_box(
                sci_form::compute_xtb(black_box(&methane_z), black_box(&methane_xyz))
                    .expect("GFN2 methane"),
            )
        })
    });

    group.finish();
}

fn configure() -> Criterion {
    Criterion::default()
        .warm_up_time(Duration::from_millis(300))
        .measurement_time(Duration::from_secs(3))
        .sample_size(50)
}

criterion_group! {
    name = benches;
    config = configure();
    targets = parse, embed, force_fields, quantum
}
criterion_main!(benches);
