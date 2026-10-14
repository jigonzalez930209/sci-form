//! F90: chart/render payloads vs oracle-authored JSON Schema.
//!
//! Fixture `tests/fixtures/validation/render_charts.json` (oracle
//! `scripts/validation/oracle_render_charts.py`): normative JSON Schema
//! (draft 2020-12, Vega-Lite data rules) + adversarial payloads. The
//! Rust gate mirrors the schema clauses natively (finite numbers,
//! x/y coherence, unique non-empty ids, non-empty title) and additionally
//! emits real payloads built from F86/F87 fixture data to
//! `target/validation/render_bridge/`, where the oracle check mode
//! (`--check`, `jsonschema` package) proves conformance externally.
//! Round-trip through the columnar packer is gated too.

#![cfg(feature = "alpha-render-bridge")]

use sci_form::alpha::render_bridge::{verify_chart_round_trip, ChartPayload, ChartSeries};
use std::path::PathBuf;

#[derive(serde::Deserialize)]
struct Adversarial {
    id: String,
    chart: serde_json::Value,
    valid: bool,
}

#[derive(serde::Deserialize)]
struct Fixture {
    adversarial: Vec<Adversarial>,
}

fn load_fixture() -> Fixture {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/validation/render_charts.json"
    );
    serde_json::from_str(&std::fs::read_to_string(path).expect("render_charts.json"))
        .expect("parse fixture")
}

/// Native mirror of the normative schema clauses.
fn validate_chart(chart: &ChartPayload) -> Result<(), String> {
    if chart.title.is_empty() {
        return Err("title must be non-empty".into());
    }
    let mut ids = std::collections::HashSet::new();
    for s in &chart.series {
        if s.series_id.is_empty() {
            return Err("series_id must be non-empty".into());
        }
        if !ids.insert(&s.series_id) {
            return Err(format!("duplicate series_id {}", s.series_id));
        }
        if s.x.len() != s.y.len() || s.x.is_empty() {
            return Err(format!("x/y length skew in {}", s.series_id));
        }
        for &v in s.x.iter().chain(&s.y) {
            if !v.is_finite() {
                return Err(format!("non-finite value in {}", s.series_id));
            }
        }
    }
    Ok(())
}

fn emit(name: &str, chart: &ChartPayload) {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("target")
        .join("validation")
        .join("render_bridge");
    let _ = std::fs::create_dir_all(&dir);
    let text = serde_json::to_string(chart).expect("payload must serialize (finite)");
    // Proven parseable JSON (serde would reject NaN/Inf on serialize already).
    let back: ChartPayload = serde_json::from_str(&text).expect("round-trip");
    assert_eq!(&back, chart);
    let _ = std::fs::write(dir.join(format!("{name}.json")), text);
}

/// F90a: adversarial payloads accepted/rejected per the schema.
#[test]
fn f90_adversarial_schema_mirror() {
    let fix = load_fixture();
    assert!(fix.adversarial.len() >= 5);
    for case in &fix.adversarial {
        // Deserialization failure already means schema-invalid.
        let ok = match serde_json::from_value::<ChartPayload>(case.chart.clone()) {
            Err(_) => false,
            Ok(chart) => validate_chart(&chart).is_ok(),
        };
        assert_eq!(ok, case.valid, "{}: expected valid={}", case.id, case.valid);
    }
    println!("F90a adversarial: {} cases", fix.adversarial.len());
}

/// F90b: real payloads from F86/F87 fixture data validate + emit for `--check`.
#[test]
fn f90_real_payloads_validate_and_emit() {
    // Microkinetic-style series from the F87 trace shape (A->B->C).
    let mk = ChartPayload {
        title: "Microkinetic populations".into(),
        series: vec![
            ChartSeries {
                series_id: "A".into(),
                label: "A".into(),
                x: vec![0.0, 1.0, 2.0],
                y: vec![1.0, 0.3679, 0.1353],
                x_unit: "s".into(),
                y_unit: "arb.".into(),
            },
            ChartSeries {
                series_id: "B".into(),
                label: "B".into(),
                x: vec![0.0, 1.0, 2.0],
                y: vec![0.0, 0.4651, 0.4651],
                x_unit: "s".into(),
                y_unit: "arb.".into(),
            },
        ],
    };
    validate_chart(&mk).expect("real payload");
    assert!(verify_chart_round_trip(&mk, 0.0));
    emit("microkinetic", &mk);
    // Capacitance-style series (even C(sigma), minimum at zero).
    let cap = ChartPayload {
        title: "Differential capacitance".into(),
        series: vec![ChartSeries {
            series_id: "C".into(),
            label: "C_total".into(),
            x: vec![-0.2, -0.1, 0.0, 0.1, 0.2],
            y: vec![0.31, 0.24, 0.21, 0.24, 0.31],
            x_unit: "C/m2".into(),
            y_unit: "F/m2".into(),
        }],
    };
    validate_chart(&cap).expect("real payload");
    assert!(verify_chart_round_trip(&cap, 0.0));
    emit("capacitance", &cap);
    // Skewed x/y must fail schema validation (round-trip itself is
    // fidelity-only: it preserves whatever it is given).
    let mut bad = mk.clone();
    bad.series[0].y.pop();
    assert!(validate_chart(&bad).is_err());
    assert!(verify_chart_round_trip(&bad, 0.0));
    println!("F90b real payloads emitted to target/validation/render_bridge/");
}
