//! Render-bridge payload contracts — core since 2026-10-14 (F90).
//!
//! Visualization-facing payloads decoupled from solver internals so Python,
//! WASM, and charting frontends share a stable format. Validated against an
//! oracle-authored JSON Schema (draft 2020-12) plus adversarial payloads.
//!
//! Transport decisions:
//! - Small results: JSON via `serde_json` (< 64 KiB).
//! - Large profile/DOS arrays: Arrow-style `RecordBatch` via `crate::transport::arrow`.
//! - WASM typed arrays: `Float64Array` / `Float32Array` via flat contiguous buffers.

use serde::{Deserialize, Serialize};

use crate::transport::arrow::RecordBatch;
/// One XY series payload for charts.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct ChartSeries {
    /// Stable series identifier.
    pub series_id: String,
    /// Human-readable label.
    pub label: String,
    /// X values.
    pub x: Vec<f64>,
    /// Y values.
    pub y: Vec<f64>,
    /// Suggested unit label for the x axis.
    pub x_unit: String,
    /// Suggested unit label for the y axis.
    pub y_unit: String,
}

/// Multi-series chart payload.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct ChartPayload {
    /// Chart title.
    pub title: String,
    /// Included series.
    pub series: Vec<ChartSeries>,
}

/// One 3D polyline frame, typically for reaction trajectories.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct TrajectoryFrame {
    /// Time-like coordinate.
    pub frame_time: f64,
    /// Flattened xyz coordinates.
    pub coords: Vec<f64>,
}

/// Generic bridge payload for render clients.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct RenderPayload {
    /// Scalar chart data.
    pub charts: Vec<ChartPayload>,
    /// Optional geometry frames.
    pub trajectory: Vec<TrajectoryFrame>,
}

/// Convert GSM path coordinates into trajectory frames for 3D visualization.
pub fn gsm_path_trajectory(
    path_coords: &[Vec<f64>],
    path_energies: &[f64],
) -> Vec<TrajectoryFrame> {
    path_coords
        .iter()
        .zip(path_energies)
        .map(|(coords, &energy)| TrajectoryFrame {
            frame_time: energy, // Use energy as the "time" coordinate for path visualization
            coords: coords.clone(),
        })
        .collect()
}

/// Pack a generic `ChartPayload` into a `RecordBatch`.
pub fn pack_chart_payload(chart: &ChartPayload) -> RecordBatch {
    let mut batch = RecordBatch::new();
    for (i, s) in chart.series.iter().enumerate() {
        let n = s.x.len();
        batch.add_float64(&format!("{}_x", s.series_id), s.x.clone(), vec![n]);
        batch.add_float64(&format!("{}_y", s.series_id), s.y.clone(), vec![n]);
        let _ = i; // series index for diagnostics
    }
    batch
}

/// Pack multiple chart payloads into `RecordBatch` in parallel using rayon.
#[cfg(feature = "parallel")]
pub fn pack_charts_parallel(charts: &[ChartPayload]) -> Vec<RecordBatch> {
    use rayon::prelude::*;
    charts.par_iter().map(pack_chart_payload).collect()
}

/// Pack multiple chart payloads into `RecordBatch` sequentially (non-parallel fallback).
#[cfg(not(feature = "parallel"))]
pub fn pack_charts_parallel(charts: &[ChartPayload]) -> Vec<RecordBatch> {
    charts.iter().map(pack_chart_payload).collect()
}

/// Verify lossless round-trip: `ChartPayload` → `RecordBatch` → reconstructed series.
/// Returns `true` when all values survive within `tol`.
pub fn verify_chart_round_trip(chart: &ChartPayload, tol: f64) -> bool {
    let batch = pack_chart_payload(chart);
    for s in &chart.series {
        let x_col = batch
            .float_columns
            .iter()
            .find(|c| c.name == format!("{}_x", s.series_id));
        let y_col = batch
            .float_columns
            .iter()
            .find(|c| c.name == format!("{}_y", s.series_id));
        match (x_col, y_col) {
            (Some(xc), Some(yc)) => {
                if xc.values.len() != s.x.len() || yc.values.len() != s.y.len() {
                    return false;
                }
                for (a, b) in xc.values.iter().zip(&s.x) {
                    if (a - b).abs() > tol {
                        return false;
                    }
                }
                for (a, b) in yc.values.iter().zip(&s.y) {
                    if (a - b).abs() > tol {
                        return false;
                    }
                }
            }
            _ => return false,
        }
    }
    true
}

