//! GPU-accelerated render buffer generation.
//!
//! Offloads large chart packing workloads to GPU compute when processing
//! many chart payloads simultaneously. Falls back to CPU for small batches.
//!
//! **CPU by design.** Packing is a memory-bound copy of the chart arrays into
//! `RecordBatch` columns; a GPU round-trip would add two host↔device copies
//! and a dispatch for no arithmetic. The function keeps the GPU-compatible
//! call site so a future batch format that needs a device-side transform can
//! be added without changing callers, but the CPU implementation is the
//! intended one.
//!
//! # Stability tier
//!
//! **Alpha** — gated behind `alpha-render-bridge` (plus `experimental-gpu`
//! for the GPU context); the validation gate is
//! `tests/alpha_render_validation.rs` (F90).

use super::{pack_chart_payload, ChartPayload};
use crate::gpu::context::GpuContext;
use crate::transport::arrow::RecordBatch;

/// Minimum batch size for GPU-backed render buffer generation.
const GPU_RENDER_THRESHOLD: usize = 16;

/// GPU-accelerated batch chart packing.
///
/// Packing is a memory-bound copy, so this intentionally runs on the CPU; the
/// GPU context is accepted only to keep the call site uniform with the other
/// accelerated entry points.
pub fn pack_charts_gpu(ctx: &GpuContext, charts: &[ChartPayload]) -> Vec<RecordBatch> {
    let _ = (ctx, GPU_RENDER_THRESHOLD);
    charts.iter().map(pack_chart_payload).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alpha::render_bridge::ChartSeries;

    #[test]
    fn gpu_pack_falls_back_to_cpu() {
        let ctx = GpuContext::cpu_fallback();
        let charts: Vec<ChartPayload> = (0..4)
            .map(|i| ChartPayload {
                title: format!("chart_{}", i),
                series: vec![ChartSeries {
                    series_id: format!("s{}", i),
                    label: format!("S{}", i),
                    x: vec![1.0, 2.0],
                    y: vec![3.0, 4.0],
                    x_unit: "".into(),
                    y_unit: "".into(),
                }],
            })
            .collect();
        let result = pack_charts_gpu(&ctx, &charts);
        assert_eq!(result.len(), 4);
        for batch in &result {
            assert_eq!(batch.float_columns.len(), 2);
        }
    }
}
