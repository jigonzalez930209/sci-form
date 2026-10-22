//! GPU-accelerated render buffer generation.
//!
//! Offloads large chart packing workloads to GPU compute when processing
//! many chart payloads simultaneously. Falls back to CPU for small batches.
//!
//! NOTE: this is currently a **CPU fallback only** — the WGSL kernel it
//! advertises does not exist yet, so the "GPU" path executes exactly the same
//! CPU code and the parity test compares CPU against CPU. The audit
//! (reports/auditoria-total.md, P2-32) flagged this; it is kept as a stable
//! call site for the kernel to be added.
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
/// Strategy: packing charts into RecordBatch is mostly memory-copy work;
/// GPU benefit comes from parallel column construction for very large payloads.
/// Current implementation: CPU fallback.
pub fn pack_charts_gpu(ctx: &GpuContext, charts: &[ChartPayload]) -> Vec<RecordBatch> {
    if !ctx.capabilities.gpu_available || charts.len() < GPU_RENDER_THRESHOLD {
        return charts.iter().map(pack_chart_payload).collect();
    }

    // GPU path: currently CPU fallback.
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
