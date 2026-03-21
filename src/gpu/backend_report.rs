//! Backend execution reporting for GPU/CPU compute dispatch.
//!
//! Every orbital grid or isosurface execution returns a report indicating
//! which backend was actually used (GPU vs CPU) and why, enabling
//! transparent benchmarking and fallback diagnostics.

/// Coarse activation state for the GPU runtime.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GpuActivationState {
    /// GPU feature is not enabled at compile time.
    FeatureDisabled,
    /// GPU feature enabled but no hardware adapter found.
    NoAdapter,
    /// GPU runtime is fully available and ready for dispatch.
    Ready,
}

/// Activation report describing the current GPU runtime state.
#[derive(Debug, Clone)]
pub struct GpuActivationReport {
    /// Backend name visible to callers (e.g. "Vulkan", "Metal", "CPU-fallback").
    pub backend: String,
    /// Whether the build includes the `experimental-gpu` feature.
    pub feature_enabled: bool,
    /// Whether a real GPU adapter was found and device created.
    pub gpu_available: bool,
    /// Whether the orbital grid GPU pipeline is ready.
    pub runtime_ready: bool,
    /// Coarse activation state.
    pub state: GpuActivationState,
    /// Human-readable explanation.
    pub reason: String,
}

/// Per-execution report for orbital grid evaluation.
#[derive(Debug, Clone)]
pub struct OrbitalGridReport {
    /// Backend that actually executed the computation.
    pub backend: String,
    /// Whether GPU was used for this execution.
    pub used_gpu: bool,
    /// Whether GPU dispatch was attempted before falling back.
    pub attempted_gpu: bool,
    /// Number of grid points evaluated.
    pub n_points: usize,
    /// Human-readable note for logs.
    pub note: String,
}

/// Per-execution report for isosurface extraction.
#[derive(Debug, Clone)]
pub struct IsosurfaceReport {
    /// Backend used for grid evaluation.
    pub grid_backend: String,
    /// Whether GPU was used for the grid evaluation step.
    pub grid_used_gpu: bool,
    /// Number of triangles produced.
    pub n_triangles: usize,
    /// Isovalue used.
    pub isovalue: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_activation_state_enum() {
        let state = GpuActivationState::FeatureDisabled;
        assert_eq!(state, GpuActivationState::FeatureDisabled);
        assert_ne!(state, GpuActivationState::Ready);
    }

    #[test]
    fn test_activation_report_cpu_fallback() {
        let report = GpuActivationReport {
            backend: "CPU-fallback".to_string(),
            feature_enabled: false,
            gpu_available: false,
            runtime_ready: false,
            state: GpuActivationState::FeatureDisabled,
            reason: "Feature not enabled".to_string(),
        };
        assert!(!report.gpu_available);
        assert!(!report.runtime_ready);
    }

    #[test]
    fn test_orbital_grid_report() {
        let report = OrbitalGridReport {
            backend: "CPU".to_string(),
            used_gpu: false,
            attempted_gpu: false,
            n_points: 1000,
            note: "CPU evaluation".to_string(),
        };
        assert_eq!(report.n_points, 1000);
        assert!(!report.used_gpu);
    }
}
