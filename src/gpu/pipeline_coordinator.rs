//! GPU Pipeline Coordinator — multi-kernel dispatch with memory-aware scheduling.
//!
//! Coordinates sequential GPU dispatches for multi-step pipelines
//! (e.g., orbital grid → marching cubes → isosurface). Enforces
//! memory budget checks before each dispatch and provides tier-aware
//! scheduling decisions.

use super::memory_budget::{GpuMemoryBudget, MemoryError};
use super::shader_registry::ShaderDescriptor;

/// A single step in a GPU computation pipeline.
#[derive(Debug, Clone)]
pub struct PipelineStep {
    /// Human-readable label for this step.
    pub label: String,
    /// Reference to the shader descriptor.
    pub shader_name: &'static str,
    /// Estimated buffer sizes for each binding (bytes).
    pub buffer_sizes: Vec<u64>,
    /// Grid dimensions to dispatch over, or element count for 1D.
    pub dispatch: DispatchShape,
}

/// How to dispatch a compute shader.
#[derive(Debug, Clone)]
pub enum DispatchShape {
    /// 1D dispatch over `n` elements (shader divides by workgroup size).
    Linear(u32),
    /// 3D dispatch over a grid [nx, ny, nz].
    Grid3D([u32; 3]),
}

impl DispatchShape {
    /// Compute workgroup count given a workgroup size.
    pub fn workgroup_count(&self, workgroup_size: [u32; 3]) -> [u32; 3] {
        match self {
            DispatchShape::Linear(n) => [
                (*n + workgroup_size[0] - 1) / workgroup_size[0],
                1,
                1,
            ],
            DispatchShape::Grid3D(dims) => [
                (dims[0] + workgroup_size[0] - 1) / workgroup_size[0],
                (dims[1] + workgroup_size[1] - 1) / workgroup_size[1],
                (dims[2] + workgroup_size[2] - 1) / workgroup_size[2],
            ],
        }
    }
}

/// Outcome of validating a pipeline step.
#[derive(Debug, Clone)]
pub enum StepDecision {
    /// GPU dispatch is feasible and recommended.
    GpuDispatch {
        workgroup_count: [u32; 3],
        total_buffer_bytes: u64,
    },
    /// GPU dispatch possible but CPU may be faster for this tier/size.
    CpuPreferred {
        reason: String,
    },
    /// Must fall back to CPU — exceeds GPU limits.
    CpuRequired {
        error: MemoryError,
    },
}

/// Coordinates multi-step GPU pipelines.
#[derive(Debug)]
pub struct PipelineCoordinator {
    budget: GpuMemoryBudget,
    gpu_available: bool,
}

impl PipelineCoordinator {
    pub fn new(budget: GpuMemoryBudget, gpu_available: bool) -> Self {
        Self { budget, gpu_available }
    }

    /// Create a coordinator that always falls back to CPU.
    pub fn cpu_only() -> Self {
        Self {
            budget: GpuMemoryBudget::webgpu_default(),
            gpu_available: false,
        }
    }

    /// Evaluate whether a pipeline step should run on GPU or CPU.
    pub fn evaluate_step(
        &self,
        step: &PipelineStep,
        shader: &ShaderDescriptor,
    ) -> StepDecision {
        if !self.gpu_available {
            return StepDecision::CpuPreferred {
                reason: "No GPU available".to_string(),
            };
        }

        // Check tier recommendation
        if !shader.tier.gpu_recommended() {
            return StepDecision::CpuPreferred {
                reason: format!("{} — {}", shader.tier, shader.description),
            };
        }

        // Check buffer sizes against budget
        let total_bytes: u64 = step.buffer_sizes.iter().sum();
        for (i, &size) in step.buffer_sizes.iter().enumerate() {
            if let Err(e) = self.budget.check_buffer(size) {
                return StepDecision::CpuRequired {
                    error: e,
                };
            }
            // Also check storage binding limit for buffer i
            if i as u32 >= self.budget.limits.max_storage_buffers_per_stage {
                return StepDecision::CpuRequired {
                    error: MemoryError::TooManyBindings {
                        current: i as u32 + 1,
                        max: self.budget.limits.max_storage_buffers_per_stage,
                    },
                };
            }
        }

        // Compute workgroup count and validate
        let wg_count = step.dispatch.workgroup_count(shader.workgroup_size);
        if let Err(e) = self.budget.check_workgroup(shader.workgroup_size, wg_count) {
            return StepDecision::CpuRequired {
                error: e,
            };
        }

        StepDecision::GpuDispatch {
            workgroup_count: wg_count,
            total_buffer_bytes: total_bytes,
        }
    }

    /// Validate an entire multi-step pipeline and return decisions for each step.
    pub fn plan_pipeline(
        &self,
        steps: &[PipelineStep],
        shaders: &[&ShaderDescriptor],
    ) -> PipelinePlan {
        assert_eq!(steps.len(), shaders.len(), "Steps and shaders must match");

        let decisions: Vec<(String, StepDecision)> = steps
            .iter()
            .zip(shaders.iter())
            .map(|(step, shader)| {
                let decision = self.evaluate_step(step, shader);
                (step.label.clone(), decision)
            })
            .collect();

        let gpu_steps = decisions
            .iter()
            .filter(|(_, d)| matches!(d, StepDecision::GpuDispatch { .. }))
            .count();
        let cpu_steps = decisions.len() - gpu_steps;

        PipelinePlan {
            decisions,
            gpu_steps,
            cpu_steps,
        }
    }

    /// Determine if a large pairwise computation should be chunked.
    ///
    /// Returns chunk sizes for atoms when the pairwise O(N²) matrix
    /// exceeds the storage buffer binding limit.
    pub fn compute_chunks(&self, n_atoms: usize, bytes_per_pair: u64) -> Vec<(usize, usize)> {
        let max_pairs = self.budget.limits.max_storage_buffer_binding_size / bytes_per_pair;
        let max_chunk = (max_pairs as f64).sqrt() as usize;

        if n_atoms <= max_chunk {
            return vec![(0, n_atoms)];
        }

        let mut chunks = Vec::new();
        let mut start = 0;
        while start < n_atoms {
            let end = (start + max_chunk).min(n_atoms);
            chunks.push((start, end));
            start = end;
        }
        chunks
    }
}

/// Result of planning an entire pipeline.
#[derive(Debug)]
pub struct PipelinePlan {
    /// Decision for each step: (label, decision).
    pub decisions: Vec<(String, StepDecision)>,
    /// Number of steps that will run on GPU.
    pub gpu_steps: usize,
    /// Number of steps that will fall back to CPU.
    pub cpu_steps: usize,
}

impl PipelinePlan {
    /// Whether all steps will run on GPU.
    pub fn fully_gpu(&self) -> bool {
        self.cpu_steps == 0
    }

    /// Generate a text report of the pipeline plan.
    pub fn report(&self) -> String {
        let mut out = format!(
            "Pipeline Plan: {} GPU / {} CPU steps\n",
            self.gpu_steps, self.cpu_steps
        );
        for (label, decision) in &self.decisions {
            let status = match decision {
                StepDecision::GpuDispatch { workgroup_count, total_buffer_bytes } => {
                    format!(
                        "GPU → wg[{},{},{}], {:.1} KB",
                        workgroup_count[0], workgroup_count[1], workgroup_count[2],
                        *total_buffer_bytes as f64 / 1024.0
                    )
                }
                StepDecision::CpuPreferred { reason } => {
                    format!("CPU (preferred) — {reason}")
                }
                StepDecision::CpuRequired { error } => {
                    format!("CPU (required) — {error}")
                }
            };
            out.push_str(&format!("  [{label}] {status}\n"));
        }
        out
    }
}

// ─── Pre-built pipeline templates ──────────────────────────────────

/// Build a standard orbital visualization pipeline.
///
/// Steps: orbital grid → marching cubes (positive lobe) → marching cubes (negative lobe).
pub fn orbital_visualization_pipeline(
    n_basis: usize,
    grid_dims: [u32; 3],
    max_primitives: usize,
) -> Vec<PipelineStep> {
    let basis_bytes = (n_basis * 32) as u64;
    let mo_bytes = (n_basis * 4) as u64;
    let prim_bytes = (n_basis * max_primitives * 8) as u64;
    let grid_points = grid_dims[0] as u64 * grid_dims[1] as u64 * grid_dims[2] as u64;
    let grid_bytes = grid_points * 4;

    // Estimate worst-case triangle output: 5 triangles per active voxel, 10% active
    let est_triangles = (grid_points / 10) * 5;
    let vertex_bytes = est_triangles * 24; // 24 bytes per Vertex (vec3 pos + vec3 normal)

    vec![
        PipelineStep {
            label: "Orbital Grid".to_string(),
            shader_name: "orbital_grid",
            buffer_sizes: vec![basis_bytes, mo_bytes, prim_bytes, 32, grid_bytes],
            dispatch: DispatchShape::Grid3D(grid_dims),
        },
        PipelineStep {
            label: "Marching Cubes (positive lobe)".to_string(),
            shader_name: "marching_cubes",
            buffer_sizes: vec![
                grid_bytes,        // scalar_field
                256 * 4,           // edge_table
                256 * 16 * 4,      // tri_table
                vertex_bytes,      // vertices output
                4,                 // tri_count atomic
                32,                // params
            ],
            dispatch: DispatchShape::Grid3D([
                grid_dims[0] - 1,
                grid_dims[1] - 1,
                grid_dims[2] - 1,
            ]),
        },
        PipelineStep {
            label: "Marching Cubes (negative lobe)".to_string(),
            shader_name: "marching_cubes",
            buffer_sizes: vec![
                grid_bytes,
                256 * 4,
                256 * 16 * 4,
                vertex_bytes,
                4,
                32,
            ],
            dispatch: DispatchShape::Grid3D([
                grid_dims[0] - 1,
                grid_dims[1] - 1,
                grid_dims[2] - 1,
            ]),
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::shader_registry;

    #[test]
    fn test_dispatch_1d() {
        let shape = DispatchShape::Linear(1000);
        let wg = shape.workgroup_count([64, 1, 1]);
        assert_eq!(wg, [16, 1, 1]); // ceil(1000/64)
    }

    #[test]
    fn test_dispatch_3d() {
        let shape = DispatchShape::Grid3D([64, 64, 64]);
        let wg = shape.workgroup_count([8, 8, 4]);
        assert_eq!(wg, [8, 8, 16]);
    }

    #[test]
    fn test_gpu_dispatch_feasible() {
        let coord = PipelineCoordinator::new(
            GpuMemoryBudget::webgpu_default(),
            true,
        );
        let step = PipelineStep {
            label: "test".to_string(),
            shader_name: "orbital_grid",
            buffer_sizes: vec![1024, 512, 2048, 32, 500_000],
            dispatch: DispatchShape::Grid3D([50, 50, 50]),
        };
        let shader = shader_registry::find_shader("orbital_grid").unwrap();
        let decision = coord.evaluate_step(&step, shader);
        assert!(matches!(decision, StepDecision::GpuDispatch { .. }));
    }

    #[test]
    fn test_cpu_fallback_when_no_gpu() {
        let coord = PipelineCoordinator::cpu_only();
        let step = PipelineStep {
            label: "test".to_string(),
            shader_name: "orbital_grid",
            buffer_sizes: vec![1024],
            dispatch: DispatchShape::Linear(100),
        };
        let shader = shader_registry::find_shader("orbital_grid").unwrap();
        let decision = coord.evaluate_step(&step, shader);
        assert!(matches!(decision, StepDecision::CpuPreferred { .. }));
    }

    #[test]
    fn test_cpu_preferred_for_tier4() {
        let coord = PipelineCoordinator::new(
            GpuMemoryBudget::webgpu_default(),
            true,
        );
        let step = PipelineStep {
            label: "smoke test".to_string(),
            shader_name: "vector_add",
            buffer_sizes: vec![1024, 1024, 32, 1024],
            dispatch: DispatchShape::Linear(256),
        };
        let shader = shader_registry::find_shader("vector_add").unwrap();
        let decision = coord.evaluate_step(&step, shader);
        assert!(matches!(decision, StepDecision::CpuPreferred { .. }));
    }

    #[test]
    fn test_buffer_too_large_forces_cpu() {
        let coord = PipelineCoordinator::new(
            GpuMemoryBudget::webgpu_default(),
            true,
        );
        let step = PipelineStep {
            label: "huge".to_string(),
            shader_name: "orbital_grid",
            buffer_sizes: vec![300_000_000], // 300 MB > 256 MB limit
            dispatch: DispatchShape::Linear(1),
        };
        let shader = shader_registry::find_shader("orbital_grid").unwrap();
        let decision = coord.evaluate_step(&step, shader);
        assert!(matches!(decision, StepDecision::CpuRequired { .. }));
    }

    #[test]
    fn test_pipeline_plan() {
        let coord = PipelineCoordinator::new(
            GpuMemoryBudget::webgpu_default(),
            true,
        );
        let steps = vec![
            PipelineStep {
                label: "grid".to_string(),
                shader_name: "orbital_grid",
                buffer_sizes: vec![1024, 256, 2048, 32, 500_000],
                dispatch: DispatchShape::Grid3D([50, 50, 50]),
            },
            PipelineStep {
                label: "mc".to_string(),
                shader_name: "marching_cubes",
                buffer_sizes: vec![500_000, 1024, 16384, 200_000, 4, 32],
                dispatch: DispatchShape::Grid3D([49, 49, 49]),
            },
        ];
        let shaders: Vec<&ShaderDescriptor> = vec![
            shader_registry::find_shader("orbital_grid").unwrap(),
            shader_registry::find_shader("marching_cubes").unwrap(),
        ];
        let plan = coord.plan_pipeline(&steps, &shaders);
        assert_eq!(plan.gpu_steps, 2);
        assert_eq!(plan.cpu_steps, 0);
        assert!(plan.fully_gpu());
    }

    #[test]
    fn test_compute_chunks_small() {
        let coord = PipelineCoordinator::new(
            GpuMemoryBudget::webgpu_default(),
            true,
        );
        let chunks = coord.compute_chunks(100, 4); // 100 atoms, 4 bytes/pair
        assert_eq!(chunks.len(), 1);
        assert_eq!(chunks[0], (0, 100));
    }

    #[test]
    fn test_compute_chunks_large() {
        let coord = PipelineCoordinator::new(
            GpuMemoryBudget::webgpu_default(),
            true,
        );
        // Force chunking: 1 MB per pair with 128 MB limit → max ~11 atoms per chunk
        let chunks = coord.compute_chunks(50, 1_000_000);
        assert!(chunks.len() > 1);
    }

    #[test]
    fn test_orbital_visualization_pipeline() {
        let steps = orbital_visualization_pipeline(7, [50, 50, 50], 3);
        assert_eq!(steps.len(), 3);
        assert_eq!(steps[0].label, "Orbital Grid");
        assert_eq!(steps[1].label, "Marching Cubes (positive lobe)");
        assert_eq!(steps[2].label, "Marching Cubes (negative lobe)");
    }

    #[test]
    fn test_pipeline_report() {
        let coord = PipelineCoordinator::new(
            GpuMemoryBudget::webgpu_default(),
            true,
        );
        let steps = orbital_visualization_pipeline(7, [30, 30, 30], 3);
        let shaders: Vec<&ShaderDescriptor> = steps
            .iter()
            .map(|s| shader_registry::find_shader(s.shader_name).unwrap())
            .collect();
        let plan = coord.plan_pipeline(&steps, &shaders);
        let report = plan.report();
        assert!(report.contains("Orbital Grid"));
        assert!(report.contains("GPU"));
    }
}
