//! GPU Memory Budget — WebGPU/WASM limit enforcement.
//!
//! WebGPU imposes strict hardware limits that vary by platform. This module
//! tracks buffer allocations against these limits and provides pre-flight
//! checks before GPU dispatch to prevent silent failures or OOM.
//!
//! # Contract with callers
//!
//! The CPU fallback path in every GPU module of this crate **always** works:
//! these budget checks are advisory gates in front of an *optional* GPU
//! dispatch, never a hard dependency. When a check fails the caller is
//! expected to fall back to the CPU implementation rather than to fail the
//! whole computation.
//!
//! # WebGPU Default Limits (W3C spec)
//!
//! | Limit | Default | Notes |
//! |-------|---------|-------|
//! | maxStorageBufferBindingSize | 128 MB | Per binding, per shader stage |
//! | maxBufferSize | 256 MB | Total per buffer object |
//! | maxUniformBufferBindingSize | 64 KB | Uniform blocks are small |
//! | maxStorageBuffersPerShaderStage | 8 | Max bindings per stage |
//! | maxComputeWorkgroupStorageSize | 16 KB | Shared workgroup memory |
//! | maxComputeInvocationsPerWorkgroup | 256 | Threads per workgroup |
//! | maxComputeWorkgroupSizeX | 256 | |
//! | maxComputeWorkgroupSizeY | 256 | |
//! | maxComputeWorkgroupSizeZ | 64 | |
//! | maxComputeWorkgroupsPerDimension | 65535 | |
//!
//! The two limits that bite most often in practice are the
//! **256 invocations per workgroup** cap and the **128 MB per storage
//! binding** cap; both are checked by [`GpuMemoryBudget::check_workgroup`]
//! and [`GpuMemoryBudget::check_buffer`] respectively.
//!
//! # Example
//!
//! ```
//! use sci_form::gpu::memory_budget::GpuMemoryBudget;
//!
//! let budget = GpuMemoryBudget::webgpu_default();
//! // 1 MB buffer is well inside the 128 MB per-binding limit.
//! assert!(budget.check_buffer(1024 * 1024).is_ok());
//! // A 300 MB buffer exceeds the 256 MB per-buffer limit.
//! assert!(budget.check_buffer(300_000_000).is_err());
//! ```

/// WebGPU hardware limits, initialized to W3C default minimums.
#[derive(Debug, Clone)]
pub struct GpuMemoryLimits {
    /// Max bytes per storage buffer binding (default: 128 MB).
    pub max_storage_buffer_binding_size: u64,
    /// Max total bytes per buffer object (default: 256 MB).
    pub max_buffer_size: u64,
    /// Max bytes per uniform buffer binding (default: 64 KB).
    pub max_uniform_buffer_binding_size: u64,
    /// Max storage buffers per shader stage (default: 8).
    pub max_storage_buffers_per_stage: u32,
    /// Max shared workgroup memory in bytes (default: 16 KB).
    pub max_workgroup_storage_size: u32,
    /// Max invocations per workgroup (default: 256).
    pub max_invocations_per_workgroup: u32,
    /// Max workgroup size X (default: 256).
    pub max_workgroup_size_x: u32,
    /// Max workgroup size Y (default: 256).
    pub max_workgroup_size_y: u32,
    /// Max workgroup size Z (default: 64).
    pub max_workgroup_size_z: u32,
    /// Max workgroups per dimension (default: 65535).
    pub max_workgroups_per_dimension: u32,
}

impl Default for GpuMemoryLimits {
    fn default() -> Self {
        Self {
            max_storage_buffer_binding_size: 134_217_728, // 128 MB
            max_buffer_size: 268_435_456,                 // 256 MB
            max_uniform_buffer_binding_size: 65_536,      // 64 KB
            max_storage_buffers_per_stage: 8,
            max_workgroup_storage_size: 16_384, // 16 KB
            max_invocations_per_workgroup: 256,
            max_workgroup_size_x: 256,
            max_workgroup_size_y: 256,
            max_workgroup_size_z: 64,
            max_workgroups_per_dimension: 65535,
        }
    }
}

/// Tracks GPU buffer allocations and enforces limits.
///
/// The budget is a plain accounting object: `allocated_bytes` and
/// `storage_bindings_used` are updated by the caller (the buffer manager)
/// after a successful allocation, and are not tied to any particular device.
#[derive(Debug, Clone)]
pub struct GpuMemoryBudget {
    /// Hardware limits this budget enforces, in bytes for the size limits.
    pub limits: GpuMemoryLimits,
    /// Total bytes currently allocated across all buffers.
    pub allocated_bytes: u64,
    /// Number of storage buffer bindings in use.
    pub storage_bindings_used: u32,
}

impl GpuMemoryBudget {
    /// Creates a budget that enforces `limits` and starts with nothing allocated.
    ///
    /// # Example
    ///
    /// ```
    /// use sci_form::gpu::memory_budget::{GpuMemoryBudget, GpuMemoryLimits};
    ///
    /// let budget = GpuMemoryBudget::new(GpuMemoryLimits::default());
    /// assert_eq!(budget.allocated_bytes, 0);
    /// assert_eq!(budget.storage_bindings_used, 0);
    /// ```
    pub fn new(limits: GpuMemoryLimits) -> Self {
        Self {
            limits,
            allocated_bytes: 0,
            storage_bindings_used: 0,
        }
    }

    /// WebGPU default limits.
    pub fn webgpu_default() -> Self {
        Self::new(GpuMemoryLimits::default())
    }

    /// Check if a buffer allocation would exceed limits.
    ///
    /// # Errors
    ///
    /// - [`MemoryError::BufferTooLarge`] when `size_bytes` exceeds
    ///   `limits.max_buffer_size` (default 256 MB).
    /// - [`MemoryError::BindingTooLarge`] when `size_bytes` exceeds
    ///   `limits.max_storage_buffer_binding_size` (default 128 MB).
    pub fn check_buffer(&self, size_bytes: u64) -> Result<(), MemoryError> {
        if size_bytes > self.limits.max_buffer_size {
            return Err(MemoryError::BufferTooLarge {
                requested: size_bytes,
                max: self.limits.max_buffer_size,
            });
        }
        if size_bytes > self.limits.max_storage_buffer_binding_size {
            return Err(MemoryError::BindingTooLarge {
                requested: size_bytes,
                max: self.limits.max_storage_buffer_binding_size,
            });
        }
        Ok(())
    }

    /// Check if a storage binding can be added.
    ///
    /// # Errors
    ///
    /// [`MemoryError::TooManyBindings`] when `storage_bindings_used` has
    /// already reached `limits.max_storage_buffers_per_stage` (default 8).
    pub fn check_storage_binding(&self) -> Result<(), MemoryError> {
        if self.storage_bindings_used >= self.limits.max_storage_buffers_per_stage {
            return Err(MemoryError::TooManyBindings {
                current: self.storage_bindings_used,
                max: self.limits.max_storage_buffers_per_stage,
            });
        }
        Ok(())
    }

    /// Check workgroup configuration validity.
    ///
    /// `size` is the workgroup size in threads per axis and `count` the
    /// number of workgroups dispatched along each axis.
    ///
    /// # Errors
    ///
    /// - [`MemoryError::WorkgroupSizeExceeded`] when any axis of `size`
    ///   exceeds the per-axis limit (defaults 256 × 256 × 64).
    /// - [`MemoryError::TooManyInvocations`] when `size[0] * size[1] * size[2]`
    ///   exceeds `limits.max_invocations_per_workgroup` (default 256).
    /// - [`MemoryError::TooManyWorkgroups`] when any entry of `count` exceeds
    ///   `limits.max_workgroups_per_dimension` (default 65535).
    pub fn check_workgroup(&self, size: [u32; 3], count: [u32; 3]) -> Result<(), MemoryError> {
        if size[0] > self.limits.max_workgroup_size_x
            || size[1] > self.limits.max_workgroup_size_y
            || size[2] > self.limits.max_workgroup_size_z
        {
            return Err(MemoryError::WorkgroupSizeExceeded {
                requested: size,
                max: [
                    self.limits.max_workgroup_size_x,
                    self.limits.max_workgroup_size_y,
                    self.limits.max_workgroup_size_z,
                ],
            });
        }
        let total = size[0] as u64 * size[1] as u64 * size[2] as u64;
        if total > self.limits.max_invocations_per_workgroup as u64 {
            return Err(MemoryError::TooManyInvocations {
                requested: total as u32,
                max: self.limits.max_invocations_per_workgroup,
            });
        }
        for (i, &c) in count.iter().enumerate() {
            if c > self.limits.max_workgroups_per_dimension {
                return Err(MemoryError::TooManyWorkgroups {
                    dimension: i as u32,
                    requested: c,
                    max: self.limits.max_workgroups_per_dimension,
                });
            }
        }
        Ok(())
    }

    /// Pre-flight check for an orbital grid computation.
    ///
    /// Returns the estimated total GPU memory needed in bytes.
    ///
    /// # Errors
    ///
    /// - [`MemoryError::BufferTooLarge`] or [`MemoryError::BindingTooLarge`]
    ///   when one of the basis, MO-coefficient, primitive or output buffers
    ///   exceeds the per-buffer / per-binding limits.
    /// - [`MemoryError::TooManyBindings`] when the pipeline needs more than
    ///   `limits.max_storage_buffers_per_stage` bindings (it needs 5).
    pub fn estimate_orbital_grid_memory(
        &self,
        n_basis: usize,
        grid_points: usize,
        max_primitives_per_basis: usize,
    ) -> Result<OrbitalGridMemoryEstimate, MemoryError> {
        // Buffers:
        // 1. basis functions: n_basis × 32 bytes (BasisFunc struct)
        // 2. MO coefficients: n_basis × 4 bytes (f32)
        // 3. primitives: n_basis × max_primitives × 8 bytes (vec2<f32>)
        // 4. grid params: 32 bytes (uniform)
        // 5. output: grid_points × 4 bytes (f32)
        let basis_bytes = (n_basis * 32) as u64;
        let mo_bytes = (n_basis * 4) as u64;
        let prim_bytes = (n_basis * max_primitives_per_basis * 8) as u64;
        let params_bytes = 32u64;
        let output_bytes = (grid_points * 4) as u64;

        let total = basis_bytes + mo_bytes + prim_bytes + params_bytes + output_bytes;

        // Check individual buffers
        self.check_buffer(basis_bytes)?;
        self.check_buffer(mo_bytes)?;
        self.check_buffer(prim_bytes)?;
        self.check_buffer(output_bytes)?;

        // Need 5 bindings: 3 storage read, 1 uniform, 1 storage read-write
        if 5 > self.limits.max_storage_buffers_per_stage {
            return Err(MemoryError::TooManyBindings {
                current: 5,
                max: self.limits.max_storage_buffers_per_stage,
            });
        }

        Ok(OrbitalGridMemoryEstimate {
            basis_bytes,
            mo_coefficients_bytes: mo_bytes,
            primitives_bytes: prim_bytes,
            params_bytes,
            output_bytes,
            total_bytes: total,
            fits_in_webgpu: total <= self.limits.max_buffer_size,
        })
    }

    /// Estimate memory for D4 dispersion pairwise computation.
    ///
    /// # Errors
    ///
    /// [`MemoryError::BufferTooLarge`] or [`MemoryError::BindingTooLarge`]
    /// when the O(N²) pairwise output buffer (`n_atoms * n_atoms * 4` bytes)
    /// exceeds the per-buffer / per-binding limits.
    pub fn estimate_d4_dispersion_memory(
        &self,
        n_atoms: usize,
    ) -> Result<PairwiseMemoryEstimate, MemoryError> {
        // Buffers:
        // 1. positions: n_atoms × 16 bytes (vec4<f32>: x,y,z,Z)
        // 2. params: n_atoms × 32 bytes (D4 params)
        // 3. config: 32 bytes (uniform)
        // 4. output pairwise: n_atoms × n_atoms × 4 bytes (f32)
        // 5. output energy: 4 bytes (f32 reduction)
        let pos_bytes = (n_atoms * 16) as u64;
        let params_bytes = (n_atoms * 32) as u64;
        let config_bytes = 32u64;
        let pairwise_bytes = (n_atoms * n_atoms * 4) as u64;
        let output_bytes = (n_atoms * 4) as u64;

        let total = pos_bytes + params_bytes + config_bytes + pairwise_bytes + output_bytes;

        self.check_buffer(pairwise_bytes)?;

        Ok(PairwiseMemoryEstimate {
            positions_bytes: pos_bytes,
            params_bytes,
            pairwise_bytes,
            total_bytes: total,
            fits_in_webgpu: total <= self.limits.max_buffer_size,
            max_atoms_for_limit: ((self.limits.max_storage_buffer_binding_size / 4) as f64).sqrt()
                as usize,
        })
    }

    /// Compute optimal workgroup dispatch for a 3D grid.
    ///
    /// Returns `(workgroup_size, workgroup_count)` for the requested grid
    /// dimensions `[nx, ny, nz]`. The workgroup size is clamped to the
    /// per-axis limits, so the result is always dispatchable.
    pub fn optimal_grid_dispatch(&self, dims: [u32; 3]) -> ([u32; 3], [u32; 3]) {
        let wg_size = [
            8u32.min(self.limits.max_workgroup_size_x),
            8u32.min(self.limits.max_workgroup_size_y),
            4u32.min(self.limits.max_workgroup_size_z),
        ];
        let wg_count = [
            dims[0].div_ceil(wg_size[0]),
            dims[1].div_ceil(wg_size[1]),
            dims[2].div_ceil(wg_size[2]),
        ];
        (wg_size, wg_count)
    }

    /// Compute optimal workgroup dispatch for a 1D array.
    ///
    /// Returns `(workgroup_size, workgroup_count)` for `n` elements. The
    /// workgroup size is clamped to `limits.max_workgroup_size_x`.
    pub fn optimal_1d_dispatch(&self, n: u32) -> (u32, u32) {
        let wg_size = 64u32.min(self.limits.max_workgroup_size_x);
        let wg_count = n.div_ceil(wg_size);
        (wg_size, wg_count)
    }
}

/// Memory estimation for orbital grid dispatch.
///
/// Every field is a **byte count** for one GPU buffer of the orbital-grid
/// pipeline (basis data, MO coefficients, primitives, uniform parameters and
/// the scalar-field output).
#[derive(Debug, Clone)]
pub struct OrbitalGridMemoryEstimate {
    /// Bytes for the basis-function buffer (`n_basis * 32` bytes).
    pub basis_bytes: u64,
    /// Bytes for the molecular-orbital coefficient buffer (`n_basis * 4` bytes, `f32`).
    pub mo_coefficients_bytes: u64,
    /// Bytes for the primitive-exponent buffer
    /// (`n_basis * max_primitives_per_basis * 8` bytes, packed `vec2<f32>`).
    pub primitives_bytes: u64,
    /// Bytes for the uniform parameter block (constant 32 bytes).
    pub params_bytes: u64,
    /// Bytes for the scalar-field output buffer (`grid_points * 4` bytes, `f32`).
    pub output_bytes: u64,
    /// Sum of all buffer sizes above, in bytes.
    pub total_bytes: u64,
    /// `true` when `total_bytes` fits in a single WebGPU buffer
    /// (i.e. is at most `limits.max_buffer_size`).
    pub fits_in_webgpu: bool,
}

/// Memory estimation for pairwise computations (D4, EEQ Coulomb).
///
/// Every field is a **byte count** for one GPU buffer of the pairwise
/// (O(N²)) kernel.
#[derive(Debug, Clone)]
pub struct PairwiseMemoryEstimate {
    /// Bytes for the atom-position buffer (`n_atoms * 16` bytes, `vec4<f32>` holding x, y, z and Z).
    pub positions_bytes: u64,
    /// Bytes for the per-atom parameter buffer (`n_atoms * 32` bytes).
    pub params_bytes: u64,
    /// Bytes for the pairwise output buffer (`n_atoms * n_atoms * 4` bytes, `f32`).
    pub pairwise_bytes: u64,
    /// Sum of all buffer sizes above, in bytes.
    pub total_bytes: u64,
    /// `true` when `total_bytes` fits in a single WebGPU buffer.
    pub fits_in_webgpu: bool,
    /// Maximum atom count that fits in one storage buffer.
    pub max_atoms_for_limit: usize,
}

/// Errors from GPU memory budget checks.
///
/// Every variant carries the offending request and the limit that was
/// violated, so a caller can decide between tiling the work and falling back
/// to the CPU path. [`MemoryError::hint`] returns actionable guidance.
#[derive(Debug, Clone)]
// The variant set may grow in a future release; downstream crates must
// match with a wildcard arm.
#[non_exhaustive]
pub enum MemoryError {
    /// A single buffer object would exceed `max_buffer_size` (default 256 MB).
    BufferTooLarge {
        /// Size of the requested buffer, in bytes.
        requested: u64,
        /// Maximum permitted buffer size, in bytes.
        max: u64,
    },
    /// A single storage-buffer binding would exceed
    /// `max_storage_buffer_binding_size` (default 128 MB).
    BindingTooLarge {
        /// Size of the requested binding, in bytes.
        requested: u64,
        /// Maximum permitted binding size, in bytes.
        max: u64,
    },
    /// The dispatch would need more storage-buffer bindings than the shader
    /// stage allows (default 8).
    TooManyBindings {
        /// Number of bindings the dispatch needs.
        current: u32,
        /// Maximum number of storage buffers per shader stage.
        max: u32,
    },
    /// The workgroup dimensions exceed the per-axis workgroup size limits
    /// (defaults 256 × 256 × 64).
    WorkgroupSizeExceeded {
        /// Requested workgroup size as `[x, y, z]`.
        requested: [u32; 3],
        /// Maximum permitted workgroup size as `[x, y, z]`.
        max: [u32; 3],
    },
    /// The product of the workgroup dimensions exceeds the per-workgroup
    /// invocation limit (default 256).
    TooManyInvocations {
        /// Number of invocations the workgroup would spawn.
        requested: u32,
        /// Maximum permitted invocations per workgroup.
        max: u32,
    },
    /// A dispatch dimension would need more workgroups than WebGPU permits
    /// (default 65535 per dimension).
    TooManyWorkgroups {
        /// Index of the offending dimension (0 = x, 1 = y, 2 = z).
        dimension: u32,
        /// Number of workgroups requested along that dimension.
        requested: u32,
        /// Maximum permitted workgroups along that dimension.
        max: u32,
    },
}

impl MemoryError {
    /// Actionable guidance for low-VRAM / limit failures.
    ///
    /// The returned string is appended to the [`std::fmt::Display`] output of
    /// the error and names the usual remedy: tile the work, pack the buffers
    /// or fall back to the CPU path.
    pub fn hint(&self) -> &'static str {
        match self {
            MemoryError::BufferTooLarge { .. } => {
                "Split the work into tiles (e.g. `mmff94_row_tiles`) or fall back to the CPU path."
            }
            MemoryError::BindingTooLarge { .. } => {
                "Reduce the per-binding payload, tile the dispatch, or fall back to the CPU path."
            }
            MemoryError::TooManyBindings { .. } => {
                "Pack several arrays into one buffer or split the kernel into passes."
            }
            MemoryError::WorkgroupSizeExceeded { .. } => {
                "Use `GpuMemoryBudget::optimal_1d_dispatch` / `optimal_grid_dispatch` for the device limits."
            }
            MemoryError::TooManyInvocations { .. } => {
                "Shrink the workgroup size (the WebGPU default is 256 invocations)."
            }
            MemoryError::TooManyWorkgroups { .. } => {
                "Split the dispatch into multiple launches; each dimension allows 65535 workgroups by default."
            }
        }
    }
}

impl std::fmt::Display for MemoryError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MemoryError::BufferTooLarge { requested, max } => {
                write!(f, "Buffer {requested} bytes exceeds max {max} bytes")
            }
            MemoryError::BindingTooLarge { requested, max } => {
                write!(f, "Binding {requested} bytes exceeds max {max} bytes")
            }
            MemoryError::TooManyBindings { current, max } => {
                write!(f, "Need {current} bindings, max {max}")
            }
            MemoryError::WorkgroupSizeExceeded { requested, max } => {
                write!(
                    f,
                    "Workgroup [{},{},{}] exceeds max [{},{},{}]",
                    requested[0], requested[1], requested[2], max[0], max[1], max[2]
                )
            }
            MemoryError::TooManyInvocations { requested, max } => {
                write!(f, "{requested} invocations exceeds max {max}")
            }
            MemoryError::TooManyWorkgroups {
                dimension,
                requested,
                max,
            } => {
                write!(
                    f,
                    "Dimension {dimension}: {requested} workgroups exceeds max {max}"
                )
            }
        }?;
        write!(f, " — {}", self.hint())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_webgpu_defaults() {
        let limits = GpuMemoryLimits::default();
        assert_eq!(limits.max_storage_buffer_binding_size, 128 * 1024 * 1024);
        assert_eq!(limits.max_buffer_size, 256 * 1024 * 1024);
        assert_eq!(limits.max_uniform_buffer_binding_size, 64 * 1024);
    }

    #[test]
    fn test_buffer_check_within_limits() {
        let budget = GpuMemoryBudget::webgpu_default();
        assert!(budget.check_buffer(1024 * 1024).is_ok()); // 1 MB
    }

    #[test]
    fn test_buffer_check_exceeds_limits() {
        let budget = GpuMemoryBudget::webgpu_default();
        assert!(budget.check_buffer(300_000_000).is_err()); // 300 MB > 256 MB
    }

    #[test]
    fn test_orbital_grid_small_molecule() {
        let budget = GpuMemoryBudget::webgpu_default();
        // H₂O: 7 basis functions, 50³ grid = 125000 points
        let est = budget.estimate_orbital_grid_memory(7, 125_000, 3).unwrap();
        assert!(est.fits_in_webgpu);
        assert!(est.total_bytes < 1_000_000); // < 1 MB
    }

    #[test]
    fn test_workgroup_check_valid() {
        let budget = GpuMemoryBudget::webgpu_default();
        assert!(budget.check_workgroup([8, 8, 4], [100, 100, 50]).is_ok());
    }

    #[test]
    fn test_workgroup_check_too_large() {
        let budget = GpuMemoryBudget::webgpu_default();
        assert!(budget.check_workgroup([512, 1, 1], [1, 1, 1]).is_err());
    }

    #[test]
    fn test_d4_memory_small_system() {
        let budget = GpuMemoryBudget::webgpu_default();
        let est = budget.estimate_d4_dispersion_memory(100).unwrap();
        assert!(est.fits_in_webgpu);
    }

    #[test]
    fn test_d4_max_atoms_calculable() {
        let budget = GpuMemoryBudget::webgpu_default();
        let est = budget.estimate_d4_dispersion_memory(10).unwrap();
        // sqrt(128MB / 4) ≈ 5792 atoms
        assert!(est.max_atoms_for_limit > 5000);
    }

    #[test]
    fn test_optimal_grid_dispatch() {
        let budget = GpuMemoryBudget::webgpu_default();
        let (wg_size, wg_count) = budget.optimal_grid_dispatch([64, 64, 64]);
        assert_eq!(wg_size, [8, 8, 4]);
        assert_eq!(wg_count, [8, 8, 16]);
    }

    #[test]
    fn test_optimal_1d_dispatch() {
        let budget = GpuMemoryBudget::webgpu_default();
        let (wg_size, wg_count) = budget.optimal_1d_dispatch(1000);
        assert_eq!(wg_size, 64);
        assert_eq!(wg_count, 16); // ceil(1000/64)
    }

    #[test]
    fn test_memory_errors_carry_actionable_hints() {
        let budget = GpuMemoryBudget::webgpu_default();
        let err = budget
            .check_buffer(budget.limits.max_buffer_size + 1)
            .unwrap_err();
        let text = err.to_string();
        assert!(text.contains("exceeds max"), "got {text:?}");
        assert!(
            text.contains("tiles") || text.contains("CPU"),
            "got {text:?}"
        );

        let err = budget.check_workgroup([4096, 1, 1], [1, 1, 1]).unwrap_err();
        assert!(err.hint().contains("optimal"), "got {:?}", err.hint());
    }
}
