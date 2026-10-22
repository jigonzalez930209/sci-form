//! GPU buffer lifecycle management for quantum chemistry workloads.
//!
//! Manages allocation, tracking, and release of GPU storage/uniform buffers,
//! enforcing WebGPU binding limits (128 MB per storage buffer, 256 MB device max).
//!
//! # Contract with callers
//!
//! The manager is a *book-keeping* layer: it owns the CPU-side view of which
//! buffers a calculation needs and how many bytes they take, and it refuses
//! allocations that would break the WebGPU limits (128 MB per storage
//! binding) or the caller-supplied budget. It never touches a device, so it
//! is usable without a GPU — the CPU fallback path of this crate uses exactly
//! the same buffer plan, just without the device-side allocation.
//!
//! # Example
//!
//! ```
//! use sci_form::gpu::buffer_manager::{BufferManager, BufferRole};
//!
//! let mut mgr = BufferManager::new(256); // 256 MB budget
//! let idx = mgr.allocate(BufferRole::OverlapMatrix, 1024, "overlap").unwrap();
//! assert_eq!(idx, 0);
//! assert_eq!(mgr.total_allocated, 1024);
//! ```

use super::aligned_types::{GpuMatrixParams, GpuScfParams};

/// Buffer role classification for slot tracking and reuse.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
// The variant set may grow in a future release; downstream crates must
// match with a wildcard arm.
#[non_exhaustive]
pub enum BufferRole {
    /// Atom position + parameter data (read-only).
    Atoms,
    /// Basis function descriptors (read-only).
    BasisFunctions,
    /// Gaussian primitive contraction table (read-only).
    Primitives,
    /// Overlap matrix S (read-write).
    OverlapMatrix,
    /// Kinetic energy matrix T (read-write).
    KineticMatrix,
    /// Nuclear attraction matrix V (read-write).
    NuclearMatrix,
    /// Core Hamiltonian H = T + V (read-write).
    CoreHamiltonian,
    /// Density matrix P (read-write).
    DensityMatrix,
    /// Fock matrix F (read-write).
    FockMatrix,
    /// Two-electron integral buffer (read-write, large).
    TwoElectronIntegrals,
    /// SCF parameters uniform (uniform).
    ScfParams,
    /// Matrix dimension parameters (uniform).
    MatrixParams,
    /// ESP grid output (read-write).
    EspGrid,
    /// Orbital grid output (read-write).
    OrbitalGrid,
}

impl BufferRole {
    /// Whether this buffer is a uniform buffer (vs storage buffer).
    pub fn is_uniform(&self) -> bool {
        matches!(self, BufferRole::ScfParams | BufferRole::MatrixParams)
    }
}

/// CPU-side representation of a managed GPU buffer.
///
/// The struct records the *plan* for a buffer (role, size, label) and
/// optionally keeps the CPU shadow copy that will be uploaded to (or was
/// read back from) the device.
#[derive(Debug, Clone)]
pub struct ManagedBuffer {
    /// Role this buffer plays in the calculation (determines uniform vs storage).
    pub role: BufferRole,
    /// Size of the buffer in bytes.
    pub size_bytes: usize,
    /// Debug label of the buffer.
    pub label: String,
    /// `true` while the buffer is mapped (i.e. readable/writable) on the device.
    pub is_mapped: bool,
    /// CPU shadow copy (for readback or initial upload).
    pub cpu_data: Option<Vec<u8>>,
}

impl ManagedBuffer {
    /// Creates an empty managed buffer of `size_bytes` bytes.
    ///
    /// The buffer starts unmapped and without a CPU shadow copy.
    pub fn new(role: BufferRole, size_bytes: usize, label: &str) -> Self {
        Self {
            role,
            size_bytes,
            label: label.to_string(),
            is_mapped: false,
            cpu_data: None,
        }
    }

    /// Creates a managed buffer pre-filled with `data`.
    ///
    /// `size_bytes` is set to `data.len()`, and the payload is kept as the CPU
    /// shadow copy that will be uploaded to the device.
    pub fn with_data(role: BufferRole, data: Vec<u8>, label: &str) -> Self {
        let size = data.len();
        Self {
            role,
            size_bytes: size,
            label: label.to_string(),
            is_mapped: false,
            cpu_data: Some(data),
        }
    }
}

/// Max storage buffer binding size (WebGPU default).
const MAX_STORAGE_BUFFER_SIZE: usize = 128 * 1024 * 1024; // 128 MB

/// Manages the lifecycle of all GPU buffers for a quantum chemistry calculation.
///
/// The manager tracks the buffers it has handed out (by index) and the total
/// number of bytes they occupy, so a caller can build a bind group without
/// re-deriving the sizes.
pub struct BufferManager {
    /// All buffers handed out so far, in allocation order.
    pub buffers: Vec<ManagedBuffer>,
    /// Total number of bytes allocated across all buffers.
    pub total_allocated: usize,
    /// Maximum number of bytes that may be allocated, in bytes.
    pub max_budget: usize,
}

impl BufferManager {
    /// Creates a manager with a budget of `max_budget_mb` megabytes.
    ///
    /// `max_budget_mb` is converted to bytes (× 1024 × 1024); a value of 0
    /// therefore forbids every allocation.
    pub fn new(max_budget_mb: usize) -> Self {
        Self {
            buffers: Vec::new(),
            total_allocated: 0,
            max_budget: max_budget_mb * 1024 * 1024,
        }
    }

    /// Allocate a new managed buffer. Returns error if it exceeds WebGPU limits.
    ///
    /// Returns the index of the new buffer inside [`BufferManager::buffers`].
    ///
    /// # Errors
    ///
    /// - `"Buffer '{label}' ({n} MB) exceeds WebGPU maxStorageBufferBindingSize (128 MB)"`
    ///   when a non-uniform buffer is larger than the 128 MB per-binding limit.
    /// - `"Buffer '{label}' ({n} bytes) would exceed GPU memory budget ({m} MB)"`
    ///   when the allocation would push `total_allocated` past `max_budget`.
    pub fn allocate(
        &mut self,
        role: BufferRole,
        size_bytes: usize,
        label: &str,
    ) -> Result<usize, String> {
        if !role.is_uniform() && size_bytes > MAX_STORAGE_BUFFER_SIZE {
            return Err(format!(
                "Buffer '{}' ({} MB) exceeds WebGPU maxStorageBufferBindingSize (128 MB)",
                label,
                size_bytes / (1024 * 1024)
            ));
        }

        if self.total_allocated + size_bytes > self.max_budget {
            return Err(format!(
                "Buffer '{}' ({} bytes) would exceed GPU memory budget ({} MB)",
                label,
                size_bytes,
                self.max_budget / (1024 * 1024)
            ));
        }

        let buf = ManagedBuffer::new(role, size_bytes, label);
        let idx = self.buffers.len();
        self.total_allocated += size_bytes;
        self.buffers.push(buf);
        Ok(idx)
    }

    /// Allocate a buffer with initial CPU data.
    ///
    /// Returns the index of the new buffer inside [`BufferManager::buffers`].
    ///
    /// # Errors
    ///
    /// - `"Buffer '{label}' ({n} MB) exceeds WebGPU maxStorageBufferBindingSize"`
    ///   when a non-uniform buffer is larger than the 128 MB per-binding limit.
    /// - `"Allocation would exceed budget: {used} + {size} > {max}"`
    ///   when the allocation would push `total_allocated` past `max_budget`.
    pub fn allocate_with_data(
        &mut self,
        role: BufferRole,
        data: Vec<u8>,
        label: &str,
    ) -> Result<usize, String> {
        let size = data.len();
        if !role.is_uniform() && size > MAX_STORAGE_BUFFER_SIZE {
            return Err(format!(
                "Buffer '{}' ({} MB) exceeds WebGPU maxStorageBufferBindingSize",
                label,
                size / (1024 * 1024)
            ));
        }
        if self.total_allocated + size > self.max_budget {
            return Err(format!(
                "Allocation would exceed budget: {} + {} > {}",
                self.total_allocated, size, self.max_budget
            ));
        }

        let buf = ManagedBuffer::with_data(role, data, label);
        let idx = self.buffers.len();
        self.total_allocated += size;
        self.buffers.push(buf);
        Ok(idx)
    }

    /// Compute buffer size for a dense n×n matrix of f32.
    ///
    /// Returns `n_basis * n_basis * 4` bytes.
    pub fn matrix_size(n_basis: usize) -> usize {
        n_basis * n_basis * std::mem::size_of::<f32>()
    }

    /// Plan all buffers needed for an SCF calculation.
    ///
    /// Reserves the atom, basis and primitive input buffers, the six
    /// `n_basis × n_basis` matrices (S, T, V, H, P, F) and the two uniform
    /// parameter blocks. The returned [`ScfBufferPlan`] holds the buffer
    /// indices and the total number of bytes reserved.
    ///
    /// # Errors
    ///
    /// Propagates the error of the first failing [`BufferManager::allocate`]
    /// call, i.e. a 128 MB per-binding violation or a budget overflow.
    pub fn plan_scf_buffers(
        &mut self,
        n_atoms: usize,
        n_basis: usize,
        n_primitives: usize,
    ) -> Result<ScfBufferPlan, String> {
        let atom_size = n_atoms * std::mem::size_of::<super::aligned_types::GpuAtom>();
        let basis_size = n_basis * std::mem::size_of::<super::aligned_types::GpuBasisFunction>();
        let prim_size =
            n_primitives * std::mem::size_of::<super::aligned_types::GpuGaussianPrimitive>();
        let mat_size = Self::matrix_size(n_basis);

        let atoms = self.allocate(BufferRole::Atoms, atom_size, "atoms")?;
        let basis = self.allocate(BufferRole::BasisFunctions, basis_size, "basis_functions")?;
        let prims = self.allocate(BufferRole::Primitives, prim_size, "primitives")?;
        let overlap = self.allocate(BufferRole::OverlapMatrix, mat_size, "overlap_matrix")?;
        let kinetic = self.allocate(BufferRole::KineticMatrix, mat_size, "kinetic_matrix")?;
        let nuclear = self.allocate(BufferRole::NuclearMatrix, mat_size, "nuclear_matrix")?;
        let core_h = self.allocate(BufferRole::CoreHamiltonian, mat_size, "core_hamiltonian")?;
        let density = self.allocate(BufferRole::DensityMatrix, mat_size, "density_matrix")?;
        let fock = self.allocate(BufferRole::FockMatrix, mat_size, "fock_matrix")?;

        let params_size = std::mem::size_of::<GpuMatrixParams>();
        let scf_params_size = std::mem::size_of::<GpuScfParams>();
        let mat_params = self.allocate(BufferRole::MatrixParams, params_size, "matrix_params")?;
        let scf_params = self.allocate(BufferRole::ScfParams, scf_params_size, "scf_params")?;

        Ok(ScfBufferPlan {
            atoms,
            basis,
            prims,
            overlap,
            kinetic,
            nuclear,
            core_h,
            density,
            fock,
            mat_params,
            scf_params,
            total_bytes: self.total_allocated,
        })
    }

    /// Release all buffers and reset allocation tracking.
    pub fn release_all(&mut self) {
        self.buffers.clear();
        self.total_allocated = 0;
    }

    /// Summary of all allocated buffers for debugging.
    ///
    /// Returns a multi-line human-readable report: one header line with the
    /// buffer count and the used/available budget in MB, followed by one
    /// indented line per buffer with its index, role, label and size in bytes.
    pub fn summary(&self) -> String {
        let mut s = format!(
            "BufferManager: {} buffers, {:.2} MB / {:.2} MB\n",
            self.buffers.len(),
            self.total_allocated as f64 / (1024.0 * 1024.0),
            self.max_budget as f64 / (1024.0 * 1024.0),
        );
        for (i, buf) in self.buffers.iter().enumerate() {
            s.push_str(&format!(
                "  [{}] {:?} '{}' — {} bytes\n",
                i, buf.role, buf.label, buf.size_bytes
            ));
        }
        s
    }
}

/// Indices into the BufferManager for all SCF-related buffers.
///
/// The six matrix buffers (overlap, kinetic, nuclear, core Hamiltonian,
/// density, Fock) each hold a dense `n_basis × n_basis` block of `f32`
/// matrix elements in row-major order; the two parameter blocks hold the
/// uniform structs consumed by the kernels.
#[derive(Debug, Clone)]
pub struct ScfBufferPlan {
    /// Index of the atom buffer (positions and per-atom parameters).
    pub atoms: usize,
    /// Index of the basis-function descriptor buffer.
    pub basis: usize,
    /// Index of the Gaussian primitive contraction table.
    pub prims: usize,
    /// Index of the overlap matrix **S** buffer.
    pub overlap: usize,
    /// Index of the kinetic energy matrix **T** buffer.
    pub kinetic: usize,
    /// Index of the nuclear attraction matrix **V** buffer.
    pub nuclear: usize,
    /// Index of the core Hamiltonian **H = T + V** buffer.
    pub core_h: usize,
    /// Index of the density matrix **P** buffer.
    pub density: usize,
    /// Index of the Fock matrix **F** buffer.
    pub fock: usize,
    /// Index of the matrix-dimension uniform buffer.
    pub mat_params: usize,
    /// Index of the SCF-parameter uniform buffer.
    pub scf_params: usize,
    /// Total number of bytes reserved for all buffers of the plan.
    pub total_bytes: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_buffer_allocation() {
        let mut mgr = BufferManager::new(256);
        let idx = mgr
            .allocate(BufferRole::OverlapMatrix, 1024, "test_overlap")
            .unwrap();
        assert_eq!(idx, 0);
        assert_eq!(mgr.total_allocated, 1024);
    }

    #[test]
    fn test_budget_overflow() {
        let mut mgr = BufferManager::new(1); // 1 MB budget
        let result = mgr.allocate(
            BufferRole::TwoElectronIntegrals,
            2 * 1024 * 1024,
            "too_large",
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_scf_plan() {
        let mut mgr = BufferManager::new(256);
        let plan = mgr.plan_scf_buffers(10, 20, 60).unwrap();
        assert!(plan.total_bytes > 0);
        assert_eq!(mgr.buffers.len(), 11); // 3 input + 6 matrices + 2 params
    }

    #[test]
    fn test_release_all() {
        let mut mgr = BufferManager::new(256);
        mgr.allocate(BufferRole::Atoms, 512, "atoms").unwrap();
        mgr.allocate(BufferRole::OverlapMatrix, 1024, "overlap")
            .unwrap();
        assert_eq!(mgr.buffers.len(), 2);
        mgr.release_all();
        assert_eq!(mgr.buffers.len(), 0);
        assert_eq!(mgr.total_allocated, 0);
    }

    #[test]
    fn test_webgpu_storage_limit() {
        let mut mgr = BufferManager::new(512);
        // 200 MB exceeds the 128 MB per-binding limit
        let result = mgr.allocate(
            BufferRole::TwoElectronIntegrals,
            200 * 1024 * 1024,
            "huge_eri",
        );
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("maxStorageBufferBindingSize"));
    }
}
