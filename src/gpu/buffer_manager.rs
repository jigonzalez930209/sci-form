//! GPU buffer lifecycle management for quantum chemistry workloads.
//!
//! Manages allocation, tracking, and release of GPU storage/uniform buffers,
//! enforcing WebGPU binding limits (128 MB per storage buffer, 256 MB device max).

use super::aligned_types::{GpuMatrixParams, GpuScfParams};

/// Buffer role classification for slot tracking and reuse.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
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
#[derive(Debug, Clone)]
pub struct ManagedBuffer {
    pub role: BufferRole,
    pub size_bytes: usize,
    pub label: String,
    pub is_mapped: bool,
    /// CPU shadow copy (for readback or initial upload).
    pub cpu_data: Option<Vec<u8>>,
}

impl ManagedBuffer {
    pub fn new(role: BufferRole, size_bytes: usize, label: &str) -> Self {
        Self {
            role,
            size_bytes,
            label: label.to_string(),
            is_mapped: false,
            cpu_data: None,
        }
    }

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
pub struct BufferManager {
    pub buffers: Vec<ManagedBuffer>,
    pub total_allocated: usize,
    pub max_budget: usize,
}

impl BufferManager {
    pub fn new(max_budget_mb: usize) -> Self {
        Self {
            buffers: Vec::new(),
            total_allocated: 0,
            max_budget: max_budget_mb * 1024 * 1024,
        }
    }

    /// Allocate a new managed buffer. Returns error if it exceeds WebGPU limits.
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
    pub fn matrix_size(n_basis: usize) -> usize {
        n_basis * n_basis * std::mem::size_of::<f32>()
    }

    /// Plan all buffers needed for an SCF calculation.
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
#[derive(Debug, Clone)]
pub struct ScfBufferPlan {
    pub atoms: usize,
    pub basis: usize,
    pub prims: usize,
    pub overlap: usize,
    pub kinetic: usize,
    pub nuclear: usize,
    pub core_h: usize,
    pub density: usize,
    pub fock: usize,
    pub mat_params: usize,
    pub scf_params: usize,
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
