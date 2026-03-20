//! GPU buffer management for input/output storage buffers.
//!
//! Manages the lifecycle of GPU-side buffers: allocation, upload (CPU→GPU),
//! dispatch, and readback (GPU→CPU). In CPU-fallback mode, buffers are
//! simply `Vec<f32>` wrappers.

use super::aligned_types::{GpuAtom, GpuBasisFunction, GpuDenseMatrix};

/// Identifies a buffer's purpose in the compute pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BufferRole {
    /// Atom positions and parameters (read-only on GPU).
    Atoms,
    /// Basis function descriptors (read-only on GPU).
    BasisFunctions,
    /// Gaussian primitives (read-only on GPU).
    Primitives,
    /// Overlap matrix S (read-write on GPU).
    OverlapMatrix,
    /// Core Hamiltonian H⁰ (read-write on GPU).
    HamiltonianMatrix,
    /// Fock matrix F (read-write on GPU).
    FockMatrix,
    /// Density matrix P (read-write on GPU).
    DensityMatrix,
    /// MO coefficient matrix C (read-write on GPU).
    CoefficientMatrix,
    /// Eigenvalues / orbital energies (read-write on GPU).
    Eigenvalues,
    /// Uniform parameters (small, read-only).
    Uniform,
}

/// A managed GPU buffer (CPU-fallback: Vec<u8>).
#[derive(Debug)]
pub struct ManagedBuffer {
    /// Role/purpose of this buffer.
    pub role: BufferRole,
    /// Raw byte data (CPU-side copy).
    pub data: Vec<u8>,
    /// Size in bytes.
    pub size: usize,
    /// Whether the buffer has been modified and needs re-upload.
    pub dirty: bool,
}

impl ManagedBuffer {
    /// Create a new buffer with the given role and byte data.
    pub fn new(role: BufferRole, data: Vec<u8>) -> Self {
        let size = data.len();
        Self {
            role,
            data,
            size,
            dirty: true,
        }
    }

    /// Create a zero buffer of the given size.
    pub fn zeros(role: BufferRole, size_bytes: usize) -> Self {
        Self {
            role,
            data: vec![0u8; size_bytes],
            size: size_bytes,
            dirty: true,
        }
    }

    /// Interpret the buffer as a slice of f32 values.
    pub fn as_f32_slice(&self) -> &[f32] {
        let ptr = self.data.as_ptr() as *const f32;
        let len = self.data.len() / std::mem::size_of::<f32>();
        unsafe { std::slice::from_raw_parts(ptr, len) }
    }

    /// Interpret the buffer as a mutable slice of f32 values.
    pub fn as_f32_slice_mut(&mut self) -> &mut [f32] {
        self.dirty = true;
        let ptr = self.data.as_mut_ptr() as *mut f32;
        let len = self.data.len() / std::mem::size_of::<f32>();
        unsafe { std::slice::from_raw_parts_mut(ptr, len) }
    }
}

/// Manages all buffers needed for a single SCF calculation.
#[derive(Debug)]
pub struct BufferManager {
    /// All allocated buffers indexed by role.
    buffers: Vec<ManagedBuffer>,
    /// Total bytes allocated.
    total_bytes: usize,
}

impl BufferManager {
    /// Create a new empty buffer manager.
    pub fn new() -> Self {
        Self {
            buffers: Vec::new(),
            total_bytes: 0,
        }
    }

    /// Allocate all buffers needed for an SCF calculation.
    ///
    /// Creates storage for: atoms, basis functions, S, H, F, P, C matrices.
    pub fn allocate_scf_buffers(&mut self, n_atoms: usize, n_basis: usize) {
        let matrix_bytes = n_basis * n_basis * std::mem::size_of::<f32>();

        // Input buffers
        self.add_buffer(ManagedBuffer::zeros(
            BufferRole::Atoms,
            n_atoms * std::mem::size_of::<GpuAtom>(),
        ));
        self.add_buffer(ManagedBuffer::zeros(
            BufferRole::BasisFunctions,
            n_basis * std::mem::size_of::<GpuBasisFunction>(),
        ));

        // Matrix buffers
        self.add_buffer(ManagedBuffer::zeros(BufferRole::OverlapMatrix, matrix_bytes));
        self.add_buffer(ManagedBuffer::zeros(BufferRole::HamiltonianMatrix, matrix_bytes));
        self.add_buffer(ManagedBuffer::zeros(BufferRole::FockMatrix, matrix_bytes));
        self.add_buffer(ManagedBuffer::zeros(BufferRole::DensityMatrix, matrix_bytes));
        self.add_buffer(ManagedBuffer::zeros(BufferRole::CoefficientMatrix, matrix_bytes));
        self.add_buffer(ManagedBuffer::zeros(
            BufferRole::Eigenvalues,
            n_basis * std::mem::size_of::<f32>(),
        ));
    }

    /// Add a buffer to the manager.
    fn add_buffer(&mut self, buffer: ManagedBuffer) {
        self.total_bytes += buffer.size;
        self.buffers.push(buffer);
    }

    /// Get a buffer by its role.
    pub fn get(&self, role: BufferRole) -> Option<&ManagedBuffer> {
        self.buffers.iter().find(|b| b.role == role)
    }

    /// Get a mutable buffer by its role.
    pub fn get_mut(&mut self, role: BufferRole) -> Option<&mut ManagedBuffer> {
        self.buffers.iter_mut().find(|b| b.role == role)
    }

    /// Upload atom data to the atoms buffer.
    pub fn upload_atoms(&mut self, atoms: &[GpuAtom]) {
        if let Some(buf) = self.get_mut(BufferRole::Atoms) {
            let src = unsafe {
                std::slice::from_raw_parts(
                    atoms.as_ptr() as *const u8,
                    atoms.len() * std::mem::size_of::<GpuAtom>(),
                )
            };
            buf.data[..src.len()].copy_from_slice(src);
            buf.dirty = true;
        }
    }

    /// Upload a dense matrix to a specific buffer role.
    pub fn upload_matrix(&mut self, role: BufferRole, matrix: &GpuDenseMatrix) {
        if let Some(buf) = self.get_mut(role) {
            let bytes = matrix.as_bytes();
            buf.data[..bytes.len()].copy_from_slice(bytes);
            buf.dirty = true;
        }
    }

    /// Read back a matrix from a buffer into a GpuDenseMatrix.
    pub fn readback_matrix(&self, role: BufferRole, rows: usize, cols: usize) -> Option<GpuDenseMatrix> {
        self.get(role).map(|buf| {
            let n = rows * cols;
            let f32_data: Vec<f32> = buf.as_f32_slice()[..n].to_vec();
            GpuDenseMatrix {
                data: f32_data,
                rows,
                cols,
            }
        })
    }

    /// Total memory allocated in bytes.
    pub fn total_allocated_bytes(&self) -> usize {
        self.total_bytes
    }

    /// Number of buffers allocated.
    pub fn num_buffers(&self) -> usize {
        self.buffers.len()
    }
}

impl Default for BufferManager {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_buffer_allocation() {
        let mut mgr = BufferManager::new();
        mgr.allocate_scf_buffers(10, 40);

        assert_eq!(mgr.num_buffers(), 8);
        assert!(mgr.get(BufferRole::OverlapMatrix).is_some());
        assert!(mgr.get(BufferRole::DensityMatrix).is_some());
        assert!(mgr.total_allocated_bytes() > 0);
    }

    #[test]
    fn test_matrix_upload_readback() {
        let mut mgr = BufferManager::new();
        mgr.allocate_scf_buffers(2, 3);

        let mut mat = GpuDenseMatrix::zeros(3, 3);
        mat.set(0, 0, 1.0);
        mat.set(1, 1, 1.0);
        mat.set(2, 2, 1.0);

        mgr.upload_matrix(BufferRole::OverlapMatrix, &mat);

        let readback = mgr.readback_matrix(BufferRole::OverlapMatrix, 3, 3).unwrap();
        assert!((readback.get(0, 0) - 1.0).abs() < 1e-6);
        assert!((readback.get(1, 1) - 1.0).abs() < 1e-6);
        assert!((readback.get(0, 1)).abs() < 1e-6);
    }
}
