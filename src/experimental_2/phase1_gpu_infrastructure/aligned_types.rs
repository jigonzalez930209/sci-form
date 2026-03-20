//! GPU-aligned data structures for zero-copy transfer to VRAM.
//!
//! All structures use `#[repr(C)]` to guarantee deterministic memory layout
//! compatible with WebGPU buffer binding requirements. Fields are padded to
//! 16-byte alignment boundaries as required by WGSL `uniform`/`storage` buffers.
//!
//! # GPU Alignment Rules (WebGPU/WGSL)
//!
//! | Type    | Alignment | Size  |
//! |---------|-----------|-------|
//! | f32     | 4 bytes   | 4 B   |
//! | vec2<f32>| 8 bytes  | 8 B   |
//! | vec3<f32>| 16 bytes | 12 B  |
//! | vec4<f32>| 16 bytes | 16 B  |
//! | mat4x4  | 16 bytes | 64 B  |
//!
//! All `vec3` fields are stored as `[f32; 4]` (vec4) to satisfy alignment.

/// Atom data packed for GPU storage buffer.
///
/// Each atom's data is exactly 32 bytes (2 × vec4<f32>) for coalesced access.
/// In WGSL:
/// ```wgsl
/// struct GpuAtom {
///     position: vec4<f32>,  // xyz = position (Bohr), w = atomic_number
///     params: vec4<f32>,    // x = zeta_s, y = zeta_p, z = eta, w = n_valence
/// };
/// ```
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct GpuAtom {
    /// [x, y, z, atomic_number] — position in Bohr, 4th component is Z as f32.
    pub position_z: [f32; 4],
    /// [zeta_s, zeta_p, eta, n_valence] — orbital exponents and chemical hardness.
    pub params: [f32; 4],
}

impl GpuAtom {
    pub fn new(x: f64, y: f64, z: f64, atomic_number: u8, zeta_s: f64, zeta_p: f64, eta: f64, n_valence: f64) -> Self {
        Self {
            position_z: [x as f32, y as f32, z as f32, atomic_number as f32],
            params: [zeta_s as f32, zeta_p as f32, eta as f32, n_valence as f32],
        }
    }

    pub fn position(&self) -> [f64; 3] {
        [self.position_z[0] as f64, self.position_z[1] as f64, self.position_z[2] as f64]
    }

    pub fn atomic_number(&self) -> u8 {
        self.position_z[3] as u8
    }
}

/// Basis function descriptor for GPU.
///
/// Each basis function is 48 bytes (3 × vec4<f32>).
/// ```wgsl
/// struct GpuBasisFunction {
///     center: vec4<f32>,       // xyz = center (Bohr), w = atom_index
///     quantum: vec4<f32>,      // x = n, y = l, z = m, w = zeta
///     energy_pad: vec4<f32>,   // x = vsip (eV), y = contraction_start, z = n_primitives, w = pad
/// };
/// ```
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct GpuBasisFunction {
    /// [x, y, z, atom_index] — center in Bohr.
    pub center_atom: [f32; 4],
    /// [n, l, m, zeta] — quantum numbers and Slater exponent.
    pub quantum_numbers: [f32; 4],
    /// [vsip, contraction_offset, n_primitives, 0.0] — energy and contraction info.
    pub energy_contraction: [f32; 4],
}

/// Gaussian primitive for GPU contraction evaluation.
///
/// 8 bytes per primitive — packed into a flat array referenced by offset.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct GpuGaussianPrimitive {
    /// Gaussian exponent α (Bohr⁻²).
    pub alpha: f32,
    /// Contraction coefficient.
    pub coefficient: f32,
}

/// Matrix dimensions and dispatch parameters for GPU compute shaders.
///
/// 16 bytes = 1 × vec4<u32>.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct GpuMatrixParams {
    /// Number of basis functions (matrix dimension).
    pub n_basis: u32,
    /// Number of atoms.
    pub n_atoms: u32,
    /// Number of occupied orbitals.
    pub n_occupied: u32,
    /// Padding to 16-byte alignment.
    pub _pad: u32,
}

/// SCF iteration parameters passed as a uniform buffer.
///
/// 32 bytes = 2 × vec4<f32>.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct GpuScfParams {
    /// [k_wh, damping_factor, convergence_threshold, iteration]
    pub scf_control: [f32; 4],
    /// [energy_current, energy_previous, max_density_change, 0.0]
    pub scf_state: [f32; 4],
}

/// Dense matrix stored in row-major order for GPU.
///
/// Wraps a flat `Vec<f32>` with dimension metadata.
#[derive(Debug, Clone)]
pub struct GpuDenseMatrix {
    /// Row-major flat data.
    pub data: Vec<f32>,
    /// Number of rows.
    pub rows: usize,
    /// Number of columns.
    pub cols: usize,
}

impl GpuDenseMatrix {
    /// Create a zero-initialized matrix.
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {
            data: vec![0.0f32; rows * cols],
            rows,
            cols,
        }
    }

    /// Get element at (i, j).
    pub fn get(&self, i: usize, j: usize) -> f32 {
        self.data[i * self.cols + j]
    }

    /// Set element at (i, j).
    pub fn set(&mut self, i: usize, j: usize, value: f32) {
        self.data[i * self.cols + j] = value;
    }

    /// Byte slice for GPU upload.
    pub fn as_bytes(&self) -> &[u8] {
        let ptr = self.data.as_ptr() as *const u8;
        let len = self.data.len() * std::mem::size_of::<f32>();
        unsafe { std::slice::from_raw_parts(ptr, len) }
    }
}

/// Conversion utility to pack a `MolecularSystem` into GPU-aligned atom array.
pub fn pack_atoms_for_gpu(
    system: &super::super::types::MolecularSystem,
    zeta_s: &[f64],
    zeta_p: &[f64],
    eta: &[f64],
    n_valence: &[f64],
) -> Vec<GpuAtom> {
    system
        .positions_bohr
        .iter()
        .enumerate()
        .map(|(i, pos)| {
            GpuAtom::new(
                pos[0],
                pos[1],
                pos[2],
                system.atomic_numbers[i],
                zeta_s.get(i).copied().unwrap_or(1.0),
                zeta_p.get(i).copied().unwrap_or(1.0),
                eta.get(i).copied().unwrap_or(0.0),
                n_valence.get(i).copied().unwrap_or(0.0),
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gpu_atom_alignment() {
        assert_eq!(std::mem::size_of::<GpuAtom>(), 32);
        assert_eq!(std::mem::align_of::<GpuAtom>(), 4);
    }

    #[test]
    fn test_gpu_basis_function_alignment() {
        assert_eq!(std::mem::size_of::<GpuBasisFunction>(), 48);
    }

    #[test]
    fn test_gpu_matrix_params_alignment() {
        assert_eq!(std::mem::size_of::<GpuMatrixParams>(), 16);
    }

    #[test]
    fn test_dense_matrix() {
        let mut mat = GpuDenseMatrix::zeros(3, 3);
        mat.set(1, 2, 3.14);
        assert!((mat.get(1, 2) - 3.14).abs() < 1e-5);
        assert_eq!(mat.as_bytes().len(), 36); // 9 * 4 bytes
    }

    #[test]
    fn test_gpu_atom_roundtrip() {
        let atom = GpuAtom::new(1.0, 2.0, 3.0, 6, 1.625, 1.625, 10.0, 4.0);
        assert_eq!(atom.atomic_number(), 6);
        let pos = atom.position();
        assert!((pos[0] - 1.0).abs() < 1e-5);
    }
}
