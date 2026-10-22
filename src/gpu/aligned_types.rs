//! GPU-aligned data structures for zero-copy transfer to VRAM.
//!
//! All structures use `#[repr(C)]` to guarantee deterministic memory layout
//! compatible with WebGPU buffer binding requirements.
//!
//! # GPU Alignment Rules (WebGPU/WGSL)
//!
//! | Type        | Alignment | Size  |
//! |-------------|-----------|-------|
//! | f32         | 4 bytes   | 4 B   |
//! | `vec2<f32>`  | 8 bytes   | 8 B   |
//! | `vec3<f32>`  | 16 bytes  | 12 B  |
//! | `vec4<f32>`  | 16 bytes  | 16 B  |
//!
//! # Units
//!
//! Positions are in **Bohr** (a₀) and Gaussian exponents in **Bohr⁻²**, which
//! is what the shaders in this crate expect. Energy-like quantities
//! (`energy_current`, `energy_previous`) are in **Hartree**; the SCF
//! convergence threshold and the density change are in Hartree and
//! electrons/bohr³ respectively, as produced by the CPU SCF path.
//!
//! # Contract with callers
//!
//! These types are plain data: they describe the memory layout of the
//! buffers that the optional GPU path uploads. The CPU fallback path always
//! works and does not need them; they are used by the GPU kernels only.
//!
//! # Feature gate
//!
//! The `Pod`/`Zeroable` derives that make the structs directly uploadable
//! are only enabled by the `experimental-gpu-rendering` feature.

#[cfg(feature = "experimental-gpu-rendering")]
use bytemuck::{Pod, Zeroable};

/// Atom data packed for GPU storage buffer (32 bytes = 2 x `vec4<f32>`).
#[repr(C)]
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "experimental-gpu-rendering", derive(Pod, Zeroable))]
pub struct GpuAtom {
    /// [x, y, z, atomic_number] — position in Bohr, 4th component is Z as f32.
    pub position_z: [f32; 4],
    /// [zeta_s, zeta_p, eta, n_valence] — orbital exponents and chemical hardness.
    ///
    /// `zeta_s`/`zeta_p` are Slater exponents in **Bohr⁻¹**, `eta` is the
    /// chemical hardness in **eV** and `n_valence` the number of valence
    /// electrons (dimensionless).
    pub params: [f32; 4],
}

impl GpuAtom {
    /// Packs one atom into the 32-byte GPU layout.
    ///
    /// `x`, `y`, `z` are the Cartesian coordinates in **Bohr**; `zeta_s` and
    /// `zeta_p` are Slater-type orbital exponents in **Bohr⁻¹** (the PM3
    /// convention); `eta` is the chemical hardness in **eV** (the unit the
    /// EEQ/CPM charge-equilibration models of this crate use); `n_valence` is
    /// the number of valence electrons (dimensionless). All values are
    /// narrowed to `f32`.
    ///
    /// # Example
    ///
    /// ```
    /// use sci_form::gpu::aligned_types::GpuAtom;
    ///
    /// let atom = GpuAtom::new(0.0, 0.0, 0.0, 6, 1.625, 1.625, 10.0, 4.0);
    /// assert_eq!(atom.atomic_number(), 6);
    /// ```
    pub fn new(
        x: f64,
        y: f64,
        z: f64,
        atomic_number: u8,
        zeta_s: f64,
        zeta_p: f64,
        eta: f64,
        n_valence: f64,
    ) -> Self {
        Self {
            position_z: [x as f32, y as f32, z as f32, atomic_number as f32],
            params: [zeta_s as f32, zeta_p as f32, eta as f32, n_valence as f32],
        }
    }

    /// Returns the Cartesian position of the atom in **Bohr**.
    pub fn position(&self) -> [f64; 3] {
        [
            self.position_z[0] as f64,
            self.position_z[1] as f64,
            self.position_z[2] as f64,
        ]
    }

    /// Returns the atomic number (proton count, dimensionless).
    pub fn atomic_number(&self) -> u8 {
        self.position_z[3] as u8
    }
}

/// Basis function descriptor for GPU (48 bytes = 3 x `vec4<f32>`).
#[repr(C)]
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "experimental-gpu-rendering", derive(Pod, Zeroable))]
pub struct GpuBasisFunction {
    /// [x, y, z, atom_index] — center in Bohr.
    pub center_atom: [f32; 4],
    /// [n, l, m, zeta] — quantum numbers and Slater exponent.
    pub quantum_numbers: [f32; 4],
    /// [vsip, contraction_offset, n_primitives, 0.0].
    pub energy_contraction: [f32; 4],
}

/// Gaussian primitive for GPU contraction evaluation (8 bytes).
#[repr(C)]
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "experimental-gpu-rendering", derive(Pod, Zeroable))]
pub struct GpuGaussianPrimitive {
    /// Gaussian exponent α (Bohr⁻²).
    pub alpha: f32,
    /// Contraction coefficient (dimensionless).
    pub coefficient: f32,
}

/// Matrix dimensions and dispatch parameters for GPU (16 bytes).
#[repr(C)]
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "experimental-gpu-rendering", derive(Pod, Zeroable))]
pub struct GpuMatrixParams {
    /// Number of basis functions (dimensionless count).
    pub n_basis: u32,
    /// Number of atoms (dimensionless count).
    pub n_atoms: u32,
    /// Number of occupied (doubly filled) spatial orbitals (dimensionless count).
    pub n_occupied: u32,
    /// Padding word, always zero; keeps the struct 16-byte aligned.
    pub _pad: u32,
}

/// SCF iteration parameters as uniform buffer (32 bytes).
///
/// Layout only: no shipped kernel consumes this block yet, but the field
/// order and units follow the CPU SCF implementation so that a GPU SCF
/// driver can reuse them verbatim.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "experimental-gpu-rendering", derive(Pod, Zeroable))]
pub struct GpuScfParams {
    /// `[k_wh, damping_factor, convergence_threshold, iteration]`.
    ///
    /// `k_wh` is the Wolfsberg–Helmholz K factor (dimensionless, typically
    /// 1.75), `damping_factor` the density damping factor (dimensionless,
    /// 0–1), `convergence_threshold` the SCF energy threshold in **Hartree**
    /// and `iteration` the current iteration counter (dimensionless).
    pub scf_control: [f32; 4],
    /// `[energy_current, energy_previous, max_density_change, 0.0]`.
    ///
    /// The two energies are total SCF energies in **Hartree** and
    /// `max_density_change` the largest density-matrix element change
    /// (dimensionless; density matrix elements are dimensionless in the
    /// atomic unit system).
    pub scf_state: [f32; 4],
}

/// Dense matrix stored in row-major order for GPU.
#[derive(Debug, Clone)]
pub struct GpuDenseMatrix {
    /// Matrix elements in row-major order (`rows * cols` values).
    pub data: Vec<f32>,
    /// Number of rows (dimensionless).
    pub rows: usize,
    /// Number of columns (dimensionless).
    pub cols: usize,
}

impl GpuDenseMatrix {
    /// Creates a `rows × cols` matrix filled with zeros.
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {
            data: vec![0.0f32; rows * cols],
            rows,
            cols,
        }
    }

    /// Returns the element at row `i`, column `j`.
    ///
    /// # Panics
    ///
    /// Panics if `i >= rows` or `j >= cols`, because the index would fall
    /// outside the backing `Vec`.
    pub fn get(&self, i: usize, j: usize) -> f32 {
        self.data[i * self.cols + j]
    }

    /// Sets the element at row `i`, column `j` to `value`.
    ///
    /// # Panics
    ///
    /// Panics if `i >= rows` or `j >= cols`, because the index would fall
    /// outside the backing `Vec`.
    pub fn set(&mut self, i: usize, j: usize, value: f32) {
        self.data[i * self.cols + j] = value;
    }

    /// Returns the raw native-endian byte image of the matrix.
    ///
    /// The slice is `rows * cols * 4` bytes long and is meant to be uploaded
    /// verbatim as a storage buffer.
    pub fn as_bytes(&self) -> &[u8] {
        let ptr = self.data.as_ptr() as *const u8;
        let len = self.data.len() * std::mem::size_of::<f32>();
        // SAFETY: `self.data` is a `Vec<f32>` that owns `self.data.len()`
        // contiguous, initialised `f32` values.  The pointer is therefore
        // valid for reads of `len` bytes, and `len` is exactly the number of
        // bytes those values occupy (`size_of::<f32>()` per element, no
        // padding).  The returned slice borrows from `self`, so the buffer
        // cannot be freed or reallocated while the slice is alive.
        unsafe { std::slice::from_raw_parts(ptr, len) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gpu_atom_alignment() {
        assert_eq!(std::mem::size_of::<GpuAtom>(), 32);
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
        mat.set(1, 2, std::f32::consts::PI);
        assert!((mat.get(1, 2) - std::f32::consts::PI).abs() < 1e-5);
        assert_eq!(mat.as_bytes().len(), 36);
    }

    #[test]
    fn test_gpu_atom_roundtrip() {
        let atom = GpuAtom::new(1.0, 2.0, 3.0, 6, 1.625, 1.625, 10.0, 4.0);
        assert_eq!(atom.atomic_number(), 6);
        let pos = atom.position();
        assert!((pos[0] - 1.0).abs() < 1e-5);
    }
}
