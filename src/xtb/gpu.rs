//! GPU-accelerated xTB gamma-matrix and Hamiltonian construction.
//!
//! Delegates to the existing `gpu::gamma_matrix_gpu` infrastructure for
//! the SCC Coulomb matrix when `experimental-gpu` is enabled.

#[cfg(feature = "experimental-gpu")]
use crate::gpu::gamma_matrix_gpu::build_gamma_gpu;

#[cfg(feature = "experimental-gpu")]
use crate::gpu::context::GpuContext;

#[cfg(feature = "experimental-gpu")]
use nalgebra::DMatrix;

/// Build the xTB damped Coulomb gamma matrix on GPU.
///
/// γ_AB = 1 / sqrt((1/η_A + 1/η_B)² + R²)
///
/// Uses the generic gamma_matrix_gpu kernel from the GPU module.
#[cfg(feature = "experimental-gpu")]
pub fn build_xtb_gamma_gpu(
    ctx: &GpuContext,
    eta: &[f64],
    positions_bohr: &[[f64; 3]],
) -> Result<DMatrix<f64>, String> {
    build_gamma_gpu(ctx, eta, positions_bohr)
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_xtb_gpu_module_compiles() {
        assert!(true);
    }
}
