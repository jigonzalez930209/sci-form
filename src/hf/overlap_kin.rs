//! Overlap (S) and kinetic energy (T) integral evaluator.
//!
//! Uses 1D Cartesian factorization and analytical derivatives.

use super::basis::{BasisSet, ShellType};
use nalgebra::DMatrix;

pub fn compute_overlap_matrix(basis: &BasisSet) -> DMatrix<f64> {
    let n = basis.n_basis();
    let mut s_mat = DMatrix::zeros(n, n);
    let mut mu = 0;
    for shell_i in &basis.shells {
        let mut nu = 0;
        for shell_j in &basis.shells {
            let (s_block, _) = compute_shell_pair_st(shell_i, shell_j);
            for i in 0..shell_i.n_functions() {
                for j in 0..shell_j.n_functions() {
                    s_mat[(mu + i, nu + j)] = s_block[i * shell_j.n_functions() + j];
                }
            }
            nu += shell_j.n_functions();
        }
        mu += shell_i.n_functions();
    }
    s_mat
}

pub fn compute_kinetic_matrix(basis: &BasisSet) -> DMatrix<f64> {
    let n = basis.n_basis();
    let mut t_mat = DMatrix::zeros(n, n);
    let mut mu = 0;
    for shell_i in &basis.shells {
        let mut nu = 0;
        for shell_j in &basis.shells {
            let (_, t_block) = compute_shell_pair_st(shell_i, shell_j);
            for i in 0..shell_i.n_functions() {
                for j in 0..shell_j.n_functions() {
                    t_mat[(mu + i, nu + j)] = t_block[i * shell_j.n_functions() + j];
                }
            }
            nu += shell_j.n_functions();
        }
        mu += shell_i.n_functions();
    }
    t_mat
}

fn compute_shell_pair_st(
    shell_a: &super::basis::Shell,
    shell_b: &super::basis::Shell,
) -> (Vec<f64>, Vec<f64>) {
    let na = shell_a.n_functions();
    let nb = shell_b.n_functions();
    let mut s_res = vec![0.0; na * nb];
    let mut t_res = vec![0.0; na * nb];

    for (&ea, &ca) in shell_a.exponents.iter().zip(&shell_a.coefficients) {
        for (&eb, &cb) in shell_b.exponents.iter().zip(&shell_b.coefficients) {
            let gamma = ea + eb;
            let (sx, tx) = compute_1d_st(shell_a.center[0], shell_b.center[0], ea, eb);
            let (sy, ty) = compute_1d_st(shell_a.center[1], shell_b.center[1], ea, eb);
            let (sz, tz) = compute_1d_st(shell_a.center[2], shell_b.center[2], ea, eb);

            let pre_s = ca * cb * (std::f64::consts::PI / gamma).powf(1.5);
            
            // For S: S3D = Sx * Sy * Sz
            // For T: T3D = Tx*Sy*Sz + Sx*Ty*Sz + Sx*Sy*Tz
            
            // Compute blocks
            let mut s_prim = vec![0.0; na * nb];
            let mut t_prim = vec![0.0; na * nb];

            build_3d_block(shell_a.shell_type, shell_b.shell_type, &sx, &sy, &sz, &tx, &ty, &tz, pre_s, &mut s_prim, &mut t_prim);

            for k in 0..(na * nb) {
                s_res[k] += s_prim[k];
                t_res[k] += t_prim[k];
            }
        }
    }
    (s_res, t_res)
}

fn compute_1d_st(a: f64, b: f64, ea: f64, eb: f64) -> ([[f64; 4]; 2], [[f64; 2]; 2]) {
    let gamma = ea + eb;
    let p = (ea * a + eb * b) / gamma;
    
    let mut s = [[0.0; 4]; 2]; // s[a][b] for a in 0..=1, b in 0..=3
    
    // Scale factor out: we pulled out (pi/gamma)^1.5 already
    s[0][0] = (-ea * eb / gamma * (a - b).powi(2)).exp();
    s[1][0] = (p - a) * s[0][0];
    s[0][1] = (p - b) * s[0][0];
    s[1][1] = (p - a) * s[0][1] + 1.0 / (2.0 * gamma) * s[0][0];
    s[0][2] = (p - b) * s[0][1] + 1.0 / (2.0 * gamma) * s[0][0];
    s[1][2] = (p - a) * s[0][2] + 2.0 / (2.0 * gamma) * s[0][1];
    s[0][3] = (p - b) * s[0][2] + 2.0 / (2.0 * gamma) * s[0][1];
    s[1][3] = (p - a) * s[0][3] + 3.0 / (2.0 * gamma) * s[0][2];

    let mut t = [[0.0; 2]; 2]; // t[a][b] for a,b in 0..=1
    // T = b*(2*l_b+1) * S_ab - 2*b^2 * S_a,b+2
    t[0][0] = eb * 1.0 * s[0][0] - 2.0 * eb * eb * s[0][2];
    t[1][0] = eb * 1.0 * s[1][0] - 2.0 * eb * eb * s[1][2];
    t[0][1] = eb * 3.0 * s[0][1] - 2.0 * eb * eb * s[0][3];
    t[1][1] = eb * 3.0 * s[1][1] - 2.0 * eb * eb * s[1][3];

    (s, t)
}

fn build_3d_block(
    ta: ShellType, tb: ShellType,
    sx: &[[f64; 4]; 2], sy: &[[f64; 4]; 2], sz: &[[f64; 4]; 2],
    tx: &[[f64; 2]; 2], ty: &[[f64; 2]; 2], tz: &[[f64; 2]; 2],
    pre: f64, s_out: &mut [f64], t_out: &mut [f64]
) {
    let nb = if tb == ShellType::P { 3 } else { 1 };
    
    // Helper closures
    let get_s = |ax: usize, ay: usize, az: usize, bx: usize, by: usize, bz: usize| -> f64 { sx[ax][bx] * sy[ay][by] * sz[az][bz] };
    let get_t = |ax: usize, ay: usize, az: usize, bx: usize, by: usize, bz: usize| -> f64 {
        tx[ax][bx] * sy[ay][by] * sz[az][bz] + sx[ax][bx] * ty[ay][by] * sz[az][bz] + sx[ax][bx] * sy[ay][by] * tz[az][bz]
    };

    match (ta, tb) {
        (ShellType::S, ShellType::S) => {
            s_out[0] = pre * get_s(0,0,0, 0,0,0);
            t_out[0] = pre * get_t(0,0,0, 0,0,0);
        }
        (ShellType::S, ShellType::P) => {
            s_out[0] = pre * get_s(0,0,0, 1,0,0); t_out[0] = pre * get_t(0,0,0, 1,0,0);
            s_out[1] = pre * get_s(0,0,0, 0,1,0); t_out[1] = pre * get_t(0,0,0, 0,1,0);
            s_out[2] = pre * get_s(0,0,0, 0,0,1); t_out[2] = pre * get_t(0,0,0, 0,0,1);
        }
        (ShellType::P, ShellType::S) => {
            s_out[0*nb] = pre * get_s(1,0,0, 0,0,0); t_out[0*nb] = pre * get_t(1,0,0, 0,0,0);
            s_out[1*nb] = pre * get_s(0,1,0, 0,0,0); t_out[1*nb] = pre * get_t(0,1,0, 0,0,0);
            s_out[2*nb] = pre * get_s(0,0,1, 0,0,0); t_out[2*nb] = pre * get_t(0,0,1, 0,0,0);
        }
        (ShellType::P, ShellType::P) => {
            let p_axes = [(1,0,0), (0,1,0), (0,0,1)];
            for i in 0..3 {
                for j in 0..3 {
                    let ax = p_axes[i];
                    let bx = p_axes[j];
                    let idx = i * nb + j;
                    s_out[idx] = pre * get_s(ax.0, ax.1, ax.2, bx.0, bx.1, bx.2);
                    t_out[idx] = pre * get_t(ax.0, ax.1, ax.2, bx.0, bx.1, bx.2);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::basis::build_sto3g_basis;

    #[test]
    fn test_overlap_diagonal_one() {
        let basis = build_sto3g_basis(&[1], &[[0.0, 0.0, 0.0]]);
        let s = compute_overlap_matrix(&basis);
        assert!((s[(0, 0)] - 1.0).abs() < 1e-10, "S_00 = {}", s[(0, 0)]);
    }

    #[test]
    fn test_overlap_symmetric() {
        let basis = build_sto3g_basis(&[8, 1, 1], &[
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.8],
            [0.0, 1.8, 0.0],
        ]);
        let s = compute_overlap_matrix(&basis);
        for i in 0..s.nrows() {
            for j in 0..s.nrows() {
                assert!((s[(i, j)] - s[(j, i)]).abs() < 1e-10);
            }
        }
    }
}
