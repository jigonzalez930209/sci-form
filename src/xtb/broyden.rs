//! Modified Broyden mixing for SCC convergence acceleration.
//!
//! Direct port of tblite's `tblite_scf_mixer_broyden` module.
//! Uses the modified Broyden method with memory-limited history
//! and LU-factorized beta matrix.

/// Broyden mixer state for SCC iteration convergence.
pub struct BroydenMixer {
    ndim: usize,
    memory: usize,
    iter: usize,
    damp: f64,
    df: Vec<Vec<f64>>,  // df[memory][ndim] — normalized delta-F vectors
    u: Vec<Vec<f64>>,   // u[memory][ndim] — Broyden update vectors
    a: Vec<Vec<f64>>,   // a[memory][memory] — dot product matrix
    omega: Vec<f64>,    // omega[memory] — iteration weights
    dq: Vec<f64>,       // dq[ndim] — current residual (output - input)
    dqlast: Vec<f64>,   // dqlast[ndim] — previous residual
    qlast_in: Vec<f64>, // qlast_in[ndim] — previous input
    q_in: Vec<f64>,     // q_in[ndim] — current input
}

impl BroydenMixer {
    /// Create a new Broyden mixer.
    ///
    /// - `ndim`: number of variables to mix (shell charges + dipoles + quadrupoles)
    /// - `memory`: maximum number of iterations to keep in history
    /// - `damp`: damping parameter (alpha), typically 0.4
    pub fn new(ndim: usize, memory: usize, damp: f64) -> Self {
        Self {
            ndim,
            memory,
            iter: 0,
            damp,
            df: vec![vec![0.0; ndim]; memory],
            u: vec![vec![0.0; ndim]; memory],
            a: vec![vec![0.0; memory]; memory],
            omega: vec![0.0; memory],
            dq: vec![0.0; ndim],
            dqlast: vec![0.0; ndim],
            qlast_in: vec![0.0; ndim],
            q_in: vec![0.0; ndim],
        }
    }

    /// Store the current input charges before SCF step.
    pub fn set(&mut self, qvec: &[f64]) {
        self.q_in[..qvec.len()].copy_from_slice(qvec);
    }

    /// Store the current input charges from multiple slices (shell charges, dipoles, quadrupoles).
    pub fn set_parts(&mut self, parts: &[&[f64]]) {
        let mut offset = 0;
        for part in parts {
            self.q_in[offset..offset + part.len()].copy_from_slice(part);
            offset += part.len();
        }
    }

    /// Compute the difference dq = new_q - q_in from a flat vector.
    pub fn diff(&mut self, qvec: &[f64]) {
        for i in 0..qvec.len() {
            self.dq[i] = qvec[i] - self.q_in[i];
        }
    }

    /// Compute the difference from multiple slices.
    pub fn diff_parts(&mut self, parts: &[&[f64]]) {
        let mut offset = 0;
        for part in parts {
            for i in 0..part.len() {
                self.dq[offset + i] = part[i] - self.q_in[offset + i];
            }
            offset += part.len();
        }
    }

    /// Apply the Broyden mixing step. Call this at the start of each iteration (after iter 0).
    /// Returns Ok(()) on success, Err on linear algebra failure.
    pub fn step(&mut self) -> Result<(), &'static str> {
        self.iter += 1;
        self.broyden_step()
    }

    /// Copy the mixed charges into a flat output vector.
    pub fn get(&self, qvec: &mut [f64]) {
        qvec.copy_from_slice(&self.q_in[..qvec.len()]);
    }

    /// Copy the mixed charges from q_in into multiple output slices.
    pub fn get_parts(&self, parts: &mut [&mut [f64]]) {
        let mut offset = 0;
        for part in parts.iter_mut() {
            part.copy_from_slice(&self.q_in[offset..offset + part.len()]);
            offset += part.len();
        }
    }

    /// Get the RMS error metric (same as tblite's get_error).
    pub fn get_error(&self) -> f64 {
        let n = self.dq.len();
        let sum_sq: f64 = self.dq.iter().map(|x| x * x).sum();
        (sum_sq / n as f64).sqrt()
    }

    /// Core Broyden mixing algorithm — direct port of tblite's `broyden` subroutine.
    fn broyden_step(&mut self) -> Result<(), &'static str> {
        let n = self.ndim;
        let itn = self.iter - 1; // 0-based iteration number (iter was already incremented)
        let it1 = itn.wrapping_sub(1) % self.memory; // mod(itn-1, memory) in 0-based

        // Parameters (matching tblite exactly)
        let alpha = self.damp;
        let omega0: f64 = 0.01;
        let minw: f64 = 1.0;
        let maxw: f64 = 100000.0;
        let wfac: f64 = 0.01;

        // First iteration: simple damping
        if self.iter == 1 {
            self.dqlast.copy_from_slice(&self.dq);
            self.qlast_in.copy_from_slice(&self.q_in);
            for i in 0..n {
                self.q_in[i] += alpha * self.dq[i];
            }
            return Ok(());
        }

        let nhistory = (itn).min(self.memory); // number of stored vectors

        // Create omega (weight) for the current iteration
        let norm_dq: f64 = self.dq.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm_dq > (wfac / maxw) {
            self.omega[it1] = wfac / norm_dq;
        } else {
            self.omega[it1] = maxw;
        }
        if self.omega[it1] < minw {
            self.omega[it1] = minw;
        }

        // Build dF(iter-1) = normalized (dq - dqlast)
        for i in 0..n {
            self.df[it1][i] = self.dq[i] - self.dqlast[i];
        }
        let norm_df: f64 = self.df[it1].iter().map(|x| x * x).sum::<f64>().sqrt();
        let inv = 1.0 / norm_df.max(f64::EPSILON);
        for i in 0..n {
            self.df[it1][i] *= inv;
        }

        // Build a, beta, c
        // a[i][it1] = df[i] . df[it1]
        let j_start = if itn > self.memory {
            itn - self.memory + 1
        } else {
            1
        };
        for j in j_start..=itn {
            let i = (j - 1) % self.memory;
            let dot: f64 = (0..n).map(|k| self.df[i][k] * self.df[it1][k]).sum();
            self.a[i][it1] = dot;
            self.a[it1][i] = dot;
        }

        // c[i] = omega[i] * (df[i] . dq)
        let mut c = vec![0.0f64; nhistory];
        for j in j_start..=itn {
            let i = (j - 1) % self.memory;
            let idx = if itn <= self.memory {
                j - 1
            } else {
                j - j_start
            };
            let dot: f64 = (0..n).map(|k| self.df[i][k] * self.dq[k]).sum();
            c[idx] = self.omega[i] * dot;
        }

        // Build beta = omega_i * omega_j * a[i][j] + omega0² * delta_ij
        let mut beta = vec![vec![0.0f64; nhistory]; nhistory];
        for jj in 0..nhistory {
            let j = j_start + jj;
            let jmod = (j - 1) % self.memory;
            for ii in 0..nhistory {
                let i = j_start + ii;
                let imod = (i - 1) % self.memory;
                beta[ii][jj] = self.omega[imod] * self.omega[jmod] * self.a[imod][jmod];
            }
            beta[jj][jj] += omega0 * omega0;
        }

        // Solve beta * gamma = c (in-place, c becomes gamma)
        Self::solve_linear(&mut beta, &mut c)?;

        // Build u[it1] = alpha * df[it1] + inv * (q_in - qlast_in)
        for i in 0..n {
            self.u[it1][i] = alpha * self.df[it1][i] + inv * (self.q_in[i] - self.qlast_in[i]);
        }

        // Save charges and deltas
        self.dqlast.copy_from_slice(&self.dq);
        self.qlast_in.copy_from_slice(&self.q_in);

        // Calculate new charges: q = q + alpha * dq
        for i in 0..n {
            self.q_in[i] += alpha * self.dq[i];
        }

        // Subtract Broyden correction: q = q - Σ omega[i] * gamma[i] * u[i]
        for jj in 0..nhistory {
            let j = j_start + jj;
            let jmod = (j - 1) % self.memory;
            let factor = self.omega[jmod] * c[jj];
            for i in 0..n {
                self.q_in[i] -= factor * self.u[jmod][i];
            }
        }

        Ok(())
    }

    /// Solve Ax = b in-place using LU decomposition with partial pivoting.
    /// On exit, `b` contains the solution vector.
    fn solve_linear(a: &mut [Vec<f64>], b: &mut [f64]) -> Result<(), &'static str> {
        let n = b.len();
        if n == 0 {
            return Ok(());
        }

        let mut ipiv = vec![0usize; n];

        // LU factorization with partial pivoting
        for k in 0..n {
            // Find pivot
            let mut max_val = a[k][k].abs();
            let mut max_idx = k;
            for i in (k + 1)..n {
                if a[i][k].abs() > max_val {
                    max_val = a[i][k].abs();
                    max_idx = i;
                }
            }
            ipiv[k] = max_idx;

            if max_val < 1e-30 {
                return Err("Singular matrix in Broyden mixer");
            }

            // Swap rows
            if max_idx != k {
                for j in 0..n {
                    let tmp = a[k][j];
                    a[k][j] = a[max_idx][j];
                    a[max_idx][j] = tmp;
                }
            }

            // Eliminate below
            for i in (k + 1)..n {
                a[i][k] /= a[k][k];
                for j in (k + 1)..n {
                    a[i][j] -= a[i][k] * a[k][j];
                }
            }
        }

        // Apply pivots to b
        for k in 0..n {
            if ipiv[k] != k {
                b.swap(k, ipiv[k]);
            }
        }

        // Forward substitution (L)
        for i in 1..n {
            for j in 0..i {
                b[i] -= a[i][j] * b[j];
            }
        }

        // Back substitution (U)
        for i in (0..n).rev() {
            for j in (i + 1)..n {
                b[i] -= a[i][j] * b[j];
            }
            b[i] /= a[i][i];
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_broyden_linear_damping() {
        // First iteration should do simple linear damping
        let mut mixer = BroydenMixer::new(3, 10, 0.4);

        // Set input = [1.0, 2.0, 3.0]
        mixer.set(&[1.0, 2.0, 3.0]);

        // New output = [1.5, 2.5, 3.5], so diff = [0.5, 0.5, 0.5]
        mixer.diff(&[1.5, 2.5, 3.5]);

        // Apply mixing
        mixer.step().unwrap();

        // Should get q_in + alpha * dq = [1.0, 2.0, 3.0] + 0.4 * [0.5, 0.5, 0.5]
        let mut out = vec![0.0; 3];
        mixer.get(&mut out);
        assert!((out[0] - 1.2).abs() < 1e-10);
        assert!((out[1] - 2.2).abs() < 1e-10);
        assert!((out[2] - 3.2).abs() < 1e-10);
    }

    #[test]
    fn test_solve_linear_simple() {
        // 2x + 1y = 5
        // 1x + 3y = 7
        // Solution: x = 1.6, y = 1.8
        let mut a = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        let mut b = vec![5.0, 7.0];
        BroydenMixer::solve_linear(&mut a, &mut b).unwrap();
        assert!((b[0] - 1.6).abs() < 1e-10);
        assert!((b[1] - 1.8).abs() < 1e-10);
    }
}
