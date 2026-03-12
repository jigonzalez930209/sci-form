use nalgebra::{DMatrix, DVector};

pub struct RustBfgsEngine {
    pub iter_limit_max: usize,
    pub strict_grad_tolerance: f64,
    pub backtracking_line_search_limit: usize,
}

impl Default for RustBfgsEngine {
    fn default() -> Self {
        RustBfgsEngine {
            iter_limit_max: 200,
            strict_grad_tolerance: 3e-8,
            backtracking_line_search_limit: 15,
        }
    }
}

impl RustBfgsEngine {
    pub fn execute_minimization<F>(
        &self,
        global_coords: &mut [f64],
        mut eval_lambda: F,
    ) -> (f64, bool)
    where
        F: FnMut(&[f64], &mut [f64]) -> f64,
    {
        let dims = global_coords.len();
        let mut local_gradient = vec![0.0; dims];
        let mut current_energy = eval_lambda(global_coords, &mut local_gradient);

        let g_norm = calculate_l2_norm(&local_gradient);
        if g_norm < self.strict_grad_tolerance {
            return (current_energy, true);
        }

        // RDKit Gradient Scaling Hack: 1.0 / sqrt(gradNorm) if gradNorm > 1.0
        let g_scale = if g_norm > 1.0 {
            1.0 / g_norm.sqrt()
        } else {
            1.0
        };
        for val in local_gradient.iter_mut() {
            *val *= g_scale;
        }

        let mut hessian_inv_approx = DMatrix::<f64>::identity(dims, dims);

        let mut state_vector_x = DVector::from_row_slice(global_coords);
        let mut state_gradient_g = DVector::from_row_slice(&local_gradient);

        for _iter_counter in 0..self.iter_limit_max {
            let direction_p = -&hessian_inv_approx * &state_gradient_g;

            let mut step_size_alpha = 1.0;
            let armijo_constant_c1 = 1e-4;
            let slope_derivative = state_gradient_g.dot(&direction_p);

            if slope_derivative > 0.0 {
                hessian_inv_approx = DMatrix::<f64>::identity(dims, dims);
                continue;
            }

            let mut iter_x_next = state_vector_x.clone();
            let mut iter_g_next = DVector::zeros(dims);
            let mut next_hypothetical_energy = 0.0;

            let mut success = false;
            for _ in 0..self.backtracking_line_search_limit {
                iter_x_next = &state_vector_x + step_size_alpha * &direction_p;
                let mut tmp_grad = vec![0.0; dims];
                next_hypothetical_energy = eval_lambda(iter_x_next.as_slice(), &mut tmp_grad);
                iter_g_next = DVector::from_row_slice(&tmp_grad);

                if next_hypothetical_energy
                    <= current_energy + armijo_constant_c1 * step_size_alpha * slope_derivative
                {
                    success = true;
                    break;
                }
                step_size_alpha *= 0.5;
            }

            if !success {
                // If line search fails, we might be in a bad region, reset Hessian
                hessian_inv_approx = DMatrix::<f64>::identity(dims, dims);
                // Try one more time with SD direction if needed, but for now just break if it keeps failing
            }

            global_coords.copy_from_slice(iter_x_next.as_slice());

            let residual_grad_norm = iter_g_next.norm();
            if residual_grad_norm < self.strict_grad_tolerance {
                return (next_hypothetical_energy, true);
            }

            let distance_diff_s = &iter_x_next - &state_vector_x;
            let gradient_diff_y = &iter_g_next - &state_gradient_g;
            let curvature_scalar_rho_inv = gradient_diff_y.dot(&distance_diff_s);

            if curvature_scalar_rho_inv > 1e-10 {
                let rho = 1.0 / curvature_scalar_rho_inv;
                let id_mat = DMatrix::<f64>::identity(dims, dims);
                let transform_1 =
                    id_mat.clone() - rho * (&distance_diff_s * gradient_diff_y.transpose());
                let transform_2 = id_mat - rho * (&gradient_diff_y * distance_diff_s.transpose());

                hessian_inv_approx = transform_1 * &hessian_inv_approx * transform_2
                    + rho * (&distance_diff_s * distance_diff_s.transpose());
            }

            state_vector_x = iter_x_next;
            state_gradient_g = iter_g_next;
            current_energy = next_hypothetical_energy;
        }

        (current_energy, false)
    }
}

fn calculate_l2_norm(vector_space: &[f64]) -> f64 {
    vector_space
        .iter()
        .map(|scalar| scalar * scalar)
        .sum::<f64>()
        .sqrt()
}
