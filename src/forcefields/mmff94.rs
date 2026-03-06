use super::traits::ForceFieldContribution;

pub struct Mmff94BufferedVanDerWaals {
    pub atom_i_idx: usize,
    pub atom_j_idx: usize,
    pub radius_star: f64,
    pub epsilon_depth: f64,
}

impl ForceFieldContribution for Mmff94BufferedVanDerWaals {
    fn evaluate_energy_and_inject_gradient(&self, coords: &[f64], grad: &mut [f64]) -> f64 {
        let root_i = self.atom_i_idx * 3;
        let root_j = self.atom_j_idx * 3;
        
        let delta_x = coords[root_i] - coords[root_j];
        let delta_y = coords[root_i+1] - coords[root_j+1];
        let delta_z = coords[root_i+2] - coords[root_j+2];
        
        let dist_squared = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z;
        let mut dist_r = dist_squared.sqrt();
        
        if dist_r < 1e-8 { dist_r = 1e-8; }
        
        let r_star_powered_7 = self.radius_star.powi(7);
        let dist_r_powered_7 = dist_r.powi(7);
        
        let repulsive_denominator = dist_r + 0.07 * self.radius_star;
        let repulsive_term = (1.07 * self.radius_star / repulsive_denominator).powi(7);
        
        let attractive_denominator = dist_r_powered_7 + 0.12 * r_star_powered_7;
        let attractive_term = (1.12 * r_star_powered_7 / attractive_denominator) - 2.0;
        
        let vdw_total_energy = self.epsilon_depth * repulsive_term * attractive_term;
        
        let gradient_rep_term = -7.0 * repulsive_term / repulsive_denominator;
        let gradient_attr_term = -7.0 * dist_r.powi(6) * (1.12 * r_star_powered_7) / (attractive_denominator * attractive_denominator);
        
        let force_scalar_magnitude = self.epsilon_depth * (gradient_rep_term * attractive_term + repulsive_term * gradient_attr_term);
        
        let vector_prefactor = force_scalar_magnitude / dist_r;
        let grad_x = vector_prefactor * delta_x;
        let grad_y = vector_prefactor * delta_y;
        let grad_z = vector_prefactor * delta_z;
        
        grad[root_i] += grad_x;
        grad[root_i+1] += grad_y;
        grad[root_i+2] += grad_z;
        
        grad[root_j] -= grad_x;
        grad[root_j+1] -= grad_y;
        grad[root_j+2] -= grad_z;
        
        vdw_total_energy
    }
}
