use super::traits::ForceFieldContribution;

pub struct UffHarmonicBondStretch {
    pub atom_i_idx: usize,
    pub atom_j_idx: usize,
    pub force_constant_kb: f64,
    pub equilibrium_r0: f64,
}

impl ForceFieldContribution for UffHarmonicBondStretch {
    fn evaluate_energy_and_inject_gradient(&self, coords: &[f64], grad: &mut [f64]) -> f64 {
        let root_i = self.atom_i_idx * 3;
        let root_j = self.atom_j_idx * 3;
        
        let mut diff_x = coords[root_i] - coords[root_j];
        let mut diff_y = coords[root_i+1] - coords[root_j+1];
        let mut diff_z = coords[root_i+2] - coords[root_j+2];
        
        let mut inter_r = (diff_x*diff_x + diff_y*diff_y + diff_z*diff_z).sqrt();
        
        if inter_r < 1e-10 {
            inter_r = 1e-10;
            diff_x = 1e-10;
        }
        
        let spatial_delta = inter_r - self.equilibrium_r0;
        let bond_energy = 0.5 * self.force_constant_kb * spatial_delta * spatial_delta;
        
        let vectorial_scalar_prefactor = self.force_constant_kb * spatial_delta / inter_r;
        let force_x = vectorial_scalar_prefactor * diff_x;
        let force_y = vectorial_scalar_prefactor * diff_y;
        let force_z = vectorial_scalar_prefactor * diff_z;
        
        grad[root_i] += force_x;
        grad[root_i+1] += force_y;
        grad[root_i+2] += force_z;
        
        grad[root_j] -= force_x;
        grad[root_j+1] -= force_y;
        grad[root_j+2] -= force_z;
        
        bond_energy
    }
}
