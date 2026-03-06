use rayon::prelude::*;

pub trait ForceFieldContribution: Send + Sync {
    fn evaluate_energy_and_inject_gradient(&self, coords: &[f64], grad: &mut [f64]) -> f64;
}

pub struct MolecularForceField {
    pub iter_terms: Vec<Box<dyn ForceFieldContribution>>,
}

impl MolecularForceField {
    pub fn new() -> Self {
        MolecularForceField { iter_terms: Vec::new() }
    }

    pub fn insert_dynamic_term(&mut self, term: Box<dyn ForceFieldContribution>) {
        self.iter_terms.push(term);
    }

    pub fn compute_system_energy_and_gradients(&self, coords: &[f64], grad: &mut [f64]) -> f64 {
        grad.fill(0.0);
        let mut local_grads = vec![0.0; grad.len()];
        let total_energy: f64 = self.iter_terms.iter().map(|term| {
            term.evaluate_energy_and_inject_gradient(coords, &mut local_grads)
        }).sum();
        
        grad.copy_from_slice(&local_grads);
        total_energy
    }
}
