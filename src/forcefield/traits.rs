/// Interfaz conductiva fundamental para cualquier operador matemático empírico molecular
pub trait ForceFieldContribution: Send + Sync {
    /// Computa algebraicamente el escalar de energía asumiendo posición local
    /// e injerta los gradientes geométricos acumulativamente a lo largo del array
    fn evaluate_energy_and_inject_gradient(&self, coords: &[f64], grad: &mut [f64]) -> f64;
}

pub struct MolecularForceField {
    pub iter_terms: Vec<Box<dyn ForceFieldContribution>>,
}

impl Default for MolecularForceField {
    fn default() -> Self {
        Self::new()
    }
}

impl MolecularForceField {
    pub fn new() -> Self {
        MolecularForceField {
            iter_terms: Vec::new(),
        }
    }

    pub fn insert_dynamic_term(&mut self, term: Box<dyn ForceFieldContribution>) {
        self.iter_terms.push(term);
    }

    /// Método de evaluación crítica masiva invocado exhaustivamente iteración por iteración.
    pub fn compute_system_energy_and_gradients(&self, coords: &[f64], grad: &mut [f64]) -> f64 {
        grad.fill(0.0);

        // Parallel evaluation using rayon
        // We use a mutex or temporary local gradients to avoid data races when writing to 'grad'
        // For simplicity and matching r.md, we use a local grad vector per term if needed,
        // or just sequential for now if the number of terms is small.
        // Actually r.md suggests:
        /*
        let mut local_grads = vec![0.0; grad.len()];
        let total_energy = self.iter_terms.iter().map(|term| {
            term.evaluate_energy_and_inject_gradient(coords, &mut local_grads)
        }).sum();
        */
        // But the above has a shared local_grads which is still a problem if parallel.

        // Let's implement it sequentially first as a baseline, or use rayon's fold/reduce if parallel.
        let mut total_energy = 0.0;
        for term in &self.iter_terms {
            total_energy += term.evaluate_energy_and_inject_gradient(coords, grad);
        }

        total_energy
    }
}
