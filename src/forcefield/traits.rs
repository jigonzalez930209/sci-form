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

        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;

            let (total_energy, total_grad) = self
                .iter_terms
                .par_iter()
                .map(|term| {
                    let mut local_grad = vec![0.0; grad.len()];
                    let energy = term.evaluate_energy_and_inject_gradient(coords, &mut local_grad);
                    (energy, local_grad)
                })
                .reduce(
                    || (0.0, vec![0.0; grad.len()]),
                    |mut left, right| {
                        left.0 += right.0;
                        for (dst, src) in left.1.iter_mut().zip(right.1.into_iter()) {
                            *dst += src;
                        }
                        left
                    },
                );

            grad.copy_from_slice(&total_grad);
            total_energy
        }

        #[cfg(not(feature = "parallel"))]
        {
            let mut total_energy = 0.0;
            for term in &self.iter_terms {
                total_energy += term.evaluate_energy_and_inject_gradient(coords, grad);
            }

            total_energy
        }
    }
}
