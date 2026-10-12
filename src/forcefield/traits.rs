/// Core interface for a molecular force-field contribution.
pub trait ForceFieldContribution: Send + Sync {
    /// Evaluate the energy term and accumulate its Cartesian gradient contribution.
    fn evaluate_energy_and_inject_gradient(&self, coords: &[f64], grad: &mut [f64]) -> f64;
    /// Short term-type label for diagnostics (validation F70).
    fn kind(&self) -> &'static str {
        "unknown"
    }
}

pub struct MolecularForceField {
    pub iter_terms: Vec<Box<dyn ForceFieldContribution>>,
}

impl Default for MolecularForceField {
    fn default() -> Self {
        Self::new()
    }
}

/// Sum the energy of a term list per term type, preserving first-seen order.
pub fn aggregate_breakdown(
    terms: &[Box<dyn ForceFieldContribution>],
    coords: &[f64],
) -> Vec<(&'static str, f64)> {
    let mut scratch = vec![0.0f64; coords.len()];
    let mut kinds: Vec<&'static str> = Vec::new();
    let mut sums: Vec<f64> = Vec::new();
    for term in terms {
        let energy = term.evaluate_energy_and_inject_gradient(coords, &mut scratch);
        match kinds.iter().position(|&k| k == term.kind()) {
            Some(idx) => sums[idx] += energy,
            None => {
                kinds.push(term.kind());
                sums.push(energy);
            }
        }
    }
    kinds.into_iter().zip(sums).collect()
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

    /// Per-term energy breakdown at `coords` (kcal/mol), in registration order.
    ///
    /// Each term is evaluated independently (its gradient contribution is
    /// discarded); use [`MolecularForceField::energy_breakdown_aggregated`]
    /// for one entry per term type.
    pub fn energy_breakdown(&self, coords: &[f64]) -> Vec<(&'static str, f64)> {
        let mut scratch = vec![0.0f64; coords.len()];
        self.iter_terms
            .iter()
            .map(|term| {
                let energy = term.evaluate_energy_and_inject_gradient(coords, &mut scratch);
                (term.kind(), energy)
            })
            .collect()
    }

    /// Energy summed per term type (bond/angle/torsion/inversion/vdw/...).
    pub fn energy_breakdown_aggregated(&self, coords: &[f64]) -> Vec<(&'static str, f64)> {
        aggregate_breakdown(&self.iter_terms, coords)
    }

    /// Evaluate all registered terms for the current geometry.
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
