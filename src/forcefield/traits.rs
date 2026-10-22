//! Core abstractions shared by every force field in the crate.
//!
//! A force field is a list of [`ForceFieldContribution`] terms held in a
//! [`MolecularForceField`]. Each term evaluates its own energy and adds its
//! own Cartesian gradient, so new terms can be mixed freely (UFF, MMFF94 and
//! the distance-geometry restraints all implement the same trait).
//!
//! Energies are in kcal/mol and gradients are `dE/dR` in kcal/mol/Å.

/// Core interface for a molecular force-field contribution.
pub trait ForceFieldContribution: Send + Sync {
    /// Evaluate the energy term and accumulate its Cartesian gradient contribution.
    fn evaluate_energy_and_inject_gradient(&self, coords: &[f64], grad: &mut [f64]) -> f64;
    /// Short term-type label for diagnostics (validation F70).
    fn kind(&self) -> &'static str {
        "unknown"
    }
}

/// Number of force-field terms accumulated into one deterministic chunk.
///
/// The chunk size is a compile-time constant: it does **not** depend on the
/// number of rayon worker threads.  Terms inside a chunk are summed
/// sequentially in registration order and the chunk partial sums are then
/// combined in chunk order, so the bit pattern of the energy and of the
/// gradient is identical for every thread count.
///
/// The chunking also bounds the number of simultaneously live per-chunk
/// gradient buffers to `ceil(n_terms / DETERMINISTIC_TERM_CHUNK)`, which keeps
/// the memory traffic of the parallel path proportional to the number of
/// chunks instead of to the number of terms (blocker D-H4).
pub const DETERMINISTIC_TERM_CHUNK: usize = 64;

/// A collection of force-field terms that together define the energy of a system.
///
/// The terms are evaluated in registration order and their energies summed;
/// the total is in kcal/mol for the UFF and MMFF94 terms. The gradient
/// convention is `dE/dR` (kcal/mol/Å), so the mechanical force is the
/// negative of the returned gradient.
pub struct MolecularForceField {
    /// The registered terms, evaluated in this order.
    pub iter_terms: Vec<Box<dyn ForceFieldContribution>>,
}

impl Default for MolecularForceField {
    fn default() -> Self {
        Self::new()
    }
}

/// Sum the energy of a term list per term type, preserving first-seen order.
///
/// Returns one `(term_type, energy)` pair per distinct term type, in the order
/// the types were first encountered. The energies are in kcal/mol and sum to
/// the total energy of the term list.
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
    /// Create an empty force field with no registered terms.
    pub fn new() -> Self {
        MolecularForceField {
            iter_terms: Vec::new(),
        }
    }

    /// Register an additional term; terms are evaluated in registration order.
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
    ///
    /// # Determinism
    ///
    /// The accumulation order is fixed and **independent of the number of
    /// rayon worker threads**:
    ///
    /// 1. the term list is split into fixed-size chunks of
    ///    [`DETERMINISTIC_TERM_CHUNK`] terms;
    /// 2. every chunk accumulates its terms sequentially, in registration
    ///    order, into a private gradient buffer;
    /// 3. the chunk results are combined in chunk order.
    ///
    /// The previous implementation used a rayon `reduce`, whose tree shape —
    /// and therefore whose floating-point rounding — depended on
    /// `rayon::current_num_threads()`.  Two runs of the same geometry with a
    /// different thread count produced different bits (blocker D-C1).
    pub fn compute_system_energy_and_gradients(&self, coords: &[f64], grad: &mut [f64]) -> f64 {
        grad.fill(0.0);

        let n_terms = self.iter_terms.len();
        if n_terms == 0 {
            return 0.0;
        }

        let grad_len = grad.len();
        let n_chunks = n_terms.div_ceil(DETERMINISTIC_TERM_CHUNK);

        // Evaluate one chunk: terms are summed sequentially in registration
        // order into a private gradient buffer.
        let evaluate_chunk = |chunk_index: usize| -> (f64, Vec<f64>) {
            let start = chunk_index * DETERMINISTIC_TERM_CHUNK;
            let end = (start + DETERMINISTIC_TERM_CHUNK).min(n_terms);
            let mut chunk_energy = 0.0f64;
            let mut chunk_grad = vec![0.0f64; grad_len];
            for term in &self.iter_terms[start..end] {
                chunk_energy += term.evaluate_energy_and_inject_gradient(coords, &mut chunk_grad);
            }
            (chunk_energy, chunk_grad)
        };

        #[cfg(feature = "parallel")]
        let chunk_results: Vec<(f64, Vec<f64>)> = {
            use rayon::prelude::*;
            // `collect` on an *indexed* parallel iterator preserves the input
            // order, so the chunk list is identical for every thread count.
            (0..n_chunks).into_par_iter().map(evaluate_chunk).collect()
        };

        #[cfg(not(feature = "parallel"))]
        let chunk_results: Vec<(f64, Vec<f64>)> = (0..n_chunks).map(evaluate_chunk).collect();

        // Combine the chunk results in chunk order (fixed, thread-count
        // independent).
        let mut total_energy = 0.0f64;
        let mut total_grad = vec![0.0f64; grad_len];
        for (chunk_energy, chunk_grad) in chunk_results {
            total_energy += chunk_energy;
            for (dst, src) in total_grad.iter_mut().zip(chunk_grad.iter()) {
                *dst += src;
            }
        }

        grad.copy_from_slice(&total_grad);
        total_energy
    }
}

#[cfg(test)]
mod tests {

    /// Environment marker used to run this test binary as a *child* process
    /// with a different `RAYON_NUM_THREADS` value.
    const CHILD_ENV: &str = "SCIFORM_DETERMINISM_CHILD";

    /// FNV-1a 64-bit hash, used to fingerprint a `f64` vector bit-exactly.
    fn fnv1a64(bytes: &[u8]) -> u64 {
        let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
        for byte in bytes {
            hash ^= u64::from(*byte);
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
        hash
    }

    /// Bit-exact fingerprint of a `f64` slice.
    fn hash_f64(values: &[f64]) -> u64 {
        let mut bytes = Vec::with_capacity(values.len() * 8);
        for value in values {
            bytes.extend_from_slice(&value.to_bits().to_le_bytes());
        }
        fnv1a64(&bytes)
    }

    /// A fixed butane-like geometry (4 C + 10 H = 14 atoms).
    ///
    /// The coordinates are hard-coded on purpose: the test must not depend on
    /// the embedding code path, which itself consumes the force-field energy
    /// and would therefore inherit the very non-determinism under test.
    fn determinism_test_system() -> (Vec<u8>, Vec<f64>) {
        let mut elements: Vec<u8> = Vec::new();
        let mut coords: Vec<f64> = Vec::new();

        let carbons: [[f64; 3]; 4] = [
            [0.000, 0.000, 0.000],
            [1.270, 0.000, 0.750],
            [2.540, 0.000, 0.000],
            [3.810, 0.000, 0.750],
        ];
        let hydrogens: [[f64; 3]; 10] = [
            [-0.363, 0.935, -0.545],
            [-0.363, -0.935, -0.545],
            [-0.726, 0.000, 0.682],
            [1.633, 0.935, 0.205],
            [1.633, -0.935, 0.205],
            [2.177, 0.935, -0.545],
            [2.177, -0.935, -0.545],
            [4.173, 0.935, 0.205],
            [4.173, -0.935, 0.205],
            [3.447, 0.000, 1.432],
        ];

        for c in carbons.iter() {
            elements.push(6);
            coords.extend_from_slice(c);
        }
        for h in hydrogens.iter() {
            elements.push(1);
            coords.extend_from_slice(h);
        }

        (elements, coords)
    }

    /// UFF energy + gradient of the fixed test system.
    fn uff_energy_and_gradient() -> (f64, Vec<f64>) {
        let (_elements, coords) = determinism_test_system();
        let mol = crate::graph::Molecule::from_smiles("CCCC").expect("butane parses");
        let ff = crate::forcefield::builder::build_uff_force_field(&mol);
        let mut grad = vec![0.0f64; coords.len()];
        let energy = ff.compute_system_energy_and_gradients(&coords, &mut grad);
        (energy, grad)
    }

    /// Generalized-Born electrostatic energy of a larger fixed test system.
    ///
    /// 40 atoms are used on purpose: the GB energy is the sum of one term
    /// per atom, so a small system does not exercise the reduction order
    /// enough to expose a thread-count dependent rounding.
    fn gb_electrostatic_energy() -> f64 {
        const N: usize = 40;
        let elements: Vec<u8> = (0..N).map(|i| if i % 3 == 0 { 6u8 } else { 1u8 }).collect();
        let positions: Vec<[f64; 3]> = (0..N)
            .map(|i| [i as f64 * 1.4, (i % 2) as f64 * 0.9, (i % 3) as f64 * 0.6])
            .collect();
        // Deterministic, non-trivial charges (alternating sign by atom index).
        let charges: Vec<f64> = (0..N)
            .map(|i| if i % 2 == 0 { -0.12 } else { 0.08 })
            .collect();
        crate::solvation::compute_gb_solvation(&elements, &positions, &charges, None, None, None)
            .electrostatic_energy_kcal_mol
    }

    /// Child half of the subprocess determinism test.
    ///
    /// It is a no-op unless the parent re-executed the test binary with
    /// `SCIFORM_DETERMINISM_CHILD=1`.
    #[test]
    fn determinism_child_reports_bits() {
        if std::env::var(CHILD_ENV).is_err() {
            return;
        }

        let (energy, grad) = uff_energy_and_gradient();
        let gb = gb_electrostatic_energy();

        println!(
            "DETERMINISM energy_bits={:016x} grad_hash={:016x} gb_bits={:016x}",
            energy.to_bits(),
            hash_f64(&grad),
            gb.to_bits()
        );
    }

    /// D-C1 regression test: the force-field energy/gradient and the
    /// Generalized-Born electrostatic energy must be **bit-identical** no
    /// matter how many rayon worker threads are used.
    ///
    /// `std::env::set_var` is not enough here — rayon reads
    /// `RAYON_NUM_THREADS` when the global pool is created — so the test
    /// re-executes its own test binary as a child process with a different
    /// thread count and compares the raw bits.
    #[test]
    fn compute_system_energy_and_gradients_is_thread_count_independent() {
        if std::env::var(CHILD_ENV).is_ok() {
            // We are the child process: the parent does the comparison.
            return;
        }

        let exe = std::env::current_exe().expect("test binary path is available");
        let child_test = "forcefield::traits::tests::determinism_child_reports_bits";

        let run_child = |threads: &str| -> String {
            let output = std::process::Command::new(&exe)
                .env(CHILD_ENV, "1")
                .env("RAYON_NUM_THREADS", threads)
                .arg("--exact")
                .arg(child_test)
                .arg("--nocapture")
                .output()
                .expect("child test process runs");
            assert!(
                output.status.success(),
                "child process with RAYON_NUM_THREADS={threads} failed: {}",
                String::from_utf8_lossy(&output.stderr)
            );
            String::from_utf8_lossy(&output.stdout).to_string()
        };

        let one_thread = run_child("1");
        let many_threads = run_child("16");

        let extract = |text: &str, key: &str| -> String {
            text.lines()
                .find(|line| line.starts_with("DETERMINISM "))
                .and_then(|line| {
                    line.split_whitespace()
                        .find(|token| token.starts_with(key))
                        .map(|token| token.to_string())
                })
                .unwrap_or_else(|| {
                    panic!("child output did not contain '{key}': {text}");
                })
        };

        assert_eq!(
            extract(&one_thread, "energy_bits="),
            extract(&many_threads, "energy_bits="),
            "UFF energy differs between 1 and 16 rayon threads"
        );
        assert_eq!(
            extract(&one_thread, "grad_hash="),
            extract(&many_threads, "grad_hash="),
            "UFF gradient differs between 1 and 16 rayon threads"
        );
        assert_eq!(
            extract(&one_thread, "gb_bits="),
            extract(&many_threads, "gb_bits="),
            "GB electrostatic energy differs between 1 and 16 rayon threads"
        );
    }

    /// The chunked accumulation must agree with a plain sequential evaluation
    /// of the same terms (the chunk boundaries must not change the result
    /// beyond the rounding of the summation grouping).
    #[test]
    fn chunked_accumulation_matches_sequential_reference() {
        let (_elements, coords) = determinism_test_system();
        let mol = crate::graph::Molecule::from_smiles("CCCC").expect("butane parses");
        let ff = crate::forcefield::builder::build_uff_force_field(&mol);

        let mut chunked_grad = vec![0.0f64; coords.len()];
        let chunked_energy = ff.compute_system_energy_and_gradients(&coords, &mut chunked_grad);

        let mut sequential_grad = vec![0.0f64; coords.len()];
        let mut sequential_energy = 0.0f64;
        for term in &ff.iter_terms {
            sequential_energy +=
                term.evaluate_energy_and_inject_gradient(&coords, &mut sequential_grad);
        }

        assert!(
            (chunked_energy - sequential_energy).abs() <= 1e-9 * sequential_energy.abs().max(1.0),
            "chunked energy {chunked_energy} differs from sequential {sequential_energy}"
        );
        for (index, (a, b)) in chunked_grad.iter().zip(sequential_grad.iter()).enumerate() {
            assert!(
                (a - b).abs() <= 1e-9 * b.abs().max(1.0),
                "gradient component {index} differs: {a} vs {b}"
            );
        }
    }
}
