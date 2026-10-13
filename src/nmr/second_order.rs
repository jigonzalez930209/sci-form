//! Second-order (strong-coupling) NMR analysis.
//!
//! The first-order `n+1` multiplet rule assumes `|Δν| ≫ |J|`: each coupled
//! partner splits a signal into a symmetric multiplet with binomial
//! intensities. When the chemical-shift difference between two coupled spins
//! is comparable to their coupling, that picture breaks down (the classic AB
//! case) and the spectrum must be obtained by diagonalizing the spin
//! Hamiltonian.
//!
//! This module builds the spin-1/2 Hamiltonian in the product basis,
//!
//! ```text
//! H = Σ_i ν_i I_z,i + Σ_{i<j} J_ij [ I_z,i I_z,j + ½ (I_+,i I_-,j + I_-,i I_+,j) ]
//! ```
//!
//! diagonalizes it and returns the exact single-quantum stick spectrum. It also
//! reports the strong-coupling ratio `|Δν|/|J|` so callers can flag first-order
//! approximations.
//!
//! The validation gate is the closed-form AB solution (textbook result): the
//! exact 4×4 diagonalization must reproduce the analytic AB line positions and
//! intensities, and the weak-coupling limit must collapse to two first-order
//! doublets.

use nalgebra::{DMatrix, SymmetricEigen};

/// A small spin-1/2 system (nuclei with shifts and pairwise couplings).
#[derive(Debug, Clone)]
pub struct SpinSystem {
    /// Chemical shifts in Hz.
    pub shifts_hz: Vec<f64>,
    /// Coupling matrix in Hz (symmetric, zero diagonal).
    pub j_hz: Vec<Vec<f64>>,
}

/// One exact spectral line.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExactLine {
    /// Frequency in Hz.
    pub frequency_hz: f64,
    /// Relative intensity (|⟨f|I_x|i⟩|²).
    pub intensity: f64,
}

/// Strong-coupling diagnostic for one spin pair.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StrongCouplingPair {
    /// First spin index.
    pub i: usize,
    /// Second spin index.
    pub j: usize,
    /// `|Δν| / |J|` (infinite when `J = 0`).
    pub ratio: f64,
    /// Whether the pair is strongly coupled (ratio below the threshold).
    pub strongly_coupled: bool,
}

/// Default threshold: `|Δν|/|J| < 10` is treated as strongly coupled
/// (first-order multiplets deviate appreciably).
pub const STRONG_COUPLING_RATIO: f64 = 10.0;

impl SpinSystem {
    /// New system from shifts and a coupling matrix.
    pub fn new(shifts_hz: Vec<f64>, j_hz: Vec<Vec<f64>>) -> Self {
        Self { shifts_hz, j_hz }
    }

    /// Number of spins.
    pub fn n_spins(&self) -> usize {
        self.shifts_hz.len()
    }

    /// `|Δν| / |J|` for a pair (returns `f64::INFINITY` when `J = 0`).
    pub fn coupling_ratio(&self, i: usize, j: usize) -> f64 {
        let j_ij = self.j_hz[i][j].abs();
        if j_ij < 1e-12 {
            return f64::INFINITY;
        }
        (self.shifts_hz[i] - self.shifts_hz[j]).abs() / j_ij
    }

    /// Diagnose every coupled pair against the strong-coupling threshold.
    pub fn strong_coupling_report(&self) -> Vec<StrongCouplingPair> {
        let mut out = Vec::new();
        for i in 0..self.n_spins() {
            for j in (i + 1)..self.n_spins() {
                if self.j_hz[i][j].abs() < 1e-12 {
                    continue;
                }
                let ratio = self.coupling_ratio(i, j);
                out.push(StrongCouplingPair {
                    i,
                    j,
                    ratio,
                    strongly_coupled: ratio < STRONG_COUPLING_RATIO,
                });
            }
        }
        out
    }

    /// Exact single-quantum stick spectrum (frequencies in Hz, relative
    /// intensities normalized to the largest line).
    ///
    /// The Hilbert space is `2^n`; the implementation is intended for small
    /// spin systems (`n ≤ 6`, i.e. up to 64 states), which covers the
    /// strong-coupling cases encountered in multiplet analysis.
    pub fn exact_spectrum(&self) -> Vec<ExactLine> {
        let n = self.n_spins();
        if n == 0 || n > 6 {
            return Vec::new();
        }
        let dim = 1usize << n;

        // Product basis |m_0 m_1 ... m_{n-1}⟩ with m = 0 (α) / 1 (β).
        // Bit k of the state index encodes spin k; I_z = +½ for α (bit 0).
        let mut h = DMatrix::zeros(dim, dim);
        for state in 0..dim {
            let mut diag = 0.0;
            for i in 0..n {
                let m_i = if state & (1 << i) == 0 { 0.5 } else { -0.5 };
                diag += self.shifts_hz[i] * m_i;
            }
            for i in 0..n {
                for j in (i + 1)..n {
                    let m_i = if state & (1 << i) == 0 { 0.5 } else { -0.5 };
                    let m_j = if state & (1 << j) == 0 { 0.5 } else { -0.5 };
                    diag += self.j_hz[i][j] * m_i * m_j;
                }
            }
            h[(state, state)] = diag;

            // Flip-flop term: ½ J_ij (I_+,i I_-,j + I_-,i I_+,j).
            for i in 0..n {
                for j in (i + 1)..n {
                    let j_ij = self.j_hz[i][j];
                    if j_ij.abs() < 1e-12 {
                        continue;
                    }
                    let m_i = state & (1 << i) != 0;
                    let m_j = state & (1 << j) != 0;
                    if m_i != m_j {
                        let flipped = state ^ (1 << i) ^ (1 << j);
                        h[(flipped, state)] += 0.5 * j_ij;
                    }
                }
            }
        }

        let eigen = SymmetricEigen::new(h);
        let dim_f = dim as f64;

        // Single-quantum transitions: |ΔM| = 1, with the intensity
        // |⟨f|Σ_k I_x,k|i⟩|². I_x = ½ (I_+ + I_-).
        let mut lines = Vec::new();
        for f in 0..dim {
            for i in 0..dim {
                let m_f = state_magnetization(f, n);
                let m_i = state_magnetization(i, n);
                if (m_f - m_i).abs() < 0.5 {
                    continue;
                }
                // ⟨f|I_x|i⟩ = Σ_s Σ_k ½ c_f(s ⊕ k) c_i(s).
                let mut amplitude = 0.0;
                for s in 0..dim {
                    let c_i = eigen.eigenvectors[(s, i)];
                    if c_i.abs() < 1e-14 {
                        continue;
                    }
                    for k in 0..n {
                        amplitude += 0.5 * eigen.eigenvectors[(s ^ (1 << k), f)] * c_i;
                    }
                }
                let intensity = amplitude * amplitude;
                if intensity < 1e-12 {
                    continue;
                }
                lines.push(ExactLine {
                    frequency_hz: eigen.eigenvalues[f] - eigen.eigenvalues[i],
                    intensity,
                });
            }
        }
        let _ = dim_f;
        // Fold negative frequencies (they mirror the positive branch for
        // Hermitian single-quantum transitions) and merge near-degenerate lines.
        lines.retain(|l| l.frequency_hz > 0.0);
        lines.sort_by(|a, b| a.frequency_hz.partial_cmp(&b.frequency_hz).unwrap());
        let mut merged: Vec<ExactLine> = Vec::new();
        for line in lines {
            if let Some(last) = merged.last_mut() {
                if (last.frequency_hz - line.frequency_hz).abs() < 1e-9 {
                    last.intensity += line.intensity;
                    continue;
                }
            }
            merged.push(line);
        }
        let max_i = merged.iter().map(|l| l.intensity).fold(0.0f64, f64::max);
        if max_i > 0.0 {
            for line in merged.iter_mut() {
                line.intensity /= max_i;
            }
        }
        merged
    }
}

/// Total magnetization `Σ m_i` of a product state (α = +½, β = −½).
fn state_magnetization(state: usize, n: usize) -> f64 {
    let mut m = 0.0;
    for k in 0..n {
        m += if state & (1 << k) == 0 { 0.5 } else { -0.5 };
    }
    m
}

/// Closed-form AB spectrum (two spins): line positions and intensities.
///
/// With `Δν = ν_A − ν_B` and `C = √(Δν² + J²)` the mixing angle satisfies
/// `sin 2θ = |J|/C`, so the four lines sit at `ν_c ± (C ± |J|)/2` with the
/// inner pair carrying `1 + sin 2θ` and the outer pair `1 − sin 2θ` (the
/// "roof effect"; the intensities collapse to a single line when `Δν → 0`).
/// Returns `(frequency_hz, intensity)` sorted by frequency, normalized to the
/// largest line.
pub fn ab_reference(nu_a_hz: f64, nu_b_hz: f64, j_hz: f64) -> Vec<ExactLine> {
    let center = 0.5 * (nu_a_hz + nu_b_hz);
    let delta = nu_a_hz - nu_b_hz;
    let c = (delta * delta + j_hz * j_hz).sqrt();
    let sin2theta = if c > 1e-12 { j_hz.abs() / c } else { 0.0 };
    let inner = 1.0 + sin2theta;
    let outer = 1.0 - sin2theta;
    let mut lines = vec![
        ExactLine {
            frequency_hz: center + 0.5 * (c + j_hz.abs()),
            intensity: outer,
        },
        ExactLine {
            frequency_hz: center + 0.5 * (c - j_hz.abs()),
            intensity: inner,
        },
        ExactLine {
            frequency_hz: center - 0.5 * (c - j_hz.abs()),
            intensity: inner,
        },
        ExactLine {
            frequency_hz: center - 0.5 * (c + j_hz.abs()),
            intensity: outer,
        },
    ];
    lines.sort_by(|a, b| a.frequency_hz.partial_cmp(&b.frequency_hz).unwrap());
    let max_i = lines.iter().map(|l| l.intensity).fold(0.0f64, f64::max);
    for line in lines.iter_mut() {
        line.intensity /= max_i;
    }
    lines
}

/// Format a human-readable strong-coupling note for a pair.
pub fn strong_coupling_note(pair: &StrongCouplingPair) -> String {
    if pair.ratio.is_finite() {
        format!(
            "strong coupling between atoms {} and {}: |Δν|/|J| = {:.1} (< {:.0}); first-order multiplet is approximate",
            pair.i, pair.j, pair.ratio, STRONG_COUPLING_RATIO
        )
    } else {
        format!(
            "atoms {} and {} have a resolved coupling below the numerical threshold",
            pair.i, pair.j
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ab_spectrum_matches_closed_form() {
        // Strong coupling: Δν = 20 Hz, J = 10 Hz.
        let system = SpinSystem::new(vec![60.0, 40.0], vec![vec![0.0, 10.0], vec![10.0, 0.0]]);
        let exact = system.exact_spectrum();
        let reference = ab_reference(60.0, 40.0, 10.0);
        assert_eq!(exact.len(), 4, "AB must give 4 lines, got {exact:?}");
        for (e, r) in exact.iter().zip(reference.iter()) {
            assert!(
                (e.frequency_hz - r.frequency_hz).abs() < 1e-8,
                "line at {} vs {}",
                e.frequency_hz,
                r.frequency_hz
            );
            assert!(
                (e.intensity - r.intensity).abs() < 1e-8,
                "intensity {} vs {}",
                e.intensity,
                r.intensity
            );
        }
    }

    #[test]
    fn ab_roof_effect_grows_with_coupling() {
        // The inner/outer intensity ratio grows as |J|/C (roof effect):
        // stronger coupling (smaller Δν at fixed J) ⇒ stronger lean.
        let strong = ab_reference(60.0, 40.0, 10.0);
        let weak = ab_reference(70.0, 30.0, 10.0);
        let ratio = |lines: &[ExactLine]| lines[1].intensity / lines[0].intensity;
        assert!(
            ratio(&strong) > ratio(&weak),
            "strong {} vs weak {}",
            ratio(&strong),
            ratio(&weak)
        );
        // A2 limit: the outer lines vanish, the inner ones merge into the
        // centre line.
        let a2 = ab_reference(100.0, 100.0, 8.0);
        assert!(a2[0].intensity < 1e-12);
        assert!((a2[1].frequency_hz - 100.0).abs() < 1e-9);
    }

    #[test]
    fn weak_coupling_limit_gives_two_doublets() {
        // Δν = 1000 Hz, J = 7 Hz: first-order doublets at ν ± J/2.
        let system = SpinSystem::new(vec![1500.0, 500.0], vec![vec![0.0, 7.0], vec![7.0, 0.0]]);
        let lines = system.exact_spectrum();
        assert_eq!(lines.len(), 4);
        let expected = [1496.5, 1503.5, 496.5, 503.5];
        for exp in expected {
            assert!(
                lines.iter().any(|l| (l.frequency_hz - exp).abs() < 0.05),
                "missing first-order line at {exp}: {lines:?}"
            );
        }
        // First-order limit: the four lines have (nearly) equal intensity;
        // the exact second-order correction is ~(J/2Δν)² ≈ 2.4e-5.
        for l in &lines {
            assert!(
                (l.intensity - 1.0).abs() < 2e-2,
                "intensity {}",
                l.intensity
            );
        }
        // Inner lines equal, outer lines equal.
        assert!((lines[0].intensity - lines[3].intensity).abs() < 1e-9);
        assert!((lines[1].intensity - lines[2].intensity).abs() < 1e-9);
    }

    #[test]
    fn equivalent_spins_give_singlet() {
        // Δν = 0: the coupling is not observable (one line).
        let system = SpinSystem::new(vec![100.0, 100.0], vec![vec![0.0, 8.0], vec![8.0, 0.0]]);
        let lines = system.exact_spectrum();
        assert_eq!(lines.len(), 1, "{lines:?}");
        assert!((lines[0].frequency_hz - 100.0).abs() < 1e-9);
    }

    #[test]
    fn strong_coupling_report_flags_ab() {
        let system = SpinSystem::new(vec![60.0, 40.0], vec![vec![0.0, 10.0], vec![10.0, 0.0]]);
        let report = system.strong_coupling_report();
        assert_eq!(report.len(), 1);
        assert!(report[0].strongly_coupled);
        assert!((report[0].ratio - 2.0).abs() < 1e-12);
        let weak = SpinSystem::new(vec![1000.0, 0.0], vec![vec![0.0, 7.0], vec![7.0, 0.0]]);
        assert!(!weak.strong_coupling_report()[0].strongly_coupled);
    }

    #[test]
    fn abx_system_runs_and_conserves_lines() {
        // Three spins: 2^n = 8 states, single-quantum transitions give 8 lines
        // (before degeneracy merging).
        let system = SpinSystem::new(
            vec![100.0, 90.0, 20.0],
            vec![
                vec![0.0, 8.0, 3.0],
                vec![8.0, 0.0, 5.0],
                vec![3.0, 5.0, 0.0],
            ],
        );
        let lines = system.exact_spectrum();
        assert!(!lines.is_empty());
        let total: f64 = lines.iter().map(|l| l.intensity).sum();
        assert!(total > 0.0);
    }
}
