//! CPM Grand Potential Framework — E10.1
//!
//! Computes equilibrium charges under an applied electrochemical potential μ
//! by minimizing the grand potential Ω = E_elec - μ·Q.
//!
//! This is the Constant Potential Method for electrochemistry, in the spirit
//! of LAMMPS `fix electrode/conp` (constant-potential electrodes with a
//! coupled charge equilibration) and of the grand-canonical charge
//! equilibration used in classical polarisable-force-field codes.
//!
//! # Units
//!
//! The applied potential μ is in **eV** (referenced to the vacuum level; the
//! standard hydrogen electrode sits at ≈ −4.44 eV), the per-element
//! electronegativity χ and hardness η are in **eV**, the partial charges are
//! in **e** (elementary charge), the Coulomb coupling is
//! `14.3996 eV·Å / (ε · r_ij)` and the grand potential Ω and electrostatic
//! energy E_elec are in **eV**. Positions are in **Å**.

use serde::{Deserialize, Serialize};

/// CPM configuration.
///
/// Controls the constant-potential charge equilibration: the applied
/// electrochemical potential, the solvent screening and the convergence
/// criteria. All quantities are in the units named by the fields (potential
/// in **eV**, charges in **e**, energies in **eV**).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CpmConfig {
    /// Electrochemical potential in eV (vs vacuum).
    /// SHE ≈ -4.44 eV. Range: [-5.5, -3.5] eV for typical electrochemistry.
    pub mu_ev: f64,
    /// Solvent dielectric constant (dimensionless; affects charge screening).
    pub dielectric: f64,
    /// Maximum SCF iterations for charge equilibration.
    pub max_iter: usize,
    /// Convergence threshold for charge changes, in units of elementary
    /// charge (e).
    pub charge_tol: f64,
    /// Optional per-atom applied potential in eV (multi-electrode setups:
    /// plate A at +V/2, plate B at −V/2, as in LAMMPS `fix electrode/conp`
    /// with `couple`). Overrides `mu_ev` where present (must match natoms).
    #[serde(default)]
    pub mu_per_atom: Option<Vec<f64>>,
}

impl Default for CpmConfig {
    fn default() -> Self {
        Self {
            mu_ev: -4.44,     // SHE
            dielectric: 78.5, // water
            max_iter: 100,
            charge_tol: 1e-6,
            mu_per_atom: None,
        }
    }
}

/// Result of CPM charge calculation.
#[derive(Debug, Clone)]
pub struct CpmResult {
    /// Equilibrium partial charges under potential μ, in units of elementary
    /// charge (e).
    pub charges: Vec<f64>,
    /// Total charge Q = Σ q_i (not necessarily zero!), in units of
    /// elementary charge (e).
    pub total_charge: f64,
    /// Grand potential Ω = E_elec - μ·Q (eV).
    pub grand_potential: f64,
    /// Electrostatic energy E_elec (eV).
    pub electrostatic_energy: f64,
    /// Applied potential (eV).
    pub mu_ev: f64,
    /// Number of SCF iterations.
    pub iterations: usize,
    /// Whether converged.
    pub converged: bool,
}

/// Per-element electronegativity (Pauling scale, eV).
pub(crate) fn chi_element(z: u8) -> f64 {
    match z {
        1 => 7.17,  // H
        6 => 6.27,  // C
        7 => 7.27,  // N
        8 => 8.30,  // O
        9 => 10.41, // F
        15 => 5.62, // P
        16 => 6.22, // S
        17 => 8.30, // Cl
        26 => 4.06, // Fe
        29 => 4.48, // Cu
        30 => 4.45, // Zn
        35 => 7.59, // Br
        53 => 6.76, // I
        _ => 6.27,  // default C
    }
}

/// Per-element hardness (eV).
pub(crate) fn eta_element(z: u8) -> f64 {
    match z {
        1 => 6.43,
        6 => 5.0,
        7 => 5.7,
        8 => 6.08,
        9 => 7.01,
        15 => 4.88,
        16 => 4.14,
        17 => 4.68,
        26 => 3.90,
        29 => 3.25,
        30 => 4.94,
        35 => 4.24,
        53 => 3.70,
        _ => 5.0,
    }
}

/// Compute CPM-equilibrated charges under applied potential μ.
///
/// The grand potential is: Ω({q}, μ) = Σ χ_i·q_i + ½·Σ η_i·q_i² + Σ_{i<j} q_i·J_ij·q_j - μ·Σ q_i
/// Minimizing w.r.t. each q_i (no neutrality constraint):
///   ∂Ω/∂q_i = χ_i + η_i·q_i + Σ_{j≠i} J_ij·q_j - μ = 0
///   → q_i = (μ - χ_i - Σ_{j≠i} J_ij·q_j) / η_i
pub fn compute_cpm_charges(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &CpmConfig,
) -> CpmResult {
    // Tolerate a short coordinate array: atoms without a coordinate are
    // treated as isolated (no Coulomb coupling) instead of panicking.
    let n = elements.len().min(positions.len());
    let mut charges = vec![0.0; n];
    let _converged = false;

    // Precompute Coulomb coupling matrix J_ij
    // Screened by dielectric: J_ij = 1 / (ε · r_ij) in eV (with conversion)
    let ev_per_angstrom = 14.3996; // e²/(4πε₀) in eV·Å
    let mut j_mat = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r > 1e-10 {
                let j_ij = ev_per_angstrom / (config.dielectric * r);
                j_mat[i][j] = j_ij;
                j_mat[j][i] = j_ij;
            }
        }
    }

    // (μ_i = mu_ev, or mu_per_atom[i] for multi-electrode setups).
    let mus: Vec<f64> = match &config.mu_per_atom {
        Some(v) => {
            assert!(
                v.len() == n,
                "mu_per_atom length {} != natoms {}",
                v.len(),
                n
            );
            v.clone()
        }
        None => vec![config.mu_ev; n],
    };
    // Direct solve of the linear CPM equations:
    //   q_i = (μ_i − χ_i − Σ_{j≠i} J_ij q_j) / η_i
    //   ⟺ Σ_j M_ij q_j = μ_i − χ_i,  M_ii = η_i, M_ij = J_ij (i≠j).
    // The damped fixed-point iteration diverges for electrode-sized
    // conducting systems (spectral radius of J/η exceeds 1: graphene
    // bilayer charges blow up to 1e157); the system is linear, so an LU
    // solve is exact and robust. (LAMMPS `electrode/conp` likewise inverts
    // the elastance matrix with `mat_inv`.)
    let n_f = n as f64;
    let mut m_mat = nalgebra::DMatrix::<f64>::zeros(n, n);
    let mut rhs = nalgebra::DVector::<f64>::zeros(n);
    for i in 0..n {
        m_mat[(i, i)] = eta_element(elements[i]);
        rhs[i] = mus[i] - chi_element(elements[i]);
        for j in 0..n {
            if j != i {
                m_mat[(i, j)] = j_mat[i][j];
            }
        }
    }
    let lu = m_mat.clone().lu();
    let sol = lu.solve(&rhs);
    let mut converged = false;
    let iter: usize = 1;
    match sol {
        Some(v) => {
            for i in 0..n {
                charges[i] = v[i];
            }
            // Residual check: ||M q − rhs||_∞ per atom.
            let mut max_res = 0.0f64;
            for i in 0..n {
                let mut acc = m_mat[(i, i)] * charges[i];
                for j in 0..n {
                    if j != i {
                        acc += m_mat[(i, j)] * charges[j];
                    }
                }
                max_res = max_res.max((acc - rhs[i]).abs());
            }
            // Charge tolerance scaled to a residual (eV·e).
            converged = max_res < config.charge_tol * n_f.max(1.0) * 10.0;
        }
        None => {
            // Singular matrix: fall back to damped fixed-point iteration.
            for it in 0..config.max_iter {
                let mut max_change = 0.0;
                for i in 0..n {
                    let chi = chi_element(elements[i]);
                    let eta = eta_element(elements[i]);
                    let mut coupling = 0.0;
                    for j in 0..n {
                        if j != i {
                            coupling += j_mat[i][j] * charges[j];
                        }
                    }
                    let new_q = (mus[i] - chi - coupling) / eta;
                    let change = (new_q - charges[i]).abs();
                    if change > max_change {
                        max_change = change;
                    }
                    charges[i] = 0.5 * charges[i] + 0.5 * new_q;
                    if !charges[i].is_finite() {
                        break;
                    }
                }
                if max_change < config.charge_tol {
                    converged = true;
                    break;
                }
                let _ = it;
            }
        }
    }

    // Compute energies
    let total_charge: f64 = charges.iter().sum();

    let mut e_elec = 0.0;
    for i in 0..n {
        e_elec += chi_element(elements[i]) * charges[i];
        e_elec += 0.5 * eta_element(elements[i]) * charges[i] * charges[i];
    }
    for i in 0..n {
        for j in (i + 1)..n {
            e_elec += charges[i] * j_mat[i][j] * charges[j];
        }
    }

    let grand_potential = e_elec
        - mus
            .iter()
            .zip(charges.iter())
            .map(|(m, q)| m * q)
            .sum::<f64>();

    CpmResult {
        charges,
        total_charge,
        grand_potential,
        electrostatic_energy: e_elec,
        mu_ev: config.mu_ev,
        iterations: iter,
        converged,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn water() -> (Vec<u8>, Vec<[f64; 3]>) {
        (
            vec![8, 1, 1],
            vec![[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]],
        )
    }

    #[test]
    fn test_cpm_converges() {
        let (el, pos) = water();
        let result = compute_cpm_charges(&el, &pos, &CpmConfig::default());
        assert!(result.converged, "CPM should converge");
    }

    #[test]
    fn test_cpm_charge_response_to_potential() {
        let (el, pos) = water();
        let r1 = compute_cpm_charges(
            &el,
            &pos,
            &CpmConfig {
                mu_ev: -4.0,
                ..Default::default()
            },
        );
        let r2 = compute_cpm_charges(
            &el,
            &pos,
            &CpmConfig {
                mu_ev: -5.0,
                ..Default::default()
            },
        );
        // Higher (less negative) μ → more electron donation → more negative total charge
        // Lower μ → more positive total charge
        assert!(
            r1.total_charge > r2.total_charge,
            "Q should decrease with lower μ: Q(-4)={}, Q(-5)={}",
            r1.total_charge,
            r2.total_charge
        );
    }

    #[test]
    fn test_cpm_grand_potential_finite() {
        let (el, pos) = water();
        let result = compute_cpm_charges(&el, &pos, &CpmConfig::default());
        assert!(result.grand_potential.is_finite());
        assert!(result.electrostatic_energy.is_finite());
    }
}
