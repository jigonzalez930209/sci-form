use super::traits::ForceFieldContribution;

// ─── MMFF94 Atom Type Assignment ─────────────────────────────────────────────

/// MMFF94 atom type (simplified subset for common organic elements).
/// Full MMFF94 has 99 types; we cover the most frequent organic types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Mmff94AtomType {
    CR = 1,    // Alkyl carbon, SP3
    CSp2 = 2,  // Vinyl carbon, SP2
    CSp = 3,   // Acetylenic carbon, SP
    CO = 4,    // Carbonyl carbon C=O
    CR4R = 20, // Carbon in 4-membered ring
    CR3R = 22, // Carbon in 3-membered ring
    CB = 37,   // Aromatic carbon
    NR = 8,    // Amine nitrogen SP3
    N2 = 9,    // Imine nitrogen SP2
    NC = 10,   // Isonitrile nitrogen SP
    NAm = 40,  // Amide nitrogen
    NR2 = 39,  // Aromatic nitrogen (pyridine)
    OR = 6,    // Alcohol/ether oxygen
    O2 = 7,    // Carbonyl oxygen
    OX = 32,   // Carboxylate oxygen
    F = 11,    // Fluorine
    Cl = 12,   // Chlorine
    Br = 13,   // Bromine
    I = 14,    // Iodine
    S = 15,    // Thiol sulfur
    SO = 17,   // Sulfoxide S=O
    SO2 = 18,  // Sulfone S(=O)2
    P = 25,    // Phosphorus SP3
    HC = 5,    // Hydrogen on carbon
    HO = 21,   // Hydrogen on oxygen
    HN = 23,   // Hydrogen on nitrogen
    HS = 71,   // Hydrogen on sulfur
    Unknown = 0,
}

/// Assign MMFF94 atom type from element number and hybridization.
pub fn assign_mmff94_type(
    element: u8,
    hyb: &crate::graph::Hybridization,
    is_aromatic: bool,
    is_amide_n: bool,
) -> Mmff94AtomType {
    use crate::graph::Hybridization::*;
    match element {
        1 => Mmff94AtomType::HC, // Default; caller should refine based on neighbor
        5 => Mmff94AtomType::CSp2, // Approximate: boron as sp2 carbon
        6 => {
            if is_aromatic {
                Mmff94AtomType::CB
            } else {
                match hyb {
                    SP => Mmff94AtomType::CSp,
                    SP2 => Mmff94AtomType::CSp2,
                    _ => Mmff94AtomType::CR,
                }
            }
        }
        7 => {
            if is_aromatic {
                Mmff94AtomType::NR2
            } else if is_amide_n {
                Mmff94AtomType::NAm
            } else {
                match hyb {
                    SP => Mmff94AtomType::NC,
                    SP2 => Mmff94AtomType::N2,
                    _ => Mmff94AtomType::NR,
                }
            }
        }
        8 => match hyb {
            SP2 => Mmff94AtomType::O2,
            _ => Mmff94AtomType::OR,
        },
        9 => Mmff94AtomType::F,
        15 => Mmff94AtomType::P,
        16 => Mmff94AtomType::S,
        17 => Mmff94AtomType::Cl,
        35 => Mmff94AtomType::Br,
        53 => Mmff94AtomType::I,
        _ => Mmff94AtomType::Unknown,
    }
}

// ─── MMFF94 Bond Stretching ──────────────────────────────────────────────────

/// MMFF94 bond stretching (quartic form):
/// E = 0.5 · k_b · Δr² · (1 + cs · Δr + 7/12 · cs² · Δr²)
/// where cs = -2.0 Å⁻¹ (cubic stretch constant)
pub struct Mmff94BondStretch {
    pub atom_i: usize,
    pub atom_j: usize,
    pub k_b: f64,  // Force constant (md/Å)
    pub r0: f64,   // Equilibrium bond length (Å)
}

const MMFF94_CUBIC_STRETCH: f64 = -2.0;

impl ForceFieldContribution for Mmff94BondStretch {
    fn evaluate_energy_and_inject_gradient(&self, coords: &[f64], grad: &mut [f64]) -> f64 {
        let ri = self.atom_i * 3;
        let rj = self.atom_j * 3;
        let dx = coords[ri] - coords[rj];
        let dy = coords[ri + 1] - coords[rj + 1];
        let dz = coords[ri + 2] - coords[rj + 2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-8);
        let dr = dist - self.r0;
        let cs = MMFF94_CUBIC_STRETCH;
        let cs2 = cs * cs;

        // Energy: 143.9325 * 0.5 * kb * dr^2 * (1 + cs*dr + 7/12*cs^2*dr^2)
        let energy = 143.9325 * 0.5 * self.k_b * dr * dr * (1.0 + cs * dr + (7.0 / 12.0) * cs2 * dr * dr);

        // Gradient: dE/dr
        let de_dr = 143.9325 * self.k_b * dr * (1.0 + 1.5 * cs * dr + (7.0 / 6.0) * cs2 * dr * dr);
        let scale = de_dr / dist;
        grad[ri] += scale * dx;
        grad[ri + 1] += scale * dy;
        grad[ri + 2] += scale * dz;
        grad[rj] -= scale * dx;
        grad[rj + 1] -= scale * dy;
        grad[rj + 2] -= scale * dz;

        energy
    }
}

// ─── MMFF94 Angle Bending ────────────────────────────────────────────────────

/// MMFF94 angle bending:
/// E = 0.5 · k_a · (Δθ)² · (1 + cb · Δθ)
/// where cb = -0.014 deg⁻¹ (cubic bend constant)
pub struct Mmff94AngleBend {
    pub atom_i: usize,
    pub atom_j: usize, // central
    pub atom_k: usize,
    pub k_a: f64,    // Force constant (md·Å/rad²)
    pub theta0: f64, // Equilibrium angle (radians)
}

const MMFF94_CUBIC_BEND: f64 = -0.014;

impl ForceFieldContribution for Mmff94AngleBend {
    fn evaluate_energy_and_inject_gradient(&self, coords: &[f64], grad: &mut [f64]) -> f64 {
        let ri = self.atom_i * 3;
        let rj = self.atom_j * 3;
        let rk = self.atom_k * 3;

        let rji = [coords[ri] - coords[rj], coords[ri + 1] - coords[rj + 1], coords[ri + 2] - coords[rj + 2]];
        let rjk = [coords[rk] - coords[rj], coords[rk + 1] - coords[rj + 1], coords[rk + 2] - coords[rj + 2]];

        let d_ji = (rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2]).sqrt().max(1e-8);
        let d_jk = (rjk[0] * rjk[0] + rjk[1] * rjk[1] + rjk[2] * rjk[2]).sqrt().max(1e-8);

        let cos_theta = (rji[0] * rjk[0] + rji[1] * rjk[1] + rji[2] * rjk[2]) / (d_ji * d_jk);
        let cos_theta_clamped = cos_theta.clamp(-1.0, 1.0);
        let theta = cos_theta_clamped.acos();
        let dt = (theta - self.theta0) * 180.0 / std::f64::consts::PI; // In degrees for MMFF94

        let cb = MMFF94_CUBIC_BEND;
        let energy = 0.043844 * 0.5 * self.k_a * dt * dt * (1.0 + cb * dt);

        // Gradient (simplified: project along angle bisector normal)
        let de_dtheta = 0.043844 * self.k_a * dt * (1.0 + 1.5 * cb * dt) * (180.0 / std::f64::consts::PI);
        let sin_theta = (1.0 - cos_theta_clamped * cos_theta_clamped).sqrt().max(1e-8);
        let dcos = -1.0 / sin_theta;
        let pref = de_dtheta * dcos;

        for dim in 0..3 {
            let term_i = pref * (rjk[dim] / (d_ji * d_jk) - cos_theta * rji[dim] / (d_ji * d_ji)) / d_ji * d_ji;
            let term_k = pref * (rji[dim] / (d_ji * d_jk) - cos_theta * rjk[dim] / (d_jk * d_jk)) / d_jk * d_jk;
            grad[ri + dim] += term_i;
            grad[rk + dim] += term_k;
            grad[rj + dim] -= term_i + term_k;
        }

        energy
    }
}

// ─── MMFF94 Torsion ──────────────────────────────────────────────────────────

/// MMFF94 torsion:
/// E = 0.5 · (V1·(1+cos φ) + V2·(1-cos 2φ) + V3·(1+cos 3φ))
pub struct Mmff94Torsion {
    pub atom_i: usize,
    pub atom_j: usize,
    pub atom_k: usize,
    pub atom_l: usize,
    pub v1: f64,
    pub v2: f64,
    pub v3: f64,
}

impl ForceFieldContribution for Mmff94Torsion {
    fn evaluate_energy_and_inject_gradient(&self, coords: &[f64], grad: &mut [f64]) -> f64 {
        let p = |idx: usize| -> [f64; 3] {
            [coords[idx * 3], coords[idx * 3 + 1], coords[idx * 3 + 2]]
        };
        let pi = p(self.atom_i);
        let pj = p(self.atom_j);
        let pk = p(self.atom_k);
        let pl = p(self.atom_l);

        // Compute dihedral using standard atan2 method
        let b1 = [pj[0] - pi[0], pj[1] - pi[1], pj[2] - pi[2]];
        let b2 = [pk[0] - pj[0], pk[1] - pj[1], pk[2] - pj[2]];
        let b3 = [pl[0] - pk[0], pl[1] - pk[1], pl[2] - pk[2]];

        let cross = |a: [f64; 3], b: [f64; 3]| -> [f64; 3] {
            [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]
        };
        let dot = |a: [f64; 3], b: [f64; 3]| -> f64 { a[0]*b[0] + a[1]*b[1] + a[2]*b[2] };

        let n1 = cross(b1, b2);
        let n2 = cross(b2, b3);
        let m1 = cross(n1, b2);

        let b2_len = dot(b2, b2).sqrt().max(1e-8);
        let x = dot(n1, n2);
        let y = dot(m1, n2) / b2_len;
        let phi = (-y).atan2(x);

        let energy = 0.5 * (self.v1 * (1.0 + phi.cos())
            + self.v2 * (1.0 - (2.0 * phi).cos())
            + self.v3 * (1.0 + (3.0 * phi).cos()));

        // Numerical gradient for torsion (analytical is complex; use central differences)
        let eps = 1e-5;
        for atom_idx in [self.atom_i, self.atom_j, self.atom_k, self.atom_l] {
            for dim in 0..3 {
                let idx = atom_idx * 3 + dim;
                let orig = coords[idx];
                let mut c_plus = coords.to_vec();
                let mut c_minus = coords.to_vec();
                c_plus[idx] = orig + eps;
                c_minus[idx] = orig - eps;

                let phi_p = {
                    let pi = [c_plus[self.atom_i*3], c_plus[self.atom_i*3+1], c_plus[self.atom_i*3+2]];
                    let pj = [c_plus[self.atom_j*3], c_plus[self.atom_j*3+1], c_plus[self.atom_j*3+2]];
                    let pk = [c_plus[self.atom_k*3], c_plus[self.atom_k*3+1], c_plus[self.atom_k*3+2]];
                    let pl = [c_plus[self.atom_l*3], c_plus[self.atom_l*3+1], c_plus[self.atom_l*3+2]];
                    let b1 = [pj[0]-pi[0], pj[1]-pi[1], pj[2]-pi[2]];
                    let b2 = [pk[0]-pj[0], pk[1]-pj[1], pk[2]-pj[2]];
                    let b3 = [pl[0]-pk[0], pl[1]-pk[1], pl[2]-pk[2]];
                    let nn1 = cross(b1, b2);
                    let nn2 = cross(b2, b3);
                    let mm1 = cross(nn1, b2);
                    let b2l = dot(b2,b2).sqrt().max(1e-8);
                    (-dot(mm1,nn2)/b2l).atan2(dot(nn1,nn2))
                };
                let phi_m = {
                    let pi = [c_minus[self.atom_i*3], c_minus[self.atom_i*3+1], c_minus[self.atom_i*3+2]];
                    let pj = [c_minus[self.atom_j*3], c_minus[self.atom_j*3+1], c_minus[self.atom_j*3+2]];
                    let pk = [c_minus[self.atom_k*3], c_minus[self.atom_k*3+1], c_minus[self.atom_k*3+2]];
                    let pl = [c_minus[self.atom_l*3], c_minus[self.atom_l*3+1], c_minus[self.atom_l*3+2]];
                    let b1 = [pj[0]-pi[0], pj[1]-pi[1], pj[2]-pi[2]];
                    let b2 = [pk[0]-pj[0], pk[1]-pj[1], pk[2]-pj[2]];
                    let b3 = [pl[0]-pk[0], pl[1]-pk[1], pl[2]-pk[2]];
                    let nn1 = cross(b1, b2);
                    let nn2 = cross(b2, b3);
                    let mm1 = cross(nn1, b2);
                    let b2l = dot(b2,b2).sqrt().max(1e-8);
                    (-dot(mm1,nn2)/b2l).atan2(dot(nn1,nn2))
                };

                let e_p = 0.5 * (self.v1 * (1.0 + phi_p.cos()) + self.v2 * (1.0 - (2.0*phi_p).cos()) + self.v3 * (1.0 + (3.0*phi_p).cos()));
                let e_m = 0.5 * (self.v1 * (1.0 + phi_m.cos()) + self.v2 * (1.0 - (2.0*phi_m).cos()) + self.v3 * (1.0 + (3.0*phi_m).cos()));
                grad[idx] += (e_p - e_m) / (2.0 * eps);
            }
        }

        energy
    }
}

// ─── MMFF94 Buffered 14-7 Van der Waals ──────────────────────────────────────

/// Dispersión estérica repulsiva/atractiva regida por el Potencial Amortiguado 14-7 (Buffered 14-7) de Halgren.
pub struct Mmff94BufferedVanDerWaals {
    pub atom_i_idx: usize,
    pub atom_j_idx: usize,
    pub radius_star: f64,   // Parámetro dimensional empírico cruzado R*ij
    pub epsilon_depth: f64, // Factor de profundidad termodinámica eps_ij
}

impl ForceFieldContribution for Mmff94BufferedVanDerWaals {
    fn evaluate_energy_and_inject_gradient(&self, coords: &[f64], grad: &mut [f64]) -> f64 {
        let root_i = self.atom_i_idx * 3;
        let root_j = self.atom_j_idx * 3;

        let delta_x = coords[root_i] - coords[root_j];
        let delta_y = coords[root_i + 1] - coords[root_j + 1];
        let delta_z = coords[root_i + 2] - coords[root_j + 2];

        let dist_squared = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
        let mut dist_r = dist_squared.sqrt();

        // Tope asintótico absoluto inferior para colisiones
        if dist_r < 1e-8 {
            dist_r = 1e-8;
        }

        // Algebra fraccionaria amortiguadora: E_vdW = eps * (1.07 R* / (R + 0.07 R*))^7 * ((1.12 R*^7 / (R^7 + 0.12 R*^7)) - 2)
        let r_star_powered_7 = self.radius_star.powi(7);
        let dist_r_powered_7 = dist_r.powi(7);

        let repulsive_denominator = dist_r + 0.07 * self.radius_star;
        let repulsive_term = (1.07 * self.radius_star / repulsive_denominator).powi(7);

        let attractive_denominator = dist_r_powered_7 + 0.12 * r_star_powered_7;
        let attractive_term = (1.12 * r_star_powered_7 / attractive_denominator) - 2.0;

        let vdw_total_energy = self.epsilon_depth * repulsive_term * attractive_term;

        // Derivación espacial analítica
        let gradient_rep_term = -7.0 * repulsive_term / repulsive_denominator;
        let gradient_attr_term = -7.0 * dist_r.powi(6) * (1.12 * r_star_powered_7)
            / (attractive_denominator * attractive_denominator);

        let force_scalar_magnitude = self.epsilon_depth
            * (gradient_rep_term * attractive_term + repulsive_term * gradient_attr_term);

        // Factorización cartesiana
        let vector_prefactor = force_scalar_magnitude / dist_r;
        let grad_x = vector_prefactor * delta_x;
        let grad_y = vector_prefactor * delta_y;
        let grad_z = vector_prefactor * delta_z;

        grad[root_i] += grad_x;
        grad[root_i + 1] += grad_y;
        grad[root_i + 2] += grad_z;

        grad[root_j] -= grad_x;
        grad[root_j + 1] -= grad_y;
        grad[root_j + 2] -= grad_z;

        vdw_total_energy
    }
}

// ─── Gradient Validation ─────────────────────────────────────────────────────

/// Validate analytical gradients against numerical (central-difference) gradients.
/// Returns max absolute error across all coordinates.
pub fn validate_gradients(
    term: &dyn ForceFieldContribution,
    coords: &[f64],
    eps: f64,
) -> f64 {
    let n = coords.len();
    let mut analytical_grad = vec![0.0; n];
    term.evaluate_energy_and_inject_gradient(coords, &mut analytical_grad);

    let mut max_err = 0.0f64;
    for i in 0..n {
        let mut c_plus = coords.to_vec();
        let mut c_minus = coords.to_vec();
        c_plus[i] += eps;
        c_minus[i] -= eps;

        let mut g_dummy = vec![0.0; n];
        let e_plus = term.evaluate_energy_and_inject_gradient(&c_plus, &mut g_dummy);
        g_dummy.fill(0.0);
        let e_minus = term.evaluate_energy_and_inject_gradient(&c_minus, &mut g_dummy);

        let numerical = (e_plus - e_minus) / (2.0 * eps);
        let err = (analytical_grad[i] - numerical).abs();
        max_err = max_err.max(err);
    }
    max_err
}

// ─── MMFF94 Builder (assemble terms for a molecule) ──────────────────────────

/// Simple bond/angle/torsion lookup parameters for building MMFF94 terms.
/// Uses fallback empirical rules when proper MMFF94 parameter tables are not available.
pub struct Mmff94Builder;

impl Mmff94Builder {
    /// Estimate MMFF94 bond stretching parameters from elements and bond order.
    fn bond_params(elem_i: u8, elem_j: u8, _bond_order: u8) -> (f64, f64) {
        // r0 (Å) from covalent radii sum, kb from empirical rule
        let r_cov = |e: u8| -> f64 {
            match e {
                1 => 0.31, 6 => 0.76, 7 => 0.71, 8 => 0.66, 9 => 0.57,
                15 => 1.07, 16 => 1.05, 17 => 1.02, 35 => 1.20, 53 => 1.39,
                _ => 1.0,
            }
        };
        let r0 = r_cov(elem_i) + r_cov(elem_j);
        let kb = 5.0; // Fallback; real MMFF94 uses typed parameters
        (kb, r0)
    }

    /// Build all MMFF94 force field terms for a parsed molecule.
    ///
    /// `elements`: atomic numbers.
    /// `bonds`: list of (atom_i, atom_j, bond_order).
    /// `coords`: flat xyz coordinates.
    ///
    /// Returns boxed force field contributions.
    pub fn build(
        elements: &[u8],
        bonds: &[(usize, usize, u8)],
    ) -> Vec<Box<dyn ForceFieldContribution>> {
        let n_atoms = elements.len();
        let mut terms: Vec<Box<dyn ForceFieldContribution>> = Vec::new();

        // Bond stretching terms
        for &(i, j, order) in bonds {
            let (kb, r0) = Self::bond_params(elements[i], elements[j], order);
            terms.push(Box::new(Mmff94BondStretch {
                atom_i: i, atom_j: j, k_b: kb, r0,
            }));
        }

        // Angle bending: find all i-j-k where (i,j) and (j,k) are bonded
        let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); n_atoms];
        for &(i, j, _) in bonds {
            neighbors[i].push(j);
            neighbors[j].push(i);
        }
        for j in 0..n_atoms {
            let nbrs = &neighbors[j];
            for a in 0..nbrs.len() {
                for b in (a + 1)..nbrs.len() {
                    let i = nbrs[a];
                    let k = nbrs[b];
                    terms.push(Box::new(Mmff94AngleBend {
                        atom_i: i, atom_j: j, atom_k: k,
                        k_a: 0.5, // Fallback force constant
                        theta0: 109.5_f64.to_radians(), // SP3 default
                    }));
                }
            }
        }

        // Torsion terms: find all i-j-k-l where (i,j), (j,k), (k,l) are bonded
        for &(j, k, _) in bonds {
            for &i in &neighbors[j] {
                if i == k { continue; }
                for &l in &neighbors[k] {
                    if l == j || l == i { continue; }
                    terms.push(Box::new(Mmff94Torsion {
                        atom_i: i, atom_j: j, atom_k: k, atom_l: l,
                        v1: 0.0, v2: 1.0, v3: 0.0, // Generic 2-fold barrier
                    }));
                }
            }
        }

        // 1-4 vdW terms (atoms separated by 3 bonds)
        for &(j, k, _) in bonds {
            for &i in &neighbors[j] {
                if i == k { continue; }
                for &l in &neighbors[k] {
                    if l == j || l == i { continue; }
                    let r_star = 3.5; // Generic
                    let eps = 0.05;
                    terms.push(Box::new(Mmff94BufferedVanDerWaals {
                        atom_i_idx: i, atom_j_idx: l,
                        radius_star: r_star, epsilon_depth: eps,
                    }));
                }
            }
        }

        terms
    }

    /// Compute total energy from all MMFF94 terms.
    pub fn total_energy(
        terms: &[Box<dyn ForceFieldContribution>],
        coords: &[f64],
    ) -> (f64, Vec<f64>) {
        let n = coords.len();
        let mut grad = vec![0.0; n];
        let mut total = 0.0;
        for term in terms {
            total += term.evaluate_energy_and_inject_gradient(coords, &mut grad);
        }
        (total, grad)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mmff94_vdw_energy() {
        let term = Mmff94BufferedVanDerWaals {
            atom_i_idx: 0,
            atom_j_idx: 1,
            radius_star: 3.6,
            epsilon_depth: 0.1,
        };
        let coords = vec![0.0, 0.0, 0.0, 3.6, 0.0, 0.0];
        let mut grad = vec![0.0; 6];
        let e = term.evaluate_energy_and_inject_gradient(&coords, &mut grad);
        // At equilibrium distance, energy should be near -epsilon
        assert!(e < 0.0 && e > -0.2, "vdW energy at R*: {e}");
    }

    #[test]
    fn test_mmff94_vdw_gradient_validation() {
        let term = Mmff94BufferedVanDerWaals {
            atom_i_idx: 0,
            atom_j_idx: 1,
            radius_star: 3.6,
            epsilon_depth: 0.1,
        };
        let coords = vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0];
        let max_err = validate_gradients(&term, &coords, 1e-5);
        assert!(max_err < 1e-4, "vdW gradient error: {max_err}");
    }

    #[test]
    fn test_mmff94_bond_stretch() {
        let term = Mmff94BondStretch {
            atom_i: 0,
            atom_j: 1,
            k_b: 5.0,
            r0: 1.5,
        };
        // At equilibrium: energy ≈ 0
        let coords_eq = vec![0.0, 0.0, 0.0, 1.5, 0.0, 0.0];
        let mut grad = vec![0.0; 6];
        let e_eq = term.evaluate_energy_and_inject_gradient(&coords_eq, &mut grad);
        assert!(e_eq.abs() < 1e-10, "bond stretch at r0: {e_eq}");

        // Stretched: energy > 0
        let coords_stretch = vec![0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        grad.fill(0.0);
        let e_str = term.evaluate_energy_and_inject_gradient(&coords_stretch, &mut grad);
        assert!(e_str > 0.0, "bond stretch energy should be positive: {e_str}");
    }

    #[test]
    fn test_mmff94_bond_stretch_gradient_validation() {
        let term = Mmff94BondStretch {
            atom_i: 0,
            atom_j: 1,
            k_b: 5.0,
            r0: 1.5,
        };
        let coords = vec![0.0, 0.0, 0.0, 2.0, 0.1, 0.0];
        let max_err = validate_gradients(&term, &coords, 1e-5);
        assert!(max_err < 1e-3, "bond stretch gradient error: {max_err}");
    }

    #[test]
    fn test_mmff94_torsion_energy() {
        let term = Mmff94Torsion {
            atom_i: 0,
            atom_j: 1,
            atom_k: 2,
            atom_l: 3,
            v1: 0.0,
            v2: 5.0,
            v3: 0.0,
        };
        // Planar trans: phi ≈ 180°
        let coords = vec![
            -1.5, 1.0, 0.0,
             0.0, 0.0, 0.0,
             1.5, 0.0, 0.0,
             3.0, 1.0, 0.0,
        ];
        let mut grad = vec![0.0; 12];
        let e = term.evaluate_energy_and_inject_gradient(&coords, &mut grad);
        assert!(e.is_finite(), "torsion energy should be finite: {e}");
    }

    #[test]
    fn test_mmff94_atom_typing() {
        use crate::graph::Hybridization;
        let t = assign_mmff94_type(6, &Hybridization::SP3, false, false);
        assert_eq!(t, Mmff94AtomType::CR);
        let t = assign_mmff94_type(6, &Hybridization::SP2, true, false);
        assert_eq!(t, Mmff94AtomType::CB);
        let t = assign_mmff94_type(7, &Hybridization::SP3, false, true);
        assert_eq!(t, Mmff94AtomType::NAm);
    }

    #[test]
    fn test_mmff94_builder_ethane() {
        // Ethane: C-C with 6 hydrogens
        let elements = vec![6, 6, 1, 1, 1, 1, 1, 1]; // C, C, H×6
        let bonds = vec![
            (0, 1, 1), // C-C
            (0, 2, 1), (0, 3, 1), (0, 4, 1), // C-H
            (1, 5, 1), (1, 6, 1), (1, 7, 1), // C-H
        ];
        // Staggered ethane coordinates (approximate)
        let coords = vec![
            0.0,  0.0,  0.0,   // C0
            1.54, 0.0,  0.0,   // C1
           -0.5,  0.9,  0.0,   // H
           -0.5, -0.9,  0.0,   // H
           -0.5,  0.0,  0.9,   // H
            2.04, 0.9,  0.0,   // H
            2.04, -0.9, 0.0,   // H
            2.04, 0.0,  0.9,   // H
        ];
        let terms = Mmff94Builder::build(&elements, &bonds);
        assert!(!terms.is_empty(), "should produce force field terms");

        let (energy, grad) = Mmff94Builder::total_energy(&terms, &coords);
        assert!(energy.is_finite(), "total energy should be finite: {energy}");
        assert!(grad.iter().all(|g| g.is_finite()), "all gradients should be finite");
    }

    #[test]
    fn test_mmff94_builder_gradient_consistency() {
        // Verify total gradient is consistent with numerical for a simple system
        let elements = vec![6, 6];
        let bonds = vec![(0, 1, 1)];
        let coords = vec![0.0, 0.0, 0.0, 1.6, 0.1, 0.0];
        let terms = Mmff94Builder::build(&elements, &bonds);
        let (_, analytical_grad) = Mmff94Builder::total_energy(&terms, &coords);

        let eps = 1e-5;
        for i in 0..coords.len() {
            let mut cp = coords.clone();
            let mut cm = coords.clone();
            cp[i] += eps;
            cm[i] -= eps;
            let (ep, _) = Mmff94Builder::total_energy(&terms, &cp);
            let (em, _) = Mmff94Builder::total_energy(&terms, &cm);
            let numerical = (ep - em) / (2.0 * eps);
            let err = (analytical_grad[i] - numerical).abs();
            assert!(
                err < 0.1,
                "gradient mismatch at coord {i}: anal={:.6} num={:.6} err={:.6}",
                analytical_grad[i], numerical, err
            );
        }
    }
}
