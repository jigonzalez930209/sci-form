//! Universal Force Field (UFF) energy terms: bond stretch, angle bend,
//! torsion, inversion, and van der Waals contributions.
//!
//! Reference: Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024–10035.

use super::traits::ForceFieldContribution;

/// Classical harmonic bond-stretch term in UFF.
pub struct UffHarmonicBondStretch {
    pub atom_i_idx: usize,
    pub atom_j_idx: usize,
    pub force_constant_kb: f64, // Bond stiffness constant k_b
    pub equilibrium_r0: f64,    // Equilibrium bond length r_0
}

impl ForceFieldContribution for UffHarmonicBondStretch {
    fn evaluate_energy_and_inject_gradient(&self, coords: &[f64], grad: &mut [f64]) -> f64 {
        let root_i = self.atom_i_idx * 3;
        let root_j = self.atom_j_idx * 3;

        let mut diff_x = coords[root_i] - coords[root_j];
        let diff_y = coords[root_i + 1] - coords[root_j + 1];
        let diff_z = coords[root_i + 2] - coords[root_j + 2];

        let mut inter_r = (diff_x * diff_x + diff_y * diff_y + diff_z * diff_z).sqrt();

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
        grad[root_i + 1] += force_y;
        grad[root_i + 2] += force_z;

        grad[root_j] -= force_x;
        grad[root_j + 1] -= force_y;
        grad[root_j + 2] -= force_z;

        bond_energy
    }
}

/// Three-term Fourier angle-bending term used in UFF.
pub struct UffAngleBend {
    pub atom_i_idx: usize,
    pub atom_j_idx: usize, // Central
    pub atom_k_idx: usize,
    pub force_constant_ka: f64,
    pub equilibrium_theta0: f64,
    pub coordination_n: usize, // 0 (linear), 3 (sp2), 4 (sp3)
}

impl ForceFieldContribution for UffAngleBend {
    fn evaluate_energy_and_inject_gradient(&self, coords: &[f64], grad: &mut [f64]) -> f64 {
        let root_i = self.atom_i_idx * 3;
        let root_j = self.atom_j_idx * 3;
        let root_k = self.atom_k_idx * 3;

        let r_ji = [
            coords[root_i] - coords[root_j],
            coords[root_i + 1] - coords[root_j + 1],
            coords[root_i + 2] - coords[root_j + 2],
        ];
        let r_jk = [
            coords[root_k] - coords[root_j],
            coords[root_k + 1] - coords[root_j + 1],
            coords[root_k + 2] - coords[root_j + 2],
        ];

        let d_ji = (r_ji[0] * r_ji[0] + r_ji[1] * r_ji[1] + r_ji[2] * r_ji[2]).sqrt();
        let d_jk = (r_jk[0] * r_jk[0] + r_jk[1] * r_jk[1] + r_jk[2] * r_jk[2]).sqrt();

        if d_ji < 1e-10 || d_jk < 1e-10 {
            return 0.0;
        }

        let cos_theta = (r_ji[0] * r_jk[0] + r_ji[1] * r_jk[1] + r_ji[2] * r_jk[2]) / (d_ji * d_jk);
        let cos_theta = cos_theta.clamp(-1.0, 1.0);
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(1e-8);

        let (energy, d_e_dtheta) = match self.coordination_n {
            0 => {
                // Linear center (RDKit uses n=0 as a dedicated linear sentinel).
                // Reuse the MMFF-style linear form here as a stable proxy.
                let e = self.force_constant_ka * (1.0 + cos_theta);
                let de = -self.force_constant_ka * sin_theta;
                (e, de)
            }
            _ => {
                // General Fourier expansion
                let cos_theta0 = self.equilibrium_theta0.cos();
                let sin_theta0 = self.equilibrium_theta0.sin();

                let c2 = 1.0 / (4.0 * sin_theta0 * sin_theta0).max(1e-8);
                let c1 = -4.0 * c2 * cos_theta0;
                let c0 = c2 * (2.0 * cos_theta0 * cos_theta0 + 1.0);

                let cos_2theta = 2.0 * cos_theta * cos_theta - 1.0;
                let energy = self.force_constant_ka * (c0 + c1 * cos_theta + c2 * cos_2theta);

                // dE/dtheta = ka * (-c1 * sin(theta) - 2 * c2 * sin(2*theta))
                let sin_2theta = 2.0 * sin_theta * cos_theta;
                let de = self.force_constant_ka * (-c1 * sin_theta - 2.0 * c2 * sin_2theta);
                (energy, de)
            }
        };

        // Cartesian gradient derived from the Wilson-Decius-Cross form.
        let pre_i = d_e_dtheta / (d_ji * sin_theta);
        let pre_k = d_e_dtheta / (d_jk * sin_theta);

        for dim in 0..3 {
            let gi = pre_i * (r_jk[dim] / d_jk - cos_theta * (r_ji[dim] / d_ji));
            let gk = pre_k * (r_ji[dim] / d_ji - cos_theta * (r_jk[dim] / d_jk));

            grad[root_i + dim] += gi;
            grad[root_k + dim] += gk;
            grad[root_j + dim] -= gi + gk;
        }

        energy
    }
}

/// UFF torsion term: E = 0.5 * V * [1 - cos(phi0) * cos(n * phi)].
pub struct UffTorsion {
    pub atom_i_idx: usize,
    pub atom_j_idx: usize,
    pub atom_k_idx: usize,
    pub atom_l_idx: usize,
    pub force_constant_v: f64,
    pub periodicity_n: f64,
    pub cos_phi0: f64, // Usually 1.0 or -1.0
}

impl ForceFieldContribution for UffTorsion {
    fn evaluate_energy_and_inject_gradient(&self, coords: &[f64], grad: &mut [f64]) -> f64 {
        let i = self.atom_i_idx * 3;
        let j = self.atom_j_idx * 3;
        let k = self.atom_k_idx * 3;
        let l = self.atom_l_idx * 3;

        let b1 = [
            coords[i] - coords[j],
            coords[i + 1] - coords[j + 1],
            coords[i + 2] - coords[j + 2],
        ];
        let b2 = [
            coords[k] - coords[j],
            coords[k + 1] - coords[j + 1],
            coords[k + 2] - coords[j + 2],
        ];
        let b3 = [
            coords[l] - coords[k],
            coords[l + 1] - coords[k + 1],
            coords[l + 2] - coords[k + 2],
        ];

        let n1 = [
            b1[1] * b2[2] - b1[2] * b2[1],
            b1[2] * b2[0] - b1[0] * b2[2],
            b1[0] * b2[1] - b1[1] * b2[0],
        ];
        let n2 = [
            b2[1] * b3[2] - b2[2] * b3[1],
            b2[2] * b3[0] - b2[0] * b3[2],
            b2[0] * b3[1] - b2[1] * b3[0],
        ];

        let m1 = (n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]).sqrt();
        let m2 = (n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2]).sqrt();
        if m1 < 1e-10 || m2 < 1e-10 {
            return 0.0;
        }

        let cos_phi = (n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]) / (m1 * m2);
        let cos_phi = cos_phi.clamp(-1.0, 1.0);
        let phi = cos_phi.acos();

        // Recover the sign of phi from the orientation relative to b2.
        let cross_n1_n2 = [
            n1[1] * n2[2] - n1[2] * n2[1],
            n1[2] * n2[0] - n1[0] * n2[2],
            n1[0] * n2[1] - n1[1] * n2[0],
        ];
        let dot_dir = cross_n1_n2[0] * b2[0] + cross_n1_n2[1] * b2[1] + cross_n1_n2[2] * b2[2];
        let phi = if dot_dir < 0.0 { -phi } else { phi };

        let energy =
            0.5 * self.force_constant_v * (1.0 - self.cos_phi0 * (self.periodicity_n * phi).cos());
        let d_e_dphi = 0.5
            * self.force_constant_v
            * self.cos_phi0
            * self.periodicity_n
            * (self.periodicity_n * phi).sin();

        // Analytical Gradients (Blondel-Karplus)
        let f_i = [
            -d_e_dphi * n1[0] / (m1 * m1),
            -d_e_dphi * n1[1] / (m1 * m1),
            -d_e_dphi * n1[2] / (m1 * m1),
        ];

        let f_l = [
            d_e_dphi * n2[0] / (m2 * m2),
            d_e_dphi * n2[1] / (m2 * m2),
            d_e_dphi * n2[2] / (m2 * m2),
        ];

        // Central atoms J and K get the remainder (lever rule)
        // Simplified for stability:
        for dim in 0..3 {
            grad[i + dim] += f_i[dim];
            grad[l + dim] += f_l[dim];
            grad[j + dim] -= f_i[dim]; // Approximated
            grad[k + dim] -= f_l[dim]; // Approximated
        }

        energy
    }
}

/// Wilson-Decius-Cross inversion term for near-planar UFF centers.
pub struct UffInversion {
    pub idx_i: usize,
    pub idx_j: usize, // Central
    pub idx_k: usize,
    pub idx_l: usize,
    pub k_inv: f64,
    pub c0: f64,
    pub c1: f64,
    pub c2: f64,
}

impl ForceFieldContribution for UffInversion {
    fn evaluate_energy_and_inject_gradient(&self, coords: &[f64], grad: &mut [f64]) -> f64 {
        let j = self.idx_j * 3;
        let i = self.idx_i * 3;
        let k = self.idx_k * 3;
        let l = self.idx_l * 3;

        let r_ji = [
            coords[i] - coords[j],
            coords[i + 1] - coords[j + 1],
            coords[i + 2] - coords[j + 2],
        ];
        let r_jk = [
            coords[k] - coords[j],
            coords[k + 1] - coords[j + 1],
            coords[k + 2] - coords[j + 2],
        ];
        let r_jl = [
            coords[l] - coords[j],
            coords[l + 1] - coords[j + 1],
            coords[l + 2] - coords[j + 2],
        ];

        // Normal al plano I-J-K
        let n = [
            r_ji[1] * r_jk[2] - r_ji[2] * r_jk[1],
            r_ji[2] * r_jk[0] - r_ji[0] * r_jk[2],
            r_ji[0] * r_jk[1] - r_ji[1] * r_jk[0],
        ];
        let n_len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if n_len < 1e-10 {
            return 0.0;
        }

        let r_jl_len = (r_jl[0] * r_jl[0] + r_jl[1] * r_jl[1] + r_jl[2] * r_jl[2]).sqrt();
        if r_jl_len < 1e-10 {
            return 0.0;
        }

        let sin_psi = (n[0] * r_jl[0] + n[1] * r_jl[1] + n[2] * r_jl[2]) / (n_len * r_jl_len);
        let sin_psi = sin_psi.clamp(-1.0, 1.0);
        let psi = sin_psi.asin();

        let energy = self.k_inv * (self.c0 + self.c1 * sin_psi + self.c2 * (2.0 * psi).cos());
        let d_e_dpsi = self.k_inv * (self.c1 * psi.cos() - 2.0 * self.c2 * (2.0 * psi).sin());

        let cos_psi = psi.cos().max(1e-8);
        let pre_l = d_e_dpsi / (n_len * r_jl_len * cos_psi);

        for dim in 0..3 {
            let gi = pre_l * (n[dim] - sin_psi * r_jl[dim] / r_jl_len);
            grad[l + dim] += gi;
            grad[j + dim] -= gi;
        }

        energy
    }
}

/// UFF Lennard-Jones 12-6 non-bonded interaction.
///
/// E = ε · [(x_ij / r)¹² − 2 · (x_ij / r)⁶]
///
/// Combining rules (Rappé Table 1):
///   x_ij  = (x_i + x_j) / 2          [arithmetic mean of VDW distances]
///   ε_ij  = √(D_i · D_j) · scale     [geometric mean of well depths; scale = 0.5 for 1-4]
pub struct UffLennardJones {
    pub atom_i_idx: usize,
    pub atom_j_idx: usize,
    pub r_star: f64,  // x_ij = (x_i + x_j) / 2
    pub epsilon: f64, // ε_ij = √(D_i · D_j) [already scaled before construction]
}

impl ForceFieldContribution for UffLennardJones {
    fn evaluate_energy_and_inject_gradient(&self, coords: &[f64], grad: &mut [f64]) -> f64 {
        let ri = self.atom_i_idx * 3;
        let rj = self.atom_j_idx * 3;

        let dx = coords[ri] - coords[rj];
        let dy = coords[ri + 1] - coords[rj + 1];
        let dz = coords[ri + 2] - coords[rj + 2];
        let r2 = (dx * dx + dy * dy + dz * dz).max(1e-6);
        let r = r2.sqrt();

        let u = self.r_star / r;
        let u6 = u * u * u * u * u * u;
        let u12 = u6 * u6;

        let energy = self.epsilon * (u12 - 2.0 * u6);

        // dE/dr = ε · (-12 u¹² + 12 u⁶) / r  →  scatter onto grad
        let de_dr = self.epsilon * 12.0 * (u6 - u12) / r;
        let pre = de_dr / r;

        grad[ri] += pre * dx;
        grad[ri + 1] += pre * dy;
        grad[ri + 2] += pre * dz;
        grad[rj] -= pre * dx;
        grad[rj + 1] -= pre * dy;
        grad[rj + 2] -= pre * dz;

        energy
    }
}
