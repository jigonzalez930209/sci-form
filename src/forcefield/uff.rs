//! Universal Force Field (UFF) energy terms: bond stretch, angle bend,
//! torsion, inversion, and van der Waals contributions.
//!
//! Reference: Rappé et al., *J. Am. Chem. Soc.* **1992**, 114, 10024–10035.
//! The functional forms and the parameter lookups follow the canonical RDKit
//! implementation (`ForceFieldHelpers/UFF`) so that energies can be
//! cross-validated against it.
//!
//! # Units
//!
//! Energies are in **kcal/mol**, gradients are `dE/dR` in **kcal/mol/Å**,
//! bond lengths and van der Waals radii in **Å**, and angles in **radians**.

use super::traits::ForceFieldContribution;

/// Classical harmonic bond-stretch term in UFF.
///
/// `E = ½ · k_b · (r - r_0)²` with `r` the distance between the two atoms
/// in Å. The energy is in kcal/mol and the injected gradient is
/// `dE/dR` in kcal/mol/Å.
///
/// Reference: Rappé et al., *J. Am. Chem. Soc.* **1992**, 114, 10024, eq. 6.
pub struct UffHarmonicBondStretch {
    /// Index of the first bonded atom.
    pub atom_i_idx: usize,
    /// Index of the second bonded atom.
    pub atom_j_idx: usize,
    /// Bond force constant `k_b` in kcal/mol/Å².
    pub force_constant_kb: f64,
    /// Equilibrium bond length `r_0` in Å.
    pub equilibrium_r0: f64,
}

impl ForceFieldContribution for UffHarmonicBondStretch {
    fn kind(&self) -> &'static str {
        "bond"
    }
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

/// Angle-bending term used in UFF.
///
/// The functional form is selected by [`Self::order`]: the general Fourier
/// polynomial (`order = 0`) or the `E = k_a/n² · (1 - cos(nθ))` family used
/// for linear, trigonal-planar and square-planar/octahedral centres. The
/// equilibrium angle is in **radians**, the force constant `k_a` in
/// kcal/mol/rad², the energy in kcal/mol and the injected gradient in
/// kcal/mol/Å.
///
/// Reference: Rappé et al., *J. Am. Chem. Soc.* **1992**, 114, 10024, eq. 13.
pub struct UffAngleBend {
    /// Index of the first terminal atom of the angle.
    pub atom_i_idx: usize,
    /// Index of the central atom of the angle.
    pub atom_j_idx: usize,
    /// Index of the second terminal atom of the angle.
    pub atom_k_idx: usize,
    /// Angle force constant `k_a` in kcal/mol/rad².
    pub force_constant_ka: f64,
    /// Equilibrium valence angle `θ_0` in radians.
    pub equilibrium_theta0: f64,
    /// Functional-form selector: `0` = Fourier polynomial, `1` = linear,
    /// `2` = double-angle, `3` = trigonal planar,
    /// `4` = square planar/octahedral (RDKit `AngleBendContrib::order`).
    pub order: usize,
}

impl ForceFieldContribution for UffAngleBend {
    fn kind(&self) -> &'static str {
        "angle"
    }
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
        let sin_sq = 1.0 - cos_theta * cos_theta;

        let (mut energy, mut d_e_dtheta) = if self.order == 0 {
            // General Fourier expansion (RDKit order 0).
            let cos_theta0 = self.equilibrium_theta0.cos();
            let sin_theta0 = self.equilibrium_theta0.sin();

            let c2 = 1.0 / (4.0 * sin_theta0 * sin_theta0).max(1e-8);
            let c1 = -4.0 * c2 * cos_theta0;
            let c0 = c2 * (2.0 * cos_theta0 * cos_theta0 + 1.0);

            let cos_2theta = 2.0 * cos_theta * cos_theta - 1.0;
            let energy = self.force_constant_ka * (c0 + c1 * cos_theta + c2 * cos_2theta);

            let sin_2theta = 2.0 * sin_theta * cos_theta;
            let de = self.force_constant_ka * (-c1 * sin_theta - 2.0 * c2 * sin_2theta);
            (energy, de)
        } else {
            // E = ka / n^2 * (1 - cos(nθ)) (RDKit orders 1..4; order 1 uses
            // `-cosθ` so the minimum sits at 180°).
            let n = self.order as f64;
            let cos_n = match self.order {
                1 => -cos_theta,
                2 => 2.0 * cos_theta * cos_theta - 1.0,
                3 => cos_theta * (cos_theta * cos_theta - 3.0 * sin_sq),
                4 => {
                    cos_theta * cos_theta * cos_theta * cos_theta
                        - 6.0 * cos_theta * cos_theta * sin_sq
                        + sin_sq * sin_sq
                }
                _ => cos_theta,
            };
            let sin_n = match self.order {
                1 => -sin_theta,
                2 => 2.0 * sin_theta * cos_theta,
                3 => sin_theta * (3.0 - 4.0 * sin_sq),
                4 => cos_theta * sin_theta * (4.0 - 8.0 * sin_sq),
                _ => sin_theta,
            };
            let e = self.force_constant_ka * (1.0 - cos_n) / (n * n);
            let de = self.force_constant_ka / n * sin_n;
            (e, de)
        };

        // Empirical penalty for near-zero angles (RDKit orders 1..4 only).
        if (1..5).contains(&self.order) && cos_theta > 0.8660 {
            let theta = cos_theta.acos();
            let pen = (-20.0 * (theta - self.equilibrium_theta0 + 0.25)).exp();
            energy += pen;
            d_e_dtheta += -20.0 * pen;
        }

        // Cartesian gradient: dE/dr = dE/dtheta * dtheta/dr with
        // dtheta/dr = -d(cos)/dr/sin (acos derivative). The minus sign is
        // essential; without it the gradient points uphill (found by F70
        // validation: 46 kcal/mol/A disagreement with finite differences).
        let pre_i = -d_e_dtheta / (d_ji * sin_theta);
        let pre_k = -d_e_dtheta / (d_jk * sin_theta);

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

/// UFF torsion term: `E = ½ · V · [1 - cos(φ₀) · cos(n·φ)]`.
///
/// `φ` is the signed dihedral angle about the j–k bond in radians, `V` is the
/// torsion barrier in kcal/mol, `n` the periodicity (dimensionless) and
/// `cos_phi0` the phase factor (dimensionless, usually `+1` or `-1`). The
/// energy is in kcal/mol and the injected gradient in kcal/mol/Å.
pub struct UffTorsion {
    /// Index of the first terminal atom of the dihedral.
    pub atom_i_idx: usize,
    /// Index of the second atom of the dihedral.
    pub atom_j_idx: usize,
    /// Index of the third atom of the dihedral.
    pub atom_k_idx: usize,
    /// Index of the fourth atom of the dihedral.
    pub atom_l_idx: usize,
    /// Torsion barrier `V` in kcal/mol.
    pub force_constant_v: f64,
    /// Torsion periodicity `n` (dimensionless, e.g. 2 or 3).
    pub periodicity_n: f64,
    /// Phase factor `cos(φ₀)` (dimensionless, usually `1.0` or `-1.0`).
    pub cos_phi0: f64,
}

impl ForceFieldContribution for UffTorsion {
    fn kind(&self) -> &'static str {
        "torsion"
    }
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

        // The canonical implementation builds the second normal from
        // r3 = p2 - p3 (opposite of our b2-based normal), which shifts the
        // unsigned dihedral by 180 degrees. For odd periodicities that flips
        // cos(n*phi); absorb it into an effective cos_term.
        let period = self.periodicity_n.round() as i64;
        let eff_cos_phi0 = if period.rem_euclid(2) == 1 {
            -self.cos_phi0
        } else {
            self.cos_phi0
        };

        let energy =
            0.5 * self.force_constant_v * (1.0 - eff_cos_phi0 * (self.periodicity_n * phi).cos());
        let d_e_dphi = 0.5
            * self.force_constant_v
            * eff_cos_phi0
            * self.periodicity_n
            * (self.periodicity_n * phi).sin();

        // Exact Cartesian gradients from the chain rule through the
        // signed dihedral (validation F70). With a = r_i-r_j, b = r_k-r_j,
        // c = r_l-r_k, unit plane normals u = n1/m1, v = n2/m2,
        // C = u.v, S = sin(phi) > 0, s = sign(((n1 x n2).b)) and
        // z = v - C*u, w = u - C*v:
        //   dC/dr_i = -(z x b)/m1,  dC/dr_l = (w x b)/m2,
        //   dC/dr_j = (z x (b-a))/m1 + (w x c)/m2,
        //   dC/dr_k = (z x a)/m1 + ((b+c) x w)/m2,
        // and dphi/dr = -s*dC/dr/S. The sum over the four atoms vanishes
        // identically (translation invariance). The previous Blondel-Karplus
        // form missed the central-bond-length factor and used a lever-rule
        // approximation for the central atoms.
        fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
            [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ]
        }
        let sin_phi = (1.0 - cos_phi * cos_phi).sqrt().max(1e-8);
        let s_sign = if dot_dir < 0.0 { -1.0 } else { 1.0 };
        let u = [n1[0] / m1, n1[1] / m1, n1[2] / m1];
        let v = [n2[0] / m2, n2[1] / m2, n2[2] / m2];
        let z = [
            v[0] - cos_phi * u[0],
            v[1] - cos_phi * u[1],
            v[2] - cos_phi * u[2],
        ];
        let w = [
            u[0] - cos_phi * v[0],
            u[1] - cos_phi * v[1],
            u[2] - cos_phi * v[2],
        ];
        let b_minus_a = [b2[0] - b1[0], b2[1] - b1[1], b2[2] - b1[2]];
        let b_plus_c = [b2[0] + b3[0], b2[1] + b3[1], b2[2] + b3[2]];
        // -s/S * dC/dr for each atom, scaled by dE/dphi below.
        let pre = -s_sign / sin_phi;
        let gi = cross(&z, &b2);
        let gl = cross(&w, &b2);
        let zbma = cross(&z, &b_minus_a);
        let wc = cross(&w, &b3);
        let za = cross(&z, &b1);
        let bpcw = cross(&b_plus_c, &w);
        for dim in 0..3 {
            grad[i + dim] += d_e_dphi * pre * (-gi[dim] / m1);
            grad[l + dim] += d_e_dphi * pre * (gl[dim] / m2);
            grad[j + dim] += d_e_dphi * pre * (zbma[dim] / m1 + wc[dim] / m2);
            grad[k + dim] += d_e_dphi * pre * (za[dim] / m1 + bpcw[dim] / m2);
        }

        energy
    }
}

/// UFF inversion (out-of-plane) term for near-planar centres.
///
/// Mirrors the canonical `Inversion.cpp`: `E = K·(C0 + C1·sinY + C2·cos2W)`
/// where `sinY` comes from the cosine between the out-of-plane bond and the
/// plane normal, and three contributions share the centre force constant
/// (already divided by three in `get_uff_inversion_params`).
///
/// `K` is in kcal/mol/rad², the three coefficients are dimensionless, the
/// energy is in kcal/mol and the injected gradient in kcal/mol/Å.
pub struct UffInversion {
    /// Index of the first in-plane atom bonded to the central atom.
    pub idx_i: usize,
    /// Index of the central atom whose planarity is measured.
    pub idx_j: usize,
    /// Index of the second in-plane atom bonded to the central atom.
    pub idx_k: usize,
    /// Index of the out-of-plane atom.
    pub idx_l: usize,
    /// Inversion force constant `K` in kcal/mol/rad², already divided by three.
    pub k_inv: f64,
    /// Constant term `C0` (dimensionless).
    pub c0: f64,
    /// `sinY` coefficient `C1` (dimensionless).
    pub c1: f64,
    /// `cos2W` coefficient `C2` (dimensionless).
    pub c2: f64,
}

impl ForceFieldContribution for UffInversion {
    fn kind(&self) -> &'static str {
        "inversion"
    }
    fn evaluate_energy_and_inject_gradient(&self, coords: &[f64], grad: &mut [f64]) -> f64 {
        let pi = self.idx_i * 3;
        let pj = self.idx_j * 3;
        let pk = self.idx_k * 3;
        let pl = self.idx_l * 3;

        let point = |base: usize| [coords[base], coords[base + 1], coords[base + 2]];
        let sub = |u: [f64; 3], v: [f64; 3]| [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
        let cross = |u: [f64; 3], v: [f64; 3]| {
            [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ]
        };
        let dot = |u: [f64; 3], v: [f64; 3]| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
        let norm = |u: [f64; 3]| (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();

        let (pa, pb, pc, pd) = (point(pi), point(pj), point(pk), point(pl));
        let r_ji = sub(pa, pb);
        let r_jk = sub(pc, pb);
        let r_jl = sub(pd, pb);
        let d_ji = norm(r_ji);
        let d_jk = norm(r_jk);
        let d_jl = norm(r_jl);
        if d_ji < 1e-16 || d_jk < 1e-16 || d_jl < 1e-16 {
            return 0.0;
        }

        // RDKit convention: n = (-rJI) x rJK (normalized).
        let n = cross([-r_ji[0], -r_ji[1], -r_ji[2]], r_jk);
        let n_len = norm(n);
        if n_len < 1e-16 {
            return 0.0;
        }
        let n_hat = [n[0] / n_len, n[1] / n_len, n[2] / n_len];
        let mut cos_y = dot(n_hat, r_jl) / d_jl;
        cos_y = cos_y.clamp(-1.0, 1.0);
        let sin_y_sq = 1.0 - cos_y * cos_y;
        let sin_y = (sin_y_sq.max(0.0)).sqrt().max(1e-8);
        let cos2w = 2.0 * sin_y * sin_y - 1.0;
        let energy = self.k_inv * (self.c0 + self.c1 * sin_y + self.c2 * cos2w);

        // Gradient (port of InversionContrib::getGrad, normalized vectors).
        let r_ji_hat = [r_ji[0] / d_ji, r_ji[1] / d_ji, r_ji[2] / d_ji];
        let r_jk_hat = [r_jk[0] / d_jk, r_jk[1] / d_jk, r_jk[2] / d_jk];
        let r_jl_hat = [r_jl[0] / d_jl, r_jl[1] / d_jl, r_jl[2] / d_jl];
        let cos_theta = dot(r_ji_hat, r_jk_hat).clamp(-1.0, 1.0);
        let sin_theta_sq = 1.0 - cos_theta * cos_theta;
        let sin_theta = (sin_theta_sq.max(0.0)).sqrt().max(1e-8);

        let d_e_dw = -self.k_inv * (self.c1 * cos_y + 4.0 * self.c2 * cos_y * sin_y);
        let t1 = cross(r_jl_hat, r_jk_hat);
        let t2 = cross(r_ji_hat, r_jl_hat);
        let t3 = cross(r_jk_hat, r_ji_hat);
        let term1 = sin_y * sin_theta;
        let term2 = cos_y / (sin_y * sin_theta_sq);

        for dim in 0..3 {
            let tg1 =
                (t1[dim] / term1 - (r_ji_hat[dim] - r_jk_hat[dim] * cos_theta) * term2) / d_ji;
            let tg3 =
                (t2[dim] / term1 - (r_jk_hat[dim] - r_ji_hat[dim] * cos_theta) * term2) / d_jk;
            let tg4 = (t3[dim] / term1 - r_jl_hat[dim] * cos_y / sin_y) / d_jl;

            grad[pi + dim] += d_e_dw * tg1;
            grad[pj + dim] -= d_e_dw * (tg1 + tg3 + tg4);
            grad[pk + dim] += d_e_dw * tg3;
            grad[pl + dim] += d_e_dw * tg4;
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
    /// Index of the first atom of the pair.
    pub atom_i_idx: usize,
    /// Index of the second atom of the pair.
    pub atom_j_idx: usize,
    /// Combined Lennard-Jones distance `x_ij = (x_i + x_j)/2` in Å.
    pub r_star: f64,
    /// Combined well depth `ε_ij = √(D_i·D_j)` in kcal/mol
    /// (already scaled, e.g. by 0.5 for 1-4 pairs, before construction).
    pub epsilon: f64,
}

impl ForceFieldContribution for UffLennardJones {
    fn kind(&self) -> &'static str {
        "vdw"
    }
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
