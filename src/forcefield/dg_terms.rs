use super::traits::ForceFieldContribution;

/// Distance Violation contribution for Distance Geometry refinement.
/// E = weight * (d/ub - 1)^2 if d > ub, or weight * (2*lb/(lb+d) - 1)^2 if d < lb.
pub struct DistanceViolation {
    pub atom_i_idx: usize,
    pub atom_j_idx: usize,
    pub lower_bound: f64,
    pub upper_bound: f64,
    pub weight: f64,
}

impl ForceFieldContribution for DistanceViolation {
    fn evaluate_energy_and_inject_gradient(&self, coords: &[f64], grad: &mut [f64]) -> f64 {
        let root_i = self.atom_i_idx * 3;
        let root_j = self.atom_j_idx * 3;

        let dx = coords[root_i] - coords[root_j];
        let dy = coords[root_i + 1] - coords[root_j + 1];
        let dz = coords[root_i + 2] - coords[root_j + 2];

        let d2 = dx * dx + dy * dy + dz * dz;
        let d = d2.sqrt();

        let ub = self.upper_bound;
        let lb = self.lower_bound;
        let ub2 = ub * ub;
        let lb2 = lb * lb;

        if d2 > ub2 {
            let val = (d2 / ub2) - 1.0;
            let energy = self.weight * val * val;
            let pre_factor = 4.0 * self.weight * val / ub2;

            grad[root_i] += pre_factor * dx;
            grad[root_i + 1] += pre_factor * dy;
            grad[root_i + 2] += pre_factor * dz;

            grad[root_j] -= pre_factor * dx;
            grad[root_j + 1] -= pre_factor * dy;
            grad[root_j + 2] -= pre_factor * dz;

            energy
        } else if d2 < lb2 && d > 1e-8 {
            let l2d2 = lb2 + d2;
            let val = (2.0 * lb2 / l2d2) - 1.0;
            let energy = self.weight * val * val;

            // dE/dd = 2 * weight * (2*lb^2/(lb^2+d^2) - 1) * (-4*lb^2*d / (lb^2+d^2)^2)
            // pre_factor = dE/dd / d = -8 * weight * lb2 * val / (l2d2 * l2d2)
            // Note: RDKit uses a slightly different pre-factor for d < lb, but this is the analytical one.
            let pre_factor = -8.0 * self.weight * lb2 * val / (l2d2 * l2d2);

            grad[root_i] += pre_factor * dx;
            grad[root_i + 1] += pre_factor * dy;
            grad[root_i + 2] += pre_factor * dz;

            grad[root_j] -= pre_factor * dx;
            grad[root_j + 1] -= pre_factor * dy;
            grad[root_j + 2] -= pre_factor * dz;

            energy
        } else {
            0.0
        }
    }
}

/// Chiral Violation contribution for Distance Geometry refinement.
pub struct ChiralViolation {
    pub atom_idx: [usize; 4], // [i, j, k, l] -> l is center? No, usually i, j, k are neighbors, l is center
    pub lower_vol: f64,
    pub upper_vol: f64,
    pub weight: f64,
}

impl ForceFieldContribution for ChiralViolation {
    fn evaluate_energy_and_inject_gradient(&self, coords: &[f64], grad: &mut [f64]) -> f64 {
        let i = self.atom_idx[0] * 3;
        let j = self.atom_idx[1] * 3;
        let k = self.atom_idx[2] * 3;
        let l = self.atom_idx[3] * 3;

        let v1 = [
            coords[i] - coords[l],
            coords[i + 1] - coords[l + 1],
            coords[i + 2] - coords[l + 2],
        ];
        let v2 = [
            coords[j] - coords[l],
            coords[j + 1] - coords[l + 1],
            coords[j + 2] - coords[l + 2],
        ];
        let v3 = [
            coords[k] - coords[l],
            coords[k + 1] - coords[l + 1],
            coords[k + 2] - coords[l + 2],
        ];

        let v2xv3 = [
            v2[1] * v3[2] - v2[2] * v3[1],
            v2[2] * v3[0] - v2[0] * v3[2],
            v2[0] * v3[1] - v2[1] * v3[0],
        ];
        let vol = v1[0] * v2xv3[0] + v1[1] * v2xv3[1] + v1[2] * v2xv3[2];

        let pre_factor;
        let energy;
        if vol < self.lower_vol {
            let diff = vol - self.lower_vol;
            energy = self.weight * diff * diff;
            pre_factor = 2.0 * self.weight * diff;
        } else if vol > self.upper_vol {
            let diff = vol - self.upper_vol;
            energy = self.weight * diff * diff;
            pre_factor = 2.0 * self.weight * diff;
        } else {
            return 0.0;
        }

        // dV/dpos1 = v2 x v3
        grad[i] += pre_factor * v2xv3[0];
        grad[i + 1] += pre_factor * v2xv3[1];
        grad[i + 2] += pre_factor * v2xv3[2];

        // dV/dpos2 = v3 x v1
        let v3xv1 = [
            v3[1] * v1[2] - v3[2] * v1[1],
            v3[2] * v1[0] - v3[0] * v1[2],
            v3[0] * v1[1] - v3[1] * v1[0],
        ];
        grad[j] += pre_factor * v3xv1[0];
        grad[j + 1] += pre_factor * v3xv1[1];
        grad[j + 2] += pre_factor * v3xv1[2];

        // dV/dpos3 = v1 x v2
        let v1xv2 = [
            v1[1] * v2[2] - v1[2] * v2[1],
            v1[2] * v2[0] - v1[0] * v2[2],
            v1[0] * v2[1] - v1[1] * v2[0],
        ];
        grad[k] += pre_factor * v1xv2[0];
        grad[k + 1] += pre_factor * v1xv2[1];
        grad[k + 2] += pre_factor * v1xv2[2];

        // dV/dpos4 = -(dV/dpos1 + dV/dpos2 + dV/dpos3)
        grad[l] -= pre_factor * (v2xv3[0] + v3xv1[0] + v1xv2[0]);
        grad[l + 1] -= pre_factor * (v2xv3[1] + v3xv1[1] + v1xv2[1]);
        grad[l + 2] -= pre_factor * (v2xv3[2] + v3xv1[2] + v1xv2[2]);

        energy
    }
}
