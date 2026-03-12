//! Energy and gradient functions for bounds violation, chiral violation, and embed torsions.

use nalgebra::{DMatrix, Vector3};
use super::*;


/// Computes the bounds-violation energy (RDKit's DistViolationContribs).
///
/// For each atom pair (i, j) with lower bound `lb` and upper bound `ub`:
/// - If d² > ub²: val = (d²/ub²) - 1.0; energy += weight * val²
/// - If d² < lb²: val = (2*lb²/(lb²+d²)) - 1.0; energy += weight * val²
/// - Otherwise: no penalty
/// Computes the bounds-violation energy (RDKit's DistViolationContribs).
pub fn bounds_violation_energy(coords: &DMatrix<f32>, bounds: &DMatrix<f64>) -> f32 {
    bounds_violation_energy_basin(coords, bounds, 1000.0)
}

/// Bounds violation energy with basin threshold filtering.
/// Only pairs with (ub - lb) <= basin_thresh are included.
/// RDKit uses basin_thresh=5.0 during embedding minimization.
pub fn bounds_violation_energy_basin(coords: &DMatrix<f32>, bounds: &DMatrix<f64>, basin_thresh: f32) -> f32 {
    let n = coords.nrows();
    let dim = coords.ncols();
    let mut energy = 0.0f32;

    for i in 1..n {
        for j in 0..i {
            let ub = bounds[(j, i)] as f32; // upper triangle = upper bound
            let lb = bounds[(i, j)] as f32; // lower triangle = lower bound

            // Basin threshold: skip wide-bounded pairs
            if ub - lb > basin_thresh {
                continue;
            }

            let mut d2 = 0.0f32;
            for d in 0..dim {
                let diff = coords[(i, d)] - coords[(j, d)];
                d2 += diff * diff;
            }

            let ub2 = ub * ub;
            let lb2 = lb * lb;

            let val = if d2 > ub2 {
                (d2 / ub2) - 1.0
            } else if d2 < lb2 {
                (2.0 * lb2 / (lb2 + d2)) - 1.0
            } else {
                0.0
            };

            if val > 0.0 {
                energy += val * val; // weight = 1.0
            }
        }
    }
    energy
}

/// Computes the analytical gradient of the bounds-violation energy.
///
/// Matches RDKit's DistViolationContribs::getGrad exactly.
pub fn bounds_violation_gradient(coords: &DMatrix<f32>, bounds: &DMatrix<f64>) -> DMatrix<f32> {
    bounds_violation_gradient_basin(coords, bounds, 1000.0)
}

/// Gradient with basin threshold filtering.
pub fn bounds_violation_gradient_basin(coords: &DMatrix<f32>, bounds: &DMatrix<f64>, basin_thresh: f32) -> DMatrix<f32> {
    let n = coords.nrows();
    let dim = coords.ncols();
    let mut grad = DMatrix::from_element(n, dim, 0.0f32);

    for i in 1..n {
        for j in 0..i {
            let ub = bounds[(j, i)] as f32;
            let lb = bounds[(i, j)] as f32;

            // Basin threshold: skip wide-bounded pairs
            if ub - lb > basin_thresh {
                continue;
            }

            let mut d2 = 0.0f32;
            let mut diffs = Vec::with_capacity(dim);
            for d_idx in 0..dim {
                let diff = coords[(i, d_idx)] - coords[(j, d_idx)];
                diffs.push(diff);
                d2 += diff * diff;
            }

            let ub2 = ub * ub;
            let lb2 = lb * lb;

            if d2 > ub2 {
                let d = d2.sqrt();
                if d < 1e-8 {
                    continue;
                }
                // RDKit: 2 * (d^2/ub^2 - 1) * (2d/ub^2) = 4d * (d^2/ub^2 - 1) / ub^2
                // dE/dd = 2 * ((d^2/ub^2) - 1) * (2d/ub^2)
                // dE/dx = dE/dd * dd/dx = dE/dd * (x_i - x_j)/d
                // pre_factor = dE/dd / d = 2 * ((d^2/ub^2) - 1) * (2/ub^2)
                let pre_factor = 4.0 * ((d2 / ub2) - 1.0) / ub2;
                for d_idx in 0..dim {
                    let g = pre_factor * diffs[d_idx];
                    grad[(i, d_idx)] += g;
                    grad[(j, d_idx)] -= g;
                }
            } else if d2 < lb2 {
                let d = d2.sqrt();
                if d < 1e-8 {
                    continue;
                }
                // RDKit: 2 * (2*lb^2/(lb^2+d^2) - 1) * (2*lb^2 * (-2d) / (lb^2+d^2)^2)
                // dE/dd = 2 * (2*lb^2/(lb^2+d^2) - 1) * (-4*lb^2*d / (lb^2+d^2)^2)
                // pre_factor = dE/dd / d = 2 * (2*lb^2/(lb^2+d^2) - 1) * (-4*lb^2 / (lb^2+d^2)^2)
                let l2d2 = lb2 + d2;
                let pre_factor = 8.0 * lb2 * (1.0 - 2.0 * lb2 / l2d2) / (l2d2 * l2d2);
                for d_idx in 0..dim {
                    let g = pre_factor * diffs[d_idx];
                    grad[(i, d_idx)] += g;
                    grad[(j, d_idx)] -= g;
                }
            }
        }
    }
    grad
}

/// L-BFGS minimization using bounds-violation energy/gradient only.
/// This matches RDKit's firstMinimization step.
/// basin_thresh: skip pairs with (ub - lb) > basin_thresh (RDKit uses 5.0)
/// weight_4d: penalty weight for 4th dimension (RDKit uses 0.1 for stage1, 1.0 for stage2)
pub fn chiral_violation_energy(coords: &DMatrix<f32>, chiral_sets: &[ChiralSet]) -> f32 {
    let mut energy = 0.0f32;
    for c in chiral_sets {
        let vol = crate::distgeom::calc_chiral_volume(
            c.neighbors[0],
            c.neighbors[1],
            c.neighbors[2],
            c.neighbors[3],
            coords,
        );
        if vol < c.lower_vol {
            let diff = vol - c.lower_vol;
            energy += diff * diff;
        } else if vol > c.upper_vol {
            let diff = vol - c.upper_vol;
            energy += diff * diff;
        }
    }
    energy
}

/// Chiral violation energy in f64, matching RDKit's f64 (Point3D) precision.
pub fn chiral_violation_energy_f64(coords_flat: &[f64], dim: usize, chiral_sets: &[ChiralSet]) -> f64 {
    let mut energy = 0.0f64;
    for c in chiral_sets {
        let (i1, i2, i3, i4) = (c.neighbors[0], c.neighbors[1], c.neighbors[2], c.neighbors[3]);
        let v1 = [coords_flat[i1*dim] - coords_flat[i4*dim],
                   coords_flat[i1*dim+1] - coords_flat[i4*dim+1],
                   coords_flat[i1*dim+2] - coords_flat[i4*dim+2]];
        let v2 = [coords_flat[i2*dim] - coords_flat[i4*dim],
                   coords_flat[i2*dim+1] - coords_flat[i4*dim+1],
                   coords_flat[i2*dim+2] - coords_flat[i4*dim+2]];
        let v3 = [coords_flat[i3*dim] - coords_flat[i4*dim],
                   coords_flat[i3*dim+1] - coords_flat[i4*dim+1],
                   coords_flat[i3*dim+2] - coords_flat[i4*dim+2]];
        // v2 x v3
        let cx = v2[1]*v3[2] - v2[2]*v3[1];
        let cy = v2[2]*v3[0] - v2[0]*v3[2];
        let cz = v2[0]*v3[1] - v2[1]*v3[0];
        let vol = v1[0]*cx + v1[1]*cy + v1[2]*cz;
        let lb = c.lower_vol as f64;
        let ub = c.upper_vol as f64;
        if vol < lb {
            let diff = vol - lb;
            energy += diff * diff;
        } else if vol > ub {
            let diff = vol - ub;
            energy += diff * diff;
        }
    }
    energy
}

/// Chiral violation gradient in f64, matching RDKit's f64 precision.
/// Adds weight * gradient to the output array.
pub fn chiral_violation_gradient_f64(coords_flat: &[f64], dim: usize, chiral_sets: &[ChiralSet], weight: f64, grad: &mut [f64]) {
    for c in chiral_sets {
        let (i1, i2, i3, i4) = (c.neighbors[0], c.neighbors[1], c.neighbors[2], c.neighbors[3]);
        let v1 = [coords_flat[i1*dim] - coords_flat[i4*dim],
                   coords_flat[i1*dim+1] - coords_flat[i4*dim+1],
                   coords_flat[i1*dim+2] - coords_flat[i4*dim+2]];
        let v2 = [coords_flat[i2*dim] - coords_flat[i4*dim],
                   coords_flat[i2*dim+1] - coords_flat[i4*dim+1],
                   coords_flat[i2*dim+2] - coords_flat[i4*dim+2]];
        let v3 = [coords_flat[i3*dim] - coords_flat[i4*dim],
                   coords_flat[i3*dim+1] - coords_flat[i4*dim+1],
                   coords_flat[i3*dim+2] - coords_flat[i4*dim+2]];
        // v2 x v3
        let v2xv3 = [v2[1]*v3[2] - v2[2]*v3[1],
                      v2[2]*v3[0] - v2[0]*v3[2],
                      v2[0]*v3[1] - v2[1]*v3[0]];
        let vol = v1[0]*v2xv3[0] + v1[1]*v2xv3[1] + v1[2]*v2xv3[2];
        let lb = c.lower_vol as f64;
        let ub = c.upper_vol as f64;
        let pre_factor = if vol < lb {
            weight * (vol - lb)
        } else if vol > ub {
            weight * (vol - ub)
        } else {
            continue;
        };
        // dV/dpos1 = v2 x v3
        grad[i1*dim]   += pre_factor * v2xv3[0];
        grad[i1*dim+1] += pre_factor * v2xv3[1];
        grad[i1*dim+2] += pre_factor * v2xv3[2];
        // dV/dpos2 = v3 x v1
        let v3xv1 = [v3[1]*v1[2] - v3[2]*v1[1],
                      v3[2]*v1[0] - v3[0]*v1[2],
                      v3[0]*v1[1] - v3[1]*v1[0]];
        grad[i2*dim]   += pre_factor * v3xv1[0];
        grad[i2*dim+1] += pre_factor * v3xv1[1];
        grad[i2*dim+2] += pre_factor * v3xv1[2];
        // dV/dpos3 = v1 x v2
        let v1xv2 = [v1[1]*v2[2] - v1[2]*v2[1],
                      v1[2]*v2[0] - v1[0]*v2[2],
                      v1[0]*v2[1] - v1[1]*v2[0]];
        grad[i3*dim]   += pre_factor * v1xv2[0];
        grad[i3*dim+1] += pre_factor * v1xv2[1];
        grad[i3*dim+2] += pre_factor * v1xv2[2];
        // dV/dpos4 = -(sum)
        grad[i4*dim]   -= pre_factor * (v2xv3[0] + v3xv1[0] + v1xv2[0]);
        grad[i4*dim+1] -= pre_factor * (v2xv3[1] + v3xv1[1] + v1xv2[1]);
        grad[i4*dim+2] -= pre_factor * (v2xv3[2] + v3xv1[2] + v1xv2[2]);
    }
}

pub fn chiral_violation_gradient(
    coords: &DMatrix<f32>,
    chiral_sets: &[ChiralSet],
    grad: &mut DMatrix<f32>,
) {
    for c in chiral_sets {
        let (idx1, idx2, idx3, idx4) = (
            c.neighbors[0],
            c.neighbors[1],
            c.neighbors[2],
            c.neighbors[3],
        );

        // v1 = pos1 - pos4, v2 = pos2 - pos4, v3 = pos3 - pos4
        let v1 = Vector3::new(
            coords[(idx1, 0)] - coords[(idx4, 0)],
            coords[(idx1, 1)] - coords[(idx4, 1)],
            coords[(idx1, 2)] - coords[(idx4, 2)],
        );
        let v2 = Vector3::new(
            coords[(idx2, 0)] - coords[(idx4, 0)],
            coords[(idx2, 1)] - coords[(idx4, 1)],
            coords[(idx2, 2)] - coords[(idx4, 2)],
        );
        let v3 = Vector3::new(
            coords[(idx3, 0)] - coords[(idx4, 0)],
            coords[(idx3, 1)] - coords[(idx4, 1)],
            coords[(idx3, 2)] - coords[(idx4, 2)],
        );

        let v2xv3 = v2.cross(&v3);
        let vol = v1.dot(&v2xv3);

        let pre_factor;
        if vol < c.lower_vol {
            pre_factor = vol - c.lower_vol;
        } else if vol > c.upper_vol {
            pre_factor = vol - c.upper_vol;
        } else {
            continue;
        }

        // dV/dpos1 = v2 x v3
        grad[(idx1, 0)] += pre_factor * v2xv3.x;
        grad[(idx1, 1)] += pre_factor * v2xv3.y;
        grad[(idx1, 2)] += pre_factor * v2xv3.z;

        // dV/dpos2 = v3 x v1
        let v3xv1 = v3.cross(&v1);
        grad[(idx2, 0)] += pre_factor * v3xv1.x;
        grad[(idx2, 1)] += pre_factor * v3xv1.y;
        grad[(idx2, 2)] += pre_factor * v3xv1.z;

        // dV/dpos3 = v1 x v2
        let v1xv2 = v1.cross(&v2);
        grad[(idx3, 0)] += pre_factor * v1xv2.x;
        grad[(idx3, 1)] += pre_factor * v1xv2.y;
        grad[(idx3, 2)] += pre_factor * v1xv2.z;

        // dV/dpos4 = -(dV/dpos1 + dV/dpos2 + dV/dpos3)
        grad[(idx4, 0)] -= pre_factor * (v2xv3.x + v3xv1.x + v1xv2.x);
        grad[(idx4, 1)] -= pre_factor * (v2xv3.y + v3xv1.y + v1xv2.y);
        grad[(idx4, 2)] -= pre_factor * (v2xv3.z + v3xv1.z + v1xv2.z);
    }
}

/// A pre-computed torsion term for use during 4D embedding minimization.
pub fn collect_embed_torsions(mol: &crate::graph::Molecule) -> Vec<EmbedTorsion> {
    use petgraph::visit::EdgeRef;
    let mut terms = Vec::new();

    for edge in mol.graph.edge_references() {
        let u = edge.source();
        let v = edge.target();
        let hyb_u = mol.graph[u].hybridization;
        let hyb_v = mol.graph[v].hybridization;

        if hyb_u == crate::graph::Hybridization::SP
            || hyb_v == crate::graph::Hybridization::SP
        {
            continue;
        }

        let neighbors_u: Vec<_> = mol.graph.neighbors(u).filter(|&x| x != v).collect();
        let neighbors_v: Vec<_> = mol.graph.neighbors(v).filter(|&x| x != u).collect();
        if neighbors_u.is_empty() || neighbors_v.is_empty() {
            continue;
        }

        let is_ring = crate::graph::min_path_excluding2(mol, u, v, u, v, 7).is_some();

        // UFF-style torsion for all bonds with non-SP atoms
        let (n_fold, gamma, weight) = crate::forcefield::energy::torsion_params(hyb_u, hyb_v);
        for &nu in &neighbors_u {
            for &nv in &neighbors_v {
                terms.push(EmbedTorsion {
                    idx: [nu.index(), u.index(), v.index(), nv.index()],
                    n_fold,
                    gamma,
                    weight: weight * 10.0,
                });
            }
        }

        // ETKDG-lite M6 for non-ring rotatable bonds
        if !is_ring {
            let m6 = crate::forcefield::etkdg_lite::infer_etkdg_parameters(mol, u.index(), v.index());
            let nu = neighbors_u[0];
            let nv = neighbors_v[0];
            for k in 0..6 {
                if m6.v[k].abs() > 1e-6 {
                    let nf = (k + 1) as f32;
                    let gam = if m6.s[k] > 0.0 { 0.0 } else { std::f32::consts::PI / nf };
                    terms.push(EmbedTorsion {
                        idx: [nu.index(), u.index(), v.index(), nv.index()],
                        n_fold: nf,
                        gamma: gam,
                        weight: m6.v[k],
                    });
                }
            }
        }
    }
    terms
}

pub(crate) fn torsion_energy_4d(coords: &DMatrix<f32>, terms: &[EmbedTorsion]) -> f32 {
    let mut energy = 0.0f32;
    for t in terms {
        let [i1, i2, i3, i4] = t.idx;
        let p1 = Vector3::new(coords[(i1, 0)], coords[(i1, 1)], coords[(i1, 2)]);
        let p2 = Vector3::new(coords[(i2, 0)], coords[(i2, 1)], coords[(i2, 2)]);
        let p3 = Vector3::new(coords[(i3, 0)], coords[(i3, 1)], coords[(i3, 2)]);
        let p4 = Vector3::new(coords[(i4, 0)], coords[(i4, 1)], coords[(i4, 2)]);
        energy += crate::forcefield::energy::torsional_energy(&p1, &p2, &p3, &p4, t.weight, t.n_fold, t.gamma);
    }
    energy
}

pub(crate) fn torsion_gradient_4d(coords: &DMatrix<f32>, terms: &[EmbedTorsion], grad: &mut DMatrix<f32>) {
    let n = coords.nrows();
    let mut grad3 = DMatrix::zeros(n, 3);
    for t in terms {
        let [i1, i2, i3, i4] = t.idx;
        let p1 = Vector3::new(coords[(i1, 0)], coords[(i1, 1)], coords[(i1, 2)]);
        let p2 = Vector3::new(coords[(i2, 0)], coords[(i2, 1)], coords[(i2, 2)]);
        let p3 = Vector3::new(coords[(i3, 0)], coords[(i3, 1)], coords[(i3, 2)]);
        let p4 = Vector3::new(coords[(i4, 0)], coords[(i4, 1)], coords[(i4, 2)]);
        crate::forcefield::gradients::analytical_grad_torsion(
            &p1, &p2, &p3, &p4, t.weight, t.n_fold, t.gamma,
            &mut grad3, i1, i2, i3, i4,
        );
    }
    for i in 0..n {
        for d in 0..3 {
            grad[(i, d)] += grad3[(i, d)];
        }
    }
}
