//! Gradient functions for the ETKDG 3D force field (f32 and f64 versions).

use nalgebra::{DMatrix, Vector3};
use petgraph::visit::EdgeRef;
use super::*;


/// Calculate gradient using the 3D ETKDG force field
pub fn etkdg_3d_gradient(
    coords: &DMatrix<f32>,
    mol: &crate::graph::Molecule,
    ff: &Etkdg3DFF,
) -> DMatrix<f32> {
    let n = mol.graph.node_count();
    let mut grad = DMatrix::zeros(n, 3);

    // Distance constraints gradient
    for c in ff.dist_12.iter().chain(ff.dist_13.iter()).chain(ff.dist_long.iter()) {
        let p1 = Vector3::new(coords[(c.i, 0)], coords[(c.i, 1)], coords[(c.i, 2)]);
        let p2 = Vector3::new(coords[(c.j, 0)], coords[(c.j, 1)], coords[(c.j, 2)]);
        crate::forcefield::gradients::analytical_grad_distance_constraint(
            &p1, &p2, c.min_len as f32, c.max_len as f32, c.k as f32, &mut grad, c.i, c.j,
        );
    }

    // Angle constraints gradient
    for ac in &ff.angle_constraints {
        super::energy::angle_constraint_gradient(coords, ac, &mut grad);
    }

    // OOP gradient (simple vol² formulation)
    if ff.oop_k.abs() > 1e-8 {
        for i in 0..n {
            let ni = petgraph::graph::NodeIndex::new(i);
            if mol.graph[ni].hybridization != crate::graph::Hybridization::SP2 { continue; }
            let nbs: Vec<_> = mol.graph.neighbors(ni).collect();
            if nbs.len() != 3 { continue; }
            let pc = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
            let p1 = Vector3::new(coords[(nbs[0].index(), 0)], coords[(nbs[0].index(), 1)], coords[(nbs[0].index(), 2)]);
            let p2 = Vector3::new(coords[(nbs[1].index(), 0)], coords[(nbs[1].index(), 1)], coords[(nbs[1].index(), 2)]);
            let p3 = Vector3::new(coords[(nbs[2].index(), 0)], coords[(nbs[2].index(), 1)], coords[(nbs[2].index(), 2)]);
            crate::forcefield::gradients::analytical_grad_oop(&pc, &p1, &p2, &p3, ff.oop_k as f32,
                &mut grad, i, nbs[0].index(), nbs[1].index(), nbs[2].index());
        }
    }

    // Pre-computed torsion gradients (flat ring + chain preferences)
    for tc in &ff.torsion_contribs {
        let p1 = Vector3::new(coords[(tc.i, 0)], coords[(tc.i, 1)], coords[(tc.i, 2)]);
        let p2 = Vector3::new(coords[(tc.j, 0)], coords[(tc.j, 1)], coords[(tc.j, 2)]);
        let p3 = Vector3::new(coords[(tc.k, 0)], coords[(tc.k, 1)], coords[(tc.k, 2)]);
        let p4 = Vector3::new(coords[(tc.l, 0)], coords[(tc.l, 1)], coords[(tc.l, 2)]);
        let m6 = crate::forcefield::etkdg_lite::M6Params {
            s: tc.signs.map(|x| x as f32),
            v: tc.v.map(|x| x as f32),
        };
        crate::forcefield::etkdg_lite::calc_torsion_grad_m6(
            &p1, &p2, &p3, &p4, &m6,
            &mut grad, tc.i, tc.j, tc.k, tc.l,
        );
    }

    // UFF torsion gradient
    if ff.torsion_k_omega.abs() > 1e-8 && n >= 4 {
        for edge in mol.graph.edge_references() {
            let u = edge.source();
            let v = edge.target();
            let hyb_u = mol.graph[u].hybridization;
            let hyb_v = mol.graph[v].hybridization;
            if hyb_u == crate::graph::Hybridization::SP || hyb_v == crate::graph::Hybridization::SP { continue; }
            let (n_fold, gamma, weight) = crate::forcefield::energy::torsion_params(hyb_u, hyb_v);
            let neighbors_u: Vec<_> = mol.graph.neighbors(u).filter(|&x| x != v).collect();
            let neighbors_v: Vec<_> = mol.graph.neighbors(v).filter(|&x| x != u).collect();
            for &nu in &neighbors_u {
                for &nv in &neighbors_v {
                    let p1 = Vector3::new(coords[(nu.index(), 0)], coords[(nu.index(), 1)], coords[(nu.index(), 2)]);
                    let p2 = Vector3::new(coords[(u.index(), 0)], coords[(u.index(), 1)], coords[(u.index(), 2)]);
                    let p3 = Vector3::new(coords[(v.index(), 0)], coords[(v.index(), 1)], coords[(v.index(), 2)]);
                    let p4 = Vector3::new(coords[(nv.index(), 0)], coords[(nv.index(), 1)], coords[(nv.index(), 2)]);
                    crate::forcefield::gradients::analytical_grad_torsion(
                        &p1, &p2, &p3, &p4, ff.torsion_k_omega as f32 * weight, n_fold, gamma,
                        &mut grad, nu.index(), u.index(), v.index(), nv.index(),
                    );
                }
            }
        }
    }

    // ETKDG-lite M6 torsion gradient
    if ff.use_m6_torsions && n >= 4 {
        for edge in mol.graph.edge_references() {
            let u = edge.source();
            let v = edge.target();
            if crate::graph::min_path_excluding2(mol, u, v, u, v, 7).is_some() { continue; }
            let m6 = crate::forcefield::etkdg_lite::infer_etkdg_parameters(mol, u.index(), v.index());
            if m6.v.iter().all(|&x| x.abs() < 1e-6) { continue; }
            let neighbors_u: Vec<_> = mol.graph.neighbors(u).filter(|&x| x != v).collect();
            let neighbors_v: Vec<_> = mol.graph.neighbors(v).filter(|&x| x != u).collect();
            if neighbors_u.is_empty() || neighbors_v.is_empty() { continue; }
            let nu = neighbors_u[0];
            let nv = neighbors_v[0];
            let p1 = Vector3::new(coords[(nu.index(), 0)], coords[(nu.index(), 1)], coords[(nu.index(), 2)]);
            let p2 = Vector3::new(coords[(u.index(), 0)], coords[(u.index(), 1)], coords[(u.index(), 2)]);
            let p3 = Vector3::new(coords[(v.index(), 0)], coords[(v.index(), 1)], coords[(v.index(), 2)]);
            let p4 = Vector3::new(coords[(nv.index(), 0)], coords[(nv.index(), 1)], coords[(nv.index(), 2)]);
            crate::forcefield::etkdg_lite::calc_torsion_grad_m6(
                &p1, &p2, &p3, &p4, &m6,
                &mut grad, nu.index(), u.index(), v.index(), nv.index(),
            );
        }
    }

    grad
}

/// L-BFGS minimizer for the 3D ETKDG force field

/// All computation done in f64 to match RDKit's double-precision force field.
pub fn etkdg_3d_gradient_f64(
    coords: &[f64],
    n: usize,
    mol: &crate::graph::Molecule,
    ff: &Etkdg3DFF,
) -> Vec<f64> {
    let c = |atom: usize, d: usize| -> f64 { coords[atom * 3 + d] };
    let dim = n * 3;
    let mut grad = vec![0.0f64; dim];

    // === RDKit accumulation order: torsions → inversions → dist_12 → angles → dist_13 → dist_long ===

    // Pre-computed M6 torsion gradient (contribs[0])
    for tc in &ff.torsion_contribs {
        let r0 = [c(tc.i,0)-c(tc.j,0), c(tc.i,1)-c(tc.j,1), c(tc.i,2)-c(tc.j,2)];
        let r1 = [c(tc.k,0)-c(tc.j,0), c(tc.k,1)-c(tc.j,1), c(tc.k,2)-c(tc.j,2)];
        let r2 = [-r1[0], -r1[1], -r1[2]];
        let r3 = [c(tc.l,0)-c(tc.k,0), c(tc.l,1)-c(tc.k,1), c(tc.l,2)-c(tc.k,2)];
        let t0 = [r0[1]*r1[2]-r0[2]*r1[1], r0[2]*r1[0]-r0[0]*r1[2], r0[0]*r1[1]-r0[1]*r1[0]];
        let d0 = (t0[0]*t0[0]+t0[1]*t0[1]+t0[2]*t0[2]).sqrt();
        let t1 = [r2[1]*r3[2]-r2[2]*r3[1], r2[2]*r3[0]-r2[0]*r3[2], r2[0]*r3[1]-r2[1]*r3[0]];
        let d1 = (t1[0]*t1[0]+t1[1]*t1[1]+t1[2]*t1[2]).sqrt();
        if d0 < 1e-10 || d1 < 1e-10 { continue; }
        let t0 = [t0[0]/d0, t0[1]/d0, t0[2]/d0];
        let t1 = [t1[0]/d1, t1[1]/d1, t1[2]/d1];
        let cos_phi = (t0[0]*t1[0]+t0[1]*t1[1]+t0[2]*t1[2]).clamp(-1.0, 1.0);
        let sin_phi_sq = 1.0 - cos_phi * cos_phi;
        let sin_phi = if sin_phi_sq > 0.0 { sin_phi_sq.sqrt() } else { 0.0 };
        let cp2 = cos_phi*cos_phi;
        let cp3 = cos_phi*cp2;
        let cp4 = cos_phi*cp3;
        let cp5 = cos_phi*cp4;
        let vv = &tc.v;
        let ss = &tc.signs;
        let de_dphi =
            -vv[0] * ss[0] * sin_phi
            - 2.0 * vv[1] * ss[1] * (2.0 * cos_phi * sin_phi)
            - 3.0 * vv[2] * ss[2] * (4.0 * cp2 * sin_phi - sin_phi)
            - 4.0 * vv[3] * ss[3] * (8.0 * cp3 * sin_phi - 4.0 * cos_phi * sin_phi)
            - 5.0 * vv[4] * ss[4] * (16.0 * cp4 * sin_phi - 12.0 * cp2 * sin_phi + sin_phi)
            - 6.0 * vv[4] * ss[4] * (32.0 * cp5 * sin_phi - 32.0 * cp3 * sin_phi + 6.0 * cos_phi * sin_phi);
        let is_zero_sin = sin_phi < 1e-10;
        let sin_term = -de_dphi * if is_zero_sin { 1.0 / cos_phi } else { 1.0 / sin_phi };
        let dct = [
            (t1[0] - cos_phi * t0[0]) / d0,
            (t1[1] - cos_phi * t0[1]) / d0,
            (t1[2] - cos_phi * t0[2]) / d0,
            (t0[0] - cos_phi * t1[0]) / d1,
            (t0[1] - cos_phi * t1[1]) / d1,
            (t0[2] - cos_phi * t1[2]) / d1,
        ];
        let g1 = tc.i * 3;
        let g2 = tc.j * 3;
        let g3 = tc.k * 3;
        let g4 = tc.l * 3;
        grad[g1]   += sin_term * (dct[2]*r1[1] - dct[1]*r1[2]);
        grad[g1+1] += sin_term * (dct[0]*r1[2] - dct[2]*r1[0]);
        grad[g1+2] += sin_term * (dct[1]*r1[0] - dct[0]*r1[1]);
        grad[g2]   += sin_term * (dct[1]*(r1[2]-r0[2]) + dct[2]*(r0[1]-r1[1]) + dct[4]*(-r3[2]) + dct[5]*r3[1]);
        grad[g2+1] += sin_term * (dct[0]*(r0[2]-r1[2]) + dct[2]*(r1[0]-r0[0]) + dct[3]*r3[2] + dct[5]*(-r3[0]));
        grad[g2+2] += sin_term * (dct[0]*(r1[1]-r0[1]) + dct[1]*(r0[0]-r1[0]) + dct[3]*(-r3[1]) + dct[4]*r3[0]);
        grad[g3]   += sin_term * (dct[1]*r0[2] + dct[2]*(-r0[1]) + dct[4]*(r3[2]-r2[2]) + dct[5]*(r2[1]-r3[1]));
        grad[g3+1] += sin_term * (dct[0]*(-r0[2]) + dct[2]*r0[0] + dct[3]*(r2[2]-r3[2]) + dct[5]*(r3[0]-r2[0]));
        grad[g3+2] += sin_term * (dct[0]*r0[1] + dct[1]*(-r0[0]) + dct[3]*(r3[1]-r2[1]) + dct[4]*(r2[0]-r3[0]));
        grad[g4]   += sin_term * (dct[4]*r2[2] - dct[5]*r2[1]);
        grad[g4+1] += sin_term * (dct[5]*r2[0] - dct[3]*r2[2]);
        grad[g4+2] += sin_term * (dct[3]*r2[1] - dct[4]*r2[0]);
    }

    // UFF Inversion gradient (contribs[1])
    for ic in &ff.inversion_contribs {
        let p1 = [c(ic.at1,0), c(ic.at1,1), c(ic.at1,2)];
        let p2 = [c(ic.at2,0), c(ic.at2,1), c(ic.at2,2)];
        let p3 = [c(ic.at3,0), c(ic.at3,1), c(ic.at3,2)];
        let p4 = [c(ic.at4,0), c(ic.at4,1), c(ic.at4,2)];
        let mut rji = [p1[0]-p2[0], p1[1]-p2[1], p1[2]-p2[2]];
        let mut rjk = [p3[0]-p2[0], p3[1]-p2[1], p3[2]-p2[2]];
        let mut rjl = [p4[0]-p2[0], p4[1]-p2[1], p4[2]-p2[2]];
        let dji = (rji[0]*rji[0]+rji[1]*rji[1]+rji[2]*rji[2]).sqrt();
        let djk = (rjk[0]*rjk[0]+rjk[1]*rjk[1]+rjk[2]*rjk[2]).sqrt();
        let djl = (rjl[0]*rjl[0]+rjl[1]*rjl[1]+rjl[2]*rjl[2]).sqrt();
        if dji < 1e-8 || djk < 1e-8 || djl < 1e-8 { continue; }
        rji[0] /= dji; rji[1] /= dji; rji[2] /= dji;
        rjk[0] /= djk; rjk[1] /= djk; rjk[2] /= djk;
        rjl[0] /= djl; rjl[1] /= djl; rjl[2] /= djl;
        let mut nx = -rji[1]*rjk[2]+rji[2]*rjk[1];
        let mut ny = -rji[2]*rjk[0]+rji[0]*rjk[2];
        let mut nz = -rji[0]*rjk[1]+rji[1]*rjk[0];
        let nl = (nx*nx+ny*ny+nz*nz).sqrt();
        if nl < 1e-8 { continue; }
        nx /= nl; ny /= nl; nz /= nl;
        let cos_y = (nx*rjl[0]+ny*rjl[1]+nz*rjl[2]).clamp(-1.0, 1.0);
        let sin_y_sq = 1.0 - cos_y * cos_y;
        let sin_y = sin_y_sq.sqrt().max(1e-8);
        let cos_theta = (rji[0]*rjk[0]+rji[1]*rjk[1]+rji[2]*rjk[2]).clamp(-1.0, 1.0);
        let sin_theta_sq = 1.0 - cos_theta * cos_theta;
        let sin_theta = sin_theta_sq.sqrt().max(1e-8);
        let de_dw = -ic.force_constant * (ic.c1 * cos_y - 4.0 * ic.c2 * cos_y * sin_y);
        let t1 = [rjl[1]*rjk[2]-rjl[2]*rjk[1], rjl[2]*rjk[0]-rjl[0]*rjk[2], rjl[0]*rjk[1]-rjl[1]*rjk[0]];
        let t2 = [rji[1]*rjl[2]-rji[2]*rjl[1], rji[2]*rjl[0]-rji[0]*rjl[2], rji[0]*rjl[1]-rji[1]*rjl[0]];
        let t3 = [rjk[1]*rji[2]-rjk[2]*rji[1], rjk[2]*rji[0]-rjk[0]*rji[2], rjk[0]*rji[1]-rjk[1]*rji[0]];
        let term1 = sin_y * sin_theta;
        let term2 = cos_y / (sin_y * sin_theta_sq);
        let tg1 = [
            (t1[0] / term1 - (rji[0] - rjk[0] * cos_theta) * term2) / dji,
            (t1[1] / term1 - (rji[1] - rjk[1] * cos_theta) * term2) / dji,
            (t1[2] / term1 - (rji[2] - rjk[2] * cos_theta) * term2) / dji,
        ];
        let tg3 = [
            (t2[0] / term1 - (rjk[0] - rji[0] * cos_theta) * term2) / djk,
            (t2[1] / term1 - (rjk[1] - rji[1] * cos_theta) * term2) / djk,
            (t2[2] / term1 - (rjk[2] - rji[2] * cos_theta) * term2) / djk,
        ];
        let tg4 = [
            (t3[0] / term1 - rjl[0] * cos_y / sin_y) / djl,
            (t3[1] / term1 - rjl[1] * cos_y / sin_y) / djl,
            (t3[2] / term1 - rjl[2] * cos_y / sin_y) / djl,
        ];
        let g1 = ic.at1 * 3;
        let g2 = ic.at2 * 3;
        let g3 = ic.at3 * 3;
        let g4 = ic.at4 * 3;
        for dd in 0..3 {
            grad[g1+dd] += de_dw * tg1[dd];
            grad[g2+dd] += -de_dw * (tg1[dd] + tg3[dd] + tg4[dd]);
            grad[g3+dd] += de_dw * tg3[dd];
            grad[g4+dd] += de_dw * tg4[dd];
        }
    }

    // Distance constraints gradient helper (flat-bottom) - used for all 3 groups
    macro_rules! dist_grad {
        ($constraints:expr) => {
            for dc in $constraints {
                let dx = c(dc.i, 0) - c(dc.j, 0);
                let dy = c(dc.i, 1) - c(dc.j, 1);
                let dz = c(dc.i, 2) - c(dc.j, 2);
                let d2 = dx*dx + dy*dy + dz*dz;
                let pre;
                if d2 < dc.min_len * dc.min_len {
                    let d = d2.sqrt().max(1e-8);
                    pre = dc.k * (d - dc.min_len) / d;
                } else if d2 > dc.max_len * dc.max_len {
                    let d = d2.sqrt().max(1e-8);
                    pre = dc.k * (d - dc.max_len) / d;
                } else {
                    continue;
                }
                let gi = dc.i * 3;
                let gj = dc.j * 3;
                grad[gi]   += pre * dx;
                grad[gi+1] += pre * dy;
                grad[gi+2] += pre * dz;
                grad[gj]   -= pre * dx;
                grad[gj+1] -= pre * dy;
                grad[gj+2] -= pre * dz;
            }
        };
    }

    // 1-2 distance constraints gradient (contribs[2])
    dist_grad!(&ff.dist_12);

    // Angle constraints gradient (contribs[3])
    for ac in &ff.angle_constraints {
        let r1x = c(ac.i, 0) - c(ac.j, 0);
        let r1y = c(ac.i, 1) - c(ac.j, 1);
        let r1z = c(ac.i, 2) - c(ac.j, 2);
        let r2x = c(ac.k, 0) - c(ac.j, 0);
        let r2y = c(ac.k, 1) - c(ac.j, 1);
        let r2z = c(ac.k, 2) - c(ac.j, 2);
        let l1 = (r1x*r1x + r1y*r1y + r1z*r1z).sqrt();
        let l2 = (r2x*r2x + r2y*r2y + r2z*r2z).sqrt();
        if l1 < 1e-8 || l2 < 1e-8 { continue; }
        let cos_theta = ((r1x*r2x + r1y*r2y + r1z*r2z) / (l1 * l2)).clamp(-1.0, 1.0);
        let theta_deg = cos_theta.acos() * 180.0 / std::f64::consts::PI;
        let angle_term = if theta_deg < ac.min_deg {
            theta_deg - ac.min_deg
        } else if theta_deg > ac.max_deg {
            theta_deg - ac.max_deg
        } else {
            continue;
        };
        let de_dtheta = 2.0 * ac.force_k * angle_term * (180.0 / std::f64::consts::PI);
        let rpx = r2y*r1z - r2z*r1y;
        let rpy = r2z*r1x - r2x*r1z;
        let rpz = r2x*r1y - r2y*r1x;
        let rp_norm = (rpx*rpx + rpy*rpy + rpz*rpz).sqrt();
        if rp_norm < 1e-8 { continue; }
        let prefactor = de_dtheta / rp_norm;
        let l1sq = l1 * l1;
        let l2sq = l2 * l2;
        let cp1x = r1y*rpz - r1z*rpy;
        let cp1y = r1z*rpx - r1x*rpz;
        let cp1z = r1x*rpy - r1y*rpx;
        let dp1x = -cp1x * prefactor / l1sq;
        let dp1y = -cp1y * prefactor / l1sq;
        let dp1z = -cp1z * prefactor / l1sq;
        let cp3x = r2y*rpz - r2z*rpy;
        let cp3y = r2z*rpx - r2x*rpz;
        let cp3z = r2x*rpy - r2y*rpx;
        let dp3x = cp3x * prefactor / l2sq;
        let dp3y = cp3y * prefactor / l2sq;
        let dp3z = cp3z * prefactor / l2sq;
        let gi = ac.i * 3;
        let gj = ac.j * 3;
        let gk = ac.k * 3;
        grad[gi]   += dp1x;
        grad[gi+1] += dp1y;
        grad[gi+2] += dp1z;
        grad[gj]   += -dp1x - dp3x;
        grad[gj+1] += -dp1y - dp3y;
        grad[gj+2] += -dp1z - dp3z;
        grad[gk]   += dp3x;
        grad[gk+1] += dp3y;
        grad[gk+2] += dp3z;
    }

    // 1-3 distance constraints gradient (contribs[4])
    dist_grad!(&ff.dist_13);

    // Long-range distance constraints gradient (contribs[5])
    dist_grad!(&ff.dist_long);

    // UFF-style torsion gradient (usually inactive: torsion_k_omega=0)
    if ff.torsion_k_omega.abs() > 1e-8 && n >= 4 {
        for edge in mol.graph.edge_references() {
            let u = edge.source();
            let v = edge.target();
            let hyb_u = mol.graph[u].hybridization;
            let hyb_v = mol.graph[v].hybridization;
            if hyb_u == crate::graph::Hybridization::SP || hyb_v == crate::graph::Hybridization::SP { continue; }
            let (n_fold, gamma, weight) = crate::forcefield::energy::torsion_params(hyb_u, hyb_v);
            let vv = ff.torsion_k_omega * weight as f64;
            let nf = n_fold as f64;
            let gm = gamma as f64;
            let neighbors_u: Vec<_> = mol.graph.neighbors(u).filter(|&x| x != v).collect();
            let neighbors_v: Vec<_> = mol.graph.neighbors(v).filter(|&x| x != u).collect();
            for &nu in &neighbors_u {
                for &nv in &neighbors_v {
                    let b1 = [c(u.index(),0)-c(nu.index(),0), c(u.index(),1)-c(nu.index(),1), c(u.index(),2)-c(nu.index(),2)];
                    let b2 = [c(v.index(),0)-c(u.index(),0), c(v.index(),1)-c(u.index(),1), c(v.index(),2)-c(u.index(),2)];
                    let b3 = [c(nv.index(),0)-c(v.index(),0), c(nv.index(),1)-c(v.index(),1), c(nv.index(),2)-c(v.index(),2)];
                    let nn1 = [b1[1]*b2[2]-b1[2]*b2[1], b1[2]*b2[0]-b1[0]*b2[2], b1[0]*b2[1]-b1[1]*b2[0]];
                    let nn2 = [b2[1]*b3[2]-b2[2]*b3[1], b2[2]*b3[0]-b2[0]*b3[2], b2[0]*b3[1]-b2[1]*b3[0]];
                    let n1sq = nn1[0]*nn1[0]+nn1[1]*nn1[1]+nn1[2]*nn1[2];
                    let n2sq = nn2[0]*nn2[0]+nn2[1]*nn2[1]+nn2[2]*nn2[2];
                    if n1sq < 1e-10 || n2sq < 1e-10 { continue; }
                    let b2_len = (b2[0]*b2[0]+b2[1]*b2[1]+b2[2]*b2[2]).sqrt();
                    if b2_len < 1e-6 { continue; }
                    let n1l = n1sq.sqrt();
                    let n2l = n2sq.sqrt();
                    let m1 = [
                        (nn1[1]*b2[2]-nn1[2]*b2[1]) / (n1l * b2_len),
                        (nn1[2]*b2[0]-nn1[0]*b2[2]) / (n1l * b2_len),
                        (nn1[0]*b2[1]-nn1[1]*b2[0]) / (n1l * b2_len),
                    ];
                    let x = (nn1[0]*nn2[0]+nn1[1]*nn2[1]+nn1[2]*nn2[2]) / (n1l * n2l);
                    let y = m1[0]*nn2[0]/n2l + m1[1]*nn2[1]/n2l + m1[2]*nn2[2]/n2l;
                    let phi = y.atan2(x);
                    let de_dphi = -vv * nf * (nf * phi - gm).sin();
                    // Blondel & Karplus gradient
                    let gg1 = [b2_len/n1sq * nn1[0], b2_len/n1sq * nn1[1], b2_len/n1sq * nn1[2]];
                    let gg4 = [-b2_len/n2sq * nn2[0], -b2_len/n2sq * nn2[1], -b2_len/n2sq * nn2[2]];
                    let b1_dot_b2 = (b1[0]*b2[0]+b1[1]*b2[1]+b1[2]*b2[2]) / (b2_len*b2_len);
                    let b3_dot_b2 = (b3[0]*b2[0]+b3[1]*b2[1]+b3[2]*b2[2]) / (b2_len*b2_len);
                    let gg2 = [(-b1_dot_b2-1.0)*gg1[0]+b3_dot_b2*gg4[0], (-b1_dot_b2-1.0)*gg1[1]+b3_dot_b2*gg4[1], (-b1_dot_b2-1.0)*gg1[2]+b3_dot_b2*gg4[2]];
                    let gg3 = [(-b3_dot_b2-1.0)*gg4[0]+b1_dot_b2*gg1[0], (-b3_dot_b2-1.0)*gg4[1]+b1_dot_b2*gg1[1], (-b3_dot_b2-1.0)*gg4[2]+b1_dot_b2*gg1[2]];
                    let gi1 = nu.index()*3;
                    let gi2 = u.index()*3;
                    let gi3 = v.index()*3;
                    let gi4 = nv.index()*3;
                    for dd in 0..3 {
                        grad[gi1+dd] += de_dphi * gg1[dd];
                        grad[gi2+dd] += de_dphi * gg2[dd];
                        grad[gi3+dd] += de_dphi * gg3[dd];
                        grad[gi4+dd] += de_dphi * gg4[dd];
                    }
                }
            }
        }
    }

    grad
}

