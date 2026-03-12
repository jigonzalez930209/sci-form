//! Energy functions for the ETKDG 3D force field (f32 and f64 versions).

use nalgebra::{DMatrix, Vector3};
use petgraph::visit::EdgeRef;
use super::*;

/// Angle computed in degrees. E = k * angleTerm^2.
pub(crate) fn angle_constraint_energy(coords: &DMatrix<f32>, ac: &AngleConstraint) -> f32 {
    let p1 = Vector3::new(coords[(ac.i, 0)], coords[(ac.i, 1)], coords[(ac.i, 2)]);
    let p2 = Vector3::new(coords[(ac.j, 0)], coords[(ac.j, 1)], coords[(ac.j, 2)]);
    let p3 = Vector3::new(coords[(ac.k, 0)], coords[(ac.k, 1)], coords[(ac.k, 2)]);
    let r1 = p1 - p2;
    let r2 = p3 - p2;
    let l1 = r1.norm();
    let l2 = r2.norm();
    if l1 < 1e-8 || l2 < 1e-8 { return 0.0; }
    let cos_theta = (r1.dot(&r2) / (l1 * l2)).clamp(-1.0, 1.0);
    let theta_deg = cos_theta.acos() * 180.0 / std::f32::consts::PI;
    let angle_term = if theta_deg < ac.min_deg as f32 {
        theta_deg - ac.min_deg as f32
    } else if theta_deg > ac.max_deg as f32 {
        theta_deg - ac.max_deg as f32
    } else {
        0.0
    };
    ac.force_k as f32 * angle_term * angle_term
}

/// Flat-bottom angle constraint gradient matching RDKit's AngleConstraintContribs.
/// Uses cross-product formulation: dE/dp_i = dE/dTheta * (r2 x (r1 x r2)) / (|r1 x r2| * |r1|^2)
pub(crate) fn angle_constraint_gradient(coords: &DMatrix<f32>, ac: &AngleConstraint, grad: &mut DMatrix<f32>) {
    let p1 = Vector3::new(coords[(ac.i, 0)], coords[(ac.i, 1)], coords[(ac.i, 2)]);
    let p2 = Vector3::new(coords[(ac.j, 0)], coords[(ac.j, 1)], coords[(ac.j, 2)]);
    let p3 = Vector3::new(coords[(ac.k, 0)], coords[(ac.k, 1)], coords[(ac.k, 2)]);
    let r1 = p1 - p2;
    let r2 = p3 - p2;
    let l1 = r1.norm();
    let l2 = r2.norm();
    if l1 < 1e-8 || l2 < 1e-8 { return; }
    let cos_theta = (r1.dot(&r2) / (l1 * l2)).clamp(-1.0, 1.0);
    let theta_deg = cos_theta.acos() * 180.0 / std::f32::consts::PI;
    let angle_term = if theta_deg < ac.min_deg as f32 {
        theta_deg - ac.min_deg as f32
    } else if theta_deg > ac.max_deg as f32 {
        theta_deg - ac.max_deg as f32
    } else {
        return; // no gradient if angle is within bounds
    };

    // dE/dTheta (in radians) = 2 * k * angleTerm * (180/PI)
    let de_dtheta = 2.0 * ac.force_k as f32 * angle_term * (180.0 / std::f32::consts::PI);

    // Cross product for gradient computation
    let rp = r2.cross(&r1);
    let rp_norm = rp.norm();
    if rp_norm < 1e-8 { return; }
    let prefactor = de_dtheta / rp_norm;

    // Gradient on atom i: dE/dp1 = -(r1 x rp) * prefactor / |r1|^2
    let dedp1 = -(r1.cross(&rp)) * prefactor / (l1 * l1);
    // Gradient on atom k: dE/dp3 =  (r2 x rp) * prefactor / |r2|^2
    let dedp3 = (r2.cross(&rp)) * prefactor / (l2 * l2);
    // Gradient on center atom j: negative sum
    let dedp2 = -dedp1 - dedp3;

    for d in 0..3 {
        grad[(ac.i, d)] += dedp1[d];
        grad[(ac.j, d)] += dedp2[d];
        grad[(ac.k, d)] += dedp3[d];
    }
}

/// Calculate cosY = cos(Wilson angle) matching RDKit's calculateCosY.
/// p1=neighbor1, p2=center, p3=neighbor2, p4=neighbor3
pub(crate) fn calculate_cos_y(p1: &Vector3<f32>, p2: &Vector3<f32>, p3: &Vector3<f32>, p4: &Vector3<f32>) -> f32 {
    let rji = p1 - p2;
    let rjk = p3 - p2;
    let rjl = p4 - p2;
    let l2ji = rji.norm_squared();
    let l2jk = rjk.norm_squared();
    let l2jl = rjl.norm_squared();
    if l2ji < 1e-16 || l2jk < 1e-16 || l2jl < 1e-16 {
        return 0.0;
    }
    let mut n = rji.cross(&rjk);
    n /= l2ji.sqrt() * l2jk.sqrt();
    let l2n = n.norm_squared();
    if l2n < 1e-16 {
        return 0.0;
    }
    n.dot(&rjl) / (l2jl.sqrt() * l2n.sqrt())
}

/// UFF Inversion energy: E = K * (C0 + C1*sinY + C2*cos2W)
pub(crate) fn uff_inversion_energy(coords: &DMatrix<f32>, ic: &UFFInversionContrib) -> f32 {
    let p1 = Vector3::new(coords[(ic.at1, 0)], coords[(ic.at1, 1)], coords[(ic.at1, 2)]);
    let p2 = Vector3::new(coords[(ic.at2, 0)], coords[(ic.at2, 1)], coords[(ic.at2, 2)]);
    let p3 = Vector3::new(coords[(ic.at3, 0)], coords[(ic.at3, 1)], coords[(ic.at3, 2)]);
    let p4 = Vector3::new(coords[(ic.at4, 0)], coords[(ic.at4, 1)], coords[(ic.at4, 2)]);
    let cos_y = calculate_cos_y(&p1, &p2, &p3, &p4);
    let sin_y_sq = (1.0 - cos_y * cos_y).max(0.0);
    let sin_y = sin_y_sq.sqrt();
    // cos(2W) = 2*sin²(Y) - 1  (since W = π/2 - Y, and cos(2W) = 2cos²(W)-1 = 2sin²(Y)-1)
    let cos_2w = 2.0 * sin_y * sin_y - 1.0;
    ic.force_constant as f32 * (ic.c0 as f32 + ic.c1 as f32 * sin_y + ic.c2 as f32 * cos_2w)
}

/// UFF Inversion gradient matching RDKit's InversionContrib::getGrad exactly
#[allow(dead_code)]
pub(crate) fn uff_inversion_gradient(coords: &DMatrix<f32>, ic: &UFFInversionContrib, grad: &mut DMatrix<f32>) {
    let p1 = Vector3::new(coords[(ic.at1, 0)], coords[(ic.at1, 1)], coords[(ic.at1, 2)]);
    let p2 = Vector3::new(coords[(ic.at2, 0)], coords[(ic.at2, 1)], coords[(ic.at2, 2)]);
    let p3 = Vector3::new(coords[(ic.at3, 0)], coords[(ic.at3, 1)], coords[(ic.at3, 2)]);
    let p4 = Vector3::new(coords[(ic.at4, 0)], coords[(ic.at4, 1)], coords[(ic.at4, 2)]);

    let mut rji = p1 - p2;
    let mut rjk = p3 - p2;
    let mut rjl = p4 - p2;
    let dji = rji.norm();
    let djk = rjk.norm();
    let djl = rjl.norm();
    if dji < 1e-8 || djk < 1e-8 || djl < 1e-8 { return; }
    rji /= dji;
    rjk /= djk;
    rjl /= djl;

    let mut n = (-rji).cross(&rjk);
    let n_len = n.norm();
    if n_len < 1e-8 { return; }
    n /= n_len;

    let mut cos_y = n.dot(&rjl);
    cos_y = cos_y.clamp(-1.0, 1.0);
    let sin_y_sq = 1.0 - cos_y * cos_y;
    let sin_y = sin_y_sq.sqrt().max(1e-8);

    let mut cos_theta = rji.dot(&rjk);
    cos_theta = cos_theta.clamp(-1.0, 1.0);
    let sin_theta_sq = 1.0 - cos_theta * cos_theta;
    let sin_theta = sin_theta_sq.sqrt().max(1e-8);

    // dE/dW = -K * (C1*cosY - 4*C2*cosY*sinY)
    let de_dw = -ic.force_constant as f32 * (ic.c1 as f32 * cos_y - 4.0 * ic.c2 as f32 * cos_y * sin_y);

    let t1 = rjl.cross(&rjk);
    let t2 = rji.cross(&rjl);
    let t3 = rjk.cross(&rji);

    let term1 = sin_y * sin_theta;
    let term2 = cos_y / (sin_y * sin_theta_sq);

    // Gradient for atom 1 (neighbor)
    let tg1 = Vector3::new(
        (t1[0] / term1 - (rji[0] - rjk[0] * cos_theta) * term2) / dji,
        (t1[1] / term1 - (rji[1] - rjk[1] * cos_theta) * term2) / dji,
        (t1[2] / term1 - (rji[2] - rjk[2] * cos_theta) * term2) / dji,
    );
    // Gradient for atom 3 (neighbor)
    let tg3 = Vector3::new(
        (t2[0] / term1 - (rjk[0] - rji[0] * cos_theta) * term2) / djk,
        (t2[1] / term1 - (rjk[1] - rji[1] * cos_theta) * term2) / djk,
        (t2[2] / term1 - (rjk[2] - rji[2] * cos_theta) * term2) / djk,
    );
    // Gradient for atom 4 (neighbor)
    let tg4 = Vector3::new(
        (t3[0] / term1 - rjl[0] * cos_y / sin_y) / djl,
        (t3[1] / term1 - rjl[1] * cos_y / sin_y) / djl,
        (t3[2] / term1 - rjl[2] * cos_y / sin_y) / djl,
    );

    for d in 0..3 {
        grad[(ic.at1, d)] += de_dw * tg1[d];
        grad[(ic.at2, d)] += -de_dw * (tg1[d] + tg3[d] + tg4[d]);
        grad[(ic.at3, d)] += de_dw * tg3[d];
        grad[(ic.at4, d)] += de_dw * tg4[d];
    }
}

/// Calculate UFF inversion energy for planarity checking (no torsion/distance terms)
pub fn uff_inversion_energy_only(
    coords: &DMatrix<f32>,
    inversion_contribs: &[UFFInversionContrib],
) -> f32 {
    let mut energy = 0.0f32;
    for ic in inversion_contribs {
        energy += uff_inversion_energy(coords, ic);
    }
    energy
}

/// Compute planarity check energy matching RDKit's construct3DImproperForceField.
/// This creates a force field with ONLY improper (UFF inversion) and SP angle terms,
/// then returns the total energy. RDKit uses oobForceScalingFactor=10.0 for angle
/// constraints in the planarity check FF (vs 1.0 in the main FF).
pub fn planarity_check_energy(coords: &DMatrix<f32>, ff: &Etkdg3DFF) -> f32 {
    let mut e = 0.0f32;
    for ic in &ff.inversion_contribs {
        e += uff_inversion_energy(coords, ic);
    }
    // SP angle constraints with force_k=10.0 (matching RDKit's planarity check FF)
    for ac in &ff.angle_constraints {
        let scaled = AngleConstraint {
            force_k: 10.0,
            ..*ac
        };
        e += angle_constraint_energy(coords, &scaled);
    }
    e
}

/// f64-precision planarity check energy matching RDKit's construct3DImproperForceField.
/// Computes UFF inversion energy + SP angle constraint energy (k=10.0) in f64.
pub fn planarity_check_energy_f64(coords_flat: &[f64], _n: usize, ff: &Etkdg3DFF) -> f64 {
    let c = |atom: usize, d: usize| -> f64 { coords_flat[atom * 3 + d] };
    let mut energy = 0.0f64;

    // UFF Inversion contribs (same as etkdg_3d_energy_f64)
    for ic in &ff.inversion_contribs {
        let rji = [c(ic.at1,0)-c(ic.at2,0), c(ic.at1,1)-c(ic.at2,1), c(ic.at1,2)-c(ic.at2,2)];
        let rjk = [c(ic.at3,0)-c(ic.at2,0), c(ic.at3,1)-c(ic.at2,1), c(ic.at3,2)-c(ic.at2,2)];
        let rjl = [c(ic.at4,0)-c(ic.at2,0), c(ic.at4,1)-c(ic.at2,1), c(ic.at4,2)-c(ic.at2,2)];
        let dji = (rji[0]*rji[0]+rji[1]*rji[1]+rji[2]*rji[2]).sqrt();
        let djk = (rjk[0]*rjk[0]+rjk[1]*rjk[1]+rjk[2]*rjk[2]).sqrt();
        let djl = (rjl[0]*rjl[0]+rjl[1]*rjl[1]+rjl[2]*rjl[2]).sqrt();
        if dji < 1e-8 || djk < 1e-8 || djl < 1e-8 { continue; }
        let nji = [rji[0]/dji, rji[1]/dji, rji[2]/dji];
        let njk = [rjk[0]/djk, rjk[1]/djk, rjk[2]/djk];
        let njl = [rjl[0]/djl, rjl[1]/djl, rjl[2]/djl];
        let nx = -nji[1]*njk[2]+nji[2]*njk[1];
        let ny = -nji[2]*njk[0]+nji[0]*njk[2];
        let nz = -nji[0]*njk[1]+nji[1]*njk[0];
        let nl = (nx*nx+ny*ny+nz*nz).sqrt();
        if nl < 1e-8 { continue; }
        let cos_y = ((nx*njl[0]+ny*njl[1]+nz*njl[2])/nl).clamp(-1.0, 1.0);
        let sin_y_sq = (1.0 - cos_y * cos_y).max(0.0);
        let sin_y = sin_y_sq.sqrt();
        let cos_2w = 2.0 * sin_y * sin_y - 1.0;
        energy += ic.force_constant * (ic.c0 + ic.c1 * sin_y + ic.c2 * cos_2w);
    }

    // SP angle constraints with force_k=10.0 (matching RDKit's planarity check FF)
    for ac in &ff.angle_constraints {
        let r1x = c(ac.i as usize, 0) - c(ac.j as usize, 0);
        let r1y = c(ac.i as usize, 1) - c(ac.j as usize, 1);
        let r1z = c(ac.i as usize, 2) - c(ac.j as usize, 2);
        let r2x = c(ac.k as usize, 0) - c(ac.j as usize, 0);
        let r2y = c(ac.k as usize, 1) - c(ac.j as usize, 1);
        let r2z = c(ac.k as usize, 2) - c(ac.j as usize, 2);
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
            0.0
        };
        energy += 10.0 * angle_term * angle_term;
    }

    energy
}

/// Calculate total energy using the 3D ETKDG force field
pub fn etkdg_3d_energy(
    coords: &DMatrix<f32>,
    mol: &crate::graph::Molecule,
    ff: &Etkdg3DFF,
) -> f32 {
    let n = mol.graph.node_count();
    let mut energy = 0.0f32;

    // Distance constraints (flat-bottom)
    for c in ff.dist_12.iter().chain(ff.dist_13.iter()).chain(ff.dist_long.iter()) {
        let p1 = Vector3::new(coords[(c.i, 0)], coords[(c.i, 1)], coords[(c.i, 2)]);
        let p2 = Vector3::new(coords[(c.j, 0)], coords[(c.j, 1)], coords[(c.j, 2)]);
        energy += crate::forcefield::energy::distance_constraint_energy(&p1, &p2, c.min_len as f32, c.max_len as f32, c.k as f32);
    }

    // Angle constraints (flat-bottom on angle in degrees)
    for ac in &ff.angle_constraints {
        energy += angle_constraint_energy(coords, ac);
    }

    // OOP for SP2 with 3 neighbors (simple vol² formulation)
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
            let v1 = p1 - pc;
            let v2 = p2 - pc;
            let v3 = p3 - pc;
            let vol = v1.dot(&v2.cross(&v3));
            energy += ff.oop_k as f32 * vol * vol;
        }
    }

    // Pre-computed torsion contributions (flat ring + chain preferences)
    for tc in &ff.torsion_contribs {
        let p1 = Vector3::new(coords[(tc.i, 0)], coords[(tc.i, 1)], coords[(tc.i, 2)]);
        let p2 = Vector3::new(coords[(tc.j, 0)], coords[(tc.j, 1)], coords[(tc.j, 2)]);
        let p3 = Vector3::new(coords[(tc.k, 0)], coords[(tc.k, 1)], coords[(tc.k, 2)]);
        let p4 = Vector3::new(coords[(tc.l, 0)], coords[(tc.l, 1)], coords[(tc.l, 2)]);
        let m6 = crate::forcefield::etkdg_lite::M6Params {
            s: tc.signs.map(|x| x as f32),
            v: tc.v.map(|x| x as f32),
        };
        energy += crate::forcefield::etkdg_lite::calc_torsion_energy_m6(&p1, &p2, &p3, &p4, &m6);
    }

    // UFF-style torsions
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
                    energy += crate::forcefield::energy::torsional_energy(&p1, &p2, &p3, &p4, ff.torsion_k_omega as f32 * weight, n_fold, gamma);
                }
            }
        }
    }

    // ETKDG-lite M6 torsion preferences
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
            energy += crate::forcefield::etkdg_lite::calc_torsion_energy_m6(&p1, &p2, &p3, &p4, &m6);
        }
    }

    energy
}


/// f64-precision version of etkdg_3d_energy.
/// All computation done in f64 to match RDKit's double-precision force field.
pub fn etkdg_3d_energy_f64(
    coords: &[f64],
    n: usize,
    mol: &crate::graph::Molecule,
    ff: &Etkdg3DFF,
) -> f64 {
    let c = |atom: usize, d: usize| -> f64 { coords[atom * 3 + d] };
    let mut energy = 0.0f64;

    // === RDKit accumulation order: torsions → inversions → dist_12 → angles → dist_13 → dist_long ===

    // Pre-computed M6 torsion contributions (contribs[0])
    for tc in &ff.torsion_contribs {
        let r1 = [c(tc.i,0)-c(tc.j,0), c(tc.i,1)-c(tc.j,1), c(tc.i,2)-c(tc.j,2)];
        let r2 = [c(tc.k,0)-c(tc.j,0), c(tc.k,1)-c(tc.j,1), c(tc.k,2)-c(tc.j,2)];
        let r3 = [c(tc.j,0)-c(tc.k,0), c(tc.j,1)-c(tc.k,1), c(tc.j,2)-c(tc.k,2)];
        let r4 = [c(tc.l,0)-c(tc.k,0), c(tc.l,1)-c(tc.k,1), c(tc.l,2)-c(tc.k,2)];
        let t1 = [r1[1]*r2[2]-r1[2]*r2[1], r1[2]*r2[0]-r1[0]*r2[2], r1[0]*r2[1]-r1[1]*r2[0]];
        let t2 = [r3[1]*r4[2]-r3[2]*r4[1], r3[2]*r4[0]-r3[0]*r4[2], r3[0]*r4[1]-r3[1]*r4[0]];
        let d1 = (t1[0]*t1[0]+t1[1]*t1[1]+t1[2]*t1[2]).sqrt();
        let d2 = (t2[0]*t2[0]+t2[1]*t2[1]+t2[2]*t2[2]).sqrt();
        if d1 < 1e-10 || d2 < 1e-10 { continue; }
        let n1 = [t1[0]/d1, t1[1]/d1, t1[2]/d1];
        let n2 = [t2[0]/d2, t2[1]/d2, t2[2]/d2];
        let cos_phi = (n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2]).clamp(-1.0, 1.0);
        let cp2 = cos_phi*cos_phi;
        let cp3 = cos_phi*cp2;
        let cp4 = cos_phi*cp3;
        let cp5 = cos_phi*cp4;
        let cp6 = cos_phi*cp5;
        let cos2 = 2.0*cp2 - 1.0;
        let cos3 = 4.0*cp3 - 3.0*cos_phi;
        let cos4 = 8.0*cp4 - 8.0*cp2 + 1.0;
        let cos5 = 16.0*cp5 - 20.0*cp3 + 5.0*cos_phi;
        let cos6 = 32.0*cp6 - 48.0*cp4 + 18.0*cp2 - 1.0;
        let v = &tc.v;
        let s = &tc.signs;
        energy += v[0] * (1.0 + s[0] * cos_phi);
        energy += v[1] * (1.0 + s[1] * cos2);
        energy += v[2] * (1.0 + s[2] * cos3);
        energy += v[3] * (1.0 + s[3] * cos4);
        energy += v[4] * (1.0 + s[4] * cos5);
        energy += v[5] * (1.0 + s[5] * cos6);
    }

    // UFF Inversion contribs (contribs[1])
    for ic in &ff.inversion_contribs {
        let p1 = [c(ic.at1,0), c(ic.at1,1), c(ic.at1,2)];
        let p2 = [c(ic.at2,0), c(ic.at2,1), c(ic.at2,2)];
        let p3 = [c(ic.at3,0), c(ic.at3,1), c(ic.at3,2)];
        let p4 = [c(ic.at4,0), c(ic.at4,1), c(ic.at4,2)];
        let rji = [p1[0]-p2[0], p1[1]-p2[1], p1[2]-p2[2]];
        let rjk = [p3[0]-p2[0], p3[1]-p2[1], p3[2]-p2[2]];
        let rjl = [p4[0]-p2[0], p4[1]-p2[1], p4[2]-p2[2]];
        let dji = (rji[0]*rji[0]+rji[1]*rji[1]+rji[2]*rji[2]).sqrt();
        let djk = (rjk[0]*rjk[0]+rjk[1]*rjk[1]+rjk[2]*rjk[2]).sqrt();
        let djl = (rjl[0]*rjl[0]+rjl[1]*rjl[1]+rjl[2]*rjl[2]).sqrt();
        if dji < 1e-8 || djk < 1e-8 || djl < 1e-8 { continue; }
        let nji = [rji[0]/dji, rji[1]/dji, rji[2]/dji];
        let njk = [rjk[0]/djk, rjk[1]/djk, rjk[2]/djk];
        let njl = [rjl[0]/djl, rjl[1]/djl, rjl[2]/djl];
        let nx = -nji[1]*njk[2]+nji[2]*njk[1];
        let ny = -nji[2]*njk[0]+nji[0]*njk[2];
        let nz = -nji[0]*njk[1]+nji[1]*njk[0];
        let nl = (nx*nx+ny*ny+nz*nz).sqrt();
        if nl < 1e-8 { continue; }
        let cos_y = ((nx*njl[0]+ny*njl[1]+nz*njl[2])/nl).clamp(-1.0, 1.0);
        let sin_y_sq = (1.0 - cos_y * cos_y).max(0.0);
        let sin_y = sin_y_sq.sqrt();
        let cos_2w = 2.0 * sin_y * sin_y - 1.0;
        let k = ic.force_constant;
        let c0 = ic.c0;
        let c1 = ic.c1;
        let c2 = ic.c2;
        energy += k * (c0 + c1 * sin_y + c2 * cos_2w);
    }

    // 1-2 distance constraints (contribs[2])
    for dc in &ff.dist_12 {
        let dx = c(dc.i, 0) - c(dc.j, 0);
        let dy = c(dc.i, 1) - c(dc.j, 1);
        let dz = c(dc.i, 2) - c(dc.j, 2);
        let d2 = dx * dx + dy * dy + dz * dz;
        if d2 < dc.min_len * dc.min_len {
            let d = d2.sqrt();
            let diff = dc.min_len - d;
            energy += 0.5 * dc.k * diff * diff;
        } else if d2 > dc.max_len * dc.max_len {
            let d = d2.sqrt();
            let diff = d - dc.max_len;
            energy += 0.5 * dc.k * diff * diff;
        }
    }

    // Angle constraints (contribs[3])
    for ac in &ff.angle_constraints {
        let r1x = c(ac.i, 0) - c(ac.j, 0);
        let r1y = c(ac.i, 1) - c(ac.j, 1);
        let r1z = c(ac.i, 2) - c(ac.j, 2);
        let r2x = c(ac.k, 0) - c(ac.j, 0);
        let r2y = c(ac.k, 1) - c(ac.j, 1);
        let r2z = c(ac.k, 2) - c(ac.j, 2);
        let l1 = (r1x * r1x + r1y * r1y + r1z * r1z).sqrt();
        let l2 = (r2x * r2x + r2y * r2y + r2z * r2z).sqrt();
        if l1 < 1e-8 || l2 < 1e-8 { continue; }
        let cos_theta = ((r1x * r2x + r1y * r2y + r1z * r2z) / (l1 * l2)).clamp(-1.0, 1.0);
        let theta_deg = cos_theta.acos() * 180.0 / std::f64::consts::PI;
        let angle_term = if theta_deg < ac.min_deg {
            theta_deg - ac.min_deg
        } else if theta_deg > ac.max_deg {
            theta_deg - ac.max_deg
        } else {
            0.0
        };
        energy += ac.force_k * angle_term * angle_term;
    }

    // 1-3 distance constraints (contribs[4])
    for dc in &ff.dist_13 {
        let dx = c(dc.i, 0) - c(dc.j, 0);
        let dy = c(dc.i, 1) - c(dc.j, 1);
        let dz = c(dc.i, 2) - c(dc.j, 2);
        let d2 = dx * dx + dy * dy + dz * dz;
        if d2 < dc.min_len * dc.min_len {
            let d = d2.sqrt();
            let diff = dc.min_len - d;
            energy += 0.5 * dc.k * diff * diff;
        } else if d2 > dc.max_len * dc.max_len {
            let d = d2.sqrt();
            let diff = d - dc.max_len;
            energy += 0.5 * dc.k * diff * diff;
        }
    }

    // Long-range distance constraints (contribs[5])
    for dc in &ff.dist_long {
        let dx = c(dc.i, 0) - c(dc.j, 0);
        let dy = c(dc.i, 1) - c(dc.j, 1);
        let dz = c(dc.i, 2) - c(dc.j, 2);
        let d2 = dx * dx + dy * dy + dz * dz;
        if d2 < dc.min_len * dc.min_len {
            let d = d2.sqrt();
            let diff = dc.min_len - d;
            energy += 0.5 * dc.k * diff * diff;
        } else if d2 > dc.max_len * dc.max_len {
            let d = d2.sqrt();
            let diff = d - dc.max_len;
            energy += 0.5 * dc.k * diff * diff;
        }
    }

    // UFF-style torsions (usually inactive: torsion_k_omega=0)
    if ff.torsion_k_omega.abs() > 1e-8 && n >= 4 {
        let tk = ff.torsion_k_omega;
        for edge in mol.graph.edge_references() {
            let u = edge.source();
            let v = edge.target();
            let hyb_u = mol.graph[u].hybridization;
            let hyb_v = mol.graph[v].hybridization;
            if hyb_u == crate::graph::Hybridization::SP || hyb_v == crate::graph::Hybridization::SP { continue; }
            let (n_fold, gamma, weight) = crate::forcefield::energy::torsion_params(hyb_u, hyb_v);
            let nf = n_fold as f64;
            let gm = gamma as f64;
            let wt = weight as f64;
            let neighbors_u: Vec<_> = mol.graph.neighbors(u).filter(|&x| x != v).collect();
            let neighbors_v: Vec<_> = mol.graph.neighbors(v).filter(|&x| x != u).collect();
            for &nu in &neighbors_u {
                for &nv in &neighbors_v {
                    let b1 = [c(u.index(),0)-c(nu.index(),0), c(u.index(),1)-c(nu.index(),1), c(u.index(),2)-c(nu.index(),2)];
                    let b2 = [c(v.index(),0)-c(u.index(),0), c(v.index(),1)-c(u.index(),1), c(v.index(),2)-c(u.index(),2)];
                    let b3 = [c(nv.index(),0)-c(v.index(),0), c(nv.index(),1)-c(v.index(),1), c(nv.index(),2)-c(v.index(),2)];
                    let nn1 = [b1[1]*b2[2]-b1[2]*b2[1], b1[2]*b2[0]-b1[0]*b2[2], b1[0]*b2[1]-b1[1]*b2[0]];
                    let nn2 = [b2[1]*b3[2]-b2[2]*b3[1], b2[2]*b3[0]-b2[0]*b3[2], b2[0]*b3[1]-b2[1]*b3[0]];
                    let nn1_l = (nn1[0]*nn1[0]+nn1[1]*nn1[1]+nn1[2]*nn1[2]).sqrt();
                    let nn2_l = (nn2[0]*nn2[0]+nn2[1]*nn2[1]+nn2[2]*nn2[2]).sqrt();
                    if nn1_l < 1e-8 || nn2_l < 1e-8 { continue; }
                    let nn1n = [nn1[0]/nn1_l, nn1[1]/nn1_l, nn1[2]/nn1_l];
                    let nn2n = [nn2[0]/nn2_l, nn2[1]/nn2_l, nn2[2]/nn2_l];
                    let b2_l = (b2[0]*b2[0]+b2[1]*b2[1]+b2[2]*b2[2]).sqrt();
                    if b2_l < 1e-8 { continue; }
                    let b2n = [b2[0]/b2_l, b2[1]/b2_l, b2[2]/b2_l];
                    let m1 = [nn1n[1]*b2n[2]-nn1n[2]*b2n[1], nn1n[2]*b2n[0]-nn1n[0]*b2n[2], nn1n[0]*b2n[1]-nn1n[1]*b2n[0]];
                    let x = nn1n[0]*nn2n[0]+nn1n[1]*nn2n[1]+nn1n[2]*nn2n[2];
                    let y = m1[0]*nn2n[0]+m1[1]*nn2n[1]+m1[2]*nn2n[2];
                    let phi = y.atan2(x);
                    energy += tk * wt * (1.0 + (nf * phi - gm).cos());
                }
            }
        }
    }

    energy
}
