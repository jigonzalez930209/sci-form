use nalgebra::Vector3;
use petgraph::visit::EdgeRef;

/// UFF VDW parameters by element: (x1 = VDW distance in Å, d1 = well depth in kcal/mol)
pub fn uff_vdw_params(element: u8) -> (f32, f32) {
    match element {
        1  => (2.886, 0.044),   // H
        5  => (3.637, 0.180),   // B
        6  => (3.851, 0.105),   // C
        7  => (3.660, 0.069),   // N
        8  => (3.500, 0.060),   // O
        9  => (3.364, 0.050),   // F
        14 => (4.295, 0.402),   // Si
        15 => (4.147, 0.305),   // P
        16 => (4.035, 0.274),   // S
        17 => (3.947, 0.227),   // Cl
        35 => (4.189, 0.251),   // Br
        53 => (4.500, 0.339),   // I
        _  => (3.851, 0.105),   // default to C
    }
}

/// LJ 12-6 VDW energy between two atoms
pub fn vdw_energy(p1: &Vector3<f32>, p2: &Vector3<f32>, r_star: f32, epsilon: f32) -> f32 {
    let r = (p1 - p2).norm();
    if r < 0.5 || r > 8.0 { return 0.0; }
    let u = r_star / r;
    let u6 = u * u * u * u * u * u;
    let u12 = u6 * u6;
    epsilon * (u12 - 2.0 * u6)
}

/// Flat-bottom distance constraint energy (matching RDKit's DistanceConstraintContribs)
/// Zero energy when minLen <= d <= maxLen, harmonic penalty outside
pub fn distance_constraint_energy(
    p1: &Vector3<f32>,
    p2: &Vector3<f32>,
    min_len: f32,
    max_len: f32,
    k: f32,
) -> f32 {
    let d2 = (p1 - p2).norm_squared();
    if d2 < min_len * min_len {
        let d = d2.sqrt();
        let diff = min_len - d;
        0.5 * k * diff * diff
    } else if d2 > max_len * max_len {
        let d = d2.sqrt();
        let diff = d - max_len;
        0.5 * k * diff * diff
    } else {
        0.0
    }
}

/// Harmonic bond stretching penalty (Hooke's Law)
pub fn bond_stretch_energy(p1: &Vector3<f32>, p2: &Vector3<f32>, k_b: f32, r_eq: f32) -> f32 {
    let r = (p1 - p2).norm();
    0.5 * k_b * (r - r_eq).powi(2)
}

/// Harmonic angle bending penalty
pub fn angle_bend_energy(
    p1: &Vector3<f32>,
    p2: &Vector3<f32>, // central atom
    p3: &Vector3<f32>,
    k_theta: f32,
    theta_eq: f32,
) -> f32 {
    let v1 = p1 - p2;
    let v2 = p3 - p2;
    let r1 = v1.norm();
    let r2 = v2.norm();
    if r1 < 1e-4 || r2 < 1e-4 {
        return 0.0;
    }
    let cos_th = (v1.dot(&v2) / (r1 * r2)).clamp(-0.999999, 0.999999);

    // Linear angle special case: use cosine potential for stability
    if (theta_eq - std::f32::consts::PI).abs() < 1e-4 {
        // Special linear potential: E = k * (1 + cos(theta))
        return k_theta * (1.0 + cos_th);
    }

    let theta = cos_th.acos();
    0.5 * k_theta * (theta - theta_eq).powi(2)
}

/// Harmonic torsional potential (Dihedral)
pub fn torsional_energy(
    p1: &Vector3<f32>,
    p2: &Vector3<f32>,
    p3: &Vector3<f32>,
    p4: &Vector3<f32>,
    v: f32,
    n: f32,
    gamma: f32,
) -> f32 {
    let b1 = p2 - p1;
    let b2 = p3 - p2;
    let b3 = p4 - p3;

    let n1 = b1.cross(&b2).normalize();
    let n2 = b2.cross(&b3).normalize();
    let m1 = n1.cross(&b2.normalize());

    let x = n1.dot(&n2);
    let y = m1.dot(&n2);
    let phi = y.atan2(x);

    v * (1.0 + (n * phi - gamma).cos())
}

/// Distance-based penalty (used in embedding refinement)
/// Matches RDKit's DistViolationContribs formulation
pub fn bounds_energy(
    p1: &Vector3<f32>,
    p2: &Vector3<f32>,
    lower: f32,
    upper: f32,
    k_bounds: f32,
) -> f32 {
    let r2 = (p1 - p2).norm_squared();
    let u2 = upper * upper;
    let l2 = lower * lower;
    if r2 > u2 && u2 > 1e-6 {
        let val = r2 / u2 - 1.0;
        k_bounds * val * val
    } else if r2 < l2 && l2 > 1e-6 {
        // RDKit formula: val = 2L²/(L²+d²) - 1
        let val = 2.0 * l2 / (l2 + r2.max(1e-6)) - 1.0;
        k_bounds * val * val
    } else {
        0.0
    }
}

/// Harmonic Out-of-Plane bending penalty
pub fn oop_energy(
    p_center: &Vector3<f32>,
    p1: &Vector3<f32>,
    p2: &Vector3<f32>,
    p3: &Vector3<f32>,
    k_oop: f32,
    phi_eq: f32,
) -> f32 {
    let v1 = p1 - p_center;
    let v2 = p2 - p_center;
    let v3 = p3 - p_center;

    let normal = v2.cross(&v3).normalize();
    let dist = v1.dot(&normal);
    let sin_phi = dist / v1.norm().max(1e-4);
    let phi = sin_phi.asin();

    0.5 * k_oop * (phi - phi_eq).powi(2)
}

/// Chiral volume penalty
pub fn chirality_energy(
    p_center: &Vector3<f32>,
    p1: &Vector3<f32>,
    p2: &Vector3<f32>,
    p3: &Vector3<f32>,
    target_vol: f32,
    k_chiral: f32,
) -> f32 {
    let v1 = p1 - p_center;
    let v2 = p2 - p_center;
    let v3 = p3 - p_center;
    let vol = v1.dot(&v2.cross(&v3));
    0.5 * k_chiral * (vol - target_vol).powi(2)
}

#[derive(Clone, Debug)]
pub struct FFParams {
    pub kb: f32,
    pub k_theta: f32,
    pub k_omega: f32,
    pub k_oop: f32,
    pub k_bounds: f32,
    pub k_chiral: f32,
    pub k_vdw: f32,
}

impl Default for FFParams {
    fn default() -> Self {
        Self {
            kb: 500.0,
            k_theta: 300.0,
            k_omega: 20.0,
            k_oop: 40.0,
            k_bounds: 200.0,
            k_chiral: 100.0,
            k_vdw: 0.0,
        }
    }
}

pub fn calculate_total_energy(
    coords: &nalgebra::DMatrix<f32>,
    mol: &crate::graph::Molecule,
    params: &FFParams,
    bounds_matrix: &nalgebra::DMatrix<f64>,
) -> f32 {
    let n = mol.graph.node_count();
    let mut energy = 0.0;

    // 1. Bond Stretch
    for edge in mol.graph.edge_references() {
        let idx1 = edge.source().index();
        let idx2 = edge.target().index();
        let p1 = Vector3::new(coords[(idx1, 0)], coords[(idx1, 1)], coords[(idx1, 2)]);
        let p2 = Vector3::new(coords[(idx2, 0)], coords[(idx2, 1)], coords[(idx2, 2)]);
        let r_eq = crate::distgeom::get_bond_length(mol, edge.source(), edge.target()) as f32;
        energy += bond_stretch_energy(&p1, &p2, params.kb, r_eq);
    }

    // 2. Angles
    for i in 0..n {
        let ni = petgraph::graph::NodeIndex::new(i);
        let nbs: Vec<_> = mol.graph.neighbors(ni).collect();
        for j in 0..nbs.len() {
            for k in (j + 1)..nbs.len() {
                let n1 = nbs[j];
                let n2 = nbs[k];
                let p1 = Vector3::new(
                    coords[(n1.index(), 0)],
                    coords[(n1.index(), 1)],
                    coords[(n1.index(), 2)],
                );
                let pc = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
                let p2 = Vector3::new(
                    coords[(n2.index(), 0)],
                    coords[(n2.index(), 1)],
                    coords[(n2.index(), 2)],
                );
                let ideal = crate::graph::get_corrected_ideal_angle(mol, ni, n1, n2) as f32;
                energy += angle_bend_energy(&p1, &pc, &p2, params.k_theta, ideal);
            }
        }
    }

    // 3. Distance Bounds
    for i in 0..n {
        for j in (i + 1)..n {
            let upper = bounds_matrix[(i, j)] as f32;
            let lower = bounds_matrix[(j, i)] as f32;
            let p1 = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
            let p2 = Vector3::new(coords[(j, 0)], coords[(j, 1)], coords[(j, 2)]);
            energy += bounds_energy(&p1, &p2, lower, upper, params.k_bounds);
        }
    }

    // 4. Out-of-Plane bending for SP2 atoms with 3 neighbors (volume-based)
    if params.k_oop.abs() > 1e-8 {
        for i in 0..n {
            let ni = petgraph::graph::NodeIndex::new(i);
            if mol.graph[ni].hybridization != crate::graph::Hybridization::SP2 {
                continue;
            }
            let nbs: Vec<_> = mol.graph.neighbors(ni).collect();
            if nbs.len() != 3 {
                continue;
            }
            let pc = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
            let p1 = Vector3::new(
                coords[(nbs[0].index(), 0)],
                coords[(nbs[0].index(), 1)],
                coords[(nbs[0].index(), 2)],
            );
            let p2 = Vector3::new(
                coords[(nbs[1].index(), 0)],
                coords[(nbs[1].index(), 1)],
                coords[(nbs[1].index(), 2)],
            );
            let p3 = Vector3::new(
                coords[(nbs[2].index(), 0)],
                coords[(nbs[2].index(), 1)],
                coords[(nbs[2].index(), 2)],
            );
            // Triple product: V = (p1-pc)·((p2-pc)×(p3-pc))
            let v1 = p1 - pc;
            let v2 = p2 - pc;
            let v3 = p3 - pc;
            let vol = v1.dot(&v2.cross(&v3));
            energy += params.k_oop * vol * vol;
        }
    }

    // 5. Torsional energy (UFF-style hybridization-dependent)
    if n >= 4 && params.k_omega.abs() > 1e-8 {
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

            // Determine torsion params based on central bond hybridization
            let (n_fold, gamma, weight) = torsion_params(hyb_u, hyb_v);

            let neighbors_u: Vec<_> = mol.graph.neighbors(u).filter(|&x| x != v).collect();
            let neighbors_v: Vec<_> = mol.graph.neighbors(v).filter(|&x| x != u).collect();

            for &nu in &neighbors_u {
                for &nv in &neighbors_v {
                    let (p1, p2, p3, p4) = (
                        Vector3::new(
                            coords[(nu.index(), 0)],
                            coords[(nu.index(), 1)],
                            coords[(nu.index(), 2)],
                        ),
                        Vector3::new(
                            coords[(u.index(), 0)],
                            coords[(u.index(), 1)],
                            coords[(u.index(), 2)],
                        ),
                        Vector3::new(
                            coords[(v.index(), 0)],
                            coords[(v.index(), 1)],
                            coords[(v.index(), 2)],
                        ),
                        Vector3::new(
                            coords[(nv.index(), 0)],
                            coords[(nv.index(), 1)],
                            coords[(nv.index(), 2)],
                        ),
                    );
                    energy += torsional_energy(&p1, &p2, &p3, &p4, params.k_omega * weight, n_fold, gamma);
                }
            }
        }
    }

    // 6. ETKDG-lite M6 torsion preferences (applied to non-ring rotatable bonds)
    if n >= 4 {
        for edge in mol.graph.edge_references() {
            let u = edge.source();
            let v = edge.target();
            // Skip ring bonds — ETKDG preferences are for rotatable bonds
            if crate::graph::min_path_excluding2(mol, u, v, u, v, 7).is_some() {
                continue;
            }
            let m6 = crate::forcefield::etkdg_lite::infer_etkdg_parameters(mol, u.index(), v.index());
            // Skip if all coefficients are zero
            if m6.v.iter().all(|&x| x.abs() < 1e-6) { continue; }

            let neighbors_u: Vec<_> = mol.graph.neighbors(u).filter(|&x| x != v).collect();
            let neighbors_v: Vec<_> = mol.graph.neighbors(v).filter(|&x| x != u).collect();
            if neighbors_u.is_empty() || neighbors_v.is_empty() { continue; }
            // Use first neighbor pair only (matching RDKit's approach for ETKDG)
            let nu = neighbors_u[0];
            let nv = neighbors_v[0];
            let (p1, p2, p3, p4) = (
                Vector3::new(coords[(nu.index(), 0)], coords[(nu.index(), 1)], coords[(nu.index(), 2)]),
                Vector3::new(coords[(u.index(), 0)], coords[(u.index(), 1)], coords[(u.index(), 2)]),
                Vector3::new(coords[(v.index(), 0)], coords[(v.index(), 1)], coords[(v.index(), 2)]),
                Vector3::new(coords[(nv.index(), 0)], coords[(nv.index(), 1)], coords[(nv.index(), 2)]),
            );
            energy += crate::forcefield::etkdg_lite::calc_torsion_energy_m6(&p1, &p2, &p3, &p4, &m6);
        }
    }

    // 7. VDW non-bonded interactions (1-4+ pairs)
    if params.k_vdw.abs() > 1e-8 {
        // Build exclusion set: 1-2 and 1-3 pairs
        let mut excluded = std::collections::HashSet::new();
        for edge in mol.graph.edge_references() {
            let a = edge.source().index();
            let b = edge.target().index();
            let (lo, hi) = if a < b { (a, b) } else { (b, a) };
            excluded.insert((lo, hi));
        }
        // 1-3 pairs: i-center-j for each angle
        for center in 0..n {
            let nc = petgraph::graph::NodeIndex::new(center);
            let nbs: Vec<_> = mol.graph.neighbors(nc).collect();
            for j in 0..nbs.len() {
                for k in (j + 1)..nbs.len() {
                    let a = nbs[j].index();
                    let b = nbs[k].index();
                    let (lo, hi) = if a < b { (a, b) } else { (b, a) };
                    excluded.insert((lo, hi));
                }
            }
        }
        // 1-4 pairs: need to identify for 0.5 scaling
        let mut is_14 = std::collections::HashSet::new();
        for edge in mol.graph.edge_references() {
            let u = edge.source();
            let v = edge.target();
            let neighbors_u: Vec<_> = mol.graph.neighbors(u).filter(|&x| x != v).collect();
            let neighbors_v: Vec<_> = mol.graph.neighbors(v).filter(|&x| x != u).collect();
            for &nu in &neighbors_u {
                for &nv in &neighbors_v {
                    let a = nu.index();
                    let b = nv.index();
                    if a == b { continue; }
                    let (lo, hi) = if a < b { (a, b) } else { (b, a) };
                    if !excluded.contains(&(lo, hi)) {
                        is_14.insert((lo, hi));
                    }
                }
            }
        }

        for i in 0..n {
            let ei = mol.graph[petgraph::graph::NodeIndex::new(i)].element;
            let (xi, di) = uff_vdw_params(ei);
            let pi = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
            for j in (i + 1)..n {
                if excluded.contains(&(i, j)) { continue; }
                let ej = mol.graph[petgraph::graph::NodeIndex::new(j)].element;
                let (xj, dj) = uff_vdw_params(ej);
                let r_star = (xi + xj) * 0.5;
                let eps_full = (di * dj).sqrt();
                let scale = if is_14.contains(&(i, j)) { 0.5 } else { 1.0 };
                let pj = Vector3::new(coords[(j, 0)], coords[(j, 1)], coords[(j, 2)]);
                energy += params.k_vdw * scale * vdw_energy(&pi, &pj, r_star, eps_full);
            }
        }
    }

    energy
}

/// Determine torsion periodicity, phase, and relative weight based on UFF rules
pub fn torsion_params(
    hyb_u: crate::graph::Hybridization,
    hyb_v: crate::graph::Hybridization,
) -> (f32, f32, f32) {
    use crate::graph::Hybridization::*;
    let pi = std::f32::consts::PI;
    match (hyb_u, hyb_v) {
        (SP3, SP3) => (3.0, 0.0, 1.0),              // staggered, normal weight
        (SP2, SP2) => (2.0, pi, 5.0),                // planar, strong weight
        (SP2, SP3) | (SP3, SP2) => (6.0, pi, 0.5),   // 6-fold weak barrier
        _ => (3.0, 0.0, 1.0),
    }
}
