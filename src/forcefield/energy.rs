use nalgebra::Vector3;
use petgraph::visit::EdgeRef;

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
        // Use a much stronger 1/r^2 based penalty to prevent collapse
        let r2_eff = r2.max(1e-4);
        let val = l2 / r2_eff - 1.0;
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
        }
    }
}

pub fn calculate_total_energy(
    coords: &nalgebra::DMatrix<f32>,
    mol: &crate::graph::Molecule,
    params: &FFParams,
    bounds_matrix: &nalgebra::DMatrix<f32>,
) -> f32 {
    let n = mol.graph.node_count();
    let mut energy = 0.0;

    // 1. Bond Stretch
    for edge in mol.graph.edge_references() {
        let idx1 = edge.source().index();
        let idx2 = edge.target().index();
        let p1 = Vector3::new(coords[(idx1, 0)], coords[(idx1, 1)], coords[(idx1, 2)]);
        let p2 = Vector3::new(coords[(idx2, 0)], coords[(idx2, 1)], coords[(idx2, 2)]);
        let r_eq = crate::distgeom::get_bond_length(mol, edge.source(), edge.target());
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
                let ideal = crate::graph::get_corrected_ideal_angle(mol, ni, n1, n2);
                energy += angle_bend_energy(&p1, &pc, &p2, params.k_theta, ideal);
            }
        }
    }

    // 3. Distance Bounds
    for i in 0..n {
        for j in (i + 1)..n {
            let upper = bounds_matrix[(i, j)];
            let lower = bounds_matrix[(j, i)];
            let p1 = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
            let p2 = Vector3::new(coords[(j, 0)], coords[(j, 1)], coords[(j, 2)]);
            energy += bounds_energy(&p1, &p2, lower, upper, params.k_bounds);
        }
    }

    // 4. Torsional energy
    if n >= 4 {
        for edge in mol.graph.edge_references() {
            let u = edge.source();
            let v = edge.target();
            if mol.graph[u].hybridization == crate::graph::Hybridization::SP
                || mol.graph[v].hybridization == crate::graph::Hybridization::SP
            {
                continue;
            }
            let mut neighbors_u = Vec::new();
            for n in mol.graph.neighbors(u) {
                if n != v {
                    neighbors_u.push(n);
                }
            }
            let mut neighbors_v = Vec::new();
            for n in mol.graph.neighbors(v) {
                if n != u {
                    neighbors_v.push(n);
                }
            }

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
                    energy += torsional_energy(&p1, &p2, &p3, &p4, params.k_omega, 3.0, 0.0);
                }
            }
        }
    }

    energy
}
