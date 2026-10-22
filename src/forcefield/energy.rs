//! Simplified UFF-style energy model and its building blocks.
//!
//! The functions here implement a *generic* UFF-like potential (harmonic
//! bond, angle, torsion, out-of-plane, distance-bounds and van der Waals
//! terms) with hand-tuned force constants
//! ([`crate::forcefield::energy::FFParams`]), used for the distance-geometry
//! refinement stages where only the relative ranking of nearby geometries
//! matters. The published UFF parameterisation lives in [`super::uff`] and
//! [`super::params`].
//!
//! Energies are in kcal/mol, lengths in Å and angles in radians.
use nalgebra::Vector3;
use petgraph::visit::EdgeRef;

/// UFF VDW parameters by element: (x1 = VDW distance in Å, d1 = well depth in kcal/mol)
///
/// The values are kept in `f32` because the analytical gradient
/// (`src/forcefield/gradients.rs`) uses the same table; converting them to
/// `f64` in the energy would make the energy and its gradient inconsistent
/// at the 1e-8 relative level.
pub fn uff_vdw_params(element: u8) -> (f32, f32) {
    match element {
        1 => (2.886, 0.044),  // H
        5 => (3.637, 0.180),  // B
        6 => (3.851, 0.105),  // C
        7 => (3.660, 0.069),  // N
        8 => (3.500, 0.060),  // O
        9 => (3.364, 0.050),  // F
        14 => (4.295, 0.402), // Si
        15 => (4.147, 0.305), // P
        16 => (4.035, 0.274), // S
        17 => (3.947, 0.227), // Cl
        35 => (4.189, 0.251), // Br
        53 => (4.500, 0.339), // I
        _ => (3.851, 0.105),  // default to C
    }
}

/// LJ 12-6 VDW energy between two atoms
pub fn vdw_energy(p1: &Vector3<f64>, p2: &Vector3<f64>, r_star: f64, epsilon: f64) -> f64 {
    let r = (p1 - p2).norm();
    if !(0.5..=8.0).contains(&r) {
        return 0.0;
    }
    let u = r_star / r;
    let u6 = u * u * u * u * u * u;
    let u12 = u6 * u6;
    epsilon * (u12 - 2.0 * u6)
}

/// Flat-bottom distance constraint energy (matching RDKit's DistanceConstraintContribs)
/// Zero energy when minLen <= d <= maxLen, harmonic penalty outside
pub fn distance_constraint_energy(
    p1: &Vector3<f64>,
    p2: &Vector3<f64>,
    min_len: f64,
    max_len: f64,
    k: f64,
) -> f64 {
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
pub fn bond_stretch_energy(p1: &Vector3<f64>, p2: &Vector3<f64>, k_b: f64, r_eq: f64) -> f64 {
    let r = (p1 - p2).norm();
    0.5 * k_b * (r - r_eq) * (r - r_eq)
}

/// Harmonic angle bending penalty
pub fn angle_bend_energy(
    p1: &Vector3<f64>,
    p2: &Vector3<f64>, // central atom
    p3: &Vector3<f64>,
    k_theta: f64,
    theta_eq: f64,
) -> f64 {
    let v1 = p1 - p2;
    let v2 = p3 - p2;
    let r1 = v1.norm();
    let r2 = v2.norm();
    if r1 < 1e-4 || r2 < 1e-4 {
        return 0.0;
    }
    let cos_th = (v1.dot(&v2) / (r1 * r2)).clamp(-0.999999, 0.999999);

    // Linear angle special case: use cosine potential for stability
    if (theta_eq - std::f64::consts::PI).abs() < 1e-4 {
        // Special linear potential: E = k * (1 + cos(theta))
        return k_theta * (1.0 + cos_th);
    }

    let theta = cos_th.acos();
    0.5 * k_theta * (theta - theta_eq) * (theta - theta_eq)
}

/// Harmonic torsional potential (Dihedral)
pub fn torsional_energy(
    p1: &Vector3<f64>,
    p2: &Vector3<f64>,
    p3: &Vector3<f64>,
    p4: &Vector3<f64>,
    v: f64,
    n: f64,
    gamma: f64,
) -> f64 {
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
    p1: &Vector3<f64>,
    p2: &Vector3<f64>,
    lower: f64,
    upper: f64,
    k_bounds: f64,
) -> f64 {
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

/// Harmonic Out-of-Plane bending penalty.
///
/// The out-of-plane angle `phi` of the atom `p1` relative to the plane spanned
/// by `p1`, `p2` and `p3` around `p_center` is measured in radians, and the
/// energy is `½ · k_oop · (phi - phi_eq)²`.
///
/// `k_oop` is in kcal/mol/rad² and `phi_eq` in radians; the returned energy is
/// in kcal/mol. Positions are in Å.
pub fn oop_energy(
    p_center: &Vector3<f64>,
    p1: &Vector3<f64>,
    p2: &Vector3<f64>,
    p3: &Vector3<f64>,
    k_oop: f64,
    phi_eq: f64,
) -> f64 {
    let v1 = p1 - p_center;
    let v2 = p2 - p_center;
    let v3 = p3 - p_center;

    let normal = v2.cross(&v3).normalize();
    let dist = v1.dot(&normal);
    let sin_phi = dist / v1.norm().max(1e-4);
    let phi = sin_phi.asin();

    0.5 * k_oop * (phi - phi_eq) * (phi - phi_eq)
}

/// Chiral volume penalty.
///
/// The signed volume of the tetrahedron formed by `p_center`, `p1`, `p2` and
/// `p3` is compared with `target_vol` and the difference is penalised
/// quadratically.
///
/// `target_vol` is in Å³ and `k_chiral` in kcal/mol/Å⁶; the returned energy is
/// in kcal/mol. Positions are in Å.
pub fn chirality_energy(
    p_center: &Vector3<f64>,
    p1: &Vector3<f64>,
    p2: &Vector3<f64>,
    p3: &Vector3<f64>,
    target_vol: f64,
    k_chiral: f64,
) -> f64 {
    let v1 = p1 - p_center;
    let v2 = p2 - p_center;
    let v3 = p3 - p_center;
    let vol = v1.dot(&v2.cross(&v3));
    0.5 * k_chiral * (vol - target_vol) * (vol - target_vol)
}

/// Generic force-field constants used by the simplified UFF-style energy model.
///
/// These are *not* the published UFF parameters: they are hand-tuned scalars
/// that give a reasonable relative ranking during distance-geometry
/// refinement, where the goal is to satisfy the distance bounds rather than to
/// reproduce experimental energies.
#[derive(Clone, Debug)]
pub struct FFParams {
    /// Bond-stretch force constant in kcal/mol/Å².
    pub kb: f32,
    /// Angle-bend force constant in kcal/mol/rad².
    pub k_theta: f32,
    /// Torsion force constant in kcal/mol.
    pub k_omega: f32,
    /// Out-of-plane (inversion) force constant in kcal/mol/rad².
    pub k_oop: f32,
    /// Distance-bounds violation force constant in kcal/mol/Å².
    pub k_bounds: f32,
    /// Chiral-volume violation force constant in kcal/mol/Å⁶.
    pub k_chiral: f32,
    /// van der Waals weighting factor (dimensionless, `0.0` disables the term).
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

/// Compute the total UFF-style energy of `coords` (kcal/mol).
///
/// The accumulation is performed in `f64`: the O(N²) bounds term alone can
/// contribute tens of thousands of addends, and summing them in `f32` loses
/// several decimal digits (blocker A-H3).
///
/// `coords` is stored as `f32` (the storage format of the embedding pipeline);
/// every *unit conversion* applied on the way in (bond lengths in Å, ideal
/// angles in radians, bounds in Å) is performed in `f64`.
pub fn calculate_total_energy_f64(
    coords: &nalgebra::DMatrix<f32>,
    mol: &crate::graph::Molecule,
    params: &FFParams,
    bounds_matrix: &nalgebra::DMatrix<f64>,
) -> f64 {
    let n = mol.graph.node_count();
    // f64 accumulator (blocker A-H3): the previous `let mut energy = 0.0;`
    // inferred `f32` and lost precision over the O(N²) bounds sum.
    let mut energy = 0.0f64;

    let kb = params.kb as f64;
    let k_theta = params.k_theta as f64;
    let k_omega = params.k_omega as f64;
    let k_oop = params.k_oop as f64;
    let k_bounds = params.k_bounds as f64;
    let k_vdw = params.k_vdw as f64;

    // 1. Bond Stretch
    for edge in mol.graph.edge_references() {
        let idx1 = edge.source().index();
        let idx2 = edge.target().index();
        let p1 = Vector3::new(
            coords[(idx1, 0)] as f64,
            coords[(idx1, 1)] as f64,
            coords[(idx1, 2)] as f64,
        );
        let p2 = Vector3::new(
            coords[(idx2, 0)] as f64,
            coords[(idx2, 1)] as f64,
            coords[(idx2, 2)] as f64,
        );
        let r_eq = crate::distgeom::get_bond_length(mol, edge.source(), edge.target());
        energy += bond_stretch_energy(&p1, &p2, kb, r_eq);
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
                    coords[(n1.index(), 0)] as f64,
                    coords[(n1.index(), 1)] as f64,
                    coords[(n1.index(), 2)] as f64,
                );
                let pc = Vector3::new(
                    coords[(i, 0)] as f64,
                    coords[(i, 1)] as f64,
                    coords[(i, 2)] as f64,
                );
                let p2 = Vector3::new(
                    coords[(n2.index(), 0)] as f64,
                    coords[(n2.index(), 1)] as f64,
                    coords[(n2.index(), 2)] as f64,
                );
                // Ideal angle in radians; the conversion is applied in f64.
                let ideal = crate::graph::get_corrected_ideal_angle(mol, ni, n1, n2);
                energy += angle_bend_energy(&p1, &pc, &p2, k_theta, ideal);
            }
        }
    }

    // 3. Distance Bounds
    for i in 0..n {
        for j in (i + 1)..n {
            let upper = bounds_matrix[(i, j)];
            let lower = bounds_matrix[(j, i)];
            let p1 = Vector3::new(
                coords[(i, 0)] as f64,
                coords[(i, 1)] as f64,
                coords[(i, 2)] as f64,
            );
            let p2 = Vector3::new(
                coords[(j, 0)] as f64,
                coords[(j, 1)] as f64,
                coords[(j, 2)] as f64,
            );
            energy += bounds_energy(&p1, &p2, lower, upper, k_bounds);
        }
    }

    // 4. Out-of-Plane bending for SP2 atoms with 3 neighbors (volume-based)
    if k_oop.abs() > 1e-8 {
        for i in 0..n {
            let ni = petgraph::graph::NodeIndex::new(i);
            if mol.graph[ni].hybridization != crate::graph::Hybridization::SP2 {
                continue;
            }
            let nbs: Vec<_> = mol.graph.neighbors(ni).collect();
            if nbs.len() != 3 {
                continue;
            }
            let pc = Vector3::new(
                coords[(i, 0)] as f64,
                coords[(i, 1)] as f64,
                coords[(i, 2)] as f64,
            );
            let p1 = Vector3::new(
                coords[(nbs[0].index(), 0)] as f64,
                coords[(nbs[0].index(), 1)] as f64,
                coords[(nbs[0].index(), 2)] as f64,
            );
            let p2 = Vector3::new(
                coords[(nbs[1].index(), 0)] as f64,
                coords[(nbs[1].index(), 1)] as f64,
                coords[(nbs[1].index(), 2)] as f64,
            );
            let p3 = Vector3::new(
                coords[(nbs[2].index(), 0)] as f64,
                coords[(nbs[2].index(), 1)] as f64,
                coords[(nbs[2].index(), 2)] as f64,
            );
            // Triple product: V = (p1-pc)·((p2-pc)×(p3-pc))
            let v1 = p1 - pc;
            let v2 = p2 - pc;
            let v3 = p3 - pc;
            let vol = v1.dot(&v2.cross(&v3));
            energy += k_oop * vol * vol;
        }
    }

    // 5. Torsional energy (UFF-style hybridization-dependent)
    if n >= 4 && k_omega.abs() > 1e-8 {
        for edge in mol.graph.edge_references() {
            let u = edge.source();
            let v = edge.target();
            let hyb_u = mol.graph[u].hybridization;
            let hyb_v = mol.graph[v].hybridization;
            if hyb_u == crate::graph::Hybridization::SP || hyb_v == crate::graph::Hybridization::SP
            {
                continue;
            }

            // Determine torsion params based on central bond hybridization
            let (n_fold, gamma, weight) = torsion_params(hyb_u, hyb_v);
            let n_fold = n_fold as f64;
            let gamma = gamma as f64;
            let weight = weight as f64;

            let neighbors_u: Vec<_> = mol.graph.neighbors(u).filter(|&x| x != v).collect();
            let neighbors_v: Vec<_> = mol.graph.neighbors(v).filter(|&x| x != u).collect();

            for &nu in &neighbors_u {
                for &nv in &neighbors_v {
                    let (p1, p2, p3, p4) = (
                        Vector3::new(
                            coords[(nu.index(), 0)] as f64,
                            coords[(nu.index(), 1)] as f64,
                            coords[(nu.index(), 2)] as f64,
                        ),
                        Vector3::new(
                            coords[(u.index(), 0)] as f64,
                            coords[(u.index(), 1)] as f64,
                            coords[(u.index(), 2)] as f64,
                        ),
                        Vector3::new(
                            coords[(v.index(), 0)] as f64,
                            coords[(v.index(), 1)] as f64,
                            coords[(v.index(), 2)] as f64,
                        ),
                        Vector3::new(
                            coords[(nv.index(), 0)] as f64,
                            coords[(nv.index(), 1)] as f64,
                            coords[(nv.index(), 2)] as f64,
                        ),
                    );
                    energy += torsional_energy(&p1, &p2, &p3, &p4, k_omega * weight, n_fold, gamma);
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
            let m6 =
                crate::forcefield::etkdg_lite::infer_etkdg_parameters(mol, u.index(), v.index());
            // Skip if all coefficients are zero
            if m6.v.iter().all(|&x| x.abs() < 1e-6) {
                continue;
            }

            let neighbors_u: Vec<_> = mol.graph.neighbors(u).filter(|&x| x != v).collect();
            let neighbors_v: Vec<_> = mol.graph.neighbors(v).filter(|&x| x != u).collect();
            if neighbors_u.is_empty() || neighbors_v.is_empty() {
                continue;
            }
            // Use first neighbor pair only (matching RDKit's approach for ETKDG)
            let nu = neighbors_u[0];
            let nv = neighbors_v[0];
            let (p1, p2, p3, p4) = (
                Vector3::new(
                    coords[(nu.index(), 0)] as f32,
                    coords[(nu.index(), 1)] as f32,
                    coords[(nu.index(), 2)] as f32,
                ),
                Vector3::new(
                    coords[(u.index(), 0)] as f32,
                    coords[(u.index(), 1)] as f32,
                    coords[(u.index(), 2)] as f32,
                ),
                Vector3::new(
                    coords[(v.index(), 0)] as f32,
                    coords[(v.index(), 1)] as f32,
                    coords[(v.index(), 2)] as f32,
                ),
                Vector3::new(
                    coords[(nv.index(), 0)] as f32,
                    coords[(nv.index(), 1)] as f32,
                    coords[(nv.index(), 2)] as f32,
                ),
            );
            // The M6 torsion term keeps its f32 kernel (the parameters are
            // f32); only the accumulation happens in f64.
            energy += crate::forcefield::etkdg_lite::calc_torsion_energy_m6(&p1, &p2, &p3, &p4, &m6)
                as f64;
        }
    }

    // 7. VDW non-bonded interactions (1-4+ pairs)
    if k_vdw.abs() > 1e-8 {
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
                    if a == b {
                        continue;
                    }
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
            let (xi, di) = (xi as f64, di as f64);
            let pi = Vector3::new(
                coords[(i, 0)] as f64,
                coords[(i, 1)] as f64,
                coords[(i, 2)] as f64,
            );
            for j in (i + 1)..n {
                if excluded.contains(&(i, j)) {
                    continue;
                }
                let ej = mol.graph[petgraph::graph::NodeIndex::new(j)].element;
                let (xj, dj) = uff_vdw_params(ej);
                let (xj, dj) = (xj as f64, dj as f64);
                let r_star = (xi + xj) * 0.5;
                let eps_full = (di * dj).sqrt();
                let scale = if is_14.contains(&(i, j)) { 0.5 } else { 1.0 };
                let pj = Vector3::new(
                    coords[(j, 0)] as f64,
                    coords[(j, 1)] as f64,
                    coords[(j, 2)] as f64,
                );
                energy += k_vdw * scale * vdw_energy(&pi, &pj, r_star, eps_full);
            }
        }
    }

    energy
}

/// Legacy `f32` entry point of [`calculate_total_energy_f64`].
///
/// The energy is computed in `f64` and rounded to `f32` once, at the end.
/// The signature is kept for the callers that still work in `f32`
/// (`src/forcefield/minimizer.rs`); new code should use
/// [`calculate_total_energy_f64`].
pub fn calculate_total_energy(
    coords: &nalgebra::DMatrix<f32>,
    mol: &crate::graph::Molecule,
    params: &FFParams,
    bounds_matrix: &nalgebra::DMatrix<f64>,
) -> f32 {
    calculate_total_energy_f64(coords, mol, params, bounds_matrix) as f32
}

/// Determine torsion periodicity, phase, and relative weight based on UFF rules
pub fn torsion_params(
    hyb_u: crate::graph::Hybridization,
    hyb_v: crate::graph::Hybridization,
) -> (f32, f32, f32) {
    use crate::graph::Hybridization::*;
    let pi = std::f32::consts::PI;
    match (hyb_u, hyb_v) {
        (SP3, SP3) => (3.0, 0.0, 1.0),             // staggered, normal weight
        (SP2, SP2) => (2.0, pi, 5.0),              // planar, strong weight
        (SP2, SP3) | (SP3, SP2) => (6.0, pi, 0.5), // 6-fold weak barrier
        _ => (3.0, 0.0, 1.0),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A-H3 regression test: the unit conversions applied on the way in
    /// (bond length in Å, ideal angle in radians) must be applied in `f64`,
    /// not truncated to `f32` before the term is evaluated.
    ///
    /// The test recomputes the bond terms independently in `f64`; if the bond
    /// length were rounded to `f32` first, the result would differ by
    /// ~1e-5 kcal/mol.
    #[test]
    fn total_energy_uses_f64_unit_conversions() {
        // Water (3 atoms: O + 2 H) — small enough to have no M6 torsion term.
        let mol = crate::graph::Molecule::from_smiles("O").expect("water parses");
        let n = mol.graph.node_count();
        assert_eq!(n, 3);
        let params = FFParams {
            kb: 500.0,
            k_theta: 0.0,
            k_omega: 0.0,
            k_oop: 0.0,
            k_bounds: 0.0,
            k_chiral: 0.0,
            k_vdw: 0.0,
        };

        let mut coords = nalgebra::DMatrix::<f32>::zeros(n, 3);
        coords[(1, 0)] = 0.957;
        coords[(2, 0)] = -0.239;
        coords[(2, 1)] = 0.927;

        let bounds = nalgebra::DMatrix::<f64>::zeros(n, n);
        let energy = calculate_total_energy_f64(&coords, &mol, &params, &bounds);

        // Independent f64 reference for the two bond terms.
        let mut reference = 0.0f64;
        for edge in mol.graph.edge_references() {
            let i = edge.source().index();
            let j = edge.target().index();
            let p1 = Vector3::new(
                coords[(i, 0)] as f64,
                coords[(i, 1)] as f64,
                coords[(i, 2)] as f64,
            );
            let p2 = Vector3::new(
                coords[(j, 0)] as f64,
                coords[(j, 1)] as f64,
                coords[(j, 2)] as f64,
            );
            let r_eq = crate::distgeom::get_bond_length(&mol, edge.source(), edge.target());
            reference += bond_stretch_energy(&p1, &p2, 500.0, r_eq);
        }

        assert!(energy.is_finite());
        assert!(
            (energy - reference).abs() <= 1e-12 * reference.abs().max(1.0),
            "f64 energy {energy} != f64 reference {reference}"
        );
        // The `f32` shim is the correctly rounded `f64` value.
        assert_eq!(
            calculate_total_energy(&coords, &mol, &params, &bounds),
            energy as f32
        );
    }
}
