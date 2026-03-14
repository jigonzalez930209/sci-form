use crate::graph::{Hybridization, Molecule};
use nalgebra::Vector3;

#[derive(Debug, Clone)]
pub struct M6Params {
    pub v: [f32; 6],
    pub s: [f32; 6],
}

impl Default for M6Params {
    fn default() -> Self {
        M6Params {
            v: [0.0; 6],
            s: [1.0; 6],
        }
    }
}

/// Calculates the RDKit M6 Empirical Torsion Energy:
/// V = V1*(1 + s1*cos(phi)) + V2*(1 + s2*cos(2*phi)) + ... + V6*(1 + s6*cos(6*phi))
pub fn calc_torsion_energy_m6(
    p1: &Vector3<f32>,
    p2: &Vector3<f32>,
    p3: &Vector3<f32>,
    p4: &Vector3<f32>,
    params: &M6Params,
) -> f32 {
    let r1 = p1 - p2;
    let r2 = p3 - p2;
    let r3 = p2 - p3;
    let r4 = p4 - p3;

    let t1 = r1.cross(&r2);
    let t2 = r3.cross(&r4);

    let d1 = t1.norm();
    let d2 = t2.norm();

    if d1 < 1e-4 || d2 < 1e-4 {
        return 0.0;
    }

    let n1 = t1 / d1;
    let n2 = t2 / d2;

    let mut cos_phi = n1.dot(&n2);
    cos_phi = cos_phi.clamp(-1.0, 1.0);

    // Using Chebyshev polynomials to compute cos(n*x) from cos(x) without trig functions
    let cos_phi2 = cos_phi * cos_phi;
    let cos_phi3 = cos_phi * cos_phi2;
    let cos_phi4 = cos_phi * cos_phi3;
    let cos_phi5 = cos_phi * cos_phi4;
    let cos_phi6 = cos_phi * cos_phi5;

    let cos2_phi = 2.0 * cos_phi2 - 1.0;
    let cos3_phi = 4.0 * cos_phi3 - 3.0 * cos_phi;
    let cos4_phi = 8.0 * cos_phi4 - 8.0 * cos_phi2 + 1.0;
    let cos5_phi = 16.0 * cos_phi5 - 20.0 * cos_phi3 + 5.0 * cos_phi;
    let cos6_phi = 32.0 * cos_phi6 - 48.0 * cos_phi4 + 18.0 * cos_phi2 - 1.0;

    let mut energy = 0.0;
    energy += params.v[0] * (1.0 + params.s[0] * cos_phi);
    energy += params.v[1] * (1.0 + params.s[1] * cos2_phi);
    energy += params.v[2] * (1.0 + params.s[2] * cos3_phi);
    energy += params.v[3] * (1.0 + params.s[3] * cos4_phi);
    energy += params.v[4] * (1.0 + params.s[4] * cos5_phi);
    energy += params.v[5] * (1.0 + params.s[5] * cos6_phi);

    energy
}

/// Computes the analytical gradient for the M6 torsion potential.
/// Matches RDKit's TorsionAngleM6::getGrad dE/dPhi computation exactly, including:
/// - Chebyshev polynomial derivatives for dE/dPhi
/// - RDKit's copy-paste bug: 6th harmonic uses V[4]/s[4] instead of V[5]/s[5]
///
/// Uses Blondel & Karplus for gradient distribution (mathematically equivalent to RDKit's Wilson B-matrix).
pub fn calc_torsion_grad_m6(
    p1: &Vector3<f32>,
    p2: &Vector3<f32>,
    p3: &Vector3<f32>,
    p4: &Vector3<f32>,
    params: &M6Params,
    grad: &mut nalgebra::DMatrix<f32>,
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
) {
    // RDKit-style dihedral computation: r[0]=p1-p2, r[1]=p3-p2, r[2]=-r[1], r[3]=p4-p3
    let r0 = p1 - p2;
    let r1 = p3 - p2;
    let r2 = -&r1;
    let r3 = p4 - p3;

    let mut t0 = r0.cross(&r1);
    let d0 = t0.norm().max(1e-5);
    t0 /= d0;
    let mut t1 = r2.cross(&r3);
    let d1 = t1.norm().max(1e-5);
    t1 /= d1;

    let cos_phi = t0.dot(&t1).clamp(-1.0, 1.0);
    let sin_phi_sq = 1.0 - cos_phi * cos_phi;
    let sin_phi = if sin_phi_sq > 0.0 {
        sin_phi_sq.sqrt()
    } else {
        0.0
    };

    let cos_phi2 = cos_phi * cos_phi;
    let cos_phi3 = cos_phi * cos_phi2;
    let cos_phi4 = cos_phi * cos_phi3;
    let cos_phi5 = cos_phi * cos_phi4;

    // dE/dPhi using Chebyshev derivatives, matching RDKit exactly.
    // NOTE: 6th term intentionally uses V[4]/s[4] to match RDKit's copy-paste bug
    let de_dphi = -params.v[0] * params.s[0] * sin_phi
        - 2.0 * params.v[1] * params.s[1] * (2.0 * cos_phi * sin_phi)
        - 3.0 * params.v[2] * params.s[2] * (4.0 * cos_phi2 * sin_phi - sin_phi)
        - 4.0 * params.v[3] * params.s[3] * (8.0 * cos_phi3 * sin_phi - 4.0 * cos_phi * sin_phi)
        - 5.0
            * params.v[4]
            * params.s[4]
            * (16.0 * cos_phi4 * sin_phi - 12.0 * cos_phi2 * sin_phi + sin_phi)
        - 6.0
            * params.v[4]
            * params.s[4]
            * (32.0 * cos_phi5 * sin_phi - 32.0 * cos_phi3 * sin_phi + 6.0 * cos_phi * sin_phi);

    // sinTerm = -dE_dPhi / sinPhi (or 1/cosPhi when sinPhi≈0), matching RDKit
    let is_zero_sin = sin_phi < 1e-5;
    let sin_term = -de_dphi
        * if is_zero_sin {
            1.0 / cos_phi
        } else {
            1.0 / sin_phi
        };

    // RDKit calcTorsionGrad: dCosPhi/dT vectors
    let d_cos_dt = [
        (t1.x - cos_phi * t0.x) / d0,
        (t1.y - cos_phi * t0.y) / d0,
        (t1.z - cos_phi * t0.z) / d0,
        (t0.x - cos_phi * t1.x) / d1,
        (t0.y - cos_phi * t1.y) / d1,
        (t0.z - cos_phi * t1.z) / d1,
    ];

    // Atom 1 (idx1) gradient
    grad[(idx1, 0)] += sin_term * (d_cos_dt[2] * r1.y - d_cos_dt[1] * r1.z);
    grad[(idx1, 1)] += sin_term * (d_cos_dt[0] * r1.z - d_cos_dt[2] * r1.x);
    grad[(idx1, 2)] += sin_term * (d_cos_dt[1] * r1.x - d_cos_dt[0] * r1.y);

    // Atom 2 (idx2) gradient
    grad[(idx2, 0)] += sin_term
        * (d_cos_dt[1] * (r1.z - r0.z)
            + d_cos_dt[2] * (r0.y - r1.y)
            + d_cos_dt[4] * (-r3.z)
            + d_cos_dt[5] * r3.y);
    grad[(idx2, 1)] += sin_term
        * (d_cos_dt[0] * (r0.z - r1.z)
            + d_cos_dt[2] * (r1.x - r0.x)
            + d_cos_dt[3] * r3.z
            + d_cos_dt[5] * (-r3.x));
    grad[(idx2, 2)] += sin_term
        * (d_cos_dt[0] * (r1.y - r0.y)
            + d_cos_dt[1] * (r0.x - r1.x)
            + d_cos_dt[3] * (-r3.y)
            + d_cos_dt[4] * r3.x);

    // Atom 3 (idx3) gradient
    grad[(idx3, 0)] += sin_term
        * (d_cos_dt[1] * r0.z
            + d_cos_dt[2] * (-r0.y)
            + d_cos_dt[4] * (r3.z - r2.z)
            + d_cos_dt[5] * (r2.y - r3.y));
    grad[(idx3, 1)] += sin_term
        * (d_cos_dt[0] * (-r0.z)
            + d_cos_dt[2] * r0.x
            + d_cos_dt[3] * (r2.z - r3.z)
            + d_cos_dt[5] * (r3.x - r2.x));
    grad[(idx3, 2)] += sin_term
        * (d_cos_dt[0] * r0.y
            + d_cos_dt[1] * (-r0.x)
            + d_cos_dt[3] * (r3.y - r2.y)
            + d_cos_dt[4] * (r2.x - r3.x));

    // Atom 4 (idx4) gradient
    grad[(idx4, 0)] += sin_term * (d_cos_dt[4] * r2.z - d_cos_dt[5] * r2.y);
    grad[(idx4, 1)] += sin_term * (d_cos_dt[5] * r2.x - d_cos_dt[3] * r2.z);
    grad[(idx4, 2)] += sin_term * (d_cos_dt[3] * r2.y - d_cos_dt[4] * r2.x);
}

fn is_partial_double_bond(mol: &Molecule, n2_idx: usize, n3_idx: usize) -> bool {
    let ni2 = petgraph::graph::NodeIndex::new(n2_idx);
    let ni3 = petgraph::graph::NodeIndex::new(n3_idx);

    // If it's already a double bond, it's planar
    if let Some(edge) = mol.graph.find_edge(ni2, ni3) {
        if mol.graph[edge].order == crate::graph::BondOrder::Double {
            return true;
        }
        if mol.graph[edge].order == crate::graph::BondOrder::Aromatic {
            return true;
        }
    }

    // Check for Amides (N-C=O)
    // Or Esters (O-C=O)
    for ni in &[ni2, ni3] {
        let other = if *ni == ni2 { ni3 } else { ni2 };
        let atom = &mol.graph[*ni];
        let other_atom = &mol.graph[other];

        if (atom.element == 7 || atom.element == 8)
            && other_atom.hybridization == Hybridization::SP2
        {
            // Check if 'other' has a double bond to O or N
            for neighbor in mol.graph.neighbors(other) {
                if neighbor == *ni {
                    continue;
                }
                if let Some(edge) = mol.graph.find_edge(other, neighbor) {
                    if mol.graph[edge].order == crate::graph::BondOrder::Double {
                        let neigh_atom = &mol.graph[neighbor];
                        if neigh_atom.element == 8 || neigh_atom.element == 7 {
                            return true;
                        }
                    }
                }
            }
        }
    }

    false
}

/// Infers a reasonable ETKDG experimental torsion profile based on hybridization for cases
/// where the 50-conformer pool has been turned off, helping the single conformer settle into
/// a biologically relevant rotamer.
pub fn infer_etkdg_parameters(mol: &Molecule, n2_idx: usize, n3_idx: usize) -> M6Params {
    let mut params = M6Params::default();

    if is_partial_double_bond(mol, n2_idx, n3_idx) {
        // V2 enforces planarity: minima at both 0° and 180° (cis and trans)
        // This allows the geometry to settle into whichever orientation is
        // energetically favorable, while preventing non-planar conformations.
        params.v[1] = 100.0;
        params.s[1] = -1.0; // Minima at 0° (cis) and 180° (trans)
        return params;
    }

    let a2 = &mol.graph[petgraph::graph::NodeIndex::new(n2_idx)];
    let a3 = &mol.graph[petgraph::graph::NodeIndex::new(n3_idx)];

    let e2 = a2.element;
    let e3 = a3.element;
    let h2 = a2.hybridization;
    let h3 = a3.hybridization;

    // SP2 - SP2 (e.g., conjugated single bonds)
    if h2 == Hybridization::SP2 && h3 == Hybridization::SP2 {
        params.v[1] = 10.0; // V2
        params.s[1] = -1.0; // Minima at 0, 180
    }
    // C(sp3) - O - C bonds (ether): V3 staggered, matching RDKit's [*:1][CX4:2]-[O:3][CX4:4] -> V3=2.5-4.0
    else if (e2 == 6 && e3 == 8 && h2 == Hybridization::SP3)
        || (e2 == 8 && e3 == 6 && h3 == Hybridization::SP3)
    {
        params.v[2] = 4.0;
        params.s[2] = 1.0;
    }
    // SP3 - SP3 (e.g., butane chain)
    else if h2 == Hybridization::SP3 && h3 == Hybridization::SP3 {
        params.v[2] = 5.0; // Staggered
        params.s[2] = 1.0;
        params.v[0] = 2.0; // Anti preference
        params.s[0] = 1.0;
    }
    // Default fallback
    else {
        params.v[2] = 3.0;
        params.s[2] = 1.0;
    }

    params
}
