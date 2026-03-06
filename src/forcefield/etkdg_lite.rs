use crate::graph::{Molecule, Hybridization};
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
pub fn calc_torsion_grad_m6(
    p1: &Vector3<f32>,
    p2: &Vector3<f32>,
    p3: &Vector3<f32>,
    p4: &Vector3<f32>,
    params: &M6Params,
    grad: &mut nalgebra::DMatrix<f32>,
    idx1: usize, idx2: usize, idx3: usize, idx4: usize,
) {
    let r12 = p1 - p2;
    let r23 = p3 - p2;
    let r34 = p4 - p3;

    let t1 = r12.cross(&r23);
    let t2 = r23.cross(&r34);

    let d1 = t1.norm();
    let d2 = t2.norm();
    let r23m = r23.norm();

    if d1 < 1e-4 || d2 < 1e-4 || r23m < 1e-4 {
        return;
    }

    let n1 = t1 / (d1 * d1);
    let n2 = t2 / (d2 * d2);

    let cos_phi = (t1.dot(&t2) / (d1 * d2)).clamp(-1.0, 1.0);
    let phi = cos_phi.acos();
    
    // dE/dphi
    let _sin_phi = (1.0 - cos_phi * cos_phi).sqrt().max(1e-4);
    let mut de_dphi = 0.0;
    
    // V = Vn * (1 + sn * cos(n * phi))
    // dV/dphi = -Vn * sn * n * sin(n * phi)
    for n in 1..=6 {
        let nf = n as f32;
        de_dphi -= params.v[n-1] * params.s[n-1] * nf * (nf * phi).sin();
    }

    // Standard torsion gradient vectors (from Case et al. or RDKit)
    let g1 = de_dphi * r23m * n1;
    let g4 = -de_dphi * r23m * n2;
    
    let g2 = ((r12.dot(&r23) / (r23m * r23m)) - 1.0) * g1 - (r34.dot(&r23) / (r23m * r23m)) * g4;
    let g3 = -(g1 + g2 + g4);

    for d in 0..3 {
        grad[(idx1, d)] += g1[d];
        grad[(idx2, d)] += g2[d];
        grad[(idx3, d)] += g3[d];
        grad[(idx4, d)] += g4[d];
    }
}

fn is_partial_double_bond(mol: &Molecule, n2_idx: usize, n3_idx: usize) -> bool {
    let ni2 = petgraph::graph::NodeIndex::new(n2_idx);
    let ni3 = petgraph::graph::NodeIndex::new(n3_idx);

    // If it's already a double bond, it's planar
    if let Some(edge) = mol.graph.find_edge(ni2, ni3) {
        if mol.graph[edge].order == crate::graph::BondOrder::Double { return true; }
        if mol.graph[edge].order == crate::graph::BondOrder::Aromatic { return true; }
    }

    // Check for Amides (N-C=O)
    // Or Esters (O-C=O)
    for ni in &[ni2, ni3] {
        let other = if *ni == ni2 { ni3 } else { ni2 };
        let atom = &mol.graph[*ni];
        let other_atom = &mol.graph[other];

        if (atom.element == 7 || atom.element == 8) && other_atom.hybridization == Hybridization::SP2 {
            // Check if 'other' has a double bond to O or N
            for neighbor in mol.graph.neighbors(other) {
                if neighbor == *ni { continue; }
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
        // High magnitude V1/V2 to enforce planarity in amides/esters
        // Matching RDKit v2 patterns like [O:1]=[C:2]!@;-[O:3]~[C:4] -> V1=100.0
        params.v[0] = 100.0; 
        params.s[0] = -1.0; // Minima at 180 (Trans)
        return params;
    }

    let a2 = &mol.graph[petgraph::graph::NodeIndex::new(n2_idx)];
    let a3 = &mol.graph[petgraph::graph::NodeIndex::new(n3_idx)];

    // SP2 - SP2 (e.g., conjugated single bonds)
    if a2.hybridization == Hybridization::SP2 && a3.hybridization == Hybridization::SP2 {
        // Matching RDKit v2 patterns like [cH1,n:1][c:2]!@;-[O:3][CH1:4] -> V1=7.2
        params.v[1] = 10.0; // V2
        params.s[1] = -1.0; // Minima at 0, 180
    }
    // SP3 - SP3 (e.g., butane chain)
    else if a2.hybridization == Hybridization::SP3 && a3.hybridization == Hybridization::SP3 {
        // V3 term (Staggered) - RDKit default is ~2.5 but pooled.
        // For single-shot, we use slightly stronger staggered guidance.
        params.v[2] = 5.0; 
        params.s[2] = 1.0;
        // V1 term (Anti preference)
        params.v[0] = 2.0; 
        params.s[0] = 1.0;
    }
    // Default fallback
    else {
        // Generic 3-fold potential
        params.v[2] = 3.0; 
        params.s[2] = 1.0;
    }

    params
}
