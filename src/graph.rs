use nalgebra::Vector3;
pub use petgraph::graph::NodeIndex;
use petgraph::graph::UnGraph;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq)]
pub enum Hybridization {
    SP,
    SP2,
    SP3,
    SP3D,
    SP3D2,
    Unknown,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub enum ChiralType {
    Unspecified,
    TetrahedralCW,  // R usually
    TetrahedralCCW, // S usually
    Other,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq)]
pub enum BondStereo {
    None,
    E,
    Z,
    Any,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
    Unknown,
}

#[derive(Clone, Debug)]
pub struct Atom {
    pub element: u8,
    pub position: Vector3<f32>,
    pub charge: f32,
    pub formal_charge: i8,
    pub hybridization: Hybridization,
    pub chiral_tag: ChiralType,
    pub explicit_h: u8,
}

#[derive(Clone, Debug)]
pub struct Bond {
    pub order: BondOrder,
    pub stereo: BondStereo,
}

impl Atom {
    pub fn new(atomic_number: u8, x: f32, y: f32, z: f32) -> Self {
        Atom {
            element: atomic_number,
            position: Vector3::new(x, y, z),
            charge: 0.0,
            formal_charge: 0,
            hybridization: Hybridization::Unknown,
            chiral_tag: ChiralType::Unspecified,
            explicit_h: 0,
        }
    }
}

// Representing Molecule utilizing petgraph
pub struct Molecule {
    pub graph: UnGraph<Atom, Bond>,
    pub name: String,
}

impl Molecule {
    pub fn new(name: &str) -> Self {
        Self {
            graph: UnGraph::new_undirected(),
            name: name.to_string(),
        }
    }

    pub fn from_smiles(smiles: &str) -> Result<Self, String> {
        let mut mol = Self::new(smiles);
        let mut parser = crate::smiles::SmilesParser::new(smiles, &mut mol);
        parser.parse()?;
        Ok(mol)
    }

    pub fn add_atom(&mut self, atom: Atom) -> NodeIndex {
        self.graph.add_node(atom)
    }

    pub fn add_bond(&mut self, a: NodeIndex, b: NodeIndex, bond: Bond) {
        self.graph.add_edge(a, b, bond);
    }
}

// =============== ATOMIC CONSTANTS ==================

/// Returns the covalent radius (in Å) for a given atomic number
/// Data derived from RDKit's atomic_data.cpp
pub fn get_covalent_radius(atomic_number: u8) -> f64 {
    match atomic_number {
        1 => 0.31,  // H
        2 => 0.28,  // He
        3 => 1.28,  // Li
        4 => 0.96,  // Be
        5 => 0.84,  // B
        6 => 0.76,  // C
        7 => 0.71,  // N
        8 => 0.66,  // O
        9 => 0.57,  // F
        10 => 0.58, // Ne
        11 => 1.66, // Na
        12 => 1.41, // Mg
        13 => 1.21, // Al
        14 => 1.11, // Si
        15 => 1.07, // P
        16 => 1.05, // S
        17 => 1.02, // Cl
        18 => 1.06, // Ar
        19 => 1.96, // K
        20 => 1.71, // Ca
        21 => 1.48, // Sc
        22 => 1.36, // Ti
        23 => 1.34, // V
        24 => 1.22, // Cr
        25 => 1.19, // Mn
        26 => 1.16, // Fe
        27 => 1.11, // Co
        28 => 1.10, // Ni
        29 => 1.12, // Cu
        30 => 1.18, // Zn
        31 => 1.24, // Ga
        32 => 1.21, // Ge
        33 => 1.21, // As
        34 => 1.16, // Se
        35 => 1.20, // Br
        36 => 1.16, // Kr
        37 => 2.10, // Rb
        38 => 1.85, // Sr
        39 => 1.63, // Y
        40 => 1.54, // Zr
        41 => 1.47, // Nb
        42 => 1.38, // Mo
        43 => 1.28, // Tc
        44 => 1.25, // Ru
        45 => 1.25, // Rh
        46 => 1.20, // Pd
        47 => 1.28, // Ag
        48 => 1.36, // Cd
        49 => 1.42, // In
        50 => 1.40, // Sn
        51 => 1.40, // Sb
        52 => 1.36, // Te
        53 => 1.39, // I
        54 => 1.40, // Xe
        55 => 2.32, // Cs
        56 => 1.96, // Ba
        72 => 1.50, // Hf
        73 => 1.38, // Ta
        74 => 1.37, // W
        75 => 1.35, // Re
        76 => 1.26, // Os
        77 => 1.27, // Ir
        78 => 1.23, // Pt
        79 => 1.24, // Au
        80 => 1.33, // Hg
        81 => 1.44, // Tl
        82 => 1.46, // Pb
        83 => 1.48, // Bi
        84 => 1.46, // Po
        85 => 1.45, // At
        86 => 1.50, // Rn
        _ => 0.76,  // Fallback to Carbon's radius if unknown mostly
    }
}

/// Returns the Van der Waals radius (in Å)
pub fn get_vdw_radius(atomic_number: u8) -> f64 {
    // RDKit's PeriodicTable::getRvdw values
    match atomic_number {
        1 => 1.20,  // H
        2 => 1.40,  // He
        3 => 1.82,  // Li
        4 => 1.53,  // Be
        5 => 1.92,  // B
        6 => 1.70,  // C
        7 => 1.60,  // N (was 1.55)
        8 => 1.55,  // O (was 1.52)
        9 => 1.50,  // F (was 1.47)
        10 => 1.54, // Ne
        11 => 2.27, // Na
        12 => 1.73, // Mg
        14 => 2.10, // Si
        15 => 1.80, // P
        16 => 1.80, // S
        17 => 1.80, // Cl (was 1.75)
        18 => 1.88, // Ar
        35 => 1.90, // Br (was 1.85)
        53 => 2.10, // I (was 1.98)
        _ => 2.0,
    }
}

/// Returns the ideal bond angle (in radians) for a given hybridization
pub fn get_ideal_angle(hybridization: &Hybridization) -> f64 {
    match hybridization {
        Hybridization::SP => std::f64::consts::PI, // 180 degrees
        Hybridization::SP2 => 2.0 * std::f64::consts::PI / 3.0, // 120 degrees
        Hybridization::SP3 => 109.5 * std::f64::consts::PI / 180.0, // 109.5 degrees
        Hybridization::SP3D => 105.0 * std::f64::consts::PI / 180.0, // 105 degrees (approx)
        Hybridization::SP3D2 => 90.0 * std::f64::consts::PI / 180.0, // 90 degrees
        Hybridization::Unknown => 109.5 * std::f64::consts::PI / 180.0, // Fallback to tetrahedral
    }
}

/// Returns the ideal bond angle accounting for 3-, 4-, and 5-membered rings topological constraints
pub fn min_path_excluding(
    mol: &Molecule,
    start: petgraph::graph::NodeIndex,
    target: petgraph::graph::NodeIndex,
    exclude: petgraph::graph::NodeIndex,
    limit: usize,
) -> Option<usize> {
    if mol.graph.contains_edge(start, target) {
        return Some(1);
    }
    let mut queue = std::collections::VecDeque::new();
    let mut visited = std::collections::HashSet::new();
    queue.push_back((start, 0));
    visited.insert(start);

    while let Some((curr, dist)) = queue.pop_front() {
        if dist >= limit {
            continue;
        }
        for n in mol.graph.neighbors(curr) {
            if n == target {
                return Some(dist + 1);
            }
            if n != exclude && !visited.contains(&n) {
                visited.insert(n);
                queue.push_back((n, dist + 1));
            }
        }
    }
    None
}

/// BFS shortest path from start to target, excluding two intermediate nodes.
pub fn min_path_excluding2(
    mol: &Molecule,
    start: petgraph::graph::NodeIndex,
    target: petgraph::graph::NodeIndex,
    excl1: petgraph::graph::NodeIndex,
    excl2: petgraph::graph::NodeIndex,
    limit: usize,
) -> Option<usize> {
    let mut queue = std::collections::VecDeque::new();
    let mut visited = std::collections::HashSet::new();
    visited.insert(excl1);
    visited.insert(excl2);
    visited.insert(start);

    // Always start from neighbors, skipping target to avoid counting the direct bond
    for n in mol.graph.neighbors(start) {
        if n == target {
            continue; // Skip direct bond — we want an alternative path
        }
        if !visited.contains(&n) {
            visited.insert(n);
            queue.push_back((n, 1usize));
        }
    }

    while let Some((curr, dist)) = queue.pop_front() {
        if dist >= limit {
            continue;
        }
        for n in mol.graph.neighbors(curr) {
            if n == target {
                return Some(dist + 1);
            }
            if !visited.contains(&n) {
                visited.insert(n);
                queue.push_back((n, dist + 1));
            }
        }
    }
    None
}

/// Returns the ideal bond angle accounting for rings topological constraints
/// Helper: compute the ring interior angle for a given hybridization and ring size.
/// Matches RDKit's _setRingAngle logic.
fn ring_angle_for(ahyb: &Hybridization, ring_size: usize) -> f64 {
    use std::f64::consts::PI;
    if ring_size == 3 || ring_size == 4 || (*ahyb == Hybridization::SP2 && ring_size <= 8) {
        PI * (1.0 - 2.0 / ring_size as f64)
    } else if *ahyb == Hybridization::SP3 {
        if ring_size == 5 {
            104.0 * PI / 180.0
        } else {
            109.5 * PI / 180.0
        }
    } else {
        get_ideal_angle(ahyb)
    }
}

pub fn get_corrected_ideal_angle(
    mol: &Molecule,
    center: petgraph::graph::NodeIndex,
    n1: petgraph::graph::NodeIndex,
    n2: petgraph::graph::NodeIndex,
) -> f64 {
    use std::f64::consts::PI;
    let ahyb = &mol.graph[center].hybridization;

    // Check if n1 and n2 are in a ring through center
    let dist = min_path_excluding(mol, n1, n2, center, 8);

    match dist {
        Some(1) => return 60.0 * PI / 180.0, // 3-ring
        Some(2) => return 90.0 * PI / 180.0, // 4-ring
        _ => {}
    }

    // RDKit special case: exocyclic angles at SP3 atoms in small rings.
    // When n1-n2 are NOT in a 3- or 4-ring through center, but center IS in a 3- or 4-ring,
    // the exocyclic angle is widened: 116° for 3-ring, 112° for 4-ring.
    // dist > 2 means the pair is not directly in any 3-ring (dist=1) or 4-ring (dist=2),
    // so it's either in a larger ring envelope (fused system) or truly exocyclic.
    if *ahyb == Hybridization::SP3 && dist.map_or(true, |d| d > 2) {
        let nbs: Vec<_> = mol.graph.neighbors(center).collect();
        let mut smallest_ring = 0usize;
        for i in 0..nbs.len() {
            for j in (i + 1)..nbs.len() {
                if let Some(p) = min_path_excluding(mol, nbs[i], nbs[j], center, 8) {
                    let rs = p + 2;
                    if rs == 3 {
                        smallest_ring = if smallest_ring == 0 { 3 } else { smallest_ring.min(3) };
                    } else if rs == 4 && smallest_ring != 3 {
                        smallest_ring = if smallest_ring == 0 { 4 } else { smallest_ring.min(4) };
                    }
                }
            }
        }
        if smallest_ring == 3 {
            return 116.0 * PI / 180.0;
        } else if smallest_ring == 4 {
            return 112.0 * PI / 180.0;
        }
    }

    // For SP2 atoms with exactly 3 neighbors: use RDKit's approach of
    // distributing 360° among all neighbor pairs based on ring membership.
    // angleTaken = sum of ring angles, exo = (2π - angleTaken) / remaining_pairs
    if *ahyb == Hybridization::SP2 {
        let nbs: Vec<_> = mol.graph.neighbors(center).collect();
        if nbs.len() == 3 {
            // Find ring size for each pair
            let pairs = [(0, 1, 2), (0, 2, 1), (1, 2, 0)];
            let mut ring_angles = [0.0f64; 3]; // angles for pairs (0,1), (0,2), (1,2)
            let mut in_ring = [false; 3];

            for (idx, &(a, b, _)) in pairs.iter().enumerate() {
                let d = min_path_excluding(mol, nbs[a], nbs[b], center, 8);
                if let Some(p) = d {
                    let rs = p + 2;
                    if rs <= 8 {
                        ring_angles[idx] = ring_angle_for(ahyb, rs);
                        in_ring[idx] = true;
                    }
                }
            }

            let ring_count = in_ring.iter().filter(|&&x| x).count();
            // Identify which pair index corresponds to (n1, n2)
            let target_idx = if (nbs[0] == n1 && nbs[1] == n2) || (nbs[1] == n1 && nbs[0] == n2)
            {
                0
            } else if (nbs[0] == n1 && nbs[2] == n2) || (nbs[2] == n1 && nbs[0] == n2) {
                1
            } else {
                2
            };

            if in_ring[target_idx] && ring_count <= 1 {
                // Only this pair in a ring — return ring angle directly
                return ring_angles[target_idx];
            } else if !in_ring[target_idx] && ring_count > 0 {
                // This pair is NOT in a ring, but other pairs are:
                // exo angle = (2π - sum_of_ring_angles) / non_ring_pairs
                let angle_taken: f64 = ring_angles.iter().sum();
                let non_ring = 3 - ring_count;
                if non_ring > 0 {
                    let exo = (2.0 * PI - angle_taken) / non_ring as f64;
                    if exo > 0.0 {
                        return exo;
                    }
                }
            } else if ring_count >= 2 && in_ring[target_idx] {
                // This pair AND other pairs are in rings.
                // For fused junctions: check if this pair's "ring" is
                // actually the fused envelope. If the other two pairs
                // account for the full angle, use ring angles directly.
                let other_angle: f64 = (0..3)
                    .filter(|&i| i != target_idx && in_ring[i])
                    .map(|i| ring_angles[i])
                    .sum();
                let this_ring_expected = 2.0 * PI - other_angle;
                // If this pair's ring angle differs significantly from expected exo,
                // it's a fused envelope — use exo angle
                if (ring_angles[target_idx] - this_ring_expected).abs() > 0.05 {
                    return this_ring_expected;
                }
                return ring_angles[target_idx];
            }
        }
    }

    // Non-SP2 or fallback: use ring angle if in a ring
    if let Some(ring_path) = dist {
        let ring_size = ring_path + 2;
        if ring_size <= 8 {
            return ring_angle_for(ahyb, ring_size);
        }
    }

    get_ideal_angle(ahyb)
}
