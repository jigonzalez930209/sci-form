use nalgebra::Vector3;
use petgraph::graph::NodeIndex;
use petgraph::graph::UnGraph;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
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

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub enum BondStereo {
    None,
    E,
    Z,
    Any,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
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
pub fn get_covalent_radius(atomic_number: u8) -> f32 {
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
pub fn get_vdw_radius(atomic_number: u8) -> f32 {
    match atomic_number {
        1 => 1.20,  // H
        2 => 1.40,  // He
        3 => 1.82,  // Li
        4 => 1.53,  // Be
        5 => 1.92,  // B
        6 => 1.70,  // C
        7 => 1.55,  // N
        8 => 1.52,  // O
        9 => 1.47,  // F
        10 => 1.54, // Ne
        11 => 2.27, // Na
        12 => 1.73, // Mg
        14 => 2.10, // Si
        15 => 1.80, // P
        16 => 1.80, // S
        17 => 1.75, // Cl
        18 => 1.88, // Ar
        35 => 1.85, // Br
        53 => 1.98, // I
        _ => 2.0,
    }
}

/// Returns the ideal bond angle (in radians) for a given hybridization
pub fn get_ideal_angle(hybridization: &Hybridization) -> f32 {
    match hybridization {
        Hybridization::SP => std::f32::consts::PI, // 180 degrees
        Hybridization::SP2 => 2.0 * std::f32::consts::PI / 3.0, // 120 degrees
        Hybridization::SP3 => 109.5 * std::f32::consts::PI / 180.0, // 109.5 degrees
        Hybridization::SP3D => 105.0 * std::f32::consts::PI / 180.0, // 105 degrees (approx)
        Hybridization::SP3D2 => 90.0 * std::f32::consts::PI / 180.0, // 90 degrees
        Hybridization::Unknown => 109.5 * std::f32::consts::PI / 180.0, // Fallback to tetrahedral
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

/// Returns the ideal bond angle accounting for rings topological constraints
pub fn get_corrected_ideal_angle(
    mol: &Molecule,
    center: petgraph::graph::NodeIndex,
    n1: petgraph::graph::NodeIndex,
    n2: petgraph::graph::NodeIndex,
) -> f32 {
    use std::f32::consts::PI;
    let ahyb = &mol.graph[center].hybridization;

    let dist = min_path_excluding(mol, n1, n2, center, 6);

    match dist {
        Some(1) => {
            return 60.0 * PI / 180.0; // 3-ring
        }
        Some(2) => {
            return 90.0 * PI / 180.0; // 4-ring
        }
        Some(3) => {
            return 108.0 * PI / 180.0; // 5-ring
        }
        Some(4) => {
            return 120.0 * PI / 180.0; // 6-ring
        }
        _ => {}
    }

    let mut in_ring_3 = false;
    let mut in_ring_4 = false;
    let mut in_ring_5 = false;

    let neighbors_of_center: Vec<_> = mol.graph.neighbors(center).collect();
    for i in 0..neighbors_of_center.len() {
        for j in (i + 1)..neighbors_of_center.len() {
            let ni = neighbors_of_center[i];
            let nj = neighbors_of_center[j];
            let d = min_path_excluding(mol, ni, nj, center, 3);
            if d == Some(1) {
                in_ring_3 = true;
            } else if d == Some(2) {
                in_ring_4 = true;
            } else if d == Some(3) {
                in_ring_5 = true;
            }
        }
    }

    if ahyb == &Hybridization::SP3 {
        if in_ring_3 {
            return 116.0 * PI / 180.0;
        } else if in_ring_4 {
            return 112.0 * PI / 180.0;
        }
    } else if ahyb == &Hybridization::SP2 {
        if in_ring_5 {
            return 126.0 * PI / 180.0;
        }
    }

    get_ideal_angle(ahyb)
}
