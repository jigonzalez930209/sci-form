//! Framework assembly: connecting SBUs into periodic crystal structures.
//!
//! The assembly process:
//! 1. Define metal-node SBUs and organic-linker SBUs
//! 2. Specify a target topology (e.g., cubic pcu for MOF-5)
//! 3. Place nodes at predefined positions within the unit cell
//! 4. Connect them via linkers along edges of the topology net
//! 5. Output the assembled periodic structure

use super::cell::UnitCell;
use super::sbu::Sbu;
use serde::{Deserialize, Serialize};

/// A placed atom in the assembled crystal structure.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CrystalAtom {
    /// Atomic number.
    pub element: u8,
    /// Fractional coordinates in the unit cell.
    pub frac_coords: [f64; 3],
}

/// An assembled periodic crystal structure.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CrystalStructure {
    /// Unit cell definition.
    pub cell: UnitCell,
    /// All atoms in the asymmetric unit (fractional coordinates).
    pub atoms: Vec<CrystalAtom>,
    /// Labels for tracking which SBU each atom came from.
    pub labels: Vec<String>,
}

/// Topology definition: positions and edges for a network.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Topology {
    /// Name of the topology (e.g., "pcu", "dia", "sql").
    pub name: String,
    /// Node positions in fractional coordinates.
    pub nodes: Vec<[f64; 3]>,
    /// Edges as (node_i, node_j, cell_offset_j).
    /// cell_offset_j = [da, db, dc] for the periodic image of node_j.
    pub edges: Vec<(usize, usize, [i32; 3])>,
}

impl Topology {
    /// Primitive cubic (pcu) topology: single node at origin with 6 edges.
    /// Common for MOF-5 type structures.
    pub fn pcu() -> Self {
        Self {
            name: "pcu".to_string(),
            nodes: vec![[0.0, 0.0, 0.0]],
            edges: vec![
                (0, 0, [1, 0, 0]), // +a
                (0, 0, [0, 1, 0]), // +b
                (0, 0, [0, 0, 1]), // +c
            ],
        }
    }

    /// Diamond (dia) topology: two nodes with tetrahedral connectivity.
    pub fn dia() -> Self {
        Self {
            name: "dia".to_string(),
            nodes: vec![[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]],
            edges: vec![
                (0, 1, [0, 0, 0]),
                (0, 1, [-1, 0, 0]),
                (0, 1, [0, -1, 0]),
                (0, 1, [0, 0, -1]),
            ],
        }
    }

    /// Square lattice (sql) for 2D frameworks.
    pub fn sql() -> Self {
        Self {
            name: "sql".to_string(),
            nodes: vec![[0.0, 0.0, 0.0]],
            edges: vec![(0, 0, [1, 0, 0]), (0, 0, [0, 1, 0])],
        }
    }
}

impl CrystalStructure {
    /// Number of atoms in the structure.
    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// Get all Cartesian coordinates (Å).
    pub fn cartesian_coords(&self) -> Vec<[f64; 3]> {
        self.atoms
            .iter()
            .map(|a| self.cell.frac_to_cart(a.frac_coords))
            .collect()
    }

    /// Get all atomic numbers.
    pub fn elements(&self) -> Vec<u8> {
        self.atoms.iter().map(|a| a.element).collect()
    }

    /// Generate a supercell by replicating the structure.
    pub fn make_supercell(&self, na: usize, nb: usize, nc: usize) -> CrystalStructure {
        let (new_cell, offsets) = self.cell.supercell(na, nb, nc);
        let mut new_atoms = Vec::new();
        let mut new_labels = Vec::new();

        let na_f = na as f64;
        let nb_f = nb as f64;
        let nc_f = nc as f64;

        for offset in &offsets {
            for (i, atom) in self.atoms.iter().enumerate() {
                new_atoms.push(CrystalAtom {
                    element: atom.element,
                    frac_coords: [
                        (atom.frac_coords[0] + offset[0]) / na_f,
                        (atom.frac_coords[1] + offset[1]) / nb_f,
                        (atom.frac_coords[2] + offset[2]) / nc_f,
                    ],
                });
                new_labels.push(self.labels[i].clone());
            }
        }

        CrystalStructure {
            cell: new_cell,
            atoms: new_atoms,
            labels: new_labels,
        }
    }
}

/// Assemble a framework by placing SBUs on a topology.
///
/// `node_sbu`: the SBU placed at each topology node
/// `linker_sbu`: the SBU placed along each topology edge
/// `topology`: the periodic net topology
/// `cell`: the unit cell for the framework
///
/// Returns the assembled crystal structure.
pub fn assemble_framework(
    node_sbu: &Sbu,
    linker_sbu: &Sbu,
    topology: &Topology,
    cell: &UnitCell,
) -> CrystalStructure {
    let mut atoms = Vec::new();
    let mut labels = Vec::new();

    // 1. Place node SBUs at each topology node position
    for (ni, node_frac) in topology.nodes.iter().enumerate() {
        let node_cart = cell.frac_to_cart(*node_frac);

        for (ai, &elem) in node_sbu.elements.iter().enumerate() {
            let atom_cart = [
                node_cart[0] + node_sbu.positions[ai][0],
                node_cart[1] + node_sbu.positions[ai][1],
                node_cart[2] + node_sbu.positions[ai][2],
            ];
            let frac = cell.cart_to_frac(atom_cart);
            atoms.push(CrystalAtom {
                element: elem,
                frac_coords: frac,
            });
            labels.push(format!("node_{}", ni));
        }
    }

    // 2. Place linker SBUs along each topology edge
    for (ei, &(ni, nj, offset)) in topology.edges.iter().enumerate() {
        let node_i_frac = topology.nodes[ni];
        let node_j_frac = [
            topology.nodes[nj][0] + offset[0] as f64,
            topology.nodes[nj][1] + offset[1] as f64,
            topology.nodes[nj][2] + offset[2] as f64,
        ];

        let node_i_cart = cell.frac_to_cart(node_i_frac);
        let node_j_cart = cell.frac_to_cart(node_j_frac);

        // Midpoint
        let mid = [
            (node_i_cart[0] + node_j_cart[0]) / 2.0,
            (node_i_cart[1] + node_j_cart[1]) / 2.0,
            (node_i_cart[2] + node_j_cart[2]) / 2.0,
        ];

        // Direction from node_i to node_j
        let dx = node_j_cart[0] - node_i_cart[0];
        let dy = node_j_cart[1] - node_i_cart[1];
        let dz = node_j_cart[2] - node_i_cart[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

        if dist < 1e-12 {
            continue;
        }

        let dir = [dx / dist, dy / dist, dz / dist];

        // Rotate linker atoms: default linker is along x-axis,
        // rotate to align with the edge direction
        let rot = rotation_from_x_to(dir);

        for (ai, &elem) in linker_sbu.elements.iter().enumerate() {
            let p = linker_sbu.positions[ai];
            // Apply rotation then translate to midpoint
            let rotated = apply_rotation(&rot, p);
            let atom_cart = [
                mid[0] + rotated[0],
                mid[1] + rotated[1],
                mid[2] + rotated[2],
            ];
            let frac = cell.cart_to_frac(atom_cart);
            atoms.push(CrystalAtom {
                element: elem,
                frac_coords: frac,
            });
            labels.push(format!("linker_{}", ei));
        }
    }

    CrystalStructure {
        cell: cell.clone(),
        atoms,
        labels,
    }
}

/// Compute the rotation matrix that maps [1,0,0] to the target direction.
fn rotation_from_x_to(target: [f64; 3]) -> [[f64; 3]; 3] {
    let x = [1.0, 0.0, 0.0];
    let dot = x[0] * target[0] + x[1] * target[1] + x[2] * target[2];

    // If nearly parallel
    if dot > 1.0 - 1e-12 {
        return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    }
    // If nearly anti-parallel, rotate 180° around z
    if dot < -1.0 + 1e-12 {
        return [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]];
    }

    // Rodrigues' rotation formula
    let cross = [
        x[1] * target[2] - x[2] * target[1],
        x[2] * target[0] - x[0] * target[2],
        x[0] * target[1] - x[1] * target[0],
    ];
    let s = (cross[0].powi(2) + cross[1].powi(2) + cross[2].powi(2)).sqrt();
    let c = dot;

    // Skew-symmetric cross-product matrix K
    // K = [  0   -v3   v2 ]
    //     [  v3   0   -v1 ]
    //     [ -v2   v1   0  ]
    let v = cross;
    let k = [[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]];

    // K²
    let k2 = mat_mul_3x3(k, k);

    // R = I + K + K² * (1 - c) / s²
    let factor = (1.0 - c) / (s * s);
    [
        [
            1.0 + k[0][0] + k2[0][0] * factor,
            k[0][1] + k2[0][1] * factor,
            k[0][2] + k2[0][2] * factor,
        ],
        [
            k[1][0] + k2[1][0] * factor,
            1.0 + k[1][1] + k2[1][1] * factor,
            k[1][2] + k2[1][2] * factor,
        ],
        [
            k[2][0] + k2[2][0] * factor,
            k[2][1] + k2[2][1] * factor,
            1.0 + k[2][2] + k2[2][2] * factor,
        ],
    ]
}

fn mat_mul_3x3(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    r
}

fn apply_rotation(rot: &[[f64; 3]; 3], p: [f64; 3]) -> [f64; 3] {
    [
        rot[0][0] * p[0] + rot[0][1] * p[1] + rot[0][2] * p[2],
        rot[1][0] * p[0] + rot[1][1] * p[1] + rot[1][2] * p[2],
        rot[2][0] * p[0] + rot[2][1] * p[1] + rot[2][2] * p[2],
    ]
}

#[cfg(test)]
mod tests {
    use super::super::sbu::{CoordinationGeometry, Sbu};
    use super::*;

    #[test]
    fn test_assemble_simple_cubic_framework() {
        let node = Sbu::metal_node(30, 0.0, CoordinationGeometry::Octahedral); // Zn at origin
        let linker = Sbu::linear_linker(&[6, 6], 1.4, "carboxylate");
        let topo = Topology::pcu();
        let cell = UnitCell::cubic(10.0);

        let structure = assemble_framework(&node, &linker, &topo, &cell);

        // 1 node with 1 atom + 3 edges with 2 atoms each = 7
        assert_eq!(structure.num_atoms(), 7);
        assert!(structure
            .atoms
            .iter()
            .all(|a| a.element == 30 || a.element == 6));
    }

    #[test]
    fn test_assemble_diamond_topology() {
        let node = Sbu::metal_node(14, 0.0, CoordinationGeometry::Tetrahedral); // Si
        let linker = Sbu::linear_linker(&[8], 1.0, "bridge"); // O bridge
        let topo = Topology::dia();
        let cell = UnitCell::cubic(5.43); // silicon lattice constant

        let structure = assemble_framework(&node, &linker, &topo, &cell);

        // 2 nodes + 4 edges = 2 × 1 + 4 × 1 = 6 atoms
        assert_eq!(structure.num_atoms(), 6);
    }

    #[test]
    fn test_supercell_atom_count() {
        let node = Sbu::metal_node(30, 0.0, CoordinationGeometry::Octahedral);
        let linker = Sbu::linear_linker(&[6, 6], 1.4, "carboxylate");
        let topo = Topology::pcu();
        let cell = UnitCell::cubic(10.0);

        let structure = assemble_framework(&node, &linker, &topo, &cell);
        let n_base = structure.num_atoms();

        let super_struct = structure.make_supercell(2, 2, 2);
        assert_eq!(super_struct.num_atoms(), n_base * 8);
    }

    #[test]
    fn test_topology_pcu() {
        let t = Topology::pcu();
        assert_eq!(t.nodes.len(), 1);
        assert_eq!(t.edges.len(), 3);
        assert_eq!(t.name, "pcu");
    }

    #[test]
    fn test_topology_dia() {
        let t = Topology::dia();
        assert_eq!(t.nodes.len(), 2);
        assert_eq!(t.edges.len(), 4);
    }

    #[test]
    fn test_rotation_identity() {
        let rot = rotation_from_x_to([1.0, 0.0, 0.0]);
        let p = [1.0, 2.0, 3.0];
        let rp = apply_rotation(&rot, p);
        for i in 0..3 {
            assert!((rp[i] - p[i]).abs() < 1e-10);
        }
    }

    #[test]
    fn test_rotation_to_y() {
        let rot = rotation_from_x_to([0.0, 1.0, 0.0]);
        let p = [1.0, 0.0, 0.0];
        let rp = apply_rotation(&rot, p);
        assert!((rp[0]).abs() < 1e-10, "x should be ~0, got {:.6}", rp[0]);
        assert!(
            (rp[1] - 1.0).abs() < 1e-10,
            "y should be ~1, got {:.6}",
            rp[1]
        );
        assert!((rp[2]).abs() < 1e-10, "z should be ~0, got {:.6}", rp[2]);
    }

    #[test]
    fn test_crystal_cartesian_coords() {
        let cell = UnitCell::cubic(10.0);
        let structure = CrystalStructure {
            cell,
            atoms: vec![
                CrystalAtom {
                    element: 6,
                    frac_coords: [0.5, 0.0, 0.0],
                },
                CrystalAtom {
                    element: 8,
                    frac_coords: [0.0, 0.5, 0.0],
                },
            ],
            labels: vec!["a".into(), "b".into()],
        };

        let coords = structure.cartesian_coords();
        assert!((coords[0][0] - 5.0).abs() < 1e-10);
        assert!((coords[1][1] - 5.0).abs() < 1e-10);
    }
}
