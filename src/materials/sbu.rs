//! Secondary Building Unit (SBU) representation for framework assembly.
//!
//! An SBU is a molecular fragment with designated connection points. Framework
//! structures (MOFs, COFs, zeolites) are built by connecting SBUs at their
//! connection sites according to geometric rules.

use serde::{Deserialize, Serialize};

/// A connection point on an SBU where another SBU can attach.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConnectionPoint {
    /// Position in Cartesian coordinates (Å) local to the SBU.
    pub position: [f64; 3],
    /// Outward-pointing direction vector (unit vector).
    pub direction: [f64; 3],
    /// Coordination type label (e.g., "carboxylate", "amine", "metal").
    pub kind: String,
}

/// Coordination geometry of a metal center or node.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum CoordinationGeometry {
    Linear,
    Trigonal,
    Tetrahedral,
    SquarePlanar,
    Octahedral,
}

impl CoordinationGeometry {
    /// Number of connection points for this geometry.
    pub fn num_connections(self) -> usize {
        match self {
            CoordinationGeometry::Linear => 2,
            CoordinationGeometry::Trigonal => 3,
            CoordinationGeometry::Tetrahedral => 4,
            CoordinationGeometry::SquarePlanar => 4,
            CoordinationGeometry::Octahedral => 6,
        }
    }

    /// Generate ideal direction vectors for this coordination geometry.
    pub fn ideal_directions(self) -> Vec<[f64; 3]> {
        match self {
            CoordinationGeometry::Linear => vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
            CoordinationGeometry::Trigonal => {
                let s3 = (3.0f64).sqrt() / 2.0;
                vec![[1.0, 0.0, 0.0], [-0.5, s3, 0.0], [-0.5, -s3, 0.0]]
            }
            CoordinationGeometry::Tetrahedral => {
                let c = 1.0 / (3.0f64).sqrt();
                vec![
                    [c, c, c],
                    [c, -c, -c],
                    [-c, c, -c],
                    [-c, -c, c],
                ]
            }
            CoordinationGeometry::SquarePlanar => {
                vec![
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [0.0, -1.0, 0.0],
                ]
            }
            CoordinationGeometry::Octahedral => {
                vec![
                    [1.0, 0.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, -1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [0.0, 0.0, -1.0],
                ]
            }
        }
    }
}

/// A Secondary Building Unit (SBU) — a molecular fragment with connection points.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Sbu {
    /// Human-readable label (e.g., "Zn4O", "BDC_linker").
    pub label: String,
    /// Atomic numbers.
    pub elements: Vec<u8>,
    /// Cartesian coordinates (Å), local frame.
    pub positions: Vec<[f64; 3]>,
    /// Connection points where other SBUs can attach.
    pub connections: Vec<ConnectionPoint>,
    /// Coordination geometry of the central metal/node (if applicable).
    pub geometry: Option<CoordinationGeometry>,
}

impl Sbu {
    /// Create a new SBU.
    pub fn new(label: &str) -> Self {
        Self {
            label: label.to_string(),
            elements: Vec::new(),
            positions: Vec::new(),
            connections: Vec::new(),
            geometry: None,
        }
    }

    /// Add an atom to the SBU.
    pub fn add_atom(&mut self, element: u8, position: [f64; 3]) {
        self.elements.push(element);
        self.positions.push(position);
    }

    /// Add a connection point.
    pub fn add_connection(&mut self, position: [f64; 3], direction: [f64; 3], kind: &str) {
        let mag = (direction[0].powi(2) + direction[1].powi(2) + direction[2].powi(2)).sqrt();
        let normalized = if mag > 1e-12 {
            [direction[0] / mag, direction[1] / mag, direction[2] / mag]
        } else {
            [0.0, 0.0, 1.0]
        };
        self.connections.push(ConnectionPoint {
            position,
            direction: normalized,
            kind: kind.to_string(),
        });
    }

    /// Number of atoms in this SBU.
    pub fn num_atoms(&self) -> usize {
        self.elements.len()
    }

    /// Number of connection points.
    pub fn num_connections(&self) -> usize {
        self.connections.len()
    }

    /// Centroid of all atoms.
    pub fn centroid(&self) -> [f64; 3] {
        let n = self.positions.len() as f64;
        if n < 1.0 {
            return [0.0; 3];
        }
        let mut c = [0.0; 3];
        for p in &self.positions {
            c[0] += p[0];
            c[1] += p[1];
            c[2] += p[2];
        }
        [c[0] / n, c[1] / n, c[2] / n]
    }

    /// Create a metal-node SBU with ideal coordination geometry at the origin.
    ///
    /// `element`: atomic number of the metal center
    /// `bond_length`: metal–connection distance (Å)
    /// `geometry`: desired coordination geometry
    pub fn metal_node(
        element: u8,
        bond_length: f64,
        geometry: CoordinationGeometry,
    ) -> Self {
        let mut sbu = Self::new("metal_node");
        sbu.add_atom(element, [0.0, 0.0, 0.0]);
        sbu.geometry = Some(geometry);

        for dir in geometry.ideal_directions() {
            let pos = [
                dir[0] * bond_length,
                dir[1] * bond_length,
                dir[2] * bond_length,
            ];
            sbu.add_connection(pos, dir, "metal");
        }

        sbu
    }

    /// Create a linear linker SBU between two connection points.
    ///
    /// Places atoms along the x-axis with connection points at each end.
    pub fn linear_linker(
        elements: &[u8],
        spacing: f64,
        kind: &str,
    ) -> Self {
        let mut sbu = Self::new("linear_linker");
        let n = elements.len();
        let total_len = (n - 1) as f64 * spacing;

        for (i, &e) in elements.iter().enumerate() {
            let x = i as f64 * spacing - total_len / 2.0;
            sbu.add_atom(e, [x, 0.0, 0.0]);
        }

        // Connection points at both ends, pointing outward
        let half = total_len / 2.0 + spacing / 2.0;
        sbu.add_connection([-half, 0.0, 0.0], [-1.0, 0.0, 0.0], kind);
        sbu.add_connection([half, 0.0, 0.0], [1.0, 0.0, 0.0], kind);

        sbu
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_metal_node_tetrahedral() {
        let sbu = Sbu::metal_node(30, 2.0, CoordinationGeometry::Tetrahedral); // Zn
        assert_eq!(sbu.num_atoms(), 1);
        assert_eq!(sbu.num_connections(), 4);
        assert_eq!(sbu.elements[0], 30);
    }

    #[test]
    fn test_metal_node_octahedral() {
        let sbu = Sbu::metal_node(26, 2.1, CoordinationGeometry::Octahedral); // Fe
        assert_eq!(sbu.num_connections(), 6);
    }

    #[test]
    fn test_linear_linker() {
        let linker = Sbu::linear_linker(&[6, 6, 6, 6, 6, 6], 1.4, "carboxylate"); // benzene-like
        assert_eq!(linker.num_atoms(), 6);
        assert_eq!(linker.num_connections(), 2);
    }

    #[test]
    fn test_connection_direction_normalized() {
        let mut sbu = Sbu::new("test");
        sbu.add_connection([0.0, 0.0, 0.0], [3.0, 4.0, 0.0], "test");
        let dir = sbu.connections[0].direction;
        let mag = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
        assert!((mag - 1.0).abs() < 1e-10, "Direction should be unit vector, got mag={mag:.6}");
    }

    #[test]
    fn test_sbu_centroid() {
        let mut sbu = Sbu::new("test");
        sbu.add_atom(6, [0.0, 0.0, 0.0]);
        sbu.add_atom(6, [2.0, 0.0, 0.0]);
        let c = sbu.centroid();
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[1]).abs() < 1e-10);
    }

    #[test]
    fn test_coordination_num_connections() {
        assert_eq!(CoordinationGeometry::Linear.num_connections(), 2);
        assert_eq!(CoordinationGeometry::Tetrahedral.num_connections(), 4);
        assert_eq!(CoordinationGeometry::Octahedral.num_connections(), 6);
    }
}
