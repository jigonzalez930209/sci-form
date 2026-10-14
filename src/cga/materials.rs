//! CGA-based Materials Assembly — E1.3
//!
//! Encodes unit cell lattice vectors as CGA translation motors,
//! places SBUs via composed motors, and assembles framework structures.

use super::motor::Motor;

/// A crystallographic frame represented as CGA translation motors.
#[derive(Debug, Clone)]
pub struct CgaFrame {
    /// Motor for lattice vector a.
    pub motor_a: Motor,
    /// Motor for lattice vector b.
    pub motor_b: Motor,
    /// Motor for lattice vector c.
    pub motor_c: Motor,
    /// Lattice vectors (Å).
    pub lattice: [[f64; 3]; 3],
}

impl CgaFrame {
    /// Build a CGA frame from three lattice vectors.
    pub fn from_lattice(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Self {
        Self {
            motor_a: Motor::translator(a),
            motor_b: Motor::translator(b),
            motor_c: Motor::translator(c),
            lattice: [a, b, c],
        }
    }

    /// Build from cell parameters (a, b, c, α, β, γ in degrees).
    pub fn from_parameters(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> Self {
        let alpha_r = alpha.to_radians();
        let beta_r = beta.to_radians();
        let gamma_r = gamma.to_radians();

        let va = [a, 0.0, 0.0];
        let vb = [b * gamma_r.cos(), b * gamma_r.sin(), 0.0];
        let cx = c * beta_r.cos();
        let cy = c * (alpha_r.cos() - beta_r.cos() * gamma_r.cos()) / gamma_r.sin();
        let cz = (c * c - cx * cx - cy * cy).max(0.0).sqrt();
        let vc = [cx, cy, cz];

        Self::from_lattice(va, vb, vc)
    }

    /// Convert fractional coordinates to Cartesian using motor composition.
    pub fn frac_to_cart(&self, frac: [f64; 3]) -> [f64; 3] {
        [
            frac[0] * self.lattice[0][0]
                + frac[1] * self.lattice[1][0]
                + frac[2] * self.lattice[2][0],
            frac[0] * self.lattice[0][1]
                + frac[1] * self.lattice[1][1]
                + frac[2] * self.lattice[2][1],
            frac[0] * self.lattice[0][2]
                + frac[1] * self.lattice[1][2]
                + frac[2] * self.lattice[2][2],
        ]
    }

    /// Convert Cartesian to fractional coordinates.
    pub fn cart_to_frac(&self, cart: [f64; 3]) -> [f64; 3] {
        // Solve lattice * frac = cart via Cramer's rule
        let l = &self.lattice;
        let det = l[0][0] * (l[1][1] * l[2][2] - l[1][2] * l[2][1])
            - l[1][0] * (l[0][1] * l[2][2] - l[0][2] * l[2][1])
            + l[2][0] * (l[0][1] * l[1][2] - l[0][2] * l[1][1]);

        if det.abs() < 1e-15 {
            return [0.0, 0.0, 0.0];
        }

        let inv_det = 1.0 / det;
        let fa = inv_det
            * (cart[0] * (l[1][1] * l[2][2] - l[1][2] * l[2][1])
                - l[1][0] * (cart[1] * l[2][2] - cart[2] * l[2][1])
                + l[2][0] * (cart[1] * l[1][2] - cart[2] * l[1][1]));
        let fb = inv_det
            * (l[0][0] * (cart[1] * l[2][2] - cart[2] * l[2][1])
                - cart[0] * (l[0][1] * l[2][2] - l[0][2] * l[2][1])
                + l[2][0] * (l[0][1] * cart[2] - l[0][2] * cart[1]));
        let fc = inv_det
            * (l[0][0] * (l[1][1] * cart[2] - l[1][2] * cart[1])
                - l[1][0] * (l[0][1] * cart[2] - l[0][2] * cart[1])
                + cart[0] * (l[0][1] * l[1][2] - l[0][2] * l[1][1]));

        [fa, fb, fc]
    }

    /// Build a composite motor for a supercell offset (na, nb, nc).
    pub fn supercell_motor(&self, na: i32, nb: i32, nc: i32) -> Motor {
        let t = [
            na as f64 * self.lattice[0][0]
                + nb as f64 * self.lattice[1][0]
                + nc as f64 * self.lattice[2][0],
            na as f64 * self.lattice[0][1]
                + nb as f64 * self.lattice[1][1]
                + nc as f64 * self.lattice[2][1],
            na as f64 * self.lattice[0][2]
                + nb as f64 * self.lattice[1][2]
                + nc as f64 * self.lattice[2][2],
        ];
        Motor::translator(t)
    }
}

/// A placed atom in the CGA-assembled crystal.
#[derive(Debug, Clone)]
pub struct CgaCrystalAtom {
    pub element: u8,
    pub cart_coords: [f64; 3],
    pub frac_coords: [f64; 3],
    pub label: String,
}

/// A CGA-assembled crystal structure.
#[derive(Debug, Clone)]
pub struct CgaCrystalStructure {
    pub frame: CgaFrame,
    pub atoms: Vec<CgaCrystalAtom>,
}

impl CgaCrystalStructure {
    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    pub fn elements(&self) -> Vec<u8> {
        self.atoms.iter().map(|a| a.element).collect()
    }

    pub fn cartesian_coords(&self) -> Vec<[f64; 3]> {
        self.atoms.iter().map(|a| a.cart_coords).collect()
    }
}

/// Place an SBU at a topology node using a CGA Motor.
///
/// `sbu_atoms`: list of (element, local_coord) pairs for the SBU.
/// `motor`: the positioning motor M = T_i R_i
pub fn place_sbu_cga(
    sbu_atoms: &[(u8, [f64; 3])],
    motor: &Motor,
    label: &str,
) -> Vec<(u8, [f64; 3], String)> {
    sbu_atoms
        .iter()
        .map(|&(el, coord)| {
            let placed = motor.transform_point(coord);
            (el, placed, label.to_string())
        })
        .collect()
}

/// Assemble a framework using CGA motors for the PCU topology.
///
/// `node_element`: atomic number for node atoms
/// `linker_element`: atomic number for linker atoms
/// `lattice_a`: cubic lattice parameter (Å)
/// `supercell`: number of unit cell repeats per direction
pub fn assemble_framework_cga(
    node_element: u8,
    linker_element: u8,
    lattice_a: f64,
    supercell: usize,
) -> CgaCrystalStructure {
    let frame = CgaFrame::from_parameters(lattice_a, lattice_a, lattice_a, 90.0, 90.0, 90.0);

    let mut atoms = Vec::new();

    // PCU topology: single node at origin, 6 edges along ±a, ±b, ±c
    let node_frac = [0.0, 0.0, 0.0];

    // Linker positions at midpoints of edges
    let linker_fracs = [[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]];

    for na in 0..(supercell as i32) {
        for nb in 0..(supercell as i32) {
            for nc in 0..(supercell as i32) {
                let sc_motor = frame.supercell_motor(na, nb, nc);

                // Place node
                let node_cart = frame.frac_to_cart(node_frac);
                let placed_node = sc_motor.transform_point(node_cart);
                let node_frac_out = [
                    (node_frac[0] + na as f64) / supercell as f64,
                    (node_frac[1] + nb as f64) / supercell as f64,
                    (node_frac[2] + nc as f64) / supercell as f64,
                ];
                atoms.push(CgaCrystalAtom {
                    element: node_element,
                    cart_coords: placed_node,
                    frac_coords: node_frac_out,
                    label: format!("node_{}_{}_{}", na, nb, nc),
                });

                // Place linkers
                for (li, &lf) in linker_fracs.iter().enumerate() {
                    let linker_cart = frame.frac_to_cart(lf);
                    let placed_linker = sc_motor.transform_point(linker_cart);
                    let l_frac_out = [
                        (lf[0] + na as f64) / supercell as f64,
                        (lf[1] + nb as f64) / supercell as f64,
                        (lf[2] + nc as f64) / supercell as f64,
                    ];
                    atoms.push(CgaCrystalAtom {
                        element: linker_element,
                        cart_coords: placed_linker,
                        frac_coords: l_frac_out,
                        label: format!("linker{}_{}_{}_{}", li, na, nb, nc),
                    });
                }
            }
        }
    }

    CgaCrystalStructure { frame, atoms }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cubic_frame_frac_to_cart() {
        let frame = CgaFrame::from_parameters(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
        let cart = frame.frac_to_cart([0.5, 0.5, 0.5]);
        assert!((cart[0] - 5.0).abs() < 1e-10);
        assert!((cart[1] - 5.0).abs() < 1e-10);
        assert!((cart[2] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_frac_cart_round_trip() {
        let frame = CgaFrame::from_parameters(10.0, 12.0, 8.0, 90.0, 90.0, 90.0);
        let frac = [0.3, 0.7, 0.1];
        let cart = frame.frac_to_cart(frac);
        let frac2 = frame.cart_to_frac(cart);
        for i in 0..3 {
            assert!(
                (frac[i] - frac2[i]).abs() < 1e-10,
                "Frac round-trip error at {}: {} → {}",
                i,
                frac[i],
                frac2[i]
            );
        }
    }

    #[test]
    fn test_pcu_assembly_atom_count() {
        // 1×1×1: 1 node + 3 linkers = 4 atoms
        let crystal = assemble_framework_cga(30, 6, 10.0, 1);
        assert_eq!(crystal.num_atoms(), 4);

        // 2×2×2: 8 nodes + 24 linkers = 32 atoms
        let crystal = assemble_framework_cga(30, 6, 10.0, 2);
        assert_eq!(crystal.num_atoms(), 32);
    }

    #[test]
    fn test_pcu_assembly_distances() {
        let crystal = assemble_framework_cga(30, 6, 10.0, 1);
        // Node at origin, linker at (5,0,0)
        let node = &crystal.atoms[0];
        let linker = &crystal.atoms[1];
        let dx = linker.cart_coords[0] - node.cart_coords[0];
        let dy = linker.cart_coords[1] - node.cart_coords[1];
        let dz = linker.cart_coords[2] - node.cart_coords[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        assert!(
            (dist - 5.0).abs() < 0.001,
            "Node-linker distance = {}, expected 5.0",
            dist
        );
    }

    #[test]
    fn test_supercell_motor_translation() {
        let frame = CgaFrame::from_parameters(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
        let motor = frame.supercell_motor(1, 2, 3);
        let result = motor.transform_point([0.0, 0.0, 0.0]);
        assert!((result[0] - 10.0).abs() < 1e-10);
        assert!((result[1] - 20.0).abs() < 1e-10);
        assert!((result[2] - 30.0).abs() < 1e-10);
    }

    #[test]
    fn test_place_sbu_cga() {
        let sbu = vec![(30u8, [0.0, 0.0, 0.0]), (8u8, [1.0, 0.0, 0.0])];
        let motor = Motor::translator([5.0, 5.0, 5.0]);
        let placed = place_sbu_cga(&sbu, &motor, "test_sbu");
        assert_eq!(placed.len(), 2);
        assert!((placed[0].1[0] - 5.0).abs() < 1e-10);
        assert!((placed[1].1[0] - 6.0).abs() < 1e-10);
    }
}
