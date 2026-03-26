//! Benchmark: run EHT on transition metal reference systems and report
//! energies and gaps for regression tracking.
//!
//! Category: benchmark

use serde::Serialize;

#[derive(Debug, Serialize)]
struct MetalReferenceCase {
    name: String,
    homo_energy: f64,
    lumo_energy: f64,
    gap: f64,
    n_orbitals: usize,
    n_electrons: usize,
    support_level: String,
}

fn ferrocene_like() -> (Vec<u8>, Vec<[f64; 3]>) {
    let mut elements = Vec::new();
    let mut positions = Vec::new();

    elements.push(26); // Fe
    positions.push([0.0, 0.0, 0.0]);

    let ring_r_c = 1.65;
    let ring_r_h = 2.70;
    let z_top = 1.66;
    let z_bot = -1.66;

    for i in 0..5 {
        let a = (i as f64) * 72.0_f64.to_radians();
        elements.push(6);
        positions.push([ring_r_c * a.cos(), ring_r_c * a.sin(), z_top]);
        elements.push(1);
        positions.push([ring_r_h * a.cos(), ring_r_h * a.sin(), z_top + 0.35]);
    }

    for i in 0..5 {
        let a = (i as f64) * 72.0_f64.to_radians() + 36.0_f64.to_radians();
        elements.push(6);
        positions.push([ring_r_c * a.cos(), ring_r_c * a.sin(), z_bot]);
        elements.push(1);
        positions.push([ring_r_h * a.cos(), ring_r_h * a.sin(), z_bot - 0.35]);
    }

    (elements, positions)
}

fn cisplatin_like() -> (Vec<u8>, Vec<[f64; 3]>) {
    let elements = vec![
        78, // Pt
        17, 17, // Cl2
        7, 7, // 2x N
        1, 1, 1, 1, 1, 1, // 2x NH3 hydrogens
    ];
    let positions = vec![
        [0.0, 0.0, 0.0],
        [2.30, 0.0, 0.0],
        [-2.30, 0.0, 0.0],
        [0.0, 2.00, 0.0],
        [0.0, -2.00, 0.0],
        [0.90, 2.60, 0.0],
        [-0.90, 2.60, 0.0],
        [0.00, 2.00, 0.95],
        [0.90, -2.60, 0.0],
        [-0.90, -2.60, 0.0],
        [0.00, -2.00, -0.95],
    ];
    (elements, positions)
}

fn fecl6_octahedral() -> (Vec<u8>, Vec<[f64; 3]>) {
    let elements = vec![26, 17, 17, 17, 17, 17, 17];
    let positions = vec![
        [0.0, 0.0, 0.0],
        [2.30, 0.0, 0.0],
        [-2.30, 0.0, 0.0],
        [0.0, 2.30, 0.0],
        [0.0, -2.30, 0.0],
        [0.0, 0.0, 2.30],
        [0.0, 0.0, -2.30],
    ];
    (elements, positions)
}

fn pdcl4_square_planar() -> (Vec<u8>, Vec<[f64; 3]>) {
    let elements = vec![46, 17, 17, 17, 17];
    let positions = vec![
        [0.0, 0.0, 0.0],
        [2.30, 0.0, 0.0],
        [-2.30, 0.0, 0.0],
        [0.0, 2.30, 0.0],
        [0.0, -2.30, 0.0],
    ];
    (elements, positions)
}

fn main() {
    let systems = vec![
        ("ferrocene_like", ferrocene_like()),
        ("cisplatin_like", cisplatin_like()),
        ("fecl6_octahedral", fecl6_octahedral()),
        ("pdcl4_square_planar", pdcl4_square_planar()),
    ];

    let mut out = Vec::new();
    for (name, (elements, positions)) in systems {
        let result = sci_form::eht::solve_eht(&elements, &positions, None)
            .unwrap_or_else(|e| panic!("{} failed: {}", name, e));
        out.push(MetalReferenceCase {
            name: name.to_string(),
            homo_energy: result.homo_energy,
            lumo_energy: result.lumo_energy,
            gap: result.gap,
            n_orbitals: result.energies.len(),
            n_electrons: result.n_electrons,
            support_level: format!("{:?}", result.support.level).to_lowercase(),
        });
    }

    println!("{}", serde_json::to_string_pretty(&out).unwrap());
}
