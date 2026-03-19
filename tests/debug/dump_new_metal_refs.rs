//! Temporary helper: dump new EHT reference values for metal systems.
//! Run with: cargo test --test dump_new_metal_refs -- --nocapture

fn ferrocene_like() -> (Vec<u8>, Vec<[f64; 3]>) {
    let mut elements = Vec::new();
    let mut positions = Vec::new();
    elements.push(26);
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
    let elements = vec![78, 17, 17, 7, 7, 1, 1, 1, 1, 1, 1];
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

#[test]
fn dump_metal_refs() {
    type System = (&'static str, Vec<u8>, Vec<[f64; 3]>);
    let systems: Vec<System> = vec![
        {
            let (e, p) = ferrocene_like();
            ("ferrocene_like", e, p)
        },
        {
            let (e, p) = cisplatin_like();
            ("cisplatin_like", e, p)
        },
        {
            let (e, p) = fecl6_octahedral();
            ("fecl6_octahedral", e, p)
        },
        {
            let (e, p) = pdcl4_square_planar();
            ("pdcl4_square_planar", e, p)
        },
    ];

    println!("\n[");
    for (i, (name, elements, positions)) in systems.iter().enumerate() {
        let result = sci_form::eht::solve_eht(elements, positions, None).unwrap();
        println!("  {{");
        println!("    \"name\": \"{}\",", name);
        println!("    \"homo_energy\": {},", result.homo_energy);
        println!("    \"lumo_energy\": {},", result.lumo_energy);
        println!("    \"gap\": {},", result.gap);
        println!("    \"n_orbitals\": {},", result.energies.len());
        println!("    \"n_electrons\": {},", result.n_electrons);
        println!("    \"support_level\": \"experimental\",");
        println!("    \"tolerance\": 0.0001,");
        println!("    \"max_variation_percent\": 0.1");
        print!("  }}");
        if i < systems.len() - 1 {
            println!(",");
        } else {
            println!();
        }

        // Verify sanity
        println!("  // HOMO-LUMO gap: {:.6} eV (non-zero = good)", result.gap);
        println!(
            "  // n_electrons: {} (even = good: {})",
            result.n_electrons,
            result.n_electrons.is_multiple_of(2)
        );
    }
    println!("]");
}
