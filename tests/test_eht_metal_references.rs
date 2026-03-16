use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct MetalReference {
    name: String,
    homo_energy: f64,
    lumo_energy: f64,
    gap: f64,
    n_orbitals: usize,
    n_electrons: usize,
    support_level: String,
    tolerance: f64,
    max_variation_percent: f64,
}

fn relative_error_percent(predicted: f64, reference: f64) -> Option<f64> {
    if reference.abs() < 1e-6 {
        return None;
    }
    Some(((predicted - reference).abs() / reference.abs()) * 100.0)
}

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

fn system_by_name(name: &str) -> (Vec<u8>, Vec<[f64; 3]>) {
    match name {
        "ferrocene_like" => ferrocene_like(),
        "cisplatin_like" => cisplatin_like(),
        "fecl6_octahedral" => fecl6_octahedral(),
        "pdcl4_square_planar" => pdcl4_square_planar(),
        _ => panic!("Unknown reference system: {}", name),
    }
}

#[test]
fn test_eht_metal_reference_regression() {
    let text = std::fs::read_to_string("tests/fixtures/eht_metal_reference.json")
        .expect("Failed to read metal reference fixture");
    let refs: Vec<MetalReference> =
        serde_json::from_str(&text).expect("Failed to parse metal reference fixture JSON");

    for reference in refs {
        let (elements, positions) = system_by_name(&reference.name);
        let result = sci_form::eht::solve_eht(&elements, &positions, None)
            .unwrap_or_else(|e| panic!("EHT solve failed for {}: {}", reference.name, e));

        assert_eq!(
            result.energies.len(),
            reference.n_orbitals,
            "Unexpected MO count for {}",
            reference.name
        );
        assert_eq!(
            result.n_electrons, reference.n_electrons,
            "Unexpected valence electron count for {}",
            reference.name
        );
        assert!(
            result.energies.windows(2).all(|w| w[0] <= w[1] + 1e-12),
            "MO energy ordering is not non-decreasing for {}",
            reference.name
        );
        assert!(
            result.homo_energy <= result.lumo_energy + 1e-12,
            "Frontier orbital ordering invalid for {}",
            reference.name
        );

        assert!(
            (result.homo_energy - reference.homo_energy).abs() <= reference.tolerance,
            "HOMO mismatch for {}: got {}, expected {}",
            reference.name,
            result.homo_energy,
            reference.homo_energy
        );
        let homo_variation = relative_error_percent(result.homo_energy, reference.homo_energy)
            .expect("HOMO reference should be non-zero for variation checks");
        assert!(
            homo_variation <= reference.max_variation_percent,
            "HOMO variation above threshold for {}: {:.6}% > {:.6}%",
            reference.name,
            homo_variation,
            reference.max_variation_percent
        );

        assert!(
            (result.lumo_energy - reference.lumo_energy).abs() <= reference.tolerance,
            "LUMO mismatch for {}: got {}, expected {}",
            reference.name,
            result.lumo_energy,
            reference.lumo_energy
        );
        let lumo_variation = relative_error_percent(result.lumo_energy, reference.lumo_energy)
            .expect("LUMO reference should be non-zero for variation checks");
        assert!(
            lumo_variation <= reference.max_variation_percent,
            "LUMO variation above threshold for {}: {:.6}% > {:.6}%",
            reference.name,
            lumo_variation,
            reference.max_variation_percent
        );

        assert!(
            (result.gap - reference.gap).abs() <= reference.tolerance,
            "Gap mismatch for {}: got {}, expected {}",
            reference.name,
            result.gap,
            reference.gap
        );
        if let Some(gap_variation) = relative_error_percent(result.gap, reference.gap) {
            assert!(
                gap_variation <= reference.max_variation_percent,
                "Gap variation above threshold for {}: {:.6}% > {:.6}%",
                reference.name,
                gap_variation,
                reference.max_variation_percent
            );
        }

        assert_eq!(
            format!("{:?}", result.support.level).to_lowercase(),
            reference.support_level,
            "Support level mismatch for {}",
            reference.name
        );
    }
}
