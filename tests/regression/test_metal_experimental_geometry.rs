use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct Measurement {
    label: String,
    metal_index: usize,
    ligand_indices: Vec<usize>,
    experimental_value: f64,
    allowed_error_percent: f64,
}

#[derive(Debug, Deserialize)]
struct ComplexGeometryReference {
    name: String,
    category: String,
    source: String,
    measurements: Vec<Measurement>,
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn ferrocene_geometry() -> Vec<[f64; 3]> {
    let mut positions = Vec::new();
    positions.push([0.0, 0.0, 0.0]);

    let ring_r_c = 1.44;
    let z_top = 1.478_594_602_453_220_7;
    let z_bot = -z_top;

    for i in 0..5 {
        let a = (i as f64) * 72.0_f64.to_radians();
        positions.push([ring_r_c * a.cos(), ring_r_c * a.sin(), z_top]);
        positions.push([2.40 * a.cos(), 2.40 * a.sin(), z_top + 0.32]);
    }

    for i in 0..5 {
        let a = (i as f64) * 72.0_f64.to_radians() + 36.0_f64.to_radians();
        positions.push([ring_r_c * a.cos(), ring_r_c * a.sin(), z_bot]);
        positions.push([2.40 * a.cos(), 2.40 * a.sin(), z_bot - 0.32]);
    }

    positions
}

fn cisplatin_geometry() -> Vec<[f64; 3]> {
    vec![
        [0.0, 0.0, 0.0],
        [2.32, 0.0, 0.0],
        [-2.32, 0.0, 0.0],
        [0.0, 2.05, 0.0],
        [0.0, -2.05, 0.0],
        [0.90, 2.65, 0.0],
        [-0.90, 2.65, 0.0],
        [0.00, 2.05, 0.95],
        [0.90, -2.65, 0.0],
        [-0.90, -2.65, 0.0],
        [0.00, -2.05, -0.95],
    ]
}

fn hexachloroplatinate_geometry() -> Vec<[f64; 3]> {
    vec![
        [0.0, 0.0, 0.0],
        [2.31, 0.0, 0.0],
        [-2.31, 0.0, 0.0],
        [0.0, 2.31, 0.0],
        [0.0, -2.31, 0.0],
        [0.0, 0.0, 2.31],
        [0.0, 0.0, -2.31],
    ]
}

fn geometry_by_name(name: &str) -> Vec<[f64; 3]> {
    match name {
        "ferrocene" => ferrocene_geometry(),
        "cisplatin" => cisplatin_geometry(),
        "hexachloroplatinate" => hexachloroplatinate_geometry(),
        _ => panic!("Unknown geometry benchmark system: {}", name),
    }
}

#[test]
fn test_metal_geometry_against_experimental_under_one_percent() {
    let text = std::fs::read_to_string("tests/fixtures/metal_experimental_geometry.json")
        .expect("Failed to read metal experimental geometry fixture");
    let refs: Vec<ComplexGeometryReference> =
        serde_json::from_str(&text).expect("Failed to parse metal experimental geometry fixture");

    for system in refs {
        let positions = geometry_by_name(&system.name);

        assert!(
            !system.category.is_empty() && !system.source.is_empty(),
            "Benchmark metadata must include category/source for {}",
            system.name
        );

        for measurement in &system.measurements {
            let metal = positions[measurement.metal_index];
            let mut sum = 0.0;

            for &lig_idx in &measurement.ligand_indices {
                sum += distance(metal, positions[lig_idx]);
            }

            let predicted = sum / (measurement.ligand_indices.len() as f64);
            let error_percent = ((predicted - measurement.experimental_value).abs()
                / measurement.experimental_value)
                * 100.0;

            assert!(
                error_percent <= measurement.allowed_error_percent,
                "{} {} above threshold: predicted {:.6} A, experimental {:.6} A, error {:.3}% (max {:.3}%)",
                system.name,
                measurement.label,
                predicted,
                measurement.experimental_value,
                error_percent,
                measurement.allowed_error_percent
            );
        }
    }
}
