use serde::Deserialize;
use std::fs;
use std::time::Instant;

#[derive(Deserialize, Debug)]
struct OracleAtom {
    element: u8,
    x: f32,
    y: f32,
    z: f32,
    formal_charge: i8,
    hybridization: String,
}

#[derive(Deserialize, Debug)]
struct OracleBond {
    start: usize,
    end: usize,
    order: String,
}

#[derive(Deserialize, Debug)]
struct OracleMolecule {
    smiles: String,
    atoms: Vec<OracleAtom>,
    bonds: Vec<OracleBond>,
}

#[derive(Deserialize, Debug)]
struct BenchmarkData {
    rdkit_time_ms: f64,
    count: usize,
    molecules: Vec<OracleMolecule>,
}

fn main() {
    let data = fs::read_to_string("tests/fixtures/benchmark_data.json")
        .expect("Error reading benchmark_data.json");
    let bench_data: BenchmarkData = serde_json::from_str(&data).expect("Error parsing JSON");

    let mut total_rmsd = 0.0;

    /*
     * Perform exactly the embedding sequence without optimization
     * Because we are measuring the RAW ETKDG coordinates generated!
     */
    let pb = Instant::now();
    for mol in &bench_data.molecules {
        let n = mol.atoms.len();
        let mut our_mol = sci_form::graph::Molecule::new(&mol.smiles);
        let mut node_indices = Vec::new();

        for atom in &mol.atoms {
            let hybridization = match atom.hybridization.as_str() {
                "SP" => sci_form::graph::Hybridization::SP,
                "SP2" => sci_form::graph::Hybridization::SP2,
                "SP3" => sci_form::graph::Hybridization::SP3,
                "SP3D" => sci_form::graph::Hybridization::SP3D,
                "SP3D2" => sci_form::graph::Hybridization::SP3D2,
                _ => sci_form::graph::Hybridization::Unknown,
            };

            let new_atom = sci_form::graph::Atom {
                element: atom.element,
                position: nalgebra::Vector3::zeros(),
                charge: 0.0,
                formal_charge: atom.formal_charge,
                hybridization,
                chiral_tag: sci_form::graph::ChiralType::Unspecified,
                explicit_h: 0,
            };
            let n_idx = our_mol.add_atom(new_atom);
            node_indices.push(n_idx);
        }

        for bond in &mol.bonds {
            let order = match bond.order.as_str() {
                "DOUBLE" => sci_form::graph::BondOrder::Double,
                "TRIPLE" => sci_form::graph::BondOrder::Triple,
                "AROMATIC" => sci_form::graph::BondOrder::Aromatic,
                _ => sci_form::graph::BondOrder::Single,
            };

            our_mol.add_bond(
                node_indices[bond.start],
                node_indices[bond.end],
                sci_form::graph::Bond {
                    order,
                    stereo: sci_form::graph::BondStereo::None,
                },
            );
        }

        let bounds = sci_form::distgeom::calculate_bounds_matrix(&our_mol);
        let smoothed = sci_form::distgeom::smooth_bounds_matrix(bounds);

        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        use rand::SeedableRng; // required for seed_from_u64

        let dists = sci_form::distgeom::pick_random_distances(&mut rng, &smoothed);
        let metric = sci_form::distgeom::compute_metric_matrix(&dists);
        let coords3d = sci_form::distgeom::generate_3d_coordinates(&mut rng, &metric);

        let mut diff_raw_sq_sum = 0.0;
        let mut pairs = 0;
        for i in 0..n {
            for j in (i + 1)..n {
                let dx_r = mol.atoms[i].x - mol.atoms[j].x;
                let dy_r = mol.atoms[i].y - mol.atoms[j].y;
                let dz_r = mol.atoms[i].z - mol.atoms[j].z;
                let rdkit_dist = (dx_r * dx_r + dy_r * dy_r + dz_r * dz_r).sqrt();

                let dx_raw = coords3d[(i, 0)] - coords3d[(j, 0)];
                let dy_raw = coords3d[(i, 1)] - coords3d[(j, 1)];
                let dz_raw = coords3d[(i, 2)] - coords3d[(j, 2)];
                let raw_dist = (dx_raw * dx_raw + dy_raw * dy_raw + dz_raw * dz_raw).sqrt();

                let diff_raw = (raw_dist - rdkit_dist).abs();
                diff_raw_sq_sum += diff_raw * diff_raw;
                pairs += 1;
            }
        }

        if pairs > 0 {
            total_rmsd += (diff_raw_sq_sum / pairs as f32).sqrt();
        }
    }
    let rust_time = pb.elapsed();

    let avg_rmsd = total_rmsd / bench_data.count as f32;
    let rust_time_ms = rust_time.as_secs_f64() * 1000.0;

    println!("=== PERFORMANCE BENCHMARK: RDKit vs Sci-Form ===");
    println!("Molecules Processed: {}", bench_data.count);
    println!("Raw Distance RMSD (vs RDKit): {:.3} Å", avg_rmsd);
    println!(
        "RDKit Total Time (Python): {:.2} ms",
        bench_data.rdkit_time_ms
    );
    println!("Sci-Form Total Time (Rust): {:.2} ms", rust_time_ms);
    println!(
        "Sci-Form is {:.2}x faster!",
        bench_data.rdkit_time_ms / rust_time_ms
    );
}
