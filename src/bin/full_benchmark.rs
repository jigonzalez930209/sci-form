use sci_form::forcefield::minimizer::minimize_energy_lbfgs;
use sci_form::forcefield::FFParams;
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

    let mut raw_total_rmsd = 0.0;
    let mut min_total_rmsd = 0.0;

    let pb = Instant::now();
    for mol in &bench_data.molecules {
        let n = mol.atoms.len();

        // NATIVE PARSING FROM RAW SMILES
        let mut our_mol =
            sci_form::graph::Molecule::from_smiles(&mol.smiles).expect("SMILES parsing failed");

        let bounds = sci_form::distgeom::calculate_bounds_matrix(&our_mol);
        let smoothed = sci_form::distgeom::smooth_bounds_matrix(bounds);

        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        use rand::SeedableRng; // required for seed_from_u64

        let dists = sci_form::distgeom::pick_random_distances(&mut rng, &smoothed);
        let metric = sci_form::distgeom::compute_metric_matrix(&dists);
        let mut coords3d_raw = sci_form::distgeom::generate_3d_coordinates(&mut rng, &metric);
        let mut coords3d_min = coords3d_raw.clone();

        // 1. Calculate Raw RMSD
        let mut raw_diff_sq_sum = 0.0;
        let mut pairs = 0;
        for i in 0..n {
            for j in (i + 1)..n {
                let dx_r = mol.atoms[i].x - mol.atoms[j].x;
                let dy_r = mol.atoms[i].y - mol.atoms[j].y;
                let dz_r = mol.atoms[i].z - mol.atoms[j].z;
                let rdkit_dist = (dx_r * dx_r + dy_r * dy_r + dz_r * dz_r).sqrt();

                let dx_our = coords3d_raw[(i, 0)] - coords3d_raw[(j, 0)];
                let dy_our = coords3d_raw[(i, 1)] - coords3d_raw[(j, 1)];
                let dz_our = coords3d_raw[(i, 2)] - coords3d_raw[(j, 2)];
                let our_dist = (dx_our * dx_our + dy_our * dy_our + dz_our * dz_our).sqrt();

                let diff = (our_dist - rdkit_dist).abs();
                raw_diff_sq_sum += diff * diff;
                pairs += 1;
            }
        }
        if pairs > 0 {
            raw_total_rmsd += (raw_diff_sq_sum / pairs as f32).sqrt();
        }

        // 2. FULL L-BFGS MINIMIZATION
        let params = FFParams {
            kb: 200.0,
            k_theta: 150.0,
            k_omega: 10.0,
            k_oop: 20.0,
            k_bounds: 200.0,
            k_chiral: 50.0,
        };
        coords3d_min =
            minimize_energy_lbfgs(&our_mol, &coords3d_min, &smoothed, &params, 100, 1e-4);

        // 3. Calculate Minimized RMSD
        let mut min_diff_sq_sum = 0.0;
        pairs = 0;
        for i in 0..n {
            for j in (i + 1)..n {
                let dx_r = mol.atoms[i].x - mol.atoms[j].x;
                let dy_r = mol.atoms[i].y - mol.atoms[j].y;
                let dz_r = mol.atoms[i].z - mol.atoms[j].z;
                let rdkit_dist = (dx_r * dx_r + dy_r * dy_r + dz_r * dz_r).sqrt();

                let dx_our = coords3d_min[(i, 0)] - coords3d_min[(j, 0)];
                let dy_our = coords3d_min[(i, 1)] - coords3d_min[(j, 1)];
                let dz_our = coords3d_min[(i, 2)] - coords3d_min[(j, 2)];
                let our_dist = (dx_our * dx_our + dy_our * dy_our + dz_our * dz_our).sqrt();

                let diff = (our_dist - rdkit_dist).abs();
                min_diff_sq_sum += diff * diff;
                pairs += 1;
            }
        }
        if pairs > 0 {
            min_total_rmsd += (min_diff_sq_sum / pairs as f32).sqrt();
        }
    }
    let rust_time = pb.elapsed();

    let avg_raw_rmsd = raw_total_rmsd / bench_data.count as f32;
    let avg_min_rmsd = min_total_rmsd / bench_data.count as f32;
    let rust_time_ms = rust_time.as_secs_f64() * 1000.0;

    println!("=== FULL PIPELINE BENCHMARK: SMILES -> 3D -> MIN ===");
    println!("Molecules Processed: {}", bench_data.count);
    println!("Raw Distance RMSD (vs RDKit): {:.3} Å", avg_raw_rmsd);
    println!("Final Distance RMSD (vs RDKit): {:.3} Å", avg_min_rmsd);
    println!("Total Time (Rust): {:.2} ms", rust_time_ms);
    println!(
        "Average time per molecule: {:.2} ms",
        rust_time_ms / bench_data.count as f64
    );
}
