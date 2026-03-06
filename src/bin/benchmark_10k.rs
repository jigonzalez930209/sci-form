use sci_form::forcefield::FFParams;
use sci_form::graph::Molecule;
use serde::Deserialize;
use std::fs;
use std::io::{BufRead, BufReader};
use std::time::Instant;

#[derive(Deserialize, Debug)]
struct OracleAtom {
    element: u8,
    x: f32,
    y: f32,
    z: f32,
}

#[derive(Deserialize, Debug)]
struct OracleMolecule {
    smiles: String,
    atoms: Vec<OracleAtom>,
}

fn compute_rmsd(coords: &nalgebra::DMatrix<f32>, oracle: &[OracleAtom], n: usize) -> f32 {
    let mut diff_sq_sum = 0.0;
    let mut pairs = 0;
    for i in 0..n {
        for j in (i + 1)..n {
            let dx_r = oracle[i].x - oracle[j].x;
            let dy_r = oracle[i].y - oracle[j].y;
            let dz_r = oracle[i].z - oracle[j].z;
            let rdkit_dist = (dx_r * dx_r + dy_r * dy_r + dz_r * dz_r).sqrt();

            let dx = coords[(i, 0)] - coords[(j, 0)];
            let dy = coords[(i, 1)] - coords[(j, 1)];
            let dz = coords[(i, 2)] - coords[(j, 2)];
            let our_dist = (dx * dx + dy * dy + dz * dz).sqrt();

            let diff = (our_dist - rdkit_dist).abs();
            diff_sq_sum += diff * diff;
            pairs += 1;
        }
    }
    if pairs > 0 {
        (diff_sq_sum / pairs as f32).sqrt()
    } else {
        f32::MAX
    }
}

fn main() {
    let smi_file = "scripts/10k_smiles.smi";
    let file = fs::File::open(smi_file).expect("Failed to open SMILES file");
    let reader = BufReader::new(file);
    let smiles_list: Vec<String> = reader
        .lines()
        .map(|l| l.unwrap().trim().to_string())
        .filter(|s| !s.is_empty())
        .collect();

    let ref_path = "tests/fixtures/rdkit_10k_reference.json";
    let reference_mols: Option<Vec<OracleMolecule>> = fs::read_to_string(ref_path)
        .ok()
        .and_then(|s| serde_json::from_str(&s).ok());

    println!("Starting Sci-Form FULL 10k benchmark (conformer pool)...");
    let start_total = Instant::now();

    let mut count = 0;
    let mut total_rmsd = 0.0;
    let mut max_rmsd: f32 = 0.0;
    let mut max_rmsd_smi = String::new();
    let mut rmsd_count = 0;

    let params = FFParams {
        kb: 300.0,
        k_theta: 200.0,
        k_omega: 15.0,
        k_oop: 50.0,
        k_bounds: 100.0,
        k_chiral: 50.0,
    };
    let num_confs = 50;
    let lbfgs_iters = 50;

    for smi in &smiles_list {
        let mol_result = Molecule::from_smiles(smi);
        if let Ok(mol) = mol_result {
            let n = mol.graph.node_count();
            let bounds = sci_form::distgeom::calculate_bounds_matrix(&mol);
            let smoothed = sci_form::distgeom::smooth_bounds_matrix(bounds);

            let mut rng = rand::rngs::StdRng::seed_from_u64(42);
            use rand::SeedableRng;

            let mut best_rmsd = f32::MAX;

            for _ in 0..num_confs {
                let dists = sci_form::distgeom::pick_random_distances(&mut rng, &smoothed);
                let metric = sci_form::distgeom::compute_metric_matrix(&dists);
                let mut coords3d = sci_form::distgeom::generate_3d_coordinates(&mut rng, &metric);

                coords3d = sci_form::forcefield::minimizer::minimize_energy_lbfgs(
                    &mol,
                    &coords3d,
                    &smoothed,
                    &params,
                    lbfgs_iters,
                    1e-4,
                );

                if let Some(ref ref_list) = reference_mols {
                    if let Some(oracle) = ref_list.iter().find(|m| m.smiles == *smi) {
                        if oracle.atoms.len() == n {
                            let rmsd = compute_rmsd(&coords3d, &oracle.atoms, n);
                            if rmsd < best_rmsd {
                                best_rmsd = rmsd;
                            }
                        }
                    }
                }
            }

            if best_rmsd < f32::MAX {
                total_rmsd += best_rmsd;
                if best_rmsd > max_rmsd {
                    max_rmsd = best_rmsd;
                    max_rmsd_smi = smi.clone();
                    println!("New MAX Min-RMSD {:.3} Å for molecule {}", max_rmsd, smi);
                }
                rmsd_count += 1;
            }
            count += 1;
            if count % 500 == 0 {
                let elapsed = start_total.elapsed().as_secs_f64();
                let avg = if rmsd_count > 0 {
                    total_rmsd / rmsd_count as f32
                } else {
                    0.0
                };
                println!(
                    "[{:>5}] Avg RMSD: {:.3} Å | Max: {:.3} Å ({}) | {:.1}s",
                    count, avg, max_rmsd, max_rmsd_smi, elapsed
                );
            }
        }
    }

    let total_time = start_total.elapsed();
    let total_ms = total_time.as_secs_f64() * 1000.0;

    println!("\n=== Sci-Form FULL 10k Benchmark Results ===");
    println!("Molecules Successfully Processed: {}", count);
    println!("Total Time: {:.2} ms", total_ms);
    println!("Average Time: {:.2} ms/mol", total_ms / count as f64);
    if rmsd_count > 0 {
        println!(
            "Average RMSD (vs RDKit reference): {:.3} Å (over {} molecules)",
            total_rmsd / rmsd_count as f32,
            rmsd_count
        );
        println!(
            "MAX RMSD (vs RDKit reference): {:.3} Å ({})",
            max_rmsd, max_rmsd_smi
        );
    }
}
