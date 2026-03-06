use sci_form::graph::Molecule;
use serde::Deserialize;
use std::fs;
use std::io::{BufRead, BufReader};
use std::time::Instant;
use nalgebra::DMatrix;

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

fn calculate_rmsd(coords: &DMatrix<f32>, reference: &DMatrix<f32>) -> f32 {
    let n = coords.nrows();
    if n == 0 { return 0.0; }
    
    // Centroid alignment
    let mut c1 = [0.0; 3];
    let mut c2 = [0.0; 3];
    for i in 0..n {
        for d in 0..3 {
            c1[d] += coords[(i, d)];
            c2[d] += reference[(i, d)];
        }
    }
    for d in 0..3 {
        c1[d] /= n as f32;
        c2[d] /= n as f32;
    }

    let mut sum_sq = 0.0;
    for i in 0..n {
        let dx = (coords[(i, 0)] - c1[0]) - (reference[(i, 0)] - c2[0]);
        let dy = (coords[(i, 1)] - c1[1]) - (reference[(i, 1)] - c2[1]);
        let dz = (coords[(i, 2)] - c1[2]) - (reference[(i, 2)] - c2[2]);
        sum_sq += dx * dx + dy * dy + dz * dz;
    }
    (sum_sq / n as f32).sqrt()
}

fn main() {
    let smi_file = "scripts/10k_smiles.smi";
    let file = fs::File::open(smi_file).expect("Failed to open SMILES file");
    let reader = BufReader::new(file);
    let mut smiles_list: Vec<String> = reader
        .lines()
        .map(|l| l.unwrap().trim().to_string())
        .filter(|s| !s.is_empty())
        .collect();
    smiles_list.truncate(100);

    let ref_path = "tests/fixtures/rdkit_10k_reference.json";
    let reference_mols: Option<Vec<OracleMolecule>> = fs::read_to_string(ref_path)
        .ok()
        .and_then(|s| serde_json::from_str(&s).ok());

    println!("Starting Sci-Form fast benchmark (single conformer, 4D->3D ETKDG)...");
    let start_total = Instant::now();

    let mut count = 0;
    let mut total_rmsd = 0.0;
    let mut max_rmsd: f32 = 0.0;
    let mut max_rmsd_smi = String::new();
    let mut rmsd_count = 0;

    for smi in &smiles_list {
        let mol_result = Molecule::from_smiles(smi);
        if let Ok(mol) = mol_result {
            let n = mol.graph.node_count();
            let mut bounds = sci_form::distgeom::calculate_bounds_matrix(&mol);
            if !sci_form::distgeom::triangle_smooth(&mut bounds) {
                println!("WARNING: triangle_smooth failed for {}", smi);
            }
            let smoothed = bounds;

            let mut rng = rand::rngs::StdRng::seed_from_u64(42);
            use rand::SeedableRng;

            let chiral_sets = sci_form::distgeom::identify_chiral_sets(&mol);
            
            let dists = sci_form::distgeom::pick_random_distances(&mut rng, &smoothed);
            let metric = sci_form::distgeom::compute_metric_matrix(&dists);
            
            // Phase 1: 4D Bounds FF Minimization
            let mut coords4d = sci_form::distgeom::generate_nd_coordinates(&mut rng, &metric, 4);
            sci_form::forcefield::bounds_ff::minimize_bounds_lbfgs(
                &mut coords4d, &smoothed, &chiral_sets, 500, 1e-4,
            );

            // Project to 3D for Phase 2
            let mut coords3d = DMatrix::from_element(n, 3, 0.0);
            for i in 0..n {
                for d in 0..3 {
                    coords3d[(i, d)] = coords4d[(i, d)];
                }
            }

            // Phase 2: ETKDG-Lite MMFF Minimization
            let params = sci_form::forcefield::FFParams {
                kb: 400.0,
                k_theta: 300.0,
                k_omega: 20.0,
                k_oop: 40.0,
                k_bounds: 10000.0,
            };
            sci_form::forcefield::minimizer::minimize_energy_lbfgs(
                &mut coords3d, &mol, &params, &smoothed, 500,
            );

            if let Some(ref ref_list) = reference_mols {
                if let Some(oracle) = ref_list.iter().find(|m| m.smiles == *smi) {
                    if oracle.atoms.len() == n {
                        let mut ref_coords = DMatrix::from_element(n, 3, 0.0);
                        for i in 0..n {
                            ref_coords[(i, 0)] = oracle.atoms[i].x;
                            ref_coords[(i, 1)] = oracle.atoms[i].y;
                            ref_coords[(i, 2)] = oracle.atoms[i].z;
                        }

                        let rmsd = calculate_rmsd(&coords3d, &ref_coords);
                        
                        if smi == "C" {
                            println!("DEBUG: C | RMSD: {:.3} Å", rmsd);
                            println!("Coords:\n{:?}", coords3d);
                            println!("Ref:\n{:?}", ref_coords);
                        }
                        
                        total_rmsd += rmsd;
                        if rmsd > max_rmsd {
                            max_rmsd = rmsd;
                            max_rmsd_smi = smi.clone();
                            println!("New MAX RMSD {:.3} Å for {}", max_rmsd, smi);
                        }
                        rmsd_count += 1;
                    }
                }
            }
            count += 1;
        }
    }

    let total_time = start_total.elapsed();
    let total_ms = total_time.as_secs_f64() * 1000.0;

    println!("\n=== Sci-Form Fast 100 Benchmark (Single Conformer) ===");
    println!("Molecules Processed: {}", count);
    println!("Average Time: {:.2} ms/mol", total_ms / count as f64);
    if rmsd_count > 0 {
        println!(
            "Average RMSD: {:.3} Å (over {} molecules)",
            total_rmsd / rmsd_count as f32,
            rmsd_count
        );
        println!("MAX RMSD: {:.3} Å ({})", max_rmsd, max_rmsd_smi);
    }
}
