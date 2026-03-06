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

    println!("Starting Sci-Form fast benchmark (single conformer, bounds FF)...");
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
                // println!("WARNING: triangle_smooth failed for {}", smi);
            }
            let smoothed = bounds;

            let mut rng = rand::rngs::StdRng::seed_from_u64(42);
            use rand::SeedableRng;

            let mut best_rmsd = f32::MAX;
            let mut found_match = false;

            let chiral_sets = sci_form::distgeom::identify_chiral_sets(&mol);
            let num_confs = 1;
            for _c in 0..num_confs {
                let dists = sci_form::distgeom::pick_random_distances(&mut rng, &smoothed);
                let metric = sci_form::distgeom::compute_metric_matrix(&dists);
                
                // Phase 1: 4D Bounds FF Minimization
                let mut coords = sci_form::distgeom::generate_nd_coordinates(&mut rng, &metric, 4);
                sci_form::forcefield::bounds_ff::minimize_bounds_lbfgs(
                    &mut coords, &smoothed, &chiral_sets, 500, 1e-4,
                );

                // Project to 3D for Phase 2
                let mut coords3d = nalgebra::DMatrix::from_element(n, 3, 0.0);
                for i in 0..n {
                    for d in 0..3 {
                        coords3d[(i, d)] = coords[(i, d)];
                    }
                }

                // Phase 2: ETKDG-Lite MMFF Minimization (M6 torsions included)
                let params = sci_form::forcefield::FFParams {
                    kb: 400.0,
                    k_theta: 300.0,
                    k_omega: 20.0, // This is now used alongside M6
                    k_oop: 40.0,
                    k_bounds: 10000.0,
                };
                sci_form::forcefield::minimizer::minimize_energy_lbfgs(
                    &mut coords3d, &mol, &params, &smoothed, 100,
                );

                if let Some(ref ref_list) = reference_mols {
                    if let Some(oracle) = ref_list.iter().find(|m| m.smiles == *smi) {
                        if oracle.atoms.len() == n {
                            let mut diff_sq_sum = 0.0;
                            let mut pairs = 0;
                            for i in 0..n {
                                for j in (i + 1)..n {
                                    let dx_r = oracle.atoms[i].x - oracle.atoms[j].x;
                                    let dy_r = oracle.atoms[i].y - oracle.atoms[j].y;
                                    let dz_r = oracle.atoms[i].z - oracle.atoms[j].z;
                                    let rdkit_dist = (dx_r * dx_r + dy_r * dy_r + dz_r * dz_r).sqrt();
                                    let dx = coords3d[(i, 0)] - coords3d[(j, 0)];
                                    let dy = coords3d[(i, 1)] - coords3d[(j, 1)];
                                    let dz = coords3d[(i, 2)] - coords3d[(j, 2)];
                                    let our_dist = (dx * dx + dy * dy + dz * dz).sqrt();
                                    
                                    let diff = (our_dist - rdkit_dist).abs();
                                    diff_sq_sum += diff * diff;
                                    pairs += 1;
                                }
                            }
                            if pairs > 0 {
                                let rmsd = (diff_sq_sum / pairs as f32).sqrt();
                                if rmsd < best_rmsd {
                                    best_rmsd = rmsd;
                                    found_match = true;
                                }
                            }
                        }
                    }
                }
            }

            if found_match {
                total_rmsd += best_rmsd;
                if best_rmsd > max_rmsd {
                    max_rmsd = best_rmsd;
                    max_rmsd_smi = smi.clone();
                    println!("New MAX RMSD {:.3} Å for {}", max_rmsd, smi);
                }
                rmsd_count += 1;
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
