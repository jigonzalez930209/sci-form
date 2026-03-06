use nalgebra::{DMatrix, Matrix3, Vector3};
use sci_form::graph::Molecule;
use serde::{Deserialize, Serialize};
use std::fs;
use std::io::{BufRead, BufReader};

#[derive(Debug, Serialize, Deserialize)]
struct OracleAtom {
    element: u8,
    x: f32,
    y: f32,
    z: f32,
}

#[derive(Debug, Serialize, Deserialize)]
struct OracleMolecule {
    smiles: String,
    atoms: Vec<OracleAtom>,
}

fn calculate_rmsd_icp_refined(
    coords: &DMatrix<f32>,
    reference: &DMatrix<f32>,
    elements: &[u8],
) -> f32 {
    let n = coords.nrows();
    if n == 0 {
        return 0.0;
    }

    let mut current_coords = coords.clone();
    let mut best_rmsd = 1e10;

    // Try a few initial rotations (Identity, and 90-degree flips) to avoid local minima
    for _trial in 0..1 {
        let mut mapping = vec![0; n];

        for _iter in 0..20 {
            // 1. Find best mapping based on elements and distance
            let mut used = vec![false; n];
            for i in 0..n {
                let p = Vector3::new(
                    current_coords[(i, 0)],
                    current_coords[(i, 1)],
                    current_coords[(i, 2)],
                );
                let mut min_d2 = 1e10;
                let mut best_j = 0;
                for j in 0..n {
                    if elements[i] == elements[j] && !used[j] {
                        let q =
                            Vector3::new(reference[(j, 0)], reference[(j, 1)], reference[(j, 2)]);
                        let d2 = (p - q).norm_squared();
                        if d2 < min_d2 {
                            min_d2 = d2;
                            best_j = j;
                        }
                    }
                }
                mapping[i] = best_j;
                used[best_j] = true;
            }

            // 2. Kabsch alignment
            let mut mapped_ref = DMatrix::from_element(n, 3, 0.0);
            for i in 0..n {
                let j = mapping[i];
                mapped_ref[(i, 0)] = reference[(j, 0)];
                mapped_ref[(i, 1)] = reference[(j, 1)];
                mapped_ref[(i, 2)] = reference[(j, 2)];
            }

            let rmsd = sci_form::forcefield::minimizer::calculate_rmsd_kabsch(
                &current_coords,
                &mapped_ref,
            );
            if rmsd < best_rmsd {
                best_rmsd = rmsd;
            }

            // For the next iteration, we need the actual coordinates to be aligned
            // But calculate_rmsd_kabsch doesn't return the rotation.
            // Let's just do a manual Kabsch to update current_coords.
            let mut c1 = Vector3::zeros();
            let mut c2 = Vector3::zeros();
            for i in 0..n {
                c1 += Vector3::new(
                    current_coords[(i, 0)],
                    current_coords[(i, 1)],
                    current_coords[(i, 2)],
                );
                c2 += Vector3::new(mapped_ref[(i, 0)], mapped_ref[(i, 1)], mapped_ref[(i, 2)]);
            }
            c1 /= n as f32;
            c2 /= n as f32;
            let mut h = Matrix3::zeros();
            for i in 0..n {
                let p = Vector3::new(
                    current_coords[(i, 0)],
                    current_coords[(i, 1)],
                    current_coords[(i, 2)],
                ) - c1;
                let q =
                    Vector3::new(mapped_ref[(i, 0)], mapped_ref[(i, 1)], mapped_ref[(i, 2)]) - c2;
                h += p * q.transpose();
            }
            let svd = h.svd(true, true);
            let u = svd.u.unwrap();
            let vt = svd.v_t.unwrap();
            let mut d_mat = Matrix3::identity();
            if (vt.transpose() * u.transpose()).determinant() < 0.0 {
                d_mat[(2, 2)] = -1.0;
            }
            let r = vt.transpose() * d_mat * u.transpose();

            for i in 0..n {
                let p = Vector3::new(
                    current_coords[(i, 0)],
                    current_coords[(i, 1)],
                    current_coords[(i, 2)],
                ) - c1;
                let new_p = (r * p) + c2;
                current_coords[(i, 0)] = new_p.x;
                current_coords[(i, 1)] = new_p.y;
                current_coords[(i, 2)] = new_p.z;
            }
        }
    }
    best_rmsd
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

    let mut max_rmsd = 0.0;
    let mut avg_rmsd = 0.0;
    let mut count = 0;

    for smi in &smiles_list {
        let mol_result = Molecule::from_smiles(smi);
        if let Ok(mol) = mol_result {
            let n = mol.graph.node_count();
            let mut bounds = sci_form::distgeom::calculate_bounds_matrix(&mol);
            sci_form::distgeom::triangle_smooth(&mut bounds);
            let smoothed = bounds;
            let mut rng = rand::rngs::StdRng::seed_from_u64(42);
            use rand::SeedableRng;
            let chiral_sets = sci_form::distgeom::identify_chiral_sets(&mol);
            let dists = sci_form::distgeom::pick_random_distances(&mut rng, &smoothed);
            let metric = sci_form::distgeom::compute_metric_matrix(&dists);
            let mut coords4d = sci_form::distgeom::generate_nd_coordinates(&mut rng, &metric, 4);
            sci_form::forcefield::bounds_ff::minimize_bounds_lbfgs(
                &mut coords4d,
                &smoothed,
                &chiral_sets,
                1000,
                1e-7,
            );

            let mut coords3d = DMatrix::from_element(n, 3, 0.0);
            for i in 0..n {
                for d in 0..3 {
                    coords3d[(i, d)] = coords4d[(i, d)];
                }
            }

            let params = sci_form::forcefield::FFParams {
                kb: 800.0,
                k_theta: 500.0,
                k_omega: 50.0,
                k_oop: 80.0,
                k_bounds: 400.0,
                k_chiral: 200.0,
            };
            coords3d = sci_form::forcefield::minimizer::minimize_energy_lbfgs(
                &mol, &coords3d, &smoothed, &params, 2000, 1e-5,
            );

            if let Some(ref ref_list) = reference_mols {
                if let Some(oracle) = ref_list.iter().find(|m| m.smiles == *smi) {
                    if oracle.atoms.len() == n {
                        let mut ref_coords = DMatrix::from_element(n, 3, 0.0);
                        let mut elements = Vec::new();
                        for i in 0..n {
                            ref_coords[(i, 0)] = oracle.atoms[i].x;
                            ref_coords[(i, 1)] = oracle.atoms[i].y;
                            ref_coords[(i, 2)] = oracle.atoms[i].z;
                            elements.push(mol.graph[petgraph::graph::NodeIndex::new(i)].element);
                        }

                        let rmsd = calculate_rmsd_icp_refined(&coords3d, &ref_coords, &elements);
                        if rmsd > max_rmsd {
                            max_rmsd = rmsd;
                            println!("New MAX RMSD {:.3} Å for {}", max_rmsd, smi);
                        }
                        avg_rmsd += rmsd;
                        count += 1;

                        if smi == "CC#C" {
                            println!("DEBUG: CC#C | RMSD: {:.3} Å", rmsd);
                            println!("DEBUG CC#C Coords generated:");
                            for i in 0..n {
                                println!("  Atom {}: ({:.3}, {:.3}, {:.3})", 
                                    i, coords3d[(i, 0)], coords3d[(i, 1)], coords3d[(i, 2)]);
                            }
                            println!("DEBUG CC#C Reference coords:");
                            for i in 0..n {
                                println!("  Atom {}: ({:.3}, {:.3}, {:.3})", 
                                    i, ref_coords[(i, 0)], ref_coords[(i, 1)], ref_coords[(i, 2)]);
                            }
                            let ana = sci_form::forcefield::gradients::compute_analytical_gradient(
                                &coords3d, &mol, &params, &smoothed,
                            );
                            println!("DEBUG CC#C Grad Norm at end: {:.6}", ana.norm());
                        }
                    }
                }
            }
        }
    }
    if count > 0 {
        println!(
            "Final results: Max RMSD: {:.3} Å, Avg RMSD: {:.3} Å",
            max_rmsd,
            avg_rmsd / count as f32
        );
    }
}
