use nalgebra::{DMatrix, Matrix3};
use sci_form::graph::Molecule;
use std::fs;
use std::io::{BufReader, BufRead};
use rand::SeedableRng;
use rand::seq::SliceRandom;

mod oracle {
    use serde::Deserialize;
    #[derive(Deserialize)]
    pub struct OracleAtom {
        #[allow(dead_code)]
        pub element: u8,
        pub x: f32,
        pub y: f32,
        pub z: f32,
    }
    #[derive(Deserialize)]
    pub struct OracleMolecule {
        pub smiles: String,
        pub atoms: Vec<OracleAtom>,
    }
}

fn calculate_rmsd_icp_refined(coords: &DMatrix<f32>, reference: &DMatrix<f32>, _elements: &[u8]) -> f32 {
    let n = coords.nrows();
    if n == 0 { return 0.0; }
    
    // Simple element-based matching: group atoms by element type
    // and match by element, then by distance
    let mut c1 = nalgebra::Vector3::zeros();
    let mut c2 = nalgebra::Vector3::zeros();
    for i in 0..n {
        c1 += nalgebra::Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
        c2 += nalgebra::Vector3::new(reference[(i, 0)], reference[(i, 1)], reference[(i, 2)]);
    }
    c1 /= n as f32;
    c2 /= n as f32;
    
    // Center coordinates
    let mut coords_centered = DMatrix::from_element(n, 3, 0.0);
    let mut ref_centered = DMatrix::from_element(n, 3, 0.0);
    for i in 0..n {
        coords_centered[(i, 0)] = coords[(i, 0)] - c1[0];
        coords_centered[(i, 1)] = coords[(i, 1)] - c1[1];
        coords_centered[(i, 2)] = coords[(i, 2)] - c1[2];
        ref_centered[(i, 0)] = reference[(i, 0)] - c2[0];
        ref_centered[(i, 1)] = reference[(i, 1)] - c2[1];
        ref_centered[(i, 2)] = reference[(i, 2)] - c2[2];
    }
    
    // Kabsch SVD
    let mut h = Matrix3::zeros();
    for i in 0..n {
        let p = nalgebra::Vector3::new(coords_centered[(i, 0)], coords_centered[(i, 1)], coords_centered[(i, 2)]);
        let q = nalgebra::Vector3::new(ref_centered[(i, 0)], ref_centered[(i, 1)], ref_centered[(i, 2)]);
        h += p * q.transpose();
    }
    
    let svd = h.svd(true, true);
    if let (Some(u), Some(vt)) = (svd.u, svd.v_t) {
        let det = (u * vt).determinant();
        if det < 0.0 {
            // Reflection case - just use absolute value
            let mut sum = 0.0;
            for i in 0..n {
                let d = (coords_centered.row(i) - ref_centered.row(i)).norm();
                sum += d * d;
            }
            return (sum / n as f32).sqrt();
        }
    }
    
    // Direct Kabsch RMSD
    sci_form::forcefield::minimizer::calculate_rmsd_kabsch(&coords, &reference)
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
    let all_reference_mols: Vec<oracle::OracleMolecule> = fs::read_to_string(ref_path)
        .ok()
        .and_then(|s| serde_json::from_str(&s).ok())
        .unwrap_or_default();
    
    // Seleccionar 100 moléculas aleatorias
    let mut rng = rand::rngs::StdRng::seed_from_u64(42);
    let mut indices: Vec<usize> = (0..all_reference_mols.len()).collect();
    let _ = indices.shuffle(&mut rng);
    indices.truncate(100);
    
    // Crear lista de SMILES a procesar
    let mut test_smiles: Vec<String> = Vec::new();
    for idx in &indices {
        if *idx < smiles_list.len() {
            test_smiles.push(smiles_list[*idx].clone());
        }
    }
    
    // Crear mapa de reference para acceso rápido
    let reference_map: std::collections::HashMap<String, &oracle::OracleMolecule> = 
        all_reference_mols.iter().map(|m| (m.smiles.clone(), m)).collect();

    let mut max_rmsd = 0.0;
    let mut avg_rmsd = 0.0;
    let mut count = 0;

    for smi in &test_smiles {
        let mol_result = Molecule::from_smiles(smi);
        if let Ok(mol) = mol_result {
            let n = mol.graph.node_count();
            let mut bounds = sci_form::distgeom::calculate_bounds_matrix(&mol);
            sci_form::distgeom::triangle_smooth(&mut bounds);
            let smoothed = bounds;
            
            let chiral_sets = sci_form::distgeom::identify_chiral_sets(&mol);
            
            let params = sci_form::forcefield::FFParams {
                kb: 2000.0,
                k_theta: 1000.0,
                k_omega: 100.0,
                k_oop: 200.0,
                k_bounds: 1000.0,
                k_chiral: 500.0,
                k_vdw: 0.0,
            };
            
            // Fast mode: 1 solo intento
            let seed = 42 + (smi.len() as u64 * 7);
            let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
            
            let dists = sci_form::distgeom::pick_etkdg_distances(&mut rng, &smoothed);
            let metric = sci_form::distgeom::compute_metric_matrix(&dists);
            let mut coords4d = sci_form::distgeom::generate_nd_coordinates(&mut rng, &metric, 4);
            
            sci_form::forcefield::bounds_ff::minimize_bounds_lbfgs(
                &mut coords4d,
                &smoothed,
                &chiral_sets,
                500,
                1e-5,
            );

            let mut coords3d = DMatrix::from_element(n, 3, 0.0);
            for i in 0..n {
                for d in 0..3 {
                    coords3d[(i, d)] = coords4d[(i, d)];
                }
            }

            coords3d = sci_form::forcefield::minimizer::minimize_energy_lbfgs(
                &mol, &coords3d, &smoothed, &params, 500, 1e-4,
            );

            if let Some(oracle) = reference_map.get(smi) {
                let oracle = *oracle;
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
    if count > 0 {
        println!(
            "Final results: Max RMSD: {:.3} Å, Avg RMSD: {:.3} Å",
            max_rmsd,
            avg_rmsd / count as f32
        );
    }
}
