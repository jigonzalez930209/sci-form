use rand::SeedableRng;
use serde::Deserialize;
use std::fs;

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

#[test]
fn test_parse_reference_data() {
    let data = fs::read_to_string("tests/fixtures/reference_coords.json")
        .expect("Should be able to read reference_coords.json");
    let mut molecules: Vec<OracleMolecule> =
        serde_json::from_str(&data).expect("JSON was not well-formatted");

    assert!(!molecules.is_empty());

    // Shuffle and pick 100 random molecules
    use rand::seq::SliceRandom;
    let mut rng = rand::thread_rng();
    molecules.shuffle(&mut rng);
    molecules.truncate(100);

    let mut total_rmsd = 0.0;
    let mut count = 0;

    for mol in molecules {
        assert!(!mol.atoms.is_empty());
        let n = mol.atoms.len();

        // Convert OracleMolecule to our format
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

            let mut new_atom = sci_form::graph::Atom {
                element: atom.element,
                position: nalgebra::Vector3::zeros(), // Ignore RDKit coords when building initial graph
                charge: 0.0,                          // Not parsed yet
                formal_charge: atom.formal_charge,
                hybridization: hybridization,
                chiral_tag: sci_form::graph::ChiralType::Unspecified,
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

        // 1. Calculate Bounds
        let bounds = sci_form::distgeom::calculate_bounds_matrix(&our_mol);
        let smoothed = sci_form::distgeom::smooth_bounds_matrix(bounds.clone());

        // 2. Embed
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let dists = sci_form::distgeom::pick_random_distances(&mut rng, &smoothed);
        let metric = sci_form::distgeom::compute_metric_matrix(&dists);

        // Project 3D
        let mut coords3d = sci_form::distgeom::generate_3d_coordinates(&metric);

        let mut minimized_coords = coords3d.clone();
        let params = sci_form::forcefield::FFParams {
            kb: 0.0,
            k_theta: 0.0,
            k_omega: 0.0,
            k_oop: 0.0,
            k_bounds: 500.0, // Force L-BFGS to rigidly enforce ONLY the distgeom bounds
        };
        sci_form::forcefield::minimize_energy_lbfgs(&mut minimized_coords, &our_mol, &params, &smoothed, 300);

        // Ensure successful generation
        assert_eq!(coords3d.nrows(), n);
        assert_eq!(coords3d.ncols(), 3);

        // Compute pairwise distances for both sets to avoid Kabsch alignment
        let mut rdkit_dmatrix = nalgebra::DMatrix::from_element(n, n, 0.0);
        let mut our_dmatrix = nalgebra::DMatrix::from_element(n, n, 0.0);

        // Compare Raw Embed
        let mut diff_raw_sq_sum = 0.0;
        let mut diff_min_sq_sum = 0.0;
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

                let dx_min = minimized_coords[(i, 0)] - minimized_coords[(j, 0)];
                let dy_min = minimized_coords[(i, 1)] - minimized_coords[(j, 1)];
                let dz_min = minimized_coords[(i, 2)] - minimized_coords[(j, 2)];
                let min_dist = (dx_min * dx_min + dy_min * dy_min + dz_min * dz_min).sqrt();

                let diff_raw = rdkit_dist - raw_dist;
                let diff_min = rdkit_dist - min_dist;

                diff_raw_sq_sum += diff_raw * diff_raw;
                diff_min_sq_sum += diff_min * diff_min;
                pairs += 1;
            }
        }

        let dist_raw_rmsd = if pairs > 0 {
            (diff_raw_sq_sum / pairs as f32).sqrt()
        } else {
            0.0
        };
        let dist_min_rmsd = if pairs > 0 {
            (diff_min_sq_sum / pairs as f32).sqrt()
        } else {
            0.0
        };

        total_rmsd += dist_min_rmsd;
        // Print 5 random molecules to visually see raw vs min
        if count % 20 == 0 {
            println!(
                "Mol {} -> Raw RMSD: {:.3} Å | Min RMSD: {:.3} Å",
                count, dist_raw_rmsd, dist_min_rmsd
            );
        }

        count += 1;

        if count % 1000 == 0 {
            println!(
                "Processed {} molecules... Current Avg Error: {:.3} Å",
                count,
                total_rmsd / count as f32
            );
        }
    }

    let final_avg_rmsd = total_rmsd / count as f32;
    println!("=== TEST COMPLETE ===");
    println!("Successfully processed {} molecules.", count);
    println!(
        "Average Distance Matrix Error (vs RDKit): {:.3} Å",
        final_avg_rmsd
    );
}

// Function to calculate RMSD between two sets of coordinates.
// Future implementation when we generate 3D coordinates.
pub fn calculate_rmsd(coords1: &[(f32, f32, f32)], coords2: &[(f32, f32, f32)]) -> f32 {
    assert_eq!(coords1.len(), coords2.len());
    let mut sum_sq_diff = 0.0;

    for (p1, p2) in coords1.iter().zip(coords2.iter()) {
        let dx = p1.0 - p2.0;
        let dy = p1.1 - p2.1;
        let dz = p1.2 - p2.2;
        sum_sq_diff += dx * dx + dy * dy + dz * dz;
    }

    (sum_sq_diff / coords1.len() as f32).sqrt()
}

#[test]
fn test_rmsd_calculation_dummy() {
    let raw1 = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)];
    let raw2 = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.1)];
    let rmsd = calculate_rmsd(&raw1, &raw2);
    assert!(rmsd < 0.1);
}
