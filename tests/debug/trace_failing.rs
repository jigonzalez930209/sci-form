#![allow(
    unused_imports,
    unused_variables,
    dead_code,
    clippy::unnecessary_cast,
    clippy::needless_range_loop,
    clippy::manual_repeat_n,
    clippy::manual_str_repeat,
    clippy::manual_is_multiple_of,
    clippy::redundant_field_names,
    clippy::useless_vec,
    clippy::single_range_in_vec_init
)]
//! Step-by-step pipeline comparison for a failing molecule.
//! Loads RDKit's bounds matrix and final coordinates, then runs our pipeline
//! and prints diagnostics at each step to find where divergence occurs.

use nalgebra::DMatrix;
use sci_form::distgeom::{
    calculate_bounds_matrix_opts, compute_initial_coords_rdkit, identify_chiral_sets,
    identify_tetrahedral_centers, pick_rdkit_distances, triangle_smooth_tol, MinstdRand,
};
use sci_form::forcefield::bounds_ff::minimize_bfgs_rdkit;
use sci_form::forcefield::etkdg_3d::{build_etkdg_3d_ff_with_torsions, minimize_etkdg_3d_bfgs};

fn load_npy_f64(path: &str) -> DMatrix<f64> {
    let data = std::fs::read(path).expect("Failed to read npy file");
    // Parse numpy .npy format: 10-byte magic, then header, then data
    // Find the end of header (terminated by \n)
    let header_end = data.iter().position(|&b| b == b'\n').unwrap();
    // Parse shape from header
    let header = std::str::from_utf8(&data[10..header_end]).unwrap();
    let shape_start = header.find("(").unwrap() + 1;
    let shape_end = header.find(")").unwrap();
    let shape_str = &header[shape_start..shape_end];
    let dims: Vec<usize> = shape_str
        .split(',')
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.trim().parse().unwrap())
        .collect();
    let nrows = dims[0];
    let ncols = if dims.len() > 1 { dims[1] } else { 1 };

    let data_start = header_end + 1;
    let f64_data: Vec<f64> = data[data_start..]
        .chunks(8)
        .map(|chunk| f64::from_le_bytes(chunk.try_into().unwrap()))
        .collect();

    // numpy stores row-major, DMatrix is column-major
    let mut mat = DMatrix::zeros(nrows, ncols);
    for i in 0..nrows {
        for j in 0..ncols {
            mat[(i, j)] = f64_data[i * ncols + j];
        }
    }
    mat
}

#[test]
fn test_trace_failing_molecule() {
    let smiles = "C#CCOC(C)CC1CC2C3CCC(C)C(O)(C3)C2O1";

    // Build molecule from reference JSON (matching RDKit's graph exactly)
    let ref_data = sci_form::fixture_io::read_text_fixture("tests/fixtures/gdb20_reference_1k.json")
        .expect("Reference file needed");
    let ref_mols: Vec<serde_json::Value> = serde_json::from_str(&ref_data).unwrap();
    let ref_mol = ref_mols
        .iter()
        .find(|m| m["smiles"].as_str().unwrap() == smiles)
        .expect("Molecule not found in reference");

    let mut mol = sci_form::graph::Molecule::new(smiles);
    let atoms = ref_mol["atoms"].as_array().unwrap();
    let bonds = ref_mol["bonds"].as_array().unwrap();
    let mut node_indices = Vec::new();
    for atom in atoms {
        let element = atom["element"].as_u64().unwrap() as u8;
        let hyb_str = atom["hybridization"].as_str().unwrap();
        let hybridization = match hyb_str {
            "SP" => sci_form::graph::Hybridization::SP,
            "SP2" => sci_form::graph::Hybridization::SP2,
            "SP3" => sci_form::graph::Hybridization::SP3,
            "SP3D" => sci_form::graph::Hybridization::SP3D,
            "SP3D2" => sci_form::graph::Hybridization::SP3D2,
            _ => sci_form::graph::Hybridization::Unknown,
        };
        let fc = atom["formal_charge"].as_i64().unwrap_or(0) as i8;
        let new_atom = sci_form::graph::Atom {
            element,
            position: nalgebra::Vector3::zeros(),
            charge: 0.0,
            formal_charge: fc,
            hybridization,
            chiral_tag: sci_form::graph::ChiralType::Unspecified,
            explicit_h: if element == 1 || element == 0 { 1 } else { 0 },
        };
        node_indices.push(mol.add_atom(new_atom));
    }
    for bond in bonds {
        let start = bond["start"].as_u64().unwrap() as usize;
        let end = bond["end"].as_u64().unwrap() as usize;
        let order_str = bond["order"].as_str().unwrap();
        let order = match order_str {
            "DOUBLE" => sci_form::graph::BondOrder::Double,
            "TRIPLE" => sci_form::graph::BondOrder::Triple,
            "AROMATIC" => sci_form::graph::BondOrder::Aromatic,
            _ => sci_form::graph::BondOrder::Single,
        };
        mol.add_bond(
            node_indices[start],
            node_indices[end],
            sci_form::graph::Bond {
                order,
                stereo: sci_form::graph::BondStereo::None,
            },
        );
    }

    let n = mol.graph.node_count();
    println!("\nMolecule: {} (built from reference JSON)", smiles);
    println!("Atoms: {}", n);

    // Load RDKit bounds
    let rdkit_bounds = load_npy_f64("/tmp/rdkit_bounds_fail.npy");
    println!(
        "RDKit bounds loaded: {}x{}",
        rdkit_bounds.nrows(),
        rdkit_bounds.ncols()
    );

    // Our bounds
    let raw = calculate_bounds_matrix_opts(&mol, true);
    let mut bounds = raw;
    let ok = triangle_smooth_tol(&mut bounds, 0.0);
    println!("Our smoothing OK: {}", ok);

    // Compare bounds
    let mut max_ub_diff = 0.0f64;
    let mut max_lb_diff = 0.0f64;
    let mut ub_diff_count = 0;
    let mut lb_diff_count = 0;
    for i in 0..n {
        for j in (i + 1)..n {
            let our_ub = bounds[(i, j)];
            let rdkit_ub = rdkit_bounds[(i, j)];
            let diff = (our_ub - rdkit_ub).abs();
            if diff > 1e-10 {
                ub_diff_count += 1;
                if diff > max_ub_diff {
                    max_ub_diff = diff;
                }
            }
            let our_lb = bounds[(j, i)];
            let rdkit_lb = rdkit_bounds[(j, i)];
            let diff = (our_lb - rdkit_lb).abs();
            if diff > 1e-10 {
                lb_diff_count += 1;
                if diff > max_lb_diff {
                    max_lb_diff = diff;
                }
            }
        }
    }
    println!("\n=== Bounds Comparison ===");
    println!(
        "UB diffs (>1e-10): {} max: {:.2e}",
        ub_diff_count, max_ub_diff
    );
    println!(
        "LB diffs (>1e-10): {} max: {:.2e}",
        lb_diff_count, max_lb_diff
    );

    // Run our pipeline
    let chiral_sets = identify_chiral_sets(&mol);
    let _tet_centers = identify_tetrahedral_centers(&mol);
    let use_4d = !chiral_sets.is_empty();
    let embed_dim = if use_4d { 4 } else { 3 };
    println!(
        "\nChiral sets: {}, use_4d: {}, embed_dim: {}",
        chiral_sets.len(),
        use_4d,
        embed_dim
    );

    let mut rng = MinstdRand::new(42);

    // Step 1: Random distances
    let dists = pick_rdkit_distances(&mut rng, &bounds);
    println!("\n=== Random Distances (first 10 pairs) ===");
    for i in 0..n.min(10) {
        for j in (i + 1)..n.min(i + 3) {
            println!("  d({},{}) = {:.10}", i, j, dists[(i, j)]);
        }
    }

    // Also compute what distances RDKit would generate using same bounds
    let mut rng2 = MinstdRand::new(42);
    let dists_from_rdkit_bounds = pick_rdkit_distances(&mut rng2, &rdkit_bounds);
    println!("\n=== Distance comparison (our bounds vs RDKit bounds, first mismatches) ===");
    let mut dist_diff_count = 0;
    let mut max_dist_diff = 0.0f64;
    for i in 0..n {
        for j in (i + 1)..n {
            let diff = (dists[(i, j)] - dists_from_rdkit_bounds[(i, j)]).abs();
            if diff > 1e-15 {
                dist_diff_count += 1;
                if diff > max_dist_diff {
                    max_dist_diff = diff;
                }
                if dist_diff_count <= 5 {
                    println!(
                        "  d({},{}) ours={:.12} rdkit_b={:.12} diff={:.2e}",
                        i,
                        j,
                        dists[(i, j)],
                        dists_from_rdkit_bounds[(i, j)],
                        diff
                    );
                }
            }
        }
    }
    println!(
        "Total distance diffs: {} max: {:.2e}",
        dist_diff_count, max_dist_diff
    );

    // Step 2: Compute metric matrix and eigen decomposition (to print eigenvalues)
    {
        let n = mol.graph.node_count();
        let d_size = n * (n + 1) / 2;
        let mut sq_packed = vec![0.0f64; d_size];
        let mut sum_sq_all = 0.0f64;
        for i in 0..n {
            let id = i * (i + 1) / 2;
            for j in 0..=i {
                let d = dists[(i, j)];
                sq_packed[id + j] = d * d;
                sum_sq_all += d * d;
            }
        }
        sum_sq_all /= (n * n) as f64;

        let mut d0 = vec![0.0f64; n];
        for i in 0..n {
            let mut row_sum = 0.0f64;
            for j in 0..n {
                let idx = if i >= j {
                    i * (i + 1) / 2 + j
                } else {
                    j * (j + 1) / 2 + i
                };
                row_sum += sq_packed[idx];
            }
            d0[i] = row_sum / n as f64 - sum_sq_all;
        }

        let mut t_packed = vec![0.0f64; d_size];
        for i in 0..n {
            let id = i * (i + 1) / 2;
            for j in 0..=i {
                let sq_val = sq_packed[if i >= j {
                    i * (i + 1) / 2 + j
                } else {
                    j * (j + 1) / 2 + i
                }];
                t_packed[id + j] = 0.5 * (d0[i] + d0[j] - sq_val);
            }
        }

        let eigen_seed = (sum_sq_all * n as f64) as i32;
        println!("\n=== Metric Matrix Diagnostics ===");
        println!("sum_sq_all = {:.10}", sum_sq_all);
        println!("eigen_seed = {}", eigen_seed);
        println!(
            "d0[0]={:.10}, d0[1]={:.10}, d0[2]={:.10}",
            d0[0], d0[1], d0[2]
        );
        println!(
            "t_packed[46]={:.10} t_packed[47]={:.10}",
            t_packed[46], t_packed[47]
        );

        let result = sci_form::distgeom::power_eigen_solver(3, &mut t_packed, n, eigen_seed);
        match result {
            Some((eigenvalues, eigenvectors)) => {
                for (i, ev) in eigenvalues.iter().enumerate() {
                    println!(
                        "  λ{} = {:.10}  v[0..3] = [{:.10}, {:.10}, {:.10}]",
                        i, ev, eigenvectors[i][0], eigenvectors[i][1], eigenvectors[i][2]
                    );
                }
            }
            None => println!("  Power iteration FAILED!"),
        }
    }

    // Step 2b: Initial coordinates
    let coords_opt = compute_initial_coords_rdkit(&mut rng, &dists, embed_dim);
    let mut coords = match coords_opt {
        Some(c) => c,
        None => {
            println!("EMBEDDING FAILED!");
            return;
        }
    };
    println!("\n=== Initial Coordinates (first 5 atoms) ===");
    for i in 0..n.min(5) {
        println!(
            "  {}: ({:.8}, {:.8}, {:.8})",
            i,
            coords[(i, 0)],
            coords[(i, 1)],
            coords[(i, 2)]
        );
    }

    // Step 3: First minimization (bounds FF)
    let basin = 5.0f32;
    let force_tol = 1e-3f64;
    let mut need_more = 1;
    let mut iters = 0;
    while need_more != 0 {
        need_more = minimize_bfgs_rdkit(
            &mut coords,
            &bounds,
            &chiral_sets,
            400,
            force_tol,
            basin,
            0.1,
            1.0,
        );
        iters += 1;
    }
    println!("\n=== After Bounds FF (stage 1, {} rounds) ===", iters);
    for i in 0..n.min(5) {
        println!(
            "  {}: ({:.8}, {:.8}, {:.8})",
            i,
            coords[(i, 0)],
            coords[(i, 1)],
            coords[(i, 2)]
        );
    }

    // Step 4: Check energy
    let e_total = sci_form::forcefield::bounds_ff::bounds_violation_energy_basin(
        &coords.map(|v| v as f32),
        &bounds,
        basin,
    );
    println!(
        "Bounds energy: {:.6}, energy/atom: {:.6}",
        e_total,
        e_total / n as f32
    );

    // Step 5: Second minimization if 4D
    if use_4d {
        let mut need_more2 = 1;
        let mut iters2 = 0;
        while need_more2 != 0 {
            need_more2 = minimize_bfgs_rdkit(
                &mut coords,
                &bounds,
                &chiral_sets,
                200,
                force_tol,
                basin,
                1.0,
                0.2,
            );
            iters2 += 1;
        }
        println!("\n=== After Bounds FF (stage 2, {} rounds) ===", iters2);
        for i in 0..n.min(5) {
            println!(
                "  {}: ({:.8}, {:.8}, {:.8})",
                i,
                coords[(i, 0)],
                coords[(i, 1)],
                coords[(i, 2)]
            );
        }
    }

    // Step 6: Drop to 3D
    let coords3d = coords.columns(0, 3).into_owned();
    println!("\n=== After Drop to 3D (first 5) ===");
    for i in 0..n.min(5) {
        println!(
            "  {}: ({:.8}, {:.8}, {:.8})",
            i,
            coords3d[(i, 0)],
            coords3d[(i, 1)],
            coords3d[(i, 2)]
        );
    }

    // Step 7: ETKDG 3D FF minimization
    let ff = build_etkdg_3d_ff_with_torsions(&mol, &coords3d, &bounds, &[]);
    let refined = minimize_etkdg_3d_bfgs(&mol, &coords3d, &ff, 300, 1e-3);
    println!("\n=== After ETKDG 3D Minimization (first 5) ===");
    for i in 0..n.min(5) {
        println!(
            "  {}: ({:.8}, {:.8}, {:.8})",
            i,
            refined[(i, 0)],
            refined[(i, 1)],
            refined[(i, 2)]
        );
    }

    // Load RDKit final coordinates for comparison
    let rdkit_coords = load_npy_f64("/tmp/rdkit_coords_fail.npy");
    println!("\n=== RDKit Final Coordinates (first 5) ===");
    for i in 0..n.min(5) {
        println!(
            "  {}: ({:.8}, {:.8}, {:.8})",
            i,
            rdkit_coords[(i, 0)],
            rdkit_coords[(i, 1)],
            rdkit_coords[(i, 2)]
        );
    }

    // Compute RMSD (with alignment)
    let our_f32 = refined.map(|v| v as f32);
    let rdkit_f32 = rdkit_coords.map(|v| v as f32);
    // Simple centroid-aligned RMSD
    let mut our_cx = 0.0f32;
    let mut our_cy = 0.0f32;
    let mut our_cz = 0.0f32;
    let mut rd_cx = 0.0f32;
    let mut rd_cy = 0.0f32;
    let mut rd_cz = 0.0f32;
    for i in 0..n {
        our_cx += our_f32[(i, 0)];
        our_cy += our_f32[(i, 1)];
        our_cz += our_f32[(i, 2)];
        rd_cx += rdkit_f32[(i, 0)];
        rd_cy += rdkit_f32[(i, 1)];
        rd_cz += rdkit_f32[(i, 2)];
    }
    our_cx /= n as f32;
    our_cy /= n as f32;
    our_cz /= n as f32;
    rd_cx /= n as f32;
    rd_cy /= n as f32;
    rd_cz /= n as f32;
    let mut sum_sq = 0.0f32;
    for i in 0..n {
        let dx = (our_f32[(i, 0)] - our_cx) - (rdkit_f32[(i, 0)] - rd_cx);
        let dy = (our_f32[(i, 1)] - our_cy) - (rdkit_f32[(i, 1)] - rd_cy);
        let dz = (our_f32[(i, 2)] - our_cz) - (rdkit_f32[(i, 2)] - rd_cz);
        sum_sq += dx * dx + dy * dy + dz * dz;
    }
    let rmsd_no_rotation = (sum_sq / n as f32).sqrt();
    println!(
        "\nRMSD (centroid-only, no rotation): {:.4}",
        rmsd_no_rotation
    );
}
