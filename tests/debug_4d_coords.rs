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
/// Debug test: Compare our 4D→3D coordinates with RDKit's pure DG output
/// to determine if the divergence happens in the 4D stage or the ETKDG 3D stage.
use nalgebra::DMatrix;

#[test]
fn debug_4d_coords_mol34() {
    // CCCNC(=O)C=O - this is Mol 34 in the test, RMSD 0.748
    let smi = "CCCNC(=O)C=O";
    let mol = sci_form::graph::Molecule::from_smiles(smi).unwrap();
    let n = mol.graph.node_count();
    println!("Molecule: {} atoms", n);

    // Build bounds matrix (same as conformer pipeline)
    let bounds = {
        let raw = sci_form::distgeom::calculate_bounds_matrix_opts(&mol, true);
        let mut b = raw;
        if sci_form::distgeom::triangle_smooth_tol(&mut b, 0.05) {
            b
        } else {
            let raw2 = sci_form::distgeom::calculate_bounds_matrix_opts(&mol, false);
            let mut b2 = raw2.clone();
            if sci_form::distgeom::triangle_smooth_tol(&mut b2, 0.0) {
                b2
            } else {
                let mut b3 = raw2;
                sci_form::distgeom::triangle_smooth_tol(&mut b3, 0.05);
                b3
            }
        }
    };

    let chiral_sets = sci_form::distgeom::identify_chiral_sets(&mol);

    // Print select bounds values for comparison with RDKit
    {
        let heavy: Vec<usize> = (0..n)
            .filter(|&i| mol.graph[petgraph::graph::NodeIndex::new(i)].element != 1)
            .collect();
        println!("Bounds matrix check (heavy atoms):");
        println!("  bounds[0,1]={:.8} (upper C0-C1)", bounds[(0, 1)]);
        println!("  bounds[1,0]={:.8} (lower C0-C1)", bounds[(1, 0)]);
        println!("  bounds[0,7]={:.8} (upper C0-O7)", bounds[(0, 7)]);
        println!("  bounds[7,0]={:.8} (lower C0-O7)", bounds[(7, 0)]);
        println!("  bounds[6,7]={:.8} (upper C6-O7)", bounds[(6, 7)]);
        println!("  bounds[7,6]={:.8} (lower C6-O7)", bounds[(7, 6)]);
        println!("  bounds[4,7]={:.8} (upper C4-O7)", bounds[(4, 7)]);
        println!("  bounds[7,4]={:.8} (lower C4-O7)", bounds[(7, 4)]);
    }

    let mut rng = sci_form::distgeom::MinstdRand::new(42);

    // Step 1: Generate initial coords
    let dists = sci_form::distgeom::pick_rdkit_distances(&mut rng, &bounds);
    let coords4d_opt = sci_form::distgeom::compute_initial_coords_rdkit(&mut rng, &dists, 4);
    let mut coords4d = coords4d_opt.expect("Should get initial coords");

    println!("Initial 4D coords (heavy atoms only):");
    for i in 0..n {
        let elem = mol.graph[petgraph::graph::NodeIndex::new(i)].element;
        if elem != 1 {
            println!(
                "  INIT atom {} (elem={}): ({:.10}, {:.10}, {:.10}, {:.10})",
                i,
                elem,
                coords4d[(i, 0)],
                coords4d[(i, 1)],
                coords4d[(i, 2)],
                coords4d[(i, 3)]
            );
        }
    }

    // Step 2: First minimization (bounds FF)
    {
        let mut need_more = 1;
        let mut iters = 0;
        while need_more != 0 {
            need_more = sci_form::forcefield::bounds_ff::minimize_bfgs_rdkit(
                &mut coords4d,
                &bounds,
                &chiral_sets,
                400,
                1e-3,
                5.0,
                0.1,
                1.0,
            );
            iters += 1;
        }
        println!("\nFirst minimization: {} passes", iters);
    }

    println!("\nPost-4D-min coords (heavy atoms):");
    for i in 0..n {
        let elem = mol.graph[petgraph::graph::NodeIndex::new(i)].element;
        if elem != 1 {
            println!(
                "  POST atom {} (elem={}): ({:.10}, {:.10}, {:.10}, 4d={:.10})",
                i,
                elem,
                coords4d[(i, 0)],
                coords4d[(i, 1)],
                coords4d[(i, 2)],
                coords4d[(i, 3)]
            );
        }
    }

    // Step 7: Drop to 3D
    let coords3d = coords4d.columns(0, 3).into_owned();

    println!("\n3D coords (after 4D drop, before ETKDG):");
    for i in 0..n {
        let elem = mol.graph[petgraph::graph::NodeIndex::new(i)].element;
        if elem != 1 {
            println!(
                "  atom {} (elem={}): ({:.8}, {:.8}, {:.8})",
                i,
                elem,
                coords3d[(i, 0)],
                coords3d[(i, 1)],
                coords3d[(i, 2)]
            );
        }
    }

    // RDKit pure DG reference coords (for CCCNC(=O)C=O with seed=42, no ETKDG)
    let rdkit_dg_heavy: Vec<(f64, f64, f64)> = vec![
        (-2.58678039, 0.57139509, 0.33073414),  // C0
        (-1.34023567, 0.21234681, -0.43757262), // C1
        (-0.41889375, -0.57745987, 0.46375821), // C2
        (0.86305662, 0.09428429, 0.65536763),   // N3
        (2.08045598, -0.65495130, 0.55725872),  // C4
        (1.96244499, -1.88274871, 0.26029674),  // O5
        (3.24934688, 0.09624547, 0.09084923),   // C6
        (3.87601159, 0.90440786, 0.81955589),   // O7
    ];

    // Compare heavy atom pairwise distances
    let heavy_idxs: Vec<usize> = (0..n)
        .filter(|&i| mol.graph[petgraph::graph::NodeIndex::new(i)].element != 1)
        .collect();

    let mut sq_sum = 0.0f64;
    let mut npairs = 0;
    for a in 0..heavy_idxs.len() {
        for b in (a + 1)..heavy_idxs.len() {
            let i = heavy_idxs[a];
            let j = heavy_idxs[b];
            let our_d = ((coords3d[(i, 0)] - coords3d[(j, 0)]).powi(2)
                + (coords3d[(i, 1)] - coords3d[(j, 1)]).powi(2)
                + (coords3d[(i, 2)] - coords3d[(j, 2)]).powi(2))
            .sqrt();
            let rdkit_d = ((rdkit_dg_heavy[a].0 - rdkit_dg_heavy[b].0).powi(2)
                + (rdkit_dg_heavy[a].1 - rdkit_dg_heavy[b].1).powi(2)
                + (rdkit_dg_heavy[a].2 - rdkit_dg_heavy[b].2).powi(2))
            .sqrt();
            let diff = (our_d - rdkit_d).abs();
            if diff > 0.01 {
                println!(
                    "  Heavy pair ({},{}) = atoms ({},{}): our={:.6} rdkit={:.6} diff={:.6}",
                    a, b, i, j, our_d, rdkit_d, diff
                );
            }
            sq_sum += (our_d - rdkit_d).powi(2);
            npairs += 1;
        }
    }
    let rmsd = (sq_sum / npairs as f64).sqrt();
    println!(
        "\nPairwise distance RMSD (heavy atoms, 4D→3D vs RDKit pure DG): {:.6}",
        rmsd
    );

    // Also compute ALL-atom pairwise RMSD in 3D (before ETKDG)
    // Full RDKit pure DG coords (all 17 atoms)
    let rdkit_dg_all: Vec<(f64, f64, f64)> = vec![
        (-2.5867803895, 0.5713950930, 0.3307341426),
        (-1.3402356730, 0.2123468065, -0.4375726167),
        (-0.4188937495, -0.5774598721, 0.4637582107),
        (0.8630566226, 0.0942842875, 0.6553676330),
        (2.0804559804, -0.6549513024, 0.5572587154),
        (1.9624449916, -1.8827487117, 0.2602967378),
        (3.2493468844, 0.0962454746, 0.0908492298),
        (3.8760115928, 0.9044078570, 0.8195558851),
        (-3.4553627869, -0.0528625078, 0.0659866434),
        (-2.3817526881, 0.4734951335, 1.4160589757),
        (-2.8335808486, 1.6329785170, 0.1521877194),
        (-0.8530144992, 1.1343616131, -0.8383510108),
        (-1.6012232886, -0.3974005493, -1.3235707096),
        (-0.2595741287, -1.5997021979, 0.0640145795),
        (-0.8569985747, -0.7209924414, 1.4826308725),
        (0.9100155321, 0.9359816233, 1.2882643157),
        (3.6460850230, -0.1693788227, -0.8877952348),
    ];

    let mut sq_sum_all = 0.0f64;
    let mut npairs_all = 0;
    for i in 0..n {
        for j in (i + 1)..n {
            let our_d = ((coords3d[(i, 0)] - coords3d[(j, 0)]).powi(2)
                + (coords3d[(i, 1)] - coords3d[(j, 1)]).powi(2)
                + (coords3d[(i, 2)] - coords3d[(j, 2)]).powi(2))
            .sqrt();
            let rdkit_d = ((rdkit_dg_all[i].0 - rdkit_dg_all[j].0).powi(2)
                + (rdkit_dg_all[i].1 - rdkit_dg_all[j].1).powi(2)
                + (rdkit_dg_all[i].2 - rdkit_dg_all[j].2).powi(2))
            .sqrt();
            sq_sum_all += (our_d - rdkit_d).powi(2);
            npairs_all += 1;
        }
    }
    let rmsd_all = (sq_sum_all / npairs_all as f64).sqrt();
    println!(
        "Pairwise distance RMSD (all atoms, 4D→3D vs RDKit pure DG): {:.6}",
        rmsd_all
    );
}

#[test]
fn debug_ethane_coords() {
    let mol = sci_form::graph::Molecule::from_smiles("CC").unwrap();
    let n = mol.graph.node_count();
    assert_eq!(n, 8);

    // Build bounds matrix
    let bounds = {
        let raw = sci_form::distgeom::calculate_bounds_matrix_opts(&mol, true);
        let mut b = raw;
        if sci_form::distgeom::triangle_smooth_tol(&mut b, 0.05) {
            b
        } else {
            let raw2 = sci_form::distgeom::calculate_bounds_matrix_opts(&mol, false);
            let mut b2 = raw2.clone();
            if sci_form::distgeom::triangle_smooth_tol(&mut b2, 0.0) {
                b2
            } else {
                let mut b3 = raw2;
                sci_form::distgeom::triangle_smooth_tol(&mut b3, 0.05);
                b3
            }
        }
    };

    // RDKit bounds (upper in [i,j] where i<j, lower in [j,i] where j>i)
    // Format: upper stored in row<col, lower stored in row>col
    // bounds[i,j] with i<j = upper bound, bounds[j,i] with j>i = lower bound
    let rdkit_bounds: Vec<(usize, usize, f64, f64)> = vec![
        // (i,j, upper, lower)
        (0, 1, 1.5240, 1.5040),
        (0, 2, 1.1194007949, 1.0994007949),
        (0, 3, 1.1194007949, 1.0994007949),
        (0, 4, 1.1194007949, 1.0994007949),
        (0, 5, 2.1950665942, 2.1150665942),
        (0, 6, 2.1950665942, 2.1150665942),
        (0, 7, 2.1950665942, 2.1150665942),
        (1, 5, 1.1194007949, 1.0994007949),
        (1, 6, 1.1194007949, 1.0994007949),
        (1, 7, 1.1194007949, 1.0994007949),
    ];

    println!("Bounds comparison (upper=b[i,j], lower=b[j,i]):");
    let mut bounds_match = true;
    for &(i, j, rdkit_upper, rdkit_lower) in &rdkit_bounds {
        let our_upper = bounds[(i, j)];
        let our_lower = bounds[(j, i)];
        let u_diff = (our_upper - rdkit_upper).abs();
        let l_diff = (our_lower - rdkit_lower).abs();
        if u_diff > 1e-6 || l_diff > 1e-6 {
            println!("  MISMATCH ({},{}): upper ours={:.10} rdkit={:.10} diff={:.2e}, lower ours={:.10} rdkit={:.10} diff={:.2e}",
                i, j, our_upper, rdkit_upper, u_diff, our_lower, rdkit_lower, l_diff);
            bounds_match = false;
        }
    }
    if bounds_match {
        println!("  All checked bounds match!");
    }

    // Sample random distances
    let mut rng = sci_form::distgeom::MinstdRand::new(42);
    let dists = sci_form::distgeom::pick_rdkit_distances(&mut rng, &bounds);

    println!("\nRandom distances (select pairs):");
    for i in 0..std::cmp::min(n, 5) {
        for j in 0..i {
            println!("  dist[{},{}] = {:.16}", i, j, dists[(i, j)]);
        }
    }

    // Initial 4D coords
    let coords4d_opt = sci_form::distgeom::compute_initial_coords_rdkit(&mut rng, &dists, 4);
    let mut coords4d = coords4d_opt.expect("Should get initial coords");

    println!("\nInitial 4D coords:");
    for i in 0..n {
        println!(
            "  atom {}: ({:.12}, {:.12}, {:.12}, {:.12})",
            i,
            coords4d[(i, 0)],
            coords4d[(i, 1)],
            coords4d[(i, 2)],
            coords4d[(i, 3)]
        );
    }

    // First minimization (4D)
    let chiral_sets = sci_form::distgeom::identify_chiral_sets(&mol);
    {
        let mut need_more = 1;
        while need_more != 0 {
            need_more = sci_form::forcefield::bounds_ff::minimize_bfgs_rdkit(
                &mut coords4d,
                &bounds,
                &chiral_sets,
                400,
                1e-3,
                5.0,
                0.1,
                1.0,
            );
        }
    }

    println!("\nPost-4D-BFGS coords:");
    for i in 0..n {
        println!(
            "  atom {}: ({:.12}, {:.12}, {:.12}, {:.12})",
            i,
            coords4d[(i, 0)],
            coords4d[(i, 1)],
            coords4d[(i, 2)],
            coords4d[(i, 3)]
        );
    }

    // Drop to 3D
    let coords3d_f64 = coords4d.columns(0, 3).into_owned();

    println!("\n3D coords (pre-ETKDG):");
    for i in 0..n {
        println!(
            "  atom {}: ({:.12}, {:.12}, {:.12})",
            i,
            coords3d_f64[(i, 0)],
            coords3d_f64[(i, 1)],
            coords3d_f64[(i, 2)]
        );
    }

    // Final coords via full pipeline
    let final_coords = sci_form::conformer::generate_3d_conformer(&mol, 42).unwrap();

    // RDKit reference
    let rdkit_ref: Vec<[f64; 3]> = vec![
        [-0.7598767677, -0.0045137979, -0.0513340444],
        [0.7378738999, -0.0346971899, 0.0053987715],
        [-1.1288615848, -0.1839044855, 0.9818889357],
        [-1.0729072847, 1.0194893703, -0.3235095125],
        [-1.2009112670, -0.7879231094, -0.6786082171],
        [1.0531716678, 0.8976810181, 0.5402779301],
        [1.2204626537, 0.0240941223, -0.9800832323],
        [1.1510486828, -0.9302259280, 0.5059693690],
    ];

    println!("\nFinal coords comparison:");
    for i in 0..n {
        let dx = final_coords[(i, 0)] as f64 - rdkit_ref[i][0];
        let dy = final_coords[(i, 1)] as f64 - rdkit_ref[i][1];
        let dz = final_coords[(i, 2)] as f64 - rdkit_ref[i][2];
        let d = (dx * dx + dy * dy + dz * dz).sqrt();
        println!(
            "  atom {}: ours=({:.10}, {:.10}, {:.10}) rdkit=({:.10}, {:.10}, {:.10}) dist={:.8}",
            i,
            final_coords[(i, 0)],
            final_coords[(i, 1)],
            final_coords[(i, 2)],
            rdkit_ref[i][0],
            rdkit_ref[i][1],
            rdkit_ref[i][2],
            d
        );
    }

    // Pairwise RMSD
    let mut sq_sum = 0.0f64;
    let mut npairs = 0;
    for i in 0..n {
        for j in (i + 1)..n {
            let our_d = (((final_coords[(i, 0)] - final_coords[(j, 0)]) as f64).powi(2)
                + ((final_coords[(i, 1)] - final_coords[(j, 1)]) as f64).powi(2)
                + ((final_coords[(i, 2)] - final_coords[(j, 2)]) as f64).powi(2))
            .sqrt();
            let rdkit_d = ((rdkit_ref[i][0] - rdkit_ref[j][0]).powi(2)
                + (rdkit_ref[i][1] - rdkit_ref[j][1]).powi(2)
                + (rdkit_ref[i][2] - rdkit_ref[j][2]).powi(2))
            .sqrt();
            sq_sum += (our_d - rdkit_d).powi(2);
            npairs += 1;
        }
    }
    println!("\nPairwise RMSD: {:.8}", (sq_sum / npairs as f64).sqrt());
}

#[test]
fn debug_rng_values() {
    let mut rng = sci_form::distgeom::MinstdRand::new(42);
    println!("First 5 RNG values:");
    for i in 0..5 {
        let v = rng.next_double();
        println!("  {}: {:.18}", i, v);
    }
}
