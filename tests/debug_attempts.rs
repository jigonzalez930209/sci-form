/// Debug test: Check which attempts fail for mol 54 and mol 94

#[test]
fn debug_attempt_tracking() {
    let test_cases = vec![
        ("CCCc1cocn1", 42u64),      // mol 54: RDKit has 1 initial_coords failure
        ("CC1(CC1)CCOC", 42u64),    // mol 94: RDKit has 2 initial_coords failures
        ("CC(C#N)(C#N)OC", 42u64),  // mol 27: RDKit has 0 failures
    ];
    
    for (smiles, seed) in test_cases {
        let mol = sci_form::graph::Molecule::from_smiles(smiles).unwrap();
        let n = mol.graph.node_count();
        
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
        
        let mut rng = sci_form::distgeom::MinstdRand::new(seed as u32);
        
        println!("\n=== {} (n={}) ===", smiles, n);
        
        let max_iter = 10;
        for iter in 0..max_iter {
            // Step 1: Pick random distances
            let dists = sci_form::distgeom::pick_rdkit_distances(&mut rng, &bounds);
            
            // Step 2: Compute initial coords
            let coords4d_opt = sci_form::distgeom::compute_initial_coords_rdkit(&mut rng, &dists, 4);
            
            match coords4d_opt {
                Some(mut coords4d) => {
                    // Step 3: First minimization
                    let mut need_more = 1;
                    let mut bfgs_passes = 0;
                    while need_more != 0 {
                        need_more = sci_form::forcefield::bounds_ff::minimize_bfgs_rdkit(
                            &mut coords4d, &bounds, &chiral_sets,
                            400, 1e-3, 5.0, 0.1, 1.0,
                        );
                        bfgs_passes += 1;
                    }
                    
                    // Step 4: Energy check
                    let coords4d_f32 = coords4d.map(|v| v as f32);
                    let coords4d_f64 = coords4d.map(|v| v as f64);
                    let mut energy = sci_form::forcefield::bounds_ff::bounds_violation_energy_basin(&coords4d_f32, &bounds, 5.0);
                    if !chiral_sets.is_empty() {
                        energy += sci_form::forcefield::bounds_ff::chiral_violation_energy(&coords4d_f32, &chiral_sets);
                    }
                    for i in 0..n {
                        if coords4d.ncols() == 4 {
                            let x4 = coords4d_f32[(i, 3)];
                            energy += 0.1 * x4 * x4;
                        }
                    }
                    let e_per_atom = energy / n as f32;
                    
                    if e_per_atom >= 0.05 {
                        println!("  Iter {}: initial_coords OK, but E/atom={:.6} >= 0.05 → reject (first_min fail)", iter, e_per_atom);
                        continue;
                    }
                    
                    // Check tetrahedral centers
                    let tet_centers = sci_form::distgeom::identify_tetrahedral_centers(&mol);
                    if !sci_form::distgeom::check_tetrahedral_centers(&coords4d_f64, &tet_centers) {
                        println!("  Iter {}: tetrahedral check failed", iter);
                        continue;
                    }
                    
                    if !chiral_sets.is_empty() && !sci_form::distgeom::check_chiral_centers(&coords4d_f64, &chiral_sets) {
                        println!("  Iter {}: chiral check failed", iter);
                        continue;
                    }
                    
                    println!("  Iter {}: SUCCESS (E/atom={:.6}, {} BFGS passes)", iter, e_per_atom, bfgs_passes);
                    break;
                }
                None => {
                    println!("  Iter {}: initial_coords FAILED (returned None)", iter);
                }
            }
        }
    }
}
