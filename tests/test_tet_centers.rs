#![allow(unused_imports, unused_variables, dead_code, clippy::unnecessary_cast, clippy::needless_range_loop, clippy::manual_repeat_n, clippy::manual_str_repeat, clippy::manual_is_multiple_of, clippy::redundant_field_names, clippy::useless_vec, clippy::single_range_in_vec_init)]
#[test]
fn test_tetrahedral_centers() {
    use sci_form::distgeom::validation::{find_sssr_pub, identify_tetrahedral_centers};
    use sci_form::graph::Molecule;

    let smiles = ["CC1(O)CN2CC12", "CC1N2CC1(O)C2", "CC1(O)C2OCC12C"];
    for smi in &smiles {
        let mol = Molecule::from_smiles(smi).unwrap();
        let rings = find_sssr_pub(&mol);
        let centers = identify_tetrahedral_centers(&mol);
        println!(
            "{}: rings={:?}, {} tetrahedral centers",
            smi,
            rings,
            centers.len()
        );
    }
}

#[test]
fn test_bounds_comparison() {
    use sci_form::distgeom::bounds::calculate_bounds_matrix_opts;
    use sci_form::distgeom::bounds::triangle_smooth_tol;
    use sci_form::graph::Molecule;

    let smi = "CC1(O)CN2CC12";
    let mol = Molecule::from_smiles(smi).unwrap();
    let n = mol.graph.node_count();

    let bounds = {
        let raw = calculate_bounds_matrix_opts(&mol, true);
        let mut b = raw;
        if triangle_smooth_tol(&mut b, 0.05) {
            b
        } else {
            let raw2 = calculate_bounds_matrix_opts(&mol, false);
            let mut b2 = raw2.clone();
            if triangle_smooth_tol(&mut b2, 0.0) {
                b2
            } else {
                let mut b3 = raw2;
                triangle_smooth_tol(&mut b3, 0.05);
                b3
            }
        }
    };

    println!("{}: bounds matrix ({}x{})", smi, n, n);
    println!("Heavy atom pairs:");
    for i in 0..7 {
        for j in (i + 1)..7 {
            // In our code: upper triangle = upper bound, lower triangle = lower bound
            let ub = bounds[(i, j)];
            let lb = bounds[(j, i)];
            println!(
                "  ({},{}) lb={:.4}  ub={:.4}  range={:.4}",
                i,
                j,
                lb,
                ub,
                ub - lb
            );
        }
    }
}
