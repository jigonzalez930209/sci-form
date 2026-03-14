#![allow(unused_imports, unused_variables, dead_code, clippy::unnecessary_cast, clippy::needless_range_loop, clippy::manual_repeat_n, clippy::manual_str_repeat, clippy::manual_is_multiple_of, clippy::redundant_field_names, clippy::useless_vec, clippy::single_range_in_vec_init)]
use sci_form::distgeom;
use sci_form::distgeom::{
    compute_initial_coords_rdkit, identify_chiral_sets, identify_tetrahedral_centers,
    pick_rdkit_distances, MinstdRand,
};
use sci_form::forcefield::bounds_ff::minimize_bfgs_rdkit;
use sci_form::forcefield::etkdg_3d::{build_etkdg_3d_ff_with_torsions, minimize_etkdg_3d_bfgs};
use sci_form::graph::Molecule;

#[test]
fn test_debug_bounds() {
    let smiles = "CC(C#N)(C#N)OC";
    let mol = Molecule::from_smiles(smiles).unwrap();
    let n = mol.graph.node_count();

    let bounds = {
        let raw = distgeom::calculate_bounds_matrix_opts(&mol, true);
        let mut b = raw;
        if distgeom::triangle_smooth_tol(&mut b, 0.05) {
            b
        } else {
            let raw2 = distgeom::calculate_bounds_matrix_opts(&mol, false);
            let mut b2 = raw2.clone();
            if distgeom::triangle_smooth_tol(&mut b2, 0.0) {
                b2
            } else {
                let mut b3 = raw2;
                distgeom::triangle_smooth_tol(&mut b3, 0.05);
                b3
            }
        }
    };
    let chiral_sets = identify_chiral_sets(&mol);

    println!("SMILES: {}", smiles);
    println!("Atoms: {}", n);

    let mut rng = MinstdRand::new(42);
    let dists = pick_rdkit_distances(&mut rng, &bounds);
    let mut coords4d = compute_initial_coords_rdkit(&mut rng, &dists, 4).expect("Embedding failed");

    println!("\n=== INITIAL 4D COORDS ===");
    for i in 0..n {
        println!(
            "  {} {:12.6} {:12.6} {:12.6} {:12.6}",
            i,
            coords4d[(i, 0)],
            coords4d[(i, 1)],
            coords4d[(i, 2)],
            coords4d[(i, 3)]
        );
    }

    // 1st BFGS: chiral_w=1.0, 4d_w=0.1
    {
        let mut nm = 1;
        let mut pass = 0;
        while nm != 0 {
            nm = minimize_bfgs_rdkit(
                &mut coords4d,
                &bounds,
                &chiral_sets,
                400,
                1e-4,
                5.0,
                0.1,
                1.0,
            );
            pass += 1;
        }
        println!("\n=== AFTER 1st BFGS ({} passes) ===", pass);
        for i in 0..n {
            println!(
                "  {} {:12.6} {:12.6} {:12.6} {:12.6}",
                i,
                coords4d[(i, 0)],
                coords4d[(i, 1)],
                coords4d[(i, 2)],
                coords4d[(i, 3)]
            );
        }
    }

    // 2nd BFGS: chiral_w=0.2, 4d_w=1.0
    {
        let mut nm = 1;
        let mut pass = 0;
        while nm != 0 {
            nm = minimize_bfgs_rdkit(
                &mut coords4d,
                &bounds,
                &chiral_sets,
                200,
                1e-4,
                5.0,
                1.0,
                0.2,
            );
            pass += 1;
        }
        println!("\n=== AFTER 2nd BFGS ({} passes) ===", pass);
        for i in 0..n {
            println!(
                "  {} {:12.6} {:12.6} {:12.6} {:12.6}",
                i,
                coords4d[(i, 0)],
                coords4d[(i, 1)],
                coords4d[(i, 2)],
                coords4d[(i, 3)]
            );
        }
    }

    // Drop to 3D
    let coords3d = coords4d.columns(0, 3).into_owned();
    println!("\n=== 3D PROJECTION ===");
    for i in 0..n {
        println!(
            "  {} {:12.6} {:12.6} {:12.6}",
            i,
            coords3d[(i, 0)],
            coords3d[(i, 1)],
            coords3d[(i, 2)]
        );
    }

    // ETKDG 3D minimization (no CSD torsions for this test)
    let ff = build_etkdg_3d_ff_with_torsions(&mol, &coords3d, &bounds, &[]);
    let refined = minimize_etkdg_3d_bfgs(&mol, &coords3d, &ff, 300, 1e-4);
    println!("\n=== AFTER 3D ETKDG MINIMIZATION ===");
    for i in 0..n {
        println!(
            "  {} {:12.6} {:12.6} {:12.6}",
            i,
            refined[(i, 0)],
            refined[(i, 1)],
            refined[(i, 2)]
        );
    }
}
