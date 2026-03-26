//! Debug: exercise SSSR ring perception on assorted ring systems.
//!
//! Category: debug

use sci_form::graph::Molecule;

fn main() {
    let molecules = [
        "C1C2CC1C2",            // bicyclo[1.1.1]pentane
        "CC12CC(C1)C2",         // norbornane-like
        "C1CC2CCC12",           // bicyclo[2.2.0]hexane
        "C1C2C3CC2C13",         // tricyclic
        "CC1C2CC1C2",           // methylbicyclopentane
        "CC1(C)C2CC1C2",        // dimethyl
        "OC1C2CC1C2",           // OH-bicyclopentane
        "C1CC2CC1C2",           // bicyclo[2.1.1]hexane
        "C12CC3CC(C1)CC(C2)C3", // adamantane
    ];

    for smi in &molecules {
        let mol = match Molecule::from_smiles(smi) {
            Ok(m) => m,
            Err(_) => {
                println!("{:30} PARSE_FAIL", smi);
                continue;
            }
        };
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
        let tet_centers = sci_form::distgeom::identify_tetrahedral_centers(&mol);

        let use_4d = !chiral_sets.is_empty();
        let embed_dim = if use_4d { 4 } else { 3 };

        let mut rng = sci_form::distgeom::MinstdRand::new(42);

        let mut step_fails = [0u32; 7]; // 0:eigen 1:energy 2:tet 3:chiral 4:planar 5:dblbond 6:smooth_fail
        let max_iter = 10 * n;
        let mut ok = false;

        for _iter in 0..max_iter {
            let dists = sci_form::distgeom::pick_rdkit_distances(&mut rng, &bounds);
            let coords_opt =
                sci_form::distgeom::compute_initial_coords_rdkit(&mut rng, &dists, embed_dim);
            let mut coords = match coords_opt {
                Some(c) => c,
                None => {
                    step_fails[0] += 1;
                    continue;
                }
            };

            {
                let mut need = 1;
                let mut restarts = 0;
                while need != 0 && restarts < 5 {
                    need = sci_form::forcefield::bounds_ff::minimize_bfgs_rdkit(
                        &mut coords,
                        &bounds,
                        &chiral_sets,
                        400,
                        1e-3,
                        5.0,
                        0.1,
                        1.0,
                    );
                    restarts += 1;
                }
            }

            let c32 = coords.map(|v| v as f32);
            let e =
                sci_form::forcefield::bounds_ff::bounds_violation_energy_basin(&c32, &bounds, 5.0)
                    + if !chiral_sets.is_empty() {
                        sci_form::forcefield::bounds_ff::chiral_violation_energy(&c32, &chiral_sets)
                    } else {
                        0.0
                    };
            if e / n as f32 >= 0.05 {
                step_fails[1] += 1;
                continue;
            }

            if !sci_form::distgeom::check_tetrahedral_centers(&coords, &tet_centers) {
                step_fails[2] += 1;
                continue;
            }

            if !chiral_sets.is_empty()
                && !sci_form::distgeom::check_chiral_centers(&coords, &chiral_sets)
            {
                step_fails[3] += 1;
                continue;
            }

            if use_4d {
                let mut need = 1;
                let mut restarts = 0;
                while need != 0 && restarts < 5 {
                    need = sci_form::forcefield::bounds_ff::minimize_bfgs_rdkit(
                        &mut coords,
                        &bounds,
                        &chiral_sets,
                        200,
                        1e-3,
                        5.0,
                        1.0,
                        0.2,
                    );
                    restarts += 1;
                }
            }

            let coords3d = coords.columns(0, 3).into_owned();
            let ff = sci_form::forcefield::etkdg_3d::build_etkdg_3d_ff_with_torsions(
                &mol,
                &coords3d,
                &bounds,
                &[],
            );
            let refined = sci_form::forcefield::etkdg_3d::minimize_etkdg_3d_bfgs(
                &mol, &coords3d, &ff, 300, 1e-3,
            );
            let ref_f32 = refined.map(|v| v as f32);

            {
                let n_improper_atoms = ff.inversion_contribs.len() / 3;
                let planarity_energy =
                    sci_form::forcefield::etkdg_3d::planarity_check_energy(&ref_f32, &ff);
                if planarity_energy > n_improper_atoms as f32 * 0.7 {
                    step_fails[4] += 1;
                    continue;
                }
            }

            if !sci_form::distgeom::check_double_bond_geometry(&mol, &refined) {
                step_fails[5] += 1;
                continue;
            }

            println!(
                "{:30} OK iter={:3}  eig={} ene={} tet={} chi={} pla={} dbl={}",
                smi,
                _iter,
                step_fails[0],
                step_fails[1],
                step_fails[2],
                step_fails[3],
                step_fails[4],
                step_fails[5]
            );
            ok = true;
            break;
        }

        if !ok {
            println!(
                "{:30} FAIL iters={}  eig={} ene={} tet={} chi={} pla={} dbl={}",
                smi,
                max_iter,
                step_fails[0],
                step_fails[1],
                step_fails[2],
                step_fails[3],
                step_fails[4],
                step_fails[5]
            );
        }
    }
}
