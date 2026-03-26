//! Debug: exercise the full distance geometry + bounds FF pipeline and
//! dump intermediate embedding attempts.
//!
//! Category: debug

use sci_form::distgeom::*;
use sci_form::forcefield::bounds_ff::*;
use sci_form::graph::Molecule;

fn main() {
    let smiles = "C#CCOC(C)CC1CC2C3CCC(C)C(O)(C3)C2O1";
    let mol = Molecule::from_smiles(smiles).unwrap();
    let n = mol.graph.node_count();

    let bounds = {
        let raw = calculate_bounds_matrix_opts(&mol, true);
        let mut b = raw;
        triangle_smooth_tol(&mut b, 0.0);
        b
    };
    let chiral_sets = identify_chiral_sets(&mol);
    let tet_centers = identify_tetrahedral_centers(&mol);
    let use_4d = !chiral_sets.is_empty();
    let embed_dim = if use_4d { 4 } else { 3 };
    let max_iterations = 10 * n;
    let mut rng = MinstdRand::new(42);

    println!(
        "chiral_sets: {}, use_4d: {}, embed_dim: {}",
        chiral_sets.len(),
        use_4d,
        embed_dim
    );

    for iter in 0..max_iterations {
        let rng_state_before = rng.get_state();
        let dists = pick_rdkit_distances(&mut rng, &bounds);
        let rng_after_dists = rng.get_state();
        let coords_opt = compute_initial_coords_rdkit(&mut rng, &dists, embed_dim);
        let rng_after_coords = rng.get_state();

        let mut coords = match coords_opt {
            Some(c) => c,
            None => {
                println!(
                    "attempt {}: embedding FAILED (rng: {} → {} → {})",
                    iter, rng_state_before, rng_after_dists, rng_after_coords
                );
                continue;
            }
        };

        // Bounds FF minimization
        {
            let mut need_more = 1;
            while need_more != 0 {
                need_more = minimize_bfgs_rdkit(
                    &mut coords,
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

        let coords_f32 = coords.map(|v| v as f32);
        let energy_per_atom = {
            let mut e = 0.0f32;
            for i in 1..n {
                for j in 0..i {
                    let ub = bounds[(j, i)] as f32;
                    let lb = bounds[(i, j)] as f32;
                    if ub - lb > 5.0 {
                        continue;
                    }
                    let mut d2 = 0.0f32;
                    for d in 0..coords_f32.ncols().min(3) {
                        let diff = coords_f32[(i, d)] - coords_f32[(j, d)];
                        d2 += diff * diff;
                    }
                    let ub2 = ub * ub;
                    let lb2 = lb * lb;
                    let val = if d2 > ub2 {
                        d2 / ub2 - 1.0
                    } else if d2 < lb2 {
                        2.0 * lb2 / (lb2 + d2) - 1.0
                    } else {
                        0.0
                    };
                    if val > 0.0 {
                        e += val * val;
                    }
                }
            }
            e / n as f32
        };

        if energy_per_atom >= 0.05 {
            println!("attempt {}: energy FAILED ({:.4})", iter, energy_per_atom);
            continue;
        }

        if !check_tetrahedral_centers(&coords, &tet_centers) {
            println!("attempt {}: tet check FAILED", iter);
            continue;
        }

        if !chiral_sets.is_empty() && !check_chiral_centers(&coords, &chiral_sets) {
            println!("attempt {}: chiral check FAILED", iter);
            continue;
        }

        println!(
            "attempt {}: PASSED all checks (rng: {} → {} → {})",
            iter, rng_state_before, rng_after_dists, rng_after_coords
        );
        println!(
            "  coords[0] = ({:.8}, {:.8}, {:.8})",
            coords[(0, 0)],
            coords[(0, 1)],
            coords[(0, 2)]
        );
        break;
    }
}
