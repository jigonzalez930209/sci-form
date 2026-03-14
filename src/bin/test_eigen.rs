use sci_form::distgeom::{
    calculate_bounds_matrix_opts, identify_chiral_sets, pick_rdkit_distances, power_eigen_solver,
    triangle_smooth_tol, MinstdRand,
};
use sci_form::graph::Molecule;

fn main() {
    let smiles = "C#CCOC(C)CC1CC2C3CCC(C)C(O)(C3)C2O1";
    let mol = Molecule::from_smiles(smiles).unwrap();
    let n = mol.graph.node_count();
    println!("Molecule: {}, N={}", smiles, n);

    let raw = calculate_bounds_matrix_opts(&mol, true);
    let mut bounds = raw;
    triangle_smooth_tol(&mut bounds, 0.0);

    let chiral_sets = identify_chiral_sets(&mol);
    let embed_dim = if !chiral_sets.is_empty() { 4 } else { 3 };
    println!("embed_dim = {}", embed_dim);

    let mut rng = MinstdRand::new(42);
    let dists = pick_rdkit_distances(&mut rng, &bounds);

    // Build metric matrix (matching compute_initial_coords_rdkit exactly)
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
    #[allow(clippy::needless_range_loop)]
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
            let sq_val = sq_packed[i * (i + 1) / 2 + j];
            t_packed[id + j] = 0.5 * (d0[i] + d0[j] - sq_val);
        }
    }

    let eigen_seed = (sum_sq_all * n as f64) as i32;
    println!("sum_sq_all = {:.15}", sum_sq_all);
    println!("eigen_seed = {}", eigen_seed);

    let result = power_eigen_solver(embed_dim, &mut t_packed, n, eigen_seed);
    match result {
        Some((eigenvalues, eigenvectors)) => {
            println!("\nEigenvalues:");
            for (i, ev) in eigenvalues.iter().enumerate() {
                println!("  λ{} = {:.15}", i, ev);
            }
            println!("\nFirst eigenvector (first 5):");
            #[allow(clippy::needless_range_loop)]
            for j in 0..5.min(n) {
                println!("  v0[{}] = {:.15}", j, eigenvectors[0][j]);
            }
            println!("\nSecond eigenvector (first 5):");
            #[allow(clippy::needless_range_loop)]
            for j in 0..5.min(n) {
                println!("  v1[{}] = {:.15}", j, eigenvectors[1][j]);
            }
        }
        None => {
            println!("Power iteration FAILED!");
        }
    }
}
