use sci_form::distgeom::bounds::{calculate_bounds_matrix_opts, triangle_smooth_tol};
use sci_form::graph::Molecule;

#[test]
fn test_compare_bounds_matrix() {
    let smiles = "C#CC(C)(COCCC1CCCNCC1)OCC(C)=O";
    let mol = Molecule::from_smiles(smiles).unwrap();
    let n = mol.graph.node_count();

    println!("Molecule: {}", smiles);
    println!("Atoms: {}", n);

    let raw = calculate_bounds_matrix_opts(&mol, true);
    let mut bounds = raw;
    let smooth_ok = triangle_smooth_tol(&mut bounds, 0.0);
    println!(
        "Triangle smooth: {}",
        if smooth_ok { "OK" } else { "FAILED" }
    );

    let rdkit_bounds_path = "/tmp/rdkit_bounds.npy";
    let data = std::fs::read(rdkit_bounds_path).expect("Run Python script first");

    let major = data[6] as u16;
    let header_len = if major >= 2 {
        u32::from_le_bytes([data[8], data[9], data[10], data[11]]) as usize
    } else {
        u16::from_le_bytes([data[8], data[9]]) as usize
    };
    let data_start = if major >= 2 {
        12 + header_len
    } else {
        10 + header_len
    };
    let float_data = &data[data_start..];

    let rdkit_bounds: Vec<f64> = float_data
        .chunks_exact(8)
        .map(|chunk| f64::from_le_bytes(chunk.try_into().unwrap()))
        .collect();

    assert_eq!(rdkit_bounds.len(), n * n, "Size mismatch");

    let mut max_ub_diff = 0.0f64;
    let mut max_lb_diff = 0.0f64;
    let mut ub_diffs = 0;
    let mut lb_diffs = 0;
    let mut total_pairs = 0;

    for i in 0..n {
        for j in (i + 1)..n {
            let our_ub = bounds[(i, j)];
            let our_lb = bounds[(j, i)];
            let rdkit_ub = rdkit_bounds[i * n + j];
            let rdkit_lb = rdkit_bounds[j * n + i];

            let ub_diff = (our_ub - rdkit_ub).abs();
            let lb_diff = (our_lb - rdkit_lb).abs();

            if ub_diff > 1e-10 {
                ub_diffs += 1;
            }
            if lb_diff > 1e-10 {
                lb_diffs += 1;
            }
            if ub_diff > max_ub_diff {
                max_ub_diff = ub_diff;
            }
            if lb_diff > max_lb_diff {
                max_lb_diff = lb_diff;
            }
            total_pairs += 1;

            if ub_diff > 0.001 || lb_diff > 0.001 {
                println!("  DIFF ({},{}): our_ub={:.10}, rdkit_ub={:.10}, d={:.2e}; our_lb={:.10}, rdkit_lb={:.10}, d={:.2e}",
                    i, j, our_ub, rdkit_ub, ub_diff, our_lb, rdkit_lb, lb_diff);
            }
        }
    }

    println!("\n=== Bounds matrix comparison ===");
    println!("Total pairs: {}", total_pairs);
    println!("UB diffs (>1e-10): {}", ub_diffs);
    println!("LB diffs (>1e-10): {}", lb_diffs);
    println!("Max UB diff: {:.2e}", max_ub_diff);
    println!("Max LB diff: {:.2e}", max_lb_diff);

    let mut our_ub_sum = 0.0f64;
    let mut our_lb_sum = 0.0f64;
    for i in 0..n {
        for j in (i + 1)..n {
            our_ub_sum += bounds[(i, j)];
            our_lb_sum += bounds[(j, i)];
        }
    }
    println!("\nChecksums:");
    println!("  Our UB sum:   {:.10}", our_ub_sum);
    println!("  RDKit UB sum: 7872.4329018951");
    println!("  Our LB sum:   {:.10}", our_lb_sum);
    println!("  RDKit LB sum: 2654.3162734026");

    println!("\n=== First row of our bounds matrix ===");
    for j in 1..std::cmp::min(n, 20) {
        println!(
            "  bm[0,{:2}] = {:.10} (UB), bm[{:2},0] = {:.10} (LB)",
            j,
            bounds[(0, j)],
            j,
            bounds[(j, 0)]
        );
    }
}
