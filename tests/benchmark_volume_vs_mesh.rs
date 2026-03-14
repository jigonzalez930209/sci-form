//! Benchmark: raw volume transfer vs marching-cubes isosurface mesh.
//!
//! Compares memory footprint and triangle count for the two rendering paths.
//!
//! Run with:  cargo test --release --test benchmark_volume_vs_mesh -- --nocapture

use sci_form::eht::{basis::build_basis, evaluate_orbital_on_grid, marching_cubes, solve_eht};
use std::time::Instant;

#[test]
#[allow(clippy::type_complexity)]
fn volume_vs_mesh_comparison() {
    let molecules: Vec<(&str, Vec<u8>, Vec<[f64; 3]>)> = vec![
        (
            "H2",
            vec![1, 1],
            vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]],
        ),
        (
            "H2O",
            vec![8, 1, 1],
            vec![
                [0.0, 0.0, 0.117],
                [-0.757, 0.0, -0.469],
                [0.757, 0.0, -0.469],
            ],
        ),
        (
            "CH4",
            vec![6, 1, 1, 1, 1],
            vec![
                [0.0, 0.0, 0.0],
                [0.629, 0.629, 0.629],
                [-0.629, -0.629, 0.629],
                [-0.629, 0.629, -0.629],
                [0.629, -0.629, -0.629],
            ],
        ),
    ];

    let spacings = [0.3, 0.2, 0.15];

    println!();
    println!(
        "{:<6} {:<8} {:>10} {:>12} {:>8} {:>12} {:>10} {:>10}",
        "Mol", "Spacing", "GridPts", "Vol(bytes)", "Tris", "Mesh(bytes)", "VolTime", "MeshTime"
    );
    println!("{}", "-".repeat(86));

    for (name, elems, positions) in &molecules {
        let result = solve_eht(elems, positions, None).unwrap();
        let basis = build_basis(elems, positions);

        for &spacing in &spacings {
            let t0 = Instant::now();
            let grid = evaluate_orbital_on_grid(
                &basis,
                &result.coefficients,
                0,
                positions,
                spacing,
                3.0,
            );
            let vol_time = t0.elapsed();
            let vol_points = grid.num_points();
            let vol_bytes = grid.values.len() * std::mem::size_of::<f32>();

            let t1 = Instant::now();
            let mesh = marching_cubes(&grid, 0.02);
            let mesh_time = t1.elapsed();
            let mesh_bytes = mesh.vertices.len() * std::mem::size_of::<f32>()
                + mesh.normals.len() * std::mem::size_of::<f32>()
                + mesh.indices.len() * std::mem::size_of::<u32>();

            println!(
                "{:<6} {:<8.2} {:>10} {:>12} {:>8} {:>12} {:>9.2}ms {:>9.2}ms",
                name,
                spacing,
                vol_points,
                format_bytes(vol_bytes),
                mesh.num_triangles,
                format_bytes(mesh_bytes),
                vol_time.as_secs_f64() * 1000.0,
                mesh_time.as_secs_f64() * 1000.0,
            );

            // Volume should always have a non-zero footprint
            assert!(vol_bytes > 0);
            // Mesh is typically much smaller than the raw volume
            if mesh.num_triangles > 0 {
                assert!(
                    mesh_bytes < vol_bytes,
                    "{} @ {}: mesh {} >= volume {}",
                    name,
                    spacing,
                    mesh_bytes,
                    vol_bytes
                );
            }
        }
    }
}

fn format_bytes(b: usize) -> String {
    if b < 1024 {
        format!("{} B", b)
    } else if b < 1024 * 1024 {
        format!("{:.1} KB", b as f64 / 1024.0)
    } else {
        format!("{:.1} MB", b as f64 / (1024.0 * 1024.0))
    }
}
