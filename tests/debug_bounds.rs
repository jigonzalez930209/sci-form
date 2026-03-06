use rand::SeedableRng;
use sci_form::distgeom;
use sci_form::graph::Molecule;

#[test]
fn test_debug_bounds() {
    let mol = Molecule::from_smiles("CCCCC").unwrap();
    let bounds = distgeom::calculate_bounds_matrix(&mol);
    let smoothed = distgeom::smooth_bounds_matrix(bounds.clone());

    println!(
        "Initial Bounds 0-4: U={}, L={}",
        bounds[(0, 4)],
        bounds[(4, 0)]
    );
    println!(
        "Smoothed Bounds 0-4: U={}, L={}",
        smoothed[(0, 4)],
        smoothed[(4, 0)]
    );

    let mut rng = rand::rngs::StdRng::seed_from_u64(42);
    let dists = distgeom::pick_random_distances(&mut rng, &smoothed);
    println!("Picked Dist 0-4: {}", dists[(0, 4)]);

    let metric = distgeom::compute_metric_matrix(&dists);
    let coords3d = distgeom::generate_3d_coordinates(&mut rng, &metric);
    let dr = (coords3d[(0, 0)] - coords3d[(4, 0)]).powi(2)
        + (coords3d[(0, 1)] - coords3d[(4, 1)]).powi(2)
        + (coords3d[(0, 2)] - coords3d[(4, 2)]).powi(2);
    println!("Coords dist 0-4: {}", dr.sqrt());
}
