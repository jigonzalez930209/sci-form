use nalgebra::DMatrix;
use sci_form::forcefield::{
    energy::calculate_total_energy, gradients::compute_analytical_gradient, FFParams,
};
use sci_form::graph::{Atom, Molecule};

fn main() {
    let mut mol = Molecule::new("test");
    // Central atom 0, neighbors 1, 2, 3
    let a0 = mol.add_atom(Atom::new(6, 0.0, 0.0, 0.0));
    let a1 = mol.add_atom(Atom::new(6, 1.0, 0.0, 0.0));
    let a2 = mol.add_atom(Atom::new(6, 0.0, 1.0, 0.0));
    let a3 = mol.add_atom(Atom::new(6, 0.02, 0.02, 0.05));

    // Add bonds so they are neighbors
    mol.add_bond(
        a0,
        a1,
        sci_form::graph::Bond {
            order: sci_form::graph::BondOrder::Single,
            stereo: sci_form::graph::BondStereo::None,
        },
    );
    mol.add_bond(
        a0,
        a2,
        sci_form::graph::Bond {
            order: sci_form::graph::BondOrder::Single,
            stereo: sci_form::graph::BondStereo::None,
        },
    );
    mol.add_bond(
        a0,
        a3,
        sci_form::graph::Bond {
            order: sci_form::graph::BondOrder::Single,
            stereo: sci_form::graph::BondStereo::None,
        },
    );

    // Planar coordinates with very small volume (~0.05)
    let mut coords = DMatrix::from_element(4, 3, 0.0);
    coords[(0, 0)] = 0.0;
    coords[(0, 1)] = 0.0;
    coords[(0, 2)] = 0.0;
    coords[(1, 0)] = 1.0;
    coords[(1, 1)] = 0.0;
    coords[(1, 2)] = 0.0;
    coords[(2, 0)] = 0.0;
    coords[(2, 1)] = 1.0;
    coords[(2, 2)] = 0.0;
    coords[(3, 0)] = 0.02;
    coords[(3, 1)] = 0.02;
    coords[(3, 2)] = 0.05;

    let params = FFParams {
        kb: 0.0,
        k_theta: 0.0,
        k_omega: 0.0,
        k_oop: 0.0,
        k_bounds: 0.0,
        k_chiral: 100.0,
        k_vdw: 0.0,
    };

    let n = 4;
    let bounds_matrix = DMatrix::from_element(n, n, 0.0);

    let e0 = calculate_total_energy(&coords, &mol, &params, &bounds_matrix);
    println!("Initial Total Energy: {}", e0);

    let analytical_grad = compute_analytical_gradient(&coords, &mol, &params, &bounds_matrix);

    let h = 1e-4;
    let mut numerical_grad = DMatrix::from_element(n, 3, 0.0);

    for i in 0..n {
        for dim in 0..3 {
            let mut c_plus = coords.clone();
            c_plus[(i, dim)] += h;
            let mut c_minus = coords.clone();
            c_minus[(i, dim)] -= h;

            let e_plus = calculate_total_energy(&c_plus, &mol, &params, &bounds_matrix);
            let e_minus = calculate_total_energy(&c_minus, &mol, &params, &bounds_matrix);

            numerical_grad[(i, dim)] = (e_plus - e_minus) / (2.0 * h);
        }
    }

    println!("Analytical vs Numerical Chiral Gradients:");
    let mut max_diff: f32 = 0.0;
    for i in 0..n {
        for dim in 0..3 {
            let a = analytical_grad[(i, dim)];
            let num = numerical_grad[(i, dim)];
            let diff = (a - num).abs();
            if diff > max_diff {
                max_diff = diff;
            }
            println!(
                "Atom {}, dim {}: An: {:.5}, Num: {:.5}, Diff: {:.5}",
                i, dim, a, num, diff
            );
        }
    }
    println!("Max diff: {}", max_diff);
}
