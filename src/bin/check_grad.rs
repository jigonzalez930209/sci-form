use sci_form::forcefield::{FFParams, energy::calculate_total_energy, gradients::compute_analytical_gradient};
use sci_form::graph::{Atom, Molecule};
use nalgebra::{DMatrix, Vector3};

fn main() {
    let mut mol = Molecule::new("test");
    // Create a 4 atom molecule
    let a1 = mol.add_atom(Atom::new(6, 0.0, 0.0, 0.0));
    let a2 = mol.add_atom(Atom::new(6, 1.5, 0.0, 0.0));
    let a3 = mol.add_atom(Atom::new(6, 1.5, 1.5, 0.0));
    let a4 = mol.add_atom(Atom::new(6, 0.0, 1.5, 0.5));
    mol.add_bond(a1, a2, sci_form::graph::Bond { order: sci_form::graph::BondOrder::Single, stereo: sci_form::graph::BondStereo::None });
    mol.add_bond(a2, a3, sci_form::graph::Bond { order: sci_form::graph::BondOrder::Single, stereo: sci_form::graph::BondStereo::None });
    mol.add_bond(a3, a4, sci_form::graph::Bond { order: sci_form::graph::BondOrder::Single, stereo: sci_form::graph::BondStereo::None });

    let mut coords = DMatrix::from_element(4, 3, 0.0);
    coords[(0, 0)] = 0.0; coords[(0, 1)] = 0.0; coords[(0, 2)] = 0.0;
    coords[(1, 0)] = 1.3; coords[(1, 1)] = 0.1; coords[(1, 2)] = 0.0;
    coords[(2, 0)] = 1.5; coords[(2, 1)] = 1.6; coords[(2, 2)] = 0.1;
    coords[(3, 0)] = 0.1; coords[(3, 1)] = 1.2; coords[(3, 2)] = 0.7;

    let params = FFParams {
        kb: 0.0,
        k_theta: 200.0,
        k_omega: 50.0,
        k_oop: 50.0,
        k_bounds: 100.0,
    };

    let n = 4;
    let mut bounds_matrix = DMatrix::from_element(n, n, 0.0);
    for i in 0..n {
        for j in 0..n {
            bounds_matrix[(i, j)] = (i as f32 - j as f32).abs() * 1.5;
        }
    }

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

    println!("Analytical vs Numerical Gradients:");
    let mut max_diff: f32 = 0.0;
    for i in 0..n {
        for dim in 0..3 {
            let a = analytical_grad[(i, dim)];
            let num = numerical_grad[(i, dim)];
            let diff = (a - num).abs();
            if diff > max_diff { max_diff = diff; }
            println!("Atom {}, dim {}: An: {:.5}, Num: {:.5}, Diff: {:.5}", i, dim, a, num, diff);
        }
    }
    println!("Max diff: {}", max_diff);
}
