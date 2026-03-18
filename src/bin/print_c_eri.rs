use sci_form::hf::basis::build_sto3g_basis;
use sci_form::hf::integrals::{compute_eris, get_eri};

fn main() {
    let elements = [6];
    let positions = [[0.0, 0.0, 0.0]];
    let basis = build_sto3g_basis(&elements, &positions);
    let eris = compute_eris(&basis);
    let n = basis.n_basis();
    
    // px is 2, py is 3
    println!("ERI[2,3,2,3] : {}", get_eri(&eris, 2, 3, 2, 3, n));
    println!("ERI[2,3,3,2] : {}", get_eri(&eris, 2, 3, 3, 2, n));
}
