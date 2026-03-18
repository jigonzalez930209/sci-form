use sci_form::hf::basis::{build_sto3g_basis};
use sci_form::hf::integrals::{compute_eris, get_eri};
use sci_form::hf::overlap_kin::compute_overlap_matrix;

fn main() {
    let elements = [6, 1];
    let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let basis = build_sto3g_basis(&elements, &positions);
    let eris = compute_eris(&basis);
    let o = compute_overlap_matrix(&basis);
    let n = basis.n_basis();
    
    // Carbon: 1s(0), 2s(1), 2px(2), 2py(3), 2pz(4)
    // Hydrogen: 1s(5)
    
    println!("ERI[2,5,0,0]: {}", get_eri(&eris, 2, 5, 0, 0, n));
    println!("ERI[2,5,2,5]: {}", get_eri(&eris, 2, 5, 2, 5, n));
    println!("ERI[2,2,0,0]: {}", get_eri(&eris, 2, 2, 0, 0, n));
    println!("ERI[2,0,2,0]: {}", get_eri(&eris, 2, 0, 2, 0, n));
    println!("S[2,5]: {}", o[(2, 5)]);
    println!("S[0,5]: {}", o[(0, 5)]);
}
