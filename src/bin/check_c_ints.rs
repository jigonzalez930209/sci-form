use sci_form::hf::basis::build_sto3g_basis;
use sci_form::hf::overlap_kin::{compute_overlap_matrix, compute_kinetic_matrix};

fn main() {
    let elements = [6];
    let positions = [[0.0, 0.0, 0.0]];
    let basis = build_sto3g_basis(&elements, &positions);
    let s = compute_overlap_matrix(&basis);
    let t = compute_kinetic_matrix(&basis);
    for i in 0..s.nrows() {
        println!("S[{},{}] = {}", i, i, s[(i, i)]);
        println!("T[{},{}] = {}", i, i, t[(i, i)]);
    }
}
