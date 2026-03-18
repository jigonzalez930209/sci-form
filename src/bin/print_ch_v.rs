use sci_form::hf::basis::{build_sto3g_basis, ANG_TO_BOHR};
use sci_form::hf::nuclear::compute_nuclear_matrix;

fn main() {
    let elements = [6, 1];
    let pos_angstrom = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let basis = build_sto3g_basis(&elements, &pos_angstrom);
    let mut pos_bohr = Vec::new();
    for p in &pos_angstrom {
        pos_bohr.push([p[0] * ANG_TO_BOHR, p[1] * ANG_TO_BOHR, p[2] * ANG_TO_BOHR]);
    }

    let v = compute_nuclear_matrix(&basis, &elements, &pos_bohr);

    println!("V[0,0]: {}", v[(0, 0)]);
    println!("V[2,2]: {}", v[(2, 2)]);
    println!("V[0,5]: {}", v[(0, 5)]);
    println!("V[2,5]: {}", v[(2, 5)]);
}
