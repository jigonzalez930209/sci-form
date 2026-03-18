use sci_form::hf::basis::build_sto3g_basis;
use sci_form::hf::integrals::{compute_eris, get_eri};

fn main() {
    let elements = [6u8, 1, 1, 1, 1];
    let positions = [
        [0.0, 0.0, 0.0],
        [0.6276, 0.6276, 0.6276],
        [-0.6276, -0.6276, 0.6276],
        [-0.6276, 0.6276, -0.6276],
        [0.6276, -0.6276, -0.6276],
    ];
    let basis = build_sto3g_basis(&elements, &positions);
    let eris = compute_eris(&basis);
    let n = basis.n_basis();
    
    // Print the specific ERIs that PySCF shows for lam=2, sig varies
    println!("lam=2 (px):");
    for sig in 5..=8 {
        let j55 = get_eri(&eris, 5, 5, 2, sig, n);
        let k55 = get_eri(&eris, 5, 2, 5, sig, n);
        let j66 = get_eri(&eris, 6, 6, 2, sig, n);
        let k66 = get_eri(&eris, 6, 2, 6, sig, n);
        println!("  sig={}: J55={:.8} K55={:.8} J66={:.8} K66={:.8}", sig, j55, k55, j66, k66);
    }
}
