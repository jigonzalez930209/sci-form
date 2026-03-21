//! Energy evaluation from density and Fock matrices.
//!
//! Electronic energy:
//!   E_elec = 0.5 · Tr[P(H⁰ + F)]
//!
//! Total energy:
//!   E_total = E_elec + E_nuc

use nalgebra::DMatrix;

/// Compute the electronic energy from density and matrices.
///
/// E_elec = 0.5 · Tr[P · (H + F)]
pub fn electronic_energy(
    density: &DMatrix<f64>,
    h_core: &DMatrix<f64>,
    fock: &DMatrix<f64>,
) -> f64 {
    let n = density.nrows();
    let hf = h_core + fock;
    let mut energy = 0.0;

    for i in 0..n {
        for j in 0..n {
            energy += density[(i, j)] * hf[(i, j)];
        }
    }

    0.5 * energy
}

/// Compute total energy = electronic + nuclear repulsion.
pub fn total_energy(
    density: &DMatrix<f64>,
    h_core: &DMatrix<f64>,
    fock: &DMatrix<f64>,
    nuclear_repulsion: f64,
) -> f64 {
    electronic_energy(density, h_core, fock) + nuclear_repulsion
}

/// Compute the one-electron energy contribution.
///
/// E_1e = Tr[P · H⁰]
pub fn one_electron_energy(density: &DMatrix<f64>, h_core: &DMatrix<f64>) -> f64 {
    let n = density.nrows();
    let mut energy = 0.0;
    for i in 0..n {
        for j in 0..n {
            energy += density[(i, j)] * h_core[(i, j)];
        }
    }
    energy
}

/// Compute the two-electron energy contribution.
///
/// E_2e = 0.5 · Tr[P · G]
/// where G = F - H⁰ is the two-electron part of the Fock matrix.
pub fn two_electron_energy(
    density: &DMatrix<f64>,
    h_core: &DMatrix<f64>,
    fock: &DMatrix<f64>,
) -> f64 {
    let g = fock - h_core;
    let n = density.nrows();
    let mut energy = 0.0;
    for i in 0..n {
        for j in 0..n {
            energy += density[(i, j)] * g[(i, j)];
        }
    }
    0.5 * energy
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_energy_decomposition() {
        let h = DMatrix::from_row_slice(2, 2, &[-1.0, -0.3, -0.3, -0.8]);
        let f = DMatrix::from_row_slice(2, 2, &[-0.9, -0.2, -0.2, -0.7]);
        let p = DMatrix::from_row_slice(2, 2, &[1.5, 0.2, 0.2, 0.5]);

        let e_1e = one_electron_energy(&p, &h);
        let e_2e = two_electron_energy(&p, &h, &f);
        let e_elec = electronic_energy(&p, &h, &f);

        assert!((e_elec - (e_1e + e_2e)).abs() < 1e-12);
    }

    #[test]
    fn test_total_energy_includes_nuclear() {
        let h = DMatrix::from_row_slice(2, 2, &[-1.0, 0.0, 0.0, -1.0]);
        let f = h.clone();
        let p = DMatrix::identity(2, 2);

        let e_nuc = 0.5;
        let e_total = total_energy(&p, &h, &f, e_nuc);
        let e_elec = electronic_energy(&p, &h, &f);

        assert!((e_total - (e_elec + e_nuc)).abs() < 1e-14);
    }
}
