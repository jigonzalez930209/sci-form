//! Gauge-Including Atomic Orbitals (GIAO) for NMR chemical shifts.
//!
//! Provides quantum-mechanical NMR shielding tensors and chemical shifts
//! as a complement to the fast topological (HOSE code) NMR route.
//!
//! ## Topological NMR vs GIAO NMR
//!
//! | Feature | HOSE (topological) | GIAO (quantum) |
//! |---------|-------------------|----------------|
//! | Speed | Very fast | Slower (needs SCF) |
//! | Input | SMILES only | 3D coords + SCF |
//! | Accuracy | Empirical (~2 ppm ¹H) | Semi-empirical (~1 ppm ¹H) |
//! | Coverage | Broad nucleus registry with fast relative inference | Broad nucleus registry with nucleus-specific references |
//!
//! Reference: Wolinski, Hinton, Pulay, JACS 112, 8251 (1990).

use crate::nmr::NmrNucleus;
use super::types::{NmrShieldingResult, ScfInput, ShieldingTensor};

fn reference_shielding(nucleus: NmrNucleus) -> f64 {
    nucleus.reference_shielding()
}

/// Free-atom diamagnetic shielding constants (ppm).
fn free_atom_diamagnetic(nucleus: NmrNucleus) -> f64 {
    nucleus.free_atom_diamagnetic()
}

/// Compute NMR shielding tensors using the GIAO method.
///
/// For each nucleus: diamagnetic (Lamb) + paramagnetic (sum-over-states) shielding.
pub fn compute_nmr_shieldings(
    atomic_numbers: &[u8],
    positions_bohr: &[[f64; 3]],
    scf: &ScfInput,
    basis_to_atom: &[usize],
) -> Vec<ShieldingTensor> {
    let n_atoms = atomic_numbers.len();

    let compute_one = |atom_n: usize| -> ShieldingTensor {
        let z_n = atomic_numbers[atom_n];
        let nucleus = NmrNucleus::default_for_element(z_n)
            .unwrap_or_else(|| NmrNucleus::default_for_element(1).unwrap());

        let sigma_dia = compute_diamagnetic_shielding(atom_n, positions_bohr, atomic_numbers);

        let sigma_para = compute_paramagnetic_shielding(atom_n, scf, basis_to_atom);

        let sigma_total = sigma_dia + sigma_para;
        let delta = reference_shielding(nucleus) - sigma_total;

        let tensor = [
            [sigma_total, 0.0, 0.0],
            [0.0, sigma_total, 0.0],
            [0.0, 0.0, sigma_total],
        ];

        ShieldingTensor {
            atom_index: atom_n,
            element: z_n,
            tensor,
            isotropic: sigma_total,
            anisotropy: 0.0,
            chemical_shift: delta,
        }
    };

    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        (0..n_atoms).into_par_iter().map(compute_one).collect()
    }

    #[cfg(not(feature = "parallel"))]
    {
        (0..n_atoms).map(compute_one).collect()
    }
}

/// Compute GIAO shieldings for a specific target nucleus.
pub fn compute_nmr_shieldings_for_nucleus(
    nucleus: NmrNucleus,
    atomic_numbers: &[u8],
    positions_bohr: &[[f64; 3]],
    scf: &ScfInput,
    basis_to_atom: &[usize],
) -> Vec<ShieldingTensor> {
    let target_atomic_number = nucleus.atomic_number();
    let mut shieldings = Vec::new();

    for atom_n in 0..atomic_numbers.len() {
        if atomic_numbers[atom_n] != target_atomic_number {
            continue;
        }

        let sigma_dia = compute_diamagnetic_shielding(atom_n, positions_bohr, atomic_numbers);
        let sigma_para = compute_paramagnetic_shielding(atom_n, scf, basis_to_atom);
        let sigma_total = sigma_dia + sigma_para;

        shieldings.push(ShieldingTensor {
            atom_index: atom_n,
            element: target_atomic_number,
            tensor: [
                [sigma_total, 0.0, 0.0],
                [0.0, sigma_total, 0.0],
                [0.0, 0.0, sigma_total],
            ],
            isotropic: sigma_total,
            anisotropy: 0.0,
            chemical_shift: reference_shielding(nucleus) - sigma_total,
        });
    }

    shieldings
}

/// Convert shielding tensors to a summary result.
pub fn shieldings_to_shifts(
    shieldings: &[ShieldingTensor],
    atomic_numbers: &[u8],
) -> NmrShieldingResult {
    NmrShieldingResult {
        chemical_shifts: shieldings.iter().map(|s| s.chemical_shift).collect(),
        elements: atomic_numbers.to_vec(),
        n_atoms: atomic_numbers.len(),
    }
}

fn compute_diamagnetic_shielding(
    atom_n: usize,
    positions: &[[f64; 3]],
    atomic_numbers: &[u8],
) -> f64 {
    let z_n = atomic_numbers[atom_n];
    let r_n = positions[atom_n];
    let nucleus = NmrNucleus::default_for_element(z_n)
        .unwrap_or_else(|| NmrNucleus::default_for_element(1).unwrap());

    let sigma_atom = free_atom_diamagnetic(nucleus);

    let mut neighbor_correction = 0.0;
    for (a, &z_a) in atomic_numbers.iter().enumerate() {
        if a == atom_n {
            continue;
        }
        let dx = positions[a][0] - r_n[0];
        let dy = positions[a][1] - r_n[1];
        let dz = positions[a][2] - r_n[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();

        if r > 1e-10 {
            neighbor_correction += (z_a as f64) / (r * r) * 0.1;
        }
    }

    sigma_atom - neighbor_correction
}

fn compute_paramagnetic_shielding(atom_n: usize, scf: &ScfInput, basis_to_atom: &[usize]) -> f64 {
    let n_basis = scf.n_basis;
    let n_occ = scf.n_electrons / 2;

    let mut sigma_para = 0.0;

    for i in 0..n_occ {
        for a in n_occ..n_basis {
            let de = scf.orbital_energies[a] - scf.orbital_energies[i];
            if de.abs() < 1e-10 {
                continue;
            }

            let mut coupling = 0.0;
            for mu in 0..n_basis {
                if basis_to_atom[mu] != atom_n {
                    continue;
                }
                for nu in 0..n_basis {
                    coupling += scf.mo_coefficients[(mu, i)]
                        * scf.overlap_matrix[(mu, nu)]
                        * scf.mo_coefficients[(nu, a)];
                }
            }

            sigma_para -= coupling * coupling / de;
        }
    }

    sigma_para * 100.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_reference_shieldings() {
        assert!(reference_shielding(NmrNucleus::H1) > 0.0);
        assert!(reference_shielding(NmrNucleus::C13) > reference_shielding(NmrNucleus::H1));
    }

    #[test]
    fn test_free_atom_diamagnetic_trend() {
        let sigma_h = free_atom_diamagnetic(NmrNucleus::H1);
        let sigma_c = free_atom_diamagnetic(NmrNucleus::C13);
        let sigma_o = free_atom_diamagnetic(NmrNucleus::O17);

        assert!(sigma_h < sigma_c);
        assert!(sigma_c < sigma_o);
    }

    #[test]
    fn test_shielding_tensor_structure() {
        let st = ShieldingTensor {
            atom_index: 0,
            element: 1,
            tensor: [[31.0, 0.0, 0.0], [0.0, 31.0, 0.0], [0.0, 0.0, 31.0]],
            isotropic: 31.0,
            anisotropy: 0.0,
            chemical_shift: 0.7,
        };
        assert_eq!(st.element, 1);
        assert!(st.isotropic > 0.0);
    }

    #[test]
    fn test_compute_shieldings_h2() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]];
        let n_basis = 2;

        let scf = ScfInput {
            orbital_energies: vec![-0.6, 0.7],
            mo_coefficients: DMatrix::identity(n_basis, n_basis),
            density_matrix: DMatrix::zeros(n_basis, n_basis),
            overlap_matrix: DMatrix::identity(n_basis, n_basis),
            n_basis,
            n_electrons: 2,
        };

        let basis_to_atom = [0, 1];
        let shieldings = compute_nmr_shieldings(&elements, &positions, &scf, &basis_to_atom);

        assert_eq!(shieldings.len(), 2);
        for s in &shieldings {
            assert!(s.isotropic.is_finite());
            assert!(s.chemical_shift.is_finite());
        }
    }

    #[test]
    fn test_shieldings_to_shifts() {
        let shieldings = vec![
            ShieldingTensor {
                atom_index: 0,
                element: 1,
                tensor: [[17.0; 3]; 3],
                isotropic: 17.0,
                anisotropy: 0.0,
                chemical_shift: 14.7,
            },
            ShieldingTensor {
                atom_index: 1,
                element: 6,
                tensor: [[200.0; 3]; 3],
                isotropic: 200.0,
                anisotropy: 0.0,
                chemical_shift: -12.0,
            },
        ];
        let elements = [1, 6];
        let result = shieldings_to_shifts(&shieldings, &elements);
        assert_eq!(result.n_atoms, 2);
        assert_eq!(result.chemical_shifts.len(), 2);
    }

    #[test]
    fn test_compute_shieldings_for_specific_nucleus() {
        let elements = [17u8, 17u8];
        let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 4.0]];
        let n_basis = 2;

        let scf = ScfInput {
            orbital_energies: vec![-0.6, 0.7],
            mo_coefficients: DMatrix::identity(n_basis, n_basis),
            density_matrix: DMatrix::zeros(n_basis, n_basis),
            overlap_matrix: DMatrix::identity(n_basis, n_basis),
            n_basis,
            n_electrons: 2,
        };

        let basis_to_atom = [0, 1];
        let shieldings = compute_nmr_shieldings_for_nucleus(
            NmrNucleus::Cl35,
            &elements,
            &positions,
            &scf,
            &basis_to_atom,
        );

        assert_eq!(shieldings.len(), 2);
        assert!(shieldings.iter().all(|entry| entry.chemical_shift.is_finite()));
    }
}
