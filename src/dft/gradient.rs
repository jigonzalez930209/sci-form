//! KS nuclear gradient: frozen-P integral derivatives plus Pulay −W:∂S/∂R.
//!
//! Integrals are differentiated by central differences of the analytic
//! H_core, ERI, XC, and overlap at the converged density. That is the
//! analytic KS force with numerically differentiated AO integrals.

use crate::scf::core_matrices::{nuclear_repulsion_energy_charges, CoreMatrices};
use crate::scf::gradients::GradientResult;
use crate::scf::overlap_matrix::build_overlap_matrix;
use crate::scf::two_electron::TwoElectronIntegrals;
use crate::scf::types::MolecularSystem;

use super::ecp;
use super::grid::MolecularGrid;
use super::ks_fock::{solve_ks_dft, DftConfig, DftResult};
use super::vxc_matrix::{build_vxc_matrix, build_vxc_uks};
use nalgebra::DMatrix;

const STEP: f64 = 1.0e-4;

/// Kohn-Sham nuclear gradient at the converged density.
///
/// The gradient is computed as a *frozen-density* derivative: the AO
/// integrals (core Hamiltonian, ERIs, exchange-correlation potential and the
/// overlap) are re-evaluated at the displaced geometry, but the converged
/// density matrix is held fixed, and the Pulay term `−W:∂S/∂R` is added
/// afterwards. Atomic integrals are differentiated by central finite
/// differences with a step of `1e-4` Bohr, so this is an analytic Kohn-Sham
/// force built from numerically differentiated integrals.
///
/// `elements` and `positions` describe the molecule in the usual convention
/// of this crate: atomic numbers and Cartesian coordinates in Ångström.
///
/// # Errors
///
/// Returns `Err` if the Kohn-Sham SCF fails (see
/// [`super::ks_fock::solve_ks_dft`]), if the geometry is invalid, or if the
/// requested basis set is unavailable for one of the elements. The error
/// string of the underlying SCF failure is propagated unchanged.
///
/// # Example
///
/// ```no_run
/// use sci_form::dft::gradient::ks_nuclear_gradient;
/// use sci_form::dft::ks_fock::DftConfig;
///
/// fn main() -> Result<(), String> {
///     let gradient = ks_nuclear_gradient(
///         &[1, 1],
///         &[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]],
///         &DftConfig::default(),
///     )?;
///     // Gradients are dE/dR in Hartree/Bohr, one [x, y, z] triple per atom.
///     println!("{:?}", gradient.gradients);
///     Ok(())
/// }
/// ```
pub fn ks_nuclear_gradient(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &DftConfig,
) -> Result<GradientResult, String> {
    let result = solve_ks_dft(elements, positions, config)?;
    frozen_p_gradient(elements, positions, config, &result)
}

fn reshape(flat: &[f64], n: usize) -> DMatrix<f64> {
    DMatrix::from_row_slice(n, n, flat)
}

fn frozen_p_gradient(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &DftConfig,
    result: &DftResult,
) -> Result<GradientResult, String> {
    let system = MolecularSystem::from_angstrom(
        elements,
        positions,
        config.charge,
        config.multiplicity.max(1),
    );
    let n = result.n_basis;
    let p = reshape(&result.density, n);
    let c = reshape(&result.coefficients, n);
    // Spin densities and the energy-weighted density. Restricted: one set of
    // doubly occupied orbitals (factor 2). Unrestricted: each spin channel
    // contributes its own singly occupied orbitals (factor 1).
    let (spin, w) = if result.unrestricted {
        let c_b = reshape(&result.coefficients_beta, n);
        let p_a = spin_density(&c, result.n_alpha);
        let p_b = spin_density(&c_b, result.n_beta);
        let w = energy_weighted(&c, &result.orbital_energies, result.n_alpha)
            + energy_weighted(&c_b, &result.orbital_energies_beta, result.n_beta);
        (Some((p_a, p_b)), w)
    } else {
        let n_occ = result.n_alpha.max(result.n_electrons / 2);
        (
            None,
            crate::scf::gradients::energy_weighted_density(&c, &result.orbital_energies, n_occ),
        )
    };

    let n_atoms = elements.len();
    let mut gradients = vec![[0.0; 3]; n_atoms];
    for atom in 0..n_atoms {
        for coord in 0..3 {
            let e_plus = frozen_energy(&system, config, &p, spin.as_ref(), atom, coord, STEP)?;
            let e_minus = frozen_energy(&system, config, &p, spin.as_ref(), atom, coord, -STEP)?;
            let s_plus = overlap_at(&system, config, atom, coord, STEP)?;
            let s_minus = overlap_at(&system, config, atom, coord, -STEP)?;
            let ds = (&s_plus - &s_minus) / (2.0 * STEP);
            let pulay = (&w * &ds).trace();
            gradients[atom][coord] = (e_plus - e_minus) / (2.0 * STEP) - pulay;
        }
    }

    let mut sum_sq = 0.0;
    let mut max_g = 0.0f64;
    for g in &gradients {
        for &gc in g {
            sum_sq += gc * gc;
            max_g = max_g.max(gc.abs());
        }
    }
    Ok(GradientResult {
        gradients,
        rms_gradient: (sum_sq / (n_atoms * 3) as f64).sqrt(),
        max_gradient: max_g,
        energy: result.total_energy,
    })
}

fn spin_density(c: &DMatrix<f64>, n_occ: usize) -> DMatrix<f64> {
    let occ = c.columns(0, n_occ);
    occ * occ.transpose()
}

/// `Σ_i ε_i c_i c_iᵀ` over the first `n_occ` orbitals of one spin channel.
fn energy_weighted(c: &DMatrix<f64>, eps: &[f64], n_occ: usize) -> DMatrix<f64> {
    let mut w = DMatrix::zeros(c.nrows(), c.nrows());
    for (k, &e) in eps.iter().enumerate().take(n_occ) {
        let col = c.column(k);
        w += (col * col.transpose()) * e;
    }
    w
}

fn displaced(system: &MolecularSystem, atom: usize, coord: usize, step: f64) -> MolecularSystem {
    let mut sys = system.clone();
    sys.positions_bohr[atom][coord] += step;
    sys
}

fn overlap_at(
    system: &MolecularSystem,
    config: &DftConfig,
    atom: usize,
    coord: usize,
    step: f64,
) -> Result<DMatrix<f64>, String> {
    let sys = displaced(system, atom, coord, step);
    let basis = super::valence_basis::ks_basis_kind(
        &sys.atomic_numbers,
        &sys.positions_bohr,
        config.basis,
    )?;
    Ok(build_overlap_matrix(&basis))
}

fn frozen_energy(
    system: &MolecularSystem,
    config: &DftConfig,
    density: &DMatrix<f64>,
    spin: Option<&(DMatrix<f64>, DMatrix<f64>)>,
    atom: usize,
    coord: usize,
    step: f64,
) -> Result<f64, String> {
    let sys = displaced(system, atom, coord, step);
    let basis = super::valence_basis::ks_basis_kind(
        &sys.atomic_numbers,
        &sys.positions_bohr,
        config.basis,
    )?;
    let z_eff = ecp::effective_charges(&sys.atomic_numbers);
    let mut core = CoreMatrices::build_with_charges(&basis, &z_eff, &sys.positions_bohr);
    ecp::add_local_ecp(&mut core.nuclear, &basis, &sys.atomic_numbers);
    core.core_hamiltonian = &core.kinetic + &core.nuclear;
    let eri = TwoElectronIntegrals::compute(&basis);
    let n = basis.n_basis;
    if n != density.nrows() {
        return Err("basis size changed under displacement".to_string());
    }
    let mut j_energy = 0.0;
    for mu in 0..n {
        for nu in 0..n {
            let mut j = 0.0;
            for lam in 0..n {
                for sig in 0..n {
                    j += density[(lam, sig)] * eri.get(mu, nu, lam, sig);
                }
            }
            j_energy += density[(mu, nu)] * j;
        }
    }
    j_energy *= 0.5;
    let grid = MolecularGrid::build(
        &sys.atomic_numbers,
        &sys.positions_bohr,
        config.grid_quality,
    );
    let exc = match spin {
        Some((p_a, p_b)) => build_vxc_uks(&basis, p_a, p_b, &grid, config.method).2,
        None => build_vxc_matrix(&basis, density, &grid, config.method).1,
    };
    let e1 = (density * &core.core_hamiltonian).trace();
    let vnn = nuclear_repulsion_energy_charges(&z_eff, &sys.positions_bohr);
    Ok(e1 + j_energy + exc + vnn)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dft::grid::GridQuality;
    use crate::dft::ks_fock::DftMethod;

    #[test]
    fn h2_ks_gradient_opposes_the_bond_stretch() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let config = DftConfig {
            method: DftMethod::Svwn,
            grid_quality: GridQuality::Coarse,
            max_iterations: 40,
            ..DftConfig::default()
        };
        let g = ks_nuclear_gradient(&elements, &positions, &config).unwrap();
        assert_eq!(g.gradients.len(), 2);
        // At 0.74 Å H2 is near equilibrium for this model; forces must be finite
        // and opposite along the bond (Newton III).
        assert!((g.gradients[0][0] + g.gradients[1][0]).abs() < 5e-3);
        assert!(g.gradients[0][0].is_finite() && g.gradients[1][0].is_finite());
    }

    /// The gradient must equal the derivative of the reported energy, in
    /// Hartree/**Bohr**. (Until v1.0 the F80 gate compared it with a
    /// derivative taken per Å and concluded the analytic gradient was off by
    /// a factor of 0.5292 = Å per Bohr; the gradient was right, the
    /// reference was in the wrong unit.)
    #[test]
    fn h2_ks_gradient_is_the_energy_derivative_in_ha_per_bohr() {
        use crate::scf::constants::BOHR_TO_ANGSTROM;
        let elements = [1u8, 1];
        let config = DftConfig {
            method: DftMethod::Svwn,
            grid_quality: GridQuality::Coarse,
            max_iterations: 200,
            energy_threshold: 1e-10,
            density_threshold: 1e-9,
            ..DftConfig::default()
        };
        let at = |x: f64| [[0.0, 0.0, 0.0], [x, 0.0, 0.0]];
        let g = ks_nuclear_gradient(&elements, &at(0.90), &config).unwrap();
        let h = 2.0e-4; // Å
        let e = |x: f64| {
            solve_ks_dft(&elements, &at(x), &config)
                .unwrap()
                .total_energy
        };
        let per_angstrom = (e(0.90 + h) - e(0.90 - h)) / (2.0 * h);
        let per_bohr = per_angstrom * BOHR_TO_ANGSTROM;
        assert!(
            (g.gradients[1][0] - per_bohr).abs() < 1e-5,
            "gradient {:.7} Ha/Bohr vs energy derivative {:.7} Ha/Bohr",
            g.gradients[1][0],
            per_bohr
        );
        assert!(g.gradients[1][0] > 0.0, "H₂ at 0.90 Å is past the minimum");
    }
}
