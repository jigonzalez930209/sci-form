//! Alpha periodic adapters needing gated beta solvers (KPM/RandNLA).
//!
//! Meshes, Bloch phases, and adapters live in core [`crate::kmesh`]
//! (re-exported below).
//!
//! # Stability tier
//!
//! **Alpha** — early research tier. Gated behind the
//! `alpha-periodic-linear` feature flag (alias: `experimental-periodic-linear`,
//! which requires the `beta-kpm` and `beta-randnla` solvers). The validation
//! gate is `tests/alpha_periodic_validation.rs` (F88): the k-mesh and band
//! structure are compared against ASE. Promotion criteria are listed in
//! [`crate::experimental_status`].
//!
//! # Units
//!
//! Energies and Fermi levels in **eV**, the density of states in
//! **states/eV**, and the k-points in **Å⁻¹**.

#[cfg(feature = "experimental-gpu")]
pub mod gpu;

pub use crate::kmesh::*;

use nalgebra::DMatrix;
/// Compute a k-averaged DOS by applying beta KPM to one operator per k-point.
///
/// The energies are in **eV** and the density of states in **states/eV**.
///
/// # Errors
///
/// - `"at least one Hamiltonian operator is required"` when `operators` is
///   empty.
/// - `"periodic KPM output requires at least 2 energy points"` when
///   `config.n_points < 2`.
/// - `"number of operators must match the number of k-points"` when the
///   operator count differs from the k-mesh size.
/// - Any error from the k-mesh builder.
pub fn compute_periodic_kpm_dos_from_operators(
    operators: &[DMatrix<f64>],
    config: &PeriodicKpmAdapterConfig,
) -> Result<PeriodicKpmDosResult, String> {
    if operators.is_empty() {
        return Err("at least one Hamiltonian operator is required".into());
    }
    if config.n_points < 2 {
        return Err("periodic KPM output requires at least 2 energy points".into());
    }

    let kmesh = monkhorst_pack_mesh(&config.kmesh)?;
    if operators.len() != kmesh.points.len() {
        return Err("number of operators must match the number of k-points".into());
    }

    let beta_config = crate::beta::kpm::KpmConfig {
        order: config.order,
        n_vectors: config.n_vectors,
        seed: config.seed,
        temperature: config.temperature_k,
    };

    let first = crate::beta::kpm::compute_kpm_dos(
        &operators[0],
        &beta_config,
        config.e_min_ev,
        config.e_max_ev,
        config.n_points,
    );
    let energies_ev = first.energies.clone();
    let mut total_dos = vec![0.0; first.total_dos.len()];

    for (operator, point) in operators.iter().zip(&kmesh.points) {
        let result = crate::beta::kpm::compute_kpm_dos(
            operator,
            &beta_config,
            config.e_min_ev,
            config.e_max_ev,
            config.n_points,
        );
        for (accumulator, value) in total_dos.iter_mut().zip(result.total_dos) {
            *accumulator += point.weight * value;
        }
    }

    Ok(PeriodicKpmDosResult {
        band_edges: estimate_band_edges_from_dos(&energies_ev, &total_dos, config.fermi_energy_ev),
        energies_ev,
        total_dos,
        kmesh: Some(kmesh.clone()),
        diagnostics: PeriodicSpectralDiagnostics {
            n_kpoints: kmesh.points.len(),
            polynomial_order: Some(config.order),
            sketch_rank: None,
            max_residual: 0.0,
            used_fallback: false,
        },
    })
}

/// Solve one generalized eigenproblem per k-point with the beta RandNLA solver.
///
/// The Hamiltonian and overlap matrices are in **eV**, the band edges and
/// the Fermi level in **eV**.
///
/// # Errors
///
/// - `"at least one Hamiltonian matrix is required"` when `hamiltonians` is
///   empty.
/// - `"Hamiltonian and overlap lists must have the same length"` when the two
///   lists differ in length.
/// - `"number of Hamiltonians must match the number of k-points"` when the
///   matrix count differs from the k-mesh size.
/// - Any error from the k-mesh builder.
pub fn solve_periodic_randnla(
    hamiltonians: &[DMatrix<f64>],
    overlaps: &[DMatrix<f64>],
    n_electrons: usize,
    config: &PeriodicRandNlaAdapterConfig,
) -> Result<PeriodicRandNlaResult, String> {
    if hamiltonians.is_empty() {
        return Err("at least one Hamiltonian matrix is required".into());
    }
    if hamiltonians.len() != overlaps.len() {
        return Err("Hamiltonian and overlap lists must have the same length".into());
    }

    let kmesh = monkhorst_pack_mesh(&config.kmesh)?;
    if hamiltonians.len() != kmesh.points.len() {
        return Err("number of Hamiltonians must match the number of k-points".into());
    }

    let n_occupied = n_electrons.div_ceil(2);
    let beta_config = crate::beta::rand_nla::RandNlaConfig {
        sketch_size: config.sketch_size,
        seed: config.seed,
        max_error: config.max_error,
        fallback_enabled: config.fallback_enabled,
    };

    let mut kpoint_results = Vec::with_capacity(kmesh.points.len());
    let mut max_homo = f64::NEG_INFINITY;
    let mut min_lumo = f64::INFINITY;
    let mut direct_gap = f64::INFINITY;
    let mut max_residual = 0.0_f64;
    let mut used_fallback = false;

    for ((hamiltonian, overlap), point) in hamiltonians.iter().zip(overlaps).zip(&kmesh.points) {
        let (energies, _coefficients, info) =
            crate::beta::rand_nla::solve_eht_randnla(hamiltonian, overlap, &beta_config);
        if n_occupied == 0 || n_occupied >= energies.len() {
            return Err("n_electrons does not map to a valid occupied/unoccupied split".into());
        }

        let homo = energies[n_occupied - 1];
        let lumo = energies[n_occupied];
        max_homo = max_homo.max(homo);
        min_lumo = min_lumo.min(lumo);
        direct_gap = direct_gap.min(lumo - homo);
        max_residual = max_residual.max(info.residual_error);
        used_fallback |= info.used_fallback;

        kpoint_results.push(PeriodicRandNlaKPointResult {
            k_fractional: point.fractional,
            weight: point.weight,
            orbital_energies_ev: energies.iter().copied().collect(),
            residual_error: info.residual_error,
            used_fallback: info.used_fallback,
        });
    }

    Ok(PeriodicRandNlaResult {
        band_edges: PeriodicBandEdgeSummary {
            homo_energy_ev: max_homo,
            lumo_energy_ev: min_lumo,
            direct_gap_ev: direct_gap.max(0.0),
            indirect_gap_ev: (min_lumo - max_homo).max(0.0),
        },
        diagnostics: PeriodicSpectralDiagnostics {
            n_kpoints: kmesh.points.len(),
            polynomial_order: None,
            sketch_rank: config.sketch_size,
            max_residual,
            used_fallback,
        },
        kmesh,
        kpoint_results,
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn periodic_kpm_adapter_matches_mesh_size() {
        let h = DMatrix::from_row_slice(2, 2, &[-1.0, -0.2, -0.2, 1.0]);
        let config = PeriodicKpmAdapterConfig {
            kmesh: KMeshConfig {
                grid: [1, 1, 2],
                centering: KMeshCentering::MonkhorstPack,
            },
            n_points: 32,
            e_min_ev: -2.0,
            e_max_ev: 2.0,
            ..Default::default()
        };

        let result = compute_periodic_kpm_dos_from_operators(&[h.clone(), h], &config).unwrap();
        assert_eq!(result.energies_ev.len(), 32);
        assert_eq!(result.kmesh.as_ref().unwrap().points.len(), 2);
        assert_eq!(result.diagnostics.n_kpoints, 2);
    }

    #[test]
    fn periodic_randnla_adapter_reports_band_edges() {
        let h = DMatrix::from_row_slice(2, 2, &[-1.0, -0.2, -0.2, 1.0]);
        let s = DMatrix::identity(2, 2);
        let result = solve_periodic_randnla(
            &[h],
            &[s],
            2,
            &PeriodicRandNlaAdapterConfig {
                kmesh: KMeshConfig::default(),
                sketch_size: Some(2),
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(result.kpoint_results.len(), 1);
        assert!(result.band_edges.direct_gap_ev >= 0.0);
        assert!(result.diagnostics.max_residual < 1e-8);
    }
}
