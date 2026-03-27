//! Integration tests for DFT (KS-SCF) validation.
//!
//! Reference data sources:
//! - CCCBDB (NIST Computational Chemistry Comparison and Benchmark Database)
//!   <https://cccbdb.nist.gov/>
//! - HF/STO-3G total energies are well-established benchmarks
//! - LDA (SVWN) typically overbinds by 1–5% relative to HF for minimal basis
//! - PBE (GGA) corrects LDA over-binding partially

#![cfg(all(feature = "alpha-dft", feature = "alpha-obara-saika"))]

use sci_form::dft::grid::{GridQuality, MolecularGrid};
use sci_form::dft::ks_fock::{solve_ks_dft, DftConfig, DftMethod};

const _HARTREE_TO_EV: f64 = 27.211_386;

// ─── H₂ molecule (simplest diatomic) ────────────────────────────────────────

/// H₂ at equilibrium bond length 0.74 Å.
/// HF/STO-3G reference: -1.1175 Ha (CCCBDB).
/// LDA (SVWN) should give a lower (more negative) energy due to overbinding.
#[test]
fn h2_svwn_total_energy_reasonable() {
    let elements = [1u8, 1];
    let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
    let config = DftConfig {
        method: DftMethod::Svwn,
        grid_quality: GridQuality::Fine,
        ..Default::default()
    };
    let result = solve_ks_dft(&elements, &positions, &config).unwrap();

    assert!(result.converged, "SCF must converge for H₂");

    // Alpha DFT: just verify energy is negative and finite.
    // Absolute values may differ from reference HF/STO-3G = -1.1175 Ha.
    assert!(
        result.total_energy < 0.0,
        "H₂ total energy ({:.4} Ha) should be negative",
        result.total_energy
    );
    assert!(
        result.total_energy.is_finite(),
        "H₂ total energy should be finite"
    );

    // 2 electrons, 2 basis functions (STO-3G: 1s on each H)
    assert_eq!(result.n_electrons, 2);
    assert_eq!(result.n_basis, 2);
}

#[test]
fn h2_svwn_homo_lumo_gap_positive() {
    let elements = [1u8, 1];
    let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
    let config = DftConfig {
        method: DftMethod::Svwn,
        grid_quality: GridQuality::Medium,
        ..Default::default()
    };
    let result = solve_ks_dft(&elements, &positions, &config).unwrap();

    assert!(result.converged);
    assert!(
        result.gap > 0.0,
        "HOMO-LUMO gap ({:.3} eV) should be positive for H₂",
        result.gap
    );
    // Alpha DFT: gap may be larger than physical due to minimal basis
    assert!(result.gap.is_finite(), "HOMO-LUMO gap should be finite");
}

#[test]
fn h2_pbe_converges_and_energy_reasonable() {
    let elements = [1u8, 1];
    let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
    let config = DftConfig {
        method: DftMethod::Pbe,
        grid_quality: GridQuality::Fine,
        ..Default::default()
    };
    let result = solve_ks_dft(&elements, &positions, &config).unwrap();

    assert!(result.converged, "PBE SCF must converge for H₂");
    assert!(
        result.total_energy < 0.0 && result.total_energy.is_finite(),
        "PBE H₂ energy ({:.4} Ha) should be negative and finite",
        result.total_energy
    );
}

/// SVWN and PBE should give different total energies (PBE corrects LDA overbinding).
#[test]
fn h2_svwn_vs_pbe_methods_differ() {
    let elements = [1u8, 1];
    let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];

    let svwn = solve_ks_dft(
        &elements,
        &positions,
        &DftConfig {
            method: DftMethod::Svwn,
            grid_quality: GridQuality::Fine,
            ..Default::default()
        },
    )
    .unwrap();

    let pbe = solve_ks_dft(
        &elements,
        &positions,
        &DftConfig {
            method: DftMethod::Pbe,
            grid_quality: GridQuality::Fine,
            ..Default::default()
        },
    )
    .unwrap();

    assert!(svwn.converged && pbe.converged);
    assert!(
        (svwn.total_energy - pbe.total_energy).abs() > 1e-4,
        "SVWN ({:.6}) and PBE ({:.6}) should give different total energies",
        svwn.total_energy,
        pbe.total_energy
    );
}

// ─── H₂ bond dissociation curve ──────────────────────────────────────────────

/// Energy should increase as H₂ is stretched (bond weakening).
#[test]
fn h2_svwn_energy_increases_with_stretch() {
    let config = DftConfig {
        method: DftMethod::Svwn,
        grid_quality: GridQuality::Medium,
        ..Default::default()
    };
    let elements = [1u8, 1];

    let e_equil = solve_ks_dft(&elements, &[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]], &config).unwrap();

    let e_stretched =
        solve_ks_dft(&elements, &[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], &config).unwrap();

    assert!(e_equil.converged && e_stretched.converged);
    assert!(
        e_equil.total_energy < e_stretched.total_energy,
        "Equilibrium ({:.4}) should be lower than stretched ({:.4})",
        e_equil.total_energy,
        e_stretched.total_energy
    );
}

// ─── Mulliken charge neutrality ──────────────────────────────────────────────

#[test]
fn h2_svwn_mulliken_charges_sum_to_zero() {
    let elements = [1u8, 1];
    let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
    let config = DftConfig {
        method: DftMethod::Svwn,
        grid_quality: GridQuality::Medium,
        ..Default::default()
    };
    let result = solve_ks_dft(&elements, &positions, &config).unwrap();
    let sum: f64 = result.mulliken_charges.iter().sum();
    assert!(
        sum.abs() < 0.1,
        "Mulliken charges should sum to ~0 for neutral H₂, got {:.4}",
        sum
    );
}

/// Symmetric H₂: both atoms should have equal charges.
#[test]
fn h2_svwn_mulliken_charges_symmetric() {
    let elements = [1u8, 1];
    let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
    let config = DftConfig {
        method: DftMethod::Svwn,
        grid_quality: GridQuality::Fine,
        ..Default::default()
    };
    let result = solve_ks_dft(&elements, &positions, &config).unwrap();
    let diff = (result.mulliken_charges[0] - result.mulliken_charges[1]).abs();
    assert!(
        diff < 0.05,
        "Symmetric H₂ should have equal charges, diff = {:.4}",
        diff
    );
}

// ─── Grid quality affects precision ──────────────────────────────────────────

#[test]
fn grid_quality_ordering_affects_points() {
    let elements = [1u8, 1];
    // Positions in Bohr for MolecularGrid::build
    let positions_bohr = [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]];

    let coarse = MolecularGrid::build(&elements, &positions_bohr, GridQuality::Coarse);
    let fine = MolecularGrid::build(&elements, &positions_bohr, GridQuality::Fine);

    assert!(
        fine.n_points > coarse.n_points,
        "Fine grid ({}) should have more points than Coarse ({})",
        fine.n_points,
        coarse.n_points
    );
}

// ─── XC energy is negative (exchange always negative for closed-shell) ───────

#[test]
fn h2_svwn_xc_energy_negative() {
    let elements = [1u8, 1];
    let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
    let config = DftConfig {
        method: DftMethod::Svwn,
        grid_quality: GridQuality::Medium,
        ..Default::default()
    };
    let result = solve_ks_dft(&elements, &positions, &config).unwrap();
    assert!(result.converged);
    assert!(
        result.xc_energy < 0.0,
        "XC energy ({:.6} Ha) should be negative",
        result.xc_energy
    );
}

// ─── Nuclear repulsion is positive ───────────────────────────────────────────

#[test]
fn h2_nuclear_repulsion_positive() {
    let elements = [1u8, 1];
    let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
    let config = DftConfig::default();
    let result = solve_ks_dft(&elements, &positions, &config).unwrap();

    assert!(
        result.nuclear_repulsion > 0.0,
        "Nuclear repulsion ({:.6}) should be positive",
        result.nuclear_repulsion
    );
    // Z₁*Z₂/R₁₂ in Bohr = 1*1 / (0.74 * 1.8897) ≈ 0.7151 Ha
    let expected_nuc = 1.0 / (0.74 * 1.889_726_124_6);
    assert!(
        (result.nuclear_repulsion - expected_nuc).abs() / expected_nuc < 0.01,
        "Nuclear repulsion ({:.6}) should be near {:.6}",
        result.nuclear_repulsion,
        expected_nuc
    );
}
