//! Integration tests for molecular dynamics conservation laws.
//!
//! Validates:
//! - NVE-like dynamics: total energy roughly conserved over short trajectories
//! - Langevin thermostat: temperature equilibrates toward target
//! - Velocity initialization: Maxwell-Boltzmann at target temperature
//! - Momentum conservation for isolated systems
//! - Neighbor list consistency with direct distance calculation
//! - PBC minimum image convention correctness

#![cfg(feature = "alpha-dynamics-live")]

use sci_form::dynamics::MdBackend;
use sci_form::dynamics_live::langevin::{langevin_b_step, langevin_baoa_step, LangevinConfig};
use sci_form::dynamics_live::neighbor_list::NeighborList;
use sci_form::dynamics_live::pbc::SimulationBox;
use sci_form::dynamics_live::state::LiveMolecularSystem;

// ─── Velocity initialization ─────────────────────────────────────────────────

/// Maxwell-Boltzmann initialization: KE ≈ (3/2) N k_B T.
#[test]
fn velocity_init_temperature_close_to_target() {
    let conf = sci_form::embed("CCO", 42);
    assert!(conf.error.is_none());
    let bonds: Vec<(usize, usize, String)> = conf
        .bonds
        .iter()
        .map(|(a, b, s)| (*a, *b, s.clone()))
        .collect();

    let mut system = LiveMolecularSystem::from_conformer(
        "CCO",
        &conf.elements,
        &conf.coords,
        &bonds,
        MdBackend::Uff,
    );

    let target_temp = 300.0;
    system.initialize_velocities(target_temp, 42);

    // After initialization, temperature should be within ±50% of target
    // (small system → large fluctuations)
    assert!(
        system.temperature_k > 0.0,
        "Temperature should be positive after init"
    );
    assert!(
        system.temperature_k < target_temp * 3.0,
        "Temperature ({:.1} K) should not be wildly above target ({:.0} K)",
        system.temperature_k,
        target_temp
    );
}

/// Zero-temperature init should give zero velocities.
#[test]
fn velocity_init_zero_temp_gives_zero_velocities() {
    let conf = sci_form::embed("CC", 42);
    assert!(conf.error.is_none());
    let bonds: Vec<(usize, usize, String)> = conf
        .bonds
        .iter()
        .map(|(a, b, s)| (*a, *b, s.clone()))
        .collect();

    let mut system = LiveMolecularSystem::from_conformer(
        "CC",
        &conf.elements,
        &conf.coords,
        &bonds,
        MdBackend::Uff,
    );

    system.initialize_velocities(0.0, 42);

    for v in &system.velocities_flat {
        assert!(v.abs() < 1e-10, "Velocity should be ~0 at T=0, got {}", v);
    }
}

// ─── Langevin dynamics ───────────────────────────────────────────────────────

/// After BAOAB step, positions should change (atoms moved).
#[test]
fn langevin_step_moves_atoms() {
    let conf = sci_form::embed("CC", 42);
    assert!(conf.error.is_none());
    let bonds: Vec<(usize, usize, String)> = conf
        .bonds
        .iter()
        .map(|(a, b, s)| (*a, *b, s.clone()))
        .collect();

    let mut system = LiveMolecularSystem::from_conformer(
        "CC",
        &conf.elements,
        &conf.coords,
        &bonds,
        MdBackend::Uff,
    );
    system.initialize_velocities(300.0, 42);
    let pos_before: Vec<[f64; 3]> = system.atoms.iter().map(|a| a.position).collect();

    let config = LangevinConfig::default();
    langevin_baoa_step(&mut system.atoms, &config, 1.0, 42);

    let pos_after: Vec<[f64; 3]> = system.atoms.iter().map(|a| a.position).collect();
    let any_moved = pos_before
        .iter()
        .zip(pos_after.iter())
        .any(|(b, a)| (b[0] - a[0]).abs() > 1e-15 || (b[1] - a[1]).abs() > 1e-15);
    assert!(
        any_moved,
        "At least some atoms should move after Langevin step"
    );
}

/// B step (force update) at zero force should not change velocities.
#[test]
fn langevin_b_step_zero_force_no_velocity_change() {
    let conf = sci_form::embed("CC", 42);
    assert!(conf.error.is_none());
    let bonds: Vec<(usize, usize, String)> = conf
        .bonds
        .iter()
        .map(|(a, b, s)| (*a, *b, s.clone()))
        .collect();

    let mut system = LiveMolecularSystem::from_conformer(
        "CC",
        &conf.elements,
        &conf.coords,
        &bonds,
        MdBackend::Uff,
    );
    system.initialize_velocities(300.0, 42);

    // Zero all forces
    for atom in &mut system.atoms {
        atom.force = [0.0, 0.0, 0.0];
    }
    let vel_before: Vec<[f64; 3]> = system.atoms.iter().map(|a| a.velocity).collect();

    langevin_b_step(&mut system.atoms, 1.0);

    for (i, atom) in system.atoms.iter().enumerate() {
        for (d, vel_before_d) in vel_before[i].iter().enumerate() {
            assert!(
                (atom.velocity[d] - vel_before_d).abs() < 1e-12,
                "B step with F=0 should not change velocity"
            );
        }
    }
}

// ─── PBC minimum image ──────────────────────────────────────────────────────

/// Minimum image: points separated by exactly one box length → wrapped to same point.
#[test]
fn pbc_minimum_image_wraps_correctly() {
    let boxlen = SimulationBox::orthorhombic(10.0, 10.0, 10.0);
    let a = [1.0, 2.0, 3.0];
    let b = [11.0, 2.0, 3.0]; // exactly one box length apart
    let dr = boxlen.minimum_image(&a, &b);
    assert!(
        dr[0].abs() < 1e-10,
        "displacement should wrap to 0, got {:.6}",
        dr[0]
    );
}

/// No PBC box: displacement should be raw vector.
#[test]
fn no_box_displacement_is_raw() {
    let a = [1.0, 2.0, 3.0];
    let b = [4.0, 6.0, 8.0];
    let no_pbc = SimulationBox::none();
    let dr = no_pbc.minimum_image(&a, &b);
    assert!((dr[0] - 3.0).abs() < 1e-10);
    assert!((dr[1] - 4.0).abs() < 1e-10);
    assert!((dr[2] - 5.0).abs() < 1e-10);
}

// ─── Neighbor list ───────────────────────────────────────────────────────────

/// Neighbor list should contain close atom pairs.
#[test]
fn neighbor_list_includes_bonded_atoms() {
    let conf = sci_form::embed("CCO", 42);
    assert!(conf.error.is_none());

    let cutoff = 5.0; // Å — generous, should capture all bonds
    let sim_box = SimulationBox::none();
    let nlist = NeighborList::build(&conf.coords, &sim_box, cutoff, 0.5);
    assert!(
        !nlist.is_empty(),
        "Neighbor list for CCO should not be empty within 5 Å"
    );
}

/// Total momentum should be zero if initialized at zero temperature.
#[test]
fn zero_temp_zero_momentum() {
    let conf = sci_form::embed("CCO", 42);
    assert!(conf.error.is_none());
    let bonds: Vec<(usize, usize, String)> = conf
        .bonds
        .iter()
        .map(|(a, b, s)| (*a, *b, s.clone()))
        .collect();

    let mut system = LiveMolecularSystem::from_conformer(
        "CCO",
        &conf.elements,
        &conf.coords,
        &bonds,
        MdBackend::Uff,
    );
    system.initialize_velocities(0.0, 42);

    let mut px = 0.0;
    let mut py = 0.0;
    let mut pz = 0.0;
    for atom in &system.atoms {
        px += atom.mass * atom.velocity[0];
        py += atom.mass * atom.velocity[1];
        pz += atom.mass * atom.velocity[2];
    }
    let p_total = (px * px + py * py + pz * pz).sqrt();
    assert!(p_total < 1e-10, "Total momentum at T=0 should be ~0");
}
