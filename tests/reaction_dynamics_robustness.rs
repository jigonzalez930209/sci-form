//! Regression tests for the reaction-dynamics robustness blockers found by
//! the v1.0 release gate.
//!
//! `tests/reaction_dynamics_study.rs` had never run in CI (the packaging
//! audit found 9 test targets that CI skipped), and one of them failed: the
//! hydration of acetaldehyde put a water oxygen 0.14 Å away from the carbonyl
//! carbon, the 1/r repulsion terms diverged, the analytic gradients overflowed
//! to `NaN`, and every remaining frame of the path was `NaN`.
//!
//! The fix has three parts:
//!
//! * a minimum inter-fragment separation, so two atoms of different
//!   fragments can never end up effectively on top of each other,
//! * a displacement cap and a descent-direction check in the constrained
//!   relaxation, so a large gradient cannot throw an atom across the
//!   molecule,
//! * a refusal to run the CI-NEB on a singular geometry, instead of
//!   propagating a `NaN` into every image.

#![cfg(feature = "alpha-reaction-dynamics")]

use sci_form::alpha::reaction_dynamics::constrained_opt::{
    enforce_min_interfragment_distance, MIN_INTERFRAGMENT_DISTANCE,
};
use sci_form::alpha::reaction_dynamics::{compute_reaction_dynamics_3d, ReactionDynamics3DConfig};

/// The configuration used by the reaction battery in
/// `tests/reaction_dynamics_study.rs`.
fn hydration_config() -> ReactionDynamics3DConfig {
    ReactionDynamics3DConfig {
        method: "xtb".to_string(),
        n_images: 7,
        neb_max_iter: 30,
        spring_k: 0.1,
        step_size: 0.005,
        use_climbing_image: true,
        ci_neb_force_threshold: 0.5,
        n_approach_frames: 9,
        n_departure_frames: 9,
        far_distance: 6.0,
        reactive_distance: 2.0,
        seed: 42,
        use_orbital_guidance: false,
        use_electrostatic_steering: false,
        optimise_complex: false,
        complex_opt_max_steps: 10,
        compute_properties: false,
        n_angular_samples: 0,
        smirks: None,
    }
}

#[test]
fn v1_reaction_dynamics_never_emits_a_non_finite_frame() {
    // Acetaldehyde + water -> geminal diol. The product-guided alignment of
    // this pair is the one that used to collide.
    let result = compute_reaction_dynamics_3d(&["CC=O", "O"], &["CC(O)O"], &hydration_config())
        .expect("the hydration path must be computable");

    assert!(
        !result.frames.is_empty(),
        "the reaction path must contain at least one frame"
    );

    for frame in &result.frames {
        assert!(
            frame.energy_kcal_mol.is_finite(),
            "frame {} ({}) has a non-finite energy: {}",
            frame.index,
            frame.phase,
            frame.energy_kcal_mol
        );
        assert!(
            frame.coords.iter().all(|c| c.is_finite()),
            "frame {} ({}) has a non-finite coordinate",
            frame.index,
            frame.phase
        );
    }
}

#[test]
fn v1_interfragment_separation_guard_pushes_clashing_fragments_apart() {
    // Two fragments with an atom pair 0.2 Å apart.
    let mut coords: Vec<f64> = vec![0.0, 0.0, 0.0, 0.2, 0.0, 0.0];
    let ranges = [(0usize, 1usize), (1usize, 2usize)];

    assert!(
        enforce_min_interfragment_distance(&mut coords, &ranges, MIN_INTERFRAGMENT_DISTANCE),
        "the guard must report that it changed the geometry"
    );

    let d = (coords[3] - coords[0]).powi(2).sqrt();
    assert!(
        d >= MIN_INTERFRAGMENT_DISTANCE - 1e-9,
        "the fragments must end up at least {MIN_INTERFRAGMENT_DISTANCE} Å apart, got {d}"
    );

    // A geometry that already satisfies the constraint is left alone.
    let mut ok: Vec<f64> = vec![0.0, 0.0, 0.0, 5.0, 0.0, 0.0];
    assert!(
        !enforce_min_interfragment_distance(&mut ok, &ranges, MIN_INTERFRAGMENT_DISTANCE),
        "the guard must not touch a geometry that already satisfies the constraint"
    );
    assert_eq!(ok, vec![0.0, 0.0, 0.0, 5.0, 0.0, 0.0]);
}
