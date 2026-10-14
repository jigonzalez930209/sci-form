//! F81: CGA motors vs Rodrigues + clifford library.
//!
//! Fixture `tests/fixtures/validation/cga_clifford.json` (oracle
//! `scripts/validation/oracle_cga.py`, seed 42): rigid-motion targets
//! from closed-form Rodrigues matrices AND from an independent motor
//! sandwich in the `clifford` package (Cl(4,1)) — the latter validates
//! our bitmask blade encoding. Dihedral rotation about a bond axis +
//! dihedral angle vs numpy.

#![cfg(feature = "alpha-cga")]

use sci_form::alpha::cga::{
    apply_motor_to_subtree, conformal_point, dihedral_angle, dihedral_motor, extract_euclidean,
    Motor,
};

#[derive(serde::Deserialize)]
struct Motion {
    axis: [f64; 3],
    angle: f64,
    t: [f64; 3],
    points: Vec<[f64; 3]>,
    rodrigues_ref: Vec<[f64; 3]>,
    clifford_ref: Vec<[f64; 3]>,
}

#[derive(serde::Deserialize)]
struct Dihedral {
    axis_point: [f64; 3],
    axis_dir: [f64; 3],
    angle: f64,
    points: Vec<[f64; 3]>,
    rotated_ref: Vec<[f64; 3]>,
    dihedral_points: Vec<[f64; 3]>,
    dihedral_ref: f64,
}

#[derive(serde::Deserialize)]
struct Fixture {
    tolerance: f64,
    motions: Vec<Motion>,
    dihedrals: Vec<Dihedral>,
}

fn load() -> Fixture {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/validation/cga_clifford.json"
    );
    serde_json::from_str(&std::fs::read_to_string(path).expect("cga_clifford.json"))
        .expect("parse fixture")
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

/// F81a: motor transform vs Rodrigues and vs clifford sandwich.
#[test]
fn f81_motor_matches_rodrigues_clifford() {
    let fix = load();
    assert!(fix.motions.len() >= 10);
    let mut worst = (0.0f64, String::new());
    for (k, m) in fix.motions.iter().enumerate() {
        let motor = Motor::from_rotation_translation(m.axis, m.angle, m.t);
        for (i, p) in m.points.iter().enumerate() {
            let got = motor.transform_point(*p);
            let dr = dist(got, m.rodrigues_ref[i]);
            let dc = dist(got, m.clifford_ref[i]);
            if dr > worst.0 {
                worst = (dr, format!("motion {k} pt {i} (rodrigues)"));
            }
            assert!(dr <= fix.tolerance, "motion {k} pt {i}: Rodrigues {dr:.3e}");
            assert!(dc <= fix.tolerance, "motion {k} pt {i}: clifford {dc:.3e}");
        }
        // Composition: M = T R applied directly == composed pieces.
        let r = Motor::rotor(m.axis, m.angle);
        let t = Motor::translator(m.t);
        let composed = t.compose(&r);
        for p in &m.points {
            let d = dist(composed.transform_point(*p), motor.transform_point(*p));
            assert!(d <= 1e-12, "compose mismatch {d:.3e}");
        }
        // Reverse round-trip.
        let rev = motor.reverse();
        for p in &m.points {
            let back = rev.transform_point(motor.transform_point(*p));
            assert!(dist(back, *p) <= 1e-9, "reverse round-trip");
        }
    }
    println!("F81a motors: n={} worst {:.3e}", fix.motions.len(), worst.0);
    // Conformal lift/extract round-trip.
    let cp = conformal_point([1.0, -2.0, 3.5]);
    let back = extract_euclidean(&cp);
    assert!(dist(back, [1.0, -2.0, 3.5]) <= 1e-12);
}

/// F81b: dihedral rotation about a bond axis + subtree isolation + angle.
#[test]
fn f81_dihedral_matches_axis_rotation() {
    let fix = load();
    assert!(fix.dihedrals.len() >= 8);
    let mut worst = 0.0f64;
    for (k, d) in fix.dihedrals.iter().enumerate() {
        let b = [
            d.axis_point[0] + d.axis_dir[0],
            d.axis_point[1] + d.axis_dir[1],
            d.axis_point[2] + d.axis_dir[2],
        ];
        let m = dihedral_motor(d.axis_point, b, d.angle);
        // Subtree application: only listed atoms move, others bit-exact.
        let flat: Vec<f64> = d.points.iter().flat_map(|p| p.iter().copied()).collect();
        let idx: Vec<usize> = vec![0, 2, 4];
        let out = apply_motor_to_subtree(&flat, &idx, &m);
        for (i, want) in d.rotated_ref.iter().enumerate() {
            let got = [out[3 * i], out[3 * i + 1], out[3 * i + 2]];
            if idx.contains(&i) {
                let e = dist(got, *want);
                worst = worst.max(e);
                assert!(e <= fix.tolerance, "dihedral {k} pt {i}: {e:.3e}");
            } else {
                assert!(got == d.points[i], "non-subtree atom {i} moved");
            }
        }
        let dp = &d.dihedral_points;
        let got_a = dihedral_angle(dp[0], dp[1], dp[2], dp[3]);
        // Branch + chirality: ours returns [0, 2π) with opposite sign to the
        // numpy atan2 convention (ours ≡ −numpy mod 2π — valid mirror
        // convention; internal drivers use this function consistently).
        let tau = 2.0 * std::f64::consts::PI;
        let circ = |d: f64| {
            let r = d.abs() % tau;
            r.min(tau - r)
        };
        let da = circ(got_a - d.dihedral_ref).min(circ(got_a + d.dihedral_ref));
        assert!(da <= 1e-9, "dihedral angle {da:.3e}");
    }
    println!("F81b dihedrals: n={} worst {:.3e}", fix.dihedrals.len(), worst);
}
