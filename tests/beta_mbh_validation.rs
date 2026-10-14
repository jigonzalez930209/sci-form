//! F92: MBH harmonic frequencies vs SI closed form + block structure.
//!
//! Fixture `tests/fixtures/validation/mbh_analytic.json` (oracle
//! `scripts/validation/oracle_mbh.py`): diatomic stretch from SI
//! constants — independent of the 108.59 conversion factor in the code.
//! Expected spectrum: 5 near-zero + 1 stretch. Block test: benzene
//! hexagon collapses to 1 rigid block with reduced DOF.
//! Documented gap: experimental +-2 cm-1 (Ghysels) unmet.

#![cfg(feature = "beta-mbh")]

use sci_form::beta::mbh::{compute_mbh_frequencies, detect_rigid_blocks, MbhConfig};

#[derive(serde::Deserialize)]
struct MbhCase {
    name: String,
    m1_amu: f64,
    m2_amu: f64,
    k_kcal_per_mol_A2: f64,
    r0_A: f64,
    nu_ref_cm: f64,
}

#[derive(serde::Deserialize)]
struct Fixture {
    stretch_tolerance_rel: f64,
    zero_tolerance_cm: f64,
    cases: Vec<MbhCase>,
}

fn load() -> Fixture {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/validation/mbh_analytic.json"
    );
    serde_json::from_str(&std::fs::read_to_string(path).expect("mbh_analytic.json"))
        .expect("parse fixture")
}

fn amu_of(m: f64) -> u8 {
    // Map oracle masses back to elements with (near-)exact masses.
    if (m - 1.008).abs() < 1e-6 {
        1
    } else if (m - 14.007).abs() < 1e-6 {
        7
    } else if (m - 12.011).abs() < 1e-6 {
        6
    } else if (m - 15.999).abs() < 1e-6 {
        8
    } else if (m - 35.453).abs() < 1e-6 {
        17
    } else {
        panic!("no element for mass {m}")
    }
}

/// F92a: diatomic stretch vs SI closed form; 5 zero modes.
#[test]
fn f92_diatomic_matches_analytic() {
    let fix = load();
    assert!(fix.cases.len() >= 4);
    for c in &fix.cases {
        let elements = [amu_of(c.m1_amu), amu_of(c.m2_amu)];
        let positions = [[0.0, 0.0, 0.0], [c.r0_A, 0.0, 0.0]];
        let k = c.k_kcal_per_mol_A2;
        let r0 = c.r0_A;
        let energy = move |coords: &[f64]| {
            let dx = coords[3] - coords[0];
            let dy = coords[4] - coords[1];
            let dz = coords[5] - coords[2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            0.5 * k * (r - r0) * (r - r0)
        };
        let cfg = MbhConfig { fd_step: 0.001 };
        let r = compute_mbh_frequencies(&elements, &positions, &[], &energy, &cfg);
        assert_eq!(r.frequencies.len(), 6, "{}: {}", c.name, r.frequencies.len());
        let stretch = r.frequencies[5];
        let d = (stretch - c.nu_ref_cm).abs() / c.nu_ref_cm;
        let zero_max = r.frequencies[..5].iter().fold(0.0f64, |m, &v| m.max(v.abs()));
        println!("F92a {}: stretch {stretch:.4} (ref {:.4}, d={d:.2e}) zero_max {zero_max:.2e}",
            c.name, c.nu_ref_cm);
        assert!(d <= fix.stretch_tolerance_rel, "{}: {d:.3e}", c.name);
        assert!(zero_max <= fix.zero_tolerance_cm, "{}: {zero_max:.3e}", c.name);
    }
}

/// F92b: benzene hexagon -> exactly 1 rigid block, reduced DOF, speedup.
#[test]
fn f92_benzene_block_structure() {
    // Ideal hexagon: C-C 1.39 A, C-H 1.09 A, planar.
    let cc = 1.39f64;
    let ch = 1.09f64;
    let mut elements = vec![6u8; 6];
    let mut positions = Vec::new();
    for i in 0..6 {
        let a = i as f64 * std::f64::consts::PI / 3.0;
        positions.push([cc * a.cos(), cc * a.sin(), 0.0]);
    }
    for i in 0..6 {
        elements.push(1);
        let a = i as f64 * std::f64::consts::PI / 3.0;
        positions.push([(cc + ch) * a.cos(), (cc + ch) * a.sin(), 0.0]);
    }
    let ring: Vec<usize> = (0..6).collect();
    let decomp = detect_rigid_blocks(12, &elements, &positions, &[(ring, true)]);
    assert_eq!(decomp.blocks.len(), 1, "benzene: 1 aromatic block");
    assert_eq!(decomp.blocks[0].atom_indices.len(), 6);
    assert_eq!(decomp.flexible_atoms.len(), 6, "6 H stay flexible");
    // Harmonic bond energy for a smoke frequency run.
    let energy = |coords: &[f64]| {
        let mut e = 0.0;
        for i in 0..12 {
            for j in (i + 1)..12 {
                let dx = coords[3 * j] - coords[3 * i];
                let dy = coords[3 * j + 1] - coords[3 * i + 1];
                let dz = coords[3 * j + 2] - coords[3 * i + 2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                let r0 = if elements[i] == 6 && elements[j] == 6 { 1.39 } else { 1.09 };
                e += 50.0 * (r - r0) * (r - r0);
            }
        }
        e
    };
    let r = compute_mbh_frequencies(&elements, &positions, &[( (0..6).collect(), true)], &energy, &MbhConfig::default());
    println!(
        "F92b benzene: blocks={} flex={} dof {}/{} speedup {:.2} nf={}",
        r.n_blocks, r.n_flexible, r.n_dof_reduced, r.n_dof_full, r.speedup, r.frequencies.len()
    );
    assert_eq!(r.n_blocks, 1);
    assert!(r.n_dof_reduced < r.n_dof_full);
    assert!((r.speedup - r.n_dof_full as f64 / r.n_dof_reduced as f64).abs() < 1e-9);
    assert!(r.frequencies.iter().all(|f| f.is_finite()));
}
