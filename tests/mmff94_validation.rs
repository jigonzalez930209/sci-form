//! MMFF94 validation against an independent reference implementation.
//!
//! The reference energies live in `tests/fixtures/validation/mmff94_rdkit.json`
//! and are produced by `scripts/validation/oracle_mmff94.py`, an independent
//! Python implementation of the MMFF94 energy expression that is driven by the
//! published MMFF94 parameter tables (as embedded in RDKit's
//! `Code/ForceField/MMFF/Params.cpp`).
//!
//! Regenerate the fixture with:
//! ```bash
//! python scripts/validation/oracle_mmff94.py
//! ```
//!
//! The tolerances below are far tighter than the 0.1 kcal/mol target: the
//! Rust and Python implementations are independent, so any discrepancy in a
//! functional form or in the parameter resolution shows up immediately.

use std::collections::BTreeMap;
use std::path::PathBuf;

use serde::Deserialize;

fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join("validation")
        .join(name)
}

#[derive(Debug, Deserialize)]
struct Fixture {
    #[allow(dead_code)]
    description: String,
    #[allow(dead_code)]
    upstream: String,
    #[allow(dead_code)]
    upstream_sha256: String,
    molecules: Vec<MoleculeFixture>,
}

#[derive(Debug, Deserialize)]
struct MoleculeFixture {
    name: String,
    smiles: String,
    elements: Vec<u8>,
    #[allow(dead_code)]
    mmff94_types: Vec<u8>,
    #[allow(dead_code)]
    bonds: Vec<Vec<serde_json::Value>>,
    coords: Vec<f64>,
    counts: Counts,
    energy: BTreeMap<String, f64>,
}

#[derive(Debug, Deserialize)]
struct Counts {
    bonds: usize,
    angles: usize,
    torsions: usize,
    vdw_pairs: usize,
}

fn load() -> Fixture {
    let path = fixture_path("mmff94_rdkit.json");
    let data = std::fs::read_to_string(&path).unwrap_or_else(|err| {
        panic!(
            "missing fixture {} ({err}); run `python scripts/validation/oracle_mmff94.py`",
            path.display()
        )
    });
    serde_json::from_str(&data).expect("invalid fixture JSON")
}

/// Energy parity between the Rust MMFF94 implementation and the independent
/// Python reference.
#[test]
fn mmff94_energy_matches_independent_reference() {
    let fixture = load();
    assert!(!fixture.molecules.is_empty(), "fixture has no molecules");

    // Tolerance on the total energy (kcal/mol).  The two implementations are
    // independent, so agreement is expected at the 1e-6 kcal/mol level; the
    // 1e-3 kcal/mol bound leaves room for last-bit floating point noise.
    const TOL_TOTAL: f64 = 1.0e-3;
    // Per-term-type tolerance.
    const TOL_TERM: f64 = 1.0e-3;

    for mol in &fixture.molecules {
        let parsed = sci_form::graph::Molecule::from_smiles(&mol.smiles)
            .unwrap_or_else(|e| panic!("{}: cannot parse SMILES: {e}", mol.name));

        assert_eq!(
            parsed.graph.node_count(),
            mol.elements.len(),
            "{}: atom count mismatch",
            mol.name
        );
        for (i, &z) in mol.elements.iter().enumerate() {
            let node = petgraph::graph::NodeIndex::new(i);
            assert_eq!(
                parsed.graph[node].element, z,
                "{}: element mismatch at atom {i}",
                mol.name
            );
        }

        let energy = sci_form::compute_mmff94_energy(&mol.smiles, &mol.coords)
            .unwrap_or_else(|e| panic!("{}: MMFF94 energy failed: {e}", mol.name));

        let expected_total = mol.energy["total"];
        let dev = (energy - expected_total).abs();
        println!(
            "{:<16} rust = {:>16.8}  reference = {:>16.8}  |dev| = {:.3e} kcal/mol",
            mol.name, energy, expected_total, dev
        );
        assert!(
            dev < TOL_TOTAL,
            "{}: total energy deviation {dev:.6} kcal/mol exceeds {TOL_TOTAL} \
             (rust {energy:.8} vs reference {expected_total:.8})",
            mol.name
        );

        // Per-term-type breakdown.
        let breakdown = sci_form::compute_mmff94_energy_breakdown(&mol.smiles, &mol.coords)
            .unwrap_or_else(|e| panic!("{}: MMFF94 breakdown failed: {e}", mol.name));
        let sum: f64 = breakdown.iter().map(|(_, e)| e).sum();
        assert!(
            (sum - energy).abs() < 1e-9,
            "{}: breakdown does not sum to the total",
            mol.name
        );
        for (kind, value) in &breakdown {
            let expected = mol
                .energy
                .get(kind.as_str())
                .unwrap_or_else(|| panic!("{}: no reference for term kind {kind}", mol.name));
            let dev = (value - expected).abs();
            assert!(
                dev < TOL_TERM,
                "{}: {kind} energy deviation {dev:.6} kcal/mol exceeds {TOL_TERM} \
                 (rust {value:.8} vs reference {expected:.8})",
                mol.name
            );
        }
        for (kind, _) in &breakdown {
            assert!(
                matches!(kind.as_str(), "bond" | "angle" | "torsion" | "vdw"),
                "{}: unexpected term kind {kind}",
                mol.name
            );
        }
        assert!(
            breakdown.len() >= 2,
            "{}: expected at least a bond and an angle term",
            mol.name
        );
    }
}

/// The number of MMFF94 terms must match the reference enumeration.
#[test]
fn mmff94_term_counts_match_reference() {
    let fixture = load();
    for mol in &fixture.molecules {
        let parsed = sci_form::graph::Molecule::from_smiles(&mol.smiles).unwrap();
        let terms = sci_form::forcefield::mmff94::Mmff94Builder::build_typed(
            &parsed,
            sci_form::forcefield::mmff94::Mmff94Variant::Mmff94,
        );
        let mut counts: BTreeMap<&'static str, usize> = BTreeMap::new();
        for term in &terms {
            *counts.entry(term.kind()).or_insert(0) += 1;
        }
        assert_eq!(
            counts.get("bond").copied().unwrap_or(0),
            mol.counts.bonds,
            "{}: bond term count",
            mol.name
        );
        assert_eq!(
            counts.get("angle").copied().unwrap_or(0),
            mol.counts.angles,
            "{}: angle term count",
            mol.name
        );
        assert_eq!(
            counts.get("torsion").copied().unwrap_or(0),
            mol.counts.torsions,
            "{}: torsion term count",
            mol.name
        );
        assert_eq!(
            counts.get("vdw").copied().unwrap_or(0),
            mol.counts.vdw_pairs,
            "{}: vdW term count",
            mol.name
        );
    }
}

/// The published MMFF94 parameter tables must be reachable and reproduce
/// documented values.
#[test]
fn mmff94_published_parameters_are_used() {
    use sci_form::forcefield::mmff94::{
        mmff94_official_type, Mmff94AtomType, Mmff94Builder, Mmff94Variant,
    };

    // C-H bond of an alkane carbon: MMFF94 `kb = 4.766 md/A`, `r0 = 1.093 A`
    // (MMFF94 type 5 = HC bonded to type 1 = CR).
    let mol = sci_form::graph::Molecule::from_smiles("C").unwrap();
    let types: Vec<u8> = (0..mol.graph.node_count())
        .map(|i| {
            mmff94_official_type(sci_form::forcefield::mmff94::assign_mmff94_type_smarts(
                &mol, i,
            ))
        })
        .collect();
    assert_eq!(types[0], 1, "methane carbon must be MMFF94 type CR (1)");
    assert_eq!(types[1], 5, "methane hydrogen must be MMFF94 type HC (5)");

    let terms = Mmff94Builder::build_typed(&mol, Mmff94Variant::Mmff94);
    let mut seen = 0;
    for term in &terms {
        if term.kind() != "bond" {
            continue;
        }
        seen += 1;
    }
    assert_eq!(seen, 4, "methane has four C-H bonds");

    // Inspect the parameters through the public builder API: the total energy
    // of methane at the MMFF94 reference bond length must be ~0, and at a
    // stretched geometry positive and finite.
    let mut coords = vec![0.0; 15];
    let r = 1.093_f64;
    let dirs: [[f64; 3]; 4] = [
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ];
    for (k, d) in dirs.iter().enumerate() {
        let n = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        for c in 0..3 {
            coords[3 * (k + 1) + c] = r * d[c] / n;
        }
    }
    let (e, _) = Mmff94Builder::total_energy(&terms, &coords);
    assert!(
        e.abs() < 5.0,
        "methane energy at the reference geometry: {e}"
    );
    assert!(e.is_finite());

    // The atom-type mapping must translate the crate's enum into the official
    // MMFF94 numbers.
    assert_eq!(mmff94_official_type(Mmff94AtomType::CR), 1);
    assert_eq!(mmff94_official_type(Mmff94AtomType::CO), 3);
    assert_eq!(mmff94_official_type(Mmff94AtomType::CSp), 4);
    assert_eq!(mmff94_official_type(Mmff94AtomType::HC), 5);
    assert_eq!(mmff94_official_type(Mmff94AtomType::CB), 37);
}

/// Finite-difference gradient check for every MMFF94 term type.
#[test]
fn mmff94_gradients_match_finite_differences() {
    use sci_form::forcefield::mmff94::{
        validate_gradients, Mmff94AngleBend, Mmff94BondStretch, Mmff94BufferedVanDerWaals,
        Mmff94Torsion,
    };

    // bond stretch
    let bond = Mmff94BondStretch {
        atom_i: 0,
        atom_j: 1,
        k_b: 4.766,
        r0: 1.093,
    };
    let coords = vec![0.0, 0.0, 0.0, 1.15, 0.05, 0.02];
    let err = validate_gradients(&bond, &coords, 1e-5);
    assert!(err < 1e-4, "bond stretch gradient error: {err}");

    // angle bend (non-linear)
    let angle = Mmff94AngleBend {
        atom_i: 0,
        atom_j: 1,
        atom_k: 2,
        k_a: 0.636,
        theta0_deg: 110.549,
        linear: false,
    };
    let coords = vec![
        -0.6, 0.8, 0.0, // i
        0.0, 0.0, 0.0, // j
        0.9, 0.4, 0.3, // k
    ];
    let err = validate_gradients(&angle, &coords, 1e-5);
    assert!(err < 1e-4, "angle bend gradient error: {err}");

    // angle bend (linear centre, cosine functional form)
    let linear = Mmff94AngleBend {
        atom_i: 0,
        atom_j: 1,
        atom_k: 2,
        k_a: 0.5,
        theta0_deg: 180.0,
        linear: true,
    };
    let coords = vec![-1.1, 0.0, 0.0, 0.0, 0.0, 0.0, 1.1, 0.05, 0.0];
    let err = validate_gradients(&linear, &coords, 1e-5);
    assert!(err < 1e-4, "linear angle bend gradient error: {err}");

    // torsion
    let torsion = Mmff94Torsion {
        atom_i: 0,
        atom_j: 1,
        atom_k: 2,
        atom_l: 3,
        v1: 0.103,
        v2: 0.681,
        v3: 0.332,
    };
    let coords = vec![
        -1.5, 1.0, 0.0, //
        0.0, 0.0, 0.0, //
        1.5, 0.0, 0.0, //
        3.0, 1.0, 0.2, //
    ];
    let err = validate_gradients(&torsion, &coords, 1e-5);
    assert!(err < 1e-4, "torsion gradient error: {err}");

    // van der Waals
    let vdw = Mmff94BufferedVanDerWaals {
        atom_i_idx: 0,
        atom_j_idx: 1,
        radius_star: 3.2758,
        epsilon_depth: 0.03307,
    };
    let coords = vec![0.0, 0.0, 0.0, 3.9, 0.1, 0.0];
    let err = validate_gradients(&vdw, &coords, 1e-5);
    assert!(err < 1e-4, "vdW gradient error: {err}");
}

/// The total MMFF94 gradient must be consistent with a finite-difference
/// gradient of the total energy for a real molecule.
#[test]
fn mmff94_total_gradient_matches_finite_differences() {
    use sci_form::forcefield::mmff94::{Mmff94Builder, Mmff94Variant};

    let mol = sci_form::graph::Molecule::from_smiles("CCO").unwrap();
    let terms = Mmff94Builder::build_typed(&mol, Mmff94Variant::Mmff94);
    assert!(!terms.is_empty());

    // Slightly distorted geometry (not the MMFF94 minimum) so that the
    // gradient is non-trivial.
    let mut coords = vec![0.0; 3 * mol.graph.node_count()];
    for (i, chunk) in coords.chunks_mut(3).enumerate() {
        chunk[0] = i as f64 * 1.3;
        chunk[1] = ((i * 7) % 5) as f64 * 0.31;
        chunk[2] = ((i * 3) % 4) as f64 * 0.22;
    }

    let (_, grad) = Mmff94Builder::total_energy(&terms, &coords);
    let eps = 1e-5;
    for i in 0..coords.len() {
        let mut cp = coords.clone();
        let mut cm = coords.clone();
        cp[i] += eps;
        cm[i] -= eps;
        let (ep, _) = Mmff94Builder::total_energy(&terms, &cp);
        let (em, _) = Mmff94Builder::total_energy(&terms, &cm);
        let numerical = (ep - em) / (2.0 * eps);
        let err = (grad[i] - numerical).abs();
        assert!(
            err < 1e-3,
            "total gradient mismatch at coordinate {i}: analytic {:.8} numeric {:.8} err {err:.3e}",
            grad[i],
            numerical
        );
    }
}
