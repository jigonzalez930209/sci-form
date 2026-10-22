//! F84: live-dynamics UFF force path vs RDKit + Verlet properties.
//!
//! Fixture `tests/fixtures/validation/uff_forces_rdkit.json` (oracle
//! `scripts/validation/oracle_dynamics.py`, F16 construction): RDKit UFF
//! `-CalcGrad` forces. F16 gates our UFF energy to 0.0; this gates the
//! force vector actually integrated by velocity Verlet, plus NVE
//! conservation, momentum conservation, reversibility, and the
//! Berendsen thermostat approach. (LAMMPS has no UFF; trajectory-vs-
//! engine comparison is documented out of reach.)

#![cfg(feature = "alpha-dynamics-live")]

use sci_form::dynamics::MdBackend;
use sci_form::dynamics_live::state::LiveMolecularSystem;

#[derive(serde::Deserialize)]
struct DynCase {
    smiles: String,
    elements: Vec<u8>,
    coords: Vec<f64>,
    reference_energy: f64,
    rdkit_forces: Vec<f64>,
}

#[derive(serde::Deserialize)]
struct Fixture {
    force_tolerance: f64,
    cases: Vec<DynCase>,
}

fn load() -> Fixture {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/validation/uff_forces_rdkit.json"
    );
    serde_json::from_str(&std::fs::read_to_string(path).expect("uff_forces_rdkit.json"))
        .expect("parse fixture")
}

fn topology(smiles: &str) -> (Vec<u8>, Vec<(usize, usize, String)>) {
    let mol = sci_form::parse(smiles).expect("parse");
    let n = mol.graph.node_count();
    let elements: Vec<u8> = (0..n)
        .map(|i| mol.graph[petgraph::graph::NodeIndex::new(i)].element)
        .collect();
    let bonds: Vec<(usize, usize, String)> = mol
        .graph
        .edge_indices()
        .map(|e| {
            let (a, b) = mol.graph.edge_endpoints(e).unwrap();
            let order = match mol.graph[e].order {
                sci_form::graph::BondOrder::Single => "SINGLE",
                sci_form::graph::BondOrder::Double => "DOUBLE",
                sci_form::graph::BondOrder::Triple => "TRIPLE",
                sci_form::graph::BondOrder::Aromatic => "AROMATIC",
                sci_form::graph::BondOrder::Unknown => "UNKNOWN",
                // `BondOrder` is `#[non_exhaustive]`: a bond order added after
                // this test was written reports itself as unknown.
                _ => "UNKNOWN",
            };
            (a.index(), b.index(), order.to_string())
        })
        .collect();
    (elements, bonds)
}

fn system_for(c: &DynCase) -> LiveMolecularSystem {
    let (elements, bonds) = topology(&c.smiles);
    assert_eq!(
        elements.len() * 3,
        c.coords.len(),
        "{}: topology/coords",
        c.smiles
    );
    assert_eq!(elements, c.elements, "{}: elements", c.smiles);
    LiveMolecularSystem::from_conformer(&c.smiles, &elements, &c.coords, &bonds, MdBackend::Uff)
}

/// F84a: UFF forces vs RDKit (-gradient), energy vs F16 reference.
#[test]
fn f84_forces_match_rdkit() {
    let fix = load();
    assert!(fix.cases.len() >= 10);
    let mut worst = (0.0f64, String::new());
    for c in &fix.cases {
        let mut sys = system_for(c);
        let e = sys.compute_forces().expect("forces");
        let de = (e - c.reference_energy).abs();
        assert!(de <= 0.01, "{}: energy {de:.3e}", c.smiles);
        let df = sys
            .forces_flat
            .iter()
            .zip(&c.rdkit_forces)
            .fold(0.0f64, |m, (&g, &w)| m.max((g - w).abs()));
        if df > worst.0 {
            worst = (df, c.smiles.clone());
        }
        assert!(df <= fix.force_tolerance, "{}: forces {df:.3e}", c.smiles);
    }
    println!(
        "F84a forces vs RDKit: n={} worst {:.3e}",
        fix.cases.len(),
        worst.0
    );
}

/// Total linear momentum (amu·Å/fs).
fn momentum(sys: &LiveMolecularSystem) -> [f64; 3] {
    let mut p = [0.0; 3];
    for a in &sys.atoms {
        for (d, p_d) in p.iter_mut().enumerate() {
            *p_d += a.mass * a.velocity[d];
        }
    }
    p
}

/// F84b: NVE Verlet — energy drift, momentum conservation, reversibility.
#[test]
fn f84_nve_verlet_properties() {
    let fix = load();
    let c = &fix.cases[0];
    let mut sys = system_for(c);
    let _pe0 = sys.compute_forces().expect("forces");
    sys.initialize_velocities(300.0, 42);
    // Zero total momentum for a clean conservation check.
    let p0 = momentum(&sys);
    let mtot: f64 = sys.atoms.iter().map(|a| a.mass).sum();
    for a in sys.atoms.iter_mut() {
        for (d, v) in a.velocity.iter_mut().enumerate() {
            *v -= p0[d] / mtot;
        }
    }
    sys.sync_flat_from_atoms();
    sys.potential_energy = sys.compute_forces().expect("forces");
    let e0 = sys.potential_energy + sys.kinetic_energy_and_temperature().0;
    sys.integrate(0.25, 2000).expect("nve");
    let e1 = sys.potential_energy + sys.kinetic_energy;
    let drift = (e1 - e0).abs() / e0.abs().max(1e-300);
    let p1 = momentum(&sys);
    let pmag: f64 = p1.iter().map(|v| v * v).sum::<f64>().sqrt();
    println!(
        "F84b NVE 1ps: drift {drift:.3e} |p| {pmag:.3e} t={}fs",
        sys.time_fs
    );
    assert!(drift <= 1e-3, "energy drift {drift:.3e}");
    assert!(pmag <= 1e-6, "momentum {pmag:.3e}");
    // Reversibility over a short horizon (long horizons diverge by chaos
    // amplification of fp noise, e.g. parallel reduction order): fresh
    // system, 50 steps forward, flip v, 50 back, recover start.
    let mut sys2 = system_for(c);
    sys2.compute_forces().expect("forces");
    sys2.initialize_velocities(300.0, 42);
    let start = sys2.positions_flat.clone();
    sys2.integrate(0.25, 50).expect("nve fwd");
    for a in sys2.atoms.iter_mut() {
        for d in 0..3 {
            a.velocity[d] = -a.velocity[d];
        }
    }
    sys2.sync_flat_from_atoms();
    sys2.integrate(0.25, 50).expect("nve back");
    let back: f64 = sys2
        .positions_flat
        .iter()
        .zip(&start)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt()
        / (start.len() as f64).sqrt();
    println!("F84b reversibility RMS {back:.3e} A");
    assert!(back <= 1e-6, "not reversible {back:.3e}");
}

/// F84c: Berendsen thermostat drives temperature to target.
#[test]
fn f84_berendsen_reaches_target() {
    let fix = load();
    let c = &fix.cases[1];
    let mut sys = system_for(c);
    sys.compute_forces().expect("forces");
    sys.initialize_velocities(300.0, 7);
    // Small-N instantaneous T fluctuates ~+-90K; gate the time average.
    let mut acc = 0.0;
    let navg = 200;
    for i in 0..600 {
        sys.verlet_step(0.5).expect("step");
        sys.berendsen_thermostat(300.0, 50.0, 0.5);
        if i >= 600 - navg {
            acc += sys.temperature_k;
        }
    }
    let tavg = acc / navg as f64;
    println!("F84c Berendsen: <T>={tavg:.2}K target=300K");
    assert!((tavg - 300.0).abs() <= 40.0);
}
