//! F77: KS-DFT SVWN/PBE vs live-generated PySCF RKS/STO-3G references.
//!
//! Fixture `tests/fixtures/validation/dft_pyscf.json` (oracle
//! `scripts/validation/generate_dft_reference.py`, PySCF 2.x,
//! grids.level=4): same Å geometries, STO-3G, xc lda,vwn / pbe,pbe.
//! Gates total energy (1 mHa) and HOMO-LUMO gap (0.05 eV) on
//! H2/H2O/CH4/NH3/N2/CO — closed-shell H/C/N/O covered by our minimal basis.

#![cfg(feature = "alpha-dft")]

use sci_form::dft::grid::GridQuality;
use sci_form::dft::ks_fock::{solve_ks_dft, DftConfig, DftMethod};

#[derive(serde::Deserialize)]
struct DftCase {
    name: String,
    elements: Vec<u8>,
    #[serde(rename = "coords_A")]
    coords: Vec<[f64; 3]>,
    xc: String,
    #[serde(rename = "total_energy_Ha")]
    e_ref: f64,
    #[serde(rename = "gap_eV")]
    gap_ref: f64,
}

#[derive(serde::Deserialize)]
struct Fixture {
    #[serde(rename = "energy_tolerance_Ha")]
    energy_tolerance_ha: f64,
    #[serde(rename = "gap_tolerance_eV")]
    gap_tolerance_ev: f64,
    cases: Vec<DftCase>,
}

fn load() -> Fixture {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/validation/dft_pyscf.json"
    );
    serde_json::from_str(&std::fs::read_to_string(path).expect("dft_pyscf.json"))
        .expect("parse fixture")
}

/// F77: total energy + gap vs PySCF for both functionals.
#[test]
fn f77_dft_matches_pyscf() {
    let fix = load();
    assert!(fix.cases.len() >= 10);
    let mut worst_e = (0.0f64, String::new());
    let mut worst_g = (0.0f64, String::new());
    for c in &fix.cases {
        let method = match c.xc.as_str() {
            "svwn" => DftMethod::Svwn,
            "pbe" => DftMethod::Pbe,
            other => panic!("unknown xc {other}"),
        };
        let r = solve_ks_dft(
            &c.elements,
            &c.coords,
            &DftConfig {
                method,
                grid_quality: GridQuality::VeryFine,
                ..Default::default()
            },
        )
        .unwrap_or_else(|err| panic!("{}/{} failed: {err}", c.name, c.xc));
        assert!(r.converged, "{}/{} did not converge", c.name, c.xc);
        let de = (r.total_energy - c.e_ref).abs();
        let dg = (r.gap - c.gap_ref).abs();
        if de > worst_e.0 {
            worst_e = (de, format!("{}/{}", c.name, c.xc));
        }
        if dg > worst_g.0 {
            worst_g = (dg, format!("{}/{}", c.name, c.xc));
        }
        println!(
            "{}/{}: E={:.8} (ref {:.8}, dE={:.2e}) gap={:.3} (ref {:.3}, dg={:.2e}) iters={}",
            c.name, c.xc, r.total_energy, c.e_ref, de, r.gap, c.gap_ref, dg, r.scf_iterations
        );
    }
    println!(
        "F77 vs PySCF: n={} worst dE={:.3e} Ha ({}) worst dgap={:.3e} eV ({})",
        fix.cases.len(),
        worst_e.0,
        worst_e.1,
        worst_g.0,
        worst_g.1
    );
    assert!(
        worst_e.0 <= fix.energy_tolerance_ha,
        "dE {:.4e} Ha at {} (tol {})",
        worst_e.0,
        worst_e.1,
        fix.energy_tolerance_ha
    );
    assert!(
        worst_g.0 <= fix.gap_tolerance_ev,
        "dgap {:.4e} eV at {} (tol {})",
        worst_g.0,
        worst_g.1,
        fix.gap_tolerance_ev
    );
}
