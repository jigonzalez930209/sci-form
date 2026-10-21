//! F79: BP symmetry functions + tiny-NN energies/forces vs numpy.
//!
//! Fixture `tests/fixtures/validation/mlff_bp.json` (oracle
//! `scripts/validation/oracle_mlff.py`, RDKit ETKDG): independent numpy
//! implementation of the Behler-Parrinello formulas (radial G2 with
//! paired shift/eta lists, angular G4, cos cutoff) + numpy forward pass
//! through tiny random nets + central-FD forces of the numpy pipeline
//! (independent of the analytic backprop under test).
//! Documented gaps: TorchANI tensor comparison (species-split vs our
//! species-agnostic layout) and MD17/SPICE energy benchmarks (need
//! trained weights).

#![cfg(feature = "alpha-mlff")]

use sci_form::ml::inference::{Activation, DenseLayer, InferenceNet};
use sci_form::ml::mlff::{compute_mlff, MlffConfig};
use sci_form::ml::symmetry_functions::{compute_aevs, SymmetryFunctionParams};
use std::collections::HashMap;

#[derive(serde::Deserialize)]
struct LayerSpec {
    #[serde(rename = "in")]
    input: usize,
    #[serde(rename = "out")]
    output: usize,
    weights: Vec<f64>,
    bias: Vec<f64>,
    act: String,
}

#[derive(serde::Deserialize)]
struct MolCase {
    smiles: String,
    elements: Vec<u8>,
    coords: Vec<[f64; 3]>,
    aev_default: Vec<Vec<f64>>,
    aev_small: Vec<Vec<f64>>,
    energy_ref: f64,
    forces_ref: Vec<[f64; 3]>,
    gelu_energy_ref: Option<f64>,
}

#[derive(serde::Deserialize)]
struct Fixture {
    aev_tolerance: f64,
    energy_tolerance: f64,
    gelu_energy_tolerance: f64,
    force_tolerance: f64,
    params: HashMap<String, SymmetryFunctionParams>,
    nets: HashMap<String, Vec<LayerSpec>>,
    ani2x_paper: SymmetryFunctionParams,
    cases: Vec<MolCase>,
}

fn load() -> Fixture {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/validation/mlff_bp.json"
    );
    serde_json::from_str(&std::fs::read_to_string(path).expect("mlff_bp.json"))
        .expect("parse fixture")
}

fn activation(name: &str) -> Activation {
    match name {
        "linear" => Activation::Linear,
        "relu" => Activation::Relu,
        "celu" => Activation::Celu,
        "gelu" => Activation::Gelu,
        other => panic!("activation {other}"),
    }
}

fn net(specs: &[LayerSpec]) -> InferenceNet {
    InferenceNet::new(
        specs
            .iter()
            .map(|s| DenseLayer {
                weights: s.weights.clone(),
                bias: s.bias.clone(),
                in_features: s.input,
                out_features: s.output,
                activation: activation(&s.act),
            })
            .collect(),
    )
}

/// F79a: AEVs vs numpy BP reference (both param sets) + layout + paper constants.
#[test]
fn f79_aev_matches_bp() {
    let fix = load();
    // Defaults match the ANI-2x paper constants (literature check).
    let d = &fix.params["default"];
    assert_eq!(d.radial_shifts.len(), 16);
    assert!((d.radial_cutoff - 5.2).abs() <= 1e-15);
    assert!((d.angular_cutoff - 3.5).abs() <= 1e-15);
    assert!(d.radial_etas.iter().all(|&e| (e - 16.0).abs() <= 1e-15));
    assert!(d.angular_zetas == vec![1.0, 2.0, 4.0, 8.0]);
    let paper = &fix.ani2x_paper;
    assert_eq!(d.radial_shifts, paper.radial_shifts);
    assert_eq!(d.angular_shifts, paper.angular_shifts);

    assert!(fix.cases.len() >= 8);
    let small = &fix.params["small"];
    for c in &fix.cases {
        for (params, want) in [(d, &c.aev_default), (small, &c.aev_small)] {
            let got = compute_aevs(&c.elements, &c.coords, params);
            assert_eq!(got.len(), want.len());
            let mut worst = 0.0f64;
            for (a, w) in got.iter().zip(want) {
                // Layout: the full (eta, shift) radial grid followed by the
                // angular block. Every (eta, shift) pair contributes its own
                // symmetry function — the radial block is not a zip of the
                // two lists.
                let n_radial = a.radial.len();
                assert_eq!(
                    n_radial,
                    params.radial_shifts.len() * params.radial_etas.len()
                );
                for (g, x) in a.radial.iter().zip(&w[..n_radial]) {
                    worst = worst.max((g - x).abs());
                }
                for (g, x) in a.angular.iter().zip(&w[n_radial..]) {
                    worst = worst.max((g - x).abs());
                }
            }
            assert!(worst <= fix.aev_tolerance, "{}: {worst:.3e}", c.smiles);
        }
        // Species-agnostic by design: labels do not enter the descriptor.
        let relabeled = vec![6u8; c.elements.len()];
        let a = compute_aevs(&c.elements, &c.coords, d);
        let b = compute_aevs(&relabeled, &c.coords, d);
        for (x, y) in a.iter().zip(&b) {
            assert_eq!(x.to_vec(), y.to_vec());
        }
    }
    println!(
        "F79a AEV vs numpy: {} molecules x 2 param sets",
        fix.cases.len()
    );
}

/// F79b: MLFF energies vs numpy forward; forces vs numpy FD.
#[test]
fn f79_mlff_energy_forces_match() {
    let fix = load();
    let params = &fix.params["default"];
    let mut worst_e = 0.0f64;
    let mut worst_f = 0.0f64;
    for c in &fix.cases {
        let mut element_nets = HashMap::new();
        for z in c
            .elements
            .iter()
            .copied()
            .collect::<std::collections::HashSet<_>>()
        {
            let key = z.to_string();
            let spec = fix
                .nets
                .get(&key)
                .unwrap_or_else(|| panic!("no net for {z}"));
            element_nets.insert(z, net(spec));
        }
        let cfg = MlffConfig {
            aev_params: params.clone(),
            element_nets,
        };
        let r = compute_mlff(&c.elements, &c.coords, &cfg).expect("mlff");
        let de = (r.energy - c.energy_ref).abs() / c.energy_ref.abs().max(1e-300);
        worst_e = worst_e.max(de);
        assert!(de <= fix.energy_tolerance, "{}: {de:.3e}", c.smiles);
        for (g, w) in r.forces.iter().zip(&c.forces_ref) {
            for d in 0..3 {
                worst_f = worst_f.max((g[d] - w[d]).abs());
            }
        }
    }
    println!("F79b energies worst {worst_e:.3e}, forces worst {worst_f:.3e}");
    assert!(worst_f <= fix.force_tolerance, "forces {worst_f:.3e}");
    // GELU-net energies (erf-approximation limited tolerance): zero nets
    // for non-C elements so the total equals the C-only reference.
    let params = &fix.params["default"];
    let dim = params.radial_shifts.len() * params.radial_etas.len()
        + params.angular_shifts.len() * params.angular_etas.len() * params.angular_zetas.len();
    let zero_net = || {
        InferenceNet::new(vec![DenseLayer {
            weights: vec![0.0; dim],
            bias: vec![0.0],
            in_features: dim,
            out_features: 1,
            activation: Activation::Linear,
        }])
    };
    let mut worst_g = 0.0f64;
    for c in fix.cases.iter().filter(|c| c.gelu_energy_ref.is_some()) {
        let mut element_nets = HashMap::new();
        element_nets.insert(6u8, net(&fix.nets["6-gelu"]));
        for z in c
            .elements
            .iter()
            .copied()
            .collect::<std::collections::HashSet<_>>()
        {
            if z != 6 {
                element_nets.insert(z, zero_net());
            }
        }
        let cfg = MlffConfig {
            aev_params: params.clone(),
            element_nets,
        };
        let r = compute_mlff(&c.elements, &c.coords, &cfg).expect("mlff gelu");
        let want = c.gelu_energy_ref.unwrap();
        let dg = (r.energy - want).abs() / want.abs().max(1e-300);
        worst_g = worst_g.max(dg);
        assert!(dg <= fix.gelu_energy_tolerance, "{}: {dg:.3e}", c.smiles);
    }
    println!("F79b GELU energies worst {worst_g:.3e}");
}
