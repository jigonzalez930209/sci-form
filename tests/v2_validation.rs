//! F101–F115: validation gates for the v2.x candidates V2-01…V2-15.
//!
//! Each ficha writes `target/validation/F1xx.json`; `scripts/validation/collect.py`
//! merges the rows into `reports/validation/summary.json`.
//!
//! Registration rule (`docs/v2-roadmap.md` §8):
//! - gates the prototype must meet (L0–L2) are asserted and recorded with
//!   `pass: true`;
//! - gates it does not meet yet (L3/L4) record the measured value with
//!   `pass: null` and `status: "gap"`. A pass is never invented.
//!
//! Oracles: `tests/fixtures/validation/v2_pyscf.json`
//! (`scripts/validation/oracle_pyscf_v2.py`, PySCF, STO-3G cartesian, all
//! electrons) and `tests/fixtures/validation/v2_pymbar.json`
//! (`scripts/validation/oracle_pymbar_v2.py`).
//!
//! Run with `cargo test --release --features alpha-v2 --test v2_validation`.

#![cfg(feature = "alpha-v2")]

use std::path::PathBuf;

use serde_json::{json, Value};

const EV_TO_KCAL: f64 = 23.060_548_867;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn fixture(name: &str) -> Value {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/validation")
        .join(name);
    let text = std::fs::read_to_string(&path).unwrap_or_else(|err| {
        panic!(
            "missing fixture {} ({err}); run the matching oracle in scripts/validation/",
            path.display()
        )
    });
    serde_json::from_str(&text).expect("invalid fixture JSON")
}

struct Molecule {
    name: String,
    elements: Vec<u8>,
    positions: Vec<[f64; 3]>,
}

fn molecule(record: &Value) -> Molecule {
    let source = record.get("molecule").unwrap_or(record);
    Molecule {
        name: record["name"].as_str().unwrap_or("").to_string(),
        elements: serde_json::from_value(source["elements"].clone()).unwrap(),
        positions: serde_json::from_value(source["positions_angstrom"].clone()).unwrap(),
    }
}

/// One validation row. `pass: None` marks a gap (measured, not met yet).
struct Row<'a> {
    ficha: &'a str,
    observable: &'a str,
    metric: &'a str,
    value: f64,
    tolerance: f64,
    pass: Option<bool>,
    n_cases: usize,
    level: &'a str,
    units: &'a str,
    reference: &'a str,
    notes: &'a str,
}

fn record(row: Row) {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("target")
        .join("validation");
    let _ = std::fs::create_dir_all(&dir);
    let mut payload = json!({
        "ficha": row.ficha,
        "observable": row.observable,
        "metric": row.metric,
        "value": row.value,
        "tolerance": row.tolerance,
        "pass": row.pass,
        "n_cases": row.n_cases,
        "source": "tests/v2_validation.rs",
        "level": row.level,
        "units": row.units,
        "reference": { "tool": row.reference },
    });
    if row.pass.is_none() {
        payload["status"] = json!("gap");
    }
    if !row.notes.is_empty() {
        payload["notes"] = json!(row.notes);
    }
    eprintln!(
        "{} {:<44} {} = {:.6e} (tol {:.1e}) {:?}",
        row.ficha, row.observable, row.metric, row.value, row.tolerance, row.pass
    );
    let path = dir.join(format!("{}.json", row.ficha));
    let mut rows: Vec<Value> = std::fs::read_to_string(&path)
        .ok()
        .and_then(|text| serde_json::from_str::<Value>(&text).ok())
        .map(|value| match value {
            Value::Array(items) => items,
            other => vec![other],
        })
        .unwrap_or_default();
    rows.retain(|r| r.get("observable").and_then(|v| v.as_str()) != Some(row.observable));
    rows.push(payload);
    let _ = std::fs::write(&path, Value::Array(rows).to_string());
}

/// Record a gate that must hold, then assert it.
fn gate(row: Row) {
    let (ficha, observable, value, tolerance) =
        (row.ficha, row.observable, row.value, row.tolerance);
    let passed = value.is_finite() && value <= tolerance;
    record(Row {
        pass: Some(passed),
        ..row
    });
    assert!(
        passed,
        "{ficha} {observable}: {value:.6e} exceeds the tolerance {tolerance:.1e}"
    );
}

/// Record a gate the prototype does not meet yet: the value, never a pass.
fn gap(row: Row) {
    record(Row { pass: None, ..row });
}

fn max_abs(values: impl IntoIterator<Item = f64>) -> f64 {
    values.into_iter().fold(0.0, |acc, v| acc.max(v.abs()))
}

/// Deterministic standard normal draws (SplitMix64 + Box–Muller).
struct Normal(u64);

impl Normal {
    fn uniform(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        ((z >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }

    fn sample(&mut self) -> f64 {
        let (u1, u2) = (self.uniform(), self.uniform());
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

// ---------------------------------------------------------------------------
// F101 (V2-01) PBE0 hybrid DFT vs PySCF
// ---------------------------------------------------------------------------

#[test]
fn f101_pbe0_matches_pyscf() {
    use sci_form::dft::grid::GridQuality;
    use sci_form::dft::{solve_ks_dft, DftConfig, DftMethod};

    let data = fixture("v2_pyscf.json");
    let cases = data["dft"].as_array().unwrap();
    let mut pbe_err = Vec::new();
    let mut pbe0_err = Vec::new();
    for case in cases {
        let mol = molecule(case);
        for (method, key, errors) in [
            (DftMethod::Pbe, "pbe_hartree", &mut pbe_err),
            (DftMethod::Pbe0, "pbe0_hartree", &mut pbe0_err),
        ] {
            let config = DftConfig {
                method,
                grid_quality: GridQuality::VeryFine,
                energy_threshold: 1e-9,
                density_threshold: 1e-7,
                max_iterations: 200,
                ..DftConfig::default()
            };
            let result = solve_ks_dft(&mol.elements, &mol.positions, &config).unwrap();
            assert!(result.converged, "{} {key} did not converge", mol.name);
            let error = result.total_energy - case[key].as_f64().unwrap();
            eprintln!("F101 {:<5} {key:<13} error {error:+.6} Ha", mol.name);
            errors.push(error);
        }
    }
    let n = cases.len();
    let reference = "PySCF RKS/STO-3G, grids.level 4";
    let pbe = max_abs(pbe_err);
    let pbe0 = max_abs(pbe0_err);
    let row = |observable, value| Row {
        ficha: "F101",
        observable,
        metric: "max_abs_error_hartree",
        value,
        tolerance: 1e-3,
        pass: None,
        n_cases: n,
        level: "L2",
        units: "Ha",
        reference,
        notes: "",
    };
    gate(row("pbe_total_energy", pbe));
    gate(row("pbe0_total_energy", pbe0));
}

// ---------------------------------------------------------------------------
// F102 (V2-02) MP2 correlation energy vs PySCF
// ---------------------------------------------------------------------------

#[test]
fn f102_mp2_matches_pyscf() {
    let data = fixture("v2_pyscf.json");
    let cases = data["mp2"].as_array().unwrap();
    let mut hf = Vec::new();
    let mut corr = Vec::new();
    for case in cases {
        let mol = molecule(case);
        let result = sci_form::compute_mp2_energy(&mol.elements, &mol.positions).unwrap();
        hf.push(result.hf_energy_hartree - case["hf_energy_hartree"].as_f64().unwrap());
        corr.push(
            result.correlation_energy_hartree
                - case["correlation_energy_hartree"].as_f64().unwrap(),
        );
    }
    let reference = "PySCF RHF+MP2/STO-3G, all electrons";
    let row = |observable, value, tolerance| Row {
        ficha: "F102",
        observable,
        metric: "max_abs_error_hartree",
        value,
        tolerance,
        pass: None,
        n_cases: cases.len(),
        level: "L2",
        units: "Ha",
        reference,
        notes: "",
    };
    gate(row("hf_reference_energy", max_abs(hf), 1e-6));
    gate(row("mp2_correlation_energy", max_abs(corr), 1e-6));
}

// ---------------------------------------------------------------------------
// F103 (V2-03) QM/MM mechanical embedding vs PySCF
// ---------------------------------------------------------------------------

#[test]
fn f103_qmmm_mechanical_embedding_matches_pyscf() {
    use sci_form::qmmm::{QmMmConfig, QmMmSystem};

    let data = fixture("v2_pyscf.json");
    let case = &data["qmmm"];
    let qm = molecule(&case["qm"]);
    let mm_positions: Vec<[f64; 3]> = serde_json::from_value(case["mm_positions"].clone()).unwrap();
    let mm_charges: Vec<f64> = serde_json::from_value(case["mm_charges"].clone()).unwrap();
    let system = QmMmSystem {
        qm_elements: qm.elements,
        qm_positions: qm.positions,
        mm_elements: vec![8, 1],
        mm_positions,
        mm_charges,
    };
    let result = sci_form::compute_qmmm_energy(&system, &QmMmConfig::default()).unwrap();
    let mechanical = case["mechanical_total_kcal_mol"].as_f64().unwrap();
    let electrostatic = case["electrostatic_total_kcal_mol"].as_f64().unwrap();
    let reference = "PySCF RHF/STO-3G + Mulliken point charges";
    gate(Row {
        ficha: "F103",
        observable: "mechanical_total_energy",
        metric: "abs_error_kcal_mol",
        value: (result.total_energy_kcal_mol - mechanical).abs(),
        tolerance: 1e-2,
        pass: None,
        n_cases: 1,
        level: "L2",
        units: "kcal/mol",
        reference,
        notes: "",
    });
    gap(Row {
        ficha: "F103",
        observable: "electrostatic_embedding_total_energy",
        metric: "abs_error_kcal_mol",
        value: (result.total_energy_kcal_mol - electrostatic).abs(),
        tolerance: 0.1,
        pass: None,
        n_cases: 1,
        level: "L3",
        units: "kcal/mol",
        reference: "PySCF qmmm.mm_charge RHF/STO-3G",
        notes: "mechanical embedding has no MM polarisation of the QM density",
    });
}

// ---------------------------------------------------------------------------
// F104 (V2-04) EPR isotropic hyperfine vs PySCF UHF spin density
// ---------------------------------------------------------------------------

#[test]
fn f104_epr_hyperfine_matches_pyscf_spin_density() {
    use sci_form::epr::EprConfig;

    let data = fixture("v2_pyscf.json");
    let cases = data["epr"].as_array().unwrap();
    // The SAD guess puts NH2 in its 2B1 ground state (the H_core guess used to
    // land on 2A1, 0.097 Ha higher).
    let mut errors = Vec::new();
    for case in cases {
        let mol = molecule(case);
        let result =
            sci_form::compute_epr_hyperfine(&mol.elements, &mol.positions, &EprConfig::default())
                .unwrap();
        for expected in case["entries"].as_array().unwrap() {
            let index = expected["atom_index"].as_u64().unwrap() as usize;
            let entry = result
                .entries
                .iter()
                .find(|e| e.atom_index == index)
                .unwrap_or_else(|| panic!("{}: no entry for atom {index}", mol.name));
            // sci-form scales A_iso by the STO-3G hydrogen 1s density (so the H
            // atom gives 1420.4 MHz by construction); PySCF's A_iso uses the
            // exact 1/π. The shared observable is the spin density at the nucleus.
            let reference = expected["spin_density_au"].as_f64().unwrap();
            eprintln!(
                "F104 {:<4} atom {index}: rho_s {:+.6} vs {reference:+.6} au",
                mol.name, entry.spin_density_au
            );
            errors.push((entry.spin_density_au - reference) / reference.abs().max(1e-3));
        }
    }
    let reference = "PySCF UHF/STO-3G spin density at the nucleus";
    gate(Row {
        ficha: "F104",
        observable: "spin_density_at_nucleus",
        metric: "max_rel_error",
        value: max_abs(errors.iter().copied()),
        tolerance: 1e-3,
        pass: None,
        n_cases: errors.len(),
        level: "L2",
        units: "1",
        reference,
        notes: "H, CH3, NH2, OH; A_iso is calibrated on the STO-3G H atom, so no L4 row on it",
    });
}

// ---------------------------------------------------------------------------
// F105 (V2-05) uncoupled polarizability vs PySCF
// ---------------------------------------------------------------------------

#[test]
fn f105_uncoupled_polarizability_matches_pyscf() {
    let data = fixture("v2_pyscf.json");
    let cases = data["polarizability"].as_array().unwrap();
    let mut uncoupled = Vec::new();
    let mut coupled = Vec::new();
    for case in cases {
        let mol = molecule(case);
        let alpha = sci_form::raman::compute_polarizability(&mol.elements, &mol.positions).unwrap();
        let reference = case["uncoupled_isotropic_ang3"].as_f64().unwrap();
        let cphf = case["coupled_isotropic_ang3"].as_f64().unwrap();
        eprintln!(
            "F105 {:<5} {:.6} vs UCHF {reference:.6} / CPHF {cphf:.6} Å³",
            mol.name, alpha.isotropic
        );
        uncoupled.push((alpha.isotropic - reference) / reference);
        coupled.push((alpha.isotropic - cphf) / cphf);
    }
    let row = |observable, value, tolerance, level, reference| Row {
        ficha: "F105",
        observable,
        metric: "max_rel_error",
        value,
        tolerance,
        pass: None,
        n_cases: cases.len(),
        level,
        units: "1",
        reference,
        notes: "",
    };
    gate(row(
        "uncoupled_isotropic_polarizability",
        max_abs(uncoupled),
        1e-3,
        "L2",
        "PySCF RHF/STO-3G sum-over-states (UCHF)",
    ));
    gap(Row {
        notes: "the prototype has no CPHF response",
        ..row(
            "coupled_isotropic_polarizability",
            max_abs(coupled),
            0.02,
            "L3",
            "PySCF finite-field RHF/STO-3G (CPHF)",
        )
    });
}

// ---------------------------------------------------------------------------
// F106 (V2-06) Koopmans core binding energies vs PySCF
// ---------------------------------------------------------------------------

#[test]
fn f106_xps_koopmans_matches_pyscf() {
    use sci_form::xps::XpsConfig;

    let data = fixture("v2_pyscf.json");
    let cases = data["xps"].as_array().unwrap();
    let mut errors = Vec::new();
    for case in cases {
        let mol = molecule(case);
        let result = sci_form::compute_xps_binding_energies(
            &mol.elements,
            &mol.positions,
            &XpsConfig::default(),
        )
        .unwrap();
        for expected in case["bindings"].as_array().unwrap() {
            let index = expected["atom_index"].as_u64().unwrap() as usize;
            let binding = result
                .bindings
                .iter()
                .find(|b| b.atom_index == index)
                .unwrap_or_else(|| panic!("{}: no binding for atom {index}", mol.name));
            errors
                .push(binding.binding_energy_ev - expected["binding_energy_ev"].as_f64().unwrap());
        }
    }
    gate(Row {
        ficha: "F106",
        observable: "koopmans_core_binding_energy",
        metric: "max_abs_error_ev",
        value: max_abs(errors.iter().copied()),
        tolerance: 0.01,
        pass: None,
        n_cases: errors.len(),
        level: "L2",
        units: "eV",
        reference: "PySCF UHF/STO-3G 1s orbital energies",
        notes: "",
    });
}

// ---------------------------------------------------------------------------
// F107 (V2-07) 3D-RISM ideal limit: chi = 0 => h = c => g = KH(-beta u)
// ---------------------------------------------------------------------------

#[test]
fn f107_rism_ideal_limit_is_exact() {
    use sci_form::rism::RismConfig;

    let config = RismConfig {
        nx: 16,
        ny: 16,
        nz: 16,
        box_length_angstrom: 12.0,
        tolerance: 1e-12,
        max_iterations: 2_000,
        ..RismConfig::default()
    };
    let n = config.nx;
    let spacing = config.box_length_angstrom / n as f64;
    let center = config.box_length_angstrom / 2.0;
    // Repulsive core (d < 0 branch) plus an attractive shell (d > 0 branch).
    let mut beta_u = Vec::with_capacity(n * n * n);
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let r = [i, j, k]
                    .map(|index| index as f64 * spacing - center)
                    .iter()
                    .map(|x| x * x)
                    .sum::<f64>()
                    .sqrt();
                beta_u.push(4.0 * (-r * r).exp() - 2.0 * (-(r - 3.0).powi(2)).exp());
            }
        }
    }
    let chi = vec![0.0; beta_u.len()];
    let result = sci_form::solve_3d_rism(&beta_u, &chi, &config).unwrap();
    assert!(result.converged);
    let error = max_abs(beta_u.iter().zip(&result.g).map(|(u, g)| {
        let d = -u;
        let exact = if d <= 0.0 { d.exp() } else { 1.0 + d };
        g - exact
    }));
    gate(Row {
        ficha: "F107",
        observable: "ideal_limit_kh_closure",
        metric: "max_abs_error_g",
        value: error,
        tolerance: 1e-8,
        pass: None,
        n_cases: beta_u.len(),
        level: "L0",
        units: "1",
        reference: "analytic: chi = 0 gives g = KH(-beta u)",
        notes: "",
    });
}

// ---------------------------------------------------------------------------
// F108 (V2-08) relative pKa vs experiment
// ---------------------------------------------------------------------------

#[test]
fn f108_relative_pka_vs_experiment() {
    use sci_form::electrochemistry::{ElectrochemistryConfig, PkaReference};

    // CRC Handbook of Chemistry and Physics, aqueous pKa at 25 °C.
    let acids = [
        ("formic", "OC=O", "[O-]C=O", 3.75),
        ("propanoic", "CCC(=O)O", "CCC(=O)[O-]", 4.87),
        ("chloroacetic", "OC(=O)CCl", "[O-]C(=O)CCl", 2.87),
        ("benzoic", "OC(=O)c1ccccc1", "[O-]C(=O)c1ccccc1", 4.20),
    ];
    let reference = PkaReference::acetic_acid();
    let config = ElectrochemistryConfig::default();
    let mut errors = Vec::new();
    for (name, acid, base, experiment) in acids {
        let result = sci_form::compute_pka(acid, base, Some(&reference), &config).unwrap();
        eprintln!(
            "F108 {name:<13} pKa {:.2} vs exp {experiment:.2}",
            result.pka
        );
        errors.push(result.pka - experiment);
    }
    let mae = errors.iter().map(|e| e.abs()).sum::<f64>() / errors.len() as f64;
    let row = Row {
        ficha: "F108",
        observable: "relative_pka_carboxylic_acids",
        metric: "mae_pka_units",
        value: mae,
        tolerance: 1.0,
        pass: None,
        n_cases: errors.len(),
        level: "L4",
        units: "pKa",
        reference: "CRC Handbook; relative to acetic acid 4.76",
        notes: "single-point UHF/STO-3G on unoptimised ETKDG geometries; with 8 conformers the pKa moves by up to 8 units",
    };
    if mae <= 1.0 {
        gate(row);
    } else {
        gap(row);
    }
}

// ---------------------------------------------------------------------------
// F109 (V2-09) FEP / BAR / MBAR vs pymbar and the exact answer
// ---------------------------------------------------------------------------

#[test]
fn f109_free_energy_estimators_match_pymbar() {
    use sci_form::free_energy::{
        estimate_free_energy, AlchemicalSample, FreeEnergyConfig, FreeEnergyEstimator,
    };

    let data = fixture("v2_pymbar.json");
    let samples: Vec<AlchemicalSample> = serde_json::from_value(data["samples"].clone()).unwrap();
    let floats = |key: &str| -> Vec<f64> { serde_json::from_value(data[key].clone()).unwrap() };
    let (mbar_ref, exp_ref, exact) = (
        floats("mbar_delta_f"),
        floats("exp_delta_f"),
        floats("exact_delta_f"),
    );
    let config = |n_states| FreeEnergyConfig {
        n_states,
        temperature_k: None,
        tolerance: 1e-12,
        ..FreeEnergyConfig::default()
    };

    let mbar = estimate_free_energy(&samples, FreeEnergyEstimator::Mbar, &config(3)).unwrap();
    assert!(mbar.converged);

    // BAR: the samples of states 0 and 1, restricted to those two states.
    let pair: Vec<AlchemicalSample> = samples
        .iter()
        .filter(|s| s.origin_state < 2)
        .map(|s| AlchemicalSample {
            origin_state: s.origin_state,
            reduced_potentials: s.reduced_potentials[..2].to_vec(),
        })
        .collect();
    let bar = estimate_free_energy(&pair, FreeEnergyEstimator::Bar, &config(2)).unwrap();

    // FEP (EXP): state-0 samples, from state 0 to each state j.
    let mut fep_errors = Vec::new();
    for (j, reference) in exp_ref.iter().enumerate().skip(1) {
        let forward: Vec<AlchemicalSample> = samples
            .iter()
            .filter(|s| s.origin_state == 0)
            .map(|s| AlchemicalSample {
                origin_state: 0,
                reduced_potentials: vec![s.reduced_potentials[0], s.reduced_potentials[j]],
            })
            .collect();
        let fep = estimate_free_energy(&forward, FreeEnergyEstimator::Fep, &config(2)).unwrap();
        fep_errors.push(fep.delta_f_kbt[1] - reference);
    }

    let row = |observable, value, tolerance, level, reference| Row {
        ficha: "F109",
        observable,
        metric: "max_abs_error_kbt",
        value,
        tolerance,
        pass: None,
        n_cases: 3,
        level,
        units: "kBT",
        reference,
        notes: "",
    };
    let pymbar = "pymbar 4 on identical samples";
    gate(row(
        "mbar_vs_pymbar",
        max_abs(mbar.delta_f_kbt.iter().zip(&mbar_ref).map(|(a, b)| a - b)),
        1e-6,
        "L2",
        pymbar,
    ));
    gate(row(
        "bar_vs_pymbar",
        (bar.delta_f_kbt[1] - data["bar_delta_f_01"].as_f64().unwrap()).abs(),
        1e-6,
        "L2",
        pymbar,
    ));
    gate(row(
        "fep_vs_pymbar_exp",
        max_abs(fep_errors),
        1e-6,
        "L2",
        pymbar,
    ));
    gate(row(
        "mbar_vs_exact",
        max_abs(mbar.delta_f_kbt.iter().zip(&exact).map(|(a, b)| a - b)),
        0.05,
        "L0",
        "analytic: Delta f = ln(K_k/K_0)/2 (500 samples per state)",
    ));
}

// ---------------------------------------------------------------------------
// F110 (V2-10) well-tempered metadynamics barrier on an analytic double well
// ---------------------------------------------------------------------------

#[test]
fn f110_metadynamics_recovers_double_well_barrier() {
    use sci_form::enhanced_sampling::{
        run_metadynamics, CollectiveVariable, DoubleWell, MetaDynamicsConfig,
    };

    let barrier = 3.0;
    let potential = DoubleWell::new(barrier, 1.0);
    let config = MetaDynamicsConfig {
        cv: CollectiveVariable::Coordinate { index: 0 },
        height_kcal_mol: 0.2,
        sigma: 0.15,
        bias_factor: 10.0,
        pace_steps: 50,
        temperature_k: 300.0,
        dt_fs: 0.5,
        // 1 ps⁻¹ cannot remove the energy the deposition injects (0.2 kcal/mol
        // every 25 fs) and the walker runs hot; 10 ps⁻¹ keeps it at 300 K.
        friction_ps_inv: 10.0,
        max_steps: 40_000,
        seed: 7,
        grid_min: -1.5,
        grid_max: 1.5,
        grid_bins: 61,
        trace_stride: 100,
    };
    let result = run_metadynamics(&potential, &[-1.0], &[1.0], &config).unwrap();
    let fes = &result.free_energy;
    let at = |x: f64| {
        let index = ((x - fes.min) / (fes.max - fes.min) * (fes.bins - 1) as f64).round() as usize;
        fes.values[index]
    };
    let estimated = at(0.0) - at(-1.0).min(at(1.0));
    eprintln!("F110 barrier {estimated:.3} vs exact {barrier:.3} kcal/mol");
    gate(Row {
        ficha: "F110",
        observable: "double_well_barrier",
        metric: "abs_error_kcal_mol",
        value: (estimated - barrier).abs(),
        tolerance: 1.0,
        pass: None,
        n_cases: 1,
        level: "L0",
        units: "kcal/mol",
        reference: "analytic: F(x) = U(x) for the coordinate CV",
        notes: "",
    });
}

// ---------------------------------------------------------------------------
// F111 (V2-11) GCMC ideal-gas limit <N> = PV/kT
// ---------------------------------------------------------------------------

#[test]
fn f111_gcmc_ideal_gas_loading() {
    use sci_form::gcmc::{GcmcConfig, GcmcFramework, GuestModel};

    let framework = GcmcFramework::empty_cubic(30.0);
    let guest = GuestModel::single_site(0.0, 3.0, 0.0);
    let (temperature, pressure) = (300.0, 50.0);
    let config = GcmcConfig {
        temperature_k: temperature,
        pressures_bar: vec![pressure],
        equilibration_steps: 20_000,
        production_steps: 200_000,
        seed: 11,
        translation_probability: 0.0,
        ..GcmcConfig::default()
    };
    let result = sci_form::compute_gcmc_isotherm(&framework, &guest, &config).unwrap();
    let volume_m3 = 30.0f64.powi(3) * 1e-30;
    let ideal = pressure * 1e5 * volume_m3 / (1.380_649e-23 * temperature);
    let loading = result.points[0].mean_loading;
    eprintln!("F111 <N> = {loading:.3} vs PV/kT = {ideal:.3}");
    gate(Row {
        ficha: "F111",
        observable: "ideal_gas_loading",
        metric: "rel_error",
        value: ((loading - ideal) / ideal).abs(),
        tolerance: 0.05,
        pass: None,
        n_cases: 1,
        level: "L0",
        units: "1",
        reference: "analytic: <N> = PV/kT",
        notes: "",
    });
}

// ---------------------------------------------------------------------------
// F112 (V2-12) kinetic Monte Carlo first-order decay
// ---------------------------------------------------------------------------

#[test]
fn f112_kmc_first_order_decay() {
    use sci_form::kinetics::{evaluate_htst_rate, ElementaryStep, ThermodynamicState};
    use sci_form::kmc::{run_kinetic_monte_carlo, KmcConfig, KmcSpecies, KmcStep};

    let step = ElementaryStep {
        step_id: "A->B".into(),
        activation_free_energy_ev: 0.6,
        reaction_free_energy_ev: -0.3,
        prefactor_s_inv: Some(1e13),
    };
    let state = ThermodynamicState::default();
    let rate = evaluate_htst_rate(&step, state).unwrap().forward_rate_s_inv;
    let (n0, seeds) = (100u64, 200u64);
    let mut survivors = 0u64;
    for seed in 0..seeds {
        let config = KmcConfig {
            species: vec![
                KmcSpecies {
                    id: "A".into(),
                    initial_count: n0,
                },
                KmcSpecies {
                    id: "B".into(),
                    initial_count: 0,
                },
            ],
            steps: vec![KmcStep {
                step: step.clone(),
                reactants: vec![(0, 1)],
                products: vec![(1, 1)],
            }],
            state,
            max_events: 10_000,
            max_time_s: 1.0 / rate,
            seed,
            frame_stride: 1_000,
        };
        survivors += run_kinetic_monte_carlo(&config).unwrap().final_counts[0];
    }
    let fraction = survivors as f64 / (n0 * seeds) as f64;
    let exact = (-1.0f64).exp();
    eprintln!("F112 N_A(1/k)/N0 = {fraction:.4} vs e^-1 = {exact:.4}");
    gate(Row {
        ficha: "F112",
        observable: "first_order_survival_at_one_lifetime",
        metric: "rel_error",
        value: ((fraction - exact) / exact).abs(),
        tolerance: 0.05,
        pass: None,
        n_cases: seeds as usize,
        level: "L0",
        units: "1",
        reference: "analytic: N_A(t)/N0 = exp(-k t)",
        notes: "",
    });
}

// ---------------------------------------------------------------------------
// F113 (V2-13) MLIP contract: units and forces
// ---------------------------------------------------------------------------

#[test]
fn f113_mlip_contract_units_and_forces() {
    use sci_form::error::SciFormError;
    use sci_form::mlip::{CallbackMlip, MlipEnergyUnit, MlipManifest, MlipModel};

    // Harmonic bond in eV: E = k (r - r0)^2 / 2, F_i = -dE/dx_i.
    let (k, r0) = (2.0, 1.0);
    let harmonic = move |_: &[u8], x: &[[f64; 3]], _: Option<&[[f64; 3]; 3]>| {
        let d = [x[1][0] - x[0][0], x[1][1] - x[0][1], x[1][2] - x[0][2]];
        let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        let energy = 0.5 * k * (r - r0).powi(2);
        let scale = k * (r - r0) / r;
        let f1 = [-scale * d[0], -scale * d[1], -scale * d[2]];
        Ok::<_, SciFormError>((energy, vec![[-f1[0], -f1[1], -f1[2]], f1], None))
    };
    let model = CallbackMlip::new(
        MlipManifest {
            schema_version: 1,
            name: "harmonic-ev".into(),
            elements: vec![1],
            cutoff_angstrom: 5.0,
            license: "test-only".into(),
            energy_unit: MlipEnergyUnit::Ev,
            supports_stress: false,
            source: None,
        },
        harmonic,
    )
    .unwrap();
    let elements = [1, 1];
    let positions = [[0.0, 0.0, 0.0], [0.9, 0.3, -0.2]];
    let result = model.evaluate(&elements, &positions, None).unwrap();
    let (energy_ev, _, _) = harmonic(&elements, &positions, None).unwrap();
    let unit_error = (result.energy_kcal_mol - energy_ev * EV_TO_KCAL).abs();

    let h = 1e-5;
    let mut force_error: f64 = 0.0;
    for atom in 0..2 {
        for axis in 0..3 {
            let mut plus = positions;
            let mut minus = positions;
            plus[atom][axis] += h;
            minus[atom][axis] -= h;
            let e_plus = model
                .evaluate(&elements, &plus, None)
                .unwrap()
                .energy_kcal_mol;
            let e_minus = model
                .evaluate(&elements, &minus, None)
                .unwrap()
                .energy_kcal_mol;
            let numerical = -(e_plus - e_minus) / (2.0 * h);
            force_error = force_error
                .max((result.forces_kcal_mol_per_angstrom[atom][axis] - numerical).abs());
        }
    }
    let row = |observable, metric, value, tolerance| Row {
        ficha: "F113",
        observable,
        metric,
        value,
        tolerance,
        pass: None,
        n_cases: 1,
        level: "L0",
        units: "kcal/mol",
        reference: "analytic harmonic bond declared in eV",
        notes: "",
    };
    gate(row(
        "ev_energy_converted_to_kcal_mol",
        "abs_error_kcal_mol",
        unit_error,
        1e-9,
    ));
    gate(row(
        "forces_are_minus_energy_gradient",
        "max_abs_error_kcal_mol_per_angstrom",
        force_error,
        1e-5,
    ));
}

// ---------------------------------------------------------------------------
// F114 (V2-14) split-conformal coverage
// ---------------------------------------------------------------------------

#[test]
fn f114_conformal_coverage_is_nominal() {
    use sci_form::active_learning::ConformalCalibrator;

    let mut normal = Normal(20_261_119);
    // y = f(x) + N(0, 1) with a perfect point predictor f.
    let mut draw = |n: usize| -> (Vec<f64>, Vec<f64>) {
        (0..n)
            .map(|i| {
                let prediction = (i as f64 * 0.01).sin();
                (prediction, prediction + normal.sample())
            })
            .unzip()
    };
    let (calibration_pred, calibration_target) = draw(2_000);
    let calibrator =
        ConformalCalibrator::fit(&calibration_pred, &calibration_target, 0.90).unwrap();
    let (test_pred, test_target) = draw(20_000);
    let covered = test_pred
        .iter()
        .zip(&test_target)
        .filter(|(p, t)| (*p - *t).abs() <= calibrator.radius)
        .count();
    let coverage = covered as f64 / test_pred.len() as f64;
    eprintln!(
        "F114 empirical coverage {coverage:.4}, radius {:.4}",
        calibrator.radius
    );
    gate(Row {
        ficha: "F114",
        observable: "split_conformal_coverage",
        metric: "abs_deviation_from_0.90",
        value: (coverage - 0.90).abs(),
        tolerance: 0.02,
        pass: None,
        n_cases: test_pred.len(),
        level: "L0",
        units: "1",
        reference: "exchangeability guarantee: coverage >= 0.90",
        notes: "",
    });
}

// ---------------------------------------------------------------------------
// F115 (V2-15) all-nuclide NMR driver == direct GIAO per nucleus
// ---------------------------------------------------------------------------

#[test]
fn f115_all_nuclide_nmr_matches_direct_giao() {
    use sci_form::GiaoNmrConfig;

    let elements = [8u8, 1, 1];
    let positions = [
        [0.0, 0.0, 0.0],
        [0.756_950, 0.0, 0.585_882],
        [-0.756_950, 0.0, 0.585_882],
    ];
    let config = GiaoNmrConfig::default();
    let all = sci_form::nmr_all::compute_all_nmr(&elements, &positions, &config).unwrap();
    assert!(all.failures.is_empty(), "failures: {:?}", all.failures);
    let mut error: f64 = 0.0;
    let mut n = 0;
    for result in &all.nuclei {
        let direct =
            sci_form::compute_giao_nmr_configured(&elements, &positions, &result.nucleus, &config)
                .unwrap();
        for (a, b) in result.shieldings.iter().zip(&direct.shieldings) {
            assert_eq!(a.atom_index, b.atom_index);
            error = error.max((a.isotropic - b.isotropic).abs());
            n += 1;
        }
    }
    assert_eq!(n, elements.len());
    gate(Row {
        ficha: "F115",
        observable: "all_nuclide_driver_vs_direct_giao",
        metric: "max_abs_error_ppm",
        value: error,
        tolerance: 1e-10,
        pass: None,
        n_cases: n,
        level: "L1",
        units: "ppm",
        reference: "sci-form compute_giao_nmr_configured",
        notes: "PySCF GIAO oracle pending (pyscf-properties)",
    });
}
