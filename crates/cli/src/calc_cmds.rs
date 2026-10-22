//! EHT, SASA, population, dipole, ESP, DOS, cell, assemble, ANI, PM3, xTB, GFN, HF-3c command handlers.

use crate::format::{
    fail, gasteiger_charges_from_geometry, parse_bond_list, parse_elems_coords, perceive_bonds,
};

pub fn cmd_eht(elements: &str, coords: &str, k: Option<f64>) {
    let (elems, flat, _) = parse_elems_coords(elements, coords);
    if flat.len() != elems.len() * 3 {
        fail(format!(
            "Invalid 'coords' length {}. Expected {} values (3 per atom).",
            flat.len(),
            elems.len() * 3
        ));
    }
    let positions: Vec<[f64; 3]> = flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    // `--k` omitted means "engine default"; an explicit value must be
    // positive (blocker C-M5: no magic-number sentinel).
    let k_opt = match k {
        None => None,
        Some(v) if v > 0.0 => Some(v),
        Some(v) => fail(format!(
            "Invalid 'k' value '{v}'. Expected a positive number (omit --k for the default 1.75)."
        )),
    };
    match sci_form::eht::solve_eht(&elems, &positions, k_opt) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("EHT error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_sasa(elements: &str, coords: &str, probe_radius: f64) {
    let (elems, flat, _) = parse_elems_coords(elements, coords);
    match sci_form::compute_sasa(&elems, &flat, Some(probe_radius)) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("SASA error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_population(elements: &str, coords: &str) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    match sci_form::compute_population(&elems, &positions) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("Population error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_dipole(elements: &str, coords: &str) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    match sci_form::compute_dipole(&elems, &positions) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("Dipole error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_esp(elements: &str, coords: &str, spacing: f64, padding: f64) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    match sci_form::compute_esp(&elems, &positions, spacing, padding) {
        Ok(grid) => println!(
            "{{\"origin\":{:?},\"spacing\":{},\"dims\":{:?},\"n_values\":{},\"values\":{:?}}}",
            grid.origin,
            grid.spacing,
            grid.dims,
            grid.values.len(),
            grid.values
        ),
        Err(e) => {
            eprintln!("ESP error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_dos(elements: &str, coords: &str, sigma: f64, e_min: f64, e_max: f64, n_points: usize) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    match sci_form::compute_dos(&elems, &positions, sigma, e_min, e_max, n_points) {
        Ok(result) => {
            let pairs: Vec<_> = result
                .energies
                .iter()
                .zip(result.total_dos.iter())
                .map(|(e, d)| format!("[{:.4},{:.6}]", e, d))
                .collect();
            println!(
                "{{\"sigma\":{},\"data\":[{}]}}",
                result.sigma,
                pairs.join(",")
            );
        }
        Err(e) => {
            eprintln!("DOS error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_cell(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) {
    let cell = sci_form::create_unit_cell(a, b, c, alpha, beta, gamma);
    let vol = cell.volume();
    let p = cell.parameters();
    println!(
        "{{\"a\":{:.4},\"b\":{:.4},\"c\":{:.4},\"alpha\":{:.2},\"beta\":{:.2},\"gamma\":{:.2},\"volume\":{:.4},\"lattice\":{:?}}}",
        p.a, p.b, p.c, p.alpha, p.beta, p.gamma, vol, cell.lattice
    );
}

pub fn cmd_assemble(topology: &str, a: f64, metal: u8, geometry: &str, supercell: usize) {
    let geom = match geometry {
        "linear" => sci_form::materials::CoordinationGeometry::Linear,
        "trigonal" => sci_form::materials::CoordinationGeometry::Trigonal,
        "tetrahedral" => sci_form::materials::CoordinationGeometry::Tetrahedral,
        "square_planar" => sci_form::materials::CoordinationGeometry::SquarePlanar,
        "octahedral" => sci_form::materials::CoordinationGeometry::Octahedral,
        _ => {
            fail(crate::format::invalid_value(
                "'geometry'",
                geometry,
                &[
                    "linear",
                    "trigonal",
                    "tetrahedral",
                    "square_planar",
                    "octahedral",
                ],
            ));
        }
    };
    let topo = match topology {
        "pcu" => sci_form::materials::Topology::pcu(),
        "dia" => sci_form::materials::Topology::dia(),
        "sql" => sci_form::materials::Topology::sql(),
        _ => {
            fail(crate::format::invalid_value(
                "'topology'",
                topology,
                &["pcu", "dia", "sql"],
            ));
        }
    };
    let node = sci_form::materials::Sbu::metal_node(metal, 0.0, geom);
    let linker = sci_form::materials::Sbu::linear_linker(&[6, 6], 1.4, "carboxylate");
    let cell = sci_form::materials::UnitCell::cubic(a);
    let mut structure = sci_form::assemble_framework(&node, &linker, &topo, &cell);
    if supercell > 1 {
        structure = structure.make_supercell(supercell, supercell, supercell);
    }
    println!("{}", serde_json::to_string_pretty(&structure).unwrap());
}

pub fn cmd_ani(elements: &str, coords: &str) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    // The network input dimension must match the real AEV length (4032 for the
    // default ANI-2x-style parameters); using a different value makes every
    // call fail with an input-dimension mismatch.
    let aev_len = sci_form::ani::aev_params::default_ani2x_params().total_aev_length();
    let mut network_map = std::collections::HashMap::new();
    for &element in &[1u8, 6, 7, 8] {
        network_map.insert(element, sci_form::ani::weights::make_test_model(aev_len));
    }
    let config = sci_form::ani::api::AniConfig {
        compute_forces: true,
        output_aevs: true,
        ..Default::default()
    };
    match sci_form::ani::api::compute_ani(&elems, &positions, &config, &network_map) {
        Ok(result) => {
            eprintln!(
                "note: using randomly-initialized test networks ({} AEV inputs); \
                 energies are not physically meaningful.",
                aev_len
            );
            println!("{}", serde_json::to_string_pretty(&result).unwrap());
        }
        Err(e) => {
            eprintln!("ANI error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_pm3(elements: &str, coords: &str) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    match sci_form::compute_pm3(&elems, &positions) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("PM3 error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_xtb(elements: &str, coords: &str) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    match sci_form::compute_xtb(&elems, &positions) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("xTB error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_gfn1(elements: &str, coords: &str) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    match sci_form::xtb::solve_gfn1(&elems, &positions) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("GFN1-xTB error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_gfn2(elements: &str, coords: &str) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    match sci_form::xtb::solve_gfn2(&elems, &positions) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("GFN2-xTB error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_hf3c(elements: &str, coords: &str) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    let config = sci_form::hf::api::HfConfig {
        n_cis_states: 5,
        corrections: true,
        ..Default::default()
    };
    match sci_form::hf::api::solve_hf3c(&elems, &positions, &config) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("HF-3c error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_stereo(smiles: &str, coords: Option<&str>) {
    // `--coords` omitted means "topology only" (blocker C-M5).
    let flat: Vec<f64> = match coords {
        None => vec![],
        Some(raw) => serde_json::from_str(raw).unwrap_or_else(|e| {
            fail(crate::format::bad_json("'coords'", raw, &e));
        }),
    };
    match sci_form::analyze_stereo(smiles, &flat) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("Stereo error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_solvation(elements: &str, coords: &str, charges: Option<&str>, probe_radius: f64) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    // `--charges` omitted means "non-polar contribution only"
    // (blocker C-M5).
    match charges {
        None => {
            let r = sci_form::compute_nonpolar_solvation(&elems, &positions, Some(probe_radius))
                .unwrap_or_else(|e| fail(format!("Solvation error: {e}")));
            println!("{}", serde_json::to_string_pretty(&r).unwrap());
        }
        Some(raw) => {
            let q: Vec<f64> = serde_json::from_str(raw).unwrap_or_else(|e| {
                fail(crate::format::bad_json("'charges'", raw, &e));
            });
            let r = sci_form::compute_gb_solvation(
                &elems,
                &positions,
                &q,
                Some(78.5),
                Some(1.0),
                Some(probe_radius),
            )
            .unwrap_or_else(|e| fail(format!("Solvation error: {e}")));
            println!("{}", serde_json::to_string_pretty(&r).unwrap());
        }
    }
}

pub fn cmd_sssr(smiles: &str) {
    match sci_form::compute_sssr(smiles) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("SSSR error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_ecfp(smiles: &str, radius: usize, n_bits: usize) {
    match sci_form::compute_ecfp(smiles, radius, n_bits) {
        Ok(fp) => println!("{}", serde_json::to_string_pretty(&fp).unwrap()),
        Err(e) => {
            eprintln!("ECFP error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_tanimoto(smiles1: &str, smiles2: &str, radius: usize, n_bits: usize) {
    let fp1 = match sci_form::compute_ecfp(smiles1, radius, n_bits) {
        Ok(fp) => fp,
        Err(e) => {
            eprintln!("ECFP error for '{}': {}", smiles1, e);
            std::process::exit(1);
        }
    };
    let fp2 = match sci_form::compute_ecfp(smiles2, radius, n_bits) {
        Ok(fp) => fp,
        Err(e) => {
            eprintln!("ECFP error for '{}': {}", smiles2, e);
            std::process::exit(1);
        }
    };
    let t = sci_form::compute_tanimoto(&fp1, &fp2);
    println!("{{\"tanimoto\": {:.6}}}", t);
}

// ─── Cross-surface parity commands ─────────────────────────────────────────

pub fn cmd_ecfp_batch(smiles_list: &str, radius: usize, n_bits: usize) {
    let smiles: Vec<String> = serde_json::from_str(smiles_list).unwrap_or_else(|e| {
        fail(crate::format::bad_json("'smiles list'", smiles_list, &e));
    });
    let refs: Vec<&str> = smiles.iter().map(|s| s.as_str()).collect();
    match sci_form::compute_ecfp_batch(&refs, radius, n_bits) {
        Ok(fps) => println!("{}", serde_json::to_string_pretty(&fps).unwrap()),
        Err(e) => {
            eprintln!("ECFP batch error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_charges_configured(
    smiles: &str,
    max_iter: usize,
    initial_damping: f64,
    convergence_threshold: f64,
) {
    let config = sci_form::charges::gasteiger::GasteigerConfig {
        max_iter,
        initial_damping,
        convergence_threshold,
    };
    match sci_form::compute_charges_configured(smiles, &config) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("Charges error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_ensemble_j_couplings(
    smiles: &str,
    conformer_coords: &str,
    energies_kcal: &str,
    temperature_k: f64,
) {
    let conformers: Vec<Vec<f64>> = serde_json::from_str(conformer_coords).unwrap_or_else(|e| {
        fail(crate::format::bad_json(
            "'conformer coords'",
            conformer_coords,
            &e,
        ));
    });
    let energies: Vec<f64> = serde_json::from_str(energies_kcal).unwrap_or_else(|e| {
        fail(crate::format::bad_json(
            "'conformer energies'",
            energies_kcal,
            &e,
        ));
    });
    match sci_form::compute_ensemble_j_couplings(smiles, &conformers, &energies, temperature_k) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("Ensemble J-coupling error: {}", e);
            std::process::exit(1);
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub fn cmd_ir_spectrum_broadened(
    elements: &str,
    coords: &str,
    method: &str,
    gamma: f64,
    wn_min: f64,
    wn_max: f64,
    n_points: usize,
    broadening: &str,
    step_size: Option<f64>,
    smiles: Option<&str>,
) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    // `--step-size` omitted means "engine default" (blocker C-M5).
    let step = step_size;
    let analysis = if method.eq_ignore_ascii_case("uff") {
        let smiles = match smiles {
            Some(s) => s,
            None => fail("Missing '--smiles': it is required for method=uff."),
        };
        sci_form::compute_vibrational_analysis_uff(smiles, &elems, &positions, step)
    } else {
        sci_form::compute_vibrational_analysis(&elems, &positions, method, step)
    };
    let analysis = match analysis {
        Ok(a) => a,
        Err(e) => {
            eprintln!("Vibrational analysis error: {}", e);
            std::process::exit(1);
        }
    };
    let spectrum = sci_form::compute_ir_spectrum_broadened(
        &analysis, gamma, wn_min, wn_max, n_points, broadening,
    );
    println!("{}", serde_json::to_string_pretty(&spectrum).unwrap());
}

pub fn cmd_cisd(elements: &str, coords: &str, n_states: usize) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    let system = sci_form::scf::types::MolecularSystem::from_angstrom(&elems, &positions, 0, 1);
    let scf_config = sci_form::scf::scf_loop::ScfConfig::default();
    let scf = sci_form::scf::scf_loop::run_scf(&system, &scf_config);
    if !scf.converged {
        eprintln!("CISD error: RHF/STO-3G SCF did not converge");
        std::process::exit(1);
    }
    let basis =
        sci_form::scf::basis::BasisSet::sto3g(&system.atomic_numbers, &system.positions_bohr);
    let n_basis = basis.n_basis;
    let engine = sci_form::scf::two_electron::TwoElectronIntegrals::compute(&basis);
    let n_pair = n_basis * (n_basis + 1) / 2;
    let mut eris = vec![0.0f64; n_pair * (n_pair + 1) / 2];
    for i in 0..n_basis {
        for j in 0..n_basis {
            for k in 0..n_basis {
                for l in 0..n_basis {
                    eris[sci_form::hf::integrals::eri_index(i, j, k, l, n_basis)] =
                        engine.get(i, j, k, l);
                }
            }
        }
    }
    let n_occupied = scf.n_electrons / 2;
    let result = sci_form::hf::compute_cisd(
        &scf.orbital_energies,
        &scf.mo_coefficients,
        &eris,
        n_basis,
        n_occupied,
        n_states,
    );
    println!("{}", serde_json::to_string_pretty(&result).unwrap());
}

pub fn cmd_dos_multimethod(
    elements: &str,
    coords: &str,
    method: &str,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    let dos_method = match method.to_ascii_lowercase().as_str() {
        "eht" => sci_form::dos::multi_method::DosMethod::Eht,
        "pm3" => sci_form::dos::multi_method::DosMethod::Pm3,
        "xtb" | "gfn0" => sci_form::dos::multi_method::DosMethod::Xtb,
        "gfn1" => sci_form::dos::multi_method::DosMethod::Gfn1,
        "gfn2" => sci_form::dos::multi_method::DosMethod::Gfn2,
        "hf3c" | "hf-3c" => sci_form::dos::multi_method::DosMethod::Hf3c,
        other => {
            eprintln!(
                "Invalid 'method' '{}'. Supported: 'eht', 'pm3', 'xtb', 'gfn1', 'gfn2', 'hf3c'.",
                other
            );
            std::process::exit(1);
        }
    };
    match sci_form::dos::multi_method::compute_dos_multimethod(
        &elems, &positions, dos_method, sigma, e_min, e_max, n_points,
    ) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("DOS error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_eht_gradient(elements: &str, coords: &str) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    match sci_form::eht::gradients::compute_eht_gradient(&elems, &positions) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("EHT gradient error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_optimize_geometry_eht(
    elements: &str,
    coords: &str,
    max_steps: usize,
    grad_threshold: f64,
    energy_threshold: f64,
    step_size: f64,
) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    let config = sci_form::eht::gradients::EhtOptConfig {
        max_steps,
        grad_threshold,
        energy_threshold,
        step_size,
    };
    match sci_form::eht::gradients::optimize_geometry_eht(&elems, &positions, Some(config)) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("EHT optimization error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_aevs_tm(elements: &str, coords: &str) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    let params = sci_form::ani::ani_tm::default_ani_tm_params();
    let aevs = sci_form::ani::ani_tm::compute_aevs_tm(&elems, &positions, &params);
    let aev_length = aevs.first().map(|a| a.len()).unwrap_or(0);
    let payload = serde_json::json!({
        "n_atoms": elems.len(),
        "aev_length": aev_length,
        "aevs": aevs,
    });
    println!("{}", serde_json::to_string_pretty(&payload).unwrap());
}

pub fn cmd_whim(elements: &str, coords: &str, weighting: &str) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    let w = match weighting.to_ascii_lowercase().as_str() {
        "unit" | "unweighted" => sci_form::ml::whim::WhimWeighting::Unit,
        "mass" => sci_form::ml::whim::WhimWeighting::Mass,
        "volume" | "vdw" => sci_form::ml::whim::WhimWeighting::Volume,
        "electronegativity" | "en" => sci_form::ml::whim::WhimWeighting::Electronegativity,
        "polarizability" => sci_form::ml::whim::WhimWeighting::Polarizability,
        other => {
            eprintln!(
                "Invalid 'weighting' '{}'. Supported: 'unit', 'unweighted', 'mass', \
                 'volume', 'vdw', 'electronegativity', 'en', 'polarizability'.",
                other
            );
            std::process::exit(1);
        }
    };
    let result = sci_form::ml::whim::compute_whim(&elems, &positions, w);
    println!("{}", serde_json::to_string_pretty(&result).unwrap());
}

pub fn cmd_rdf(
    elements: &str,
    coords: &str,
    charges: Option<&str>,
    beta: f64,
    r_max: f64,
    n_points: usize,
) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    // `--charges` omitted means "derive them with Gasteiger" (blocker C-M5).
    let q = match charges {
        None => match gasteiger_charges_from_geometry(&elems, &positions) {
            Ok(q) => q,
            Err(e) => {
                fail(format!("Gasteiger charge error: {e}"));
            }
        },
        Some(raw) => serde_json::from_str(raw).unwrap_or_else(|e| {
            fail(crate::format::bad_json("'charges'", raw, &e));
        }),
    };
    let result =
        sci_form::ml::rdf_descriptors::compute_rdf(&elems, &positions, &q, beta, r_max, n_points);
    println!("{}", serde_json::to_string_pretty(&result).unwrap());
}

pub fn cmd_getaway(elements: &str, coords: &str, bonds: Option<&str>, max_lag: usize) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    let bond_list = match parse_bond_list(bonds) {
        Some(b) => b,
        None => perceive_bonds(&elems, &positions, 0.3),
    };
    let result = sci_form::ml::getaway::compute_getaway(&elems, &positions, &bond_list, max_lag);
    println!("{}", serde_json::to_string_pretty(&result).unwrap());
}

fn rmse(predictions: &[f64], targets: &[f64]) -> f64 {
    if predictions.is_empty() {
        return 0.0;
    }
    let sum: f64 = predictions
        .iter()
        .zip(targets.iter())
        .map(|(p, t)| (p - t) * (p - t))
        .sum();
    (sum / predictions.len() as f64).sqrt()
}

fn parse_training_data(features: &str, targets: &str) -> (Vec<Vec<f64>>, Vec<f64>) {
    let feats: Vec<Vec<f64>> = serde_json::from_str(features).unwrap_or_else(|e| {
        fail(crate::format::bad_json("'features'", features, &e));
    });
    let tgts: Vec<f64> = serde_json::from_str(targets).unwrap_or_else(|e| {
        fail(crate::format::bad_json("'targets'", targets, &e));
    });
    if feats.len() != tgts.len() {
        eprintln!(
            "features length {} != targets length {}",
            feats.len(),
            tgts.len()
        );
        std::process::exit(1);
    }
    (feats, tgts)
}

#[allow(clippy::too_many_arguments)]
pub fn cmd_train_random_forest(
    features: &str,
    targets: &str,
    n_trees: usize,
    max_depth: usize,
    min_samples_split: usize,
    min_samples_leaf: usize,
    sample_fraction: f64,
    seed: u64,
) {
    let (feats, tgts) = parse_training_data(features, targets);
    let config = sci_form::ml::advanced_models::RandomForestConfig {
        n_trees,
        tree_config: sci_form::ml::advanced_models::TreeConfig {
            max_depth,
            min_samples_split,
            min_samples_leaf,
            max_features: 0,
        },
        sample_fraction,
        seed,
    };
    let model = sci_form::ml::advanced_models::train_random_forest(&feats, &tgts, &config);
    let predictions: Vec<f64> = feats.iter().map(|f| model.predict(f)).collect();
    let payload = serde_json::json!({
        "model_type": "random_forest",
        "n_estimators": model.n_trees,
        "n_samples": tgts.len(),
        "train_rmse": rmse(&predictions, &tgts),
        "predictions": predictions,
    });
    println!("{}", serde_json::to_string_pretty(&payload).unwrap());
}

#[allow(clippy::too_many_arguments)]
pub fn cmd_train_gradient_boosting(
    features: &str,
    targets: &str,
    n_estimators: usize,
    learning_rate: f64,
    max_depth: usize,
    min_samples_split: usize,
    min_samples_leaf: usize,
    subsample: f64,
    seed: u64,
) {
    let (feats, tgts) = parse_training_data(features, targets);
    let config = sci_form::ml::advanced_models::GradientBoostingConfig {
        n_estimators,
        learning_rate,
        tree_config: sci_form::ml::advanced_models::TreeConfig {
            max_depth,
            min_samples_split,
            min_samples_leaf,
            max_features: 0,
        },
        subsample,
        seed,
    };
    let model = sci_form::ml::advanced_models::train_gradient_boosting(&feats, &tgts, &config);
    let predictions: Vec<f64> = feats.iter().map(|f| model.predict(f)).collect();
    let payload = serde_json::json!({
        "model_type": "gradient_boosting",
        "n_estimators": model.n_estimators,
        "n_samples": tgts.len(),
        "train_rmse": rmse(&predictions, &tgts),
        "predictions": predictions,
    });
    println!("{}", serde_json::to_string_pretty(&payload).unwrap());
}

pub fn cmd_assign_peaks(
    frequencies: &str,
    intensities: &str,
    elements: &str,
    mode_displacements: Option<&str>,
) {
    let freqs: Vec<f64> = serde_json::from_str(frequencies).unwrap_or_else(|e| {
        fail(crate::format::bad_json("'frequencies'", frequencies, &e));
    });
    let ints: Vec<f64> = serde_json::from_str(intensities).unwrap_or_else(|e| {
        fail(crate::format::bad_json("'intensities'", intensities, &e));
    });
    let elems: Vec<u8> = serde_json::from_str(elements).unwrap_or_else(|e| {
        fail(crate::format::bad_json("'elements'", elements, &e));
    });
    // `--mode-displacements` omitted means "the modes are unknown"
    // (blocker C-M5).
    let displacements: Option<Vec<Vec<[f64; 3]>>> = mode_displacements.map(|raw| {
        serde_json::from_str(raw).unwrap_or_else(|e| {
            fail(crate::format::bad_json("'mode displacements'", raw, &e));
        })
    });
    let result = sci_form::ir::peak_assignment::assign_peaks(
        &freqs,
        &ints,
        &elems,
        displacements.as_deref(),
    );
    println!("{}", serde_json::to_string_pretty(&result).unwrap());
}

pub fn cmd_space_group(number: Option<u16>, symbol: Option<&str>) {
    let group = match (number, symbol) {
        (Some(n), _) => sci_form::materials::space_groups::space_group_by_number(n),
        (None, Some(s)) => sci_form::materials::space_groups::space_group_by_symbol(s),
        (None, None) => {
            fail("Missing '--number' or '--symbol'. Supported: '--number' (1-230) or '--symbol'.");
        }
    };
    match group {
        Some(sg) => println!("{}", serde_json::to_string_pretty(&sg).unwrap()),
        None => {
            let message = match number {
                Some(n) => format!("Invalid 'number' '{n}'. Supported: an ITC number in 1-230."),
                None => format!(
                    "Invalid 'symbol' '{}'. Supported: a Hermann-Mauguin symbol, e.g. 'P 21/c'.",
                    symbol.unwrap_or_default()
                ),
            };
            fail(message);
        }
    }
}

pub fn cmd_all_space_groups() {
    let groups: Vec<(u16, String)> = sci_form::materials::space_groups::all_space_groups()
        .into_iter()
        .map(|(n, s)| (n, s.to_string()))
        .collect();
    println!("{}", serde_json::to_string_pretty(&groups).unwrap());
}

#[allow(clippy::too_many_arguments)]
pub fn cmd_periodic_molecule(
    elements: &str,
    frac_coords: &str,
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
    bond_tolerance: Option<f64>,
) {
    let elems: Vec<u8> = serde_json::from_str(elements).unwrap_or_else(|e| {
        fail(crate::format::bad_json("'elements'", elements, &e));
    });
    let frac: Vec<[f64; 3]> = serde_json::from_str(frac_coords).unwrap_or_else(|e| {
        fail(crate::format::bad_json(
            "'fractional coords'",
            frac_coords,
            &e,
        ));
    });
    let cell = sci_form::create_unit_cell(a, b, c, alpha, beta, gamma);
    // `--bond-tolerance` omitted means "engine default" (blocker C-M5).
    let tolerance = bond_tolerance;
    let pm = sci_form::periodic::build_periodic_molecule(&elems, &frac, &cell, tolerance);
    println!("{}", serde_json::to_string_pretty(&pm).unwrap());
}

pub fn cmd_hapticity(elements: &str, coords: &str, bonds: Option<&str>) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    let bond_list = match parse_bond_list(bonds) {
        Some(b) => b,
        None => perceive_bonds(&elems, &positions, 0.3),
    };
    let result = sci_form::periodic::detect_hapticity(&elems, &positions, &bond_list);
    println!("{}", serde_json::to_string_pretty(&result).unwrap());
}

#[allow(clippy::too_many_arguments)]
pub fn cmd_optimize_framework(
    smiles: &str,
    coords: &str,
    method: &str,
    max_iter: usize,
    force_tol: f64,
    energy_tol: f64,
    max_step: f64,
    fixed_atoms: Option<&str>,
    cell: Option<[f64; 6]>,
) {
    let mol = match sci_form::parse(smiles) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("Parse error: {}", e);
            std::process::exit(1);
        }
    };
    let n = mol.graph.node_count();
    let flat: Vec<f64> = serde_json::from_str(coords).unwrap_or_else(|e| {
        fail(crate::format::bad_json("'coords'", coords, &e));
    });
    if flat.len() != n * 3 {
        eprintln!("coords length {} != 3 * atoms {}", flat.len(), n);
        std::process::exit(1);
    }
    let elems: Vec<u8> = (0..n)
        .map(|i| mol.graph[sci_form::graph::NodeIndex::new(i)].element)
        .collect();
    let initial_positions: Vec<[f64; 3]> =
        flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();

    let opt_method = match method.to_ascii_lowercase().as_str() {
        "bfgs" => sci_form::materials::geometry_opt::OptMethod::Bfgs,
        "steepest_descent" | "sd" => sci_form::materials::geometry_opt::OptMethod::SteepestDescent,
        other => {
            fail(crate::format::invalid_value(
                "'method'",
                other,
                &["bfgs", "steepest_descent", "sd"],
            ));
        }
    };
    // `--fixed-atoms` omitted means "fix no atom" (blocker C-M5).
    let fixed: Vec<usize> =
        crate::format::parse_usize_list(fixed_atoms, "'fixed atoms'").unwrap_or_default();
    let lattice =
        cell.map(|c| sci_form::create_unit_cell(c[0], c[1], c[2], c[3], c[4], c[5]).lattice);

    let config = sci_form::materials::geometry_opt::FrameworkOptConfig {
        method: opt_method,
        max_iter,
        force_tol,
        energy_tol,
        max_step,
        fixed_atoms: fixed,
        lattice,
    };

    // UFF-backed energy/force closure (UFF works in kcal/mol and kcal/mol/Å).
    let ff = sci_form::forcefield::builder::build_uff_force_field(&mol);
    const KCAL_PER_EV: f64 = 23.0605;
    let energy_force = |_els: &[u8], pos: &[[f64; 3]]| -> Result<(f64, Vec<[f64; 3]>), String> {
        let flat: Vec<f64> = pos.iter().flat_map(|p| p.iter().copied()).collect();
        let mut grad = vec![0.0f64; flat.len()];
        let energy_kcal = ff.compute_system_energy_and_gradients(&flat, &mut grad);
        let mut forces = vec![[0.0f64; 3]; pos.len()];
        for (i, f) in forces.iter_mut().enumerate() {
            f[0] = -grad[3 * i] / KCAL_PER_EV;
            f[1] = -grad[3 * i + 1] / KCAL_PER_EV;
            f[2] = -grad[3 * i + 2] / KCAL_PER_EV;
        }
        Ok((energy_kcal / KCAL_PER_EV, forces))
    };

    match sci_form::materials::geometry_opt::optimize_framework(
        &elems,
        &initial_positions,
        &energy_force,
        &config,
    ) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("Framework optimization error: {}", e);
            std::process::exit(1);
        }
    }
}

// ─── Capabilities that used to be Rust-only (blocker C-H2) ─────────────────

pub fn cmd_uhf(elements: &str, coords: &str, charge: i32, multiplicity: u32) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    let r = sci_form::compute_uhf(&elems, &positions, charge, multiplicity)
        .unwrap_or_else(|e| fail(format!("UHF error: {e}")));
    println!("{}", serde_json::to_string_pretty(&r.summary()).unwrap());
}

pub fn cmd_rohf(elements: &str, coords: &str, charge: i32, multiplicity: u32) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    let r = sci_form::compute_rohf(&elems, &positions, charge, multiplicity)
        .unwrap_or_else(|e| fail(format!("ROHF error: {e}")));
    println!("{}", serde_json::to_string_pretty(&r.summary()).unwrap());
}

pub fn cmd_pm3_gradient(elements: &str, coords: &str) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    let r = sci_form::compute_pm3_gradient(&elems, &positions)
        .unwrap_or_else(|e| fail(format!("PM3 gradient error: {e}")));
    println!("{}", serde_json::to_string_pretty(&r).unwrap());
}

pub fn cmd_xtb_gradient(elements: &str, coords: &str) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    let r = sci_form::compute_xtb_gradient(&elems, &positions)
        .unwrap_or_else(|e| fail(format!("xTB gradient error: {e}")));
    println!("{}", serde_json::to_string_pretty(&r).unwrap());
}

pub fn cmd_uff_breakdown(smiles: &str, coords: &str) {
    let flat: Vec<f64> = serde_json::from_str(coords)
        .unwrap_or_else(|e| fail(crate::format::bad_json("'coords'", coords, &e)));
    let terms = sci_form::compute_uff_energy_breakdown(smiles, &flat)
        .unwrap_or_else(|e| fail(format!("UFF breakdown error: {e}")));
    println!("{}", serde_json::to_string_pretty(&terms).unwrap());
}

pub fn cmd_mmff94_breakdown(smiles: &str, coords: &str) {
    let flat: Vec<f64> = serde_json::from_str(coords)
        .unwrap_or_else(|e| fail(crate::format::bad_json("'coords'", coords, &e)));
    let terms = sci_form::compute_mmff94_energy_breakdown(smiles, &flat)
        .unwrap_or_else(|e| fail(format!("MMFF94 breakdown error: {e}")));
    println!("{}", serde_json::to_string_pretty(&terms).unwrap());
}

pub fn cmd_tpsa(smiles: &str) {
    let tpsa = sci_form::compute_tpsa_from_smiles(smiles)
        .unwrap_or_else(|e| fail(format!("TPSA error: {e}")));
    println!("{}", serde_json::to_string_pretty(&tpsa).unwrap());
}

pub fn cmd_ensemble(smiles: &str) {
    let r = sci_form::predict_ensemble_from_smiles(smiles)
        .unwrap_or_else(|e| fail(format!("Ensemble error: {e}")));
    println!("{}", serde_json::to_string_pretty(&r).unwrap());
}

pub fn cmd_parse_cif(path: &str) {
    let text = read_input(path);
    let structure = sci_form::parse_cif(&text).unwrap_or_else(|e| fail(format!("CIF error: {e}")));
    println!("{}", serde_json::to_string_pretty(&structure).unwrap());
}

pub fn cmd_write_cif(path: &str) {
    let text = read_input(path);
    let structure: sci_form::materials::CifStructure = serde_json::from_str(&text)
        .unwrap_or_else(|e| fail(crate::format::bad_json("'structure'", &text, &e)));
    print!("{}", sci_form::write_cif(&structure));
}

/// Read a file, or standard input when the path is "-".
fn read_input(path: &str) -> String {
    use std::io::Read;
    if path == "-" {
        let mut buffer = String::new();
        std::io::stdin()
            .read_to_string(&mut buffer)
            .unwrap_or_else(|e| fail(format!("cannot read standard input: {e}")));
        return buffer;
    }
    std::fs::read_to_string(path).unwrap_or_else(|e| fail(format!("cannot read {path}: {e}")))
}
