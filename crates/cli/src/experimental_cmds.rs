//! CLI handlers for experimental modules.
//!
//! Every handler either delegates to a real core implementation in
//! [`sci_form`] or has been removed — no handler returns hard-coded
//! placeholder data.

#[allow(unused_imports)]
use crate::format::parse_elems_coords;

fn parse_flat_coords_json(coords: &str) -> Vec<f64> {
    serde_json::from_str(coords).unwrap_or_else(|e| {
        eprintln!("Bad coords JSON: {}", e);
        std::process::exit(1);
    })
}

#[cfg(feature = "alpha-gsm")]
fn parse_gsm_methods_json(methods: &str) -> Vec<sci_form::alpha::gsm::GsmEnergyBackend> {
    let raw_methods: Vec<String> = serde_json::from_str(methods).unwrap_or_else(|e| {
        eprintln!("Bad methods JSON: {}", e);
        std::process::exit(1);
    });
    raw_methods
        .into_iter()
        .map(|method| {
            method.parse().unwrap_or_else(|e| {
                eprintln!("{}", e);
                std::process::exit(1);
            })
        })
        .collect()
}

/// EEQ geometry-dependent charges.
pub fn cmd_eeq(elements: &str, coords: &str, total_charge: f64) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    let config = sci_form::charges_eeq::EeqConfig {
        total_charge,
        regularization: 1e-10,
    };
    let r = sci_form::charges_eeq::compute_eeq_charges(&elems, &positions, &config);
    println!(
        "{}",
        serde_json::json!({
            "charges": r.charges,
            "coordination_numbers": r.coordination_numbers,
            "total_charge": r.total_charge
        })
    );
}

/// ALPB implicit solvation.
#[cfg(feature = "experimental-alpb")]
pub fn cmd_alpb(elements: &str, coords: &str, charges_str: &str, dielectric: f64) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    let charges: Vec<f64> = if charges_str.is_empty() {
        vec![0.0; elems.len()]
    } else {
        serde_json::from_str(charges_str).unwrap_or_else(|e| {
            eprintln!("bad charges JSON: {}", e);
            std::process::exit(1);
        })
    };
    let config = sci_form::solvation_alpb::AlpbConfig {
        solvent_dielectric: dielectric,
        probe_radius: 1.4,
        surface_tension: 0.005,
    };
    let r = sci_form::solvation_alpb::compute_alpb_solvation(&elems, &positions, &charges, &config);
    println!(
        "{}",
        serde_json::json!({
            "electrostatic_energy": r.electrostatic_energy,
            "nonpolar_energy": r.nonpolar_energy,
            "total_energy": r.total_energy,
            "born_radii": r.born_radii,
            "alpb_factor": r.alpb_factor
        })
    );
}

/// D4 dispersion correction.
#[cfg(feature = "experimental-d4")]
pub fn cmd_d4(elements: &str, coords: &str, three_body: bool) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    let config = sci_form::dispersion::D4Config {
        s6: 1.0,
        s8: 0.95,
        a1: 0.45,
        a2: 4.0,
        three_body,
        s9: 1.0,
    };
    let r = sci_form::dispersion::compute_d4_energy(&elems, &positions, &config);
    println!(
        "{}",
        serde_json::json!({
            "e2_body": r.e2_body,
            "e3_body": r.e3_body,
            "total_energy": r.total_energy,
            "total_kcal_mol": r.total_kcal_mol,
            "coordination_numbers": r.coordination_numbers
        })
    );
}

/// CPM charges at a given potential.
#[cfg(feature = "experimental-cpm")]
pub fn cmd_cpm(elements: &str, coords: &str, mu_ev: f64, dielectric: f64) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    let config = sci_form::beta::cpm::CpmConfig {
        mu_ev,
        dielectric,
        max_iter: 100,
        charge_tol: 1e-6,
        mu_per_atom: None,
    };
    let r = sci_form::beta::cpm::compute_cpm_charges(&elems, &positions, &config);
    println!(
        "{}",
        serde_json::json!({
            "charges": r.charges,
            "total_charge": r.total_charge,
            "grand_potential": r.grand_potential,
            "electrostatic_energy": r.electrostatic_energy,
            "mu_ev": r.mu_ev,
            "iterations": r.iterations,
            "converged": r.converged
        })
    );
}

// ─── Alpha: EDL ────────────────────────────────────────────────────────────

/// Compute an EDL profile.
#[cfg(feature = "alpha-edl")]
pub fn cmd_edl_profile(
    surface_potential_v: f64,
    model: &str,
    ionic_strength_m: f64,
    temperature_k: f64,
    n_points: usize,
    extent_angstrom: f64,
) {
    use sci_form::alpha::edl::*;
    let edl_model = match model {
        "helmholtz" => EdlModel::Helmholtz,
        "gouy-chapman" => EdlModel::GouyChapman,
        "gouy-chapman-stern" | "gcs" => EdlModel::GouyChapmanStern,
        _ => {
            eprintln!("Unknown EDL model: {}", model);
            std::process::exit(1);
        }
    };
    let config = EdlConfig {
        model: edl_model,
        temperature_k,
        ionic_strength_m,
        numerics: EdlNumerics {
            n_points,
            extent_angstrom,
        },
        ..Default::default()
    };
    match compute_edl_profile(surface_potential_v, &config) {
        Ok(r) => {
            println!(
                "{}",
                serde_json::json!({
                    "distance_axis_angstrom": r.distance_axis_angstrom,
                    "electrostatic_potential_v": r.electrostatic_potential_v,
                    "charge_density_c_per_m3": r.charge_density_c_per_m3,
                    "compact_layer_drop_v": r.compact_layer_drop_v,
                    "diffuse_layer_drop_v": r.diffuse_layer_drop_v,
                    "total_interfacial_drop_v": r.total_interfacial_drop_v,
                    "capacitance_total_f_per_m2": r.differential_capacitance.total_f_per_m2,
                    "model_name": r.model_name,
                    "converged": r.converged
                })
            );
        }
        Err(e) => {
            eprintln!("EDL error: {}", e);
            std::process::exit(1);
        }
    }
}

/// Scan EDL capacitance over a potential range.
#[cfg(feature = "alpha-edl")]
pub fn cmd_edl_capacitance_scan(
    v_min: f64,
    v_max: f64,
    n_points: usize,
    ionic_strength_m: f64,
    model: &str,
) {
    use sci_form::alpha::edl::*;
    let edl_model = match model {
        "helmholtz" => EdlModel::Helmholtz,
        "gouy-chapman" => EdlModel::GouyChapman,
        "gouy-chapman-stern" | "gcs" => EdlModel::GouyChapmanStern,
        _ => {
            eprintln!("Unknown EDL model: {}", model);
            std::process::exit(1);
        }
    };
    let config = EdlConfig {
        model: edl_model,
        ionic_strength_m,
        ..Default::default()
    };
    match scan_edl_capacitance(v_min, v_max, n_points, &config) {
        Ok(scan) => println!("{}", serde_json::to_string(&scan).unwrap()),
        Err(e) => {
            eprintln!("EDL scan error: {}", e);
            std::process::exit(1);
        }
    }
}

// ─── Alpha: Kinetics ──────────────────────────────────────────────────────

/// Evaluate an HTST transition rate.
#[cfg(feature = "alpha-kinetics")]
pub fn cmd_htst_rate(
    activation_free_energy_ev: f64,
    reaction_free_energy_ev: f64,
    temperature_k: f64,
) {
    use sci_form::alpha::kinetics::*;
    let step = ElementaryStep {
        step_id: "cli".into(),
        activation_free_energy_ev,
        reaction_free_energy_ev,
        prefactor_s_inv: None,
    };
    let state = ThermodynamicState {
        temperature_k,
        pressure_bar: 1.0,
    };
    match evaluate_htst_rate(&step, state) {
        Ok(r) => {
            println!(
                "{}",
                serde_json::json!({
                    "step_id": r.step_id,
                    "forward_rate_s_inv": r.forward_rate_s_inv,
                    "reverse_rate_s_inv": r.reverse_rate_s_inv,
                    "equilibrium_constant": r.equilibrium_constant,
                    "temperature_k": r.state.temperature_k
                })
            );
        }
        Err(e) => {
            eprintln!("HTST error: {}", e);
            std::process::exit(1);
        }
    }
}

/// Show available GSM backends for a molecule.
#[cfg(feature = "alpha-gsm")]
pub fn cmd_gsm_backend_plan(smiles: &str) {
    match sci_form::alpha::gsm::plan_gsm_backends_for_smiles(smiles) {
        Ok(plan) => println!("{}", serde_json::to_string_pretty(&plan).unwrap()),
        Err(e) => {
            eprintln!("GSM backend plan error: {}", e);
            std::process::exit(1);
        }
    }
}

/// Compare GSM backends on the same geometry.
#[cfg(feature = "alpha-gsm")]
pub fn cmd_gsm_compare_backends(smiles: &str, coords: &str, methods: &str) {
    let flat = parse_flat_coords_json(coords);
    let parsed_methods = parse_gsm_methods_json(methods);
    match sci_form::alpha::gsm::compare_gsm_backends(smiles, &flat, &parsed_methods) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("GSM backend comparison error: {}", e);
            std::process::exit(1);
        }
    }
}

/// Compute a simplified NEB-like path between two geometries.
#[allow(clippy::too_many_arguments)]
pub fn cmd_simplified_neb_path(
    smiles: &str,
    start_coords: &str,
    end_coords: &str,
    n_images: usize,
    n_iter: usize,
    spring_k: f64,
    step_size: f64,
    method: &str,
) {
    let start = parse_flat_coords_json(start_coords);
    let end = parse_flat_coords_json(end_coords);
    match sci_form::compute_simplified_neb_path_configurable(
        smiles, &start, &end, n_images, n_iter, spring_k, step_size, method,
    ) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("Simplified NEB error: {}", e);
            std::process::exit(1);
        }
    }
}

/// Compute single-point energy with any NEB-capable backend.
pub fn cmd_neb_energy(smiles: &str, coords: &str, method: &str) {
    let flat = parse_flat_coords_json(coords);
    match sci_form::neb_backend_energy_kcal(method, smiles, &flat) {
        Ok(energy) => {
            println!(
                "{}",
                serde_json::json!({
                    "method": method,
                    "energy_kcal_mol": energy,
                })
            );
        }
        Err(e) => {
            eprintln!("NEB energy error: {}", e);
            std::process::exit(1);
        }
    }
}

/// Compute energy + gradient with any NEB-capable backend.
pub fn cmd_neb_gradient(smiles: &str, coords: &str, method: &str) {
    let flat = parse_flat_coords_json(coords);
    match sci_form::neb_backend_energy_and_gradient(method, smiles, &flat) {
        Ok((energy, gradient)) => {
            println!(
                "{}",
                serde_json::json!({
                    "method": method,
                    "energy_kcal_mol": energy,
                    "gradient_kcal_mol_ang": gradient,
                })
            );
        }
        Err(e) => {
            eprintln!("NEB gradient error: {}", e);
            std::process::exit(1);
        }
    }
}

// ─── Alpha: Periodic ──────────────────────────────────────────────────────

/// Generate and print a k-mesh.
#[cfg(feature = "alpha-periodic-linear")]
pub fn cmd_kmesh(grid: [usize; 3], centering: &str) {
    use sci_form::alpha::periodic_linear::*;
    let c = match centering {
        "monkhorst-pack" | "mp" => KMeshCentering::MonkhorstPack,
        "gamma" | "gamma-centered" => KMeshCentering::GammaCentered,
        _ => {
            eprintln!("Unknown centering: {}", centering);
            std::process::exit(1);
        }
    };
    match monkhorst_pack_mesh(&KMeshConfig { grid, centering: c }) {
        Ok(mesh) => {
            println!(
                "{}",
                serde_json::json!({
                    "n_points": mesh.points.len(),
                    "grid": grid,
                    "weights_sum": mesh.points.iter().map(|p| p.weight).sum::<f64>()
                })
            );
        }
        Err(e) => {
            eprintln!("k-mesh error: {}", e);
            std::process::exit(1);
        }
    }
}

/// Compute Helmholtz-only profile.
#[cfg(feature = "alpha-edl")]
pub fn cmd_edl_helmholtz(surface_potential_v: f64, ionic_strength_m: f64, temperature_k: f64) {
    use sci_form::alpha::edl::*;
    let config = EdlConfig {
        model: EdlModel::Helmholtz,
        ionic_strength_m,
        temperature_k,
        ..Default::default()
    };
    match compute_helmholtz_profile(surface_potential_v, &config) {
        Ok(r) => {
            println!(
                "{}",
                serde_json::json!({
                    "model": "helmholtz",
                    "surface_potential_v": surface_potential_v,
                    "capacitance_f_per_m2": r.differential_capacitance.total_f_per_m2,
                    "compact_layer_drop_v": r.compact_layer_drop_v,
                    "diffuse_layer_drop_v": r.diffuse_layer_drop_v,
                    "total_drop_v": r.total_interfacial_drop_v
                })
            );
        }
        Err(e) => {
            eprintln!("Helmholtz error: {}", e);
            std::process::exit(1);
        }
    }
}

/// Compute Gouy-Chapman profile.
#[cfg(feature = "alpha-edl")]
pub fn cmd_edl_gouy_chapman(surface_potential_v: f64, ionic_strength_m: f64, temperature_k: f64) {
    use sci_form::alpha::edl::*;
    let config = EdlConfig {
        model: EdlModel::GouyChapman,
        ionic_strength_m,
        temperature_k,
        ..Default::default()
    };
    match compute_gouy_chapman_profile(surface_potential_v, &config) {
        Ok(r) => {
            println!(
                "{}",
                serde_json::json!({
                    "model": "gouy-chapman",
                    "surface_potential_v": surface_potential_v,
                    "capacitance_f_per_m2": r.differential_capacitance.total_f_per_m2,
                    "compact_layer_drop_v": r.compact_layer_drop_v,
                    "diffuse_layer_drop_v": r.diffuse_layer_drop_v,
                    "total_drop_v": r.total_interfacial_drop_v
                })
            );
        }
        Err(e) => {
            eprintln!("Gouy-Chapman error: {}", e);
            std::process::exit(1);
        }
    }
}

/// Compute Gouy-Chapman-Stern profile.
#[cfg(feature = "alpha-edl")]
pub fn cmd_edl_gcs(
    surface_potential_v: f64,
    ionic_strength_m: f64,
    temperature_k: f64,
    stern_thickness_a: f64,
) {
    use sci_form::alpha::edl::*;
    let config = EdlConfig {
        model: EdlModel::GouyChapmanStern,
        ionic_strength_m,
        temperature_k,
        stern_thickness_angstrom: stern_thickness_a,
        ..Default::default()
    };
    match compute_gcs_profile(surface_potential_v, &config) {
        Ok(r) => {
            println!(
                "{}",
                serde_json::json!({
                    "model": "gouy-chapman-stern",
                    "surface_potential_v": surface_potential_v,
                    "stern_thickness_a": stern_thickness_a,
                    "capacitance_f_per_m2": r.differential_capacitance.total_f_per_m2,
                    "compact_layer_drop_v": r.compact_layer_drop_v,
                    "diffuse_layer_drop_v": r.diffuse_layer_drop_v,
                    "total_drop_v": r.total_interfacial_drop_v
                })
            );
        }
        Err(e) => {
            eprintln!("GCS error: {}", e);
            std::process::exit(1);
        }
    }
}

// ─── Render-bridge chart output ───────────────────────────────────────────

/// Print a chart payload as JSON.
#[cfg(all(
    feature = "alpha-render-bridge",
    any(
        feature = "alpha-edl",
        feature = "alpha-kinetics",
        feature = "alpha-periodic-linear",
        feature = "alpha-gsm"
    )
))]
fn print_chart(chart: &sci_form::render_bridge::ChartPayload) {
    println!(
        "{}",
        serde_json::json!({
            "title": chart.title,
            "series_count": chart.series.len(),
            "series": chart.series.iter().map(|s| serde_json::json!({
                "series_id": s.series_id,
                "label": s.label,
                "x": s.x,
                "y": s.y,
                "x_unit": s.x_unit,
                "y_unit": s.y_unit
            })).collect::<Vec<_>>()
        })
    );
}

/// Integrate the charge density of an EDL profile into a surface charge
/// density in coulombs per square metre.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-edl"))]
fn integrated_charge_c_per_m2(profile: &sci_form::alpha::edl::EdlProfileResult) -> f64 {
    const ANGSTROM_TO_M: f64 = 1.0e-10;
    let mut integral = 0.0_f64;
    for i in 1..profile.distance_axis_angstrom.len() {
        let dx_m = (profile.distance_axis_angstrom[i] - profile.distance_axis_angstrom[i - 1])
            * ANGSTROM_TO_M;
        integral += 0.5
            * (profile.charge_density_c_per_m3[i] + profile.charge_density_c_per_m3[i - 1])
            * dx_m;
    }
    integral
}

/// Generate an EDL profile chart.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-edl"))]
pub fn cmd_edl_chart(surface_potential_v: f64, ionic_strength_m: f64) {
    use sci_form::alpha::edl::*;
    use sci_form::alpha::render_bridge::edl_profile_chart;
    let config = EdlConfig {
        model: EdlModel::GouyChapman,
        ionic_strength_m,
        ..Default::default()
    };
    match compute_gouy_chapman_profile(surface_potential_v, &config) {
        Ok(profile) => print_chart(&edl_profile_chart(&profile)),
        Err(e) => {
            eprintln!("EDL chart error: {}", e);
            std::process::exit(1);
        }
    }
}

/// Generate a capacitance scan chart.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-edl"))]
pub fn cmd_capacitance_chart(
    v_min: f64,
    v_max: f64,
    n_points: usize,
    ionic_strength_m: f64,
    model: &str,
) {
    use sci_form::alpha::edl::*;
    use sci_form::render_bridge::{ChartPayload, ChartSeries};
    let edl_model = match model {
        "helmholtz" => EdlModel::Helmholtz,
        "gouy-chapman" => EdlModel::GouyChapman,
        "gouy-chapman-stern" | "gcs" => EdlModel::GouyChapmanStern,
        _ => {
            eprintln!("Unknown EDL model: {}", model);
            std::process::exit(1);
        }
    };
    let config = EdlConfig {
        model: edl_model,
        ionic_strength_m,
        ..Default::default()
    };
    let scan = match scan_edl_capacitance(v_min, v_max, n_points, &config) {
        Ok(scan) => scan,
        Err(e) => {
            eprintln!("Capacitance chart error: {}", e);
            std::process::exit(1);
        }
    };
    let mut potentials = Vec::with_capacity(scan.len());
    let mut capacitances = Vec::with_capacity(scan.len());
    let mut charges = Vec::with_capacity(scan.len());
    for (potential, capacitance) in &scan {
        let profile = match compute_edl_profile(*potential, &config) {
            Ok(profile) => profile,
            Err(e) => {
                eprintln!("Capacitance chart error: {}", e);
                std::process::exit(1);
            }
        };
        potentials.push(*potential);
        capacitances.push(*capacitance);
        charges.push(integrated_charge_c_per_m2(&profile));
    }
    let chart = ChartPayload {
        title: "EDL Capacitance Scan".to_string(),
        series: vec![
            ChartSeries {
                series_id: "charge_vs_potential".into(),
                label: "Q(μ)".into(),
                x: potentials.clone(),
                y: charges,
                x_unit: "V".into(),
                y_unit: "C/m^2".into(),
            },
            ChartSeries {
                series_id: "capacitance".into(),
                label: "C(μ)".into(),
                x: potentials,
                y: capacitances,
                x_unit: "V".into(),
                y_unit: "F/m^2".into(),
            },
        ],
    };
    print_chart(&chart);
}

/// Generate an Arrhenius chart from temperature-dependent rates.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-kinetics"))]
pub fn cmd_arrhenius_chart(temperatures_str: &str, rates_str: &str) {
    use sci_form::render_bridge::{ChartPayload, ChartSeries};
    let temps: Vec<f64> = serde_json::from_str(temperatures_str).unwrap_or_else(|e| {
        eprintln!("Failed to parse temperatures: {}", e);
        std::process::exit(1);
    });
    let rates: Vec<f64> = serde_json::from_str(rates_str).unwrap_or_else(|e| {
        eprintln!("Failed to parse rates: {}", e);
        std::process::exit(1);
    });
    if temps.len() != rates.len() {
        eprintln!("Temperature and rate arrays must have same length");
        std::process::exit(1);
    }
    let chart = ChartPayload {
        title: "Arrhenius Plot".to_string(),
        series: vec![ChartSeries {
            series_id: "arrhenius".into(),
            label: "ln(k) vs 1000/T".into(),
            x: temps.iter().map(|t| 1000.0 / t).collect(),
            y: rates.iter().map(|r| r.ln()).collect(),
            x_unit: "1000/K".into(),
            y_unit: "ln(s⁻¹)".into(),
        }],
    };
    print_chart(&chart);
}

/// Generate a thermal-properties chart from an elementary step.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-kinetics"))]
pub fn cmd_thermal_prop_chart(
    activation_free_energy_ev: f64,
    reaction_free_energy_ev: f64,
    temperatures_str: &str,
) {
    use sci_form::alpha::kinetics::{evaluate_htst_temperature_sweep, ElementaryStep};
    use sci_form::alpha::render_bridge::arrhenius_chart;
    use sci_form::render_bridge::{ChartPayload, ChartSeries};
    let temperatures: Vec<f64> = serde_json::from_str(temperatures_str).unwrap_or_else(|e| {
        eprintln!("Failed to parse temperatures: {}", e);
        std::process::exit(1);
    });
    let step = ElementaryStep {
        step_id: "thermal".to_string(),
        activation_free_energy_ev,
        reaction_free_energy_ev,
        prefactor_s_inv: None,
    };
    let results = match evaluate_htst_temperature_sweep(&step, &temperatures, 1.0) {
        Ok(results) => results,
        Err(e) => {
            eprintln!("Thermal property chart error: {}", e);
            std::process::exit(1);
        }
    };
    let base = arrhenius_chart(&results);
    let mut series = base.series.into_iter().collect::<Vec<_>>();
    series.push(ChartSeries {
        series_id: "equilibrium_constant".into(),
        label: "ln(K_eq) vs 1000/T".into(),
        x: results
            .iter()
            .map(|r| 1000.0 / r.state.temperature_k)
            .collect(),
        y: results
            .iter()
            .map(|r| r.equilibrium_constant.ln())
            .collect(),
        x_unit: "1000/K".into(),
        y_unit: "ln".into(),
    });
    let chart = ChartPayload {
        title: "Thermal Properties".to_string(),
        series,
    };
    print_chart(&chart);
}

/// Generate a phase portrait chart from a solved microkinetic network.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-kinetics"))]
pub fn cmd_phase_portrait_chart(network: Option<&str>, n_steps: usize) {
    use sci_form::alpha::kinetics::solve_microkinetic_network;
    use sci_form::render_bridge::{ChartPayload, ChartSeries};
    let config = microkinetic_config(network, n_steps);
    if config.species.len() < 2 {
        eprintln!("A phase portrait needs a network with at least two species");
        std::process::exit(1);
    }
    let trace = match solve_microkinetic_network(&config) {
        Ok(trace) => trace,
        Err(e) => {
            eprintln!("microkinetic solve failed: {e}");
            std::process::exit(1);
        }
    };
    let mut x = Vec::with_capacity(trace.frames.len());
    let mut y = Vec::with_capacity(trace.frames.len());
    for frame in &trace.frames {
        let first = frame
            .species
            .iter()
            .find(|s| s.species_id == config.species[0])
            .map(|s| s.population)
            .unwrap_or_default();
        let second = frame
            .species
            .iter()
            .find(|s| s.species_id == config.species[1])
            .map(|s| s.population)
            .unwrap_or_default();
        x.push(first);
        y.push(second);
    }
    let chart = ChartPayload {
        title: "Phase Portrait".to_string(),
        series: vec![ChartSeries {
            series_id: "phase_portrait".into(),
            label: format!("[{}] vs [{}]", config.species[0], config.species[1]),
            x,
            y,
            x_unit: "arb.".into(),
            y_unit: "arb.".into(),
        }],
    };
    print_chart(&chart);
}

/// Generate a band structure chart from a periodic EHT calculation.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
pub fn cmd_band_structure_chart(
    elements: &str,
    coords: &str,
    lattice: &str,
    n_kpoints_per_segment: usize,
) {
    use sci_form::alpha::periodic_linear::*;
    use sci_form::alpha::render_bridge::band_structure_chart;
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    let lattice_values: Vec<f64> = serde_json::from_str(lattice).unwrap_or_else(|e| {
        eprintln!("Bad lattice JSON: {}", e);
        std::process::exit(1);
    });
    if lattice_values.len() != 9 {
        eprintln!("Lattice needs exactly nine components");
        std::process::exit(1);
    }
    let mut lattice_matrix = [[0.0_f64; 3]; 3];
    for (row, chunk) in lattice_values.chunks_exact(3).enumerate() {
        lattice_matrix[row] = [chunk[0], chunk[1], chunk[2]];
    }
    let n_electrons = match sci_form::eht::solve_eht(&elems, &positions, None) {
        Ok(result) => result.n_electrons,
        Err(e) => {
            eprintln!("EHT failed: {e}");
            std::process::exit(1);
        }
    };
    let config = BandStructureAdapterConfig {
        n_kpoints_per_segment,
        path: BandStructureAdapterConfig::default().path,
    };
    match compute_band_structure_adapter(&elems, &positions, &lattice_matrix, n_electrons, &config)
    {
        Ok(bands) => print_chart(&band_structure_chart(&bands)),
        Err(e) => {
            eprintln!("Band structure chart error: {}", e);
            std::process::exit(1);
        }
    }
}

/// Generate a DOS chart.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
pub fn cmd_dos_chart(energies_str: &str, dos_str: &str) {
    use sci_form::alpha::periodic_linear::{
        PeriodicBandEdgeSummary, PeriodicKpmDosResult, PeriodicSpectralDiagnostics,
    };
    use sci_form::alpha::render_bridge::periodic_dos_chart;
    let energies: Vec<f64> = serde_json::from_str(energies_str).unwrap_or_else(|e| {
        eprintln!("Failed to parse energies: {}", e);
        std::process::exit(1);
    });
    let dos_values: Vec<f64> = serde_json::from_str(dos_str).unwrap_or_else(|e| {
        eprintln!("Failed to parse DOS values: {}", e);
        std::process::exit(1);
    });
    let dos_result = PeriodicKpmDosResult {
        energies_ev: energies,
        total_dos: dos_values,
        kmesh: None,
        band_edges: PeriodicBandEdgeSummary::default(),
        diagnostics: PeriodicSpectralDiagnostics::default(),
    };
    print_chart(&periodic_dos_chart(&dos_result));
}

// ─── Periodic chart helpers ───────────────────────────────────────────────

/// Parse a JSON array of real-space operator blocks.
///
/// A block is `{"matrix": [[...], ...], "translation": [i, j, k]}` where the
/// matrix is square and the translation is expressed in lattice-vector units.
#[cfg(feature = "alpha-periodic-linear")]
fn parse_blocks(json_str: &str) -> Vec<(nalgebra::DMatrix<f64>, [i32; 3])> {
    let value: serde_json::Value = match serde_json::from_str(json_str) {
        Ok(value) => value,
        Err(e) => {
            eprintln!("Invalid real-space block JSON: {e}");
            std::process::exit(1);
        }
    };
    let array = match value.as_array() {
        Some(array) => array,
        None => {
            eprintln!("Real-space blocks must be a JSON array");
            std::process::exit(1);
        }
    };
    let mut blocks = Vec::with_capacity(array.len());
    for item in array {
        let matrix = match item.get("matrix").and_then(|m| m.as_array()) {
            Some(matrix) => matrix,
            None => {
                eprintln!("Each real-space block needs a 'matrix' array");
                std::process::exit(1);
            }
        };
        let translation = match item.get("translation").and_then(|t| t.as_array()) {
            Some(translation) => translation,
            None => {
                eprintln!("Each real-space block needs a 'translation' array");
                std::process::exit(1);
            }
        };
        if translation.len() != 3 {
            eprintln!("Lattice translations need exactly three components");
            std::process::exit(1);
        }
        let mut trans = [0i32; 3];
        for (slot, value) in trans.iter_mut().zip(translation.iter()) {
            *slot = match value.as_i64() {
                Some(v) => v as i32,
                None => {
                    eprintln!("Lattice translations must be integers");
                    std::process::exit(1);
                }
            };
        }
        let n = matrix.len();
        if n == 0 {
            eprintln!("Real-space matrices must be non-empty");
            std::process::exit(1);
        }
        let mut flat = Vec::with_capacity(n * n);
        for row in matrix {
            let row = match row.as_array() {
                Some(row) => row,
                None => {
                    eprintln!("Real-space matrix rows must be arrays");
                    std::process::exit(1);
                }
            };
            if row.len() != n {
                eprintln!("Real-space matrices must be square");
                std::process::exit(1);
            }
            for entry in row {
                match entry.as_f64() {
                    Some(value) => flat.push(value),
                    None => {
                        eprintln!("Real-space matrix entries must be numbers");
                        std::process::exit(1);
                    }
                }
            }
        }
        blocks.push((nalgebra::DMatrix::from_row_slice(n, n, &flat), trans));
    }
    if blocks.is_empty() {
        eprintln!("At least one real-space block is required");
        std::process::exit(1);
    }
    blocks
}

/// Borrow parsed blocks in the shape expected by the core k-mesh API.
#[cfg(feature = "alpha-periodic-linear")]
fn as_core_refs(
    blocks: &[(nalgebra::DMatrix<f64>, [i32; 3])],
) -> Vec<(&nalgebra::DMatrix<f64>, [i32; 3])> {
    blocks
        .iter()
        .map(|(matrix, trans)| (matrix, *trans))
        .collect()
}

/// Parse a k-mesh centering label.
#[cfg(feature = "alpha-periodic-linear")]
fn parse_centering(centering: &str) -> sci_form::kmesh::KMeshCentering {
    use sci_form::kmesh::KMeshCentering;
    match centering {
        "monkhorst-pack" | "mp" => KMeshCentering::MonkhorstPack,
        "gamma" | "gamma-centered" => KMeshCentering::GammaCentered,
        _ => {
            eprintln!("Unknown centering: {centering}");
            std::process::exit(1);
        }
    }
}

/// Row-major representation of a square matrix.
#[cfg(feature = "alpha-periodic-linear")]
fn matrix_rows(matrix: &nalgebra::DMatrix<f64>) -> Vec<Vec<f64>> {
    let n = matrix.nrows();
    (0..n)
        .map(|i| (0..n).map(|j| matrix[(i, j)]).collect())
        .collect()
}

/// Ascending eigenvalues of a symmetric matrix.
#[cfg(feature = "alpha-periodic-linear")]
fn symmetric_eigenvalues(matrix: &nalgebra::DMatrix<f64>) -> Vec<f64> {
    let eigen = nalgebra::SymmetricEigen::new(matrix.clone());
    let mut values: Vec<f64> = eigen.eigenvalues.iter().copied().collect();
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    values
}

/// Interpolate a high-symmetry k-path and return the cumulative path distance
/// together with the sampled fractional k-points.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
fn interpolate_path(
    path: &[([f64; 3], String)],
    n_points_per_segment: usize,
) -> (Vec<f64>, Vec<[f64; 3]>) {
    let n = n_points_per_segment.max(1);
    let mut distance = Vec::new();
    let mut points = Vec::new();
    let mut accumulated = 0.0_f64;
    let mut previous: Option<[f64; 3]> = None;
    for index in 0..path.len().saturating_sub(1) {
        let start = path[index].0;
        let end = path[index + 1].0;
        for step in 0..=n {
            if index > 0 && step == 0 {
                continue;
            }
            let t = step as f64 / n as f64;
            let mut k = [0.0_f64; 3];
            for axis in 0..3 {
                k[axis] = start[axis] + t * (end[axis] - start[axis]);
            }
            if let Some(prev) = previous {
                accumulated += ((k[0] - prev[0]).powi(2)
                    + (k[1] - prev[1]).powi(2)
                    + (k[2] - prev[2]).powi(2))
                .sqrt();
            }
            previous = Some(k);
            distance.push(accumulated);
            points.push(k);
        }
    }
    (distance, points)
}

/// Compute the Fermi-surface contour of the selected band in the `kz = 0` plane.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
fn fermi_contour(
    refs: &[(&nalgebra::DMatrix<f64>, [i32; 3])],
    fermi_energy_ev: Option<f64>,
    n_kx: usize,
    n_ky: usize,
) -> (Vec<f64>, Vec<f64>, f64) {
    let n_kx = n_kx.max(2);
    let n_ky = n_ky.max(2);
    let kx_values: Vec<f64> = (0..n_kx).map(|i| -0.5 + i as f64 / n_kx as f64).collect();
    let ky_values: Vec<f64> = (0..n_ky).map(|j| -0.5 + j as f64 / n_ky as f64).collect();

    let mut samples: Vec<Vec<Vec<f64>>> = Vec::with_capacity(n_ky);
    let mut all_energies: Vec<f64> = Vec::new();
    for &ky in &ky_values {
        let mut row = Vec::with_capacity(n_kx);
        for &kx in &kx_values {
            let hamiltonian = sci_form::kmesh::build_bloch_hamiltonian(refs, [kx, ky, 0.0]);
            let eigenvalues = symmetric_eigenvalues(&hamiltonian);
            all_energies.extend(eigenvalues.iter().copied());
            row.push(eigenvalues);
        }
        samples.push(row);
    }

    let fermi = match fermi_energy_ev {
        Some(value) => value,
        None => {
            let mut sorted = all_energies.clone();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            if sorted.is_empty() {
                0.0
            } else {
                sorted[sorted.len() / 2]
            }
        }
    };

    let band = samples
        .first()
        .and_then(|row| row.first())
        .map(|energies| {
            energies
                .iter()
                .filter(|&&e| e < fermi)
                .count()
                .min(energies.len().saturating_sub(1))
        })
        .unwrap_or(0);

    let mut contour_x = Vec::new();
    let mut contour_y = Vec::new();
    for (j, &ky) in ky_values.iter().enumerate() {
        for i in 0..n_kx.saturating_sub(1) {
            let e0 = samples[j][i].get(band).copied().unwrap_or(f64::NAN);
            let e1 = samples[j][i + 1].get(band).copied().unwrap_or(f64::NAN);
            if !e0.is_finite() || !e1.is_finite() {
                continue;
            }
            let f0 = e0 - fermi;
            let f1 = e1 - fermi;
            if f0 == 0.0 {
                contour_x.push(kx_values[i]);
                contour_y.push(ky);
                continue;
            }
            if f0 * f1 > 0.0 {
                continue;
            }
            let span = kx_values[i + 1] - kx_values[i];
            let energy_at = |x: f64| -> f64 {
                if span == 0.0 {
                    e0
                } else {
                    e0 + (x - kx_values[i]) / span * (e1 - e0)
                }
            };
            let (mut lo, mut hi) = (kx_values[i], kx_values[i + 1]);
            for _ in 0..60 {
                let mid = 0.5 * (lo + hi);
                if (energy_at(mid) - fermi) * f0 <= 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            contour_x.push(0.5 * (lo + hi));
            contour_y.push(ky);
        }
    }
    (contour_x, contour_y, fermi)
}

/// Generate a k-point path chart from the core high-symmetry path.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
pub fn cmd_kpoint_path_chart(n_points_per_segment: usize) {
    use sci_form::alpha::periodic_linear::BandStructureAdapterConfig;
    use sci_form::render_bridge::{ChartPayload, ChartSeries};
    let path = BandStructureAdapterConfig::default().path;
    let (distance, points) = interpolate_path(&path, n_points_per_segment);
    let chart = ChartPayload {
        title: "K-Point Path".to_string(),
        series: vec![
            ChartSeries {
                series_id: "k_x".into(),
                label: "k_x".into(),
                x: distance.clone(),
                y: points.iter().map(|k| k[0]).collect(),
                x_unit: "1/B".into(),
                y_unit: "1/B".into(),
            },
            ChartSeries {
                series_id: "k_y".into(),
                label: "k_y".into(),
                x: distance.clone(),
                y: points.iter().map(|k| k[1]).collect(),
                x_unit: "1/B".into(),
                y_unit: "1/B".into(),
            },
            ChartSeries {
                series_id: "k_z".into(),
                label: "k_z".into(),
                x: distance,
                y: points.iter().map(|k| k[2]).collect(),
                x_unit: "1/B".into(),
                y_unit: "1/B".into(),
            },
        ],
    };
    print_chart(&chart);
}

/// Generate a Fermi-surface slice chart from real-space operator blocks.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
pub fn cmd_fermi_surface_chart(
    blocks: &str,
    fermi_energy_ev: Option<f64>,
    n_kx: usize,
    n_ky: usize,
) {
    use sci_form::render_bridge::{ChartPayload, ChartSeries};
    let parsed = parse_blocks(blocks);
    let refs = as_core_refs(&parsed);
    let (contour_x, contour_y, fermi) = fermi_contour(&refs, fermi_energy_ev, n_kx, n_ky);
    let chart = ChartPayload {
        title: "Fermi Surface".to_string(),
        series: vec![ChartSeries {
            series_id: "fermi_contour".into(),
            label: format!("E_F = {fermi:.4} eV"),
            x: contour_x,
            y: contour_y,
            x_unit: "1/B".into(),
            y_unit: "1/B".into(),
        }],
    };
    print_chart(&chart);
}

// ─── GSM chart helpers ────────────────────────────────────────────────────

/// Run the core GSM search with the requested energy backend.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-gsm"))]
fn gsm_path(
    smiles: &str,
    reactant: &[f64],
    product: &[f64],
    method: &str,
) -> sci_form::alpha::gsm::GsmResult {
    use sci_form::alpha::gsm::{find_transition_state_with_backend, GsmConfig, GsmEnergyBackend};
    let backend: GsmEnergyBackend = match method.parse() {
        Ok(backend) => backend,
        Err(e) => {
            eprintln!("Unknown GSM backend '{method}': {e}");
            std::process::exit(1);
        }
    };
    match find_transition_state_with_backend(
        smiles,
        reactant,
        product,
        backend,
        &GsmConfig::default(),
    ) {
        Ok(result) => result,
        Err(e) => {
            eprintln!("GSM path error: {e}");
            std::process::exit(1);
        }
    }
}

/// Accumulated path length of the string nodes.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-gsm"))]
fn path_distances(path_coords: &[Vec<f64>]) -> Vec<f64> {
    let mut distances = Vec::with_capacity(path_coords.len());
    let mut accumulated = 0.0_f64;
    let mut previous: Option<&Vec<f64>> = None;
    for coords in path_coords {
        if let Some(prev) = previous {
            accumulated += prev
                .iter()
                .zip(coords.iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt();
        }
        previous = Some(coords);
        distances.push(accumulated);
    }
    distances
}

/// Generate a trajectory chart from a GSM reaction path.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-gsm"))]
pub fn cmd_trajectory_chart(
    smiles: &str,
    reactant_coords: &str,
    product_coords: &str,
    method: &str,
) {
    use sci_form::render_bridge::{ChartPayload, ChartSeries};
    let reactant = parse_flat_coords_json(reactant_coords);
    let product = parse_flat_coords_json(product_coords);
    let result = gsm_path(smiles, &reactant, &product, method);
    let distances = path_distances(&result.path_coords);
    let index: Vec<f64> = (0..result.path_energies.len()).map(|i| i as f64).collect();
    let chart = ChartPayload {
        title: "GSM Trajectory".to_string(),
        series: vec![
            ChartSeries {
                series_id: "path_energy".into(),
                label: "E(path)".into(),
                x: index.clone(),
                y: result.path_energies.clone(),
                x_unit: "node".into(),
                y_unit: "kcal/mol".into(),
            },
            ChartSeries {
                series_id: "path_distance".into(),
                label: "path length".into(),
                x: index,
                y: distances,
                x_unit: "node".into(),
                y_unit: "A".into(),
            },
        ],
    };
    print_chart(&chart);
}

/// Generate a reaction-coordinate chart from a GSM reaction path.
#[cfg(all(feature = "alpha-render-bridge", feature = "alpha-gsm"))]
pub fn cmd_reaction_coordinate_chart(
    smiles: &str,
    reactant_coords: &str,
    product_coords: &str,
    method: &str,
) {
    use sci_form::render_bridge::{ChartPayload, ChartSeries};
    let reactant = parse_flat_coords_json(reactant_coords);
    let product = parse_flat_coords_json(product_coords);
    let result = gsm_path(smiles, &reactant, &product, method);
    let chart = ChartPayload {
        title: "Reaction Coordinate".to_string(),
        series: vec![ChartSeries {
            series_id: "reaction_coordinate".into(),
            label: "E(rc)".into(),
            x: path_distances(&result.path_coords),
            y: result.path_energies.clone(),
            x_unit: "A".into(),
            y_unit: "kcal/mol".into(),
        }],
    };
    print_chart(&chart);
}

// ─── PERIODIC LINEAR COMMANDS ──────────────────────────────────────

/// Compute the k-averaged density of states with the core KPM solver.
#[cfg(feature = "alpha-periodic-linear")]
#[allow(clippy::too_many_arguments)]
pub fn cmd_compute_periodic_dos(
    hamiltonian_blocks: &str,
    grid: [usize; 3],
    centering: &str,
    order: usize,
    e_min: f64,
    e_max: f64,
    n_points: usize,
    fermi_energy_ev: Option<f64>,
) {
    use sci_form::alpha::periodic_linear::*;
    let blocks = parse_blocks(hamiltonian_blocks);
    let centering = parse_centering(centering);
    let kmesh_config = KMeshConfig { grid, centering };
    let mesh = match monkhorst_pack_mesh(&kmesh_config) {
        Ok(mesh) => mesh,
        Err(e) => {
            eprintln!("k-mesh error: {e}");
            std::process::exit(1);
        }
    };
    let hamiltonian_refs = as_core_refs(&blocks);
    let (hamiltonians, _overlaps) =
        sci_form::kmesh::assemble_periodic_operators(&hamiltonian_refs, &[], &mesh);
    let config = PeriodicKpmAdapterConfig {
        kmesh: kmesh_config,
        order,
        e_min_ev: e_min,
        e_max_ev: e_max,
        n_points,
        fermi_energy_ev,
        ..Default::default()
    };
    match compute_periodic_kpm_dos_from_operators(&hamiltonians, &config) {
        Ok(result) => {
            println!(
                "{}",
                serde_json::json!({
                    "energies_ev": result.energies_ev,
                    "total_dos": result.total_dos,
                    "n_kpoints": result.diagnostics.n_kpoints,
                    "order": order,
                    "e_min": e_min,
                    "e_max": e_max,
                    "fermi_energy_ev": fermi_energy_ev,
                    "band_edges": result.band_edges,
                    "diagnostics": result.diagnostics,
                })
            );
        }
        Err(e) => {
            eprintln!("periodic DOS error: {e}");
            std::process::exit(1);
        }
    }
}

/// Solve the periodic generalized eigenproblem with the core RandNLA solver.
#[cfg(feature = "alpha-periodic-linear")]
pub fn cmd_solve_periodic_randnla(
    hamiltonian_blocks: &str,
    overlap_blocks: &str,
    n_electrons: usize,
    grid: [usize; 3],
    centering: &str,
    sketch_size: Option<usize>,
) {
    use sci_form::alpha::periodic_linear::*;
    let hamiltonian_blocks = parse_blocks(hamiltonian_blocks);
    let overlap_blocks = parse_blocks(overlap_blocks);
    let centering = parse_centering(centering);
    let kmesh_config = KMeshConfig { grid, centering };
    let mesh = match monkhorst_pack_mesh(&kmesh_config) {
        Ok(mesh) => mesh,
        Err(e) => {
            eprintln!("k-mesh error: {e}");
            std::process::exit(1);
        }
    };
    let hamiltonian_refs = as_core_refs(&hamiltonian_blocks);
    let overlap_refs = as_core_refs(&overlap_blocks);
    let (hamiltonians, overlaps) =
        sci_form::kmesh::assemble_periodic_operators(&hamiltonian_refs, &overlap_refs, &mesh);
    let config = PeriodicRandNlaAdapterConfig {
        kmesh: kmesh_config,
        sketch_size,
        ..Default::default()
    };
    match sci_form::alpha::periodic_linear::solve_periodic_randnla(
        &hamiltonians,
        &overlaps,
        n_electrons,
        &config,
    ) {
        Ok(result) => {
            println!(
                "{}",
                serde_json::json!({
                    "n_kpoints": result.diagnostics.n_kpoints,
                    "sketch_size": sketch_size,
                    "homo_energy_ev": result.band_edges.homo_energy_ev,
                    "lumo_energy_ev": result.band_edges.lumo_energy_ev,
                    "band_gap_ev": result.band_edges.direct_gap_ev,
                    "direct_gap_ev": result.band_edges.direct_gap_ev,
                    "indirect_gap_ev": result.band_edges.indirect_gap_ev,
                    "max_residual": result.diagnostics.max_residual,
                    "used_fallback": result.diagnostics.used_fallback,
                    "kpoint_results": result.kpoint_results,
                })
            );
        }
        Err(e) => {
            eprintln!("periodic RandNLA error: {e}");
            std::process::exit(1);
        }
    }
}

/// Compute the Bloch phase factor for a lattice translation.
#[cfg(feature = "alpha-periodic-linear")]
pub fn cmd_bloch_phase(k_x: f64, k_y: f64, k_z: f64, t_x: i32, t_y: i32, t_z: i32) {
    let (cos_theta, sin_theta) = sci_form::kmesh::bloch_phase([k_x, k_y, k_z], [t_x, t_y, t_z]);
    println!(
        "{}",
        serde_json::json!({
            "k": [k_x, k_y, k_z],
            "translation": [t_x, t_y, t_z],
            "cos": cos_theta,
            "sin": sin_theta
        })
    );
}

/// Build the Bloch Hamiltonian `H(k)` from real-space operator blocks.
#[cfg(feature = "alpha-periodic-linear")]
pub fn cmd_build_bloch_hamiltonian(blocks: &str, k_x: f64, k_y: f64, k_z: f64) {
    let parsed = parse_blocks(blocks);
    let refs = as_core_refs(&parsed);
    let hamiltonian = sci_form::kmesh::build_bloch_hamiltonian(&refs, [k_x, k_y, k_z]);
    let n = hamiltonian.nrows();
    println!(
        "{}",
        serde_json::json!({
            "n_basis": n,
            "k": [k_x, k_y, k_z],
            "matrix_shape": [n, n],
            "matrix": matrix_rows(&hamiltonian),
            "eigenvalues": symmetric_eigenvalues(&hamiltonian),
        })
    );
}

/// Assemble periodic Hamiltonians and overlaps for every k-point of a mesh.
#[cfg(feature = "alpha-periodic-linear")]
pub fn cmd_assemble_periodic_operators(
    hamiltonian_blocks: &str,
    overlap_blocks: &str,
    grid: [usize; 3],
    centering: &str,
) {
    use sci_form::alpha::periodic_linear::{monkhorst_pack_mesh, KMeshConfig};
    let hamiltonian_blocks = parse_blocks(hamiltonian_blocks);
    let overlap_blocks = if overlap_blocks.trim().is_empty() {
        Vec::new()
    } else {
        parse_blocks(overlap_blocks)
    };
    let centering = parse_centering(centering);
    let mesh = match monkhorst_pack_mesh(&KMeshConfig { grid, centering }) {
        Ok(mesh) => mesh,
        Err(e) => {
            eprintln!("k-mesh error: {e}");
            std::process::exit(1);
        }
    };
    let hamiltonian_refs = as_core_refs(&hamiltonian_blocks);
    let overlap_refs = as_core_refs(&overlap_blocks);
    let (hamiltonians, overlaps) =
        sci_form::kmesh::assemble_periodic_operators(&hamiltonian_refs, &overlap_refs, &mesh);
    let n_basis = hamiltonians.first().map(|m| m.nrows()).unwrap_or_default();
    println!(
        "{}",
        serde_json::json!({
            "n_kpoints": hamiltonians.len(),
            "n_basis": n_basis,
            "assembled": hamiltonians.len() == mesh.points.len(),
            "kpoints": mesh.points.iter().map(|p| serde_json::json!({
                "fractional": p.fractional,
                "index": p.index,
                "weight": p.weight,
            })).collect::<Vec<_>>(),
            "hamiltonians": hamiltonians.iter().map(matrix_rows).collect::<Vec<_>>(),
            "overlaps": overlaps.iter().map(matrix_rows).collect::<Vec<_>>(),
        })
    );
}

/// Check that the electron count fits into the available band/k-point states.
#[cfg(feature = "alpha-periodic-linear")]
pub fn cmd_validate_electron_count(n_electrons: usize, n_bands: usize, n_kpoints: usize) {
    let valid = n_electrons <= n_bands * n_kpoints;
    println!(
        "{}",
        serde_json::json!({
            "n_electrons": n_electrons,
            "n_bands": n_bands,
            "n_kpoints": n_kpoints,
            "valid": valid
        })
    );
}

// ─── KINETICS COMMANDS ────────────────────────────────────────────

/// Parse a `MicrokineticNetworkConfig` from a file path or inline JSON, or
/// fall back to a small built-in demo network (A ⇌ B with a 0.5 eV barrier).
#[cfg(feature = "alpha-kinetics")]
fn microkinetic_config(
    network: Option<&str>,
    n_steps: usize,
) -> sci_form::alpha::kinetics::MicrokineticNetworkConfig {
    use sci_form::alpha::kinetics::{ElementaryStep, MicrokineticNetworkConfig};
    if let Some(spec) = network {
        let text = if std::path::Path::new(spec).is_file() {
            std::fs::read_to_string(spec).unwrap_or_else(|e| {
                eprintln!("cannot read network file {spec}: {e}");
                std::process::exit(1);
            })
        } else {
            spec.to_string()
        };
        match serde_json::from_str::<MicrokineticNetworkConfig>(&text) {
            Ok(mut cfg) => {
                if cfg.n_frames < 2 {
                    cfg.n_frames = n_steps.max(2);
                }
                cfg
            }
            Err(e) => {
                eprintln!("invalid network JSON: {e}");
                std::process::exit(1);
            }
        }
    } else {
        MicrokineticNetworkConfig {
            species: vec!["A".into(), "B".into()],
            steps: vec![ElementaryStep {
                step_id: "A->B".into(),
                activation_free_energy_ev: 0.5,
                reaction_free_energy_ev: -0.1,
                prefactor_s_inv: None,
            }],
            stoichiometry: vec![vec![(0, -1.0), (1, 1.0)]],
            initial_populations: vec![1.0, 0.0],
            total_time_s: 1e-3,
            n_frames: n_steps.max(2),
            state: sci_form::alpha::kinetics::ThermodynamicState {
                temperature_k: 298.15,
                pressure_bar: 1.0,
            },
        }
    }
}

#[cfg(feature = "alpha-kinetics")]
pub fn cmd_extract_kinetics_diagnostics(network: Option<&str>) {
    let config = microkinetic_config(network, 20);
    match sci_form::alpha::kinetics::solve_microkinetic_network(&config) {
        Ok(trace) => {
            let diag = sci_form::alpha::kinetics::extract_kinetics_diagnostics(&trace);
            println!(
                "{}",
                serde_json::to_value(&diag).unwrap_or_else(|_| serde_json::json!({}))
            );
        }
        Err(e) => {
            eprintln!("microkinetic solve failed: {e}");
            std::process::exit(1);
        }
    }
}

#[cfg(feature = "alpha-kinetics")]
pub fn cmd_solve_microkinetic_network(network: Option<&str>, n_steps: usize) {
    let config = microkinetic_config(network, n_steps);
    match sci_form::alpha::kinetics::solve_microkinetic_network(&config) {
        Ok(trace) => {
            let diag = sci_form::alpha::kinetics::extract_kinetics_diagnostics(&trace);
            let last = trace.frames.last();
            println!(
                "{}",
                serde_json::json!({
                    "n_steps": trace.frames.len(),
                    "n_steps_requested": n_steps,
                    "n_rates": trace.rates.len(),
                    "converged": diag.is_steady_state,
                    "final_time_s": last.map(|f| f.time_s).unwrap_or(0.0),
                    "final_populations": last.map(|f| f.species.iter().map(|s| s.population).collect::<Vec<_>>()).unwrap_or_default(),
                    "diagnostics": serde_json::to_value(&diag).unwrap_or_else(|_| serde_json::json!({})),
                })
            );
        }
        Err(e) => {
            eprintln!("microkinetic solve failed: {e}");
            std::process::exit(1);
        }
    }
}

#[cfg(feature = "alpha-kinetics")]
pub fn cmd_solve_microkinetic_steady_state(network: Option<&str>, n_steps: usize, tol: f64) {
    let config = microkinetic_config(network, n_steps);
    match sci_form::alpha::kinetics::solve_microkinetic_steady_state(&config, tol) {
        Ok(trace) => {
            let diag = sci_form::alpha::kinetics::extract_kinetics_diagnostics(&trace);
            let last = trace.frames.last();
            println!(
                "{}",
                serde_json::json!({
                    "n_steps": trace.frames.len(),
                    "n_steps_requested": n_steps,
                    "converged": diag.is_steady_state,
                    "residual": diag.max_rate_of_change,
                    "tol": tol,
                    "final_populations": last.map(|f| f.species.iter().map(|s| s.population).collect::<Vec<_>>()).unwrap_or_default(),
                    "diagnostics": serde_json::to_value(&diag).unwrap_or_else(|_| serde_json::json!({})),
                })
            );
        }
        Err(e) => {
            eprintln!("microkinetic steady-state solve failed: {e}");
            std::process::exit(1);
        }
    }
}

#[cfg(feature = "alpha-kinetics")]
#[allow(clippy::too_many_arguments)]
pub fn cmd_analyze_gsm_mbh_htst_step(
    smiles: &str,
    reactant_coords: &str,
    product_coords: &str,
    method: &str,
    step_id: &str,
    temperature_k: f64,
    pressure_bar: f64,
    n_nodes: usize,
    mbh_fd_step: f64,
) {
    let reactant = parse_flat_coords_json(reactant_coords);
    let product = parse_flat_coords_json(product_coords);
    let backend: sci_form::alpha::gsm::GsmEnergyBackend = method.parse().unwrap_or_else(|e| {
        eprintln!("{}", e);
        std::process::exit(1);
    });
    let config = sci_form::alpha::kinetics::HtstAdapterConfig {
        step_id: step_id.to_string(),
        state: sci_form::alpha::kinetics::ThermodynamicState {
            temperature_k,
            pressure_bar,
        },
        gsm: sci_form::alpha::gsm::GsmConfig {
            max_nodes: n_nodes,
            ..Default::default()
        },
        mbh_fd_step,
    };

    match sci_form::alpha::kinetics::analyze_gsm_mbh_htst_step_with_backend(
        smiles, &reactant, &product, backend, &config,
    ) {
        Ok(result) => {
            println!(
                "{}",
                serde_json::json!({
                    "backend": backend.as_str(),
                    "step_id": result.rate.step_id,
                    "temperature_k": result.rate.state.temperature_k,
                    "pressure_bar": result.rate.state.pressure_bar,
                    "forward_rate_s_inv": result.rate.forward_rate_s_inv,
                    "reverse_rate_s_inv": result.rate.reverse_rate_s_inv,
                    "equilibrium_constant": result.rate.equilibrium_constant,
                    "ts_energy_kcal_mol": result.ts_energy_kcal_mol,
                    "activation_energy_kcal_mol": result.activation_energy_kcal_mol,
                    "reverse_barrier_kcal_mol": result.reverse_barrier_kcal_mol,
                    "imaginary_frequency_cm1": result.imaginary_frequency_cm1,
                    "mbh_frequencies_cm1": result.mbh_frequencies_cm1,
                    "n_blocks": result.n_blocks,
                    "n_flexible": result.n_flexible,
                    "path_energies_kcal_mol": result.path_energies_kcal_mol,
                    "path_coords": result.path_coords,
                    "ts_coords": result.ts_coords,
                    "n_nodes": result.n_nodes,
                    "energy_evaluations": result.energy_evaluations
                })
            );
        }
        Err(e) => {
            eprintln!("GSM+MBH+HTST analysis error: {}", e);
            std::process::exit(1);
        }
    }
}

// ──── RENDER BRIDGE: chart utilities ──────────────────────────────────────

/// Parse either the chart wrapper envelope or a bare chart payload.
fn parse_chart(json_str: &str) -> Result<sci_form::render_bridge::ChartPayload, String> {
    let value: serde_json::Value =
        serde_json::from_str(json_str).map_err(|e| format!("invalid JSON: {e}"))?;
    let chart_value = value.get("data").cloned().unwrap_or(value);
    let chart: sci_form::render_bridge::ChartPayload =
        serde_json::from_value(chart_value).map_err(|e| format!("invalid chart payload: {e}"))?;
    if chart.series.iter().any(|s| s.x.len() != s.y.len()) {
        return Err("chart series must have matching x and y lengths".into());
    }
    Ok(chart)
}

/// Serialize a chart payload to the canonical wrapper JSON.
pub fn cmd_chart_to_json(json_input: String) {
    match parse_chart(&json_input) {
        Ok(chart) => {
            println!(
                "{}",
                serde_json::json!({
                    "format": "json",
                    "version": "1.0",
                    "data": chart,
                })
            );
        }
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(1);
        }
    }
}

/// Deserialize a chart payload (wrapper or bare chart) to chart JSON.
pub fn cmd_chart_from_json(json_input: String) {
    match parse_chart(&json_input) {
        Ok(chart) => println!("{}", serde_json::to_string(&chart).unwrap()),
        Err(e) => println!("{}", serde_json::json!({"error": e})),
    }
}

/// Validate the structure of a chart payload (wrapper or bare chart).
pub fn cmd_validate_chart_schema(json_input: String) {
    println!(
        "{}",
        serde_json::json!({"valid": parse_chart(&json_input).is_ok()})
    );
}
