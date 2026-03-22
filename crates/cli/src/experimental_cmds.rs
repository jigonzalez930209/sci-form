//! CLI handlers for experimental modules.

#[allow(unused_imports)]
use crate::format::parse_elems_coords;

/// EEQ geometry-dependent charges.
#[cfg(feature = "experimental-eeq")]
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
