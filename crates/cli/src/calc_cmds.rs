//! EHT, SASA, population, dipole, ESP, DOS, cell, assemble, ANI, PM3, xTB, GFN, HF-3c command handlers.

use crate::format::parse_elems_coords;

pub fn cmd_eht(elements: &str, coords: &str, k: f64) {
    let (elems, flat, _) = parse_elems_coords(elements, coords);
    if flat.len() != elems.len() * 3 {
        eprintln!(
            "coords length {} != 3 * elements {}",
            flat.len(),
            elems.len()
        );
        std::process::exit(1);
    }
    let positions: Vec<[f64; 3]> = flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    let k_opt = if k <= 0.0 { None } else { Some(k) };
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
            "{{\"origin\":{:?},\"spacing\":{},\"dims\":{:?},\"n_values\":{}}}",
            grid.origin,
            grid.spacing,
            grid.dims,
            grid.values.len()
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
            eprintln!("Unknown geometry: {}", geometry);
            std::process::exit(1);
        }
    };
    let topo = match topology {
        "pcu" => sci_form::materials::Topology::pcu(),
        "dia" => sci_form::materials::Topology::dia(),
        "sql" => sci_form::materials::Topology::sql(),
        _ => {
            eprintln!("Unknown topology: {}", topology);
            std::process::exit(1);
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
    let mut network_map = std::collections::HashMap::new();
    for &element in &[1u8, 6, 7, 8] {
        network_map.insert(element, sci_form::ani::weights::make_test_model(384));
    }
    let config = sci_form::ani::api::AniConfig {
        compute_forces: true,
        output_aevs: true,
        ..Default::default()
    };
    match sci_form::ani::api::compute_ani(&elems, &positions, &config, &network_map) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
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
        n_cis_states: 0,
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

pub fn cmd_stereo(smiles: &str, coords: &str) {
    let flat: Vec<f64> = if coords == "[]" || coords.is_empty() {
        vec![]
    } else {
        serde_json::from_str(coords).unwrap_or_else(|e| {
            eprintln!("bad coords: {}", e);
            std::process::exit(1);
        })
    };
    match sci_form::analyze_stereo(smiles, &flat) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(e) => {
            eprintln!("Stereo error: {}", e);
            std::process::exit(1);
        }
    }
}

pub fn cmd_solvation(elements: &str, coords: &str, charges: &str, probe_radius: f64) {
    let (elems, _, positions) = parse_elems_coords(elements, coords);
    if charges.is_empty() || charges == "[]" {
        // Non-polar only
        let r = sci_form::compute_nonpolar_solvation(&elems, &positions, Some(probe_radius));
        println!("{}", serde_json::to_string_pretty(&r).unwrap());
    } else {
        let q: Vec<f64> = serde_json::from_str(charges).unwrap_or_else(|e| {
            eprintln!("bad charges: {}", e);
            std::process::exit(1);
        });
        let r = sci_form::compute_gb_solvation(
            &elems,
            &positions,
            &q,
            Some(78.5),
            Some(1.0),
            Some(probe_radius),
        );
        println!("{}", serde_json::to_string_pretty(&r).unwrap());
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
