//! Materials assembly WASM bindings: unit cell, framework, space groups,
//! periodic molecules, hapticity and framework optimization.

use crate::helpers::{bad_json, json_error, perceive_bonds, serialize_or_error, smiles_error};
use wasm_bindgen::prelude::*;

/// Create a unit cell and return parameters + volume as JSON.
#[wasm_bindgen]
pub fn create_unit_cell(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> String {
    let cell = sci_form::create_unit_cell(a, b, c, alpha, beta, gamma);
    let vol = cell.volume();
    let p = cell.parameters();
    format!(
        "{{\"a\":{:.4},\"b\":{:.4},\"c\":{:.4},\"alpha\":{:.2},\"beta\":{:.2},\"gamma\":{:.2},\"volume\":{:.4}}}",
        p.a, p.b, p.c, p.alpha, p.beta, p.gamma, vol
    )
}

/// Assemble a framework crystal structure.
///
/// `topology`: "pcu", "dia", or "sql".
/// `metal`: atomic number of the metal center.
/// `geometry`: "linear", "tetrahedral", "octahedral", "square_planar", "trigonal".
/// `lattice_a`: cubic lattice parameter in Å.
/// `supercell`: replication factor (1 = no replication).
#[wasm_bindgen]
pub fn assemble_framework(
    topology: &str,
    metal: u8,
    geometry: &str,
    lattice_a: f64,
    supercell: usize,
) -> String {
    let geom = match geometry {
        "linear" => sci_form::materials::CoordinationGeometry::Linear,
        "trigonal" => sci_form::materials::CoordinationGeometry::Trigonal,
        "tetrahedral" => sci_form::materials::CoordinationGeometry::Tetrahedral,
        "square_planar" => sci_form::materials::CoordinationGeometry::SquarePlanar,
        "octahedral" => sci_form::materials::CoordinationGeometry::Octahedral,
        _ => {
            return crate::helpers::json_error(&crate::helpers::invalid_value(
                "'geometry'",
                geometry,
                &["linear", "trigonal", "tetrahedral", "square_planar", "octahedral"],
            ))
        }
    };
    let topo = match topology {
        "pcu" => sci_form::materials::Topology::pcu(),
        "dia" => sci_form::materials::Topology::dia(),
        "sql" => sci_form::materials::Topology::sql(),
        _ => {
            return crate::helpers::json_error(&crate::helpers::invalid_value(
                "'topology'",
                topology,
                &["pcu", "dia", "sql"],
            ))
        }
    };
    let node = sci_form::materials::Sbu::metal_node(metal, 0.0, geom);
    let linker = sci_form::materials::Sbu::linear_linker(&[6, 6], 1.4, "carboxylate");
    let cell = sci_form::materials::UnitCell::cubic(lattice_a);
    let mut structure = sci_form::assemble_framework(&node, &linker, &topo, &cell);
    if supercell > 1 {
        structure = structure.make_supercell(supercell, supercell, supercell);
    }
    crate::helpers::serialize_or_error(&structure)
}

/// Look up a crystallographic space group by ITC number (1–230).
#[wasm_bindgen]
pub fn space_group_by_number(number: u16) -> String {
    match sci_form::materials::space_groups::space_group_by_number(number) {
        Some(sg) => serialize_or_error(&sg),
        None => json_error(&crate::helpers::invalid_value(
            "'number'",
            &number.to_string(),
            &["1-230"],
        )),
    }
}

/// Look up a crystallographic space group by Hermann-Mauguin symbol.
#[wasm_bindgen]
pub fn space_group_by_symbol(symbol: &str) -> String {
    match sci_form::materials::space_groups::space_group_by_symbol(symbol) {
        Some(sg) => serialize_or_error(&sg),
        None => json_error(&crate::helpers::invalid_value(
            "'symbol'",
            symbol,
            &["a Hermann-Mauguin symbol, e.g. 'P 21/c'"],
        )),
    }
}

/// List all 230 space groups as `[number, symbol]` pairs.
#[wasm_bindgen]
pub fn all_space_groups() -> String {
    let groups: Vec<(u16, String)> = sci_form::materials::space_groups::all_space_groups()
        .into_iter()
        .map(|(n, s)| (n, s.to_string()))
        .collect();
    serialize_or_error(&groups)
}

/// Build a periodic molecular graph from fractional coordinates and a unit cell.
///
/// `cell_json`: JSON array `[a, b, c, alpha, beta, gamma]` (Å, degrees).
/// `bond_tolerance`: extra bonding tolerance in Å (0 = default 0.3).
#[wasm_bindgen]
pub fn build_periodic_molecule(
    elements_json: &str,
    frac_coords_json: &str,
    cell_json: &str,
    bond_tolerance: f64,
) -> String {
    let elements: Vec<u8> = match serde_json::from_str(elements_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'elements'", elements_json, &e)),
    };
    let frac_coords: Vec<[f64; 3]> = match serde_json::from_str(frac_coords_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'fractional coords'", frac_coords_json, &e)),
    };
    let cell: [f64; 6] = match serde_json::from_str(cell_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'cell'", cell_json, &e)),
    };
    let unit = sci_form::create_unit_cell(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]);
    let tolerance = if bond_tolerance > 0.0 {
        Some(bond_tolerance)
    } else {
        None
    };
    let pm = sci_form::periodic::build_periodic_molecule(&elements, &frac_coords, &unit, tolerance);
    serialize_or_error(&pm)
}

/// Detect haptic (η) metal–ring interactions and metallocene motifs.
///
/// `bonds_json`: JSON array of `[i, j]` atom-index pairs; pass `""` to
/// perceive the bonds from the geometry.
#[wasm_bindgen]
pub fn detect_hapticity(elements_json: &str, coords_flat_json: &str, bonds_json: &str) -> String {
    let (elements, positions) =
        match crate::helpers::parse_elements_and_positions(elements_json, coords_flat_json) {
            Ok(v) => v,
            Err(e) => return json_error(&e),
        };
    let bonds = if bonds_json.is_empty() {
        perceive_bonds(&elements, &positions, 0.3)
    } else {
        match serde_json::from_str::<Vec<(usize, usize)>>(bonds_json) {
            Ok(v) => v,
            Err(e) => return json_error(&bad_json("'bonds'", bonds_json, &e)),
        }
    };
    serialize_or_error(&sci_form::periodic::detect_hapticity(
        &elements, &positions, &bonds,
    ))
}

/// Optimize a molecular/framework geometry with the shared BFGS engine.
///
/// The energy/force backend is the UFF force field built from `smiles`;
/// energies are reported in eV and forces in eV/Å. `cell_json` optionally
/// enables periodic (minimum-image) optimization and is the JSON array
/// `[a, b, c, alpha, beta, gamma]`.
#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn optimize_framework(
    smiles: &str,
    coords_flat_json: &str,
    method: &str,
    max_iter: usize,
    force_tol: f64,
    energy_tol: f64,
    max_step: f64,
    fixed_atoms_json: &str,
    cell_json: &str,
) -> String {
    let mol = match sci_form::parse(smiles) {
        Ok(m) => m,
        Err(e) => return json_error(&smiles_error(&e, smiles)),
    };
    let n = mol.graph.node_count();
    let flat: Vec<f64> = match serde_json::from_str(coords_flat_json) {
        Ok(v) => v,
        Err(e) => return json_error(&bad_json("'coords'", coords_flat_json, &e)),
    };
    if flat.len() != n * 3 {
        return json_error(&format!("coords length {} != 3 * atoms {}", flat.len(), n));
    }
    let elements: Vec<u8> = (0..n)
        .map(|i| mol.graph[sci_form::graph::NodeIndex::new(i)].element)
        .collect();
    let initial_positions: Vec<[f64; 3]> =
        flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();

    let opt_method = match method.to_ascii_lowercase().as_str() {
        "bfgs" => sci_form::materials::geometry_opt::OptMethod::Bfgs,
        "steepest_descent" | "sd" => sci_form::materials::geometry_opt::OptMethod::SteepestDescent,
        other => {
            return json_error(&crate::helpers::invalid_value(
                "'method'",
                other,
                &["bfgs", "steepest_descent", "sd"],
            ))
        }
    };

    let fixed_atoms: Vec<usize> = if fixed_atoms_json.is_empty() {
        Vec::new()
    } else {
        match serde_json::from_str::<Vec<usize>>(fixed_atoms_json) {
            Ok(v) => v,
            Err(e) => return json_error(&bad_json("'fixed atoms'", fixed_atoms_json, &e)),
        }
    };

    let lattice = if cell_json.is_empty() {
        None
    } else {
        let cell: [f64; 6] = match serde_json::from_str(cell_json) {
            Ok(v) => v,
            Err(e) => return json_error(&bad_json("'cell'", cell_json, &e)),
        };
        Some(
            sci_form::create_unit_cell(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5])
                .lattice,
        )
    };

    let config = sci_form::materials::geometry_opt::FrameworkOptConfig {
        method: opt_method,
        max_iter,
        force_tol,
        energy_tol,
        max_step,
        fixed_atoms,
        lattice,
    };

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
        &elements,
        &initial_positions,
        &energy_force,
        &config,
    ) {
        Ok(result) => serialize_or_error(&result),
        Err(e) => json_error(&e),
    }
}
