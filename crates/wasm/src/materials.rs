//! Materials assembly WASM bindings: unit cell, framework.

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
        _ => return crate::helpers::json_error(&format!("unknown geometry: {}", geometry)),
    };
    let topo = match topology {
        "pcu" => sci_form::materials::Topology::pcu(),
        "dia" => sci_form::materials::Topology::dia(),
        "sql" => sci_form::materials::Topology::sql(),
        _ => return crate::helpers::json_error(&format!("unknown topology: {}", topology)),
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
