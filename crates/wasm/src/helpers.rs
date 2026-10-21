//! Shared helper utilities for WASM bindings.

pub(crate) fn parse_elements_and_positions(
    elements: &str,
    coords_flat: &str,
) -> Result<(Vec<u8>, Vec<[f64; 3]>), String> {
    let elems: Vec<u8> =
        serde_json::from_str(elements).map_err(|e| format!("bad elements: {}", e))?;
    let flat: Vec<f64> =
        serde_json::from_str(coords_flat).map_err(|e| format!("bad coords: {}", e))?;

    if flat.len() != elems.len() * 3 {
        return Err(format!(
            "coords length {} != 3 * elements {}",
            flat.len(),
            elems.len()
        ));
    }

    let positions: Vec<[f64; 3]> = flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    Ok((elems, positions))
}

pub(crate) fn parse_flat_coords(coords_flat: &str) -> Result<Vec<f64>, String> {
    serde_json::from_str(coords_flat).map_err(|e| format!("bad coords: {}", e))
}

pub(crate) fn parse_elements(elements: &str) -> Result<Vec<u8>, String> {
    serde_json::from_str(elements).map_err(|e| format!("bad elements: {}", e))
}

pub(crate) fn json_error(msg: &str) -> String {
    format!("{{\"error\":\"{}\"}}", msg)
}

/// Perceive covalent bonds from 3D coordinates using covalent radii.
///
/// Used by the bindings that need a bond list (GETAWAY descriptors, hapticity
/// detection) or partial charges (RDF descriptors) but only receive a
/// geometry. When the caller supplies an explicit bond list it is used
/// verbatim instead.
pub(crate) fn perceive_bonds(
    elements: &[u8],
    positions: &[[f64; 3]],
    tolerance: f64,
) -> Vec<(usize, usize)> {
    let n = elements.len().min(positions.len());
    let mut bonds = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            let r_sum = sci_form::graph::get_covalent_radius(elements[i])
                + sci_form::graph::get_covalent_radius(elements[j]);
            if d > 0.4 && d < r_sum + tolerance {
                bonds.push((i, j));
            }
        }
    }
    bonds
}

/// Gasteiger-Marsili charges derived from geometry-perceived bonds.
pub(crate) fn gasteiger_charges_from_geometry(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<Vec<f64>, String> {
    let bonds = perceive_bonds(elements, positions, 0.3);
    let formal_charges = vec![0i8; elements.len()];
    let result = sci_form::charges::gasteiger::gasteiger_marsili_charges(
        elements,
        &bonds,
        &formal_charges,
        6,
    )?;
    Ok(result.charges)
}

pub(crate) fn serialize_or_error<T: serde::Serialize>(val: &T) -> String {
    serde_json::to_string(val).unwrap_or_else(|e| json_error(&e.to_string()))
}
