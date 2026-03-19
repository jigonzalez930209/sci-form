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

pub(crate) fn serialize_or_error<T: serde::Serialize>(val: &T) -> String {
    serde_json::to_string(val).unwrap_or_else(|e| json_error(&e.to_string()))
}
