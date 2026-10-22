//! Shared helper utilities for WASM bindings.
//!
//! # Error message style
//!
//! Every user-facing error follows one template (blocker C-H9):
//!
//! ```text
//! <What failed> '<actual value>'. Supported: 'a', 'b', 'c'.
//! ```
//!
//! All literals are quoted, the offending input is always included and the
//! valid values are always listed when there is a fixed set.  Rust internal
//! names (`positions`, `unwrap`, `None`, a `file:line`) never leak into a
//! message.

/// Format a canonical "invalid value" error message.
///
/// `<what> '<value>'. Supported: 'a', 'b', 'c'.`
pub(crate) fn invalid_value(what: &str, value: &str, supported: &[&str]) -> String {
    if supported.is_empty() {
        return format!("Invalid {what} '{value}'.");
    }
    let list = supported
        .iter()
        .map(|s| format!("'{s}'"))
        .collect::<Vec<_>>()
        .join(", ");
    format!("Invalid {what} '{value}'. Supported: {list}.")
}

/// Format a canonical "malformed JSON argument" error message.
///
/// The offending fragment is quoted (truncated) so the user can see *what*
/// was wrong, and the `serde_json` diagnostic is appended.
pub(crate) fn bad_json(what: &str, value: &str, err: &serde_json::Error) -> String {
    let fragment = truncate(value, 80);
    format!("Invalid {what} JSON '{fragment}'. ({err})")
}

/// Format a canonical "malformed JSON array argument" error message.
///
/// Same as [`bad_json`] but states the expected shape.
pub(crate) fn bad_json_array(what: &str, value: &str, err: &serde_json::Error) -> String {
    let fragment = truncate(value, 80);
    format!(
        "Invalid {what} JSON '{fragment}'. Expected a JSON array of numbers. ({err})"
    )
}

/// Truncate a long string for inclusion in an error message.
pub(crate) fn truncate(value: &str, max_chars: usize) -> String {
    let trimmed = value.trim();
    if trimmed.chars().count() <= max_chars {
        return trimmed.to_string();
    }
    let head: String = trimmed.chars().take(max_chars).collect();
    format!("{head}…")
}

/// Position of the last unclosed `open` bracket in `text`, if any.
fn unclosed_position(text: &str, open: u8, close: u8) -> Option<usize> {
    let mut stack: Vec<usize> = Vec::new();
    for (i, b) in text.bytes().enumerate() {
        if b == open {
            stack.push(i);
        } else if b == close {
            stack.pop();
        }
    }
    stack.last().copied()
}

/// Enrich a SMILES parser error with the offending input and the position.
///
/// The core parser reports bare messages such as `Expected ')'`.  The
/// bindings rewrite them into the canonical style and add the position of
/// the offending bracket, which is the single most common user error in the
/// library (blocker C-H9).
pub(crate) fn smiles_error(message: &str, smiles: &str) -> String {
    let fragment = truncate(smiles, 120);
    let detail = match message {
        "Expected ')'" => match unclosed_position(smiles, b'(', b')') {
            Some(pos) => format!("expected ')' to close the '(' at position {pos}"),
            None => "expected ')'".to_string(),
        },
        "Expected ']'" => match unclosed_position(smiles, b'[', b']') {
            Some(pos) => format!("expected ']' to close the '[' at position {pos}"),
            None => "expected ']'".to_string(),
        },
        "Unexpected end of SMILES input" => "unexpected end of input".to_string(),
        other => other.to_string(),
    };
    format!("Invalid 'smiles' '{fragment}'. {detail}")
}

pub(crate) fn parse_elements_and_positions(
    elements: &str,
    coords_flat: &str,
) -> Result<(Vec<u8>, Vec<[f64; 3]>), String> {
    let elems: Vec<u8> =
        serde_json::from_str(elements).map_err(|e| bad_json_array("'elements'", elements, &e))?;
    let flat: Vec<f64> =
        serde_json::from_str(coords_flat).map_err(|e| bad_json_array("'coords'", coords_flat, &e))?;

    if flat.len() != elems.len() * 3 {
        return Err(format!(
            "Invalid 'coords' length {} for {} elements. Expected 3 values per atom ({} values).",
            flat.len(),
            elems.len(),
            elems.len() * 3
        ));
    }

    let positions: Vec<[f64; 3]> = flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    Ok((elems, positions))
}

pub(crate) fn parse_flat_coords(coords_flat: &str) -> Result<Vec<f64>, String> {
    serde_json::from_str(coords_flat).map_err(|e| bad_json_array("'coords'", coords_flat, &e))
}

pub(crate) fn parse_elements(elements: &str) -> Result<Vec<u8>, String> {
    serde_json::from_str(elements).map_err(|e| bad_json_array("'elements'", elements, &e))
}

/// Build the canonical WASM error payload: `{"error": "<message>"}`.
///
/// The message is escaped by `serde_json`, so the result **always** parses
/// with `JSON.parse()`.  The previous hand-rolled
/// `format!("{{\"error\":\"{}\"}}", msg)` produced invalid JSON whenever the
/// message contained a double quote — and `serde_json`'s own diagnostics
/// embed the offending value in quotes, so *every* parse failure of a JSON
/// argument produced a payload that could not be parsed back (blocker C-C1).
pub(crate) fn json_error(msg: &str) -> String {
    serde_json::json!({ "error": msg }).to_string()
}

/// Build the canonical WASM error payload for a list of errors.
///
/// Used where a function returns a JSON *array* of results: the array shape
/// must be preserved so that `JSON.parse()` keeps working.
pub(crate) fn json_error_list(msg: &str) -> String {
    serde_json::json!([{ "error": msg }]).to_string()
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

#[cfg(test)]
mod tests {
    use super::*;

    /// C-C1 regression: every WASM error path must emit a payload that
    /// `JSON.parse()` can read back, whatever the message contains.
    #[test]
    fn json_error_is_always_valid_json() {
        let messages = [
            "bad elements: invalid type: string \"abc\", expected a sequence at line 1 column 5",
            "plain message",
            "quote \" backslash \\ newline \n tab \t",
            "unicode ✓ Å",
            "",
        ];
        for msg in messages {
            let payload = json_error(msg);
            let parsed: serde_json::Value =
                serde_json::from_str(&payload).unwrap_or_else(|e| {
                    panic!("json_error produced invalid JSON for {msg:?}: {e}\n{payload}")
                });
            assert_eq!(
                parsed.get("error").and_then(|v| v.as_str()),
                Some(msg),
                "the error message must round-trip unchanged"
            );
        }
    }

    #[test]
    fn json_error_list_is_always_valid_json() {
        let payload = json_error_list("bad coords: \"oops\"");
        let parsed: serde_json::Value = serde_json::from_str(&payload).expect("valid JSON");
        let arr = parsed.as_array().expect("the payload must stay an array");
        assert_eq!(arr.len(), 1);
        assert_eq!(arr[0].get("error").and_then(|v| v.as_str()), Some("bad coords: \"oops\""));
    }

    /// The parser errors that the WASM helpers produce are the ones that used
    /// to break the JSON envelope; make sure the exact serde message keeps
    /// the payload parseable.
    #[test]
    fn serde_parse_failures_stay_parseable() {
        let err = parse_elements("\"abc\"").expect_err("a JSON string is not an element list");
        let payload = json_error(&err);
        let parsed: serde_json::Value =
            serde_json::from_str(&payload).expect("the serde message must not break the JSON");
        assert_eq!(parsed.get("error").and_then(|v| v.as_str()), Some(err.as_str()));
    }
}
