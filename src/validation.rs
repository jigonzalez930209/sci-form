//! Shared input validation for geometry-carrying public entry points.
//!
//! A single `NaN` or `±Inf` coordinate is enough to break a whole family of
//! solvers: sorting eigenvalues with `partial_cmp().unwrap()` panics,
//! matrix inversions silently produce `NaN`, and the resulting "energies"
//! look plausible but are wrong.  Every public entry point that accepts
//! Cartesian coordinates therefore validates them *before* doing any work.
//!
//! The helpers below are intentionally small and dependency-free so that they
//! can be used from every layer (SCF, xTB, PM3, EHT, DFT, spectroscopy, …).
//!
//! Error-message style (kept consistent across the crate):
//!
//! * what failed,
//! * the offending value,
//! * what was expected.

/// Validate a molecular geometry given as one `[x, y, z]` triple per atom.
///
/// Checks, in order:
///
/// 1. the number of positions matches the number of elements,
/// 2. the system is not empty,
/// 3. every coordinate is finite (no `NaN`, no `±Inf`).
///
/// `elements` and `positions` may describe the same system in different
/// units (Ångström, bohr, fractional) — only the *shape* and the finiteness
/// of the coordinates are checked here.
pub(crate) fn validate_geometry(elements: &[u8], positions: &[[f64; 3]]) -> Result<(), String> {
    let n_atoms = elements.len();
    if positions.len() != n_atoms {
        return Err(format!(
            "positions has {} entries but {} atoms require {} ({} atoms x 3 coordinates)",
            positions.len(),
            n_atoms,
            n_atoms,
            n_atoms
        ));
    }
    if n_atoms == 0 {
        return Err("empty geometry: 0 atoms (at least one atom is required)".to_string());
    }
    validate_coords_finite(positions)
}

/// Validate a flat coordinate array (`[x0, y0, z0, x1, …]`) against the atom
/// count.
///
/// Mirrors [`validate_geometry`] but for the flat layout used by the
/// top-level `compute_*` helpers in `crate::lib`.
pub(crate) fn validate_flat_geometry(elements: &[u8], coords: &[f64]) -> Result<(), String> {
    let n_atoms = elements.len();
    let expected = n_atoms * 3;
    if coords.len() != expected {
        return Err(format!(
            "coords has {} values but {} atoms require {} ({} atoms x 3 coordinates)",
            coords.len(),
            n_atoms,
            expected,
            n_atoms
        ));
    }
    if n_atoms == 0 {
        return Err("empty geometry: 0 atoms (at least one atom is required)".to_string());
    }
    for (i, &v) in coords.iter().enumerate() {
        if !v.is_finite() {
            let atom = i / 3;
            let axis = ["x", "y", "z"][i % 3];
            return Err(format!(
                "atom {} has a non-finite {}-coordinate ({}) (NaN and Inf are not valid Cartesian coordinates)",
                atom, axis, v
            ));
        }
    }
    Ok(())
}

/// Validate that every component of every coordinate triple is finite.
pub(crate) fn validate_coords_finite(positions: &[[f64; 3]]) -> Result<(), String> {
    for (atom, p) in positions.iter().enumerate() {
        for (axis, &v) in p.iter().enumerate() {
            if !v.is_finite() {
                let axis_name = ["x", "y", "z"][axis];
                return Err(format!(
                    "atom {} has a non-finite {}-coordinate ({}) (NaN and Inf are not valid Cartesian coordinates)",
                    atom, axis_name, v
                ));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_validate_geometry_accepts_consistent_input() {
        assert!(validate_geometry(&[1, 8], &[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]).is_ok());
    }

    #[test]
    fn test_validate_geometry_rejects_length_mismatch() {
        let err = validate_geometry(&[1, 8], &[[0.0, 0.0, 0.0]]).unwrap_err();
        assert!(err.contains("positions has 1 entries"), "{}", err);
        assert!(err.contains("2 atoms require 2"), "{}", err);
    }

    #[test]
    fn test_validate_geometry_rejects_empty_system() {
        let err = validate_geometry(&[], &[]).unwrap_err();
        assert!(err.contains("empty geometry"), "{}", err);
    }

    #[test]
    fn test_validate_geometry_rejects_nan() {
        let err = validate_geometry(&[1, 8], &[[0.0, 0.0, 0.0], [f64::NAN, 0.0, 0.0]]).unwrap_err();
        assert!(err.contains("non-finite"), "{}", err);
        assert!(err.contains("atom 1"), "{}", err);
    }

    #[test]
    fn test_validate_geometry_rejects_inf() {
        let err =
            validate_geometry(&[1, 8], &[[0.0, 0.0, 0.0], [0.0, f64::INFINITY, 0.0]]).unwrap_err();
        assert!(err.contains("non-finite"), "{}", err);
        assert!(err.contains("y-coordinate"), "{}", err);
    }

    #[test]
    fn test_validate_flat_geometry_message_matches_documented_style() {
        let err = validate_flat_geometry(&[1, 1, 1], &[0.0; 8]).unwrap_err();
        assert_eq!(
            err,
            "coords has 8 values but 3 atoms require 9 (3 atoms x 3 coordinates)"
        );
    }

    #[test]
    fn test_validate_flat_geometry_rejects_nan() {
        let err =
            validate_flat_geometry(&[1, 1], &[0.0, f64::NAN, 0.0, 0.0, 0.0, 0.0]).unwrap_err();
        assert!(err.contains("atom 0"), "{}", err);
        assert!(err.contains("y-coordinate"), "{}", err);
    }
}
