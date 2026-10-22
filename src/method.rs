//! Canonical parsing of the method names that the public API accepts as strings.
//!
//! # Why (blocker C-H3)
//!
//! Several public entry points take a `method: &str`. Every binding used to
//! parse that string with its own `match`, and the copies drifted apart:
//!
//! * the orbital mesh accepted nine aliases but its error message listed
//!   four,
//! * the DOS multimethod had the same `match` three times (Python, WASM,
//!   CLI),
//! * the IR broadening fell back to **Lorentzian** for anything that was not
//!   `"gaussian"`, so the typo `gauassian` silently changed the line shape
//!   instead of reporting it.
//!
//! The parsing now lives in the core, next to the enum it produces, and
//! every binding calls [`str::parse`]. The error is the same everywhere:
//!
//! ```text
//! Invalid method 'gauassian'. Supported: 'gaussian', 'lorentzian'.
//! ```
//!
//! The message style is the one the crate uses everywhere else (blocker
//! C-H9): what failed, the offending value, and the full list of accepted
//! values.

/// Build the canonical error message for an unknown method name.
///
/// `<what> '<actual value>'. Supported: 'a', 'b', 'c'.`
pub fn unknown_method(what: &str, value: &str, supported: &[&str]) -> String {
    let list = supported
        .iter()
        .map(|s| format!("'{s}'"))
        .collect::<Vec<_>>()
        .join(", ");
    format!("Invalid {what} '{value}'. Supported: {list}.")
}

/// The canonical error for an unknown method name.
pub(crate) fn invalid_method(value: &str, supported: &[&str]) -> String {
    unknown_method("method", value, supported)
}

/// Match a method name against a list of accepted spellings.
///
/// The comparison is case-insensitive; the first spelling that matches wins.
/// The error lists the *canonical* names, not the aliases.
pub(crate) fn match_method<'a>(
    value: &str,
    table: &[(&'a str, &'a str)],
    supported: &[&str],
) -> Result<&'a str, String> {
    let lowered = value.to_ascii_lowercase();
    for (alias, canonical) in table {
        if *alias == lowered {
            return Ok(canonical);
        }
    }
    Err(invalid_method(value, supported))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_error_names_the_value_and_every_supported_name() {
        let err = invalid_method("gauassian", &["gaussian", "lorentzian"]);
        assert_eq!(
            err,
            "Invalid method 'gauassian'. Supported: 'gaussian', 'lorentzian'."
        );
    }

    #[test]
    fn the_match_is_case_insensitive_and_returns_the_canonical_name() {
        let table = [("eht", "eht"), ("huckel", "eht")];
        assert_eq!(
            match_method("Huckel", &table, &["eht", "pm3"]).unwrap(),
            "eht"
        );
        assert!(match_method("dft", &table, &["eht", "pm3"]).is_err());
    }

    /// Every method enum parses its canonical names, accepts its aliases and
    /// rejects everything else with the canonical error.
    #[test]
    fn every_method_enum_rejects_an_unknown_name() {
        macro_rules! check {
            ($ty:ty, $good:expr, $alias:expr) => {{
                let good: $ty = $good.parse().expect("the canonical name must parse");
                let alias: $ty = $alias.parse().expect("the alias must parse");
                assert_eq!(good, alias);
                let err = "definitely-not-a-method"
                    .parse::<$ty>()
                    .expect_err("an unknown name must be rejected");
                assert!(
                    err.starts_with("Invalid method 'definitely-not-a-method'."),
                    "unexpected error: {err}"
                );
                assert!(
                    err.contains("Supported:"),
                    "the error must list the supported names: {err}"
                );
            }};
        }

        check!(crate::mesh::MeshMethod, "eht", "Huckel");
        check!(crate::ir::HessianMethod, "pm3", "PM3");
        check!(crate::dos::multi_method::DosMethod, "gfn2", "GFN2");
        check!(crate::reactivity::BroadeningType, "gaussian", "Gaussian");
        check!(crate::ir::BroadeningType, "lorentzian", "Lorentzian");
        check!(
            crate::materials::geometry_opt::OptMethod,
            "steepest_descent",
            "steepest-descent"
        );
        check!(crate::ml::whim::WhimWeighting, "mass", "MASS");
        check!(crate::edl::EdlModel, "gouy-chapman-stern", "gcs");
        check!(crate::dynamics::MdBackend, "xtb", "GFN0");
        check!(crate::kmesh::KMeshCentering, "monkhorst-pack", "mp");
    }

    /// The DFT module is behind the `alpha-dft` feature.
    #[cfg(feature = "alpha-dft")]
    #[test]
    fn the_dft_functional_rejects_an_unknown_name() {
        let err = "lda-plus"
            .parse::<crate::dft::DftMethod>()
            .expect_err("an unknown functional must be rejected");
        assert!(err.contains("'svwn'"), "{}", err);
        assert!(err.contains("'pbe'"), "{}", err);
    }

    /// The broadening used to fall back to the Lorentzian shape, so the typo
    /// `gauassian` silently changed the line shape.
    #[test]
    fn a_misspelled_line_shape_is_an_error_not_a_silent_fallback() {
        let err = "gauassian"
            .parse::<crate::ir::BroadeningType>()
            .expect_err("the typo must be rejected");
        assert!(err.contains("'gaussian'"), "{}", err);
        assert!(err.contains("'lorentzian'"), "{}", err);
    }
}
