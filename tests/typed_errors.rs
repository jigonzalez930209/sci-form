//! Public error taxonomy (2.0 migration path).

use sci_form::SciFormError;

#[test]
fn smiles_parse_failures_are_classified() {
    let error = match sci_form::parse("Q") {
        Ok(_) => panic!("Q is not an element"),
        Err(error) => error,
    };
    assert!(error.is_parse(), "expected Parse, got {error:?}");
    assert!(error.to_string().contains('Q'));
    let as_str: &str = &error;
    assert!(!as_str.is_empty());
}

#[test]
fn coordinate_length_failures_are_geometry_errors() {
    let error = sci_form::compute_mmff94_energy("CCO", &[0.0; 8])
        .expect_err("mismatched coordinates must fail");
    assert!(error.is_geometry(), "expected Geometry, got {error:?}");
    assert!(error.to_string().contains("coords length"));
}

#[test]
fn validation_and_unsupported_variants_are_available() {
    let validation = SciFormError::validation("n_samples must be > 0");
    assert!(validation.is_validation());
    let unsupported = SciFormError::unsupported("element Z=113 has no parameters");
    assert!(unsupported.is_unsupported());
    let numerical = SciFormError::numerical("SCF did not converge");
    assert!(numerical.is_numerical());
}

#[test]
fn legacy_strings_convert_and_back() {
    let error = SciFormError::from("SCF did not converge");
    assert!(error.is_numerical());
    let message: String = error.into();
    assert_eq!(message, "SCF did not converge");
}
