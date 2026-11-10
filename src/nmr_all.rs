//! All-nuclide NMR registry and multi-nucleus driver.
//!
//! The existing [`crate::nmr::NmrNucleus`] registry enumerates 85 NMR-active
//! isotopes. This module exposes the per-element query and runs the existing
//! GIAO route once per element present in the molecule. It is a coverage and
//! orchestration layer, not a new shielding theory: elements whose GIAO basis
//! or SCF path is unavailable are reported as failures instead of silently
//! returning a zero.
//!
//! # Units
//!
//! Chemical shifts are in ppm; positions are Å. The GIAO result carries the
//! isotope label and reference convention.

use serde::{Deserialize, Serialize};

use crate::error::SciFormError;
use crate::nmr::NmrNucleus;

/// One NMR-active isotope known to the registry.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NmrIsotopeInfo {
    /// Atomic number.
    pub element: u8,
    /// Mass number.
    pub mass_number: u16,
    /// Canonical label, e.g. `"13C"`.
    pub label: String,
}

/// One element whose GIAO calculation failed.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NmrFailure {
    /// Atomic number.
    pub element: u8,
    /// Isotope label that was attempted.
    pub label: String,
    /// Error message.
    pub error: String,
}

/// Result of the multi-nucleus driver.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AllNmrResult {
    /// Successful GIAO results, one per element.
    pub nuclei: Vec<crate::GiaoNmrResult>,
    /// Elements that could not be evaluated.
    pub failures: Vec<NmrFailure>,
    /// Notes and caveats.
    pub notes: Vec<String>,
}

/// Return every NMR-active isotope known for an element.
pub fn active_isotopes(element: u8) -> Vec<NmrIsotopeInfo> {
    NmrNucleus::ALL
        .iter()
        .filter(|nucleus| nucleus.atomic_number() == element)
        .map(|nucleus| NmrIsotopeInfo {
            element,
            mass_number: nucleus.mass_number(),
            label: nucleus.canonical().to_string(),
        })
        .collect()
}

/// Return every element with at least one registered NMR-active isotope.
pub fn active_elements() -> Vec<u8> {
    let mut elements: Vec<u8> = NmrNucleus::ALL
        .iter()
        .map(|nucleus| nucleus.atomic_number())
        .collect();
    elements.sort_unstable();
    elements.dedup();
    elements
}

/// Run the existing GIAO route for the default isotope of every element present.
///
/// The current GIAO route is closed-shell and STO-3G. An element without a
/// supported basis is recorded in `failures`; the function only returns an error
/// when no element can be evaluated at all.
pub fn compute_all_nmr(
    elements: &[u8],
    positions_angstrom: &[[f64; 3]],
    config: &crate::GiaoNmrConfig,
) -> Result<AllNmrResult, SciFormError> {
    crate::validation::validate_geometry(elements, positions_angstrom)?;
    let mut unique: Vec<u8> = elements.to_vec();
    unique.sort_unstable();
    unique.dedup();
    let mut nuclei = Vec::new();
    let mut failures = Vec::new();

    for element in unique {
        let Some(nucleus) = NmrNucleus::default_for_element(element) else {
            failures.push(NmrFailure {
                element,
                label: format!("Z={element}"),
                error: "no default isotope registered".to_string(),
            });
            continue;
        };
        let label = nucleus.canonical().to_string();
        match crate::compute_giao_nmr_configured(elements, positions_angstrom, &label, config) {
            Ok(result) => nuclei.push(result),
            Err(error) => failures.push(NmrFailure {
                element,
                label,
                error: error.to_string(),
            }),
        }
    }

    if nuclei.is_empty() {
        return Err(SciFormError::unsupported(format!(
            "no element could be evaluated by the current GIAO route: {failures:?}"
        )));
    }

    let mut notes = vec![
        "One GIAO calculation per element; the current route is closed-shell STO-3G.".to_string(),
        "Failures are explicit; an unsupported element is never reported as zero shift."
            .to_string(),
    ];
    if !failures.is_empty() {
        notes.push(format!(
            "{} element(s) could not be evaluated; see `failures`.",
            failures.len()
        ));
    }

    Ok(AllNmrResult {
        nuclei,
        failures,
        notes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn carbon_registry_contains_c13() {
        let isotopes = active_isotopes(6);
        assert!(isotopes.iter().any(|isotope| isotope.label == "13C"));
    }

    #[test]
    fn active_elements_are_sorted_and_unique() {
        let elements = active_elements();
        assert!(elements.windows(2).all(|pair| pair[0] < pair[1]));
        assert!(elements.contains(&1));
        assert!(elements.contains(&6));
    }

    #[test]
    fn water_driver_reports_hydrogen_or_explicit_failure() {
        let water = [
            [0.0, 0.0, 0.0],
            [0.9572, 0.0, 0.0],
            [-0.239_987, 0.927_297, 0.0],
        ];
        let result = compute_all_nmr(&[8, 1, 1], &water, &crate::GiaoNmrConfig::default());
        match result {
            Ok(result) => assert!(!result.nuclei.is_empty()),
            Err(error) => assert!(error.is_unsupported()),
        }
    }

    #[test]
    fn invalid_geometry_is_rejected() {
        assert!(compute_all_nmr(
            &[1],
            &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            &crate::GiaoNmrConfig::default()
        )
        .is_err());
    }
}
