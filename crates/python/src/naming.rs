//! Canonical function names and deprecated aliases for the Python surface.
//!
//! # Naming rule (blocker C-H1)
//!
//! The **canonical** name of a capability is the name of the Rust root
//! function, and every action is a `compute_*` (or `analyze_*` /
//! `predict_*`) verb.  The `*_calculate` suffix that the Python surface
//! used to carry has no precedent in RDKit, PySCF, ASE or xtb and is gone.
//!
//! The old Python names are kept as **deprecated aliases** so that no
//! existing script breaks; calling one emits a `DeprecationWarning` that
//! names the replacement.
//!
//! The Rust core is the source of truth: it does **not** grow aliases.

use pyo3::prelude::*;

/// Which module a name lives in.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) enum Scope {
    /// The `sci_form` package root.
    Root,
    /// The `sci_form.alpha` submodule.
    Alpha,
    /// The `sci_form.beta` submodule.
    Beta,
}

/// `(scope, canonical name, deprecated alias)`.
pub(crate) const ALIASES: &[(Scope, &str, &str)] = &[
    // ── Charges, surface, electrostatics ──────────────────────────────────
    (Scope::Root, "compute_charges", "charges"),
    (Scope::Root, "compute_sasa", "sasa"),
    (Scope::Root, "compute_population", "population"),
    (Scope::Root, "compute_dipole", "dipole"),
    (Scope::Root, "compute_dos", "dos"),
    (Scope::Root, "compute_bond_orders", "bond_orders"),
    (Scope::Root, "compute_frontier_descriptors", "frontier_descriptors"),
    (Scope::Root, "compute_fukui_descriptors", "fukui_descriptors"),
    (Scope::Root, "compute_reactivity_ranking", "reactivity_ranking"),
    (Scope::Root, "compute_rmsd", "rmsd"),
    (Scope::Root, "compute_esp", "esp"),
    // ── Force fields ───────────────────────────────────────────────────────
    (Scope::Root, "compute_uff_energy", "uff_energy"),
    (
        Scope::Root,
        "compute_uff_energy_with_aromatic_heuristics",
        "uff_energy_with_aromatic_heuristics",
    ),
    (Scope::Root, "compute_mmff94_energy", "mmff94_energy"),
    (Scope::Root, "compute_empirical_pka", "empirical_pka"),
    // ── Electronic structure ───────────────────────────────────────────────
    (Scope::Root, "compute_pm3", "pm3_calculate"),
    (Scope::Root, "compute_xtb", "xtb_calculate"),
    (Scope::Root, "compute_gfn1", "gfn1_calculate"),
    (Scope::Root, "compute_gfn2", "gfn2_calculate"),
    (Scope::Root, "compute_hf3c", "hf3c_calculate"),
    (Scope::Root, "compute_ani", "ani_calculate"),
    // ── Analysis ───────────────────────────────────────────────────────────
    (Scope::Root, "compute_uv_vis_spectrum", "uvvis_spectrum"),
    (Scope::Root, "analyze_graph_features", "graph_features"),
    (Scope::Root, "compute_topology", "topology_analysis"),
    (Scope::Root, "analyze_stereo", "stereo_analysis"),
    (Scope::Root, "compute_orbital_mesh", "orbital_mesh"),
    // ── Rings and fingerprints ─────────────────────────────────────────────
    (Scope::Root, "compute_sssr", "sssr"),
    (Scope::Root, "compute_ecfp", "ecfp"),
    (Scope::Root, "compute_tanimoto", "tanimoto"),
    (Scope::Root, "compute_rmsd_matrix", "rmsd_matrix"),
    // ── Machine learning ───────────────────────────────────────────────────
    (Scope::Root, "compute_ml_descriptors", "ml_descriptors"),
    (Scope::Root, "predict_ml_properties", "ml_predict"),
    // ── NMR ────────────────────────────────────────────────────────────────
    (Scope::Root, "predict_nmr_shifts", "nmr_shifts"),
    (
        Scope::Root,
        "predict_nmr_shifts_for_nucleus",
        "nmr_shifts_for_nucleus",
    ),
    (Scope::Root, "predict_nmr_couplings", "nmr_couplings"),
    (Scope::Root, "compute_nmr_spectrum", "nmr_spectrum"),
    (Scope::Root, "compute_giao_nmr", "giao_nmr"),
    (Scope::Root, "compute_hose_codes", "hose_codes"),
    // ── Spectroscopy ───────────────────────────────────────────────────────
    (Scope::Root, "compute_vibrational_analysis", "vibrational_analysis"),
    (Scope::Root, "compute_ir_spectrum", "ir_spectrum"),
    (Scope::Root, "compute_stda_uvvis", "stda_uvvis"),
    (Scope::Root, "compute_fluorescence", "fluorescence_spectrum"),
    // ── Solvation ──────────────────────────────────────────────────────────
    (Scope::Root, "compute_gb_solvation", "gb_solvation"),
    (Scope::Root, "compute_alpb_solvation", "alpb_solvation"),
    (Scope::Root, "compute_nonpolar_solvation", "nonpolar_solvation"),
    // ── EHT ────────────────────────────────────────────────────────────────
    (Scope::Root, "compute_eht_or_uff_fallback", "eht_or_uff_fallback"),
    (Scope::Root, "solve_eht", "eht_calculate"),
    // ── Materials ──────────────────────────────────────────────────────────
    (Scope::Root, "create_unit_cell", "unit_cell"),
    // ── Dynamics ───────────────────────────────────────────────────────────
    (Scope::Root, "compute_md_trajectory", "md_trajectory"),
    (Scope::Root, "compute_md_trajectory_nvt", "md_trajectory_nvt"),
    (
        Scope::Root,
        "compute_simplified_neb_path",
        "simplified_neb_path",
    ),
    (
        Scope::Root,
        "search_conformers_with_uff",
        "search_conformers",
    ),
    // ── Experimental tier ──────────────────────────────────────────────────
    (Scope::Root, "compute_eeq_charges", "eeq_charges"),
    (Scope::Root, "compute_eeq_energy", "eeq_energy"),
    (Scope::Root, "compute_alpb_born_radii", "alpb_born_radii"),
    (Scope::Root, "compute_d4_energy", "d4_energy"),
    (Scope::Root, "compute_cpm_charges", "cpm_charges"),
    (Scope::Root, "compute_cpm_surface", "cpm_surface"),
    (Scope::Root, "compute_edl_profile", "edl_profile"),
    (
        Scope::Root,
        "compute_edl_capacitance_scan",
        "edl_capacitance_scan",
    ),
    (Scope::Root, "compute_helmholtz_profile", "edl_helmholtz"),
    (Scope::Root, "compute_gouy_chapman_profile", "edl_gouy_chapman"),
    (Scope::Root, "compute_gcs_profile", "edl_gcs"),
    (Scope::Root, "compute_kmesh", "kmesh"),
    (Scope::Root, "compute_htst_rate", "htst_rate"),
    (Scope::Root, "compute_htst_sweep", "htst_temperature_sweep"),
    (Scope::Root, "compute_edl_chart", "edl_chart"),
    (Scope::Root, "compute_capacitance_chart", "capacitance_chart"),
    (Scope::Root, "compute_arrhenius_chart", "arrhenius_chart"),
    (Scope::Root, "compute_thermal_prop_chart", "thermal_prop_chart"),
    (Scope::Root, "compute_phase_portrait_chart", "phase_portrait_chart"),
    (Scope::Root, "compute_band_structure_chart", "band_structure_chart"),
    (Scope::Root, "compute_dos_chart", "density_of_states_chart"),
    (Scope::Root, "compute_kpoint_path_chart", "kpoint_path_chart"),
    (Scope::Root, "compute_fermi_surface_chart", "fermi_surface_chart"),
    (Scope::Root, "compute_trajectory_chart", "trajectory_chart"),
    (
        Scope::Root,
        "compute_reaction_coordinate_chart",
        "reaction_coordinate_chart",
    ),
    // ── Alpha tier ─────────────────────────────────────────────────────────
    (Scope::Alpha, "alpha_compute_dft", "dft_calculate"),
    (
        Scope::Alpha,
        "alpha_compute_reaxff_gradient",
        "reaxff_gradient",
    ),
    (Scope::Alpha, "alpha_compute_reaxff_energy", "reaxff_energy"),
    (Scope::Alpha, "alpha_compute_eem_charges", "eem_charges"),
    (
        Scope::Alpha,
        "alpha_rotate_dihedral_cga",
        "rotate_dihedral_cga",
    ),
    (
        Scope::Alpha,
        "alpha_refine_torsion_cga",
        "refine_torsion_cga",
    ),
    (Scope::Alpha, "alpha_gsm_interpolate", "gsm_interpolate"),
    (Scope::Alpha, "alpha_gsm_backend_plan", "gsm_backend_plan"),
    (Scope::Alpha, "alpha_gsm_compare_backends", "gsm_compare_backends"),
    (Scope::Alpha, "alpha_gsm_find_ts", "gsm_find_ts"),
    (Scope::Alpha, "alpha_sdr_embed", "sdr_embed_py"),
    // ── Beta tier ──────────────────────────────────────────────────────────
    (Scope::Beta, "beta_compute_kpm_dos", "kpm_dos"),
    (Scope::Beta, "beta_compute_kpm_mulliken", "kpm_mulliken"),
    (Scope::Beta, "beta_compute_mbh_frequencies", "mbh_frequencies"),
    (Scope::Beta, "beta_solve_eht_randnla", "eht_randnla"),
    (
        Scope::Beta,
        "beta_compute_cpm_charges",
        "cpm_charges",
    ),
];

/// Build the deprecated wrapper for `canonical` and register it as `alias`.
///
/// The wrapper forwards every argument and emits a `DeprecationWarning`
/// that names the replacement.
///
/// The wrapper is a plain Python function (not a `PyCFunction`) so that
/// `__wrapped__` and `__doc__` can be set on it: `inspect.signature()` then
/// reports the signature of the canonical function, and `help()` shows its
/// documentation.
fn add_alias(module: &Bound<'_, PyAny>, canonical: &str, alias: &str) -> PyResult<()> {
    let py = module.py();
    let message = format!(
        "sci_form.{alias} is deprecated and will be removed in sci-form 2.0; \
         use sci_form.{canonical} instead."
    );
    let code = format!(
        "def {alias}(*args, **kwargs):\n\
         \x20   import warnings\n\
         \x20   warnings.warn({message:?}, DeprecationWarning, stacklevel=2)\n\
         \x20   return {canonical}(*args, **kwargs)\n"
    );
    let dict = module.getattr("__dict__")?;
    let globals = dict.downcast::<pyo3::types::PyDict>()?;
    py.run(
        std::ffi::CString::new(code).expect("alias code has no interior nul").as_c_str(),
        Some(globals),
        Some(globals),
    )?;

    let func = module.getattr(alias)?;
    func.setattr("__doc__", format!("Deprecated alias of `sci_form.{canonical}`."))?;
    // `__wrapped__` makes `inspect.signature()` follow the alias to the
    // canonical function.
    func.setattr("__wrapped__", module.getattr(canonical)?)?;
    Ok(())
}

/// Register every deprecated alias on the package and its submodules.
///
/// `root` must be the `sci_form` module *after* every canonical function has
/// been registered, including the `alpha` and `beta` submodules.
pub(crate) fn register_aliases(root: &Bound<'_, PyModule>) -> PyResult<()> {
    let alpha = root.getattr("alpha").ok();
    let beta = root.getattr("beta").ok();
    for (scope, canonical, alias) in ALIASES {
        let module: Bound<'_, PyAny> = match scope {
            Scope::Root => root.clone().into_any(),
            Scope::Alpha => match alpha.clone() {
                Some(m) => m,
                None => continue,
            },
            Scope::Beta => match beta.clone() {
                Some(m) => m,
                None => continue,
            },
        };
        // A canonical function may be missing when the matching feature flag
        // is off; that is not an error.
        if module.hasattr(canonical)? {
            add_alias(&module, canonical, alias)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aliases_are_unique_and_point_at_canonical_names() {
        let mut seen = std::collections::HashSet::new();
        for (scope, canonical, alias) in ALIASES {
            assert!(
                seen.insert((*scope, *alias)),
                "duplicate alias {alias} in {scope:?}"
            );
            assert_ne!(*canonical, *alias, "an alias must differ from the canonical name");
            assert!(
                canonical.chars().all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_'),
                "canonical name {canonical} is not a plain identifier"
            );
        }
    }
}
