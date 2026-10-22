// PyO3 binding functions uniformly return `Result<T, PyErr>` (i.e. `PyResult<T>`)
// and convert the core crate's `String` errors with the taxonomy in
// `errors.rs` (`ValueError` for malformed input, `RuntimeError` for engine
// failures, `NotConvergedError` for non-convergence).
// Clippy's `useless_conversion` flags the identity `Result<T, PyErr>` ->
// `PyResult<T>` conversion that this pattern produces on every function, and
// several scientific entry points legitimately take more than seven arguments.
#![allow(clippy::useless_conversion)]
#![allow(clippy::too_many_arguments)]

use pyo3::prelude::*;
use pyo3::types::PyList;

mod alpha;
mod analysis;
mod beta;
mod buffer;
mod compare;
mod dynamics;
mod eht;
mod electronic;
mod embed;
mod errors;
mod experimental;
mod forcefield;
mod hf_ani_esp;
mod materials;
mod mesh;
mod ml;
mod naming;
mod nmr;
mod options;
mod population;
mod potentials;
mod properties;
mod reactivity;
mod rings;
mod smirks;
mod solvation;
mod spectroscopy;
mod stereo;
mod system;
mod transport;

/// The public API of the `sci_form` package, grouped by capability area.
///
/// `__all__` makes `dir(sci_form)` and `help(sci_form)` meaningful instead of
/// dumping an alphabetical list of 220+ names (blocker C-H5).  The
/// deprecated aliases are listed in their own group so that they are still
/// discoverable, but they are clearly marked as deprecated.
pub(crate) const __ALL__: &[&str] = &[
    // ── Package metadata ──────────────────────────────────────────────────
    "__version__",
    "version",
    "NotConvergedError",
    // ── Conformer generation ──────────────────────────────────────────────
    "embed",
    "embed_batch",
    "embed_diverse",
    "embed_coords_into",
    "parse",
    "to_canonical_smiles",
    "ConformerResult",
    "ParseResult",
    "ParsedAtom",
    // ── Charges, surface and electrostatics ───────────────────────────────
    "compute_charges",
    "compute_charges_configured",
    "compute_sasa",
    "compute_population",
    "compute_dipole",
    "compute_esp",
    "compute_esp_grid_into",
    "compute_dos",
    "compute_dos_multimethod",
    "compute_rmsd",
    "compute_rmsd_matrix",
    "ChargeResult",
    "SasaResult",
    "DipoleResult",
    "DosResult",
    "MultiMethodDosResult",
    "AlignmentResult",
    "EspGrid",
    // ── Electronic structure ──────────────────────────────────────────────
    "compute_pm3",
    "compute_xtb",
    "compute_gfn1",
    "compute_gfn2",
    "compute_hf3c",
    "compute_ani",
    "compute_cisd",
    "compute_aevs_tm",
    "compute_eht_gradient",
    "optimize_geometry_eht",
    "compute_eht_or_uff_fallback",
    "eht_support",
    "eht_calculate",
    "solve_eht",
    "eht_orbital_mesh",
    "Pm3Result",
    "XtbResult",
    "Gfn1Result",
    "Gfn2Result",
    "Hf3cResult",
    "AniResult",
    "CisdResult",
    "CisdExcitation",
    "AniTmAevs",
    "EhtResult",
    "EhtSupport",
    "EhtGradient",
    "EhtOptResult",
    // ── Force fields ──────────────────────────────────────────────────────
    "compute_uff_energy",
    "compute_uff_energy_with_aromatic_heuristics",
    "compute_mmff94_energy",
    "compute_empirical_pka",
    "EmpiricalPkaResult",
    "EmpiricalPkaSite",
    "UffHeuristicEnergy",
    // ── Molecular descriptors and machine learning ─────────────────────────
    "compute_ml_descriptors",
    "predict_ml_properties",
    "compute_whim",
    "compute_rdf",
    "compute_getaway",
    "train_random_forest",
    "train_gradient_boosting",
    "MolecularDescriptors",
    "MlPropertyResult",
    "WhimDescriptors",
    "RdfDescriptors",
    "GetawayDescriptors",
    "RegressionModel",
    // ── Spectroscopy ──────────────────────────────────────────────────────
    "compute_stda_uvvis",
    "compute_fluorescence",
    "compute_vibrational_analysis",
    "compute_ir_spectrum",
    "compute_ir_spectrum_broadened",
    "assign_peaks",
    "StdaUvVisSpectrum",
    "StdaExcitation",
    "FluorescenceSpectrum",
    "VibrationalAnalysis",
    "VibrationalMode",
    "IrSpectrum",
    "IrPeak",
    "PeakAssignment",
    "AssignmentResult",
    "FunctionalGroupMatch",
    // ── NMR ────────────────────────────────────────────────────────────────
    "predict_nmr_shifts",
    "predict_nmr_shifts_for_nucleus",
    "predict_nmr_couplings",
    "compute_ensemble_j_couplings",
    "compute_nmr_spectrum",
    "compute_giao_nmr",
    "compute_hose_codes",
    "NmrShiftResult",
    "NucleusShiftSeries",
    "ChemicalShift",
    "JCoupling",
    "NmrPeak",
    "NmrSpectrum",
    "GiaoNmrResult",
    "GiaoNmrEntry",
    // ── Solvation ──────────────────────────────────────────────────────────
    "compute_gb_solvation",
    "compute_alpb_solvation",
    "compute_nonpolar_solvation",
    "GbSolvation",
    "NonPolarSolvation",
    // ── Rings, fingerprints and clustering ────────────────────────────────
    "compute_sssr",
    "compute_ecfp",
    "compute_ecfp_batch",
    "compute_tanimoto",
    "butina_cluster",
    "filter_diverse",
    "SssrResult",
    "RingInfo",
    "ECFPFingerprint",
    "ClusterResult",
    // ── Stereochemistry, topology and reactivity ───────────────────────────
    "analyze_stereo",
    "analyze_graph_features",
    "compute_topology",
    "compute_uv_vis_spectrum",
    "compute_frontier_descriptors",
    "compute_fukui_descriptors",
    "compute_reactivity_ranking",
    "compute_bond_orders",
    "StereoAnalysis",
    "Stereocenter",
    "DoubleBondStereo",
    "GraphFeatureAnalysis",
    "MetalCoordinationCenter",
    "TopologyAnalysisResult",
    "UvVisSpectrum",
    "UvVisPeak",
    "FrontierDescriptors",
    "BondOrderResult",
    "PopulationResult",
    "FukuiDescriptors",
    "ReactivityRanking",
    "ReactivitySiteScore",
    // ── Orbital meshes ────────────────────────────────────────────────────
    "compute_orbital_mesh",
    "OrbitalMeshResult",
    // ── Materials and periodic systems ─────────────────────────────────────
    "create_unit_cell",
    "assemble_framework",
    "space_group_by_number",
    "space_group_by_symbol",
    "all_space_groups",
    "build_periodic_molecule",
    "detect_hapticity",
    "optimize_framework",
    "UnitCell",
    "CrystalStructure",
    "SymmetryOp",
    "SpaceGroup",
    "SpaceGroupEntry",
    "PeriodicBond",
    "PeriodicMolecule",
    "HapticInteraction",
    "HapticityAnalysis",
    "FrameworkOptResult",
    // ── Dynamics, reactions and transport ─────────────────────────────────
    "compute_md_trajectory",
    "compute_md_trajectory_nvt",
    "compute_simplified_neb_path",
    "search_conformers_with_uff",
    "parse_smirks",
    "apply_smirks",
    "pack_conformers",
    "split_worker_tasks",
    "estimate_workers",
    "MdTrajectory",
    "MdFrame",
    "NebPathResult",
    "NebImage",
    "ConformerSearchResult",
    "ConformerEnsembleMember",
    "ConformerClusterSummary",
    "SmirksTransform",
    "SmirksResult",
    "RecordBatch",
    // ── System planning ───────────────────────────────────────────────────
    "system_capabilities",
    "system_method_plan",
    "compare_methods",
    "SystemCapabilities",
    "MethodCapability",
    "MethodMetadata",
    "PropertyMethodPlan",
    "SystemMethodPlan",
    "MethodComparisonResult",
    "MethodComparisonEntry",
    // ── Submodules ────────────────────────────────────────────────────────
    "alpha",
    "beta",
    "potentials",
];

/// Names of the feature-gated experimental functions.
///
/// They only exist when the matching Cargo feature is enabled, so they are
/// added to `__all__` only when they are actually present in the module.
pub(crate) const __ALL_EXPERIMENTAL__: &[&str] = &[
    "compute_eeq_charges",
    "compute_eeq_energy",
    "compute_alpb_solvation",
    "compute_alpb_born_radii",
    "compute_d4_energy",
    "compute_cpm_charges",
    "compute_cpm_surface",
    "beta_compute_cpm_charges",
    "compute_edl_profile",
    "compute_edl_capacitance_scan",
    "compute_helmholtz_profile",
    "compute_gouy_chapman_profile",
    "compute_gcs_profile",
    "compute_kmesh",
    "compute_periodic_dos",
    "solve_periodic_randnla",
    "bloch_phase",
    "build_bloch_hamiltonian",
    "assemble_periodic_operators",
    "bz_integrate_scalar",
    "bz_integrate_vector",
    "validate_electron_count",
    "validate_kmesh_weights",
    "compute_htst_rate",
    "compute_htst_sweep",
    "extract_kinetics_diagnostics",
    "solve_microkinetic_network",
    "solve_microkinetic_steady_state",
    "analyze_gsm_mbh_htst_step",
    "chart_to_json",
    "chart_from_json",
    "validate_chart_schema",
    "pack_chart_payload",
    "compute_edl_chart",
    "compute_capacitance_chart",
    "compute_arrhenius_chart",
    "compute_thermal_prop_chart",
    "compute_phase_portrait_chart",
    "compute_band_structure_chart",
    "compute_dos_chart",
    "compute_kpoint_path_chart",
    "compute_fermi_surface_chart",
    "compute_trajectory_chart",
    "compute_reaction_coordinate_chart",
    "EeqChargeResult",
    "EeqEnergyResult",
    "AlpbResult",
    "AlpbBornRadii",
    "D4Result",
    "CpmResult",
    "CpmSurface",
    "EdlProfileResult",
    "KMeshResult",
    "HtstRateResult",
    "ChartPayload",
    "ChartSeries",
    "TrajectoryFrame",
];

/// sci_form Python module
#[pymodule]
fn sci_form(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__package__", "sci_form")?;
    m.add("__path__", PyList::empty(m.py()))?;
    if let Ok(spec) = m.getattr("__spec__") {
        if !spec.is_none() {
            spec.setattr("submodule_search_locations", PyList::empty(m.py()))?;
        }
    }
    // The crate version is the single source of truth for `__version__`.
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    errors::register(m)?;

    embed::register(m)?;
    system::register(m)?;
    compare::register(m)?;
    eht::register(m)?;
    properties::register(m)?;
    population::register(m)?;
    reactivity::register(m)?;
    analysis::register(m)?;
    forcefield::register(m)?;
    materials::register(m)?;
    transport::register(m)?;
    electronic::register(m)?;
    ml::register(m)?;
    spectroscopy::register(m)?;
    nmr::register(m)?;
    hf_ani_esp::register(m)?;
    mesh::register(m)?;
    dynamics::register(m)?;
    stereo::register(m)?;
    solvation::register(m)?;
    rings::register(m)?;
    smirks::register(m)?;
    experimental::register(m)?;
    alpha::register(m)?;
    beta::register(m)?;
    potentials::register(m)?;

    // The deprecated aliases are registered last so that they always point
    // at the canonical implementation.
    naming::register_aliases(m)?;

    // `__all__` drives `dir()` and `help()`; it is built from the static
    // tables above plus the deprecated aliases.  Every entry is checked
    // against the module, so a name in `__all__` always exists.
    let mut all: Vec<String> = Vec::new();
    for name in __ALL__.iter().chain(__ALL_EXPERIMENTAL__.iter()) {
        if m.hasattr(name)? {
            all.push((*name).to_string());
        }
    }
    for (scope, _, alias) in naming::ALIASES {
        if *scope == naming::Scope::Root && m.hasattr(alias)? {
            all.push((*alias).to_string());
        }
    }
    m.add("__all__", all.clone())?;

    // `__dir__` makes `dir(sci_form)` show the grouped public API.
    // It is defined in the module namespace so that it sees `__all__`.
    let globals = m.dict();
    let dir_code = "def __dir__():\n    return list(__all__)\n";
    m.py().run(
        std::ffi::CString::new(dir_code)
            .expect("dir code has no interior nul")
            .as_c_str(),
        Some(&globals),
        Some(&globals),
    )?;

    Ok(())
}
