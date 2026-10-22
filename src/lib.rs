//! **sci-form** — computational chemistry library for conformer generation, quantum methods,
//! spectroscopy, molecular properties, and materials science.
//!
//! # What this library is (and is not)
//!
//! sci-form is a **fast, approximate** computational-chemistry toolkit. Its centre of gravity is
//! ETKDG-style 3D conformer generation (see [`conformer`], [`distgeom`], [`etkdg`]) plus a
//! collection of cheap electronic-structure and property models that run in milliseconds.
//!
//! It is **not** a general-purpose quantum-chemistry package:
//!
//! - **EHT** ([`eht`]) is a Hoffmann extended-Hückel treatment: valence-orbital, no electron
//!   correlation, no basis-set superposition, parameterised per element. It is a qualitative
//!   orbital picture (HOMO/LUMO gaps, orbital shapes), not a quantitative electronic-structure
//!   method.
//! - **PM3** ([`pm3`]) and **xTB** ([`xtb`], GFN0/GFN1/GFN2) are **semi-empirical / tight-binding
//!   approximations** with fitted parameters. Their energies are on their own internal scales
//!   (eV, Hartree) and are not directly comparable to experiment or to a correlated ab initio
//!   method.
//! - **UFF** ([`forcefield`]) and **MMFF94** ([`forcefield::mmff94`]) are **classical force
//!   fields** for geometry and energy ranking. They have no electronic structure at all and
//!   cannot produce orbitals, charges or spectra.
//! - **ANI** ([`ani`]) is a neural-network potential whose bundled weights are test weights, not
//!   a trained production model.
//!
//! For quantitative energetics, use a correlated wave-function or DFT code. sci-form's value is
//! fast screening, geometry generation, and property estimation at scale.
//!
//! # Coordinate convention
//!
//! All coordinate arrays are **flat** `[x0, y0, z0, x1, y1, z1, ...]` in ångströms (Å).
//! The alternative "matrix" form is one `[x, y, z]` triple per atom, in Å, and is used by the
//! electronic-structure entry points.
//!
//! Element arrays are atomic numbers (`u8`): 1 = H, 6 = C, 7 = N, 8 = O, etc.
//! Atom indices are 0-based and refer to the order of atoms in the parsed molecule / element
//! array.
//!
//! # Unit conventions by method
//!
//! Units are not uniform across the crate, so check the doc of the specific function:
//!
//! | Quantity | Unit |
//! |---|---|
//! | Geometry (coordinates, bond lengths, cell parameters) | Å (cell angles in degrees) |
//! | UFF / MMFF94 energies | kcal/mol |
//! | UFF / MMFF94 gradients | kcal/mol/Å |
//! | EHT, PM3, xTB orbital energies and gaps | eV |
//! | PM3 / xTB / HF total energies | Hartree (eV where stated per field) |
//! | PM3 / xTB gradients | eV/Å |
//! | Vibrational frequencies, IR wavenumbers | cm⁻¹ |
//! | NMR chemical shifts and shieldings | ppm |
//! | NMR J-couplings | Hz |
//! | Dipole moment | Debye |
//! | Partial charges | elementary charge (e) |
//! | Electrostatic potential | Hartree/e |
//! | DOS/PDOS | 1/eV |
//! | Surface area | Å² |
//! | Pressure | bar |
//!
//! Energies are always **total system energies**, not per-atom normalised, unless the field
//! name says otherwise. Gradients are `dE/dR` with `R` in Å, so they carry the sign convention
//! "force = −gradient".
//!
//! # Architecture
//!
//! The crate is organized in layers from molecular construction up to derived properties:
//!
//! - **Parsing & topology** — [`smiles`], [`graph`], [`smarts`], [`smirks`], [`rings`], [`topology`]
//! - **Geometry** — [`distgeom`], [`etkdg`], [`conformer`], [`alignment`], [`forcefield`], [`optimization`]
//! - **Charges & surface** — [`charges`], [`charges_eeq`], [`dipole`], [`surface`], [`esp`], [`population`]
//! - **Electrochemistry** — [`cpm`]
//! - **Solvation** — [`solvation`], [`solvation_alpb`]
//! - **Electronic structure** — [`scf`], [`eht`], [`pm3`], [`xtb`], [`hf`], [`dos`]
//! - **Spectroscopy** — [`spectroscopy`], [`ir`], [`nmr`]
//! - **Neural potentials** — [`ani`]
//! - **ML & descriptors** — [`ml`], [`reactivity`]
//! - **Materials** — [`materials`], [`periodic`]
//! - **Transport & infra** — [`transport`], [`gpu`], [`mesh`]
//!
//! # How do I compute a property?
//!
//! | I want | Call | Result |
//! |---|---|---|
//! | A 3D conformer from SMILES | [`embed`] | flat Å coordinates |
//! | Many conformers, deduplicated | [`embed_diverse`], [`search_conformers_with_uff`] | ranked conformers + RMSD clusters |
//! | A batch of conformers in parallel | [`embed_batch`] | one result per SMILES |
//! | Partial charges | [`compute_charges`] | Gasteiger charges (e) |
//! | Force-field energy | [`compute_uff_energy`], [`compute_mmff94_energy`] | kcal/mol |
//! | Per-term energy breakdown | [`compute_uff_energy_breakdown`] | kcal/mol per term kind |
//! | Orbital energies / HOMO-LUMO gap | [`compute_population`], [`compute_gfn1`], [`compute_gfn2`] | eV |
//! | Dipole moment | [`compute_dipole`] | Debye |
//! | Frontier / Fukui reactivity | [`compute_frontier_descriptors`], [`compute_fukui_descriptors`] | per-atom descriptors |
//! | Electrostatic potential grid | [`compute_esp`] | Hartree/e on a grid |
//! | Density of states | [`compute_dos`] | 1/eV vs eV |
//! | Vibrational frequencies + IR | [`compute_vibrational_analysis`], [`compute_ir_spectrum`] | cm⁻¹ |
//! | NMR shifts and spectrum | [`predict_nmr_shifts`], [`compute_nmr_spectrum`], [`compute_giao_nmr`] | ppm |
//! | Solvation free energy | [`compute_gb_solvation`], [`compute_alpb_solvation`] | kcal/mol |
//! | Molecular volume / surface | [`compute_sasa`], [`materials`] | Å², Å³ |
//! | Crystallographic structure | [`parse_cif`], [`create_unit_cell`] | cell + sites |
//! | RMSD between two structures | [`compute_rmsd`] | Å |
//! | A molecular dynamics trajectory | [`compute_md_trajectory`] | frames + energies |
//! | A reaction path | [`compute_simplified_neb_path`] | images + energies |
//! | Is this method valid for my elements? | [`system_capabilities`], [`system_method_plan`] | support metadata |
//!
//! # Feature flags
//!
//! The default feature set is `parallel` (rayon-backed batch paths). The remaining flags gate
//! **experimental** modules that are not part of the stable surface:
//!
//! | Feature | Gates |
//! |---|---|
//! | `parallel` | rayon-parallel batch APIs (default **on**) |
//! | `experimental-gpu`, `experimental-gpu-rendering` | the `gpu` compute / mesh rendering backends |
//! | `alpha-cga`, `alpha-gsm`, `alpha-sdr`, `alpha-edl`, `alpha-periodic-linear`, `alpha-kinetics`, `alpha-render-bridge`, `alpha-reaction-dynamics` | the `alpha` research modules |
//! | `beta-kpm`, `beta-mbh`, `beta-cpm`, `beta-randnla`, `beta-riemannian` | the `beta` modules |
//! | `alpha-dft` | the `dft` module |
//! | `alpha-dynamics-live` | the `dynamics_live` module |
//! | `alpha-mlff` | the `mlff` re-export module |
//!
//! The legacy `experimental-*` flags (e.g. `experimental-riemannian`) are kept as no-op
//! aliases of the corresponding `beta-*` / `alpha-*` flags so existing build scripts keep
//! resolving.
//!
//! # Targets
//!
//! The public surface is mirrored across four targets: Rust (this crate), Python (`crates/python/`),
//! WASM/TypeScript (`crates/wasm/`), and CLI (`crates/cli/`). Changes to public types or functions
//! must be propagated to all targets.
//!
//! # Example
//!
//! ```
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // Embed ethanol and report the C-C bond length in Å.
//! let result = sci_form::embed("CCO", 42);
//! assert!(result.error.is_none());
//!
//! let coords = &result.coords;
//! let dx = coords[3] - coords[0];
//! let dy = coords[4] - coords[1];
//! let dz = coords[5] - coords[2];
//! let bond = (dx * dx + dy * dy + dz * dz).sqrt();
//! assert!(bond > 1.0 && bond < 2.0);
//! # Ok(())
//! # }
//! ```

// Scientific/numerical code patterns that are idiomatic in this domain
#![allow(clippy::too_many_arguments)]
#![allow(clippy::needless_range_loop)]

pub mod alignment;
pub mod ani;
pub mod canonical_smiles;
pub mod cga;
pub mod charges;
pub mod charges_eeq;
pub mod clustering;
/// End-to-end 3D conformer generation (ETKDG), the crate's headline entry point.
///
/// [`conformer::generate_3d_conformer`] takes a [`graph::Molecule`] and a seed and returns an
/// `n × 3` coordinate matrix in Å, using the RDKit ETKDG recipe: distance-bounds matrix →
/// Floyd–Warshall smoothing → random-distance sampling → metric-matrix embedding (3D, or 4D when
/// the molecule has chiral centres) → bounds-force-field minimisation → ETKDG-3D minimisation →
/// planarity and double-bond geometry checks.
pub mod conformer;
pub mod cpm;
pub mod dipole;
pub mod dispersion;
/// Distance geometry: bounds matrices, chirality, embedding and validation.
///
/// This is the geometric core of ETKDG. It follows the classic distance-geometry pipeline
/// described in Blainey, Todeschini & van Gunsteren, *Chem. Rev.* **2000**, 100, 2557, and
/// mirrors RDKit's `DistGeom` / `DistGeomHelpers` structure:
///
/// 1. [`distgeom::calculate_bounds_matrix`] builds a lower/upper bound matrix from the molecular
///    graph and the UFF-derived reference bond lengths.
/// 2. [`distgeom::smooth_bounds_matrix`] / [`distgeom::triangle_smooth`] close the bounds under
///    the triangle inequality (Floyd–Warshall).
/// 3. [`distgeom::embedding`] samples random distances inside those bounds, builds a metric
///    matrix and diagonalises it into 3D (or 4D) coordinates.
/// 4. [`distgeom::validation`] checks tetrahedral centres, chiral volumes and planarity.
pub mod distgeom;
pub mod dos;
/// Molecular dynamics: velocity-Verlet integrators, thermostats and NEB paths.
///
/// Energies are in kcal/mol, times in femtoseconds (fs), temperatures in kelvin (K), masses in
/// atomic mass units (amu) and forces in kcal/mol/Å.
pub mod dynamics;
pub mod edl;
pub mod eht;
pub mod esp;
/// ETKDG internals: dihedral angles, torsion energy terms and their gradients.
///
/// All angles are in **radians** in the range `[-π, π]`; torsion energies are in kcal/mol.
pub mod etkdg;
pub mod experimental_status;
#[doc(hidden)]
pub mod fixture_io;
pub mod forcefield;
pub mod gpu;
pub mod graph;
pub mod hf;
pub mod imd;
pub mod ir;
pub mod kmesh;
pub mod materials;
pub mod mesh;
pub mod ml;
pub mod nmr;
pub mod optimization;
pub mod periodic;
pub mod pm3;
pub mod population;
pub mod potentials;
pub mod reactivity;
pub mod render_bridge;
pub mod riemannian;
pub mod rings;
pub mod scf;
pub mod sdr;
pub mod smarts;
/// SMILES parsing into a molecular graph.
///
/// A recursive-descent parser covering the Daylight SMILES grammar subset used by the crate:
/// aliphatic and aromatic atoms, bracket atoms with isotopes and formal charges, ring-closure
/// numbers, branches, `@`/`@@` tetrahedral stereo, `/`/`\` double-bond stereo, and implicit
/// hydrogen addition.
pub mod smiles;
pub mod smirks;
pub mod solvation;
pub mod solvation_alpb;
pub mod spectroscopy;
pub mod stereo;
pub mod surface;
/// Topology analysis: metal coordination environments from 3D coordinates.
///
/// Given atomic numbers and Cartesian positions (Å), detects metal centres, their ligand
/// spheres, coordination numbers and a geometric guess of the coordination polyhedron
/// (linear, trigonal, tetrahedral, square-planar, trigonal-bipyramidal, octahedral).
pub mod topology;
pub mod transport;
pub mod xtb;

#[cfg(any(
    feature = "alpha-cga",
    feature = "alpha-gsm",
    feature = "alpha-sdr",
    feature = "alpha-edl",
    feature = "alpha-periodic-linear",
    feature = "alpha-kinetics",
    feature = "alpha-render-bridge",
    feature = "alpha-reaction-dynamics"
))]
pub mod alpha;

#[cfg(any(
    feature = "beta-kpm",
    feature = "beta-mbh",
    feature = "beta-cpm",
    feature = "beta-randnla",
    feature = "beta-riemannian"
))]
pub mod beta;

// ─── ALPHA modules (feature-gated) ──────────────────────────────────────────

#[cfg(feature = "alpha-dynamics-live")]
pub mod dynamics_live;

#[cfg(feature = "alpha-dft")]
pub mod dft;

#[cfg(feature = "alpha-mlff")]
pub mod mlff {
    //! Re-export MLFF types from the ml module.
    pub use crate::ml::inference::*;
    pub use crate::ml::mlff::*;
    pub use crate::ml::symmetry_functions::*;
}

// ─── Internal helpers (not part of the public API) ───────────────────────────

/// Shared input validation (geometry, finiteness, parameter ranges).
pub(crate) mod validation;

use serde::{Deserialize, Serialize};

// ─── Public API Types ────────────────────────────────────────────────────────

/// Result of a 3D conformer generation for a single molecule.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConformerResult {
    /// Input SMILES string.
    pub smiles: String,
    /// Number of atoms in the molecule.
    pub num_atoms: usize,
    /// Flat xyz coordinates: [x0, y0, z0, x1, y1, z1, ...].
    /// Empty on failure.
    pub coords: Vec<f64>,
    /// Atom elements (atomic numbers) in the same order as coords.
    pub elements: Vec<u8>,
    /// Bond list as (start_atom, end_atom, order_string).
    pub bonds: Vec<(usize, usize, String)>,
    /// Error message if generation failed.
    pub error: Option<String>,
    /// Generation time in milliseconds.
    pub time_ms: f64,
}

/// Configuration for conformer generation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConformerConfig {
    /// RNG seed (same seed = reproducible output).
    pub seed: u64,
    /// Number of threads for batch processing (0 = auto-detect).
    pub num_threads: usize,
}

/// One conformer entry in a ranked conformer-search ensemble.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConformerEnsembleMember {
    /// RNG seed used for this embedding attempt.
    pub seed: u64,
    /// Cluster identifier assigned after RMSD clustering.
    pub cluster_id: Option<usize>,
    /// Flat xyz coordinates: [x0, y0, z0, x1, y1, z1, ...].
    pub coords: Vec<f64>,
    /// UFF energy in kcal/mol.
    pub energy_kcal_mol: f64,
}

/// One RMSD-based cluster summary in a conformer-search ensemble.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConformerClusterSummary {
    /// Cluster identifier.
    pub cluster_id: usize,
    /// Seed of the representative (lowest-energy) conformer in this cluster.
    pub representative_seed: u64,
    /// Number of conformers assigned to this cluster.
    pub size: usize,
    /// Seeds of all conformers assigned to this cluster.
    pub member_seeds: Vec<u64>,
}

/// Result of a conformer-search workflow with UFF ranking and RMSD filtering.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConformerSearchResult {
    /// Number of successful embedded conformers before filtering.
    pub generated: usize,
    /// Number of conformers retained after RMSD duplicate filtering.
    pub unique: usize,
    /// Number of rotatable bonds detected in the molecule.
    pub rotatable_bonds: usize,
    /// Ranked conformers (lowest UFF energy first).
    pub conformers: Vec<ConformerEnsembleMember>,
    /// RMSD-based clusters with representative/member mapping.
    pub clusters: Vec<ConformerClusterSummary>,
    /// Non-fatal notes and warnings.
    pub notes: Vec<String>,
}

/// Capability status for one operation on a given element set.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MethodCapability {
    /// Whether the operation is available for the provided element set.
    pub available: bool,
    /// Confidence level for the operation.
    pub confidence: eht::SupportLevel,
    /// Elements that block operation support.
    pub unsupported_elements: Vec<u8>,
    /// Human-readable warnings.
    pub warnings: Vec<String>,
}

/// Explicit computational method exposed by top-level planning APIs.
///
/// The variant names map to the concrete implementations available in the crate; see
/// [`system_method_plan`] for how a method is chosen for a given element set.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
// The variant set may grow in a future release; downstream crates must
// match with a wildcard arm.
#[non_exhaustive]
pub enum ScientificMethod {
    /// ETKDG distance-geometry embedding ([`conformer`]) — geometry only, no energy.
    Embed,
    /// UFF classical force field ([`forcefield::uff`]) — energy and gradients in kcal/mol.
    Uff,
    /// Extended Hückel theory ([`eht`]) — orbital energies in eV, qualitative electronic structure.
    Eht,
    /// PM3 semi-empirical method ([`pm3`]) — total energy in Hartree, orbital energies in eV.
    Pm3,
    /// xTB tight binding ([`xtb`]) — GFN2-xTB by default; total energy in Hartree.
    Xtb,
    /// MMFF94 classical force field ([`forcefield::mmff94`]) — energy in kcal/mol.
    Mmff94,
    /// ANI neural-network potential ([`ani`]) — energy in kcal/mol with bundled test weights.
    Ani,
    /// HF-3c minimal-basis Hartree–Fock with dispersion and gCP corrections ([`hf`]).
    Hf3c,
}

/// Property domain used when choosing a recommended computational method.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
// The variant set may grow in a future release; downstream crates must
// match with a wildcard arm.
#[non_exhaustive]
pub enum PropertyRequest {
    /// 3D molecular geometry generation.
    Geometry,
    /// Force-field total energy (UFF or MMFF94).
    ForceFieldEnergy,
    /// Molecular orbital energies and the HOMO–LUMO gap.
    Orbitals,
    /// Atomic partial charges (Mulliken / Löwdin).
    Population,
    /// 3D molecular-orbital isosurface grid for visualisation.
    OrbitalGrid,
}

/// Structured metadata for one method on a specific element set.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MethodMetadata {
    /// The method this row describes.
    pub method: ScientificMethod,
    /// Whether the method can run at all for the inspected element set.
    pub available: bool,
    /// Qualitative confidence in the method for this element set.
    pub confidence: eht::SupportLevel,
    /// Numeric confidence in the range `0.0` (unusable) to `1.0` (fully trusted).
    ///
    /// The score is a fixed lookup table, not a calibrated probability: UFF/MMFF94 score 0.95
    /// when supported and 0.75 when experimental, EHT/PM3/xTB score 0.95/0.60, ANI 0.90/0.70
    /// and HF-3c 0.85/0.65.
    pub confidence_score: f64,
    /// Human-readable statements of what the method cannot do here (e.g. unsupported elements).
    pub limitations: Vec<String>,
    /// Human-readable warnings surfaced by the underlying support analysis.
    pub warnings: Vec<String>,
}

/// Recommended method plan for one requested property domain.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PropertyMethodPlan {
    /// The property domain this plan addresses.
    pub property: PropertyRequest,
    /// The recommended method, or `None` when no method supports the element set.
    pub recommended: Option<ScientificMethod>,
    /// The fallback method to use when the recommended one is unavailable, if any.
    pub fallback: Option<ScientificMethod>,
    /// Free-text reasons behind the recommendation.
    pub rationale: Vec<String>,
    /// Per-method metadata rows considered while building the plan.
    pub methods: Vec<MethodMetadata>,
}

/// Full method-planning summary for a molecule element set.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SystemMethodPlan {
    /// Capability summary for the core operations on this element set.
    pub capabilities: SystemCapabilities,
    /// Plan for the [`PropertyRequest::Geometry`] domain.
    pub geometry: PropertyMethodPlan,
    /// Plan for the [`PropertyRequest::ForceFieldEnergy`] domain.
    pub force_field_energy: PropertyMethodPlan,
    /// Plan for the [`PropertyRequest::Orbitals`] domain.
    pub orbitals: PropertyMethodPlan,
    /// Plan for the [`PropertyRequest::Population`] domain.
    pub population: PropertyMethodPlan,
    /// Plan for the [`PropertyRequest::OrbitalGrid`] domain.
    pub orbital_grid: PropertyMethodPlan,
}

/// Capability summary for core operations on a given element set.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SystemCapabilities {
    /// Capability of the ETKDG embedding path (element presence only, always "experimental").
    pub embed: MethodCapability,
    /// Capability of the UFF force field (element coverage of the UFF parameter tables).
    pub uff: MethodCapability,
    /// Capability of the EHT electronic-structure path.
    pub eht: MethodCapability,
    /// Capability of the population analysis path (tracks EHT).
    pub population: MethodCapability,
    /// Capability of the orbital-grid rendering path (tracks EHT).
    pub orbital_grid: MethodCapability,
}

/// Result of an electronic workflow with optional UFF fallback.
///
/// The serde representation is externally tagged with the `mode` key, so the JSON form is
/// `{"mode": "eht", "result": {...}}` or `{"mode": "uff_fallback", ...}`.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "mode", rename_all = "snake_case")]
// The variant set may grow in a future release; downstream crates must
// match with a wildcard arm.
#[non_exhaustive]
pub enum ElectronicWorkflowResult {
    /// EHT ran successfully and produced a full electronic-structure result.
    Eht {
        /// The EHT solution: orbital energies in eV, MO coefficients, HOMO/LUMO levels and gap in eV.
        result: eht::EhtResult,
    },
    /// EHT was not applicable and the workflow fell back to a UFF force-field energy.
    UffFallback {
        /// UFF total energy in kcal/mol.
        energy_kcal_mol: f64,
        /// Human-readable explanation of why the fallback was taken.
        reason: String,
        /// The EHT support metadata that triggered the fallback.
        support: eht::EhtSupport,
    },
}

/// Execution status for one method inside a multi-method comparison run.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
// The variant set may grow in a future release; downstream crates must
// match with a wildcard arm.
#[non_exhaustive]
pub enum MethodComparisonStatus {
    /// The method ran and produced a payload.
    Success,
    /// The method was not run because it is unavailable or not permitted for this system.
    Unavailable,
    /// The method was run but failed; the message is in [`MethodComparisonEntry::error`].
    Error,
}

/// Compact per-method payload for multi-method comparison.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
// The variant set may grow in a future release; downstream crates must
// match with a wildcard arm.
#[non_exhaustive]
pub enum MethodComparisonPayload {
    /// Frontier-orbital observables from the EHT method.
    Eht {
        /// HOMO orbital energy in eV.
        homo_energy: f64,
        /// LUMO orbital energy in eV.
        lumo_energy: f64,
        /// HOMO–LUMO gap in eV.
        gap: f64,
        /// EHT capability metadata for the inspected element set.
        support: eht::EhtSupport,
    },
    /// Force-field energy observable.
    Uff {
        /// UFF total energy in kcal/mol.
        energy_kcal_mol: f64,
    },
}

/// One method row in the multi-method comparison workflow.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MethodComparisonEntry {
    /// The method this row describes.
    pub method: ScientificMethod,
    /// Whether the method ran successfully.
    pub status: MethodComparisonStatus,
    /// Whether the method is usable at all for this element set.
    pub available: bool,
    /// Qualitative confidence in the method for this element set.
    pub confidence: eht::SupportLevel,
    /// Numeric confidence in the range `0.0`–`1.0` (see [`MethodMetadata::confidence_score`]).
    pub confidence_score: f64,
    /// Human-readable warnings surfaced by the support analysis.
    pub warnings: Vec<String>,
    /// Human-readable statements of what the method cannot do here.
    pub limitations: Vec<String>,
    /// The method's numeric output, present only when [`Self::status`] is
    /// [`MethodComparisonStatus::Success`].
    pub payload: Option<MethodComparisonPayload>,
    /// Error message when the method failed or was skipped; `None` on success.
    pub error: Option<String>,
}

/// Structured comparison result for multiple methods on the same geometry/system.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MethodComparisonResult {
    /// The method plan that was applied to this system.
    pub plan: SystemMethodPlan,
    /// One row per method that was considered (UFF and EHT).
    pub comparisons: Vec<MethodComparisonEntry>,
}

/// Aromaticity analysis summary from graph-level bond annotations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AromaticityAnalysis {
    /// Per-atom aromatic flag; `true` when the atom belongs to at least one aromatic bond.
    pub aromatic_atoms: Vec<bool>,
    /// Aromatic bonds as `(atom_i, atom_j)` index pairs.
    pub aromatic_bonds: Vec<(usize, usize)>,
}

/// Graph-level stereocenter analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StereocenterAnalysis {
    /// Atom indices whose stereochemistry was explicitly tagged in the SMILES (`@`/`@@`).
    pub tagged_tetrahedral_centers: Vec<usize>,
    /// Atom indices inferred to be tetrahedral stereocentres: four neighbours, sp³ hybridisation
    /// and at least three distinct neighbour elements.
    pub inferred_tetrahedral_centers: Vec<usize>,
}

/// Combined structural graph feature analysis for UI/API consumers.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GraphFeatureAnalysis {
    /// Aromaticity results.
    pub aromaticity: AromaticityAnalysis,
    /// Stereocentre results.
    pub stereocenters: StereocenterAnalysis,
}

/// Configuration for the high-level GIAO NMR route.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GiaoNmrConfig {
    /// Total molecular charge.
    pub charge: i32,
    /// Spin multiplicity. The current public SCF route supports singlet closed-shell systems only.
    pub multiplicity: u32,
    /// Maximum number of SCF iterations.
    pub max_scf_iterations: usize,
    /// Enable rayon-parallel ERI evaluation when the feature is available.
    pub use_parallel_eri: bool,
    /// Allow the internal hydrogen-like fallback basis for unsupported elements.
    pub allow_basis_fallback: bool,
}

/// One isotope-specific shielding entry returned by the high-level GIAO API.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GiaoNmrEntry {
    /// Index of the nucleus in the input element/position arrays.
    pub atom_index: usize,
    /// Atomic number of the nucleus.
    pub element: u8,
    /// 3×3 chemical-shift tensor in ppm (rows = x, y, z components).
    pub tensor: [[f64; 3]; 3],
    /// Isotropic shielding σ_iso = Tr(σ)/3, in ppm.
    pub isotropic: f64,
    /// Anisotropic shielding in ppm.
    pub anisotropy: f64,
    /// Chemical shift δ = σ_ref − σ, in ppm.
    pub chemical_shift: f64,
}

/// Result of a high-level SCF-backed GIAO NMR calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GiaoNmrResult {
    /// Canonical nucleus label, e.g. `"1H"` or `"13C"`.
    pub nucleus: String,
    /// Atomic number of the requested nucleus.
    pub target_atomic_number: u8,
    /// Electronic-structure method label, always `"RHF/GIAO"` for this route.
    pub method: String,
    /// Basis set label, always `"STO-3G"` for this route.
    pub basis_set: String,
    /// Total molecular charge used in the SCF calculation.
    pub charge: i32,
    /// Spin multiplicity used in the SCF calculation (1 = singlet).
    pub multiplicity: u32,
    /// Whether the SCF reached its convergence threshold.
    pub scf_converged: bool,
    /// Number of SCF iterations performed.
    pub scf_iterations: usize,
    /// Converged SCF total energy in Hartree.
    pub total_energy_hartree: f64,
    /// Number of basis functions in the STO-3G expansion.
    pub n_basis: usize,
    /// Number of nuclei of the requested element found in the system.
    pub n_target_atoms: usize,
    /// Chemical shifts in ppm, in the same order as [`Self::shieldings`].
    pub chemical_shifts: Vec<f64>,
    /// Per-nucleus shielding entries.
    pub shieldings: Vec<GiaoNmrEntry>,
    /// Atomic numbers that fell back to the hydrogen-like minimal basis instead of explicit
    /// STO-3G data.
    pub fallback_elements: Vec<u8>,
    /// Human-readable notes about convergence, basis fallbacks and nucleus caveats.
    pub notes: Vec<String>,
}

impl Default for ConformerConfig {
    fn default() -> Self {
        Self {
            seed: 42,
            num_threads: 0,
        }
    }
}

impl Default for GiaoNmrConfig {
    fn default() -> Self {
        Self {
            charge: 0,
            multiplicity: 1,
            max_scf_iterations: 100,
            use_parallel_eri: false,
            allow_basis_fallback: false,
        }
    }
}

// ─── Public API Functions ────────────────────────────────────────────────────

/// Library version string.
pub fn version() -> String {
    format!("sci-form {}", env!("CARGO_PKG_VERSION"))
}

/// Report EHT capability metadata for a given element list.
pub fn eht_support(elements: &[u8]) -> eht::EhtSupport {
    eht::analyze_eht_support(elements)
}

fn is_uff_element_supported(z: u8) -> bool {
    matches!(
        z,
        1 | 5
            | 6
            | 7
            | 8
            | 9
            | 14
            | 15
            | 16
            | 17
            | 22
            | 23
            | 24
            | 25
            | 26
            | 27
            | 28
            | 29
            | 30
            | 32
            | 33
            | 34
            | 35
            | 42
            | 46
            | 47
            | 50
            | 51
            | 52
            | 53
            | 78
            | 79
    )
}

fn unique_sorted_unsupported(elements: &[u8], pred: impl Fn(u8) -> bool) -> Vec<u8> {
    let mut out: Vec<u8> = elements.iter().copied().filter(|&z| !pred(z)).collect();
    out.sort_unstable();
    out.dedup();
    out
}

/// Report operation capability metadata for a given element list.
pub fn system_capabilities(elements: &[u8]) -> SystemCapabilities {
    let eht_support = eht_support(elements);
    let uff_unsupported = unique_sorted_unsupported(elements, is_uff_element_supported);

    let embed = MethodCapability {
        available: !elements.is_empty(),
        confidence: eht::SupportLevel::Experimental,
        unsupported_elements: Vec::new(),
        warnings: vec![
            "Embed capability is inferred from element presence only; final success still depends on full molecular graph and geometry constraints.".to_string(),
        ],
    };

    let uff = if uff_unsupported.is_empty() {
        MethodCapability {
            available: true,
            confidence: eht::SupportLevel::Supported,
            unsupported_elements: Vec::new(),
            warnings: Vec::new(),
        }
    } else {
        MethodCapability {
            available: false,
            confidence: eht::SupportLevel::Unsupported,
            unsupported_elements: uff_unsupported.clone(),
            warnings: vec![format!(
                "UFF atom typing is unavailable for elements {:?}.",
                uff_unsupported
            )],
        }
    };

    let eht = MethodCapability {
        available: eht_support.unsupported_elements.is_empty(),
        confidence: eht_support.level,
        unsupported_elements: eht_support.unsupported_elements.clone(),
        warnings: eht_support.warnings.clone(),
    };

    let population = MethodCapability {
        available: eht.available,
        confidence: eht.confidence,
        unsupported_elements: eht.unsupported_elements.clone(),
        warnings: eht.warnings.clone(),
    };

    let orbital_grid = MethodCapability {
        available: eht.available,
        confidence: eht.confidence,
        unsupported_elements: eht.unsupported_elements.clone(),
        warnings: eht.warnings.clone(),
    };

    SystemCapabilities {
        embed,
        uff,
        eht,
        population,
        orbital_grid,
    }
}

fn confidence_score_for_method(method: ScientificMethod, capability: &MethodCapability) -> f64 {
    if !capability.available {
        return 0.0;
    }

    match method {
        ScientificMethod::Embed => 0.8,
        ScientificMethod::Uff | ScientificMethod::Mmff94 => match capability.confidence {
            eht::SupportLevel::Supported => 0.95,
            eht::SupportLevel::Experimental => 0.75,
            eht::SupportLevel::Unsupported => 0.0,
        },
        ScientificMethod::Eht | ScientificMethod::Pm3 | ScientificMethod::Xtb => {
            match capability.confidence {
                eht::SupportLevel::Supported => 0.95,
                eht::SupportLevel::Experimental => 0.6,
                eht::SupportLevel::Unsupported => 0.0,
            }
        }
        ScientificMethod::Ani => match capability.confidence {
            eht::SupportLevel::Supported => 0.90,
            eht::SupportLevel::Experimental => 0.7,
            eht::SupportLevel::Unsupported => 0.0,
        },
        ScientificMethod::Hf3c => match capability.confidence {
            eht::SupportLevel::Supported => 0.85,
            eht::SupportLevel::Experimental => 0.65,
            eht::SupportLevel::Unsupported => 0.0,
        },
    }
}

fn build_method_metadata(
    method: ScientificMethod,
    capability: &MethodCapability,
    extra_limitations: &[&str],
) -> MethodMetadata {
    let mut limitations: Vec<String> = extra_limitations.iter().map(|s| s.to_string()).collect();

    if !capability.unsupported_elements.is_empty() {
        limitations.push(format!(
            "Unsupported elements for this method: {:?}.",
            capability.unsupported_elements
        ));
    }

    if matches!(method, ScientificMethod::Eht)
        && matches!(capability.confidence, eht::SupportLevel::Experimental)
    {
        limitations.push(
            "Transition-metal EHT parameters remain provisional and should be treated as experimental."
                .to_string(),
        );
    }

    MethodMetadata {
        method,
        available: capability.available,
        confidence: capability.confidence,
        confidence_score: confidence_score_for_method(method, capability),
        limitations,
        warnings: capability.warnings.clone(),
    }
}

fn build_property_plan(
    property: PropertyRequest,
    recommended: Option<ScientificMethod>,
    fallback: Option<ScientificMethod>,
    rationale: Vec<String>,
    methods: Vec<MethodMetadata>,
) -> PropertyMethodPlan {
    PropertyMethodPlan {
        property,
        recommended,
        fallback,
        rationale,
        methods,
    }
}

/// Build a structured method plan with recommendations, fallback paths, and confidence scores.
pub fn system_method_plan(elements: &[u8]) -> SystemMethodPlan {
    let capabilities = system_capabilities(elements);

    let geometry_method = build_method_metadata(
        ScientificMethod::Embed,
        &capabilities.embed,
        &["Geometry generation still depends on graph topology, stereochemistry, and embedding constraints."],
    );
    let geometry = build_property_plan(
        PropertyRequest::Geometry,
        geometry_method.available.then_some(ScientificMethod::Embed),
        None,
        vec!["Embedding is the top-level geometry generation path in sci-form.".to_string()],
        vec![geometry_method],
    );

    let uff_method = build_method_metadata(
        ScientificMethod::Uff,
        &capabilities.uff,
        &["This recommendation applies to force-field energy evaluation, not molecular orbital analysis."],
    );
    let force_field_energy = build_property_plan(
        PropertyRequest::ForceFieldEnergy,
        uff_method.available.then_some(ScientificMethod::Uff),
        None,
        vec![
            "UFF is the top-level force-field energy path when atom typing is available."
                .to_string(),
        ],
        vec![uff_method],
    );

    let eht_method = build_method_metadata(
        ScientificMethod::Eht,
        &capabilities.eht,
        &["EHT is the only current top-level orbital method in sci-form."],
    );
    let orbitals = build_property_plan(
        PropertyRequest::Orbitals,
        eht_method.available.then_some(ScientificMethod::Eht),
        None,
        vec!["Orbital energies and MO coefficients are produced by the EHT workflow.".to_string()],
        vec![eht_method.clone()],
    );

    let population_method = build_method_metadata(
        ScientificMethod::Eht,
        &capabilities.population,
        &["Population analysis is derived from the EHT density and overlap matrices."],
    );
    let population = build_property_plan(
        PropertyRequest::Population,
        population_method.available.then_some(ScientificMethod::Eht),
        None,
        vec!["Population analysis currently requires a successful EHT calculation.".to_string()],
        vec![population_method],
    );

    let orbital_grid_method = build_method_metadata(
        ScientificMethod::Eht,
        &capabilities.orbital_grid,
        &["Orbital-grid rendering currently depends on EHT molecular-orbital coefficients."],
    );
    let orbital_grid = build_property_plan(
        PropertyRequest::OrbitalGrid,
        orbital_grid_method
            .available
            .then_some(ScientificMethod::Eht),
        None,
        vec![
            "Orbital-grid generation currently requires a successful EHT calculation.".to_string(),
        ],
        vec![orbital_grid_method],
    );

    SystemMethodPlan {
        capabilities,
        geometry,
        force_field_energy,
        orbitals,
        population,
        orbital_grid,
    }
}

/// Deprecated alias for [`eht_support`].
#[deprecated(
    since = "0.16.0",
    note = "renamed to `eht_support`; the `get_` prefix is not idiomatic in Rust"
)]
pub fn get_eht_support(elements: &[u8]) -> eht::EhtSupport {
    eht_support(elements)
}

/// Deprecated alias for [`system_capabilities`].
#[deprecated(
    since = "0.16.0",
    note = "renamed to `system_capabilities`; the `get_` prefix is not idiomatic in Rust"
)]
pub fn get_system_capabilities(elements: &[u8]) -> SystemCapabilities {
    system_capabilities(elements)
}

/// Deprecated alias for [`system_method_plan`].
#[deprecated(
    since = "0.16.0",
    note = "renamed to `system_method_plan`; the `get_` prefix is not idiomatic in Rust"
)]
pub fn get_system_method_plan(elements: &[u8]) -> SystemMethodPlan {
    system_method_plan(elements)
}

/// Generate a 3D conformer from a SMILES string.
pub fn embed(smiles: &str, seed: u64) -> ConformerResult {
    #[cfg(not(target_arch = "wasm32"))]
    let start = std::time::Instant::now();

    let mol = match graph::Molecule::from_smiles(smiles) {
        Ok(m) => m,
        Err(e) => {
            return ConformerResult {
                smiles: smiles.to_string(),
                num_atoms: 0,
                coords: vec![],
                elements: vec![],
                bonds: vec![],
                error: Some(e),
                #[cfg(not(target_arch = "wasm32"))]
                time_ms: start.elapsed().as_secs_f64() * 1000.0,
                #[cfg(target_arch = "wasm32")]
                time_ms: 0.0,
            };
        }
    };

    let n = mol.graph.node_count();
    let elements: Vec<u8> = (0..n)
        .map(|i| mol.graph[petgraph::graph::NodeIndex::new(i)].element)
        .collect();
    let bonds: Vec<(usize, usize, String)> = mol
        .graph
        .edge_indices()
        .map(|e| {
            let (a, b) = mol.graph.edge_endpoints(e).unwrap();
            let order = match mol.graph[e].order {
                graph::BondOrder::Single => "SINGLE",
                graph::BondOrder::Double => "DOUBLE",
                graph::BondOrder::Triple => "TRIPLE",
                graph::BondOrder::Aromatic => "AROMATIC",
                graph::BondOrder::Unknown => "UNKNOWN",
            };
            (a.index(), b.index(), order.to_string())
        })
        .collect();

    // Try multiple seeds for difficult molecules — if the first seed fails,
    // retry with different seeds before giving up.
    let max_seed_retries = 3u64;
    let mut last_err = String::new();
    for retry in 0..max_seed_retries {
        let current_seed = seed.wrapping_add(retry.wrapping_mul(997));
        match conformer::generate_3d_conformer(&mol, current_seed) {
            Ok(coords) => {
                let mut flat = Vec::with_capacity(n * 3);
                for i in 0..n {
                    flat.push(coords[(i, 0)] as f64);
                    flat.push(coords[(i, 1)] as f64);
                    flat.push(coords[(i, 2)] as f64);
                }
                return ConformerResult {
                    smiles: smiles.to_string(),
                    num_atoms: n,
                    coords: flat,
                    elements,
                    bonds,
                    error: None,
                    #[cfg(not(target_arch = "wasm32"))]
                    time_ms: start.elapsed().as_secs_f64() * 1000.0,
                    #[cfg(target_arch = "wasm32")]
                    time_ms: 0.0,
                };
            }
            Err(e) => {
                last_err = e;
            }
        }
    }
    ConformerResult {
        smiles: smiles.to_string(),
        num_atoms: n,
        coords: vec![],
        elements,
        bonds,
        error: Some(last_err),
        #[cfg(not(target_arch = "wasm32"))]
        time_ms: start.elapsed().as_secs_f64() * 1000.0,
        #[cfg(target_arch = "wasm32")]
        time_ms: 0.0,
    }
}

/// Batch-embed multiple SMILES in parallel.
///
/// Uses rayon thread pool for parallel processing.
/// `config.num_threads = 0` means auto-detect CPU count.
///
/// The requested thread count is clamped to the number of available logical
/// CPUs: asking for more threads than the machine has (e.g. `1 << 20`) used
/// to make the pool creation hang or exhaust the process address space.
/// A pool that still cannot be created is reported as an error instead of
/// panicking.
#[cfg(feature = "parallel")]
pub fn embed_batch(
    smiles_list: &[&str],
    config: &ConformerConfig,
) -> Result<Vec<ConformerResult>, String> {
    use rayon::prelude::*;

    if config.num_threads > 0 {
        let available = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1);
        let threads = config.num_threads.min(available);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|e| {
                format!(
                    "failed to create a rayon thread pool with {} threads: {}",
                    threads, e
                )
            })?;
        Ok(pool.install(|| {
            smiles_list
                .par_iter()
                .map(|smi| embed(smi, config.seed))
                .collect()
        }))
    } else {
        Ok(smiles_list
            .par_iter()
            .map(|smi| embed(smi, config.seed))
            .collect())
    }
}

/// Batch-embed multiple SMILES sequentially (no rayon dependency).
#[cfg(not(feature = "parallel"))]
pub fn embed_batch(
    smiles_list: &[&str],
    config: &ConformerConfig,
) -> Result<Vec<ConformerResult>, String> {
    Ok(smiles_list
        .iter()
        .map(|smi| embed(smi, config.seed))
        .collect())
}

/// Generate multiple conformers for a SMILES and filter by Butina RMSD clustering.
///
/// Generates `n_conformers` embeddings with different seeds, then clusters the
/// successful ones by RMSD and returns only the cluster centroids (diverse set).
///
/// # Arguments
/// - `smiles`: SMILES string
/// - `n_conformers`: number of embedding attempts (different seeds)
/// - `rmsd_cutoff`: RMSD threshold for clustering (Å), typically 1.0
/// - `base_seed`: base seed for reproducibility (seeds = base_seed..base_seed+n_conformers)
pub fn embed_diverse(
    smiles: &str,
    n_conformers: usize,
    rmsd_cutoff: f64,
    base_seed: u64,
) -> Vec<ConformerResult> {
    let all_results: Vec<ConformerResult> = (0..n_conformers as u64)
        .map(|i| embed(smiles, base_seed.wrapping_add(i)))
        .collect();

    let successful: Vec<(usize, &ConformerResult)> = all_results
        .iter()
        .enumerate()
        .filter(|(_, r)| r.error.is_none() && !r.coords.is_empty())
        .collect();

    if successful.len() <= 1 {
        return all_results
            .into_iter()
            .filter(|r| r.error.is_none())
            .collect();
    }

    let coords_vecs: Vec<Vec<f64>> = successful.iter().map(|(_, r)| r.coords.clone()).collect();
    let cluster_result = clustering::butina_cluster(&coords_vecs, rmsd_cutoff);

    cluster_result
        .centroid_indices
        .iter()
        .map(|&ci| {
            let (orig_idx, _) = successful[ci];
            all_results[orig_idx].clone()
        })
        .collect()
}

/// Parse a SMILES string and return molecular structure (no 3D generation).
///
/// # Errors
///
/// Returns `Err` with a human-readable message when the SMILES string is
/// syntactically invalid or uses a feature the parser does not support.
pub fn parse(smiles: &str) -> Result<graph::Molecule, String> {
    graph::Molecule::from_smiles(smiles)
}

/// Compute Gasteiger-Marsili partial charges from a SMILES string.
///
/// Parses the molecule, extracts bonds and elements, then runs 6 iterations
/// of electronegativity equalization.
///
/// The returned charges are in elementary charge (e) and sum to the total
/// formal charge of the molecule.
///
/// # Errors
///
/// Returns `Err` if the SMILES string cannot be parsed.
pub fn compute_charges(smiles: &str) -> Result<charges::gasteiger::ChargeResult, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    let n = mol.graph.node_count();
    let elements: Vec<u8> = (0..n)
        .map(|i| mol.graph[petgraph::graph::NodeIndex::new(i)].element)
        .collect();
    let formal_charges: Vec<i8> = (0..n)
        .map(|i| mol.graph[petgraph::graph::NodeIndex::new(i)].formal_charge)
        .collect();
    let bonds: Vec<(usize, usize, u8)> = mol
        .graph
        .edge_indices()
        .map(|e| {
            let (a, b) = mol.graph.edge_endpoints(e).unwrap();
            let order = match mol.graph[e].order {
                graph::BondOrder::Single => 1u8,
                graph::BondOrder::Double => 2,
                graph::BondOrder::Triple => 3,
                graph::BondOrder::Aromatic => 4,
                graph::BondOrder::Unknown => 1,
            };
            (a.index(), b.index(), order)
        })
        .collect();
    let modes: Vec<charges::gasteiger::GasteigerMode> = (0..n)
        .map(
            |i| match mol.graph[petgraph::graph::NodeIndex::new(i)].hybridization {
                graph::Hybridization::SP => charges::gasteiger::GasteigerMode::Sp,
                graph::Hybridization::SP2 => charges::gasteiger::GasteigerMode::Sp2,
                graph::Hybridization::SP3 => charges::gasteiger::GasteigerMode::Sp3,
                _ => charges::gasteiger::GasteigerMode::Default,
            },
        )
        .collect();
    charges::gasteiger::gasteiger_marsili_charges_full(
        &elements,
        &bonds,
        &formal_charges,
        Some(&modes),
        &charges::gasteiger::GasteigerConfig::default(),
    )
}

/// Compute solvent-accessible surface area from SMILES + 3D coordinates.
///
/// The `coords_flat` parameter is a flat [x0,y0,z0, x1,y1,z1,...] array in Å.
///
/// The result carries the total and per-atom SASA in Å² and the probe radius
/// that was used (Å). A `None` probe radius selects the module default of
/// 1.4 Å (the solvent-accessible convention).
///
/// # Errors
///
/// Returns `Err` when the element list and the coordinate array are
/// inconsistent (mismatched lengths, empty input, non-finite coordinates).
pub fn compute_sasa(
    elements: &[u8],
    coords_flat: &[f64],
    probe_radius: Option<f64>,
) -> Result<surface::sasa::SasaResult, String> {
    validation::validate_flat_geometry(elements, coords_flat)?;
    let positions: Vec<[f64; 3]> = coords_flat
        .chunks_exact(3)
        .map(|c| [c[0], c[1], c[2]])
        .collect();

    #[cfg(feature = "parallel")]
    {
        Ok(surface::sasa::compute_sasa_parallel(
            elements,
            &positions,
            probe_radius,
            None,
        ))
    }

    #[cfg(not(feature = "parallel"))]
    {
        Ok(surface::sasa::compute_sasa(
            elements,
            &positions,
            probe_radius,
            None,
        ))
    }
}

/// Compute Mulliken & Löwdin population analysis from atomic elements and positions.
///
/// `positions` is one `[x, y, z]` triple per atom, in Å.
///
/// The result also carries the EHT HOMO/LUMO levels and gap (all in eV), plus
/// coverage warnings when an element has no valence model or EHT basis
/// parameters. Partial charges are in elementary charge (e).
///
/// # Errors
///
/// Returns `Err` when the geometry is invalid or when the EHT solve fails
/// (for example because an element has no EHT valence-orbital parameters).
pub fn compute_population(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<population::PopulationResult, String> {
    validation::validate_geometry(elements, positions)?;
    let eht_result = eht::solve_eht(elements, positions, None)?;
    Ok(population::compute_population_with_orbitals(
        elements,
        positions,
        &eht_result.coefficients,
        &eht_result.energies,
        eht_result.n_electrons,
    ))
}

/// Compute molecular dipole moment (Debye) from atomic elements and positions.
pub fn compute_dipole(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<dipole::DipoleResult, String> {
    validation::validate_geometry(elements, positions)?;
    let eht_result = eht::solve_eht(elements, positions, None)?;
    Ok(dipole::compute_dipole_from_eht(
        elements,
        positions,
        &eht_result.coefficients,
        eht_result.n_electrons,
    ))
}

/// Compute atom-resolved HOMO/LUMO frontier descriptors from EHT.
///
/// Frontier orbital energies and the gap are in eV; per-atom frontier
/// coefficients are dimensionless weights that sum to one over the atoms
/// that carry the frontier orbital.
pub fn compute_frontier_descriptors(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<reactivity::FrontierDescriptors, String> {
    validation::validate_geometry(elements, positions)?;
    let eht_result = eht::solve_eht(elements, positions, None)?;
    Ok(reactivity::compute_frontier_descriptors(
        elements,
        positions,
        &eht_result,
    ))
}

/// Compute Fukui-function workflows and condensed per-atom descriptors from EHT.
///
/// The condensed per-atom Fukui indices are dimensionless; the frontier
/// orbital energies they are derived from are in eV.
pub fn compute_fukui_descriptors(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<reactivity::FukuiDescriptors, String> {
    validation::validate_geometry(elements, positions)?;
    let eht_result = eht::solve_eht(elements, positions, None)?;
    Ok(reactivity::compute_fukui_descriptors(
        elements,
        positions,
        &eht_result,
    ))
}

/// Build empirical reactivity rankings using condensed Fukui descriptors and Mulliken charges.
///
/// Ranking scores are dimensionless; the underlying charges are in
/// elementary charge (e) and the frontier energies in eV.
pub fn compute_reactivity_ranking(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<reactivity::ReactivityRanking, String> {
    validation::validate_geometry(elements, positions)?;
    let eht_result = eht::solve_eht(elements, positions, None)?;
    let fukui = reactivity::compute_fukui_descriptors(elements, positions, &eht_result);
    let pop = population::compute_population(
        elements,
        positions,
        &eht_result.coefficients,
        eht_result.n_electrons,
    );
    Ok(reactivity::rank_reactivity_sites(
        &fukui,
        &pop.mulliken_charges,
    ))
}

/// Build an exploratory UV-Vis-like spectrum from low-cost EHT transitions.
///
/// `sigma` is the Gaussian broadening width in eV, `e_min`/`e_max` bound
/// the scanned transition-energy window in eV and `n_points` is the number of
/// grid points. The returned energies are in eV and the intensities are
/// dimensionless relative oscillator strengths.
///
/// This is an exploratory approximation built on EHT level spacings, not a
/// TDDFT or CIS spectrum.
pub fn compute_uv_vis_spectrum(
    elements: &[u8],
    positions: &[[f64; 3]],
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> Result<reactivity::UvVisSpectrum, String> {
    validation::validate_geometry(elements, positions)?;
    let eht_result = eht::solve_eht(elements, positions, None)?;
    Ok(reactivity::compute_uv_vis_like_spectrum(
        &eht_result,
        sigma,
        e_min,
        e_max,
        n_points,
    ))
}

/// Compute Wiberg-like and Mayer-like bond orders from EHT.
///
/// The returned bond orders are dimensionless indices in roughly the 0-3
/// range (1 = single, 2 = double, 3 = triple), derived from EHT density
/// matrix elements rather than from a correlated wave function.
pub fn compute_bond_orders(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<population::BondOrderResult, String> {
    validation::validate_geometry(elements, positions)?;
    let eht_result = eht::solve_eht(elements, positions, None)?;
    Ok(population::compute_bond_orders(
        elements,
        positions,
        &eht_result.coefficients,
        eht_result.n_electrons,
    ))
}

/// Compute structured topology analysis for transition-metal coordination environments.
pub fn compute_topology(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> topology::TopologyAnalysisResult {
    topology::analyze_topology(elements, positions)
}

/// Analyze aromaticity and graph-level stereocenters from a SMILES string.
///
/// # Errors
///
/// Returns `Err` when the SMILES string cannot be parsed.
pub fn analyze_graph_features(smiles: &str) -> Result<GraphFeatureAnalysis, String> {
    use petgraph::visit::EdgeRef;

    let mol = parse(smiles)?;
    let n_atoms = mol.graph.node_count();
    let mut aromatic_atoms = vec![false; n_atoms];
    let mut aromatic_bonds = Vec::new();

    for edge in mol.graph.edge_references() {
        if matches!(edge.weight().order, graph::BondOrder::Aromatic) {
            let i = edge.source().index();
            let j = edge.target().index();
            aromatic_atoms[i] = true;
            aromatic_atoms[j] = true;
            aromatic_bonds.push((i, j));
        }
    }

    let mut tagged_tetrahedral_centers = Vec::new();
    let mut inferred_tetrahedral_centers = Vec::new();
    for i in 0..n_atoms {
        let idx = petgraph::graph::NodeIndex::new(i);
        let atom = &mol.graph[idx];
        if matches!(
            atom.chiral_tag,
            graph::ChiralType::TetrahedralCW | graph::ChiralType::TetrahedralCCW
        ) {
            tagged_tetrahedral_centers.push(i);
        }

        let neighbors: Vec<_> = mol.graph.neighbors(idx).collect();
        if neighbors.len() == 4 && matches!(atom.hybridization, graph::Hybridization::SP3) {
            let mut signature: Vec<u8> = neighbors.iter().map(|n| mol.graph[*n].element).collect();
            signature.sort_unstable();
            signature.dedup();
            if signature.len() >= 3 {
                inferred_tetrahedral_centers.push(i);
            }
        }
    }

    Ok(GraphFeatureAnalysis {
        aromaticity: AromaticityAnalysis {
            aromatic_atoms,
            aromatic_bonds,
        },
        stereocenters: StereocenterAnalysis {
            tagged_tetrahedral_centers,
            inferred_tetrahedral_centers,
        },
    })
}

/// Compute electronic properties with automatic fallback to UFF energy.
///
/// If EHT is unsupported for the element set, this function routes directly to UFF.
/// If EHT is experimental and `allow_experimental_eht` is false, it also routes to UFF.
pub fn compute_eht_or_uff_fallback(
    smiles: &str,
    elements: &[u8],
    positions: &[[f64; 3]],
    allow_experimental_eht: bool,
) -> Result<ElectronicWorkflowResult, String> {
    let support = eht_support(elements);
    let should_fallback = match support.level {
        eht::SupportLevel::Unsupported => true,
        eht::SupportLevel::Experimental => !allow_experimental_eht,
        eht::SupportLevel::Supported => false,
    };

    if should_fallback {
        let coords_flat: Vec<f64> = positions.iter().flat_map(|p| p.iter().copied()).collect();
        let energy = compute_uff_energy(smiles, &coords_flat).map_err(|e| {
            format!(
                "EHT is not appropriate for this system (support: {:?}) and UFF fallback failed: {}",
                support.level, e
            )
        })?;
        return Ok(ElectronicWorkflowResult::UffFallback {
            energy_kcal_mol: energy,
            reason: if matches!(support.level, eht::SupportLevel::Unsupported) {
                "EHT unsupported for one or more elements; routed to UFF-only workflow.".to_string()
            } else {
                "EHT confidence is experimental and experimental mode is disabled; routed to UFF-only workflow."
                    .to_string()
            },
            support,
        });
    }

    let result = eht::solve_eht(elements, positions, None)?;
    Ok(ElectronicWorkflowResult::Eht { result })
}

/// Compare multiple supported methods on the same geometry/system.
///
/// This workflow executes available methods independently and returns per-method
/// status, confidence metadata, warnings, limitations, and compact outputs.
pub fn compare_methods(
    smiles: &str,
    elements: &[u8],
    positions: &[[f64; 3]],
    allow_experimental_eht: bool,
) -> MethodComparisonResult {
    let plan = system_method_plan(elements);
    let mut comparisons = Vec::new();

    let coords_flat: Vec<f64> = positions.iter().flat_map(|p| p.iter().copied()).collect();

    {
        let meta = build_method_metadata(
            ScientificMethod::Uff,
            &plan.capabilities.uff,
            &["Comparison uses UFF force-field energy as the UFF observable."],
        );
        if !meta.available {
            comparisons.push(MethodComparisonEntry {
                method: ScientificMethod::Uff,
                status: MethodComparisonStatus::Unavailable,
                available: false,
                confidence: meta.confidence,
                confidence_score: meta.confidence_score,
                warnings: meta.warnings,
                limitations: meta.limitations,
                payload: None,
                error: Some("UFF is unavailable for this element set.".to_string()),
            });
        } else {
            match compute_uff_energy(smiles, &coords_flat) {
                Ok(energy) => comparisons.push(MethodComparisonEntry {
                    method: ScientificMethod::Uff,
                    status: MethodComparisonStatus::Success,
                    available: true,
                    confidence: meta.confidence,
                    confidence_score: meta.confidence_score,
                    warnings: meta.warnings,
                    limitations: meta.limitations,
                    payload: Some(MethodComparisonPayload::Uff {
                        energy_kcal_mol: energy,
                    }),
                    error: None,
                }),
                Err(err) => comparisons.push(MethodComparisonEntry {
                    method: ScientificMethod::Uff,
                    status: MethodComparisonStatus::Error,
                    available: true,
                    confidence: meta.confidence,
                    confidence_score: meta.confidence_score,
                    warnings: meta.warnings,
                    limitations: meta.limitations,
                    payload: None,
                    error: Some(err),
                }),
            }
        }
    }

    {
        let meta = build_method_metadata(
            ScientificMethod::Eht,
            &plan.capabilities.eht,
            &["Comparison uses frontier orbital energies and gap as the EHT observable."],
        );

        if !meta.available {
            comparisons.push(MethodComparisonEntry {
                method: ScientificMethod::Eht,
                status: MethodComparisonStatus::Unavailable,
                available: false,
                confidence: meta.confidence,
                confidence_score: meta.confidence_score,
                warnings: meta.warnings,
                limitations: meta.limitations,
                payload: None,
                error: Some("EHT is unavailable for this element set.".to_string()),
            });
        } else if matches!(meta.confidence, eht::SupportLevel::Experimental)
            && !allow_experimental_eht
        {
            comparisons.push(MethodComparisonEntry {
                method: ScientificMethod::Eht,
                status: MethodComparisonStatus::Unavailable,
                available: true,
                confidence: meta.confidence,
                confidence_score: meta.confidence_score,
                warnings: meta.warnings,
                limitations: meta.limitations,
                payload: None,
                error: Some(
                    "EHT confidence is experimental and allow_experimental_eht=false.".to_string(),
                ),
            });
        } else {
            match eht::solve_eht(elements, positions, None) {
                Ok(result) => comparisons.push(MethodComparisonEntry {
                    method: ScientificMethod::Eht,
                    status: MethodComparisonStatus::Success,
                    available: true,
                    confidence: meta.confidence,
                    confidence_score: meta.confidence_score,
                    warnings: meta.warnings,
                    limitations: meta.limitations,
                    payload: Some(MethodComparisonPayload::Eht {
                        homo_energy: result.homo_energy,
                        lumo_energy: result.lumo_energy,
                        gap: result.gap,
                        support: result.support,
                    }),
                    error: None,
                }),
                Err(err) => comparisons.push(MethodComparisonEntry {
                    method: ScientificMethod::Eht,
                    status: MethodComparisonStatus::Error,
                    available: true,
                    confidence: meta.confidence,
                    confidence_score: meta.confidence_score,
                    warnings: meta.warnings,
                    limitations: meta.limitations,
                    payload: None,
                    error: Some(err),
                }),
            }
        }
    }

    MethodComparisonResult { plan, comparisons }
}

/// Compute ESP grid from atomic elements, positions and Mulliken charges.
///
/// Invalid `spacing`/`padding` values are rejected with an error; see
/// [`compute_esp_with_config`] for explicit evaluation settings.
pub fn compute_esp(
    elements: &[u8],
    positions: &[[f64; 3]],
    spacing: f64,
    padding: f64,
) -> Result<esp::EspGrid, String> {
    compute_esp_with_config(
        elements,
        positions,
        esp::EspGridConfig {
            spacing,
            padding,
            min_distance: esp::DEFAULT_ESP_MIN_DISTANCE,
        },
    )
}

/// Compute the ESP grid with explicit evaluation settings
/// (`spacing`, `padding` and the nuclear `min_distance` clamp).
pub fn compute_esp_with_config(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: esp::EspGridConfig,
) -> Result<esp::EspGrid, String> {
    config.validate()?;
    let pop = compute_population(elements, positions)?;
    let (grid, _report) = gpu::esp_grid_gpu::compute_esp_grid_with_report(
        elements,
        positions,
        &pop.mulliken_charges,
        config,
    );
    Ok(grid)
}

/// Compute DOS/PDOS from EHT orbital energies.
///
/// `sigma` is the Gaussian smearing width in eV, `e_min`/`e_max` bound the
/// energy window in eV and `n_points` is the grid resolution. The result
/// carries the energy grid in eV and the total DOS and per-atom PDOS in
/// states/eV.
///
/// # Errors
///
/// Returns `Err` when `n_points < 2`, when `sigma` is not a positive
/// finite number, or when the EHT solve fails.
pub fn compute_dos(
    elements: &[u8],
    positions: &[[f64; 3]],
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> Result<dos::DosResult, String> {
    let eht_result = eht::solve_eht(elements, positions, None)?;
    let flat: Vec<f64> = positions.iter().flat_map(|p| p.iter().copied()).collect();

    #[cfg(feature = "parallel")]
    {
        dos::compute_pdos_parallel(
            elements,
            &flat,
            &eht_result.energies,
            &eht_result.coefficients,
            eht_result.n_electrons,
            sigma,
            e_min,
            e_max,
            n_points,
        )
    }

    #[cfg(not(feature = "parallel"))]
    {
        dos::compute_pdos(
            elements,
            &flat,
            &eht_result.energies,
            &eht_result.coefficients,
            eht_result.n_electrons,
            sigma,
            e_min,
            e_max,
            n_points,
        )
    }
}

/// Compute RMSD between two coordinate sets after Kabsch alignment.
pub fn compute_rmsd(coords: &[f64], reference: &[f64]) -> f64 {
    alignment::compute_rmsd(coords, reference)
}

/// Compute UFF force field energy for a molecule.
///
/// `coords` is a flat `[x0, y0, z0, x1, y1, z1, ...]` array in Å with one
/// triplet per atom. The returned energy is the UFF total energy in kcal/mol.
///
/// # Errors
///
/// Returns `Err` when the SMILES string cannot be parsed, or when `coords`
/// does not have exactly `3 × n_atoms` entries.
pub fn compute_uff_energy(smiles: &str, coords: &[f64]) -> Result<f64, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    let n = mol.graph.node_count();
    if coords.len() != n * 3 {
        return Err(format!("coords length {} != 3 * atoms {}", coords.len(), n));
    }
    let ff = forcefield::builder::build_uff_force_field(&mol);
    let mut gradient = vec![0.0f64; n * 3];
    let energy = ff.compute_system_energy_and_gradients(coords, &mut gradient);
    Ok(energy)
}

/// Compute MMFF94 force field energy for a molecule.
///
/// `smiles`: SMILES string for bond/topology information.
/// `coords`: flat xyz coordinates `[x0,y0,z0, x1,y1,z1, ...]`.
///
/// The energy is the MMFF94 valence + dispersion part of the total MMFF94
/// energy: bond stretching, angle bending, torsion and van der Waals
/// (buffered 14-7) terms, evaluated with the published MMFF94 parameter
/// tables.  MMFF94's electrostatic, stretch-bend and out-of-plane bending
/// terms are not implemented, so this is *not* the complete MMFF94 energy.
///
/// Returns total MMFF94 energy in kcal/mol.
///
/// # Errors
///
/// Returns `Err` when the SMILES string cannot be parsed, or when `coords`
/// does not have exactly `3 × n_atoms` entries.
pub fn compute_mmff94_energy(smiles: &str, coords: &[f64]) -> Result<f64, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    let n = mol.graph.node_count();
    if coords.len() != n * 3 {
        return Err(format!("coords length {} != 3 * atoms {}", coords.len(), n));
    }
    let terms = forcefield::mmff94::Mmff94Builder::build_typed(
        &mol,
        forcefield::mmff94::Mmff94Variant::Mmff94,
    );
    let (energy, _grad) = forcefield::mmff94::Mmff94Builder::total_energy(&terms, coords);
    Ok(energy)
}

/// Per-term-type UFF energy breakdown (kcal/mol) for a molecule.
///
/// Types are `bond`, `angle`, `torsion`, `inversion` and `vdw` (plus any
/// diagnostics term), in first-seen order; the values sum to
/// [`compute_uff_energy`].
pub fn compute_uff_energy_breakdown(
    smiles: &str,
    coords: &[f64],
) -> Result<Vec<(String, f64)>, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    let n = mol.graph.node_count();
    if coords.len() != n * 3 {
        return Err(format!("coords length {} != 3 * atoms {}", coords.len(), n));
    }
    let ff = forcefield::builder::build_uff_force_field(&mol);
    Ok(ff
        .energy_breakdown_aggregated(coords)
        .into_iter()
        .map(|(kind, energy)| (kind.to_string(), energy))
        .collect())
}

/// Per-term-type MMFF94 energy breakdown (kcal/mol) for a molecule.
///
/// The values sum to [`compute_mmff94_energy`].
pub fn compute_mmff94_energy_breakdown(
    smiles: &str,
    coords: &[f64],
) -> Result<Vec<(String, f64)>, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    let n = mol.graph.node_count();
    if coords.len() != n * 3 {
        return Err(format!("coords length {} != 3 * atoms {}", coords.len(), n));
    }
    let terms = forcefield::mmff94::Mmff94Builder::build_typed(
        &mol,
        forcefield::mmff94::Mmff94Variant::Mmff94,
    );
    Ok(forcefield::traits::aggregate_breakdown(&terms, coords)
        .into_iter()
        .map(|(kind, energy)| (kind.to_string(), energy))
        .collect())
}

/// Run a PM3 semi-empirical calculation.
///
/// `elements`: atomic numbers.
/// `positions`: Cartesian coordinates in Å, one `[x,y,z]` per atom.
///
/// Returns orbital energies, total energy, heat of formation, and Mulliken charges.
/// Compute PM3 semi-empirical energy and orbital information.
///
/// Main-group elements use Stewart's PM3 parameters and transition metals the
/// PM3(tm) set ([`pm3::pm3_method`] reports which one applies to an element).
/// MOPAC parity is a documented gap (validation F42).
///
/// # Units
///
/// The total energy is in Hartree, the orbital energies and the heat of
/// formation in eV, and the Mulliken charges in elementary charge (e).
///
/// # Errors
///
/// Returns `Err` when the geometry fails validation (mismatched element and
/// position counts, empty input, non-finite coordinates), when an element has
/// no PM3 parameter set, or when the SCF cycle does not converge.
pub fn compute_pm3(elements: &[u8], positions: &[[f64; 3]]) -> Result<pm3::Pm3Result, String> {
    validation::validate_geometry(elements, positions)?;
    pm3::solve_pm3(elements, positions)
}

/// Compute PM3 analytical energy gradients.
///
/// Returns dE/dR in eV/Å for each atom, computed from the converged
/// SCF density matrix using Hellmann-Feynman + Pulay corrections.
pub fn compute_pm3_gradient(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<pm3::Pm3GradientResult, String> {
    validation::validate_geometry(elements, positions)?;
    pm3::compute_pm3_gradient(elements, positions)
}

/// Run an xTB tight-binding calculation.
///
/// `elements`: atomic numbers.
/// `positions`: Cartesian coordinates in Å, one `[x,y,z]` per atom.
///
/// Returns orbital energies, total energy, gap, and Mulliken charges.
/// Uses GFN2-xTB internally for reliable SCC convergence and accurate energetics.
/// Compute a GFN1-xTB single point (electrical structure + total energy).
///
/// This is the root-level entry point for the GFN1-xTB method; the
/// implementation lives in [`xtb::gfn1`].
///
/// # Units
///
/// Orbital energies, HOMO/LUMO levels and the gap are in eV; the electronic,
/// repulsive and total energies are in Hartree; Mulliken charges are in
/// elementary charge (e).
///
/// # Errors
///
/// Returns `Err` when the geometry is invalid, when an element has no GFN1
/// parameter set, or when the SCC cycle does not converge.
pub fn compute_gfn1(elements: &[u8], positions: &[[f64; 3]]) -> Result<xtb::Gfn1Result, String> {
    xtb::gfn1::solve_gfn1(elements, positions)
}

/// Compute a GFN2-xTB single point (electrical structure + total energy).
///
/// This is the root-level entry point for the GFN2-xTB method; the
/// implementation lives in [`xtb::gfn2`].
///
/// # Units
///
/// Orbital energies, HOMO/LUMO levels and the gap are in eV; the electronic,
/// repulsive and total energies are in Hartree; Mulliken charges are in
/// elementary charge (e).
///
/// # Errors
///
/// Returns `Err` when the geometry is invalid, when an element has no GFN2
/// parameter set, or when the SCC cycle does not converge.
pub fn compute_gfn2(elements: &[u8], positions: &[[f64; 3]]) -> Result<xtb::Gfn2Result, String> {
    xtb::gfn2::solve_gfn2(elements, positions)
}

/// Compute a generic xTB single point (electrical structure + total energy).
///
/// This is a convenience wrapper that runs **GFN2-xTB** internally (the most reliable of the
/// GFN family for SCC convergence and energetics) and re-packages the result into the common
/// [`xtb::XtbResult`] type. Use [`compute_gfn1`] or [`compute_gfn2`] directly when the specific
/// GFN parameterisation matters.
///
/// # Units
///
/// Orbital energies, HOMO/LUMO levels and the gap are in eV; electronic, repulsive and total
/// energies are in Hartree; Mulliken charges are in elementary charge (e).
///
/// # Errors
///
/// Returns `Err` when the geometry fails validation (mismatched element/position counts, empty
/// input, non-finite coordinates) or when the GFN2-xTB SCC cycle fails to converge.
pub fn compute_xtb(elements: &[u8], positions: &[[f64; 3]]) -> Result<xtb::XtbResult, String> {
    validation::validate_geometry(elements, positions)?;
    let gfn2 = xtb::gfn2::solve_gfn2(elements, positions)?;
    Ok(xtb::XtbResult {
        orbital_energies: gfn2.orbital_energies,
        electronic_energy: gfn2.electronic_energy,
        repulsive_energy: gfn2.repulsive_energy,
        total_energy: gfn2.total_energy,
        n_basis: gfn2.n_basis,
        n_electrons: gfn2.n_electrons,
        homo_energy: gfn2.homo_energy,
        lumo_energy: gfn2.lumo_energy,
        gap: gfn2.gap,
        mulliken_charges: gfn2.mulliken_charges,
        scc_iterations: gfn2.scc_iterations,
        converged: gfn2.converged,
        notes: Vec::new(),
    })
}

/// Compute xTB analytical energy gradients.
///
/// Returns dE/dR in eV/Å for each atom, computed from the converged
/// SCC density using Hellmann-Feynman + Pulay + repulsive corrections.
pub fn compute_xtb_gradient(
    elements: &[u8],
    positions: &[[f64; 3]],
) -> Result<xtb::XtbGradientResult, String> {
    validation::validate_geometry(elements, positions)?;
    xtb::compute_xtb_gradient(elements, positions)
}

// ─── Spectroscopy API (Track D) ─────────────────────────────────────────────

/// Compute an sTDA UV-Vis spectrum with proper oscillator strengths.
///
/// Uses EHT/xTB molecular orbital transitions with transition dipole moment
/// evaluation and configurable broadening (Gaussian or Lorentzian).
///
/// `sigma` is the broadening width in eV, `e_min`/`e_max` bound the
/// scanned energy window in eV and `n_points` is the grid resolution.
/// The returned transition energies are in eV and the oscillator strengths
/// are dimensionless.
pub fn compute_stda_uvvis(
    elements: &[u8],
    positions: &[[f64; 3]],
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
    broadening: reactivity::BroadeningType,
) -> Result<reactivity::StdaUvVisSpectrum, String> {
    validation::validate_geometry(elements, positions)?;
    reactivity::compute_stda_uvvis_spectrum(
        elements, positions, sigma, e_min, e_max, n_points, broadening,
    )
}

/// Fluorescence spectrum from the same sTDA singlets as UV-Vis absorption.
///
/// Selects the emitting state by Kasha's rule, builds the spontaneous-emission
/// lineshape, and reports the Einstein A radiative rate. A Stokes shift is
/// applied only when `config.stokes_ev` is set; this path does not optimize S1.
pub fn compute_fluorescence(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: spectroscopy::FluorescenceConfig,
) -> Result<spectroscopy::FluorescenceSpectrum, String> {
    validation::validate_geometry(elements, positions)?;
    spectroscopy::compute_fluorescence_spectrum(elements, positions, &config)
}

/// Perform vibrational analysis via numerical Hessian.
///
/// Computes vibrational frequencies (cm⁻¹), normal modes, IR intensities,
/// and zero-point vibrational energy from a semiempirical Hessian.
///
/// `method`: "eht", "pm3", "xtb", or "uff"
pub fn compute_vibrational_analysis(
    elements: &[u8],
    positions: &[[f64; 3]],
    method: &str,
    step_size: Option<f64>,
) -> Result<ir::VibrationalAnalysis, String> {
    validation::validate_geometry(elements, positions)?;
    let hessian_method =
        match method {
            "eht" => ir::HessianMethod::Eht,
            "pm3" => ir::HessianMethod::Pm3,
            "xtb" => ir::HessianMethod::Xtb,
            "uff" => return Err(
                "UFF vibrational analysis requires SMILES; use compute_vibrational_analysis_uff"
                    .to_string(),
            ),
            _ => return Err(format!("Unknown method '{}', use eht/pm3/xtb/uff", method)),
        };
    ir::compute_vibrational_analysis(elements, positions, hessian_method, step_size)
}

/// Perform vibrational analysis using UFF analytical Hessian.
///
/// Uses gradient-difference method (O(6N) gradient evaluations) instead of
/// the standard O(9N²) energy evaluations, leveraging UFF's analytical gradients.
pub fn compute_vibrational_analysis_uff(
    smiles: &str,
    elements: &[u8],
    positions: &[[f64; 3]],
    step_size: Option<f64>,
) -> Result<ir::VibrationalAnalysis, String> {
    validation::validate_geometry(elements, positions)?;
    ir::vibrations::compute_vibrational_analysis_uff(smiles, elements, positions, step_size)
}

/// Generate a Lorentzian-broadened IR spectrum from vibrational analysis.
///
/// `gamma`: line width in cm⁻¹ (typically 10–30)
/// `wn_min`, `wn_max`: wavenumber range in cm⁻¹
/// `n_points`: grid resolution
pub fn compute_ir_spectrum(
    analysis: &ir::VibrationalAnalysis,
    gamma: f64,
    wn_min: f64,
    wn_max: f64,
    n_points: usize,
) -> ir::IrSpectrum {
    ir::compute_ir_spectrum(analysis, gamma, wn_min, wn_max, n_points)
}

/// Predict representative NMR chemical shifts from SMILES.
///
/// The shifts are in ppm and come from an empirical, topology-based
/// increment scheme (not from a GIAO calculation) - use
/// [`compute_giao_nmr`] for a quantum-chemical estimate.
///
/// # Errors
///
/// Returns `Err` if the SMILES string cannot be parsed.
pub fn predict_nmr_shifts(smiles: &str) -> Result<nmr::NmrShiftResult, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    Ok(nmr::predict_chemical_shifts(&mol))
}

/// Predict chemical shifts for a specific NMR nucleus from SMILES.
pub fn predict_nmr_shifts_for_nucleus(
    smiles: &str,
    nucleus: &str,
) -> Result<Vec<nmr::ChemicalShift>, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    let nucleus = parse_nmr_nucleus(nucleus)?;
    Ok(nmr::predict_chemical_shifts_for_nucleus(&mol, nucleus))
}

/// Predict J-coupling constants for a molecule.
///
/// If positions are provided, uses Karplus equation for 3D-dependent ³J values.
/// Otherwise uses topological estimates.
pub fn predict_nmr_couplings(
    smiles: &str,
    positions: &[[f64; 3]],
) -> Result<Vec<nmr::JCoupling>, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    if !positions.is_empty() {
        let elements: Vec<u8> = (0..mol.graph.node_count())
            .map(|i| mol.graph[petgraph::graph::NodeIndex::new(i)].element)
            .collect();
        validation::validate_geometry(&elements, positions)?;
    }
    Ok(nmr::predict_j_couplings(&mol, positions))
}

fn parse_nmr_nucleus(nucleus: &str) -> Result<nmr::NmrNucleus, String> {
    nmr::NmrNucleus::parse(nucleus)
}

struct GiaoNmrPreflight {
    nucleus: nmr::NmrNucleus,
    system: scf::types::MolecularSystem,
    basis: scf::basis::BasisSet,
    fallback_elements: Vec<u8>,
}

const EXPLICIT_GIAO_STO3G_ELEMENTS: &[u8] =
    &[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 35, 53];

fn has_explicit_giao_sto3g_support(z: u8) -> bool {
    EXPLICIT_GIAO_STO3G_ELEMENTS.contains(&z)
}

fn format_element_label(z: u8) -> String {
    nmr::NmrNucleus::default_for_element(z)
        .map(|nucleus| format!("{}({})", nucleus.element_symbol(), z))
        .unwrap_or_else(|| format!("Z{}", z))
}

fn format_element_set(elements: &[u8]) -> String {
    let mut unique = elements.to_vec();
    unique.sort_unstable();
    unique.dedup();
    unique
        .into_iter()
        .map(format_element_label)
        .collect::<Vec<_>>()
        .join(", ")
}

fn preflight_giao_nmr(
    elements: &[u8],
    positions: &[[f64; 3]],
    nucleus: &str,
    config: &GiaoNmrConfig,
) -> Result<GiaoNmrPreflight, String> {
    if elements.is_empty() {
        return Err("GIAO NMR requires at least one atom.".to_string());
    }
    validation::validate_geometry(elements, positions)?;
    if config.max_scf_iterations == 0 {
        return Err("GIAO NMR requires max_scf_iterations > 0.".to_string());
    }
    if config.multiplicity != 1 {
        return Err(format!(
            "The public GIAO NMR API currently supports singlet closed-shell systems only; got multiplicity {}.",
            config.multiplicity
        ));
    }

    let parsed_nucleus = parse_nmr_nucleus(nucleus)?;
    if !elements
        .iter()
        .any(|&z| z == parsed_nucleus.atomic_number())
    {
        return Err(format!(
            "Requested nucleus {} is not present in the provided element list.",
            parsed_nucleus.canonical()
        ));
    }

    let system = scf::types::MolecularSystem::from_angstrom(
        elements,
        positions,
        config.charge,
        config.multiplicity,
    );
    if !system.n_electrons().is_multiple_of(2) {
        return Err(format!(
            "The public GIAO NMR API currently supports even-electron closed-shell systems only; got {} electrons.",
            system.n_electrons()
        ));
    }

    let fallback_elements: Vec<u8> = elements
        .iter()
        .copied()
        .filter(|&z| !has_explicit_giao_sto3g_support(z))
        .collect();

    if !config.allow_basis_fallback && !fallback_elements.is_empty() {
        return Err(format!(
            "The public GIAO NMR SCF route only has explicit STO-3G basis data for {}. This system includes fallback-only elements: {}. Re-run with GiaoNmrConfig.allow_basis_fallback=true for a screening-only attempt.",
            format_element_set(EXPLICIT_GIAO_STO3G_ELEMENTS),
            format_element_set(&fallback_elements)
        ));
    }

    let basis = scf::basis::BasisSet::sto3g(&system.atomic_numbers, &system.positions_bohr);
    let n_occupied = system.n_electrons() / 2;
    if n_occupied > basis.n_basis {
        return Err(format!(
            "GIAO NMR preflight failed for {}: STO-3G supplies {} basis functions for {} occupied orbitals. This usually means one or more elements only have the hydrogen-like fallback basis in the current SCF path.",
            parsed_nucleus.canonical(),
            basis.n_basis,
            n_occupied
        ));
    }

    Ok(GiaoNmrPreflight {
        nucleus: parsed_nucleus,
        system,
        basis,
        fallback_elements,
    })
}

/// Compute isotope-specific GIAO NMR shieldings using the public SCF route.
pub fn compute_giao_nmr(
    elements: &[u8],
    positions: &[[f64; 3]],
    nucleus: &str,
) -> Result<GiaoNmrResult, String> {
    compute_giao_nmr_configured(elements, positions, nucleus, &GiaoNmrConfig::default())
}

/// Compute isotope-specific GIAO NMR shieldings using the public SCF route and explicit settings.
pub fn compute_giao_nmr_configured(
    elements: &[u8],
    positions: &[[f64; 3]],
    nucleus: &str,
    config: &GiaoNmrConfig,
) -> Result<GiaoNmrResult, String> {
    let request = preflight_giao_nmr(elements, positions, nucleus, config)?;

    let scf_config = scf::scf_loop::ScfConfig {
        max_iterations: config.max_scf_iterations,
        use_parallel_eri: config.use_parallel_eri,
        parallel_threshold: if config.use_parallel_eri { 0 } else { 20 },
        ..scf::scf_loop::ScfConfig::default()
    };

    let scf = scf::scf_loop::run_scf(&request.system, &scf_config);
    let shieldings = spectroscopy::compute_nmr_shieldings_for_nucleus(
        request.nucleus,
        &request.system.atomic_numbers,
        &request.system.positions_bohr,
        &spectroscopy::ScfInput::from(&scf),
        &request.basis.function_to_atom,
    );
    if shieldings.is_empty() {
        return Err(format!(
            "Requested nucleus {} is not present in the provided element list.",
            request.nucleus.canonical()
        ));
    }

    let entries: Vec<GiaoNmrEntry> = shieldings
        .iter()
        .map(|entry| GiaoNmrEntry {
            atom_index: entry.atom_index,
            element: entry.element,
            tensor: entry.tensor,
            isotropic: entry.isotropic,
            anisotropy: entry.anisotropy,
            chemical_shift: entry.chemical_shift,
        })
        .collect();

    let mut notes = vec![format!(
        "GIAO NMR for {} via the public SCF route (RHF/STO-3G).",
        request.nucleus.unicode_label()
    )];
    if scf.converged {
        notes.push(format!(
            "SCF converged in {} iterations with total energy {:.6} Hartree.",
            scf.scf_iterations, scf.total_energy
        ));
    } else {
        notes.push(format!(
            "SCF did not reach the requested threshold after {} iterations; treat the returned shieldings as screening-level only.",
            scf.scf_iterations
        ));
    }
    if request.fallback_elements.is_empty() {
        notes.push(
            "All elements in this system use explicit STO-3G basis data in the current SCF path."
                .to_string(),
        );
    } else {
        notes.push(format!(
            "Fallback basis enabled for {}. Heavy-element shieldings are qualitative in this mode.",
            format_element_set(&request.fallback_elements)
        ));
    }
    if request.nucleus.is_quadrupolar() {
        notes.push(
            "Quadrupolar nucleus selected: this API returns isotropic shieldings and chemical shifts only; relaxation-driven broadening still needs experimental lineshape modeling.".to_string(),
        );
    }
    if !request.nucleus.is_primary_target() {
        notes.push(
            "Non-1H/13C nuclei remain sensitive to the reference shielding and basis choice; benchmark the target family before quantitative use.".to_string(),
        );
    }

    Ok(GiaoNmrResult {
        nucleus: request.nucleus.canonical().to_string(),
        target_atomic_number: request.nucleus.atomic_number(),
        method: "RHF/GIAO".to_string(),
        basis_set: "STO-3G".to_string(),
        charge: config.charge,
        multiplicity: config.multiplicity,
        scf_converged: scf.converged,
        scf_iterations: scf.scf_iterations,
        total_energy_hartree: scf.total_energy,
        n_basis: scf.n_basis,
        n_target_atoms: entries.len(),
        chemical_shifts: entries.iter().map(|entry| entry.chemical_shift).collect(),
        shieldings: entries,
        fallback_elements: request.fallback_elements,
        notes,
    })
}

/// Generate a complete NMR spectrum from SMILES.
///
/// `nucleus`: any supported nucleus alias, for example "1H", "13C", "35Cl", or "195Pt"
/// `gamma`: Lorentzian line width in ppm
/// `ppm_min`, `ppm_max`: spectral window
/// `n_points`: grid resolution
pub fn compute_nmr_spectrum(
    smiles: &str,
    nucleus: &str,
    gamma: f64,
    ppm_min: f64,
    ppm_max: f64,
    n_points: usize,
) -> Result<nmr::NmrSpectrum, String> {
    compute_nmr_spectrum_with_coords(smiles, &[], nucleus, gamma, ppm_min, ppm_max, n_points)
}

/// Generate a complete NMR spectrum from SMILES and optional 3D coordinates.
/// When coordinates are provided, vicinal ³J couplings use the Karplus equation.
pub fn compute_nmr_spectrum_with_coords(
    smiles: &str,
    positions: &[[f64; 3]],
    nucleus: &str,
    gamma: f64,
    ppm_min: f64,
    ppm_max: f64,
    n_points: usize,
) -> Result<nmr::NmrSpectrum, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    if !positions.is_empty() {
        let elements: Vec<u8> = (0..mol.graph.node_count())
            .map(|i| mol.graph[petgraph::graph::NodeIndex::new(i)].element)
            .collect();
        validation::validate_geometry(&elements, positions)?;
    }
    let couplings = nmr::predict_j_couplings(&mol, positions);
    let nuc = parse_nmr_nucleus(nucleus)?;
    let shifts = nmr::predict_chemical_shifts_for_nucleus(&mol, nuc);
    Ok(nmr::spectrum::compute_nmr_spectrum_for_shifts(
        &shifts, &couplings, nuc, gamma, ppm_min, ppm_max, n_points,
    ))
}

/// Generate HOSE codes for all atoms in a molecule.
///
/// HOSE (homotopic environment for each atom) codes are the canonical
/// strings used by the nmrdb/ACD shift-prediction machinery to identify
/// chemically equivalent environments. `max_radius` is the maximum number
/// of bonds explored around each atom.
///
/// # Errors
///
/// Returns `Err` if the SMILES string cannot be parsed.
pub fn compute_hose_codes(smiles: &str, max_radius: usize) -> Result<Vec<nmr::HoseCode>, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    Ok(nmr::hose::generate_hose_codes(&mol, max_radius))
}

/// Compute orbital isosurface mesh for visualization.
///
/// Works with all implemented electronic-structure methods: EHT, PM3, xTB, HF-3c.
///
/// - `elements`: atomic numbers
/// - `positions`: Cartesian coordinates `[x,y,z]` in Å
/// - `method`: "eht", "pm3", "xtb", or "hf3c"
/// - `mo_index`: which molecular orbital to visualize
/// - `spacing`: grid spacing in Å (0.2 typical)
/// - `padding`: padding around molecule in Å (3.0 typical)
/// - `isovalue`: isosurface threshold (0.02 typical)
pub fn compute_orbital_mesh(
    elements: &[u8],
    positions: &[[f64; 3]],
    method: &str,
    mo_index: usize,
    spacing: f64,
    padding: f64,
    isovalue: f32,
) -> Result<mesh::OrbitalMeshResult, String> {
    let m = match method.to_lowercase().as_str() {
        "eht" | "huckel" => mesh::MeshMethod::Eht,
        "pm3" => mesh::MeshMethod::Pm3,
        "xtb" | "gfn0" | "gfn-xtb" => mesh::MeshMethod::Xtb,
        "hf3c" | "hf-3c" | "hf" => mesh::MeshMethod::Hf3c,
        _ => {
            return Err(format!(
                "Unknown method '{}'. Supported: eht, pm3, xtb, hf3c",
                method
            ))
        }
    };
    mesh::compute_orbital_mesh(elements, positions, m, mo_index, spacing, padding, isovalue)
}

/// Compute molecular descriptors for ML property prediction.
///
/// `elements`: atomic numbers.
/// `bonds`: (atom_i, atom_j, bond_order) list.
/// `charges`: partial charges (or empty slice for default).
/// `aromatic_atoms`: aromatic flags per atom (or empty slice).
pub fn compute_ml_descriptors(
    elements: &[u8],
    bonds: &[(usize, usize, u8)],
    charges: &[f64],
    aromatic_atoms: &[bool],
) -> ml::MolecularDescriptors {
    ml::compute_descriptors(elements, bonds, charges, aromatic_atoms)
}

/// Predict molecular properties using ML proxy models.
///
/// Returns LogP, molar refractivity, solubility, Lipinski flags,
/// and druglikeness score from molecular descriptors.
pub fn predict_ml_properties(desc: &ml::MolecularDescriptors) -> ml::MlPropertyResult {
    ml::predict_properties(desc)
}

/// Predict molecular properties using ensemble ML models.
///
/// Returns consensus LogP (with uncertainty), TPSA, pKa estimates,
/// Veber bioavailability rules, BBB permeability, and confidence score.
pub fn predict_ensemble(
    elements: &[u8],
    bonds: &[(usize, usize, u8)],
    charges: &[f64],
    aromatic_atoms: &[bool],
) -> ml::EnsembleResult {
    let desc = ml::compute_descriptors(elements, bonds, charges, aromatic_atoms);
    ml::predict_ensemble(&desc, elements, bonds)
}

/// Compute Topological Polar Surface Area (TPSA) in Å².
///
/// This is the Ertl fragment-contribution TPSA: polar fragments contribute
/// fixed area increments, so the result is a topological quantity and does
/// not depend on a 3D geometry.
pub fn compute_tpsa(elements: &[u8], bonds: &[(usize, usize, u8)]) -> f64 {
    ml::compute_tpsa(elements, bonds)
}

/// Compute the Topological Polar Surface Area (TPSA) in Å², from a SMILES.
///
/// The bond list is derived from the molecular graph, so no 3D geometry is
/// needed.
///
/// # Errors
///
/// Returns `Err` when the SMILES string cannot be parsed.
pub fn compute_tpsa_from_smiles(smiles: &str) -> Result<f64, String> {
    let (elements, bonds) = elements_and_bonds(smiles)?;
    Ok(ml::compute_tpsa(&elements, &bonds))
}

/// The elements and the bond list of a molecule, derived from its graph.
///
/// The bond order is encoded as 1 = single, 2 = double, 3 = triple and
/// 4 = aromatic.
type ElementsAndBonds = (Vec<u8>, Vec<(usize, usize, u8)>);

fn elements_and_bonds(smiles: &str) -> Result<ElementsAndBonds, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    let n = mol.graph.node_count();
    let elements: Vec<u8> = (0..n)
        .map(|i| mol.graph[graph::NodeIndex::new(i)].element)
        .collect();
    let bonds: Vec<(usize, usize, u8)> = mol
        .graph
        .edge_indices()
        .map(|e| {
            let (a, b) = mol.graph.edge_endpoints(e).expect("edge endpoints exist");
            let order = match mol.graph[e].order {
                graph::BondOrder::Single => 1u8,
                graph::BondOrder::Double => 2,
                graph::BondOrder::Triple => 3,
                graph::BondOrder::Aromatic => 4,
                graph::BondOrder::Unknown => 1,
            };
            (a.index(), b.index(), order)
        })
        .collect();
    Ok((elements, bonds))
}

/// Predict molecular properties using ensemble ML models, from a SMILES.
///
/// The bond list is derived from the molecular graph.
///
/// # Errors
///
/// Returns `Err` when the SMILES string cannot be parsed.
pub fn predict_ensemble_from_smiles(smiles: &str) -> Result<ml::EnsembleResult, String> {
    let (elements, bonds) = elements_and_bonds(smiles)?;
    Ok(predict_ensemble(&elements, &bonds, &[], &[]))
}

/// Run short exploratory molecular dynamics with Velocity Verlet (NVE-like).
///
/// `coords` is a flat `[x0,y0,z0, x1,y1,z1, ...]` array in Å, `n_steps` is
/// the number of integration steps, `dt_fs` the timestep in femtoseconds and
/// `seed` the RNG seed for the initial Maxwell-Boltzmann velocities.
///
/// The trajectory stores positions in Å, times in fs, velocities in Å/fs,
/// kinetic and potential energies in kcal/mol and temperatures in K.
///
/// # Errors
///
/// Returns `Err` if the SMILES cannot be parsed, the coordinate array does
/// not match the atom count, or the UFF force field cannot be built for the
/// molecule.
pub fn compute_md_trajectory(
    smiles: &str,
    coords: &[f64],
    n_steps: usize,
    dt_fs: f64,
    seed: u64,
) -> Result<dynamics::MdTrajectory, String> {
    dynamics::simulate_velocity_verlet_uff(smiles, coords, n_steps, dt_fs, seed, None)
}

/// Run short exploratory molecular dynamics with Velocity Verlet + Berendsen NVT thermostat.
///
/// In addition to the arguments of [`compute_md_trajectory`], `target_temp_k`
/// is the thermostat set point in K and `thermostat_tau_fs` the coupling
/// time in fs. The Berendsen thermostat does not produce a canonical
/// ensemble, so the reported temperature fluctuates around the set point
/// rather than sampling it exactly.
pub fn compute_md_trajectory_nvt(
    smiles: &str,
    coords: &[f64],
    n_steps: usize,
    dt_fs: f64,
    seed: u64,
    target_temp_k: f64,
    thermostat_tau_fs: f64,
) -> Result<dynamics::MdTrajectory, String> {
    dynamics::simulate_velocity_verlet_uff(
        smiles,
        coords,
        n_steps,
        dt_fs,
        seed,
        Some((target_temp_k, thermostat_tau_fs)),
    )
}

/// Build a simplified NEB path between two geometries.
///
/// `start_coords` and `end_coords` are flat `[x,y,z,...]` arrays in Å,
/// `n_images` is the number of images between the two endpoints (including
/// them), `n_iter` the number of relaxation iterations per image,
/// `spring_k` the harmonic spring constant in kcal/mol/Å² and `step_size`
/// the gradient step in Å.
///
/// The returned path energies are in kcal/mol relative to the start image.
pub fn compute_simplified_neb_path(
    smiles: &str,
    start_coords: &[f64],
    end_coords: &[f64],
    n_images: usize,
    n_iter: usize,
    spring_k: f64,
    step_size: f64,
) -> Result<dynamics::NebPathResult, String> {
    dynamics::compute_simplified_neb_path(
        smiles,
        start_coords,
        end_coords,
        n_images,
        n_iter,
        spring_k,
        step_size,
    )
}

/// Build a simplified NEB path using a configurable energy backend.
///
/// `method`: one of `"uff"`, `"mmff94"`, `"pm3"`, `"xtb"`, `"gfn1"`, `"gfn2"`, `"hf3c"`.
pub fn compute_simplified_neb_path_configurable(
    smiles: &str,
    start_coords: &[f64],
    end_coords: &[f64],
    n_images: usize,
    n_iter: usize,
    spring_k: f64,
    step_size: f64,
    method: &str,
) -> Result<dynamics::NebPathResult, String> {
    dynamics::compute_simplified_neb_path_configurable(
        smiles,
        start_coords,
        end_coords,
        n_images,
        n_iter,
        spring_k,
        step_size,
        method,
    )
}

/// Compute single-point energy (kcal/mol) with any NEB-capable backend.
///
/// `method`: one of `"uff"`, `"mmff94"`, `"pm3"`, `"xtb"`, `"gfn1"`, `"gfn2"`, `"hf3c"`.
///
/// # Errors
///
/// Returns `Err` when the method name is not recognised, when the SMILES
/// string cannot be parsed, when the coordinate array does not match the atom
/// count, or when the selected backend fails.
pub fn neb_backend_energy_kcal(method: &str, smiles: &str, coords: &[f64]) -> Result<f64, String> {
    dynamics::neb_backend_energy_kcal(method, smiles, coords)
}

/// Compute energy (kcal/mol) and gradient (kcal/mol/Å) with any NEB-capable backend.
///
/// Returns `(energy, gradient_flat)`.
pub fn neb_backend_energy_and_gradient(
    method: &str,
    smiles: &str,
    coords: &[f64],
) -> Result<(f64, Vec<f64>), String> {
    dynamics::neb_backend_energy_and_gradient(method, smiles, coords)
}

/// Compute a full reaction dynamics path: embed reactant + product fragments,
/// build oriented complexes, run NEB for the reactive zone, and generate
/// approach/departure frames — all with energies computed in Rust.
///
/// See [`dynamics::ReactionDynamicsConfig`] for tuning parameters.
pub fn compute_reaction_dynamics(
    reactant_smiles: &[&str],
    product_smiles: &[&str],
    method: &str,
    config: &dynamics::ReactionDynamicsConfig,
) -> Result<dynamics::ReactionDynamicsResult, String> {
    dynamics::compute_reaction_dynamics(reactant_smiles, product_smiles, method, config)
}

fn coords_flat_to_matrix_f32(coords: &[f64], n_atoms: usize) -> nalgebra::DMatrix<f32> {
    let mut m = nalgebra::DMatrix::<f32>::zeros(n_atoms, 3);
    for i in 0..n_atoms {
        m[(i, 0)] = coords[3 * i] as f32;
        m[(i, 1)] = coords[3 * i + 1] as f32;
        m[(i, 2)] = coords[3 * i + 2] as f32;
    }
    m
}

fn coords_matrix_f32_to_flat(m: &nalgebra::DMatrix<f32>) -> Vec<f64> {
    let mut out = Vec::with_capacity(m.nrows() * 3);
    for i in 0..m.nrows() {
        out.push(m[(i, 0)] as f64);
        out.push(m[(i, 1)] as f64);
        out.push(m[(i, 2)] as f64);
    }
    out
}

/// Search conformers by sampling multiple embeddings, optimizing torsions, and ranking with UFF.
///
/// This workflow:
/// 1. Generates up to `n_samples` conformers with different seeds.
/// 2. Applies Monte Carlo torsion sampling plus a greedy rotatable-bond refinement pass.
/// 3. Scores each conformer with UFF energy (kcal/mol).
/// 4. Filters near-duplicates using RMSD thresholding.
/// 5. Builds explicit RMSD clusters and returns representative/member mapping.
pub fn search_conformers_with_uff(
    smiles: &str,
    n_samples: usize,
    seed: u64,
    rmsd_threshold: f64,
) -> Result<ConformerSearchResult, String> {
    if n_samples == 0 {
        return Err("n_samples must be > 0".to_string());
    }
    if rmsd_threshold <= 0.0 {
        return Err("rmsd_threshold must be > 0".to_string());
    }

    let mol = graph::Molecule::from_smiles(smiles)?;
    let n_atoms = mol.graph.node_count();
    let bounds = distgeom::smooth_bounds_matrix(distgeom::calculate_bounds_matrix(&mol));

    let mut generated = Vec::new();
    let mut notes = Vec::new();
    let mut rotatable_bonds = 0usize;

    for i in 0..n_samples {
        let sample_seed = seed.wrapping_add(i as u64);
        let conf = embed(smiles, sample_seed);

        if conf.error.is_some() || conf.coords.len() != n_atoms * 3 {
            continue;
        }

        let mut coords = coords_flat_to_matrix_f32(&conf.coords, n_atoms);
        let rot_mc = forcefield::optimize_torsions_monte_carlo_bounds(
            &mut coords,
            &mol,
            &bounds,
            sample_seed ^ 0x9E37_79B9_7F4A_7C15,
            64,
            0.4,
        );
        let rot_greedy = forcefield::optimize_torsions_bounds(&mut coords, &mol, &bounds, 2);
        let rot = rot_mc.max(rot_greedy);
        rotatable_bonds = rot;
        let coords_flat = coords_matrix_f32_to_flat(&coords);

        match compute_uff_energy(smiles, &coords_flat) {
            Ok(energy_kcal_mol) => generated.push(ConformerEnsembleMember {
                seed: sample_seed,
                cluster_id: None,
                coords: coords_flat,
                energy_kcal_mol,
            }),
            Err(_) => continue,
        }
    }

    if generated.is_empty() {
        return Err("failed to generate any valid conformers".to_string());
    }

    generated.sort_by(|a, b| {
        a.energy_kcal_mol
            .partial_cmp(&b.energy_kcal_mol)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let generated_count = generated.len();

    let mut unique = Vec::new();
    let mut cluster_members: Vec<Vec<u64>> = Vec::new();
    for candidate in generated {
        let existing_cluster = unique.iter().position(|u: &ConformerEnsembleMember| {
            compute_rmsd(&candidate.coords, &u.coords) < rmsd_threshold
        });

        if let Some(cluster_id) = existing_cluster {
            cluster_members[cluster_id].push(candidate.seed);
        } else {
            unique.push(candidate.clone());
            cluster_members.push(vec![candidate.seed]);
        }
    }

    for (cluster_id, member) in unique.iter_mut().enumerate() {
        member.cluster_id = Some(cluster_id);
    }

    let clusters: Vec<ConformerClusterSummary> = unique
        .iter()
        .enumerate()
        .map(|(cluster_id, representative)| ConformerClusterSummary {
            cluster_id,
            representative_seed: representative.seed,
            size: cluster_members[cluster_id].len(),
            member_seeds: cluster_members[cluster_id].clone(),
        })
        .collect();

    notes.push(
        "Conformers are preconditioned with Monte Carlo torsion sampling + greedy torsion refinement, ranked by UFF energy, deduplicated by Kabsch-aligned RMSD threshold, and summarized as explicit RMSD clusters."
            .to_string(),
    );

    Ok(ConformerSearchResult {
        generated: generated_count,
        unique: unique.len(),
        rotatable_bonds,
        conformers: unique,
        clusters,
        notes,
    })
}

/// Compute UFF energy and apply an aromaticity-informed heuristic correction.
pub fn compute_uff_energy_with_aromatic_heuristics(
    smiles: &str,
    coords: &[f64],
) -> Result<reactivity::UffHeuristicEnergy, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    let n = mol.graph.node_count();
    if coords.len() != n * 3 {
        return Err(format!("coords length {} != 3 * atoms {}", coords.len(), n));
    }

    let ff = forcefield::builder::build_uff_force_field(&mol);
    let mut gradient = vec![0.0f64; n * 3];
    let raw = ff.compute_system_energy_and_gradients(coords, &mut gradient);
    Ok(reactivity::apply_aromatic_uff_correction(&mol, raw))
}

/// Estimate acidic/basic pKa sites from graph environments and Gasteiger-charge heuristics.
///
/// The returned pKa values are in pKa units (dimensionless, the negative log
/// of the acid dissociation constant) and are heuristic estimates, not
/// measured or computed values.
///
/// # Errors
///
/// Returns `Err` when the SMILES string cannot be parsed or when the
/// Gasteiger charge calculation fails.
pub fn compute_empirical_pka(smiles: &str) -> Result<reactivity::EmpiricalPkaResult, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    let charges = compute_charges(smiles)?;
    Ok(reactivity::estimate_empirical_pka(&mol, &charges.charges))
}

/// Create a periodic unit cell from lattice parameters (a, b, c in Å; α, β, γ in degrees).
pub fn create_unit_cell(
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
) -> materials::UnitCell {
    materials::UnitCell::from_parameters(&materials::CellParameters {
        a,
        b,
        c,
        alpha,
        beta,
        gamma,
    })
}

/// Assemble a framework crystal structure from node/linker SBUs on a topology.
///
/// Returns the assembled crystal structure as a JSON-serializable `CrystalStructure`.
pub fn assemble_framework(
    node: &materials::Sbu,
    linker: &materials::Sbu,
    topology: &materials::Topology,
    cell: &materials::UnitCell,
) -> materials::CrystalStructure {
    materials::assemble_framework(node, linker, topology, cell)
}

// ─── CIF Import/Export ────────────────────────────────────────────────────────

/// Parse a CIF (Crystallographic Information File) string into a crystal structure.
///
/// Returns cell parameters, space group info, and atom site positions.
///
/// # Errors
///
/// Returns `Err` when the CIF text cannot be parsed, when the cell parameters
/// are missing or inconsistent, or when the atom-site loop is malformed.
pub fn parse_cif(input: &str) -> Result<materials::CifStructure, String> {
    materials::cif::parse_cif(input)
}

/// Write a crystal structure to CIF format string.
pub fn write_cif(structure: &materials::CifStructure) -> String {
    materials::cif::write_cif(structure)
}

// ─── UHF / ROHF — Open-Shell SCF ──────────────────────────────────────────────

/// Run an Unrestricted Hartree-Fock (UHF) calculation for open-shell systems.
///
/// `elements`: atomic numbers.
/// `positions`: Cartesian coordinates in Å, one `[x,y,z]` per atom.
/// `charge`: total molecular charge (0 = neutral).
/// `multiplicity`: spin multiplicity (1 = singlet, 2 = doublet, 3 = triplet, …).
///
/// Returns separate alpha/beta orbitals, <S²> expectation, spin contamination.
pub fn compute_uhf(
    elements: &[u8],
    positions: &[[f64; 3]],
    charge: i32,
    multiplicity: u32,
) -> Result<scf::uhf::UhfResult, String> {
    let system =
        scf::types::MolecularSystem::try_from_angstrom(elements, positions, charge, multiplicity)?;
    scf::uhf::run_uhf_parallel(&system)
}

/// Run a Restricted Open-shell Hartree-Fock (ROHF) calculation.
///
/// Uses canonical ROHF projection: same spatial orbitals for both spins,
/// eliminating spin contamination. Suited for clean-spin-state calculations.
///
/// `elements`: atomic numbers.
/// `positions`: Cartesian coordinates in Å, one `[x,y,z]` per atom.
/// `charge`: total molecular charge.
/// `multiplicity`: spin multiplicity.
pub fn compute_rohf(
    elements: &[u8],
    positions: &[[f64; 3]],
    charge: i32,
    multiplicity: u32,
) -> Result<scf::uhf::UhfResult, String> {
    let system =
        scf::types::MolecularSystem::try_from_angstrom(elements, positions, charge, multiplicity)?;
    scf::uhf::run_rohf(&system)
}

/// Run UHF/ROHF with full configuration control.
pub fn compute_uhf_configured(
    elements: &[u8],
    positions: &[[f64; 3]],
    charge: i32,
    multiplicity: u32,
    config: &scf::uhf::UhfConfig,
) -> Result<scf::uhf::UhfResult, String> {
    let system =
        scf::types::MolecularSystem::try_from_angstrom(elements, positions, charge, multiplicity)?;
    scf::uhf::run_uhf(&system, config)
}

/// Run a complete HF-3c calculation (Hartree-Fock with D3, gCP, SRB corrections).
///
/// `elements`: atomic numbers.
/// `positions`: Cartesian coordinates in Å, one `[x,y,z]` per atom.
/// `config`: calculation parameters (or use `HfConfig::default()`).
///
/// Returns orbital energies, total energy, correction energies, and optional CIS states.
pub fn compute_hf3c(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &hf::HfConfig,
) -> Result<hf::Hf3cResult, String> {
    validation::validate_geometry(elements, positions)?;
    hf::solve_hf3c(elements, positions, config)
}

/// Compute ANI neural-network potential energy (and optionally forces).
///
/// Uses internally-generated test weights — suitable for testing and demonstration.
/// For physically meaningful results, use `ani::compute_ani()` with trained weights.
///
/// `elements`: atomic numbers (supported: H, C, N, O, F, S, Cl).
/// `positions`: Cartesian coordinates in Å, one `[x,y,z]` per atom.
///
/// The energy is in kcal/mol and the forces (when requested) in kcal/mol/Å.
///
/// # Errors
///
/// Returns `Err` when the geometry fails validation or when an element has
/// no ANI network.
pub fn compute_ani(elements: &[u8], positions: &[[f64; 3]]) -> Result<ani::AniResult, String> {
    validation::validate_geometry(elements, positions)?;
    ani::api::compute_ani_test(elements, positions)
}

/// Compute ANI neural-network potential energy with custom config and pre-loaded models.
///
/// `elements`: atomic numbers.
/// `positions`: Cartesian coordinates in Å.
/// `config`: ANI configuration (cutoff, force computation flag).
/// `models`: pre-loaded element→network map from `ani::weights::load_weights()`.
pub fn compute_ani_with_models(
    elements: &[u8],
    positions: &[[f64; 3]],
    config: &ani::AniConfig,
    models: &std::collections::HashMap<u8, ani::nn::FeedForwardNet>,
) -> Result<ani::AniResult, String> {
    validation::validate_geometry(elements, positions)?;
    ani::compute_ani(elements, positions, config, models)
}

/// Compute ESP grid and return full result with values, origin, spacing and dimensions.
///
/// This is a convenience alias for `compute_esp()` returning the `EspGrid` struct.
pub fn compute_esp_grid(
    elements: &[u8],
    positions: &[[f64; 3]],
    spacing: f64,
    padding: f64,
) -> Result<esp::EspGrid, String> {
    compute_esp(elements, positions, spacing, padding)
}

// ─── Canonical SMILES API ─────────────────────────────────────────────────────

/// Generate a canonical SMILES string from an input SMILES.
///
/// Produces a deterministic SMILES representation regardless of
/// input atom ordering (e.g., "CCO" and "OCC" produce the same output).
///
/// # Errors
///
/// Returns `Err` when the SMILES string cannot be parsed.
pub fn to_canonical_smiles(smiles: &str) -> Result<String, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    Ok(canonical_smiles::to_canonical_smiles(&mol))
}

// ─── Stereochemistry API ─────────────────────────────────────────────────────

/// Analyze stereochemistry: detect chiral centers (R/S) and E/Z double bonds.
///
/// `coords` is a flat `[x0, y0, z0, x1, y1, z1, ...]` array in Å. When the
/// coordinates are empty the analysis is graph-based only and the 3D
/// configuration of double bonds cannot be assigned.
///
/// # Errors
///
/// Returns `Err` when the SMILES string cannot be parsed.
pub fn analyze_stereo(smiles: &str, coords: &[f64]) -> Result<stereo::StereoAnalysis, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    let positions: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    Ok(stereo::analyze_stereo(&mol, &positions))
}

// ─── Solvation API ───────────────────────────────────────────────────────────

/// Compute non-polar solvation energy from SASA and atomic solvation parameters.
///
/// The returned energy is in kcal/mol and is a sum of `gamma_i * SASA_i` terms
/// with atomic surface tensions `gamma_i` in kcal/(mol·Å²).
///
/// # Errors
///
/// Returns `Err` when the geometry is invalid (mismatched lengths, empty
/// system, non-finite coordinates) or when the probe radius is negative or
/// non-finite.
pub fn compute_nonpolar_solvation(
    elements: &[u8],
    positions: &[[f64; 3]],
    probe_radius: Option<f64>,
) -> Result<solvation::NonPolarSolvation, String> {
    solvation::compute_nonpolar_solvation(elements, positions, probe_radius)
}

/// Compute Generalized Born (Still-HCT) electrostatic + non-polar solvation.
///
/// The non-polar term is `surface_tension · SASA` with the shared
/// Shrake-Rupley implementation ([`surface::sasa`]). For the ALPB-corrected
/// variant (xTB-compatible) use [`compute_alpb_solvation`].
///
/// # Errors
///
/// Returns `Err` when the geometry is invalid, when the charge vector does
/// not have one finite charge per atom, when a dielectric constant is
/// non-finite or below 1.0, or when the probe radius is negative or
/// non-finite.
pub fn compute_gb_solvation(
    elements: &[u8],
    positions: &[[f64; 3]],
    charges: &[f64],
    solvent_dielectric: Option<f64>,
    solute_dielectric: Option<f64>,
    probe_radius: Option<f64>,
) -> Result<solvation::GbSolvation, String> {
    solvation::compute_gb_solvation(
        elements,
        positions,
        charges,
        solvent_dielectric,
        solute_dielectric,
        probe_radius,
    )
}

/// Compute ALPB (Analytical Linearized Poisson-Boltzmann) solvation energy.
///
/// ALPB extends the same Generalized Born electrostatics with the Klamt
/// correction `(ε-1)/(ε + A·x)` and the shared SASA non-polar term. Use it
/// instead of [`compute_gb_solvation`] when matching xTB/ALPB references.
///
/// # Errors
///
/// Returns `Err` when the geometry is invalid, when the charge vector does
/// not have one finite charge per atom, when the solvent dielectric constant
/// is non-finite or below 1.0, when the probe radius is negative or
/// non-finite, or when the surface tension is not finite.
pub fn compute_alpb_solvation(
    elements: &[u8],
    positions: &[[f64; 3]],
    charges: &[f64],
    config: solvation_alpb::AlpbConfig,
) -> Result<solvation_alpb::AlpbResult, String> {
    solvation_alpb::compute_alpb_solvation(elements, positions, charges, &config)
}

// ─── Clustering API ──────────────────────────────────────────────────────────

/// Cluster conformers using Butina (Taylor-Butina) algorithm based on RMSD.
pub fn butina_cluster(conformers: &[Vec<f64>], rmsd_cutoff: f64) -> clustering::ClusterResult {
    clustering::butina_cluster(conformers, rmsd_cutoff)
}

/// Compute all-pairs RMSD matrix for a set of conformers.
pub fn compute_rmsd_matrix(conformers: &[Vec<f64>]) -> Vec<Vec<f64>> {
    clustering::compute_rmsd_matrix(conformers)
}

// ─── Rings & Fingerprints API ────────────────────────────────────────────────

/// Compute the Smallest Set of Smallest Rings (SSSR).
///
/// The returned rings are lists of 0-based atom indices; no units are involved.
///
/// # Errors
///
/// Returns `Err` when the SMILES string cannot be parsed.
pub fn compute_sssr(smiles: &str) -> Result<rings::sssr::SssrResult, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    Ok(rings::sssr::compute_sssr(&mol))
}

/// Compute Extended-Connectivity Fingerprint (ECFP/Morgan).
pub fn compute_ecfp(
    smiles: &str,
    radius: usize,
    n_bits: usize,
) -> Result<rings::ecfp::ECFPFingerprint, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    Ok(rings::ecfp::compute_ecfp(&mol, radius, n_bits))
}

/// Compute ECFP fingerprints for a batch of SMILES (parallelized with rayon when enabled).
pub fn compute_ecfp_batch(
    smiles_list: &[&str],
    radius: usize,
    n_bits: usize,
) -> Result<Vec<rings::ecfp::ECFPFingerprint>, String> {
    let mols: Result<Vec<_>, _> = smiles_list
        .iter()
        .map(|s| graph::Molecule::from_smiles(s))
        .collect();
    let mols = mols?;
    Ok(rings::ecfp::compute_ecfp_batch(&mols, radius, n_bits))
}

/// Compute Tanimoto similarity between two ECFP fingerprints.
///
/// The Tanimoto (Jaccard) coefficient of two bit vectors is
/// `|A ∩ B| / |A ∪ B|`; the returned similarity is dimensionless and lies in
/// `[0, 1]`, with `1.0` meaning identical fingerprints.
pub fn compute_tanimoto(
    fp1: &rings::ecfp::ECFPFingerprint,
    fp2: &rings::ecfp::ECFPFingerprint,
) -> f64 {
    rings::ecfp::compute_tanimoto(fp1, fp2)
}

// ─── Charges configured API ─────────────────────────────────────────────────

/// Compute Gasteiger-Marsili charges with configurable parameters.
pub fn compute_charges_configured(
    smiles: &str,
    config: &charges::gasteiger::GasteigerConfig,
) -> Result<charges::gasteiger::ChargeResult, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    let n = mol.graph.node_count();
    let elements: Vec<u8> = (0..n)
        .map(|i| mol.graph[petgraph::graph::NodeIndex::new(i)].element)
        .collect();
    let bonds: Vec<(usize, usize, u8)> = mol
        .graph
        .edge_indices()
        .map(|e| {
            let (a, b) = mol.graph.edge_endpoints(e).unwrap();
            let order = match mol.graph[e].order {
                graph::BondOrder::Single => 1u8,
                graph::BondOrder::Double => 2,
                graph::BondOrder::Triple => 3,
                graph::BondOrder::Aromatic => 4,
                graph::BondOrder::Unknown => 1,
            };
            (a.index(), b.index(), order)
        })
        .collect();
    let formal_charges: Vec<i8> = (0..n)
        .map(|i| mol.graph[petgraph::graph::NodeIndex::new(i)].formal_charge)
        .collect();
    let modes: Vec<charges::gasteiger::GasteigerMode> = (0..n)
        .map(
            |i| match mol.graph[petgraph::graph::NodeIndex::new(i)].hybridization {
                graph::Hybridization::SP => charges::gasteiger::GasteigerMode::Sp,
                graph::Hybridization::SP2 => charges::gasteiger::GasteigerMode::Sp2,
                graph::Hybridization::SP3 => charges::gasteiger::GasteigerMode::Sp3,
                _ => charges::gasteiger::GasteigerMode::Default,
            },
        )
        .collect();
    charges::gasteiger::gasteiger_marsili_charges_full(
        &elements,
        &bonds,
        &formal_charges,
        Some(&modes),
        config,
    )
}

// ─── NMR Enhanced API ────────────────────────────────────────────────────────

/// Predict J-couplings averaged over a conformer ensemble using Boltzmann weighting.
pub fn compute_ensemble_j_couplings(
    smiles: &str,
    conformer_coords: &[Vec<f64>],
    energies_kcal: &[f64],
    temperature_k: f64,
) -> Result<Vec<nmr::JCoupling>, String> {
    let mol = graph::Molecule::from_smiles(smiles)?;
    let positions: Vec<Vec<[f64; 3]>> = conformer_coords
        .iter()
        .map(|c| c.chunks(3).map(|p| [p[0], p[1], p[2]]).collect())
        .collect();
    Ok(nmr::coupling::ensemble_averaged_j_couplings(
        &mol,
        &positions,
        energies_kcal,
        temperature_k,
    ))
}

// ─── IR Enhanced API ─────────────────────────────────────────────────────────

/// Compute IR spectrum with configurable broadening type.
pub fn compute_ir_spectrum_broadened(
    analysis: &ir::VibrationalAnalysis,
    gamma: f64,
    wn_min: f64,
    wn_max: f64,
    n_points: usize,
    broadening: &str,
) -> ir::IrSpectrum {
    let bt = match broadening {
        "gaussian" | "Gaussian" => ir::BroadeningType::Gaussian,
        _ => ir::BroadeningType::Lorentzian,
    };
    ir::vibrations::compute_ir_spectrum_with_broadening(
        analysis, gamma, wn_min, wn_max, n_points, bt,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cisplatin_style_support_metadata() {
        let smiles = "[Pt](Cl)(Cl)([NH3])[NH3]";
        let mol = parse(smiles).expect("Cisplatin-style example should parse");
        let elements: Vec<u8> = (0..mol.graph.node_count())
            .map(|i| mol.graph[petgraph::graph::NodeIndex::new(i)].element)
            .collect();

        let caps = system_capabilities(&elements);
        assert!(caps.eht.available);
        assert!(matches!(
            caps.eht.confidence,
            eht::SupportLevel::Experimental
        ));
    }

    #[test]
    fn test_pt_system_routes_to_uff_when_experimental_disabled() {
        let smiles = "[Pt](Cl)(Cl)([NH3])[NH3]";
        let mol = parse(smiles).expect("Pt example should parse");
        let n = mol.graph.node_count();
        let elements: Vec<u8> = (0..n)
            .map(|i| mol.graph[petgraph::graph::NodeIndex::new(i)].element)
            .collect();
        let positions = vec![[0.0, 0.0, 0.0]; n];

        let result = compute_eht_or_uff_fallback(smiles, &elements, &positions, false)
            .expect("Fallback workflow should return a result");

        assert!(matches!(
            result,
            ElectronicWorkflowResult::UffFallback { .. }
        ));
    }

    #[test]
    fn test_method_plan_prefers_supported_workflows_for_organic_systems() {
        let plan = system_method_plan(&[6, 1, 1, 1, 1]);

        assert_eq!(plan.geometry.recommended, Some(ScientificMethod::Embed));
        assert_eq!(
            plan.force_field_energy.recommended,
            Some(ScientificMethod::Uff)
        );
        assert_eq!(plan.orbitals.recommended, Some(ScientificMethod::Eht));
        assert_eq!(plan.population.recommended, Some(ScientificMethod::Eht));
        assert_eq!(plan.orbital_grid.recommended, Some(ScientificMethod::Eht));
        assert!(plan.orbitals.methods[0].confidence_score > 0.9);
    }

    #[test]
    fn test_method_plan_marks_metal_orbital_workflow_experimental() {
        let plan = system_method_plan(&[78, 17, 17, 7, 7]);

        assert_eq!(plan.orbitals.recommended, Some(ScientificMethod::Eht));
        assert_eq!(
            plan.force_field_energy.recommended,
            Some(ScientificMethod::Uff)
        );
        assert!(matches!(
            plan.orbitals.methods[0].confidence,
            eht::SupportLevel::Experimental
        ));
        assert!(plan.orbitals.methods[0].confidence_score < 0.9);
        assert!(!plan.orbitals.methods[0].warnings.is_empty());
    }

    #[test]
    fn test_method_plan_reports_unavailable_workflows_for_unsupported_elements() {
        let plan = system_method_plan(&[104]);

        assert_eq!(plan.force_field_energy.recommended, None);
        assert_eq!(plan.orbitals.recommended, None);
        assert_eq!(plan.population.recommended, None);
        assert_eq!(plan.orbital_grid.recommended, None);
        assert!(!plan.orbitals.methods[0].limitations.is_empty());
    }

    #[test]
    fn test_compare_methods_supported_system_returns_success_rows() {
        let result = compare_methods("CC", &[6, 6], &[[0.0, 0.0, 0.0], [1.54, 0.0, 0.0]], true);
        assert_eq!(result.comparisons.len(), 2);
        assert!(result
            .comparisons
            .iter()
            .any(|entry| matches!(entry.method, ScientificMethod::Uff) && entry.available));
        assert!(result.comparisons.iter().any(|entry| matches!(
            entry.method,
            ScientificMethod::Eht
        ) && matches!(
            entry.status,
            MethodComparisonStatus::Success
        )));
    }

    #[test]
    fn test_compare_methods_blocks_experimental_eht_when_disabled() {
        let result = compare_methods("[O]", &[78], &[[0.0, 0.0, 0.0]], false);
        let eht_row = result
            .comparisons
            .iter()
            .find(|entry| matches!(entry.method, ScientificMethod::Eht))
            .expect("EHT row must exist");
        assert!(matches!(
            eht_row.status,
            MethodComparisonStatus::Unavailable
        ));
    }

    #[test]
    fn test_compare_methods_reports_unavailable_for_unsupported_elements() {
        let result = compare_methods("[O]", &[104], &[[0.0, 0.0, 0.0]], true);
        assert!(result
            .comparisons
            .iter()
            .all(|entry| matches!(entry.status, MethodComparisonStatus::Unavailable)));
    }

    #[test]
    fn test_compute_fukui_descriptors_returns_atomwise_output() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let result = compute_fukui_descriptors(&elements, &positions).unwrap();
        assert_eq!(result.num_atoms, 3);
        assert_eq!(result.condensed.len(), 3);
    }

    #[test]
    fn test_compute_uv_vis_spectrum_returns_requested_grid_size() {
        let elements = [6u8, 6, 1, 1, 1, 1];
        let positions = [
            [0.0, 0.0, 0.0],
            [1.34, 0.0, 0.0],
            [-0.6, 0.92, 0.0],
            [-0.6, -0.92, 0.0],
            [1.94, 0.92, 0.0],
            [1.94, -0.92, 0.0],
        ];
        let spectrum = compute_uv_vis_spectrum(&elements, &positions, 0.2, 0.5, 8.0, 256).unwrap();
        assert_eq!(spectrum.energies_ev.len(), 256);
        assert_eq!(spectrum.intensities.len(), 256);
    }

    #[test]
    fn test_analyze_graph_features_reports_benzene_aromaticity() {
        let analysis = analyze_graph_features("c1ccccc1").unwrap();
        assert_eq!(analysis.aromaticity.aromatic_atoms.len(), 12);
        assert_eq!(
            analysis
                .aromaticity
                .aromatic_atoms
                .iter()
                .filter(|v| **v)
                .count(),
            6
        );
        assert_eq!(analysis.aromaticity.aromatic_bonds.len(), 6);
    }

    #[test]
    fn test_compute_empirical_pka_finds_acidic_site_for_acetic_acid() {
        let result = compute_empirical_pka("CC(=O)O").unwrap();
        assert!(!result.acidic_sites.is_empty());
    }

    #[test]
    fn test_compute_uff_energy_with_aromatic_heuristics_applies_correction() {
        let conf = embed("c1ccccc1", 42);
        assert!(conf.error.is_none());

        let result = compute_uff_energy_with_aromatic_heuristics("c1ccccc1", &conf.coords).unwrap();
        assert!(result.aromatic_bond_count >= 6);
        assert!(result.corrected_energy_kcal_mol <= result.raw_energy_kcal_mol);
    }

    #[test]
    fn test_uff_energy_breakdown_sums_to_total() {
        let conf = embed("CCO", 42);
        assert!(conf.error.is_none());
        let breakdown = compute_uff_energy_breakdown("CCO", &conf.coords).unwrap();
        let total = compute_uff_energy("CCO", &conf.coords).unwrap();
        let sum: f64 = breakdown.iter().map(|(_, e)| e).sum();
        assert!(
            (sum - total).abs() < 1e-9,
            "breakdown {sum} vs total {total}: {breakdown:?}"
        );
        let kinds: Vec<&str> = breakdown.iter().map(|(k, _)| k.as_str()).collect();
        for expected in ["bond", "angle", "torsion", "vdw"] {
            assert!(kinds.contains(&expected), "missing {expected}: {kinds:?}");
        }
    }

    #[test]
    fn test_mmff94_energy_breakdown_sums_to_total() {
        let conf = embed("CCO", 42);
        assert!(conf.error.is_none());
        let breakdown = compute_mmff94_energy_breakdown("CCO", &conf.coords).unwrap();
        let total = compute_mmff94_energy("CCO", &conf.coords).unwrap();
        let sum: f64 = breakdown.iter().map(|(_, e)| e).sum();
        assert!(
            (sum - total).abs() < 1e-9,
            "breakdown {sum} vs total {total}: {breakdown:?}"
        );
        assert!(!breakdown.is_empty());
    }

    #[test]
    fn test_search_conformers_with_uff_returns_ranked_unique_ensemble() {
        let result = search_conformers_with_uff("CCCC", 10, 42, 0.2).unwrap();
        assert!(result.generated >= 1);
        assert!(result.unique >= 1);
        assert_eq!(result.unique, result.conformers.len());
        assert_eq!(result.unique, result.clusters.len());
        assert!(result.rotatable_bonds >= 1);

        let mut total_members = 0usize;
        for (i, cluster) in result.clusters.iter().enumerate() {
            assert_eq!(cluster.cluster_id, i);
            assert!(cluster.size >= 1);
            total_members += cluster.size;
            assert_eq!(result.conformers[i].cluster_id, Some(i));
            assert_eq!(result.conformers[i].seed, cluster.representative_seed);
        }
        assert_eq!(total_members, result.generated);

        for i in 1..result.conformers.len() {
            assert!(
                result.conformers[i - 1].energy_kcal_mol <= result.conformers[i].energy_kcal_mol
            );
        }
    }

    #[test]
    fn test_search_conformers_with_uff_large_rmsd_threshold_collapses_duplicates() {
        let result = search_conformers_with_uff("CCCC", 8, 123, 10.0).unwrap();
        assert_eq!(result.unique, 1);
        assert_eq!(result.clusters.len(), 1);
        assert_eq!(result.clusters[0].size, result.generated);
    }

    #[test]
    fn test_compute_md_trajectory_velocity_verlet_runs() {
        let conf = embed("CC", 42);
        assert!(conf.error.is_none());

        let trj = compute_md_trajectory("CC", &conf.coords, 10, 0.25, 7).unwrap();
        assert_eq!(trj.frames.len(), 11);
        assert!(trj
            .frames
            .iter()
            .all(|f| f.coords.iter().all(|v| v.is_finite())));
    }

    #[test]
    fn test_compute_md_trajectory_nvt_runs() {
        let conf = embed("CCO", 42);
        assert!(conf.error.is_none());

        let trj =
            compute_md_trajectory_nvt("CCO", &conf.coords, 12, 0.25, 17, 300.0, 10.0).unwrap();
        assert_eq!(trj.frames.len(), 13);
        assert!(trj.frames.iter().all(|f| f.temperature_k.is_finite()));
    }

    #[test]
    fn test_compute_simplified_neb_path_runs() {
        let c1 = embed("CC", 42);
        let c2 = embed("CC", 43);
        assert!(c1.error.is_none());
        assert!(c2.error.is_none());

        let path =
            compute_simplified_neb_path("CC", &c1.coords, &c2.coords, 6, 20, 0.01, 1e-5).unwrap();
        assert_eq!(path.images.len(), 6);
        assert!(path
            .images
            .iter()
            .all(|img| img.potential_energy_kcal_mol.is_finite()));
    }

    #[test]
    fn test_compute_hf3c_water() {
        let conf = embed("O", 42);
        assert!(conf.error.is_none());
        let pos: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

        let result = compute_hf3c(&conf.elements, &pos, &hf::HfConfig::default());
        assert!(result.is_ok(), "HF-3c should succeed for water");
        let r = result.unwrap();
        assert!(r.energy.is_finite());
        assert!(!r.orbital_energies.is_empty());
    }

    #[test]
    fn test_compute_ani_water() {
        let conf = embed("O", 42);
        assert!(conf.error.is_none());
        let pos: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

        let result = compute_ani(&conf.elements, &pos);
        assert!(result.is_ok(), "ANI should succeed for water");
        let r = result.unwrap();
        assert!(r.energy.is_finite());
        assert_eq!(r.forces.len(), 3); // 3 atoms in water
        assert_eq!(r.species.len(), 3);
    }

    #[test]
    fn test_compute_esp_grid_water() {
        let conf = embed("O", 42);
        assert!(conf.error.is_none());
        let pos: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

        let result = compute_esp_grid(&conf.elements, &pos, 0.5, 3.0);
        assert!(result.is_ok(), "ESP grid should succeed for water");
        let g = result.unwrap();
        assert!(!g.values.is_empty());
        assert!(g.spacing > 0.0);
        assert!(g.dims[0] > 0 && g.dims[1] > 0 && g.dims[2] > 0);
    }
}
