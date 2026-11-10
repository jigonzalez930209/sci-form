//! Thermodynamic-cycle electrochemistry: solution pKa and redox potentials.
//!
//! The workflow is explicit: the caller supplies the two molecular states
//! (acid/conjugate base or oxidised/reduced form). Each state is embedded into
//! conformers, evaluated with a closed-shell SCF, and solvated with ALPB. The
//! resulting electronic-plus-solvation energies enter a thermodynamic cycle.
//!
//! # Relative versus absolute values
//!
//! Absolute pKa and redox potentials require a solvation free energy for the
//! proton or electron and inherit the systematic error of the electronic
//! structure method. The recommended workflow is therefore a **relative**
//! cycle against a reference couple with a known experimental value: the
//! systematic error largely cancels. Absolute values are still available, but
//! the result carries the convention and a warning in [`PkaResult::notes`] or
//! [`RedoxResult::notes`].
//!
//! # Units and conventions
//!
//! * energies are in **kcal/mol**;
//! * temperature is in **kelvin**;
//! * pKa is dimensionless;
//! * potentials are in **volt**;
//! * the default proton solvation free energy is **-265.9 kcal/mol**;
//! * the default absolute potential of the standard hydrogen electrode is
//!   **4.44 V**.
//!
//! # Accuracy
//!
//! The default electronic method is HF/STO-3G because it supports charged
//! closed-shell states and the existing GFN2-xTB path is neutral-only. This is
//! a fast first-principles screening model, not a claim of chemical accuracy.
//! The result reports the conformer spread, the method and the exact
//! conventions used.

use serde::{Deserialize, Serialize};

use crate::error::SciFormError;
use crate::solvation_alpb::AlpbConfig;

/// Conversion factor: one electron-volt in kcal/mol.
pub const EV_TO_KCAL_MOL: f64 = 23.060_548_867;

/// Conversion factor: one Hartree in kcal/mol.
pub const HARTREE_TO_KCAL_MOL: f64 = 627.509_474_063_1;

/// Molar gas constant in kcal/(mol·K).
pub const R_KCAL_PER_MOL_K: f64 = 0.001_987_204_258;

/// Faraday constant in kcal/(mol·V).
pub const FARADAY_KCAL_PER_MOL_V: f64 = 23.060_548_867;

/// Default solvation free energy of a proton in water, in kcal/mol.
pub const DEFAULT_PROTON_SOLVATION_KCAL_MOL: f64 = -265.9;

/// Default absolute potential of the standard hydrogen electrode, in volt.
pub const DEFAULT_SHE_ABSOLUTE_POTENTIAL_V: f64 = 4.44;

/// Electronic-structure backend for the electrochemistry workflow.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ElectrochemistryMethod {
    /// Closed-shell UHF/STO-3G. Supports charged states; the default.
    Uhf,
    /// GFN2-xTB. Faster, but only neutral closed-shell states are supported.
    Xtb,
}

/// How the conformer energies are condensed into one state energy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ConformerStatistic {
    /// Use the lowest-energy conformer.
    Minimum,
    /// Use the Boltzmann average at the configured temperature.
    Boltzmann,
}

/// Configuration for the pKa and redox workflows.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ElectrochemistryConfig {
    /// Electronic-structure backend.
    pub method: ElectrochemistryMethod,
    /// Solvent dielectric constant (water is 78.5).
    pub solvent_dielectric: f64,
    /// SASA probe radius in Å.
    pub probe_radius: f64,
    /// Non-polar surface tension in kcal/(mol·Å²).
    pub surface_tension: f64,
    /// Temperature in K.
    pub temperature_k: f64,
    /// Number of embedding attempts per state.
    pub n_conformers: usize,
    /// Base seed for conformer generation.
    pub seed: u64,
    /// Solvation free energy of a proton in the chosen solvent, kcal/mol.
    pub proton_solvation_kcal_mol: f64,
    /// Absolute potential of the reference hydrogen electrode, volt.
    pub she_absolute_potential_v: f64,
    /// Conformer condensation rule.
    pub conformer_statistic: ConformerStatistic,
}

impl Default for ElectrochemistryConfig {
    fn default() -> Self {
        Self {
            method: ElectrochemistryMethod::Uhf,
            solvent_dielectric: 78.5,
            probe_radius: 1.4,
            surface_tension: 0.005,
            temperature_k: 298.15,
            n_conformers: 1,
            seed: 42,
            proton_solvation_kcal_mol: DEFAULT_PROTON_SOLVATION_KCAL_MOL,
            she_absolute_potential_v: DEFAULT_SHE_ABSOLUTE_POTENTIAL_V,
            conformer_statistic: ConformerStatistic::Minimum,
        }
    }
}

impl ElectrochemistryConfig {
    fn alpb(&self) -> AlpbConfig {
        AlpbConfig {
            solvent_dielectric: self.solvent_dielectric,
            probe_radius: self.probe_radius,
            surface_tension: self.surface_tension,
        }
    }

    fn validate(&self) -> Result<(), SciFormError> {
        if !self.solvent_dielectric.is_finite() || self.solvent_dielectric < 1.0 {
            return Err(SciFormError::validation(
                "solvent_dielectric must be finite and >= 1.0",
            ));
        }
        if !self.probe_radius.is_finite() || self.probe_radius < 0.0 {
            return Err(SciFormError::validation(
                "probe_radius must be finite and >= 0",
            ));
        }
        if !self.surface_tension.is_finite() {
            return Err(SciFormError::validation("surface_tension must be finite"));
        }
        if !self.temperature_k.is_finite() || self.temperature_k <= 0.0 {
            return Err(SciFormError::validation(
                "temperature_k must be finite and > 0",
            ));
        }
        if !self.proton_solvation_kcal_mol.is_finite() {
            return Err(SciFormError::validation(
                "proton_solvation_kcal_mol must be finite",
            ));
        }
        if !self.she_absolute_potential_v.is_finite() {
            return Err(SciFormError::validation(
                "she_absolute_potential_v must be finite",
            ));
        }
        if self.n_conformers == 0 {
            return Err(SciFormError::validation("n_conformers must be >= 1"));
        }
        Ok(())
    }
}

/// A reference acid with a known experimental pKa.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PkaReference {
    /// Reference acid SMILES.
    pub acid_smiles: String,
    /// Reference conjugate-base SMILES.
    pub base_smiles: String,
    /// Experimental pKa of the reference acid.
    pub pka: f64,
}

impl PkaReference {
    /// Acetic acid in water, pKa 4.76 at 298.15 K.
    pub fn acetic_acid() -> Self {
        Self {
            acid_smiles: "CC(=O)O".to_string(),
            base_smiles: "CC(=O)[O-]".to_string(),
            pka: 4.76,
        }
    }
}

/// A reference redox couple with a known experimental potential.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RedoxReference {
    /// Oxidised SMILES.
    pub oxidized_smiles: String,
    /// Reduced SMILES.
    pub reduced_smiles: String,
    /// Number of electrons in the reference couple.
    pub n_electrons: u32,
    /// Experimental potential versus the configured reference electrode, volt.
    pub potential_v: f64,
}

/// One conformer of a molecular state with its energy decomposition.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConformerEnergy {
    /// Electronic energy in kcal/mol.
    pub electronic_kcal_mol: f64,
    /// ALPB total solvation energy in kcal/mol.
    pub solvation_kcal_mol: f64,
    /// Sum of the two terms in kcal/mol.
    pub total_kcal_mol: f64,
    /// Mulliken charges of the conformer, in elementary charge.
    pub mulliken_charges: Vec<f64>,
    /// Flat xyz coordinates in Å.
    pub coords: Vec<f64>,
}

/// Energy summary for one molecular state.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StateEnergy {
    /// Input SMILES.
    pub smiles: String,
    /// Total formal charge parsed from the SMILES.
    pub charge: i32,
    /// Element atomic numbers.
    pub elements: Vec<u8>,
    /// Number of successful conformers.
    pub n_conformers: usize,
    /// All successful conformer energies, in kcal/mol.
    pub conformers: Vec<ConformerEnergy>,
    /// Lowest conformer energy in kcal/mol.
    pub minimum_kcal_mol: f64,
    /// Boltzmann-average conformer energy in kcal/mol.
    pub boltzmann_kcal_mol: f64,
    /// Energy selected by [`ConformerStatistic`], in kcal/mol.
    pub selected_kcal_mol: f64,
    /// Spread between the lowest and highest successful conformer, kcal/mol.
    pub spread_kcal_mol: f64,
    /// Mulliken charges of the selected conformer.
    pub selected_charges: Vec<f64>,
}

/// Summary of the reference acid used by a relative pKa.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PkaReferenceResult {
    /// Reference acid SMILES.
    pub acid_smiles: String,
    /// Reference base SMILES.
    pub base_smiles: String,
    /// Reference pKa supplied by the caller.
    pub pka: f64,
    /// Reference acid state energy.
    pub acid: StateEnergy,
    /// Reference base state energy.
    pub base: StateEnergy,
    /// Reference deprotonation free energy in kcal/mol.
    pub delta_g_kcal_mol: f64,
}

/// Result of a pKa thermodynamic cycle.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PkaResult {
    /// Acid SMILES (HA).
    pub acid_smiles: String,
    /// Conjugate-base SMILES (A-).
    pub base_smiles: String,
    /// Acid state energy.
    pub acid: StateEnergy,
    /// Base state energy.
    pub base: StateEnergy,
    /// Target deprotonation free energy in kcal/mol.
    pub delta_g_kcal_mol: f64,
    /// Predicted pKa.
    pub pka: f64,
    /// True when the pKa is relative to a reference acid.
    pub relative: bool,
    /// Reference acid, when the relative workflow was used.
    pub reference: Option<PkaReferenceResult>,
    /// Temperature in K.
    pub temperature_k: f64,
    /// Proton solvation free energy used in the absolute cycle, kcal/mol.
    pub proton_solvation_kcal_mol: f64,
    /// Human-readable notes and caveats.
    pub notes: Vec<String>,
}

/// Summary of the reference couple used by a relative redox potential.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RedoxReferenceResult {
    /// Oxidised SMILES.
    pub oxidized_smiles: String,
    /// Reduced SMILES.
    pub reduced_smiles: String,
    /// Reference potential supplied by the caller, volt.
    pub potential_v: f64,
    /// Reference oxidised state.
    pub oxidized: StateEnergy,
    /// Reference reduced state.
    pub reduced: StateEnergy,
    /// Reference reaction free energy in kcal/mol.
    pub delta_g_kcal_mol: f64,
}

/// Result of a redox thermodynamic cycle.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RedoxResult {
    /// Oxidised SMILES.
    pub oxidized_smiles: String,
    /// Reduced SMILES.
    pub reduced_smiles: String,
    /// Oxidised state energy.
    pub oxidized: StateEnergy,
    /// Reduced state energy.
    pub reduced: StateEnergy,
    /// Number of electrons transferred.
    pub n_electrons: u32,
    /// Reaction free energy in kcal/mol.
    pub delta_g_kcal_mol: f64,
    /// Absolute electrode potential in volt.
    pub absolute_potential_v: f64,
    /// Potential versus the configured reference electrode in volt.
    pub potential_vs_reference_v: f64,
    /// True when the potential is relative to a reference couple.
    pub relative: bool,
    /// Reference couple, when the relative workflow was used.
    pub reference: Option<RedoxReferenceResult>,
    /// Absolute potential of the reference electrode used, volt.
    pub reference_absolute_potential_v: f64,
    /// Human-readable notes and caveats.
    pub notes: Vec<String>,
}

fn total_formal_charge(smiles: &str) -> Result<i32, SciFormError> {
    let molecule = crate::parse(smiles)?;
    Ok(molecule
        .graph
        .node_weights()
        .map(|atom| atom.formal_charge as i32)
        .sum())
}

fn electron_count(elements: &[u8], charge: i32) -> Result<usize, SciFormError> {
    let neutral: i32 = elements
        .iter()
        .map(|&z| crate::scf::basis::valence_electrons(z) as i32)
        .sum();
    let count = neutral - charge;
    if count <= 0 {
        return Err(SciFormError::validation(
            "charge leaves a non-positive electron count",
        ));
    }
    Ok(count as usize)
}

/// Evaluate one molecular state with the configured electronic method and ALPB.
///
/// `smiles` is parsed and embedded with the existing ETKDG pipeline. The total
/// formal charge comes from the parsed SMILES, so `[O-]` and `[NH4+]` are
/// handled without a separate charge argument.
pub fn evaluate_state(
    smiles: &str,
    config: &ElectrochemistryConfig,
) -> Result<StateEnergy, SciFormError> {
    config.validate()?;
    let charge = total_formal_charge(smiles)?;
    if config.method == ElectrochemistryMethod::Xtb && charge != 0 {
        return Err(SciFormError::unsupported(
            "GFN2-xTB backend currently supports neutral closed-shell states only; use the UHF backend for charged states",
        ));
    }

    let conformers = crate::embed_diverse(smiles, config.n_conformers, 1.0, config.seed);
    let first = conformers
        .iter()
        .find(|conformer| conformer.error.is_none() && !conformer.coords.is_empty())
        .ok_or_else(|| {
            SciFormError::numerical(format!(
                "conformer generation failed for {smiles}: {}",
                conformers
                    .first()
                    .and_then(|conformer| conformer.error.clone())
                    .unwrap_or_else(|| "no successful conformer".to_string())
            ))
        })?;
    let elements = first.elements.clone();

    let mut evaluated = Vec::new();
    for conformer in &conformers {
        if conformer.error.is_some() || conformer.coords.len() != elements.len() * 3 {
            continue;
        }
        let positions = flat_to_positions(&conformer.coords);
        let (electronic_kcal, mulliken) = match config.method {
            ElectrochemistryMethod::Uhf => {
                let n_electrons = electron_count(&conformer.elements, charge)?;
                if n_electrons % 2 != 0 {
                    return Err(SciFormError::unsupported(
                        "the electrochemistry workflow currently requires a closed-shell state with an even electron count",
                    ));
                }
                let scf = crate::compute_uhf(&conformer.elements, &positions, charge, 1)?;
                if !scf.converged {
                    return Err(SciFormError::numerical(format!(
                        "UHF did not converge for {smiles}"
                    )));
                }
                (scf.total_energy * HARTREE_TO_KCAL_MOL, scf.mulliken_charges)
            }
            ElectrochemistryMethod::Xtb => {
                let xtb = crate::compute_xtb(&conformer.elements, &positions)?;
                (xtb.total_energy * EV_TO_KCAL_MOL, xtb.mulliken_charges)
            }
        };
        let alpb = crate::solvation_alpb::compute_alpb_solvation(
            &conformer.elements,
            &positions,
            &mulliken,
            &config.alpb(),
        )
        .map_err(SciFormError::from)?;
        evaluated.push(ConformerEnergy {
            electronic_kcal_mol: electronic_kcal,
            solvation_kcal_mol: alpb.total_energy,
            total_kcal_mol: electronic_kcal + alpb.total_energy,
            mulliken_charges: mulliken,
            coords: conformer.coords.clone(),
        });
    }

    if evaluated.is_empty() {
        return Err(SciFormError::numerical(format!(
            "no conformer of {smiles} could be evaluated"
        )));
    }

    let minimum = evaluated
        .iter()
        .map(|conformer| conformer.total_kcal_mol)
        .fold(f64::INFINITY, f64::min);
    let maximum = evaluated
        .iter()
        .map(|conformer| conformer.total_kcal_mol)
        .fold(f64::NEG_INFINITY, f64::max);
    let boltzmann = boltzmann_average(&evaluated, config.temperature_k);
    let selected = match config.conformer_statistic {
        ConformerStatistic::Minimum => minimum,
        ConformerStatistic::Boltzmann => boltzmann,
    };
    let selected_index = evaluated
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            a.total_kcal_mol
                .partial_cmp(&b.total_kcal_mol)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(index, _)| index)
        .expect("evaluated is non-empty");
    let selected_charges = evaluated[selected_index].mulliken_charges.clone();

    Ok(StateEnergy {
        smiles: smiles.to_string(),
        charge,
        elements,
        n_conformers: evaluated.len(),
        conformers: evaluated,
        minimum_kcal_mol: minimum,
        boltzmann_kcal_mol: boltzmann,
        selected_kcal_mol: selected,
        spread_kcal_mol: maximum - minimum,
        selected_charges,
    })
}

fn deprotonation_energy(
    acid: &StateEnergy,
    base: &StateEnergy,
    proton_solvation_kcal_mol: f64,
) -> f64 {
    base.selected_kcal_mol + proton_solvation_kcal_mol - acid.selected_kcal_mol
}

/// Compute a solution pKa from an explicit acid/conjugate-base pair.
///
/// With `reference = Some(...)`, the result is a relative pKa:
///
/// ```text
/// pKa(target) = pKa(ref) + (ΔG(target) - ΔG(ref)) / (RT ln 10)
/// ```
///
/// The reference cancels the proton-solvation convention and most of the
/// method's systematic error. With `reference = None`, the absolute cycle uses
/// the configured proton solvation free energy; that value is reported and
/// flagged in the notes.
///
/// The caller supplies both SMILES because automatic protonation-site
/// enumeration is a separate problem. `base_smiles` must be the conjugate base.
pub fn compute_pka(
    acid_smiles: &str,
    base_smiles: &str,
    reference: Option<&PkaReference>,
    config: &ElectrochemistryConfig,
) -> Result<PkaResult, SciFormError> {
    config.validate()?;
    let acid = evaluate_state(acid_smiles, config)?;
    let base = evaluate_state(base_smiles, config)?;
    let delta_g = deprotonation_energy(&acid, &base, config.proton_solvation_kcal_mol);
    if !delta_g.is_finite() {
        return Err(SciFormError::numerical(
            "deprotonation free energy is not finite",
        ));
    }

    let (pka, relative, reference_result) = if let Some(reference) = reference {
        let ref_acid = evaluate_state(&reference.acid_smiles, config)?;
        let ref_base = evaluate_state(&reference.base_smiles, config)?;
        let ref_delta_g =
            deprotonation_energy(&ref_acid, &ref_base, config.proton_solvation_kcal_mol);
        let pka = reference.pka
            + (delta_g - ref_delta_g)
                / (R_KCAL_PER_MOL_K * config.temperature_k * std::f64::consts::LN_10);
        (
            pka,
            true,
            Some(PkaReferenceResult {
                acid_smiles: reference.acid_smiles.clone(),
                base_smiles: reference.base_smiles.clone(),
                pka: reference.pka,
                acid: ref_acid,
                base: ref_base,
                delta_g_kcal_mol: ref_delta_g,
            }),
        )
    } else {
        (
            delta_g / (R_KCAL_PER_MOL_K * config.temperature_k * std::f64::consts::LN_10),
            false,
            None,
        )
    };
    if !pka.is_finite() {
        return Err(SciFormError::numerical("pKa is not finite"));
    }

    let mut notes = vec![
        "pKa uses an explicit acid/conjugate-base pair supplied by the caller.".to_string(),
        "Electronic energy plus ALPB solvation; thermal and vibrational corrections are not included."
            .to_string(),
        format!("Electronic method: {:?}.", config.method),
    ];
    if relative {
        notes.push(
            "Relative pKa: the reference cancels the proton-solvation convention and most systematic error."
                .to_string(),
        );
    } else {
        notes.push(format!(
            "Absolute pKa: proton solvation fixed at {:.3} kcal/mol; treat the absolute value as a screening estimate.",
            config.proton_solvation_kcal_mol
        ));
    }
    if acid.spread_kcal_mol > 2.0 || base.spread_kcal_mol > 2.0 {
        notes.push(
            "The conformer spread exceeds 2 kcal/mol; inspect the ensemble before using the absolute pKa."
                .to_string(),
        );
    }

    Ok(PkaResult {
        acid_smiles: acid_smiles.to_string(),
        base_smiles: base_smiles.to_string(),
        acid,
        base,
        delta_g_kcal_mol: delta_g,
        pka,
        relative,
        reference: reference_result,
        temperature_k: config.temperature_k,
        proton_solvation_kcal_mol: config.proton_solvation_kcal_mol,
        notes,
    })
}

/// Convert a reaction free energy into absolute and reference-relative potentials.
///
/// `delta_g_kcal_mol` is the free energy of `Ox + n e- -> Red`. The electron
/// free energy is zero by convention; `reference_absolute_potential_v` is the
/// absolute potential of the reference electrode (4.44 V for SHE).
pub fn redox_potentials_from_delta_g(
    delta_g_kcal_mol: f64,
    n_electrons: u32,
    reference_absolute_potential_v: f64,
) -> Result<(f64, f64), SciFormError> {
    if !delta_g_kcal_mol.is_finite() {
        return Err(SciFormError::validation("delta_g_kcal_mol must be finite"));
    }
    if n_electrons == 0 || n_electrons > 4 {
        return Err(SciFormError::validation(
            "n_electrons must be between 1 and 4",
        ));
    }
    if !reference_absolute_potential_v.is_finite() {
        return Err(SciFormError::validation(
            "reference_absolute_potential_v must be finite",
        ));
    }
    let absolute = -delta_g_kcal_mol / (n_electrons as f64 * FARADAY_KCAL_PER_MOL_V);
    Ok((absolute, absolute - reference_absolute_potential_v))
}

/// Compute a redox potential from an explicit oxidised/reduced pair.
///
/// With `reference = Some(...)`, the potential is relative to a known couple:
///
/// ```text
/// E(target) = E(ref) - (ΔG(target) - ΔG(ref)) / (n F)
/// ```
///
/// Without a reference, the absolute cycle subtracts the configured reference
/// electrode potential (SHE at 4.44 V by default).
pub fn compute_redox_potential(
    oxidized_smiles: &str,
    reduced_smiles: &str,
    n_electrons: u32,
    reference: Option<&RedoxReference>,
    config: &ElectrochemistryConfig,
) -> Result<RedoxResult, SciFormError> {
    config.validate()?;
    if n_electrons == 0 || n_electrons > 4 {
        return Err(SciFormError::validation(
            "n_electrons must be between 1 and 4",
        ));
    }
    if let Some(reference) = reference {
        if reference.n_electrons != n_electrons {
            return Err(SciFormError::validation(
                "the reference couple must transfer the same number of electrons",
            ));
        }
    }

    let oxidized = evaluate_state(oxidized_smiles, config)?;
    let reduced = evaluate_state(reduced_smiles, config)?;
    let delta_g = reduced.selected_kcal_mol - oxidized.selected_kcal_mol;
    if !delta_g.is_finite() {
        return Err(SciFormError::numerical("redox free energy is not finite"));
    }
    let (absolute, absolute_vs_reference) =
        redox_potentials_from_delta_g(delta_g, n_electrons, config.she_absolute_potential_v)?;

    let (potential, relative, reference_result) = if let Some(reference) = reference {
        let ref_oxidized = evaluate_state(&reference.oxidized_smiles, config)?;
        let ref_reduced = evaluate_state(&reference.reduced_smiles, config)?;
        let ref_delta_g = ref_reduced.selected_kcal_mol - ref_oxidized.selected_kcal_mol;
        let potential = reference.potential_v
            - (delta_g - ref_delta_g) / (n_electrons as f64 * FARADAY_KCAL_PER_MOL_V);
        (
            potential,
            true,
            Some(RedoxReferenceResult {
                oxidized_smiles: reference.oxidized_smiles.clone(),
                reduced_smiles: reference.reduced_smiles.clone(),
                potential_v: reference.potential_v,
                oxidized: ref_oxidized,
                reduced: ref_reduced,
                delta_g_kcal_mol: ref_delta_g,
            }),
        )
    } else {
        (absolute_vs_reference, false, None)
    };
    if !potential.is_finite() {
        return Err(SciFormError::numerical("redox potential is not finite"));
    }

    let mut notes = vec![
        "Redox potential uses an explicit oxidised/reduced pair supplied by the caller."
            .to_string(),
        "Electron free energy is zero by convention.".to_string(),
        "Electronic energy plus ALPB solvation; thermal and vibrational corrections are not included."
            .to_string(),
        format!("Electronic method: {:?}.", config.method),
    ];
    if relative {
        notes.push(
            "Relative redox potential: the reference cancels the electron convention and most systematic error."
                .to_string(),
        );
    } else {
        notes.push(format!(
            "Absolute redox potential: reference electrode fixed at {:.3} V; treat it as a screening estimate.",
            config.she_absolute_potential_v
        ));
    }
    if oxidized.spread_kcal_mol > 2.0 || reduced.spread_kcal_mol > 2.0 {
        notes.push(
            "The conformer spread exceeds 2 kcal/mol; inspect the ensemble before using the absolute potential."
                .to_string(),
        );
    }

    Ok(RedoxResult {
        oxidized_smiles: oxidized_smiles.to_string(),
        reduced_smiles: reduced_smiles.to_string(),
        oxidized,
        reduced,
        n_electrons,
        delta_g_kcal_mol: delta_g,
        absolute_potential_v: absolute,
        potential_vs_reference_v: potential,
        relative,
        reference: reference_result,
        reference_absolute_potential_v: config.she_absolute_potential_v,
        notes,
    })
}

fn flat_to_positions(coords: &[f64]) -> Vec<[f64; 3]> {
    let (chunks, _) = coords.as_chunks::<3>();
    chunks
        .iter()
        .map(|chunk| [chunk[0], chunk[1], chunk[2]])
        .collect()
}

fn boltzmann_average(conformers: &[ConformerEnergy], temperature_k: f64) -> f64 {
    let minimum = conformers
        .iter()
        .map(|conformer| conformer.total_kcal_mol)
        .fold(f64::INFINITY, f64::min);
    let beta = 1.0 / (R_KCAL_PER_MOL_K * temperature_k);
    let mut weights = Vec::with_capacity(conformers.len());
    let mut partition = 0.0;
    for conformer in conformers {
        let weight = (-(conformer.total_kcal_mol - minimum) * beta).exp();
        partition += weight;
        weights.push(weight);
    }
    if partition <= 0.0 || !partition.is_finite() {
        return minimum;
    }
    conformers
        .iter()
        .zip(weights)
        .map(|(conformer, weight)| conformer.total_kcal_mol * weight)
        .sum::<f64>()
        / partition
}

/// Convenience helper for the pKa formula, exposed for tests and downstream
/// workflows that already have the deprotonation free energy.
pub fn pka_from_delta_g(delta_g_kcal_mol: f64, temperature_k: f64) -> Result<f64, SciFormError> {
    if !delta_g_kcal_mol.is_finite() {
        return Err(SciFormError::validation("delta_g_kcal_mol must be finite"));
    }
    if !temperature_k.is_finite() || temperature_k <= 0.0 {
        return Err(SciFormError::validation("temperature_k must be > 0"));
    }
    Ok(delta_g_kcal_mol / (R_KCAL_PER_MOL_K * temperature_k * std::f64::consts::LN_10))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fast_config() -> ElectrochemistryConfig {
        ElectrochemistryConfig {
            n_conformers: 1,
            ..ElectrochemistryConfig::default()
        }
    }

    #[test]
    fn pka_formula_matches_acetic_acid_scale() {
        let pka = pka_from_delta_g(6.49, 298.15).unwrap();
        assert!((pka - 4.76).abs() < 0.05, "pka = {pka}");
    }

    #[test]
    fn boltzmann_average_is_between_minimum_and_maximum() {
        let conformers = vec![
            ConformerEnergy {
                electronic_kcal_mol: 0.0,
                solvation_kcal_mol: 0.0,
                total_kcal_mol: 0.0,
                mulliken_charges: vec![],
                coords: vec![],
            },
            ConformerEnergy {
                electronic_kcal_mol: 1.0,
                solvation_kcal_mol: 0.0,
                total_kcal_mol: 1.0,
                mulliken_charges: vec![],
                coords: vec![],
            },
        ];
        let average = boltzmann_average(&conformers, 298.15);
        assert!(average > 0.0 && average < 1.0);
    }

    #[test]
    fn config_rejects_bad_temperature() {
        let mut config = ElectrochemistryConfig::default();
        config.temperature_k = 0.0;
        assert!(config.validate().is_err());
    }

    #[test]
    fn pka_rejects_non_finite_energy() {
        assert!(pka_from_delta_g(f64::NAN, 298.15).is_err());
    }

    #[test]
    fn formal_charge_is_read_from_smiles() {
        assert_eq!(total_formal_charge("CC(=O)[O-]").unwrap(), -1);
        assert_eq!(total_formal_charge("[NH4+]").unwrap(), 1);
        assert_eq!(total_formal_charge("CCO").unwrap(), 0);
    }

    #[test]
    fn reference_acid_reproduces_its_known_pka() {
        let reference = PkaReference::acetic_acid();
        let result = compute_pka(
            &reference.acid_smiles,
            &reference.base_smiles,
            Some(&reference),
            &fast_config(),
        )
        .unwrap();
        assert!(result.relative);
        assert!((result.pka - reference.pka).abs() < 1e-9);
    }

    #[test]
    fn redox_potential_conversion_is_consistent() {
        // A one-electron reduction with ΔG = -23.06 kcal/mol is 1 V absolute.
        let (absolute, vs_she) =
            redox_potentials_from_delta_g(-FARADAY_KCAL_PER_MOL_V, 1, 4.44).unwrap();
        assert!((absolute - 1.0).abs() < 1e-12);
        assert!((vs_she + 3.44).abs() < 1e-12);
    }

    #[test]
    #[ignore = "heavy: UHF/STO-3G on six states"]
    fn trifluoroacetic_acid_is_more_acidic_than_acetic_acid() {
        let reference = PkaReference::acetic_acid();
        let result = compute_pka(
            "FC(F)(F)C(=O)O",
            "FC(F)(F)C(=O)[O-]",
            Some(&reference),
            &fast_config(),
        )
        .unwrap();
        eprintln!("TFA relative pKa = {:.3}", result.pka);
        assert!(result.pka.is_finite());
        assert!(result.pka < reference.pka);
    }

    #[test]
    fn xtb_rejects_charged_state() {
        let config = ElectrochemistryConfig {
            method: ElectrochemistryMethod::Xtb,
            ..fast_config()
        };
        let error = evaluate_state("CC(=O)[O-]", &config).unwrap_err();
        assert!(error.is_unsupported());
    }
}
