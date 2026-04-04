use super::saddle::{find_transition_state, GsmResult};
use super::string::GsmConfig;
use crate::eht::{analyze_eht_support, SupportLevel};
use crate::graph::Molecule;
use crate::hf::{solve_hf3c, HfConfig};
use crate::scf::constants::{EV_TO_HARTREE, HARTREE_TO_KCAL};
use crate::xtb::{gfn1::solve_gfn1, gfn2::solve_gfn2, is_xtb_supported};
use serde::{Deserialize, Serialize};
use std::cell::RefCell;
use std::str::FromStr;

const EV_TO_KCAL: f64 = EV_TO_HARTREE * HARTREE_TO_KCAL;
const HF_EXPLICIT_BASIS_ELEMENTS: [u8; 8] = [1, 6, 7, 8, 9, 15, 16, 17];
const REAXFF_CHON_ELEMENTS: [u8; 4] = [1, 6, 7, 8];

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GsmEnergyBackend {
    Uff,
    Mmff94,
    Eht,
    Pm3,
    Xtb,
    Gfn1,
    Gfn2,
    Hf3c,
    Reaxff,
}

impl GsmEnergyBackend {
    pub const ALL: [Self; 9] = [
        Self::Uff,
        Self::Mmff94,
        Self::Eht,
        Self::Pm3,
        Self::Xtb,
        Self::Gfn1,
        Self::Gfn2,
        Self::Hf3c,
        Self::Reaxff,
    ];

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Uff => "uff",
            Self::Mmff94 => "mmff94",
            Self::Eht => "eht",
            Self::Pm3 => "pm3",
            Self::Xtb => "xtb",
            Self::Gfn1 => "gfn1",
            Self::Gfn2 => "gfn2",
            Self::Hf3c => "hf3c",
            Self::Reaxff => "reaxff",
        }
    }
}

impl FromStr for GsmEnergyBackend {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.trim().to_ascii_lowercase().as_str() {
            "uff" => Ok(Self::Uff),
            "mmff94" | "mmff" => Ok(Self::Mmff94),
            "eht" | "he" => Ok(Self::Eht),
            "pm3" => Ok(Self::Pm3),
            "xtb" | "gfn0" | "gfn0_xtb" | "gfn0-xtb" => Ok(Self::Xtb),
            "gfn1" | "gfn1_xtb" | "gfn1-xtb" => Ok(Self::Gfn1),
            "gfn2" | "gfn2_xtb" | "gfn2-xtb" => Ok(Self::Gfn2),
            "hf3c" | "hf-3c" => Ok(Self::Hf3c),
            "reaxff" => Ok(Self::Reaxff),
            other => Err(format!(
                "Unknown GSM backend '{}'. Expected one of: uff, mmff94, eht, pm3, xtb, gfn1, gfn2, hf3c, reaxff",
                other
            )),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GsmBackendStatus {
    Success,
    Unavailable,
    Error,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GsmBackendCapability {
    pub backend: GsmEnergyBackend,
    pub available: bool,
    pub suitable_for_reaction_path: bool,
    pub model_family: String,
    pub energy_unit: String,
    pub support_level: Option<SupportLevel>,
    pub warnings: Vec<String>,
    pub limitations: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GsmBackendEvaluation {
    pub backend: GsmEnergyBackend,
    pub status: GsmBackendStatus,
    pub available: bool,
    pub suitable_for_reaction_path: bool,
    pub energy_kcal_mol: Option<f64>,
    pub raw_energy: Option<f64>,
    pub raw_energy_unit: Option<String>,
    pub converged: Option<bool>,
    pub warnings: Vec<String>,
    pub limitations: Vec<String>,
    pub error: Option<String>,
}

pub fn plan_gsm_backends(elements: &[u8]) -> Vec<GsmBackendCapability> {
    GsmEnergyBackend::ALL
        .into_iter()
        .map(|backend| capability_for_backend(backend, elements))
        .collect()
}

pub fn plan_gsm_backends_for_smiles(smiles: &str) -> Result<Vec<GsmBackendCapability>, String> {
    let elements = elements_from_smiles(smiles)?;
    Ok(plan_gsm_backends(&elements))
}

pub fn evaluate_gsm_backend(
    smiles: &str,
    coords: &[f64],
    backend: GsmEnergyBackend,
) -> Result<GsmBackendEvaluation, String> {
    let elements = elements_from_smiles(smiles)?;
    evaluate_gsm_backend_with_elements(smiles, &elements, coords, backend)
}

pub fn compare_gsm_backends(
    smiles: &str,
    coords: &[f64],
    backends: &[GsmEnergyBackend],
) -> Result<Vec<GsmBackendEvaluation>, String> {
    let elements = elements_from_smiles(smiles)?;
    let selected = if backends.is_empty() {
        GsmEnergyBackend::ALL.to_vec()
    } else {
        backends.to_vec()
    };
    selected
        .into_iter()
        .map(|backend| evaluate_gsm_backend_with_elements(smiles, &elements, coords, backend))
        .collect()
}

#[allow(dead_code)]
pub(crate) fn evaluate_gsm_backend_energy_kcal(
    smiles: &str,
    elements: &[u8],
    coords: &[f64],
    backend: GsmEnergyBackend,
) -> Result<f64, String> {
    backend_energy_kcal_mol(backend, smiles, elements, coords)
}

#[allow(dead_code)]
pub(crate) fn rings_from_smiles(smiles: &str) -> Result<Vec<(Vec<usize>, bool)>, String> {
    let sssr = crate::compute_sssr(smiles)?;
    Ok(sssr
        .rings
        .into_iter()
        .map(|ring| (ring.atoms, ring.is_aromatic))
        .collect())
}

#[allow(dead_code)]
pub(crate) fn elements_for_smiles(smiles: &str) -> Result<Vec<u8>, String> {
    elements_from_smiles(smiles)
}

pub fn find_transition_state_with_backend(
    smiles: &str,
    reactant: &[f64],
    product: &[f64],
    backend: GsmEnergyBackend,
    config: &GsmConfig,
) -> Result<GsmResult, String> {
    let elements = elements_from_smiles(smiles)?;
    validate_flat_coords(reactant, elements.len())?;
    validate_flat_coords(product, elements.len())?;
    if reactant.len() != product.len() {
        return Err("reactant/product coordinate lengths must match".to_string());
    }

    let capability = capability_for_backend(backend, &elements);
    if !capability.available || !capability.suitable_for_reaction_path {
        return Err(capability_message(&capability));
    }

    let owned_smiles = smiles.to_string();
    let owned_elements = elements.clone();
    let first_error = RefCell::new(None::<String>);
    let energy_fn = |coords: &[f64]| -> f64 {
        match backend_energy_kcal_mol(backend, &owned_smiles, &owned_elements, coords) {
            Ok(energy) => energy,
            Err(err) => {
                let mut slot = first_error.borrow_mut();
                if slot.is_none() {
                    *slot = Some(err);
                }
                f64::INFINITY
            }
        }
    };

    let result = find_transition_state(reactant, product, &energy_fn, config);
    if let Some(err) = first_error.into_inner() {
        return Err(err);
    }
    if !result.ts_energy.is_finite() || result.path_energies.iter().any(|e| !e.is_finite()) {
        return Err(format!(
            "GSM backend '{}' produced a non-finite reaction path",
            backend.as_str()
        ));
    }
    Ok(result)
}

fn evaluate_gsm_backend_with_elements(
    smiles: &str,
    elements: &[u8],
    coords: &[f64],
    backend: GsmEnergyBackend,
) -> Result<GsmBackendEvaluation, String> {
    validate_flat_coords(coords, elements.len())?;
    let capability = capability_for_backend(backend, elements);
    if !capability.available || !capability.suitable_for_reaction_path {
        let error_message = capability_message(&capability);
        return Ok(GsmBackendEvaluation {
            backend,
            status: GsmBackendStatus::Unavailable,
            available: capability.available,
            suitable_for_reaction_path: capability.suitable_for_reaction_path,
            energy_kcal_mol: None,
            raw_energy: None,
            raw_energy_unit: None,
            converged: None,
            warnings: capability.warnings,
            limitations: capability.limitations,
            error: Some(error_message),
        });
    }

    let evaluation = match backend_raw_energy(backend, smiles, elements, coords) {
        Ok((raw_energy, raw_energy_unit, energy_kcal_mol, converged)) => GsmBackendEvaluation {
            backend,
            status: GsmBackendStatus::Success,
            available: capability.available,
            suitable_for_reaction_path: capability.suitable_for_reaction_path,
            energy_kcal_mol: Some(energy_kcal_mol),
            raw_energy: Some(raw_energy),
            raw_energy_unit: Some(raw_energy_unit.to_string()),
            converged,
            warnings: capability.warnings,
            limitations: capability.limitations,
            error: None,
        },
        Err(err) => GsmBackendEvaluation {
            backend,
            status: GsmBackendStatus::Error,
            available: capability.available,
            suitable_for_reaction_path: capability.suitable_for_reaction_path,
            energy_kcal_mol: None,
            raw_energy: None,
            raw_energy_unit: None,
            converged: None,
            warnings: capability.warnings,
            limitations: capability.limitations,
            error: Some(err),
        },
    };

    Ok(evaluation)
}

fn capability_for_backend(backend: GsmEnergyBackend, elements: &[u8]) -> GsmBackendCapability {
    match backend {
        GsmEnergyBackend::Uff => GsmBackendCapability {
            backend,
            available: true,
            suitable_for_reaction_path: true,
            model_family: "force_field".to_string(),
            energy_unit: "kcal/mol".to_string(),
            support_level: None,
            warnings: vec![
                "UFF is broadly covered but still topology- and typing-dependent.".to_string(),
            ],
            limitations: vec![
                "Useful for exploratory paths, not a quantitative electronic barrier method."
                    .to_string(),
            ],
        },
        GsmEnergyBackend::Mmff94 => GsmBackendCapability {
            backend,
            available: true,
            suitable_for_reaction_path: true,
            model_family: "force_field".to_string(),
            energy_unit: "kcal/mol".to_string(),
            support_level: None,
            warnings: vec![
                "MMFF94 is most reliable for organic-like valence states and conventional atom typing."
                    .to_string(),
            ],
            limitations: vec![
                "Reactive bond breaking/forming can exceed the intended MMFF94 domain.".to_string(),
            ],
        },
        GsmEnergyBackend::Eht => {
            let support = analyze_eht_support(elements);
            let mut limitations = vec![
                "EHT exposes orbital trends, not a calibrated total-energy surface for barrier heights."
                    .to_string(),
                "Use PM3, xTB/GFN, HF-3c, or a force field backend for GSM path searches."
                    .to_string(),
            ];
            if !support.unsupported_elements.is_empty() {
                limitations.push(format!(
                    "Unsupported EHT elements: {:?}",
                    support.unsupported_elements
                ));
            }
            GsmBackendCapability {
                backend,
                available: support.level != SupportLevel::Unsupported,
                suitable_for_reaction_path: false,
                model_family: "orbital_method".to_string(),
                energy_unit: "eV".to_string(),
                support_level: Some(support.level),
                warnings: support.warnings,
                limitations,
            }
        }
        GsmEnergyBackend::Pm3 => {
            let unsupported = unsupported_elements(elements, crate::pm3::is_pm3_supported);
            GsmBackendCapability {
                backend,
                available: unsupported.is_empty(),
                suitable_for_reaction_path: unsupported.is_empty(),
                model_family: "semiempirical_scf".to_string(),
                energy_unit: "eV".to_string(),
                support_level: None,
                warnings: Vec::new(),
                limitations: unsupported_message("PM3", &unsupported),
            }
        }
        GsmEnergyBackend::Xtb | GsmEnergyBackend::Gfn1 | GsmEnergyBackend::Gfn2 => {
            let unsupported = unsupported_elements(elements, is_xtb_supported);
            let label = backend.as_str().to_ascii_uppercase();
            GsmBackendCapability {
                backend,
                available: unsupported.is_empty(),
                suitable_for_reaction_path: unsupported.is_empty(),
                model_family: "tight_binding".to_string(),
                energy_unit: "eV".to_string(),
                support_level: None,
                warnings: Vec::new(),
                limitations: unsupported_message(&label, &unsupported),
            }
        }
        GsmEnergyBackend::Hf3c => {
            let fallback: Vec<u8> = elements
                .iter()
                .copied()
                .filter(|z| !HF_EXPLICIT_BASIS_ELEMENTS.contains(z))
                .collect();
            let warnings = if fallback.is_empty() {
                Vec::new()
            } else {
                vec![format!(
                    "HF-3c will use the current hydrogen-like STO-3G fallback for elements {:?}.",
                    fallback
                )]
            };
            GsmBackendCapability {
                backend,
                available: true,
                suitable_for_reaction_path: true,
                model_family: "ab_initio_composite".to_string(),
                energy_unit: "hartree".to_string(),
                support_level: None,
                warnings,
                limitations: vec![
                    "HF-3c is much more expensive inside GSM because the path search uses many finite-difference single points."
                        .to_string(),
                ],
            }
        }
        GsmEnergyBackend::Reaxff => {
            let unsupported = unsupported_elements(elements, |z| REAXFF_CHON_ELEMENTS.contains(&z));
            let limitations = unsupported_message("ReaxFF default CHON parameter set", &unsupported);
            #[cfg(not(feature = "alpha-reaxff"))]
            {
                limitations.push(
                    "This build does not enable the 'alpha-reaxff' feature flag.".to_string(),
                );
            }
            GsmBackendCapability {
                backend,
                available: unsupported.is_empty() && cfg!(feature = "alpha-reaxff"),
                suitable_for_reaction_path: unsupported.is_empty() && cfg!(feature = "alpha-reaxff"),
                model_family: "reactive_force_field".to_string(),
                energy_unit: "kcal/mol".to_string(),
                support_level: None,
                warnings: vec![
                    "Current alpha ReaxFF exposure uses the default CHON parameter set.".to_string(),
                ],
                limitations,
            }
        }
    }
}

fn backend_energy_kcal_mol(
    backend: GsmEnergyBackend,
    smiles: &str,
    elements: &[u8],
    coords: &[f64],
) -> Result<f64, String> {
    let (_, _, energy_kcal_mol, _) = backend_raw_energy(backend, smiles, elements, coords)?;
    Ok(energy_kcal_mol)
}

fn backend_raw_energy(
    backend: GsmEnergyBackend,
    smiles: &str,
    elements: &[u8],
    coords: &[f64],
) -> Result<(f64, &'static str, f64, Option<bool>), String> {
    validate_flat_coords(coords, elements.len())?;
    let positions = flat_to_positions(coords);
    match backend {
        GsmEnergyBackend::Uff => {
            let energy = crate::compute_uff_energy(smiles, coords)?;
            Ok((energy, "kcal/mol", energy, None))
        }
        GsmEnergyBackend::Mmff94 => {
            let energy = crate::compute_mmff94_energy(smiles, coords)?;
            Ok((energy, "kcal/mol", energy, None))
        }
        GsmEnergyBackend::Eht => Err(
            "EHT is not exposed as a GSM backend because it does not provide a barrier-ready total-energy surface."
                .to_string(),
        ),
        GsmEnergyBackend::Pm3 => {
            let result = crate::compute_pm3(elements, &positions)?;
            if !result.converged {
                return Err("PM3 did not converge at this geometry".to_string());
            }
            Ok((
                result.total_energy,
                "eV",
                result.total_energy * EV_TO_KCAL,
                Some(result.converged),
            ))
        }
        GsmEnergyBackend::Xtb => {
            let result = crate::compute_xtb(elements, &positions)?;
            if !result.converged {
                return Err("GFN0-xTB did not converge at this geometry".to_string());
            }
            Ok((
                result.total_energy,
                "eV",
                result.total_energy * EV_TO_KCAL,
                Some(result.converged),
            ))
        }
        GsmEnergyBackend::Gfn1 => {
            let result = solve_gfn1(elements, &positions)?;
            if !result.converged {
                return Err("GFN1-xTB did not converge at this geometry".to_string());
            }
            Ok((
                result.total_energy,
                "eV",
                result.total_energy * EV_TO_KCAL,
                Some(result.converged),
            ))
        }
        GsmEnergyBackend::Gfn2 => {
            let result = solve_gfn2(elements, &positions)?;
            if !result.converged {
                return Err("GFN2-xTB did not converge at this geometry".to_string());
            }
            Ok((
                result.total_energy,
                "eV",
                result.total_energy * EV_TO_KCAL,
                Some(result.converged),
            ))
        }
        GsmEnergyBackend::Hf3c => {
            let result = solve_hf3c(elements, &positions, &HfConfig::default())?;
            if !result.converged {
                return Err("HF-3c SCF did not converge at this geometry".to_string());
            }
            Ok((
                result.energy,
                "hartree",
                result.energy * HARTREE_TO_KCAL,
                Some(result.converged),
            ))
        }
        GsmEnergyBackend::Reaxff => {
            #[cfg(feature = "alpha-reaxff")]
            {
                let params = crate::forcefield::reaxff::params::ReaxffParams::default_chon();
                let atom_params = build_reaxff_atom_params(elements, &params);
                let eem_params: Vec<_> = elements
                    .iter()
                    .map(|&z| crate::forcefield::reaxff::eem::default_eem_params(z))
                    .collect();
                let charges = crate::forcefield::reaxff::eem::solve_eem(coords, &eem_params, 0.0)?;
                let bond_orders = crate::forcefield::reaxff::bond_order::compute_bond_orders(
                    coords,
                    &atom_params,
                    params.cutoff,
                );
                let bonded = crate::forcefield::reaxff::energy::compute_bonded_energy(
                    coords,
                    &bond_orders,
                    &params,
                );
                let (van_der_waals, coulomb) =
                    crate::forcefield::reaxff::nonbonded::compute_nonbonded_energy(
                        coords,
                        &charges,
                        elements,
                        params.cutoff,
                    );
                let total = bonded.total + van_der_waals + coulomb;
                Ok((total, "kcal/mol", total, None))
            }
            #[cfg(not(feature = "alpha-reaxff"))]
            {
                Err("ReaxFF backend requires the 'alpha-reaxff' feature flag".to_string())
            }
        }
    }
}

#[cfg(feature = "alpha-reaxff")]
fn build_reaxff_atom_params(
    elements: &[u8],
    params: &crate::forcefield::reaxff::params::ReaxffParams,
) -> Vec<crate::forcefield::reaxff::params::ReaxffAtomParams> {
    let fallback = params.atom_params.last().cloned().unwrap_or(
        crate::forcefield::reaxff::params::ReaxffAtomParams {
            element: 1,
            r_sigma: 1.0,
            r_pi: 0.0,
            r_pipi: 0.0,
            p_bo1: 0.0,
            p_bo2: 1.0,
            p_bo3: 0.0,
            p_bo4: 1.0,
            p_bo5: 0.0,
            p_bo6: 1.0,
            valence: 1.0,
        },
    );
    elements
        .iter()
        .map(|&z| {
            params
                .element_index(z)
                .map(|idx| params.atom_params[idx].clone())
                .unwrap_or_else(|| fallback.clone())
        })
        .collect()
}

fn unsupported_elements<F>(elements: &[u8], mut supported: F) -> Vec<u8>
where
    F: FnMut(u8) -> bool,
{
    let mut unsupported = Vec::new();
    for &z in elements {
        if !supported(z) && !unsupported.contains(&z) {
            unsupported.push(z);
        }
    }
    unsupported
}

fn unsupported_message(label: &str, unsupported: &[u8]) -> Vec<String> {
    if unsupported.is_empty() {
        Vec::new()
    } else {
        vec![format!(
            "{} does not support elements {:?}",
            label, unsupported
        )]
    }
}

fn flat_to_positions(coords: &[f64]) -> Vec<[f64; 3]> {
    coords
        .chunks(3)
        .map(|chunk| [chunk[0], chunk[1], chunk[2]])
        .collect()
}

fn validate_flat_coords(coords: &[f64], n_atoms: usize) -> Result<(), String> {
    if coords.len() != n_atoms * 3 {
        return Err(format!(
            "coords length {} must equal 3 * atom count {}",
            coords.len(),
            n_atoms
        ));
    }
    Ok(())
}

fn elements_from_smiles(smiles: &str) -> Result<Vec<u8>, String> {
    let mol = Molecule::from_smiles(smiles)?;
    Ok(mol
        .graph
        .node_indices()
        .map(|idx| mol.graph[idx].element)
        .collect())
}

fn capability_message(capability: &GsmBackendCapability) -> String {
    if let Some(message) = capability.limitations.first() {
        return message.clone();
    }
    if let Some(message) = capability.warnings.first() {
        return message.clone();
    }
    format!(
        "GSM backend '{}' is unavailable for this system",
        capability.backend.as_str()
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eht_is_listed_but_not_barrier_ready() {
        let capability = plan_gsm_backends(&[6, 1, 1, 1, 1])
            .into_iter()
            .find(|entry| entry.backend == GsmEnergyBackend::Eht)
            .unwrap();
        assert!(capability.available);
        assert!(!capability.suitable_for_reaction_path);
    }

    #[test]
    fn pm3_marks_unsupported_elements() {
        let capability = plan_gsm_backends(&[78])
            .into_iter()
            .find(|entry| entry.backend == GsmEnergyBackend::Pm3)
            .unwrap();
        assert!(!capability.available);
        assert!(capability
            .limitations
            .iter()
            .any(|message| message.contains("78")));
    }

    #[test]
    fn evaluate_uff_backend_on_water_geometry() {
        let evaluation = evaluate_gsm_backend(
            "[H][H]",
            &[0.0, 0.0, 0.0, 0.0, 0.0, 0.74],
            GsmEnergyBackend::Uff,
        )
        .unwrap();
        assert_eq!(evaluation.status, GsmBackendStatus::Success);
        assert!(evaluation.energy_kcal_mol.unwrap().is_finite());
    }

    #[test]
    fn compare_defaults_to_all_backends() {
        let evaluations =
            compare_gsm_backends("[H][H]", &[0.0, 0.0, 0.0, 0.0, 0.0, 0.74], &[]).unwrap();
        assert_eq!(evaluations.len(), GsmEnergyBackend::ALL.len());
    }

    #[test]
    fn ts_search_rejects_eht_backend() {
        let err = find_transition_state_with_backend(
            "[H][H]",
            &[0.0, 0.0, 0.0, 0.0, 0.0, 0.74],
            &[0.0, 0.0, 0.0, 0.0, 0.0, 1.2],
            GsmEnergyBackend::Eht,
            &GsmConfig::default(),
        )
        .unwrap_err();
        assert!(err.contains("EHT"));
    }
}
