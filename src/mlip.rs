//! Universal MLIP model contract and registry.
//!
//! The library does not bundle a trained universal potential. This module
//! defines the contract that a user-supplied model (MACE, CHGNet, Orb, eSEN or
//! a private checkpoint) must satisfy, plus a callback-backed implementation
//! for bindings and tests. The ONNX adapter is a follow-up; it is intentionally
//! not hidden behind a hard-coded model.
//!
//! # Contract
//!
//! A model evaluates
//!
//! ```text
//! (elements, positions, cell?) -> (energy, forces, stress?)
//! ```
//!
//! with positions in Å, energy in kcal/mol, forces in kcal/mol/Å and stress in
//! kcal/mol/Å³. The manifest declares the element coverage, cutoff, licence and
//! whether the model predicts stress.
//!
//! # Units
//!
//! The contract uses kcal/mol because that is the native unit of the classical
//! force fields and the MD integrators. An ONNX model trained in eV must convert
//! at the adapter boundary; the conversion is not hidden in the registry.

use serde::{Deserialize, Serialize};

use crate::error::SciFormError;

/// Native energy unit of a model checkpoint.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum MlipEnergyUnit {
    /// kcal/mol; the native contract unit.
    #[default]
    KcalMol,
    /// electron-volt; the adapter converts at the boundary.
    Ev,
}

/// Manifest that accompanies a model checkpoint.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MlipManifest {
    /// Contract version. The current version is `1`.
    pub schema_version: u32,
    /// Human-readable model name.
    pub name: String,
    /// Atomic numbers the model was trained on.
    pub elements: Vec<u8>,
    /// Neighbour cutoff in Å.
    pub cutoff_angstrom: f64,
    /// Licence identifier or a short licence note.
    pub license: String,
    /// Native energy unit of the checkpoint.
    #[serde(default)]
    pub energy_unit: MlipEnergyUnit,
    /// Whether the model returns a stress tensor.
    pub supports_stress: bool,
    /// Optional source or checkpoint hash.
    pub source: Option<String>,
}

impl MlipManifest {
    /// Validate the manifest before a model is registered.
    pub fn validate(&self) -> Result<(), SciFormError> {
        if self.schema_version != 1 {
            return Err(SciFormError::unsupported(format!(
                "unsupported MLIP schema version {}; expected 1",
                self.schema_version
            )));
        }
        if self.name.trim().is_empty() {
            return Err(SciFormError::validation("MLIP name must not be empty"));
        }
        if self.elements.is_empty() {
            return Err(SciFormError::validation(
                "MLIP manifest must declare at least one element",
            ));
        }
        if self.elements.iter().any(|z| *z == 0 || *z > 118) {
            return Err(SciFormError::validation(
                "MLIP element list contains an invalid atomic number",
            ));
        }
        if !self.cutoff_angstrom.is_finite() || self.cutoff_angstrom <= 0.0 {
            return Err(SciFormError::validation(
                "MLIP cutoff must be finite and > 0",
            ));
        }
        if self.license.trim().is_empty() {
            return Err(SciFormError::validation(
                "MLIP licence must be declared; model weights are never bundled silently",
            ));
        }
        Ok(())
    }

    /// True when every element is covered by the model.
    pub fn covers(&self, elements: &[u8]) -> bool {
        elements.iter().all(|z| self.elements.contains(z))
    }
}

/// Energy, forces and optional stress from one model evaluation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MlipResult {
    /// Total energy in kcal/mol.
    pub energy_kcal_mol: f64,
    /// Forces in kcal/mol/Å, one vector per atom.
    pub forces_kcal_mol_per_angstrom: Vec<[f64; 3]>,
    /// Optional stress tensor in kcal/mol/Å³.
    pub stress_kcal_mol_per_angstrom3: Option<[[f64; 3]; 3]>,
    /// Model name copied from the manifest.
    pub model: String,
}

/// A model that satisfies the MLIP contract.
pub trait MlipModel {
    /// Model manifest.
    fn manifest(&self) -> &MlipManifest;

    /// Evaluate energy, forces and optional stress.
    fn evaluate(
        &self,
        elements: &[u8],
        positions_angstrom: &[[f64; 3]],
        cell: Option<&[[f64; 3]; 3]>,
    ) -> Result<MlipResult, SciFormError>;
}

/// A model backed by a user closure.
pub struct CallbackMlip<F>
where
    F: Fn(
        &[u8],
        &[[f64; 3]],
        Option<&[[f64; 3]; 3]>,
    ) -> Result<(f64, Vec<[f64; 3]>, Option<[[f64; 3]; 3]>), SciFormError>,
{
    manifest: MlipManifest,
    evaluator: F,
}

impl<F> CallbackMlip<F>
where
    F: Fn(
        &[u8],
        &[[f64; 3]],
        Option<&[[f64; 3]; 3]>,
    ) -> Result<(f64, Vec<[f64; 3]>, Option<[[f64; 3]; 3]>), SciFormError>,
{
    /// Create a callback-backed model after validating its manifest.
    pub fn new(manifest: MlipManifest, evaluator: F) -> Result<Self, SciFormError> {
        manifest.validate()?;
        Ok(Self {
            manifest,
            evaluator,
        })
    }
}

impl<F> MlipModel for CallbackMlip<F>
where
    F: Fn(
        &[u8],
        &[[f64; 3]],
        Option<&[[f64; 3]; 3]>,
    ) -> Result<(f64, Vec<[f64; 3]>, Option<[[f64; 3]; 3]>), SciFormError>,
{
    fn manifest(&self) -> &MlipManifest {
        &self.manifest
    }

    fn evaluate(
        &self,
        elements: &[u8],
        positions_angstrom: &[[f64; 3]],
        cell: Option<&[[f64; 3]; 3]>,
    ) -> Result<MlipResult, SciFormError> {
        crate::validation::validate_geometry(elements, positions_angstrom)?;
        if !self.manifest.covers(elements) {
            let missing: Vec<u8> = elements
                .iter()
                .copied()
                .filter(|z| !self.manifest.elements.contains(z))
                .collect();
            return Err(SciFormError::unsupported(format!(
                "model '{}' does not cover element(s) {missing:?}",
                self.manifest.name
            )));
        }
        if self.manifest.supports_stress && cell.is_none() {
            return Err(SciFormError::validation(format!(
                "model '{}' predicts stress but no cell was supplied",
                self.manifest.name
            )));
        }
        let (energy, forces, stress) = (self.evaluator)(elements, positions_angstrom, cell)?;
        if !energy.is_finite() {
            return Err(SciFormError::numerical("MLIP energy is not finite"));
        }
        if forces.len() != elements.len()
            || forces
                .iter()
                .any(|force| force.iter().any(|value| !value.is_finite()))
        {
            return Err(SciFormError::numerical(
                "MLIP forces must contain one finite vector per atom",
            ));
        }
        if let Some(stress) = stress {
            if stress
                .iter()
                .any(|row| row.iter().any(|value| !value.is_finite()))
            {
                return Err(SciFormError::numerical("MLIP stress is not finite"));
            }
        } else if self.manifest.supports_stress {
            return Err(SciFormError::numerical(
                "MLIP manifest declares stress but the model returned none",
            ));
        }
        Ok(MlipResult {
            energy_kcal_mol: energy,
            forces_kcal_mol_per_angstrom: forces,
            stress_kcal_mol_per_angstrom3: stress,
            model: self.manifest.name.clone(),
        })
    }
}

/// A registry of named models.
#[derive(Default)]
pub struct MlipRegistry {
    models: Vec<Box<dyn MlipModel>>,
}

impl MlipRegistry {
    /// Create an empty registry.
    pub fn new() -> Self {
        Self { models: Vec::new() }
    }

    /// Register a model. Names must be unique.
    pub fn register<M: MlipModel + 'static>(&mut self, model: M) -> Result<(), SciFormError> {
        if self
            .models
            .iter()
            .any(|existing| existing.manifest().name == model.manifest().name)
        {
            return Err(SciFormError::validation(format!(
                "an MLIP named '{}' is already registered",
                model.manifest().name
            )));
        }
        model.manifest().validate()?;
        self.models.push(Box::new(model));
        Ok(())
    }

    /// Look up a model by name.
    pub fn get(&self, name: &str) -> Option<&dyn MlipModel> {
        self.models
            .iter()
            .find(|model| model.manifest().name == name)
            .map(|model| model.as_ref())
    }

    /// Manifest of every registered model.
    pub fn manifests(&self) -> Vec<&MlipManifest> {
        self.models.iter().map(|model| model.manifest()).collect()
    }

    /// Number of registered models.
    pub fn len(&self) -> usize {
        self.models.len()
    }

    /// True when no model is registered.
    pub fn is_empty(&self) -> bool {
        self.models.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn manifest() -> MlipManifest {
        MlipManifest {
            schema_version: 1,
            name: "test-lj".into(),
            elements: vec![1, 8],
            cutoff_angstrom: 6.0,
            license: "test-only".into(),
            energy_unit: MlipEnergyUnit::KcalMol,
            supports_stress: false,
            source: None,
        }
    }

    fn lj_model() -> CallbackMlip<
        impl Fn(
            &[u8],
            &[[f64; 3]],
            Option<&[[f64; 3]; 3]>,
        ) -> Result<(f64, Vec<[f64; 3]>, Option<[[f64; 3]; 3]>), SciFormError>,
    > {
        CallbackMlip::new(
            manifest(),
            |elements: &[u8], positions: &[[f64; 3]], _cell: Option<&[[f64; 3]; 3]>| {
                let mut energy = 0.0;
                for i in 0..positions.len() {
                    for j in (i + 1)..positions.len() {
                        let dx = positions[i][0] - positions[j][0];
                        let dy = positions[i][1] - positions[j][1];
                        let dz = positions[i][2] - positions[j][2];
                        let r2 = dx * dx + dy * dy + dz * dz;
                        energy += 1.0 / r2.max(0.25);
                    }
                }
                let forces = vec![[0.0, 0.0, 0.0]; elements.len()];
                Ok((energy, forces, None))
            },
        )
        .unwrap()
    }

    #[test]
    fn callback_model_satisfies_the_contract() {
        let model = lj_model();
        let result = model
            .evaluate(&[8, 1], &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], None)
            .unwrap();
        assert!(result.energy_kcal_mol > 0.0);
        assert_eq!(result.forces_kcal_mol_per_angstrom.len(), 2);
        assert_eq!(result.model, "test-lj");
    }

    #[test]
    fn missing_element_is_rejected() {
        let model = lj_model();
        assert!(model
            .evaluate(&[6, 1], &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], None)
            .is_err());
    }

    #[test]
    fn registry_rejects_duplicate_names() {
        let mut registry = MlipRegistry::new();
        registry.register(lj_model()).unwrap();
        assert!(registry.register(lj_model()).is_err());
        assert_eq!(registry.len(), 1);
        assert!(registry.get("test-lj").is_some());
    }

    #[test]
    fn invalid_manifest_is_rejected() {
        let mut bad = manifest();
        bad.elements.clear();
        assert!(bad.validate().is_err());
    }
}
