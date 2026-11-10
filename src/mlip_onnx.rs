//! Optional ONNX adapter for the MLIP model contract.
//!
//! This module is compiled only with the `mlip-onnx` feature. It loads a
//! user-supplied ONNX checkpoint and executes the fixed two-input contract:
//!
//! ```text
//! positions       float32 or float64  [N, 3]  Å
//! atomic_numbers  int64               [N]
//! outputs:
//! energy          float32 or float64  [1]     manifest.energy_unit
//! forces          float32 or float64  [N, 3]  manifest.energy_unit / Å
//! ```
//!
//! No model weights are bundled. The caller is responsible for the checkpoint
//! licence and for exporting it with this input/output contract.

use std::path::Path;
use std::sync::Arc;

use crate::error::SciFormError;
use crate::mlip::{MlipEnergyUnit, MlipManifest, MlipModel, MlipResult};
use tract_onnx::prelude::{Framework, InferenceModelExt, IntoRunnable, TypedFact, TypedOp};

const EV_TO_KCAL_MOL: f64 = 23.060_548_867;

/// ONNX-backed MLIP.
pub struct MlipOnnx {
    manifest: MlipManifest,
    model: Arc<tract_onnx::prelude::RunnableModel<TypedFact, Box<dyn TypedOp>>>,
}

impl MlipOnnx {
    /// Load an ONNX checkpoint and validate the manifest.
    ///
    /// The model must expose exactly two outputs: energy and forces. Inputs are
    /// `positions` and `atomic_numbers`; the order is the ONNX input order.
    pub fn from_path(path: &Path, manifest: MlipManifest) -> Result<Self, SciFormError> {
        manifest.validate()?;
        if manifest.supports_stress {
            return Err(SciFormError::unsupported(
                "the ONNX adapter currently supports energy/forces only; stress is a follow-up",
            ));
        }
        let model = tract_onnx::onnx()
            .model_for_path(path)
            .map_err(|error| {
                SciFormError::validation(format!("could not load ONNX model: {error}"))
            })?
            .into_optimized()
            .map_err(|error| {
                SciFormError::numerical(format!("could not optimize ONNX model: {error}"))
            })?
            .into_runnable()
            .map_err(|error| {
                SciFormError::numerical(format!("could not make ONNX model runnable: {error}"))
            })?;
        Ok(Self { manifest, model })
    }

    fn parse_f64_tensor(
        value: &tract_onnx::prelude::TValue,
        what: &str,
    ) -> Result<Vec<f64>, SciFormError> {
        if let Ok(view) = value.to_plain_array_view::<f32>() {
            return Ok(view.iter().map(|value| *value as f64).collect());
        }
        if let Ok(view) = value.to_plain_array_view::<f64>() {
            return Ok(view.iter().copied().collect());
        }
        Err(SciFormError::numerical(format!(
            "ONNX output '{what}' is not a float32/float64 tensor"
        )))
    }
}

impl MlipModel for MlipOnnx {
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
        if cell.is_some() {
            return Err(SciFormError::unsupported(
                "the current ONNX adapter is non-periodic; cell input is a follow-up",
            ));
        }
        if !self.manifest.covers(elements) {
            return Err(SciFormError::unsupported(
                "the ONNX model does not cover every element in the system",
            ));
        }
        if self.model.model().inputs.len() != 2 {
            return Err(SciFormError::validation(format!(
                "the ONNX MLIP contract requires two inputs, found {}",
                self.model.model().inputs.len()
            )));
        }

        use tract_onnx::prelude::*;
        let n = elements.len();
        let positions: Vec<f32> = positions_angstrom
            .iter()
            .flat_map(|position| position.iter().map(|value| *value as f32))
            .collect();
        let positions = tract_ndarray::Array2::from_shape_vec((n, 3), positions)
            .map_err(|error| SciFormError::validation(format!("positions shape: {error}")))?
            .into_tensor();
        let atomic_numbers: Vec<i64> = elements.iter().map(|value| *value as i64).collect();
        let atomic_numbers = tract_ndarray::Array1::from_vec(atomic_numbers).into_tensor();

        let outputs = self
            .model
            .run(tvec!(positions.into(), atomic_numbers.into()))
            .map_err(|error| SciFormError::numerical(format!("ONNX inference failed: {error}")))?;
        if outputs.len() < 2 {
            return Err(SciFormError::numerical(
                "the ONNX MLIP must return energy and forces",
            ));
        }
        let energy_values = Self::parse_f64_tensor(&outputs[0], "energy")?;
        let force_values = Self::parse_f64_tensor(&outputs[1], "forces")?;
        if energy_values.len() != 1 {
            return Err(SciFormError::numerical(
                "the ONNX energy output must contain exactly one value",
            ));
        }
        if force_values.len() != n * 3 {
            return Err(SciFormError::numerical(format!(
                "the ONNX forces output has {} values, expected {}",
                force_values.len(),
                n * 3
            )));
        }

        let factor = match self.manifest.energy_unit {
            MlipEnergyUnit::KcalMol => 1.0,
            MlipEnergyUnit::Ev => EV_TO_KCAL_MOL,
        };
        let (force_chunks, _) = force_values.as_chunks::<3>();
        let forces = force_chunks
            .iter()
            .map(|chunk| [chunk[0] * factor, chunk[1] * factor, chunk[2] * factor])
            .collect();
        Ok(MlipResult {
            energy_kcal_mol: energy_values[0] * factor,
            forces_kcal_mol_per_angstrom: forces,
            stress_kcal_mol_per_angstrom3: None,
            model: self.manifest.name.clone(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn toy_checkpoint_runs_energy_and_forces() {
        let manifest = MlipManifest {
            schema_version: 1,
            name: "toy".into(),
            elements: vec![1],
            cutoff_angstrom: 5.0,
            license: "test-only".into(),
            energy_unit: MlipEnergyUnit::KcalMol,
            supports_stress: false,
            source: None,
        };
        let path =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/mlip_toy.onnx");
        let model = MlipOnnx::from_path(&path, manifest).unwrap();
        let result = model
            .evaluate(&[1, 1], &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], None)
            .unwrap();
        assert!((result.energy_kcal_mol - 1.0).abs() < 1e-5);
        assert_eq!(result.forces_kcal_mol_per_angstrom.len(), 2);
        assert!(result
            .forces_kcal_mol_per_angstrom
            .iter()
            .all(|force| force.iter().all(|value| value.abs() < 1e-6)));
    }

    #[test]
    fn missing_checkpoint_is_reported() {
        let manifest = MlipManifest {
            schema_version: 1,
            name: "missing".into(),
            elements: vec![1, 8],
            cutoff_angstrom: 5.0,
            license: "user-supplied".into(),
            energy_unit: MlipEnergyUnit::KcalMol,
            supports_stress: false,
            source: None,
        };
        let error = MlipOnnx::from_path(Path::new("/nonexistent/model.onnx"), manifest)
            .err()
            .expect("missing checkpoint must fail");
        assert!(error.to_string().contains("ONNX"));
    }
}
