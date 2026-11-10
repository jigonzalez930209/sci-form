//! Uncertainty calibration and active learning for regression models.
//!
//! The module adds a distribution-free conformal interval around any point
//! prediction and a query-by-committee acquisition score. It is deliberately
//! model-agnostic: the only requirement is a prediction and an uncertainty
//! estimate. The existing Random Forest already exposes both through
//! `predict_with_variance`.
//!
//! # Units
//!
//! Predictions and targets are in the unit of the regression target. The
//! conformal radius has the same unit; coverage is dimensionless.

use serde::{Deserialize, Serialize};

use crate::error::SciFormError;
use crate::ml::advanced_models::RandomForest;

/// A model that returns a point prediction and a standard-deviation estimate.
pub trait UncertaintyModel {
    /// Point prediction in the target unit.
    fn predict(&self, features: &[f64]) -> f64;

    /// Standard deviation in the target unit. It may be zero for a
    /// deterministic model, but the caller should then provide a conformal
    /// calibration set instead of trusting the interval.
    fn uncertainty(&self, features: &[f64]) -> f64;
}

impl UncertaintyModel for RandomForest {
    fn predict(&self, features: &[f64]) -> f64 {
        self.predict(features)
    }

    fn uncertainty(&self, features: &[f64]) -> f64 {
        let (_, variance) = self.predict_with_variance(features);
        variance.max(0.0).sqrt()
    }
}

/// A calibrated conformal interval.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConformalCalibrator {
    /// Absolute-residual quantile in the target unit.
    pub radius: f64,
    /// Requested coverage in `(0, 1)`.
    pub coverage: f64,
    /// Number of calibration points.
    pub n_calibration: usize,
}

impl ConformalCalibrator {
    /// Fit the conformal radius from predictions and targets.
    ///
    /// The finite-sample quantile uses the `ceil((n+1) * coverage) / n`
    /// order statistic, capped at the largest residual. With fewer than two
    /// points the radius is the largest absolute residual.
    pub fn fit(predictions: &[f64], targets: &[f64], coverage: f64) -> Result<Self, SciFormError> {
        if predictions.len() != targets.len() {
            return Err(SciFormError::validation(
                "predictions and targets must have the same length",
            ));
        }
        if predictions.is_empty() {
            return Err(SciFormError::validation(
                "at least one calibration point is required",
            ));
        }
        if !coverage.is_finite() || coverage <= 0.0 || coverage >= 1.0 {
            return Err(SciFormError::validation(
                "coverage must be finite and in (0, 1)",
            ));
        }
        let mut residuals: Vec<f64> = predictions
            .iter()
            .zip(targets)
            .map(|(prediction, target)| (prediction - target).abs())
            .collect();
        if residuals.iter().any(|value| !value.is_finite()) {
            return Err(SciFormError::validation(
                "predictions and targets must be finite",
            ));
        }
        residuals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let n = residuals.len();
        let rank = ((n + 1) as f64 * coverage).ceil() as usize;
        let radius = residuals[rank.saturating_sub(1).min(n - 1)];
        Ok(Self {
            radius,
            coverage,
            n_calibration: n,
        })
    }
}

/// A prediction with its calibrated interval and domain flag.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PredictionInterval {
    /// Point prediction.
    pub mean: f64,
    /// Model standard deviation.
    pub std: f64,
    /// Lower interval bound.
    pub lower: f64,
    /// Upper interval bound.
    pub upper: f64,
    /// Requested coverage.
    pub coverage: f64,
    /// True when the prediction is inside the labelled feature envelope.
    pub in_domain: bool,
}

/// Acquisition function for the next label.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Acquisition {
    /// Rank by model uncertainty.
    Uncertainty,
    /// Rank by distance from the labelled set.
    Distance,
    /// Rank by a weighted combination.
    Hybrid { uncertainty_weight: f64 },
}

/// One labelled point in an active-learning session.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LabeledPoint {
    /// Feature vector.
    pub features: Vec<f64>,
    /// Target value.
    pub target: f64,
}

/// A model plus its labelled set and conformal calibration.
pub struct ActiveLearningSession<M: UncertaintyModel> {
    model: M,
    labeled: Vec<LabeledPoint>,
    calibrator: Option<ConformalCalibrator>,
    acquisition: Acquisition,
    /// Maximum distance to the labelled set still considered in domain. If
    /// `None`, the 95th percentile of the labelled pairwise distance is used.
    pub domain_radius: Option<f64>,
}

impl<M: UncertaintyModel> ActiveLearningSession<M> {
    /// Create a session around a model.
    pub fn new(model: M, acquisition: Acquisition) -> Self {
        Self {
            model,
            labeled: Vec::new(),
            calibrator: None,
            acquisition,
            domain_radius: None,
        }
    }

    /// Add a labelled point and invalidate the calibration.
    pub fn add_label(&mut self, features: Vec<f64>, target: f64) {
        self.labeled.push(LabeledPoint { features, target });
        self.calibrator = None;
    }

    /// Number of labelled points.
    pub fn n_labeled(&self) -> usize {
        self.labeled.len()
    }

    /// Fit the conformal interval on the current labelled set.
    ///
    /// The model is used as-is; retraining is the caller's responsibility.
    pub fn calibrate(&mut self, coverage: f64) -> Result<ConformalCalibrator, SciFormError> {
        if self.labeled.len() < 2 {
            return Err(SciFormError::validation(
                "at least two labelled points are required for calibration",
            ));
        }
        let predictions: Vec<f64> = self
            .labeled
            .iter()
            .map(|point| self.model.predict(&point.features))
            .collect();
        let targets: Vec<f64> = self.labeled.iter().map(|point| point.target).collect();
        let calibrator = ConformalCalibrator::fit(&predictions, &targets, coverage)?;
        self.calibrator = Some(calibrator.clone());
        Ok(calibrator)
    }

    /// Predict with an interval.
    pub fn predict(&self, features: &[f64]) -> PredictionInterval {
        let mean = self.model.predict(features);
        let std = self.model.uncertainty(features).max(0.0);
        let radius = self
            .calibrator
            .as_ref()
            .map(|calibrator| calibrator.radius)
            .unwrap_or(1.96 * std);
        let distance = self.distance_to_labeled(features);
        let in_domain = match self.domain_radius {
            Some(radius) => distance <= radius,
            None => true,
        };
        PredictionInterval {
            mean,
            std,
            lower: mean - radius,
            upper: mean + radius,
            coverage: self
                .calibrator
                .as_ref()
                .map(|calibrator| calibrator.coverage)
                .unwrap_or(0.95),
            in_domain,
        }
    }

    /// Rank candidate indices for the next labels.
    pub fn propose(&self, candidates: &[Vec<f64>], n: usize) -> Vec<usize> {
        let mut scored: Vec<(usize, f64)> = candidates
            .iter()
            .enumerate()
            .map(|(index, features)| {
                let uncertainty = self.model.uncertainty(features).max(0.0);
                let distance = self.distance_to_labeled(features);
                let score = match self.acquisition {
                    Acquisition::Uncertainty => uncertainty,
                    Acquisition::Distance => distance,
                    Acquisition::Hybrid { uncertainty_weight } => {
                        uncertainty_weight * uncertainty + (1.0 - uncertainty_weight) * distance
                    }
                };
                (index, score)
            })
            .collect();
        scored.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        scored.into_iter().take(n).map(|(index, _)| index).collect()
    }

    /// Empirical coverage of the current calibrator on a test set.
    pub fn coverage(&self, predictions: &[f64], targets: &[f64]) -> Result<f64, SciFormError> {
        let calibrator = self.calibrator.as_ref().ok_or_else(|| {
            SciFormError::validation("calibrate the session before measuring coverage")
        })?;
        if predictions.len() != targets.len() || predictions.is_empty() {
            return Err(SciFormError::validation(
                "predictions and targets must be non-empty and the same length",
            ));
        }
        let covered = predictions
            .iter()
            .zip(targets)
            .filter(|(prediction, target)| (**prediction - **target).abs() <= calibrator.radius)
            .count();
        Ok(covered as f64 / predictions.len() as f64)
    }

    fn distance_to_labeled(&self, features: &[f64]) -> f64 {
        self.labeled
            .iter()
            .map(|point| euclidean(features, &point.features))
            .fold(f64::INFINITY, f64::min)
    }
}

fn euclidean(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(left, right)| (left - right) * (left - right))
        .sum::<f64>()
        .sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    struct ConstantModel {
        value: f64,
        uncertainty: f64,
    }

    impl UncertaintyModel for ConstantModel {
        fn predict(&self, _features: &[f64]) -> f64 {
            self.value
        }

        fn uncertainty(&self, _features: &[f64]) -> f64 {
            self.uncertainty
        }
    }

    #[test]
    fn conformal_radius_has_requested_coverage() {
        let predictions: Vec<f64> = (0..100).map(|value| value as f64).collect();
        let targets: Vec<f64> = predictions.iter().map(|value| value + 1.0).collect();
        let calibrator = ConformalCalibrator::fit(&predictions, &targets, 0.9).unwrap();
        assert!((calibrator.radius - 1.0).abs() < 1e-12);
        assert_eq!(calibrator.n_calibration, 100);
    }

    #[test]
    fn interval_uses_calibrated_radius() {
        let model = ConstantModel {
            value: 0.0,
            uncertainty: 0.1,
        };
        let mut session = ActiveLearningSession::new(model, Acquisition::Uncertainty);
        session.add_label(vec![0.0], 0.5);
        session.add_label(vec![1.0], -0.5);
        session.calibrate(0.9).unwrap();
        let interval = session.predict(&[2.0]);
        assert!((interval.mean - 0.0).abs() < 1e-12);
        assert!(interval.upper >= interval.lower);
        assert!((interval.coverage - 0.9).abs() < 1e-12);
    }

    #[test]
    fn acquisition_prefers_high_uncertainty() {
        struct PositionModel;
        impl UncertaintyModel for PositionModel {
            fn predict(&self, features: &[f64]) -> f64 {
                features[0]
            }
            fn uncertainty(&self, features: &[f64]) -> f64 {
                features[0].abs()
            }
        }
        let session = ActiveLearningSession::new(PositionModel, Acquisition::Uncertainty);
        let candidates = vec![vec![0.1], vec![5.0], vec![-2.0]];
        let proposed = session.propose(&candidates, 1);
        assert_eq!(proposed, vec![1]);
    }

    #[test]
    fn invalid_coverage_is_rejected() {
        assert!(ConformalCalibrator::fit(&[0.0], &[0.0], 1.0).is_err());
        assert!(ConformalCalibrator::fit(&[0.0], &[0.0], 0.0).is_err());
        assert!(ConformalCalibrator::fit(&[0.0], &[0.0], f64::NAN).is_err());
    }
}
