//! Advanced ML models: Random Forest, Gradient Boosting, and cross-validation.
//!
//! Pure-Rust implementations for molecular property prediction beyond
//! simple linear models. Includes:
//! - Random Forest (bagged decision trees)
//! - Gradient Boosted Trees (GBM / GBRT)
//! - K-fold cross-validation
//! - Model recalibration via isotonic regression

use serde::{Deserialize, Serialize};

// ─── Decision Tree ───────────────────────────────────────────────────────────

/// A single decision tree node (binary split).
///
/// The tree is a binary regression tree: every internal node splits the
/// samples on one feature at one threshold and sends the lower branch
/// (`value <= threshold`) left and the upper branch right.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TreeNode {
    /// Terminal node holding the constant prediction for its samples.
    Leaf {
        /// Predicted value, in the units of the regression target.
        value: f64,
    },
    /// Internal node that splits on one feature.
    Split {
        /// Index of the split feature in the feature vector.
        feature: usize,
        /// Split threshold, in the units of the split feature.
        threshold: f64,
        /// Branch taken when `feature <= threshold`.
        left: Box<TreeNode>,
        /// Branch taken when `feature > threshold`.
        right: Box<TreeNode>,
    },
}

impl TreeNode {
    /// Predict a single sample.
    pub fn predict(&self, features: &[f64]) -> f64 {
        match self {
            TreeNode::Leaf { value } => *value,
            TreeNode::Split {
                feature,
                threshold,
                left,
                right,
            } => {
                if features.get(*feature).copied().unwrap_or(0.0) <= *threshold {
                    left.predict(features)
                } else {
                    right.predict(features)
                }
            }
        }
    }
}

/// Configuration for tree building.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct TreeConfig {
    /// Maximum depth of the tree.
    pub max_depth: usize,
    /// Minimum samples to split a node.
    pub min_samples_split: usize,
    /// Minimum samples in a leaf.
    pub min_samples_leaf: usize,
    /// Number of features to consider at each split (0 = sqrt(n_features)).
    pub max_features: usize,
}

impl Default for TreeConfig {
    fn default() -> Self {
        Self {
            max_depth: 10,
            min_samples_split: 5,
            min_samples_leaf: 2,
            max_features: 0,
        }
    }
}

/// Build a decision tree from data.
pub fn build_tree(
    features: &[Vec<f64>],
    targets: &[f64],
    config: &TreeConfig,
    rng_seed: u64,
) -> TreeNode {
    let indices: Vec<usize> = (0..targets.len()).collect();
    let n_features = features.first().map_or(0, |f| f.len());
    let max_feat = if config.max_features == 0 {
        (n_features as f64).sqrt().ceil() as usize
    } else {
        config.max_features.min(n_features)
    }
    .max(1);
    build_tree_recursive(features, targets, &indices, config, max_feat, 0, rng_seed)
}

fn build_tree_recursive(
    features: &[Vec<f64>],
    targets: &[f64],
    indices: &[usize],
    config: &TreeConfig,
    max_features: usize,
    depth: usize,
    seed: u64,
) -> TreeNode {
    let n = indices.len();

    // Leaf conditions
    if n < config.min_samples_split || depth >= config.max_depth || n < 2 * config.min_samples_leaf
    {
        let mean = indices.iter().map(|&i| targets[i]).sum::<f64>() / n.max(1) as f64;
        return TreeNode::Leaf { value: mean };
    }

    let n_features = features.first().map_or(0, |f| f.len());
    if n_features == 0 {
        let mean = indices.iter().map(|&i| targets[i]).sum::<f64>() / n as f64;
        return TreeNode::Leaf { value: mean };
    }

    // Select random subset of features
    let feature_subset = select_features(n_features, max_features, seed);

    // Find best split
    let mut best_score = f64::INFINITY;
    let mut best_feature = 0;
    let mut best_threshold = 0.0;
    let mut best_left = vec![];
    let mut best_right = vec![];

    for &feat in &feature_subset {
        // Get unique sorted values for this feature
        let mut vals: Vec<f64> = indices.iter().map(|&i| features[i][feat]).collect();
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        vals.dedup();

        if vals.len() < 2 {
            continue;
        }

        // Try midpoints as thresholds
        for w in vals.windows(2) {
            let threshold = (w[0] + w[1]) / 2.0;
            let (left, right): (Vec<usize>, Vec<usize>) = indices
                .iter()
                .partition(|&&i| features[i][feat] <= threshold);

            if left.len() < config.min_samples_leaf || right.len() < config.min_samples_leaf {
                continue;
            }

            let score = mse_score(&left, targets) + mse_score(&right, targets);
            if score < best_score {
                best_score = score;
                best_feature = feat;
                best_threshold = threshold;
                best_left = left;
                best_right = right;
            }
        }
    }

    if best_left.is_empty() || best_right.is_empty() {
        let mean = indices.iter().map(|&i| targets[i]).sum::<f64>() / n as f64;
        return TreeNode::Leaf { value: mean };
    }

    let left_node = build_tree_recursive(
        features,
        targets,
        &best_left,
        config,
        max_features,
        depth + 1,
        seed.wrapping_mul(6364136223846793005).wrapping_add(1),
    );
    let right_node = build_tree_recursive(
        features,
        targets,
        &best_right,
        config,
        max_features,
        depth + 1,
        seed.wrapping_mul(6364136223846793005).wrapping_add(3),
    );

    TreeNode::Split {
        feature: best_feature,
        threshold: best_threshold,
        left: Box::new(left_node),
        right: Box::new(right_node),
    }
}

fn mse_score(indices: &[usize], targets: &[f64]) -> f64 {
    let n = indices.len() as f64;
    if n < 1.0 {
        return 0.0;
    }
    let mean = indices.iter().map(|&i| targets[i]).sum::<f64>() / n;
    indices
        .iter()
        .map(|&i| {
            let d = targets[i] - mean;
            d * d
        })
        .sum::<f64>()
}

fn select_features(n_features: usize, max_features: usize, seed: u64) -> Vec<usize> {
    if max_features >= n_features {
        return (0..n_features).collect();
    }

    // Simple LCG-based selection without replacement
    let mut selected = Vec::with_capacity(max_features);
    let mut available: Vec<usize> = (0..n_features).collect();
    let mut s = seed;
    for _ in 0..max_features {
        if available.is_empty() {
            break;
        }
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let idx = (s >> 33) as usize % available.len();
        selected.push(available.swap_remove(idx));
    }
    selected
}

// ─── Random Forest ───────────────────────────────────────────────────────────

/// Random Forest model: ensemble of bagged decision trees.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RandomForest {
    /// The bagged trees; their predictions are averaged.
    pub trees: Vec<TreeNode>,
    /// Number of trees in the forest (equal to the number requested in the
    /// configuration).
    pub n_trees: usize,
    /// Out-of-bag mean squared error, computed from the samples each tree's
    /// bootstrap did not draw. `None` when no sample was ever out of bag
    /// (e.g. a bootstrap that covers the whole training set).
    ///
    /// The value has the unit of the regression target squared (e.g. eV² for
    /// an energy target).
    pub oob_score: Option<f64>,
}

/// Configuration for Random Forest.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct RandomForestConfig {
    /// Number of trees to grow.
    pub n_trees: usize,
    /// Tree-growth parameters applied to every tree of the forest.
    pub tree_config: TreeConfig,
    /// Fraction of samples for each bootstrap (0.0–1.0). Default: 1.0
    pub sample_fraction: f64,
    /// Seed of the deterministic pseudo-random generator that drives the
    /// bootstrap sampling and the feature selection.
    pub seed: u64,
}

impl Default for RandomForestConfig {
    fn default() -> Self {
        Self {
            n_trees: 100,
            tree_config: TreeConfig::default(),
            sample_fraction: 1.0,
            seed: 42,
        }
    }
}

/// Train a Random Forest regressor.
pub fn train_random_forest(
    features: &[Vec<f64>],
    targets: &[f64],
    config: &RandomForestConfig,
) -> RandomForest {
    let n = targets.len();
    let n_sample = ((n as f64 * config.sample_fraction).ceil() as usize).max(1);
    let mut trees = Vec::with_capacity(config.n_trees);

    // Out-of-bag accumulation: every tree predicts the samples that its
    // bootstrap did not draw, which gives an unbiased estimate of the
    // generalisation error (the "OOB score").
    let mut oob_sq_error = 0.0f64;
    let mut oob_count = 0usize;

    for t in 0..config.n_trees {
        let seed = config.seed.wrapping_add(t as u64 * 1000003);
        // Bootstrap sample
        let bootstrap = bootstrap_indices(n, n_sample, seed);
        let boot_features: Vec<Vec<f64>> = bootstrap.iter().map(|&i| features[i].clone()).collect();
        let boot_targets: Vec<f64> = bootstrap.iter().map(|&i| targets[i]).collect();
        let tree = build_tree(&boot_features, &boot_targets, &config.tree_config, seed);

        // Samples that are out of the bag for this tree.
        let mut in_bag = vec![false; n];
        for &i in &bootstrap {
            in_bag[i] = true;
        }
        for i in 0..n {
            if in_bag[i] {
                continue;
            }
            let pred = tree.predict(&features[i]);
            let err = pred - targets[i];
            oob_sq_error += err * err;
            oob_count += 1;
        }

        trees.push(tree);
    }

    let oob_score = if oob_count > 0 {
        Some(oob_sq_error / oob_count as f64)
    } else {
        None
    };

    RandomForest {
        n_trees: config.n_trees,
        trees,
        oob_score,
    }
}

impl RandomForest {
    /// Predict a single sample (mean of all trees).
    pub fn predict(&self, features: &[f64]) -> f64 {
        let sum: f64 = self.trees.iter().map(|t| t.predict(features)).sum();
        sum / self.trees.len().max(1) as f64
    }

    /// Predict and return individual tree predictions (for uncertainty).
    pub fn predict_with_variance(&self, features: &[f64]) -> (f64, f64) {
        let preds: Vec<f64> = self.trees.iter().map(|t| t.predict(features)).collect();
        let mean = preds.iter().sum::<f64>() / preds.len().max(1) as f64;
        let var =
            preds.iter().map(|p| (p - mean) * (p - mean)).sum::<f64>() / preds.len().max(1) as f64;
        (mean, var)
    }
}

fn bootstrap_indices(n: usize, n_sample: usize, seed: u64) -> Vec<usize> {
    let mut indices = Vec::with_capacity(n_sample);
    let mut s = seed;
    for _ in 0..n_sample {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        indices.push((s >> 33) as usize % n);
    }
    indices
}

// ─── Gradient Boosting ───────────────────────────────────────────────────────

/// Gradient Boosted Trees regressor.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GradientBoosting {
    /// The sequentially fitted regression trees; their contributions are
    /// summed with the learning rate.
    pub trees: Vec<TreeNode>,
    /// Shrinkage factor applied to every tree contribution (dimensionless).
    pub learning_rate: f64,
    /// Constant offset the model starts from, i.e. the mean of the training
    /// targets, in the units of the regression target.
    pub initial_value: f64,
    /// Number of boosting stages (trees) that were fitted.
    pub n_estimators: usize,
}

/// Configuration for Gradient Boosting.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct GradientBoostingConfig {
    /// Number of boosting stages to fit.
    pub n_estimators: usize,
    /// Shrinkage factor applied to each tree contribution (dimensionless).
    pub learning_rate: f64,
    /// Tree-growth parameters applied to every stage.
    pub tree_config: TreeConfig,
    /// Fraction of the training set used at each stage (0.0–1.0).
    pub subsample: f64,
    /// Seed of the deterministic pseudo-random generator that drives the
    /// stage subsampling.
    pub seed: u64,
}

impl Default for GradientBoostingConfig {
    fn default() -> Self {
        Self {
            n_estimators: 100,
            learning_rate: 0.1,
            tree_config: TreeConfig {
                max_depth: 5,
                min_samples_split: 5,
                min_samples_leaf: 2,
                max_features: 0,
            },
            subsample: 0.8,
            seed: 42,
        }
    }
}

/// Train a Gradient Boosting regressor (L2 loss / least squares).
pub fn train_gradient_boosting(
    features: &[Vec<f64>],
    targets: &[f64],
    config: &GradientBoostingConfig,
) -> GradientBoosting {
    let n = targets.len();
    let initial_value = targets.iter().sum::<f64>() / n.max(1) as f64;
    let mut predictions = vec![initial_value; n];
    let mut trees = Vec::with_capacity(config.n_estimators);

    for t in 0..config.n_estimators {
        // Compute negative gradient (residuals for L2 loss)
        let residuals: Vec<f64> = (0..n).map(|i| targets[i] - predictions[i]).collect();

        // Subsample
        let seed = config.seed.wrapping_add(t as u64 * 999983);
        let n_sub = ((n as f64 * config.subsample).ceil() as usize).max(1);
        let sub_idx = bootstrap_indices(n, n_sub, seed);

        let sub_features: Vec<Vec<f64>> = sub_idx.iter().map(|&i| features[i].clone()).collect();
        let sub_residuals: Vec<f64> = sub_idx.iter().map(|&i| residuals[i]).collect();

        let tree = build_tree(&sub_features, &sub_residuals, &config.tree_config, seed);

        // Update predictions
        for i in 0..n {
            predictions[i] += config.learning_rate * tree.predict(&features[i]);
        }

        trees.push(tree);
    }

    GradientBoosting {
        trees,
        learning_rate: config.learning_rate,
        initial_value,
        n_estimators: config.n_estimators,
    }
}

impl GradientBoosting {
    /// Predict a single sample.
    pub fn predict(&self, features: &[f64]) -> f64 {
        self.initial_value
            + self.learning_rate * self.trees.iter().map(|t| t.predict(features)).sum::<f64>()
    }
}

// ─── Cross-Validation ────────────────────────────────────────────────────────

/// Result of k-fold cross-validation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CrossValidationResult {
    /// Mean absolute error per fold.
    pub fold_mae: Vec<f64>,
    /// Root mean squared error per fold.
    pub fold_rmse: Vec<f64>,
    /// R² per fold.
    pub fold_r2: Vec<f64>,
    /// Overall MAE (mean across folds).
    pub mean_mae: f64,
    /// Overall RMSE (mean across folds).
    pub mean_rmse: f64,
    /// Overall R² (mean across folds).
    pub mean_r2: f64,
    /// Number of folds.
    pub k: usize,
}

/// Model type for cross-validation.
///
/// Selects which regressor [`cross_validate`] trains inside each fold.
pub enum ModelType {
    /// Train a [`RandomForest`] in every fold.
    RandomForest(RandomForestConfig),
    /// Train a [`GradientBoosting`] regressor in every fold.
    GradientBoosting(GradientBoostingConfig),
}

/// Perform k-fold cross-validation.
pub fn cross_validate(
    features: &[Vec<f64>],
    targets: &[f64],
    model_type: &ModelType,
    k: usize,
    seed: u64,
) -> CrossValidationResult {
    let n = targets.len();
    let k = k.max(2).min(n);

    // Create fold assignments (shuffled)
    let mut indices: Vec<usize> = (0..n).collect();
    // Fisher-Yates shuffle
    let mut s = seed;
    for i in (1..n).rev() {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let j = (s >> 33) as usize % (i + 1);
        indices.swap(i, j);
    }

    let fold_size = n / k;
    let mut fold_mae = Vec::with_capacity(k);
    let mut fold_rmse = Vec::with_capacity(k);
    let mut fold_r2 = Vec::with_capacity(k);

    for fold in 0..k {
        let test_start = fold * fold_size;
        let test_end = if fold == k - 1 {
            n
        } else {
            (fold + 1) * fold_size
        };

        let test_indices: Vec<usize> = indices[test_start..test_end].to_vec();
        let train_indices: Vec<usize> = indices[..test_start]
            .iter()
            .chain(indices[test_end..].iter())
            .copied()
            .collect();

        let train_features: Vec<Vec<f64>> =
            train_indices.iter().map(|&i| features[i].clone()).collect();
        let train_targets: Vec<f64> = train_indices.iter().map(|&i| targets[i]).collect();

        // Train model
        let predictions: Vec<f64> = match model_type {
            ModelType::RandomForest(config) => {
                let model = train_random_forest(&train_features, &train_targets, config);
                test_indices
                    .iter()
                    .map(|&i| model.predict(&features[i]))
                    .collect()
            }
            ModelType::GradientBoosting(config) => {
                let model = train_gradient_boosting(&train_features, &train_targets, config);
                test_indices
                    .iter()
                    .map(|&i| model.predict(&features[i]))
                    .collect()
            }
        };

        let test_targets: Vec<f64> = test_indices.iter().map(|&i| targets[i]).collect();
        let (mae, rmse, r2) = compute_metrics(&predictions, &test_targets);
        fold_mae.push(mae);
        fold_rmse.push(rmse);
        fold_r2.push(r2);
    }

    let mean_mae = fold_mae.iter().sum::<f64>() / k as f64;
    let mean_rmse = fold_rmse.iter().sum::<f64>() / k as f64;
    let mean_r2 = fold_r2.iter().sum::<f64>() / k as f64;

    CrossValidationResult {
        fold_mae,
        fold_rmse,
        fold_r2,
        mean_mae,
        mean_rmse,
        mean_r2,
        k,
    }
}

fn compute_metrics(predictions: &[f64], targets: &[f64]) -> (f64, f64, f64) {
    let n = predictions.len() as f64;
    if n < 1.0 {
        return (0.0, 0.0, 0.0);
    }

    let mut sum_ae = 0.0;
    let mut sum_se = 0.0;
    let mean_t = targets.iter().sum::<f64>() / n;
    let mut ss_tot = 0.0;

    for i in 0..predictions.len() {
        let err = predictions[i] - targets[i];
        sum_ae += err.abs();
        sum_se += err * err;
        ss_tot += (targets[i] - mean_t) * (targets[i] - mean_t);
    }

    let mae = sum_ae / n;
    let rmse = (sum_se / n).sqrt();
    let r2 = if ss_tot > 1e-12 {
        1.0 - sum_se / ss_tot
    } else {
        0.0
    };

    (mae, rmse, r2)
}

// ─── Isotonic Regression Recalibration ───────────────────────────────────────

/// Isotonic regression for model recalibration.
/// Fits a monotone non-decreasing function to prediction-target pairs.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsotonicCalibrator {
    /// The sorted prediction values that act as knots, in the units of the
    /// model prediction.
    pub knots_x: Vec<f64>,
    /// The fitted target value at each knot, in the units of the regression
    /// target.
    pub knots_y: Vec<f64>,
}

impl IsotonicCalibrator {
    /// Fit isotonic regression from (prediction, target) pairs.
    pub fn fit(predictions: &[f64], targets: &[f64]) -> Self {
        let n = predictions.len().min(targets.len());
        if n == 0 {
            return Self {
                knots_x: vec![],
                knots_y: vec![],
            };
        }

        // Sort by prediction
        let mut pairs: Vec<(f64, f64)> = predictions[..n]
            .iter()
            .zip(targets[..n].iter())
            .map(|(&p, &t)| (p, t))
            .collect();
        pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        // Pool adjacent violators (PAVA)
        let mut blocks: Vec<(f64, f64, usize)> = pairs.iter().map(|&(x, y)| (x, y, 1)).collect(); // (x, y_mean, count)

        let mut i = 0;
        while i < blocks.len() - 1 {
            if blocks[i].1 > blocks[i + 1].1 {
                // Merge blocks
                let n1 = blocks[i].2 as f64;
                let n2 = blocks[i + 1].2 as f64;
                let merged_y = (n1 * blocks[i].1 + n2 * blocks[i + 1].1) / (n1 + n2);
                blocks[i].1 = merged_y;
                blocks[i].2 += blocks[i + 1].2;
                blocks.remove(i + 1);
                if i > 0 {
                    i = i.saturating_sub(1);
                }
            } else {
                i += 1;
            }
        }

        let knots_x: Vec<f64> = blocks.iter().map(|b| b.0).collect();
        let knots_y: Vec<f64> = blocks.iter().map(|b| b.1).collect();

        Self { knots_x, knots_y }
    }

    /// Calibrate a prediction using the fitted isotonic function.
    pub fn calibrate(&self, prediction: f64) -> f64 {
        if self.knots_x.is_empty() {
            return prediction;
        }
        if prediction <= self.knots_x[0] {
            return self.knots_y[0];
        }
        if prediction >= *self.knots_x.last().unwrap() {
            return *self.knots_y.last().unwrap();
        }

        // Linear interpolation between knots
        for i in 0..self.knots_x.len() - 1 {
            if prediction >= self.knots_x[i] && prediction <= self.knots_x[i + 1] {
                let t = (prediction - self.knots_x[i]) / (self.knots_x[i + 1] - self.knots_x[i]);
                return self.knots_y[i] + t * (self.knots_y[i + 1] - self.knots_y[i]);
            }
        }

        prediction
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decision_tree() {
        let features = vec![
            vec![1.0],
            vec![2.0],
            vec![3.0],
            vec![4.0],
            vec![5.0],
            vec![6.0],
            vec![7.0],
            vec![8.0],
        ];
        let targets = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];

        let config = TreeConfig {
            max_depth: 5,
            min_samples_split: 2,
            min_samples_leaf: 1,
            max_features: 0,
        };
        let tree = build_tree(&features, &targets, &config, 42);

        // Should predict close to the training values
        let pred = tree.predict(&[4.0]);
        assert!((pred - 4.0).abs() < 2.0);
    }

    #[test]
    fn test_random_forest() {
        let features: Vec<Vec<f64>> = (0..20).map(|i| vec![i as f64]).collect();
        let targets: Vec<f64> = (0..20).map(|i| (i as f64) * 2.0 + 1.0).collect();

        let config = RandomForestConfig {
            n_trees: 10,
            seed: 42,
            ..Default::default()
        };
        let model = train_random_forest(&features, &targets, &config);
        assert_eq!(model.trees.len(), 10);

        let (pred, var) = model.predict_with_variance(&[10.0]);
        assert!(pred > 0.0);
        assert!(var >= 0.0);
    }

    /// The out-of-bag score must be a real number and must be close to the
    /// error the forest makes on unseen data (it used to be always `None`).
    #[test]
    fn test_random_forest_oob_score() {
        // A deterministic linear target: the forest should learn it well.
        let features: Vec<Vec<f64>> = (0..40).map(|i| vec![i as f64 * 0.1]).collect();
        let targets: Vec<f64> = features.iter().map(|f| 2.0 * f[0] + 1.0).collect();
        let config = RandomForestConfig {
            n_trees: 12,
            seed: 7,
            ..RandomForestConfig::default()
        };
        let rf = train_random_forest(&features, &targets, &config);
        let oob = rf.oob_score.expect("an OOB score should be available");
        assert!(oob.is_finite(), "OOB score must be finite: {oob}");
        assert!(
            oob < 1.0,
            "OOB error should be small for a linear target: {oob}"
        );
    }

    #[test]
    fn test_gradient_boosting() {
        let features: Vec<Vec<f64>> = (0..20).map(|i| vec![i as f64]).collect();
        let targets: Vec<f64> = (0..20).map(|i| (i as f64) * 2.0 + 1.0).collect();

        let config = GradientBoostingConfig {
            n_estimators: 20,
            learning_rate: 0.3,
            seed: 42,
            ..Default::default()
        };
        let model = train_gradient_boosting(&features, &targets, &config);
        let pred = model.predict(&[10.0]);
        assert!(pred > 0.0);
    }

    #[test]
    fn test_cross_validation() {
        let features: Vec<Vec<f64>> = (0..30).map(|i| vec![i as f64]).collect();
        let targets: Vec<f64> = (0..30).map(|i| (i as f64) * 2.0).collect();

        let config = RandomForestConfig {
            n_trees: 5,
            seed: 42,
            ..Default::default()
        };
        let cv = cross_validate(&features, &targets, &ModelType::RandomForest(config), 5, 42);
        assert_eq!(cv.k, 5);
        assert_eq!(cv.fold_mae.len(), 5);
    }

    #[test]
    fn test_isotonic_calibrator() {
        let preds = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let targets = vec![1.1, 2.5, 2.8, 4.2, 5.1];

        let cal = IsotonicCalibrator::fit(&preds, &targets);
        let calibrated = cal.calibrate(3.0);
        // Should return a value near 2.8 (the fitted value at x=3)
        assert!((calibrated - 2.8).abs() < 1.0);
    }
}
