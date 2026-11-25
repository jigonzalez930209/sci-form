# Uncertainty-aware active learning

The first v2.2 machine-learning prototype adds a conformal interval around a
regression prediction and a query-by-committee acquisition score. It is
model-agnostic: any model that exposes a prediction and an uncertainty estimate
can use it. The existing Random Forest already exposes both through
`predict_with_variance`.

> **Status.** Implemented in `feat/wave-b`. Rust, Python, WASM and CLI
> (`conformal`) expose the calibration; model-bound bindings are pending. See the
> [v2.x roadmap](/v2-roadmap).
>
> **Feature.** Core since 2.2.0; `alpha-active-learning` remains as an empty
> compatibility alias.

## Why conformal calibration

An ensemble standard deviation is not a coverage guarantee. Conformal prediction
turns any point predictor into an interval with a finite-sample coverage
statement under exchangeability:

$$
\hat{C}(x) = \left[\hat{\mu}(x)-q,\ \hat{\mu}(x)+q\right],
$$

where `q` is the empirical quantile of the absolute residuals on a held-out
calibration set. The implementation uses the finite-sample order statistic

$$
q = \operatorname{Quantile}_{1-\alpha}
\left(\left\{|y_i-\hat{\mu}(x_i)|\right\}\right).
$$

The output reports the requested coverage and the number of calibration points,
so a user can see whether the interval is meaningful.

## Acquisition

Given a candidate pool, the session ranks points by one of:

| Acquisition | Score |
|---|---|
| `Uncertainty` | Model standard deviation |
| `Distance` | Euclidean distance to the nearest labelled point |
| `Hybrid { w }` | `w · uncertainty + (1-w) · distance` |

Distance is a simple domain signal, not a chemistry-aware similarity. It is a
placeholder for fingerprints or learned embeddings.

## API

```rust
use sci_form::active_learning::{
    Acquisition, ActiveLearningSession,
};
use sci_form::ml::advanced_models::{
    train_random_forest, RandomForestConfig,
};

let model = train_random_forest(&features, &targets, &RandomForestConfig::default());
let mut session = ActiveLearningSession::new(model, Acquisition::Hybrid {
    uncertainty_weight: 0.7,
});

for (features, target) in labelled {
    session.add_label(features, target);
}

let calibrator = session.calibrate(0.9)?;
println!("conformal radius = {:.3}", calibrator.radius);

let candidates = vec![vec![0.1, 0.2], vec![4.0, 1.0]];
let next = session.propose(&candidates, 1);
let interval = session.predict(&candidates[next[0]]);
println!("{:.3} ± {:.3} (in domain: {})",
    interval.mean, interval.std, interval.in_domain);
```

### Types

| Type | Meaning |
|---|---|
| `ConformalCalibrator` | Radius, requested coverage and calibration size |
| `PredictionInterval` | Mean, standard deviation, bounds, coverage and domain flag |
| `Acquisition` | Uncertainty, distance or hybrid score |
| `ActiveLearningSession` | Model, labelled points, calibrator and acquisition |

Predictions, targets and the conformal radius share the target unit. Coverage
is dimensionless.

## Random Forest integration

`RandomForest` implements the `UncertaintyModel` trait:

- `predict` returns the forest mean;
- `uncertainty` returns the square root of the tree-variance.

The variance is the spread of the individual tree predictions, not a Bayesian
posterior. It is useful for ranking, but the conformal radius is the quantity
with a coverage interpretation.

## Validation

| Gate | Observable | Result |
|---|---|---|
| Conformal quantile | 100 residuals of 1.0, coverage 0.9 | radius = 1.0 |
| Interval | Calibrated session | bounds contain the mean, coverage recorded |
| Acquisition | Candidate with larger uncertainty | ranked first |
| Invalid coverage | 0, 1, NaN | rejected |
| Domain coverage | Distance acquisition vs 100 random selections | gap ratio 0.398 (F135) |

A real coverage test needs a held-out chemical dataset and exchangeable
calibration/test splits. F135 validates the acquisition loop on a synthetic
1D oracle; a chemical benchmark is still pending.

## Bindings

| Surface | Entry point |
|---|---|
| CLI | `sci-form conformal '<predictions-json>' '<targets-json>' --coverage 0.9` |
| CLI | `sci-form active-learning '<ActiveLearningRequest-json>'` |
| WASM | `conformal_calibrate(predictions_json, targets_json, coverage)` |
| WASM | `active_learning_propose(request_json)` |
| Python | `conformal_calibrate(predictions, targets, coverage)` |
| Python | `active_learning_propose(request_json)` |

`ActiveLearningRequest` carries a serialized Random Forest, the labelled
points, the candidate features, the acquisition and the number of labels. The
generic `ActiveLearningSession<M: UncertaintyModel>` remains the Rust API for
custom models.
## Limitations

- Conformal coverage assumes exchangeability. A new chemical series can break
  it even when the marginal coverage looked correct.
- The distance acquisition is a raw Euclidean metric in feature space; it is
  not a chemistry-aware novelty detector.
- The session does not retrain the model automatically. The caller trains the
  next model and creates a new session.
- The interval is marginal, not conditional on a specific feature region.
- The Random Forest tree variance can be zero when all trees agree; use the
  conformal interval, not the raw standard deviation, for coverage.
- The JSON surface uses a serialized Random Forest; custom models stay in Rust.

## Next steps

1. Add a public calibration/coverage benchmark on a chemical dataset.
2. Add fingerprint or learned-embedding distances.
3. Add a retraining loop with a trainer callback.
4. Add Mondrian/conditional conformal intervals.
5. Add a retraining loop benchmark on a chemical dataset.
6. Connect the acquisition score to the universal MLIP label budget.
