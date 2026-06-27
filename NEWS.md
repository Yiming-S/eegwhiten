# eegwhiten 0.3.0

## New features

* **Riemannian tangent space mapping.** `tangent_space()` and
  `untangent_space()` project SPD covariance matrices to/from the tangent space
  at a reference mean (Riemann / log-Euclidean / Euclidean), with an isometric
  vectorization. `epoch_covariances()` builds per-trial covariances from an EEG
  tensor as input.
* **Per-domain recentering.** `recenter()` aligns a matrix by its own
  covariance (true Euclidean Alignment semantics), and `whiten_batch()` gains
  `mode = "recenter"`. `predict.whiten_model()` gains `recenter = TRUE` to
  center new data on its own mean for cross-session / cross-subject transfer.
* **Relative (generalized-eigenvalue) whitening.** `whiten_relative()` whitens
  one covariance by a reference covariance and diagonalizes the target — the
  CSP / xDAWN building block — returning filters, variance ratios, and patterns.
* **Online updates with forgetting.** `whiten_model_update()` gains a `decay`
  argument for exponential down-weighting of past data to track non-stationary
  drift.
* **Diagonal shrinkage target.** `whiten_model()` (and wrappers) gain
  `shrink_target = "diagonal"` to shrink correlations toward zero while
  preserving per-channel variances.
* **Non-finite handling.** `whiten_model()` (and wrappers) gain
  `na_action = "omit"` to drop trials containing non-finite values.
* **Spatial filters and patterns.** `whitening_patterns()` returns the
  whitening filters and the corresponding forward patterns (Haufe et al., 2014).
* **Richer diagnostics.** `check_whitening()` now also reports `cov_dev_frob`,
  `logdet`, and a bounded `whiteness` score.

## Performance

* **Single eigendecomposition per fit.** Covariance-based methods
  (`PCA`, `SVD`, `ZCA`) no longer eigendecompose the covariance twice; the
  decomposition is computed once and reused.
* **Analytic shrinkage path in tuning.** `auto_tune_whitening(fast_tune = TRUE)`
  eigendecomposes each fold's covariance once and reuses it across the entire
  `lambda` / `n_comp` / `var_threshold` grid (identity target, covariance-based
  methods), which is substantially faster and numerically equivalent to a full
  refit per candidate.

## Bug fixes

* `auto_tune_whitening()` no longer errors when `var_threshold_grid` is
  supplied.
