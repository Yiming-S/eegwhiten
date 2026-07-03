# eegwhiten 1.0.0

First stable release. The public API is now considered stable: the six whitening
transforms, the `whiten_model()` fit/`predict()` workflow (including dimension
reduction, shrinkage, and robust covariance estimators), the covariance-alignment
tools (Euclidean alignment, recentering, barycenter whitening), the Riemannian
tangent-space helpers, and the diagnostics/reporting functions. No breaking
changes from 0.5.2; this release marks the transition out of the experimental
lifecycle.

# eegwhiten 0.5.2

* Added a gallery of demonstration figures to the README and pkgdown site
  (whitening before/after, the six transforms, tangent-space alignment, and
  covariance-estimator comparisons), plus live illustrations in the vignette.
  All figures are reproducible with `tools/make-figures.R`.

# eegwhiten 0.5.1

## Bug fixes

* `whiten_batch(mode = "ea", center = TRUE)` now reports the per-matrix
  centering vector actually subtracted (the column means) instead of a zero
  vector, so reconstructing from `(Z, W, center)` recovers the original data.
  `euclidean_alignment()` now returns these vectors in a `centers` list.
* A model fitted with `lambda = "auto"` and `shrink_target = "diagonal"` now
  reports `lambda_method = "ss"` (Schaefer-Strimmer, the estimator actually
  used) instead of the ignored `"oas"`/`"lw"` input, in the object and in
  `print()` / `summary()`.

# eegwhiten 0.5.0

## New features

* `riemann_distance()` and `spd_mean()` expose the affine-invariant Riemannian
  (and log-Euclidean / Euclidean) distance and mean of SPD covariance matrices.
* `barycenter_whitener()` fits a whitener that maps the geometric mean of a set
  of covariances to the identity (align to the Riemannian barycenter).
* `lambda = "auto"` now works with `shrink_target = "diagonal"`, using the
  scale-invariant Schaefer-Strimmer shrinkage intensity (verified against
  `corpcor::estimate.lambda`).
* `epoch_covariances()` accepts `lambda = "auto"` for a per-epoch Ledoit-Wolf
  shrinkage intensity.

# eegwhiten 0.4.1

* `predict()` now works directly on an `auto_tune_whitening()` result, applying
  the selected best model.
* The scale-invariant positive-definiteness fix now also covers `whiten_fit()`
  and the standalone `ZCA_cor()` / `PCA_cor()` functions (correlation-based
  methods check the correlation matrix), and the relative tolerance is a single
  named constant.
* The introductory vignette now covers the alignment, Riemannian tangent-space,
  relative-whitening, shrinkage-target, and online-update features.
* Refreshed the package description; removed dead internal helpers.

# eegwhiten 0.4.0

## Correctness

* `cov_estimator = "tyler"` now restores the covariance scale (Tyler's iteration
  estimates only the trace-normalized shape), so whitening with it produces
  approximately unit-variance data instead of decorrelated-but-mis-scaled data.
* The positive-definiteness check is now scale-invariant (relative to the
  largest eigenvalue), so well-conditioned data at any magnitude -- e.g. EEG in
  volts -- is no longer falsely rejected. Correlation-based methods test
  positive-definiteness on the correlation matrix.

## API consistency

* Functions taking a model and data now use a consistent model-first argument
  order: `unwhiten(model, Z)`, `unwhiten_fast(model, Z)`,
  `unwhiten_tensor(model, Z_tensor)`, and `inverse_ea(model, Z)`. The old
  data-first order still works with a deprecation warning.
* `epoch_covariances()` shrinkage argument is renamed `lambda` (matching the
  rest of the package); `shrinkage` is deprecated.
* `predict.whiten_model()` argument `recenter` is renamed `self_center` to avoid
  clashing with the `recenter()` function; `recenter` is deprecated.
* `recenter()`, `whiten_relative()`, the tuning result (`whiten_tune`), and
  `ea_model` now have `print()` methods (and `predict()` where applicable),
  instead of returning bare lists.

## Packaging

* Added pkgdown configuration, a coverage workflow, R-devel/oldrel CI, a
  `CITATION`, README badges, `cran-comments.md`, `CONTRIBUTING.md`, and
  `CODE_OF_CONDUCT.md`.

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
