# Changelog

## eegwhiten 0.4.1

- [`predict()`](https://rdrr.io/r/stats/predict.html) now works directly
  on an
  [`auto_tune_whitening()`](https://yiming-s.github.io/eegwhiten/reference/auto_tune_whitening.md)
  result, applying the selected best model.
- The scale-invariant positive-definiteness fix now also covers
  [`whiten_fit()`](https://yiming-s.github.io/eegwhiten/reference/whiten_fit.md)
  and the standalone
  [`ZCA_cor()`](https://yiming-s.github.io/eegwhiten/reference/ZCA_cor.md)
  /
  [`PCA_cor()`](https://yiming-s.github.io/eegwhiten/reference/PCA_cor.md)
  functions (correlation-based methods check the correlation matrix),
  and the relative tolerance is a single named constant.
- The introductory vignette now covers the alignment, Riemannian
  tangent-space, relative-whitening, shrinkage-target, and online-update
  features.
- Refreshed the package description; removed dead internal helpers.

## eegwhiten 0.4.0

### Correctness

- `cov_estimator = "tyler"` now restores the covariance scale (Tyler’s
  iteration estimates only the trace-normalized shape), so whitening
  with it produces approximately unit-variance data instead of
  decorrelated-but-mis-scaled data.
- The positive-definiteness check is now scale-invariant (relative to
  the largest eigenvalue), so well-conditioned data at any magnitude –
  e.g. EEG in volts – is no longer falsely rejected. Correlation-based
  methods test positive-definiteness on the correlation matrix.

### API consistency

- Functions taking a model and data now use a consistent model-first
  argument order: `unwhiten(model, Z)`, `unwhiten_fast(model, Z)`,
  `unwhiten_tensor(model, Z_tensor)`, and `inverse_ea(model, Z)`. The
  old data-first order still works with a deprecation warning.
- [`epoch_covariances()`](https://yiming-s.github.io/eegwhiten/reference/epoch_covariances.md)
  shrinkage argument is renamed `lambda` (matching the rest of the
  package); `shrinkage` is deprecated.
- [`predict.whiten_model()`](https://yiming-s.github.io/eegwhiten/reference/predict.whiten_model.md)
  argument `recenter` is renamed `self_center` to avoid clashing with
  the
  [`recenter()`](https://yiming-s.github.io/eegwhiten/reference/recenter.md)
  function; `recenter` is deprecated.
- [`recenter()`](https://yiming-s.github.io/eegwhiten/reference/recenter.md),
  [`whiten_relative()`](https://yiming-s.github.io/eegwhiten/reference/whiten_relative.md),
  the tuning result (`whiten_tune`), and `ea_model` now have
  [`print()`](https://rdrr.io/r/base/print.html) methods (and
  [`predict()`](https://rdrr.io/r/stats/predict.html) where applicable),
  instead of returning bare lists.

### Packaging

- Added pkgdown configuration, a coverage workflow, R-devel/oldrel CI, a
  `CITATION`, README badges, `cran-comments.md`, `CONTRIBUTING.md`, and
  `CODE_OF_CONDUCT.md`.

## eegwhiten 0.3.0

### New features

- **Riemannian tangent space mapping.**
  [`tangent_space()`](https://yiming-s.github.io/eegwhiten/reference/tangent_space.md)
  and
  [`untangent_space()`](https://yiming-s.github.io/eegwhiten/reference/untangent_space.md)
  project SPD covariance matrices to/from the tangent space at a
  reference mean (Riemann / log-Euclidean / Euclidean), with an
  isometric vectorization.
  [`epoch_covariances()`](https://yiming-s.github.io/eegwhiten/reference/epoch_covariances.md)
  builds per-trial covariances from an EEG tensor as input.
- **Per-domain recentering.**
  [`recenter()`](https://yiming-s.github.io/eegwhiten/reference/recenter.md)
  aligns a matrix by its own covariance (true Euclidean Alignment
  semantics), and
  [`whiten_batch()`](https://yiming-s.github.io/eegwhiten/reference/whiten_batch.md)
  gains `mode = "recenter"`.
  [`predict.whiten_model()`](https://yiming-s.github.io/eegwhiten/reference/predict.whiten_model.md)
  gains `recenter = TRUE` to center new data on its own mean for
  cross-session / cross-subject transfer.
- **Relative (generalized-eigenvalue) whitening.**
  [`whiten_relative()`](https://yiming-s.github.io/eegwhiten/reference/whiten_relative.md)
  whitens one covariance by a reference covariance and diagonalizes the
  target — the CSP / xDAWN building block — returning filters, variance
  ratios, and patterns.
- **Online updates with forgetting.**
  [`whiten_model_update()`](https://yiming-s.github.io/eegwhiten/reference/whiten_model_update.md)
  gains a `decay` argument for exponential down-weighting of past data
  to track non-stationary drift.
- **Diagonal shrinkage target.**
  [`whiten_model()`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md)
  (and wrappers) gain `shrink_target = "diagonal"` to shrink
  correlations toward zero while preserving per-channel variances.
- **Non-finite handling.**
  [`whiten_model()`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md)
  (and wrappers) gain `na_action = "omit"` to drop trials containing
  non-finite values.
- **Spatial filters and patterns.**
  [`whitening_patterns()`](https://yiming-s.github.io/eegwhiten/reference/whitening_patterns.md)
  returns the whitening filters and the corresponding forward patterns
  (Haufe et al., 2014).
- **Richer diagnostics.**
  [`check_whitening()`](https://yiming-s.github.io/eegwhiten/reference/check_whitening.md)
  now also reports `cov_dev_frob`, `logdet`, and a bounded `whiteness`
  score.

### Performance

- **Single eigendecomposition per fit.** Covariance-based methods
  (`PCA`, `SVD`, `ZCA`) no longer eigendecompose the covariance twice;
  the decomposition is computed once and reused.
- **Analytic shrinkage path in tuning.**
  `auto_tune_whitening(fast_tune = TRUE)` eigendecomposes each fold’s
  covariance once and reuses it across the entire `lambda` / `n_comp` /
  `var_threshold` grid (identity target, covariance-based methods),
  which is substantially faster and numerically equivalent to a full
  refit per candidate.

### Bug fixes

- [`auto_tune_whitening()`](https://yiming-s.github.io/eegwhiten/reference/auto_tune_whitening.md)
  no longer errors when `var_threshold_grid` is supplied.
