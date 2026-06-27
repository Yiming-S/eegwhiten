# Automatically Tune Whitening Hyperparameters

Performs cross-validation across candidate whitening configurations and
returns the best model fit on full data.

## Usage

``` r
auto_tune_whitening(
  X,
  y = NULL,
  methods = c("SVD", "ZCA", "ZCA-cor", "PCA", "PCA-cor", "Cholesky"),
  n_comp_grid = NULL,
  var_threshold_grid = NULL,
  lambda_grid = list("auto", 0, 0.05, 0.1),
  lambda_method = c("oas", "lw"),
  shrink_target = c("identity", "diagonal"),
  cov_estimator_grid = "empirical",
  sample_weight = NULL,
  cv_folds = 5L,
  seed = NULL,
  scoring = c("auto", "unsupervised", "accuracy"),
  score_fn = NULL,
  tyler_tol = 1e-06,
  tyler_max_iter = 200,
  tyler_eps = 1e-06,
  eig_method = c("auto", "base", "rspectra"),
  fast = FALSE,
  sign_reference = NULL,
  top_n = 5L,
  fast_tune = TRUE
)
```

## Arguments

- X:

  Numeric matrix `[n_samples x n_features]`.

- y:

  Optional labels for supervised tuning.

- methods:

  Candidate whitening methods.

- n_comp_grid:

  Optional integer grid for `n_comp`.

- var_threshold_grid:

  Optional numeric grid in `(0, 1]` for `var_threshold`.

- lambda_grid:

  Candidate shrinkage values. Can include numeric values in `[0, 1]` and
  `"auto"`.

- lambda_method:

  Method used when `lambda = "auto"`.

- shrink_target:

  Covariance shrinkage target passed to
  [`whiten_model`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md);
  one of `"identity"` or `"diagonal"`.

- cov_estimator_grid:

  Candidate covariance estimators.

- sample_weight:

  Optional sample weights (used in training folds).

- cv_folds:

  Number of folds.

- seed:

  Optional random seed.

- scoring:

  Scoring mode: `"auto"`, `"unsupervised"`, or `"accuracy"`.

- score_fn:

  Optional custom scoring function. For unsupervised tuning, called as
  `score_fn(Z_valid)`. For supervised tuning, called as
  `score_fn(Z_train, y_train, Z_valid, y_valid)`.

- tyler_tol:

  Convergence tolerance for Tyler covariance.

- tyler_max_iter:

  Maximum iterations for Tyler covariance.

- tyler_eps:

  Numerical stabilization floor for Tyler covariance.

- eig_method:

  Eigen solver backend.

- fast:

  Logical; fast mode for eigensolver.

- sign_reference:

  Optional reference vectors for sign stabilization.

- top_n:

  Number of top configurations retained in the summary table.

- fast_tune:

  Logical; if `TRUE` (default), use an analytic shrinkage path for
  covariance-based methods (`"PCA"`, `"SVD"`, `"ZCA"`) with
  `shrink_target = "identity"`, no `sample_weight`, and no
  `sign_reference`. Each fold's covariance is eigendecomposed once (with
  an exact base solver, so `eig_method` and `fast` are not used for
  these candidates) and reused across the entire
  `lambda`/`n_comp`/`var_threshold` grid, which can be dramatically
  faster. Set to `FALSE` to force a full model refit per candidate.
  Results are equivalent up to floating point; candidates whose resolved
  `lambda` reaches 1 (an isotropic, degenerate target) are automatically
  routed through the full refit so scoring and the returned model stay
  consistent.

## Value

A list of class `whiten_tune` containing the best model, best
parameters, and cross-validation ranking.

## See also

[`whiten_model`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md),
[`report_whitening`](https://yiming-s.github.io/eegwhiten/reference/report_whitening.md)

## Examples

``` r
set.seed(42)
X <- matrix(rnorm(150 * 6), 150, 6)
tuned <- auto_tune_whitening(
  X,
  methods = c("PCA", "ZCA"),
  n_comp_grid = c(3),
  lambda_grid = list(0, 0.1),
  cv_folds = 3, seed = 1, top_n = 4
)
head(tuned$ranking)
#>   method n_comp var_threshold lambda cov_estimator mean_score   sd_score
#> 1    PCA      3            NA    0.1     empirical -0.4338521 0.10042868
#> 2    PCA      3            NA      0     empirical -0.4416311 0.09964211
#> 3    PCA     NA            NA    0.1     empirical -0.7476605 0.09974017
#> 4    PCA     NA            NA      0     empirical -0.7568348 0.09621823
#>   n_success
#> 1         3
#> 2         3
#> 3         3
#> 4         3
```
