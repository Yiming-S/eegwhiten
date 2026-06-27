# Batch whitening for a list of EEG matrices

Each matrix in `X_list` is whitened independently using `whiten_matrix`.
For proper train/test usage, fit a model on training data with
[`whiten_model()`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md)
and call [`predict()`](https://rdrr.io/r/stats/predict.html) on test
data instead of using this function.

## Usage

``` r
whiten_batch(
  X_list,
  center = TRUE,
  method = c("SVD", "ZCA", "ZCA-cor", "PCA", "PCA-cor", "Cholesky"),
  n_comp = NULL,
  var_threshold = NULL,
  lambda = "auto",
  lambda_method = c("oas", "lw"),
  sample_weight = NULL,
  cov_estimator = c("empirical", "mcd", "tyler"),
  tyler_tol = 1e-06,
  tyler_max_iter = 200,
  tyler_eps = 1e-06,
  eig_method = c("auto", "base", "rspectra"),
  fast = FALSE,
  sign_reference = NULL,
  shared_model = FALSE,
  mode = c("independent", "shared_model", "ea", "recenter"),
  ea_mean = c("logeuclid", "riemann", "euclid"),
  ea_input = c("auto", "raw", "cov"),
  ea_tol = 1e-06,
  ea_max_iter = 50,
  parallel = FALSE,
  n_cores = 1L,
  shrink_target = c("identity", "diagonal"),
  na_action = c("error", "omit")
)
```

## Arguments

- X_list:

  List of numeric matrices, each \[n_trials x n_channels\].

- center:

  Logical; whether to subtract the column means in each matrix.

- method:

  Whitening method; one of `"SVD"`, `"ZCA"`, `"ZCA-cor"`, `"PCA"`,
  `"PCA-cor"`, `"Cholesky"`. Underscore aliases `"ZCA_cor"` and
  `"PCA_cor"` are also accepted.

- n_comp:

  Integer (optional); number of components to keep for `"PCA"`,
  `"PCA-cor"`, and `"SVD"`.

- var_threshold:

  Numeric in (0, 1\]; cumulative explained-variance threshold for
  automatic `n_comp` selection.

- lambda:

  Numeric in \[0, 1\] or `"auto"`; covariance shrinkage strength.

- lambda_method:

  Method used when `lambda = "auto"`. One of `"oas"` or `"lw"`.

- sample_weight:

  Optional sample weights. For a single matrix, provide one numeric
  vector. For multiple matrices in `mode = "independent"`, provide a
  list with one vector per matrix.

- cov_estimator:

  Covariance estimator; one of `"empirical"`, `"mcd"`, or `"tyler"`.

- tyler_tol:

  Convergence tolerance for `cov_estimator = "tyler"`.

- tyler_max_iter:

  Maximum iterations for Tyler covariance estimation.

- tyler_eps:

  Numerical stabilization floor for Tyler estimation.

- eig_method:

  Eigen solver backend; one of `"auto"`, `"base"`, or `"rspectra"`.

- fast:

  Logical; if `TRUE`, allow faster approximate settings for iterative
  eigensolvers.

- sign_reference:

  Optional reference vectors used for stable component sign orientation
  in `"PCA"`, `"PCA-cor"`, and `"SVD"`.

- shared_model:

  Logical; if `TRUE`, fit one model on the first matrix and apply it to
  all matrices in `X_list`. If `FALSE`, each matrix is whitened
  independently.

- mode:

  Batch whitening strategy: `"independent"` (fit each matrix
  separately), `"shared_model"` (fit once on first matrix), `"ea"`
  (Euclidean Alignment with a single global reference covariance), or
  `"recenter"` (per-domain recentering: align each matrix by its own
  covariance so every domain is whitened to the identity; see
  [`recenter`](https://yiming-s.github.io/eegwhiten/reference/recenter.md)).

- ea_mean:

  Mean type used in `"ea"` mode; one of `"riemann"`, `"logeuclid"`, or
  `"euclid"`.

- ea_input:

  Input interpretation in `"ea"` mode; one of `"auto"`, `"raw"`, or
  `"cov"`.

- ea_tol:

  Tolerance for Riemannian mean convergence in `"ea"` mode.

- ea_max_iter:

  Maximum iterations for Riemannian mean in `"ea"` mode.

- parallel:

  Logical; if `TRUE`, use multicore processing for per-matrix operations
  on Unix-like systems.

- n_cores:

  Integer; number of worker processes when `parallel = TRUE`.

- shrink_target:

  Covariance shrinkage target; one of `"identity"` or `"diagonal"`. See
  [`whiten_model`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md).

- na_action:

  How to handle non-finite values; one of `"error"` or `"omit"` (drop
  offending rows). Applied for the `"independent"`, `"shared_model"`,
  and `"recenter"` modes; for `mode = "ea"` pre-clean the inputs
  instead.

## Value

A list of results as returned by
[`whiten_matrix()`](https://yiming-s.github.io/eegwhiten/reference/whiten_matrix.md).

## See also

[`whiten_matrix`](https://yiming-s.github.io/eegwhiten/reference/whiten_matrix.md),
[`whiten_model`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md),
[`euclidean_alignment`](https://yiming-s.github.io/eegwhiten/reference/euclidean_alignment.md)
for direct EA access.

## Examples

``` r
set.seed(1)
X_list <- list(
  matrix(rnorm(200 * 8), 200, 8),
  matrix(rnorm(180 * 8), 180, 8)
)
out <- whiten_batch(X_list, mode = "shared_model",
                    method = "PCA", n_comp = 4, lambda = 0.1)
length(out)
#> [1] 2
```
