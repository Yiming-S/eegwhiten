# One-shot whitening for a single EEG matrix

Convenience wrapper around
[`whiten_model()`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md)
that fits a whitening model on `X` and immediately applies it to `X`,
returning both the whitened data and the underlying model parameters.

## Usage

``` r
whiten_matrix(
  X,
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
  shrink_target = c("identity", "diagonal"),
  na_action = c("error", "omit")
)
```

## Arguments

- X:

  Numeric matrix \[n_trials x n_channels\].

- center:

  Logical; whether to subtract the column means.

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

  Optional non-negative sample weights with length `nrow(X)`.

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

- shrink_target:

  Covariance shrinkage target; one of `"identity"` or `"diagonal"`. See
  [`whiten_model`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md).

- na_action:

  How to handle non-finite values in `X`; one of `"error"` or `"omit"`.

## Value

A list with components:

- Z:

  Whitened data matrix.

- W:

  Whitening matrix.

- center:

  Centering vector used.

- method:

  Method name.

- n_comp:

  Number of retained components (if applicable).

- var_threshold:

  Variance threshold used for automatic selection.

- explained_var:

  Cumulative explained variance of retained components.

- lambda:

  Shrinkage value used during fitting.

- lambda_method:

  Automatic shrinkage estimator, if used.

- model:

  The fitted `whiten_model` object.

## See also

[`whiten_model`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md)
for the fit/transform workflow,
[`whiten_batch`](https://yiming-s.github.io/eegwhiten/reference/whiten_batch.md)
for batch whitening.

## Examples

``` r
set.seed(42)
X <- matrix(rnorm(200 * 16), 200, 16)
colnames(X) <- paste0("Ch", seq_len(16))
res <- whiten_matrix(X, center = TRUE, method = "PCA")
str(res$Z)
#>  num [1:200, 1:16] -1.1185 0.4954 0.0181 3.4823 -1.733 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr [1:16] "PC1" "PC2" "PC3" "PC4" ...
#>  - attr(*, "method")= chr "PCA"
```
