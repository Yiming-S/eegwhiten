# Fit a whitening model for EEG data

The input matrix `X` is assumed to have trials/epochs in rows and
channels or features in columns. This function estimates a centering
vector and a whitening matrix, optionally using regularization
(shrinkage) and dimensionality reduction.

## Usage

``` r
whiten_model(
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

  Logical; whether to subtract column means before computing the
  covariance.

- method:

  Whitening method; one of `"SVD"`, `"ZCA"`, `"ZCA-cor"`, `"PCA"`,
  `"PCA-cor"`, `"Cholesky"`. Underscore aliases `"ZCA_cor"` and
  `"PCA_cor"` are also accepted.

- n_comp:

  Integer (optional); Number of components to keep. Only applicable for
  "PCA", "PCA-cor", and "SVD".

- var_threshold:

  Numeric in (0, 1\]; cumulative explained variance threshold for
  automatic component selection. Ignored by non-reduction methods.

- lambda:

  Numeric (0 to 1) or `"auto"`; regularization parameter for covariance
  shrinkage. 0 = No regularization (Empirical covariance). 1 = Full
  shrinkage to target. Values between 0 and 1 mix the two. Useful for
  high-dimensional/low-sample EEG data.

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

  Covariance shrinkage target; one of `"identity"` (shrink toward a
  scaled identity `(tr(S)/d) I`, the default) or `"diagonal"` (shrink
  correlations toward zero while preserving the per-channel variances
  `diag(S)`). `lambda = "auto"` requires `"identity"`.

- na_action:

  How to handle non-finite values in `X`; one of `"error"` (default) or
  `"omit"` (drop offending rows/trials with a warning).

## Value

An object of class `"whiten_model"`.

## See also

[`predict.whiten_model`](https://yiming-s.github.io/eegwhiten/reference/predict.whiten_model.md)
to apply the model,
[`unwhiten`](https://yiming-s.github.io/eegwhiten/reference/unwhiten.md)
for inverse transforms,
[`whiten_fit`](https://yiming-s.github.io/eegwhiten/reference/whiten_fit.md)
for the low-level API,
[`whiten_matrix`](https://yiming-s.github.io/eegwhiten/reference/whiten_matrix.md)
for one-shot whitening,
[`check_whitening`](https://yiming-s.github.io/eegwhiten/reference/check_whitening.md)
for diagnostics.

## Examples

``` r
set.seed(42)
X_train <- matrix(rnorm(300 * 16), 300, 16)
colnames(X_train) <- paste0("Ch", 1:16)

m <- whiten_model(X_train, method = "PCA", var_threshold = 0.95,
                  lambda = "auto")
Z <- predict(m, X_train)
print(m)
#> <whiten_model>
#>   method      : PCA 
#>   dimensions  : 16 -> 16 
#>   n_comp      : 16 
#>   var_threshold: 0.95 
#>   explained_var: 1 
#>   lambda      : 1 
#>   shrink_target: identity 
#>   eig_method  : auto 
#>   fast_mode   : FALSE 
#>   cov_estimator: empirical 
#>   weighted_fit: FALSE 
#>   n_obs       : 300 
#>   updates     : 0 
#>   lambda_mode :auto (oas)
```
