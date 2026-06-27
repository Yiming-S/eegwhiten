# Introduction to eegwhiten

## Overview

`eegwhiten` provides covariance whitening and covariance-based
transfer-learning tools for EEG and other multichannel signals. The core
is a model-based fit/transform workflow:

- **Rows** are trials/epochs/samples; **columns** are channels/features.
- Fit once on training data, then apply to new data.
- Six whitening methods, covariance shrinkage, robust estimators,
  dimensionality reduction, and cross-validated tuning.

On top of whitening it adds the pieces needed for cross-session /
cross-subject work: **Euclidean alignment** and **per-domain
recentering**, **Riemannian tangent-space** mapping, and
**generalized-eigenvalue (“relative”) whitening**.

``` r

library(eegwhiten)
```

## Quick start: fit, transform, invert

``` r

X_train <- matrix(rnorm(300 * 16), 300, 16)
X_test  <- matrix(rnorm(120 * 16), 120, 16)
colnames(X_train) <- paste0("Ch", 1:16)
colnames(X_test)  <- colnames(X_train)

m <- whiten_model(X_train, method = "PCA", var_threshold = 0.95, lambda = "auto")
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

Z_test <- predict(m, X_test)
dim(Z_test)
#> [1] 120  16

# Inverse transform (note: model-first argument order)
X_rec <- unwhiten(m, Z_test)
dim(X_rec)
#> [1] 120  16
```

## Diagnostics

After whitening, `cov(Z)` should be close to the identity.
[`check_whitening()`](https://yiming-s.github.io/eegwhiten/reference/check_whitening.md)
reports several measures, including a bounded `whiteness` score in
`[0, 1]`.

``` r

d <- check_whitening(predict(m, X_train))
c(diag_mean = d$diag_mean, offdiag_frob = d$offdiag_frob,
  whiteness = d$whiteness, logdet = d$logdet)
#>    diag_mean offdiag_frob    whiteness       logdet 
#>    1.0000000    0.8378828    0.7824632   -0.3792011
```

## Whitening methods

| Method   | Description                          | Dim. reduction |
|----------|--------------------------------------|:--------------:|
| PCA      | Eigen-decomposition based            |      Yes       |
| PCA-cor  | PCA on correlation (scale-invariant) |      Yes       |
| SVD      | Singular value decomposition         |      Yes       |
| ZCA      | Zero-phase, closest to original data |       No       |
| ZCA-cor  | ZCA on correlation                   |       No       |
| Cholesky | Cholesky factorization               |       No       |

## Regularization and shrinkage targets

Covariance shrinkage stabilizes whitening for high-dimensional /
low-sample data. `lambda = "auto"` chooses the intensity with OAS or
Ledoit-Wolf. The target is a scaled identity by default, or `"diagonal"`
to shrink correlations toward zero while preserving each channel’s
variance.

``` r

m_auto <- whiten_model(X_train, method = "ZCA", lambda = "auto", lambda_method = "oas")
m_auto$lambda
#> [1] 1

m_diag <- whiten_model(X_train, method = "ZCA", lambda = 0.3, shrink_target = "diagonal")
m_diag$shrink_target
#> [1] "diagonal"
```

Robust covariance estimators (`"tyler"`, and `"mcd"` with
**robustbase**) are available via `cov_estimator`.

## Spatial filters and patterns

The whitening filters (rows of `W`) and the corresponding forward
patterns are both available. Patterns – not raw filter weights – are the
right object to interpret as scalp topographies (Haufe et al., 2014).

``` r

fp <- whitening_patterns(whiten_model(X_train, method = "PCA", n_comp = 4))
dim(fp$filters)   # components x channels
#> [1]  4 16
dim(fp$patterns)  # channels x components
#> [1] 16  4
```

## Hyperparameter tuning

[`auto_tune_whitening()`](https://yiming-s.github.io/eegwhiten/reference/auto_tune_whitening.md)
cross-validates over methods and parameters and returns a tuned object
you can predict with directly.

``` r

tuned <- auto_tune_whitening(
  X_train, methods = c("PCA", "ZCA"),
  n_comp_grid = c(8, 12), lambda_grid = list("auto", 0.1),
  cv_folds = 3, seed = 123, top_n = 4
)
print(tuned)
#> <whiten_tune>
#>   scoring     : unsupervised 
#>   cv_folds    : 3 
#>   best method : PCA 
#>   best n_comp : 8 
#>   best lambda : 0.1 
#>   best cov_est: empirical 
#>   top 4 configurations:
#>  method n_comp var_threshold lambda cov_estimator mean_score  sd_score
#>     PCA      8            NA    0.1     empirical -0.8503463 0.1056367
#>     PCA      8            NA   auto     empirical -0.9261697 0.2005509
#>     PCA     12            NA    0.1     empirical -1.1336164 0.2216557
#>     PCA     NA            NA   auto     empirical -1.6602069 0.3187931
#>  n_success
#>          3
#>          3
#>          3
#>          3

Z <- predict(tuned, X_test)   # uses the best model
dim(Z)
#> [1] 120   8
```

## Cross-session / cross-subject transfer

Independent per-matrix whitening removes the very between-session
covariance structure transfer learning needs. `eegwhiten` offers
alignment instead.

**Per-domain recentering**
([`recenter()`](https://yiming-s.github.io/eegwhiten/reference/recenter.md))
aligns each session by its *own* covariance, mapping every domain to the
identity – the Euclidean Alignment idea of He and Wu (2019):

``` r

X_list <- list(
  matrix(rnorm(200 * 8), 200, 8),
  matrix(rnorm(180 * 8), 180, 8) %*% diag(1:8)  # a shifted session
)

out <- whiten_batch(X_list, mode = "recenter")
round(diag(cov(out[[2]]$Z)), 3)   # each domain -> identity
#> [1] 1 1 1 1 1 1 1 1
```

**Euclidean alignment**
([`euclidean_alignment()`](https://yiming-s.github.io/eegwhiten/reference/euclidean_alignment.md))
instead builds one shared transform from a grand-mean covariance:

``` r

ea <- euclidean_alignment(X_list, input = "raw", mean_method = "logeuclid")
str(ea$aligned, max.level = 1)
#> List of 2
#>  $ : num [1:200, 1:8] -1.196 -0.572 0.591 1.129 0.689 ...
#>  $ : num [1:180, 1:8] 0.566 -1.6467 1.492 -0.0236 -1.0221 ...
```

At predict time, `self_center = TRUE` centers test data on its own mean,
which helps when the test session has drifted relative to training:

``` r

m_tr <- whiten_model(X_list[[1]], method = "ZCA", lambda = 0.1)
Z_drift <- predict(m_tr, X_list[[2]], self_center = TRUE)
dim(Z_drift)
#> [1] 180   8
```

## Riemannian tangent space

The standard Riemannian BCI pipeline projects per-trial covariance
matrices to the tangent space at a reference mean and vectorizes them,
so ordinary linear classifiers can be used.
[`epoch_covariances()`](https://yiming-s.github.io/eegwhiten/reference/epoch_covariances.md)
builds the covariances from a trial tensor;
[`tangent_space()`](https://yiming-s.github.io/eegwhiten/reference/tangent_space.md)
maps them.

``` r

X_tensor <- array(rnorm(40 * 8 * 128), dim = c(40, 8, 128))
covs <- epoch_covariances(X_tensor, lambda = 0.05)
V <- tangent_space(covs, mean_method = "riemann")
dim(V)                       # 40 trials x p(p+1)/2 features
#> [1] 40 36

covs_rec <- untangent_space(V)   # inverse map
max(abs(covs_rec[[1]] - covs[[1]]))
#> [1] 4.440892e-15
```

## Relative / generalized-eigenvalue whitening

[`whiten_relative()`](https://yiming-s.github.io/eegwhiten/reference/whiten_relative.md)
whitens one covariance by a reference covariance and diagonalizes the
target – the building block of CSP and xDAWN. The returned `values` are
the per-component variance ratios.

``` r

Sigma_signal <- cov(matrix(rnorm(200 * 6), 200, 6) %*% diag(c(4, 3, 2, 1, 1, 1)))
Sigma_noise  <- cov(matrix(rnorm(200 * 6), 200, 6))
gr <- whiten_relative(Sigma_signal, reference = Sigma_noise, n_comp = 3)
round(gr$values, 3)
#> [1] 14.154  8.414  3.884
```

## Online updates with forgetting

[`whiten_model_update()`](https://yiming-s.github.io/eegwhiten/reference/whiten_model_update.md)
refreshes a model from summary moments. A `decay < 1` down-weights past
data to track non-stationary drift in streaming settings.

``` r

X_new <- matrix(rnorm(80 * 16), 80, 16)
m2 <- whiten_model_update(m, X_new, decay = 0.9)
m2$n_obs
#> [1] 380
```

## Tensor data

3D EEG tensors (trials x channels x time) are handled directly.

``` r

m_t <- whiten_model_tensor(X_tensor, method = "ZCA", lambda = 0)
Z_t <- predict_tensor(m_t, X_tensor)
dim(Z_t)
#> [1]  40   8 128
X_t_rec <- unwhiten_tensor(m_t, Z_t)
dim(X_t_rec)
#> [1]  40   8 128
```

## Reporting

``` r

txt <- report_whitening(m, data = X_train, file = NULL)
cat(substr(txt, 1, 280), "...\n")
#> # Whitening Report
#> 
#> - Generated at: 2026-06-27 20:19:36 UTC
#> - Method: PCA
#> - Dimensions: 16 -> 16
#> - n_comp: 16
#> - var_threshold: 0.95
#> - explained_var: 1
#> - lambda: 1
#> - lambda_input: auto
#> - lambda_method: oas
#> - shrink_target: identity
#> - cov_estimator: empirical
#> - sample_weighted: FAL ...
```

## References

- He, H., and Wu, D. (2019). Transfer learning for brain-computer
  interfaces: A Euclidean space data alignment approach. *IEEE TBME*.
- Barachant, A., et al. (2012). Multiclass brain-computer interface
  classification by Riemannian geometry. *IEEE TBME*.
- Haufe, S., et al. (2014). On the interpretation of weight vectors of
  linear models in multivariate neuroimaging. *NeuroImage*.
