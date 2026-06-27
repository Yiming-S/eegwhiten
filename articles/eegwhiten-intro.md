# Introduction to eegwhiten

## Overview

`eegwhiten` provides whitening transforms for EEG and multichannel
signal data. It supports six methods (PCA, PCA-cor, ZCA, ZCA-cor, SVD,
Cholesky) with a model-based API: fit once on training data, then apply
to new data.

## Quick Start

``` r

library(eegwhiten)
set.seed(42)

# Simulate EEG data: 300 trials x 16 channels
X_train <- matrix(rnorm(300 * 16), 300, 16)
X_test  <- matrix(rnorm(120 * 16), 120, 16)
colnames(X_train) <- paste0("Ch", 1:16)
colnames(X_test)  <- colnames(X_train)

# Fit a whitening model on training data
m <- whiten_model(
  X_train,
  method       = "PCA",
  var_threshold = 0.95,
  lambda       = "auto"
)
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

# Apply to test data
Z_test <- predict(m, X_test)
dim(Z_test)
#> [1] 120  16
```

## Diagnostics

After whitening, the covariance of the transformed data should be close
to the identity matrix. Use
[`check_whitening()`](https://yiming-s.github.io/eegwhiten/reference/check_whitening.md)
to verify:

``` r

d <- check_whitening(Z_test)
d$diag_mean       # should be close to 1
#> [1] 1.006482
d$offdiag_frob    # should be close to 0
#> [1] 1.491929
```

## Inverse Transform

Recover an approximation of the original data from whitened data:

``` r

X_rec <- unwhiten(Z_test, m)
#> Warning: unwhiten(data, model) is deprecated; use unwhiten(model, data).
dim(X_rec)
#> [1] 120  16
```

For full-rank methods (no dimensionality reduction), the reconstruction
is exact.

## One-Shot Whitening

For simple workflows where you don’t need the model object:

``` r

res <- whiten_matrix(X_train, method = "ZCA", lambda = 0.1)
str(res$Z)
#>  num [1:300, 1:16] 1.471 -0.412 0.359 0.642 0.457 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr [1:16] "L1" "L2" "L3" "L4" ...
#>  - attr(*, "method")= chr "ZCA"
```

## Whitening Methods

| Method   | Description                          | Dim. reduction |
|----------|--------------------------------------|:--------------:|
| PCA      | Eigen-decomposition based            |      Yes       |
| PCA-cor  | PCA on correlation (scale-invariant) |      Yes       |
| SVD      | Singular value decomposition         |      Yes       |
| ZCA      | Zero-phase, closest to original data |       No       |
| ZCA-cor  | ZCA on correlation                   |       No       |
| Cholesky | Cholesky factorization               |       No       |

## Regularization

Covariance shrinkage toward a scaled identity target prevents
ill-conditioning: - `lambda = 0`: no shrinkage - `lambda = "auto"`:
automatic selection via OAS or Ledoit-Wolf - `lambda` between 0 and 1:
manual shrinkage

``` r

m_reg <- whiten_model(X_train, method = "PCA", lambda = 0.2)
m_auto <- whiten_model(X_train, method = "PCA", lambda = "auto",
                       lambda_method = "oas")
m_auto$lambda
#> [1] 1
```

## Batch Whitening

Process multiple sessions or subjects at once:

``` r

X_list <- list(
  matrix(rnorm(200 * 8), 200, 8),
  matrix(rnorm(180 * 8), 180, 8)
)

# Shared model from first matrix
out <- whiten_batch(X_list, mode = "shared_model",
                    method = "PCA", n_comp = 4, lambda = 0.1)
length(out)
#> [1] 2
```

## Euclidean Alignment

Align multiple sessions to a global reference covariance:

``` r

ea <- euclidean_alignment(X_list, input = "raw",
                          mean_method = "logeuclid")
str(ea$aligned, max.level = 1)
#> List of 2
#>  $ : num [1:200, 1:8] 0.0512 -0.1365 -1.2633 0.6798 0.5009 ...
#>  $ : num [1:180, 1:8] -0.794 -1.117 -1.221 0.134 -0.351 ...
```

## Incremental Updates

Update a model with new data without refitting from scratch:

``` r

X_new <- matrix(rnorm(80 * 16), 80, 16)
m_updated <- whiten_model_update(m, X_new)
m_updated$n_obs
#> [1] 380
```

## Hyperparameter Tuning

Automatically select the best method and parameters via
cross-validation:

``` r

tuned <- auto_tune_whitening(
  X_train,
  methods      = c("PCA", "ZCA"),
  n_comp_grid  = c(8, 12),
  lambda_grid  = list("auto", 0.1),
  cv_folds     = 3,
  seed         = 123,
  top_n        = 4
)
head(tuned$ranking)
#>   method n_comp var_threshold lambda cov_estimator mean_score  sd_score
#> 1    PCA      8            NA    0.1     empirical -0.8503463 0.1056367
#> 2    PCA     12            NA    0.1     empirical -1.1336164 0.2216557
#> 3    PCA      8            NA   auto     empirical -1.1809458 0.6216179
#> 4    PCA     12            NA   auto     empirical -1.4801419 0.6928170
#>   n_success
#> 1         3
#> 2         3
#> 3         3
#> 4         3
```

## Tensor Data

Handle 3D EEG tensors (trials x channels x time):

``` r

X_tensor <- array(rnorm(20 * 8 * 50), dim = c(20, 8, 50))
m_t <- whiten_model_tensor(X_tensor, method = "ZCA", lambda = 0)
Z_tensor <- predict_tensor(m_t, X_tensor)
dim(Z_tensor)
#> [1] 20  8 50
```

## Reporting

Generate a markdown report:

``` r

txt <- report_whitening(m, data = X_train, file = NULL)
cat(substr(txt, 1, 300), "...\n")
#> # Whitening Report
#> 
#> - Generated at: 2026-06-27 19:52:26 UTC
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
#> - sample_weighted: FALSE
#> - eig_method: aut ...
```
