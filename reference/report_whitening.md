# One-Click Whitening Report

Generates a markdown report with model settings and whitening
diagnostics.

## Usage

``` r
report_whitening(
  model,
  data = NULL,
  Z = NULL,
  tune_result = NULL,
  file = "whitening-report.md",
  digits = 4
)
```

## Arguments

- model:

  A `whiten_model` or `whiten_tune` object.

- data:

  Optional original-space data matrix used to compute diagnostics.

- Z:

  Optional whitened matrix. If missing and `data` is provided,
  `predict(model, data)` is used.

- tune_result:

  Optional tuning result from
  [`auto_tune_whitening()`](https://yiming-s.github.io/eegwhiten/reference/auto_tune_whitening.md).

- file:

  Output markdown path. If `NULL`, returns report text without writing
  to disk.

- digits:

  Numeric formatting precision.

## Value

If `file` is non-NULL, invisibly returns the output path. Otherwise
returns the markdown text.

## See also

[`whiten_model`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md),
[`auto_tune_whitening`](https://yiming-s.github.io/eegwhiten/reference/auto_tune_whitening.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(200 * 6), 200, 6)
m <- whiten_model(X, method = "ZCA", lambda = 0.1)
txt <- report_whitening(m, data = X, file = NULL)
cat(txt)
#> # Whitening Report
#> 
#> - Generated at: 2026-07-02 04:30:41 UTC
#> - Method: ZCA
#> - Dimensions: 6 -> 6
#> - n_comp: NULL
#> - var_threshold: NULL
#> - explained_var: NA
#> - lambda: 0.1
#> - lambda_input: 0.1
#> - lambda_method: NA
#> - shrink_target: identity
#> - cov_estimator: empirical
#> - sample_weighted: FALSE
#> - eig_method: auto
#> - fast: FALSE
#> - n_obs: 200
#> - update_count: 0
#> 
#> ## Whitening Diagnostics
#> 
#> - diag_mean: 0.9981
#> - diag_min: 0.9757
#> - diag_max: 1.01
#> - offdiag_frob: 0.02447
```
