# Incrementally Update a Whitening Model

Updates a fitted `whiten_model` using only summary moments from new
observations, avoiding a full refit from all historical data. Currently
supports models fitted with `cov_estimator = "empirical"`.

## Usage

``` r
whiten_model_update(model, X_new, sample_weight = NULL, decay = 1)
```

## Arguments

- model:

  A `whiten_model` object.

- X_new:

  Numeric matrix of new samples with `ncol(X_new) = model$dim_in`.

- sample_weight:

  Optional non-negative sample weights for `X_new`.

- decay:

  Forgetting factor in `(0, 1]` applied to the stored moments before
  adding the new batch. `decay = 1` (default) gives the standard
  cumulative update that weights all observations equally; `decay < 1`
  exponentially down-weights past data, which lets the model track
  non-stationary drift in streaming / online settings.

## Value

Updated `whiten_model` object.

## See also

[`whiten_model`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md)

## Examples

``` r
set.seed(1)
X1 <- matrix(rnorm(200 * 8), 200, 8)
m <- whiten_model(X1, method = "ZCA", lambda = 0.1)
X2 <- matrix(rnorm(100 * 8), 100, 8)
m_updated <- whiten_model_update(m, X2)
m_updated$n_obs
#> [1] 300

# Online tracking with exponential forgetting
m_drift <- whiten_model_update(m, X2, decay = 0.9)
```
