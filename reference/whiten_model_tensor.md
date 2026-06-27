# Fit a whitening model from a 3D EEG tensor

Flattens `[n_trials, n_channels, n_time]` into
`[n_trials * n_time, n_channels]` and fits
[`whiten_model()`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md).

## Usage

``` r
whiten_model_tensor(X_tensor, ...)
```

## Arguments

- X_tensor:

  Numeric 3D array `[n_trials, n_channels, n_time]`.

- ...:

  Passed to
  [`whiten_model()`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md).

## Value

A `whiten_model` object.

## See also

[`predict_tensor`](https://yiming-s.github.io/eegwhiten/reference/predict_tensor.md),
[`unwhiten_tensor`](https://yiming-s.github.io/eegwhiten/reference/unwhiten_tensor.md),
[`whiten_model`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md)

## Examples

``` r
set.seed(1)
X_t <- array(rnorm(20 * 6 * 40), dim = c(20, 6, 40))
m <- whiten_model_tensor(X_t, method = "ZCA", lambda = 0)
```
