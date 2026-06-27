# Apply a whitening model to a 3D EEG tensor

Apply a whitening model to a 3D EEG tensor

## Usage

``` r
predict_tensor(model, X_tensor)
```

## Arguments

- model:

  A `whiten_model` object.

- X_tensor:

  Numeric 3D array `[n_trials, n_channels, n_time]`.

## Value

Whitened 3D array `[n_trials, n_components, n_time]`.

## See also

[`whiten_model_tensor`](https://yiming-s.github.io/eegwhiten/reference/whiten_model_tensor.md),
[`unwhiten_tensor`](https://yiming-s.github.io/eegwhiten/reference/unwhiten_tensor.md)

## Examples

``` r
set.seed(1)
X_t <- array(rnorm(20 * 6 * 40), dim = c(20, 6, 40))
m <- whiten_model_tensor(X_t, method = "ZCA", lambda = 0)
Z_t <- predict_tensor(m, X_t)
dim(Z_t)
#> [1] 20  6 40
```
