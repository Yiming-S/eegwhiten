# Inverse transform a whitened 3D EEG tensor

Inverse transform a whitened 3D EEG tensor

## Usage

``` r
unwhiten_tensor(model, Z_tensor)
```

## Arguments

- model:

  A `whiten_model` object.

- Z_tensor:

  Numeric 3D array `[n_trials, n_components, n_time]`.

## Value

Reconstructed 3D array `[n_trials, n_channels, n_time]`.

## Details

The deprecated argument order `unwhiten_tensor(Z_tensor, model)` is
still accepted (with a warning) for backward compatibility.

## See also

[`whiten_model_tensor`](https://yiming-s.github.io/eegwhiten/reference/whiten_model_tensor.md),
[`predict_tensor`](https://yiming-s.github.io/eegwhiten/reference/predict_tensor.md),
[`unwhiten`](https://yiming-s.github.io/eegwhiten/reference/unwhiten.md)

Other inverse transforms:
[`inverse_ea()`](https://yiming-s.github.io/eegwhiten/reference/inverse_ea.md),
[`unwhiten()`](https://yiming-s.github.io/eegwhiten/reference/unwhiten.md),
[`unwhiten_fast()`](https://yiming-s.github.io/eegwhiten/reference/unwhiten_fast.md)

## Examples

``` r
set.seed(1)
X_t <- array(rnorm(20 * 6 * 40), dim = c(20, 6, 40))
m <- whiten_model_tensor(X_t, method = "ZCA", lambda = 0)
Z_t <- predict_tensor(m, X_t)
X_rec <- unwhiten_tensor(m, Z_t)
```
