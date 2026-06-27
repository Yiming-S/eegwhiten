# Optimized inverse transform from whitened space

Uses cached decomposition terms from `whiten_model` when available.

## Usage

``` r
unwhiten_fast(model, Z)
```

## Arguments

- model:

  A `whiten_model` object.

- Z:

  Whitened data matrix.

## Value

Approximate reconstruction of the original data matrix.

## Details

The deprecated argument order `unwhiten_fast(Z, model)` is still
accepted (with a warning) for backward compatibility.

## See also

[`unwhiten`](https://yiming-s.github.io/eegwhiten/reference/unwhiten.md),
[`predict.whiten_model`](https://yiming-s.github.io/eegwhiten/reference/predict.whiten_model.md)

Other inverse transforms:
[`inverse_ea()`](https://yiming-s.github.io/eegwhiten/reference/inverse_ea.md),
[`unwhiten()`](https://yiming-s.github.io/eegwhiten/reference/unwhiten.md),
[`unwhiten_tensor()`](https://yiming-s.github.io/eegwhiten/reference/unwhiten_tensor.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(200 * 6), 200, 6)
m <- whiten_model(X, method = "PCA", n_comp = 4, lambda = 0.1)
Z <- predict(m, X)
X_rec <- unwhiten_fast(m, Z)
```
