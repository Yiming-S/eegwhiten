# Inverse transform from whitened space back to original EEG space

Inverse transform from whitened space back to original EEG space

## Usage

``` r
unwhiten(model, Z)
```

## Arguments

- model:

  A `whiten_model` object.

- Z:

  Whitened data matrix.

## Value

Approximate reconstruction of the original data matrix.

## Details

The deprecated argument order `unwhiten(Z, model)` is still accepted
(with a warning) for backward compatibility.

## See also

[`unwhiten_fast`](https://yiming-s.github.io/eegwhiten/reference/unwhiten_fast.md),
[`predict.whiten_model`](https://yiming-s.github.io/eegwhiten/reference/predict.whiten_model.md)

Other inverse transforms:
[`inverse_ea()`](https://yiming-s.github.io/eegwhiten/reference/inverse_ea.md),
[`unwhiten_fast()`](https://yiming-s.github.io/eegwhiten/reference/unwhiten_fast.md),
[`unwhiten_tensor()`](https://yiming-s.github.io/eegwhiten/reference/unwhiten_tensor.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(200 * 6), 200, 6)
m <- whiten_model(X, method = "ZCA", lambda = 0)
Z <- predict(m, X)
X_rec <- unwhiten(m, Z)
max(abs(X_rec - X))
#> [1] 1.887379e-14
```
