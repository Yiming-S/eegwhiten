# Inverse Transform for EA-Aligned Data

Inverse Transform for EA-Aligned Data

## Usage

``` r
inverse_ea(model, Z, center = NULL)
```

## Arguments

- model:

  An `ea_model` object.

- Z:

  An aligned matrix (or list of aligned matrices).

- center:

  Optional center vector added back for raw-data inverse mode.

## Value

Matrix (or list of matrices) mapped back from EA-aligned space.

## Details

The deprecated argument order `inverse_ea(Z, model)` is still accepted
(with a warning) for backward compatibility.

## See also

[`euclidean_alignment`](https://yiming-s.github.io/eegwhiten/reference/euclidean_alignment.md),
[`predict.ea_model`](https://yiming-s.github.io/eegwhiten/reference/predict.ea_model.md)

Other inverse transforms:
[`unwhiten()`](https://yiming-s.github.io/eegwhiten/reference/unwhiten.md),
[`unwhiten_fast()`](https://yiming-s.github.io/eegwhiten/reference/unwhiten_fast.md),
[`unwhiten_tensor()`](https://yiming-s.github.io/eegwhiten/reference/unwhiten_tensor.md)

## Examples

``` r
set.seed(1)
X1 <- matrix(rnorm(200 * 5), 200, 5)
X2 <- matrix(rnorm(180 * 5), 180, 5)
ea <- euclidean_alignment(list(X1, X2), input = "raw", center = TRUE)
Z1 <- predict(ea$model, X1)
X1_rec <- inverse_ea(ea$model, Z1)
```
