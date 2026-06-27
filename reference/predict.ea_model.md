# Apply an EA Model to New Data

Apply an EA Model to New Data

## Usage

``` r
# S3 method for class 'ea_model'
predict(object, newdata, ...)
```

## Arguments

- object:

  An `ea_model` object returned by
  [`euclidean_alignment()`](https://yiming-s.github.io/eegwhiten/reference/euclidean_alignment.md).

- newdata:

  A numeric matrix or a list of numeric matrices.

- ...:

  Unused.

## Value

Aligned matrix or list of aligned matrices.

## See also

[`euclidean_alignment`](https://yiming-s.github.io/eegwhiten/reference/euclidean_alignment.md),
[`inverse_ea`](https://yiming-s.github.io/eegwhiten/reference/inverse_ea.md)

## Examples

``` r
set.seed(1)
X1 <- matrix(rnorm(200 * 5), 200, 5)
X2 <- matrix(rnorm(180 * 5), 180, 5)
ea <- euclidean_alignment(list(X1, X2), input = "raw")
Z_new <- predict(ea$model, matrix(rnorm(100 * 5), 100, 5))
```
