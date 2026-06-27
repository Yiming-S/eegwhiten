# Apply a whitening model to new EEG data

Apply a whitening model to new EEG data

## Usage

``` r
# S3 method for class 'whiten_model'
predict(object, newdata, self_center = FALSE, ..., recenter = NULL)
```

## Arguments

- object:

  A `whiten_model` object.

- newdata:

  Numeric matrix with the same number of columns as the data used to fit
  the model.

- self_center:

  Logical; if `TRUE`, center `newdata` by its own column means instead
  of the stored training mean. This is useful for cross-session /
  cross-subject transfer where the test data has shifted relative to the
  training distribution. Ignored when the model was fitted with
  `center = FALSE`.

- ...:

  Unused.

- recenter:

  Deprecated; use `self_center`.

## Value

Whitened data matrix (possibly with reduced dimensions).

## See also

[`whiten_model`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md),
[`unwhiten`](https://yiming-s.github.io/eegwhiten/reference/unwhiten.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(200 * 8), 200, 8)
m <- whiten_model(X, method = "ZCA", lambda = 0.1)
Z <- predict(m, X)
dim(Z)
#> [1] 200   8
```
