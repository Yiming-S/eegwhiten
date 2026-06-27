# Plot covariance deviation from identity for whitened data

Plot covariance deviation from identity for whitened data

## Usage

``` r
# S3 method for class 'whiten_model'
plot(x, data = NULL, Z = NULL, main = NULL, ...)
```

## Arguments

- x:

  A `whiten_model` object.

- data:

  Optional numeric matrix in original feature space.

- Z:

  Optional whitened data matrix.

- main:

  Optional plot title.

- ...:

  Passed to
  [`graphics::image()`](https://rdrr.io/r/graphics/image.html).

## Value

Invisibly returns `cov(Z) - I`.
