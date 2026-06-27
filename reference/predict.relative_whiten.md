# Apply a relative-whitening model to new data

Apply a relative-whitening model to new data

## Usage

``` r
# S3 method for class 'relative_whiten'
predict(object, newdata, ...)
```

## Arguments

- object:

  A `relative_whiten` object.

- newdata:

  Numeric matrix with the same number of channels.

- ...:

  Unused.

## Value

Filtered data matrix `newdata %*% t(W)`.
