# Apply the best tuned whitening model to new data

Convenience method that whitens `newdata` with the best model selected
by
[`auto_tune_whitening`](https://yiming-s.github.io/eegwhiten/reference/auto_tune_whitening.md),
so a tuning result can be used directly as a transformer without
extracting `$best_model`.

## Usage

``` r
# S3 method for class 'whiten_tune'
predict(object, newdata, ...)
```

## Arguments

- object:

  A `whiten_tune` object from
  [`auto_tune_whitening`](https://yiming-s.github.io/eegwhiten/reference/auto_tune_whitening.md).

- newdata:

  Numeric matrix to whiten.

- ...:

  Passed to
  [`predict.whiten_model`](https://yiming-s.github.io/eegwhiten/reference/predict.whiten_model.md)
  (e.g. `self_center`).

## Value

Whitened data matrix.
