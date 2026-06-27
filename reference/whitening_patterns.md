# Spatial filters and forward patterns of a whitening model

Returns the linear spatial filters (the rows of the whitening matrix
`W`, mapping channels to whitened components) and the corresponding
forward patterns (mapping components back to channel space). Patterns
are the appropriate object to visualize or interpret as scalp
topographies (Haufe et al., 2014); filter weights alone can be
misleading.

## Usage

``` r
whitening_patterns(model)
```

## Arguments

- model:

  A `whiten_model` object.

## Value

A list with:

- filters:

  Matrix `[n_components x n_channels]`; the whitening filters `W` (each
  row extracts one whitened component).

- patterns:

  Matrix `[n_channels x n_components]`; the forward patterns (columns),
  i.e. how each component projects back onto channels.

## References

Haufe, S., et al. (2014). On the interpretation of weight vectors of
linear models in multivariate neuroimaging. NeuroImage.

## See also

[`whiten_model`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md),
[`unwhiten`](https://yiming-s.github.io/eegwhiten/reference/unwhiten.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(200 * 8), 200, 8)
m <- whiten_model(X, method = "PCA", n_comp = 4, lambda = 0.1)
fp <- whitening_patterns(m)
dim(fp$patterns)
#> [1] 8 4
```
