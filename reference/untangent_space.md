# Map tangent-space vectors back to SPD covariance matrices

Inverse of
[`tangent_space`](https://yiming-s.github.io/eegwhiten/reference/tangent_space.md):
reconstructs SPD matrices from their tangent-space feature vectors and a
reference point.

## Usage

``` r
untangent_space(V, ref = NULL)
```

## Arguments

- V:

  A numeric matrix `[n x p(p+1)/2]` of tangent-space features, or a
  single numeric vector. If `ref` is `NULL`, the reference stored in
  `attr(V, "reference")` (as returned by `tangent_space`) is used.

- ref:

  Optional reference SPD matrix `[p x p]`.

## Value

A list of `[p x p]` SPD matrices.

## See also

[`tangent_space`](https://yiming-s.github.io/eegwhiten/reference/tangent_space.md)

Other riemannian:
[`epoch_covariances()`](https://yiming-s.github.io/eegwhiten/reference/epoch_covariances.md),
[`riemann_distance()`](https://yiming-s.github.io/eegwhiten/reference/riemann_distance.md),
[`spd_mean()`](https://yiming-s.github.io/eegwhiten/reference/spd_mean.md),
[`tangent_space()`](https://yiming-s.github.io/eegwhiten/reference/tangent_space.md)

## Examples

``` r
set.seed(1)
X <- array(rnorm(10 * 5 * 64), dim = c(10, 5, 64))
covs <- epoch_covariances(X)
V <- tangent_space(covs, mean_method = "logeuclid")
covs_rec <- untangent_space(V)
max(abs(covs_rec[[1]] - covs[[1]]))
#> [1] 1.776357e-15
```
