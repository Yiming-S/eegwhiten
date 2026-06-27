# Mean of a set of SPD covariance matrices

Computes the average of several symmetric positive-definite covariance
matrices under one of three geometries: the affine-invariant Riemannian
(geometric / Karcher) mean, the log-Euclidean mean, or the arithmetic
(Euclidean) mean.

## Usage

``` r
spd_mean(
  covs,
  metric = c("riemann", "logeuclid", "euclid"),
  tol = 1e-06,
  max_iter = 50,
  eps = 1e-12
)
```

## Arguments

- covs:

  A list of `[p x p]` SPD matrices, or a 3D array `[n, p, p]`.

- metric:

  One of `"riemann"`, `"logeuclid"`, or `"euclid"`.

- tol:

  Convergence tolerance for the Riemannian mean.

- max_iter:

  Maximum iterations for the Riemannian mean.

- eps:

  Numeric stability floor for eigenvalues.

## Value

A `[p x p]` SPD matrix.

## See also

[`riemann_distance`](https://yiming-s.github.io/eegwhiten/reference/riemann_distance.md),
[`barycenter_whitener`](https://yiming-s.github.io/eegwhiten/reference/barycenter_whitener.md),
[`tangent_space`](https://yiming-s.github.io/eegwhiten/reference/tangent_space.md)

Other riemannian:
[`epoch_covariances()`](https://yiming-s.github.io/eegwhiten/reference/epoch_covariances.md),
[`riemann_distance()`](https://yiming-s.github.io/eegwhiten/reference/riemann_distance.md),
[`tangent_space()`](https://yiming-s.github.io/eegwhiten/reference/tangent_space.md),
[`untangent_space()`](https://yiming-s.github.io/eegwhiten/reference/untangent_space.md)

## Examples

``` r
set.seed(1)
covs <- lapply(1:5, function(i) cov(matrix(rnorm(200 * 4), 200, 4)))
M <- spd_mean(covs, metric = "riemann")
dim(M)
#> [1] 4 4
```
