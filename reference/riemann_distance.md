# Distance between two SPD covariance matrices

Computes the distance between two symmetric positive-definite (SPD)
covariance matrices under one of three geometries. The affine-invariant
Riemannian metric is the natural distance on the SPD manifold and
underlies covariance-based alignment and classification.

## Usage

``` r
riemann_distance(
  A,
  B,
  metric = c("riemann", "logeuclid", "euclid"),
  eps = 1e-12
)
```

## Arguments

- A, B:

  Symmetric positive-definite matrices of the same dimension.

- metric:

  One of `"riemann"` (affine-invariant Riemannian / AIRM), `"logeuclid"`
  (log-Euclidean), or `"euclid"` (Frobenius).

- eps:

  Numeric stability floor for eigenvalues.

## Value

A non-negative scalar distance.

## Details

For the Riemannian metric the distance is
`sqrt(sum(log(eigvals(A^{-1} B))^2))`; for log-Euclidean it is
`||logm(A) - logm(B)||_F`; for Euclidean it is `||A - B||_F`.

## References

Pennec, X., Fillard, P., and Ayache, N. (2006). A Riemannian framework
for tensor computing. International Journal of Computer Vision.

## See also

[`spd_mean`](https://yiming-s.github.io/eegwhiten/reference/spd_mean.md),
[`tangent_space`](https://yiming-s.github.io/eegwhiten/reference/tangent_space.md)

Other riemannian:
[`epoch_covariances()`](https://yiming-s.github.io/eegwhiten/reference/epoch_covariances.md),
[`spd_mean()`](https://yiming-s.github.io/eegwhiten/reference/spd_mean.md),
[`tangent_space()`](https://yiming-s.github.io/eegwhiten/reference/tangent_space.md),
[`untangent_space()`](https://yiming-s.github.io/eegwhiten/reference/untangent_space.md)

## Examples

``` r
set.seed(1)
A <- cov(matrix(rnorm(200 * 5), 200, 5))
B <- cov(matrix(rnorm(200 * 5), 200, 5))
riemann_distance(A, B)
#> [1] 0.3470935
```
