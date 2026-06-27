# Fit a whitener that maps the barycenter of a set of covariances to identity

Builds a single linear whitener `W = M^{-1/2}`, where `M` is the mean
covariance
([`spd_mean`](https://yiming-s.github.io/eegwhiten/reference/spd_mean.md))
of the supplied covariance matrices under the chosen geometry. Applying
`W` maps data whose covariance is the barycenter `M` to the identity – a
covariance-space alternative to
[`euclidean_alignment`](https://yiming-s.github.io/eegwhiten/reference/euclidean_alignment.md)
that lets you align to the Riemannian (geometric) mean of a set of
session/epoch covariances.

## Usage

``` r
barycenter_whitener(
  covs,
  metric = c("riemann", "logeuclid", "euclid"),
  tol = 1e-06,
  max_iter = 50,
  eps = 1e-12
)
```

## Arguments

- covs:

  A list of `[p x p]` SPD matrices, or a 3D array `[n, p, p]` (e.g. from
  [`epoch_covariances`](https://yiming-s.github.io/eegwhiten/reference/epoch_covariances.md)).

- metric:

  Geometry of the barycenter; one of `"riemann"`, `"logeuclid"`, or
  `"euclid"`.

- tol, max_iter:

  Convergence controls for the Riemannian mean.

- eps:

  Numeric stability floor for eigenvalues.

## Value

An `ea_model` object (see
[`euclidean_alignment`](https://yiming-s.github.io/eegwhiten/reference/euclidean_alignment.md))
with the whitening matrix `W` and reference covariance `M`; apply it
with [`predict()`](https://rdrr.io/r/stats/predict.html) and invert with
[`inverse_ea`](https://yiming-s.github.io/eegwhiten/reference/inverse_ea.md).

## See also

[`spd_mean`](https://yiming-s.github.io/eegwhiten/reference/spd_mean.md),
[`euclidean_alignment`](https://yiming-s.github.io/eegwhiten/reference/euclidean_alignment.md),
[`recenter`](https://yiming-s.github.io/eegwhiten/reference/recenter.md)

Other alignment:
[`euclidean_alignment()`](https://yiming-s.github.io/eegwhiten/reference/euclidean_alignment.md),
[`recenter()`](https://yiming-s.github.io/eegwhiten/reference/recenter.md),
[`whiten_relative()`](https://yiming-s.github.io/eegwhiten/reference/whiten_relative.md)

## Examples

``` r
set.seed(1)
X <- array(rnorm(20 * 6 * 128), dim = c(20, 6, 128))
covs <- epoch_covariances(X, lambda = 0.05)
bw <- barycenter_whitener(covs, metric = "riemann")
Z <- predict(bw, matrix(rnorm(100 * 6), 100, 6))
```
