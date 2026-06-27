# Per-epoch channel covariance matrices from an EEG tensor

Computes one channel-by-channel covariance matrix per trial/epoch from a
3D EEG tensor. These covariances are the natural input for Riemannian
pipelines such as
[`tangent_space`](https://yiming-s.github.io/eegwhiten/reference/tangent_space.md)
and Euclidean / Riemannian recentering.

## Usage

``` r
epoch_covariances(X_tensor, lambda = 0, center = TRUE, shrinkage = NULL)
```

## Arguments

- X_tensor:

  Numeric 3D array `[n_trials, n_channels, n_time]`.

- lambda:

  Numeric in `[0, 1]` or the string `"auto"`; optional shrinkage of each
  covariance toward a scaled identity
  (`(1 - lambda) C + lambda (tr(C)/p) I`). `"auto"` selects a per-epoch
  Ledoit-Wolf intensity from that epoch's time samples. Useful when the
  number of time samples is small relative to the number of channels.
  Defaults to 0 (no shrinkage). When `n_time <= n_channels`, covariances
  are rank-deficient (not positive-definite) unless `lambda > 0`; a
  warning is emitted in that case.

- center:

  Logical; whether to subtract each channel's mean over time before
  computing the covariance.

- shrinkage:

  Deprecated; use `lambda`.

## Value

A list of `[n_channels x n_channels]` covariance matrices, one per
trial. They are positive-definite (suitable as
[`tangent_space`](https://yiming-s.github.io/eegwhiten/reference/tangent_space.md)
input) only when `n_time > n_channels` or `lambda > 0`.

## See also

[`tangent_space`](https://yiming-s.github.io/eegwhiten/reference/tangent_space.md),
[`recenter`](https://yiming-s.github.io/eegwhiten/reference/recenter.md)

Other riemannian:
[`riemann_distance()`](https://yiming-s.github.io/eegwhiten/reference/riemann_distance.md),
[`spd_mean()`](https://yiming-s.github.io/eegwhiten/reference/spd_mean.md),
[`tangent_space()`](https://yiming-s.github.io/eegwhiten/reference/tangent_space.md),
[`untangent_space()`](https://yiming-s.github.io/eegwhiten/reference/untangent_space.md)

## Examples

``` r
set.seed(1)
X <- array(rnorm(20 * 8 * 64), dim = c(20, 8, 64))
covs <- epoch_covariances(X, lambda = 0.05)
length(covs)
#> [1] 20
```
