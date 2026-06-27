# Relative whitening via a generalized eigenvalue decomposition

Computes a linear transform that simultaneously whitens a reference
covariance to the identity and diagonalizes a target covariance. This is
the generalized eigenvalue problem `Sigma w = gamma * reference w` that
underlies spatial filtering methods such as CSP and xDAWN: the resulting
filters maximize the ratio of `Sigma` to `reference` variance.

## Usage

``` r
whiten_relative(Sigma, reference, n_comp = NULL, eps = 1e-10)
```

## Arguments

- Sigma:

  Target symmetric positive-definite covariance `[p x p]` (for example a
  signal- or class-conditional covariance).

- reference:

  Reference symmetric positive-definite covariance `[p x p]` (for
  example a noise or baseline covariance) that is mapped to the
  identity.

- n_comp:

  Optional integer; number of leading components (largest variance
  ratios) to keep. If `NULL`, all `p` are returned.

- eps:

  Numeric stability floor for eigenvalues.

## Value

A list with:

- W:

  Filter matrix `[n_comp x p]`; apply as `Z = X %*% t(W)`.

- values:

  Generalized eigenvalues (variance ratios), descending.

- patterns:

  Forward patterns `[p x n_comp]` (columns), the spatial activity
  pattern of each component.

- reference_inv_sqrt:

  The reference inverse square root `R`.

## Details

Given the reference inverse square root `R = reference^{-1/2}` and the
eigendecomposition `R Sigma R = V Gamma V^T`, the filters are
`W = V^T R`. They satisfy `W reference W^T = I` (reference is whitened)
and `W Sigma W^T = Gamma` (target is diagonal), with the diagonal
`Gamma` giving the per-component variance ratio.

## References

Blankertz, B., Tomioka, R., Lemm, S., Kawanabe, M., and Mueller, K.-R.
(2008). Optimizing spatial filters for robust EEG single-trial analysis.
IEEE Signal Processing Magazine.

## See also

[`whiten_model`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md),
[`euclidean_alignment`](https://yiming-s.github.io/eegwhiten/reference/euclidean_alignment.md)

Other alignment:
[`barycenter_whitener()`](https://yiming-s.github.io/eegwhiten/reference/barycenter_whitener.md),
[`euclidean_alignment()`](https://yiming-s.github.io/eegwhiten/reference/euclidean_alignment.md),
[`recenter()`](https://yiming-s.github.io/eegwhiten/reference/recenter.md)

## Examples

``` r
set.seed(1)
p <- 6
Xa <- matrix(rnorm(200 * p), 200, p)
Xb <- matrix(rnorm(200 * p), 200, p) %*% diag(seq_len(p))
res <- whiten_relative(cov(Xb), cov(Xa))
round(res$values, 3)
#> [1] 35.606 25.780 14.661  9.248  4.259  1.126
```
