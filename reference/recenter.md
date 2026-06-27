# Recenter a single EEG matrix to the identity

Whitens one data matrix by the inverse square root of its own covariance
(or a supplied reference covariance), mapping that covariance to the
identity. Applied independently to each session/subject, this is the
per-domain recentering at the heart of Euclidean Alignment (He and Wu,
2019) and Riemannian recentering (Zanini et al., 2018): every domain's
mean covariance is moved to a common point, which removes the bulk of
between-session / between-subject covariance shift before
classification.

## Usage

``` r
recenter(X, reference = NULL, center = TRUE, lambda = 0, eps = 1e-10)
```

## Arguments

- X:

  Numeric matrix `[n_samples x n_channels]`.

- reference:

  Optional reference covariance `[p x p]`. If `NULL` (default), the
  covariance of `X` itself is used, mapping `X` to the identity.

- center:

  Logical; whether to subtract the column means before computing the
  covariance and before whitening.

- lambda:

  Numeric in `[0, 1]`; optional shrinkage of the covariance toward a
  scaled identity for stability with short epochs.

- eps:

  Numeric stability floor for eigenvalues.

## Value

A list with:

- Z:

  Recentered data matrix.

- W:

  Alignment matrix `reference^{-1/2}`.

- reference_cov:

  Covariance that was mapped to the identity.

- center:

  Column means subtracted (zeros if `center = FALSE`).

## Details

Unlike
[`euclidean_alignment`](https://yiming-s.github.io/eegwhiten/reference/euclidean_alignment.md),
which applies a single shared transform built from the grand-mean
covariance, `recenter` aligns each matrix by its own covariance, so each
domain is individually whitened to the identity.

## References

He, H., and Wu, D. (2019). Transfer learning for brain-computer
interfaces: A Euclidean space data alignment approach. IEEE Transactions
on Biomedical Engineering.

## See also

[`euclidean_alignment`](https://yiming-s.github.io/eegwhiten/reference/euclidean_alignment.md),
[`whiten_batch`](https://yiming-s.github.io/eegwhiten/reference/whiten_batch.md)
with `mode = "recenter"`.

Other alignment:
[`euclidean_alignment()`](https://yiming-s.github.io/eegwhiten/reference/euclidean_alignment.md),
[`whiten_relative()`](https://yiming-s.github.io/eegwhiten/reference/whiten_relative.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(200 * 6), 200, 6)
rc <- recenter(X)
diag(round(cov(rc$Z), 6))
#> [1] 1 1 1 1 1 1
```
