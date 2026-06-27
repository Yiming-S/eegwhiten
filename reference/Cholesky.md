# Cholesky-based whitening on a covariance matrix

Compute a whitening matrix and factor loadings based on the Cholesky
factorization of a symmetric positive-definite covariance matrix. This
is useful for EEG data where the covariance is estimated across channels
or features.

## Usage

``` r
Cholesky(Sigma, returnW = TRUE, PhiPsi = TRUE, return_decomp = FALSE)
```

## Arguments

- Sigma:

  Symmetric positive-definite covariance matrix.

- returnW:

  Logical; if TRUE, return the whitening matrix `W`.

- PhiPsi:

  Logical; if TRUE, return factor loadings `Phi` and standardized
  loadings `Psi`.

- return_decomp:

  Logical; if TRUE, return decomposition terms used for fast inverse
  transformation.

## Value

A list with some of the elements:

- W:

  Whitening matrix (if `returnW = TRUE`).

- Phi:

  Factor loadings in the original space (if `PhiPsi = TRUE`).

- Psi:

  Standardized loadings (if `PhiPsi = TRUE`).

## See also

[`PCA`](https://yiming-s.github.io/eegwhiten/reference/PCA.md),
[`ZCA`](https://yiming-s.github.io/eegwhiten/reference/ZCA.md),
[`whiten_fit`](https://yiming-s.github.io/eegwhiten/reference/whiten_fit.md)

## Examples

``` r
S <- cov(matrix(rnorm(200 * 6), 200, 6))
res <- Cholesky(S)
dim(res$W)
#> [1] 6 6
```
