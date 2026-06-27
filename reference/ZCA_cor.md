# ZCA whitening on the correlation matrix

Perform ZCA whitening using the correlation matrix implied by `Sigma`.
This is a scale-invariant whitening transform that equalizes variance
and removes linear correlations.

## Usage

``` r
ZCA_cor(
  Sigma,
  returnW = TRUE,
  PhiPsi = TRUE,
  return_decomp = FALSE,
  eig_method = c("auto", "base", "rspectra"),
  fast = FALSE
)
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

- eig_method:

  Eigen solver backend; one of `"auto"`, `"base"`, or `"rspectra"`.

- fast:

  Logical; if `TRUE`, use faster iterative eigensolver settings when
  available.

## Value

A list with some of the elements:

- W:

  Whitening matrix based on the correlation structure.

- Phi:

  Factor loadings in the original space.

- Psi:

  Standardized loadings.

## See also

[`ZCA`](https://yiming-s.github.io/eegwhiten/reference/ZCA.md),
[`PCA_cor`](https://yiming-s.github.io/eegwhiten/reference/PCA_cor.md),
[`whiten_fit`](https://yiming-s.github.io/eegwhiten/reference/whiten_fit.md)

## Examples

``` r
S <- cov(matrix(rnorm(200 * 6), 200, 6))
res <- ZCA_cor(S)
dim(res$W)
#> [1] 6 6
```
