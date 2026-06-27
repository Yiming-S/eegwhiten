# ZCA whitening on the covariance matrix

Perform ZCA (Zero-phase Component Analysis) whitening on a symmetric
positive-definite covariance matrix. This transform keeps the whitened
data as close as possible (in L2 sense) to the original data while
decorrelating components.

## Usage

``` r
ZCA(
  Sigma,
  returnW = TRUE,
  PhiPsi = TRUE,
  return_decomp = FALSE,
  eig_method = c("auto", "base", "rspectra"),
  fast = FALSE,
  precomp = NULL
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

  Logical; if `TRUE`, allow faster approximate settings for iterative
  eigensolvers.

- precomp:

  Optional precomputed eigendecomposition of `Sigma` (a list with
  `values` and `vectors`) used internally to avoid recomputing the
  decomposition.

## Value

A list with some of the elements:

- W:

  ZCA whitening matrix.

- Phi:

  Loadings in the original space.

- Psi:

  Standardized loadings.

## See also

[`ZCA_cor`](https://yiming-s.github.io/eegwhiten/reference/ZCA_cor.md),
[`PCA`](https://yiming-s.github.io/eegwhiten/reference/PCA.md),
[`whiten_fit`](https://yiming-s.github.io/eegwhiten/reference/whiten_fit.md)

## Examples

``` r
S <- cov(matrix(rnorm(200 * 6), 200, 6))
res <- ZCA(S)
dim(res$W)
#> [1] 6 6
```
