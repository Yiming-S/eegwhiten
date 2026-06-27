# PCA-based whitening on the covariance matrix

Perform whitening based on the eigen-decomposition. Supports
dimensionality reduction via `n_comp`.

## Usage

``` r
PCA(
  Sigma,
  n_comp = NULL,
  sign_ref = NULL,
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

- n_comp:

  Integer; number of components to keep.

- sign_ref:

  Optional reference vectors used to stabilize component signs across
  runs. Must have compatible dimensions with eigenvectors.

- returnW:

  Logical; if TRUE, return the whitening matrix `W`.

- PhiPsi:

  Logical; if TRUE, return factor loadings.

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
  decomposition. Reused only when it holds at least `n_comp` components.

## See also

[`whiten_fit`](https://yiming-s.github.io/eegwhiten/reference/whiten_fit.md),
[`PCA_cor`](https://yiming-s.github.io/eegwhiten/reference/PCA_cor.md),
[`ZCA`](https://yiming-s.github.io/eegwhiten/reference/ZCA.md)

## Examples

``` r
S <- cov(matrix(rnorm(200 * 6), 200, 6))
res <- PCA(S, n_comp = 3)
dim(res$W)
#> [1] 3 6
```
