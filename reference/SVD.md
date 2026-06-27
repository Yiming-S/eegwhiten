# SVD-based whitening on the covariance matrix

Perform whitening using the singular value decomposition (SVD) of a
symmetric positive-definite covariance matrix. For SPD matrices, SVD and
eigen-decomposition coincide up to ordering, but this function may be
convenient in some workflows.

## Usage

``` r
SVD(
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

  Integer; number of components to keep. If NULL, keeps all.

- sign_ref:

  Optional reference vectors used to stabilize component signs across
  runs.

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
  decomposition. Reused only when it holds at least `n_comp` components.

## Value

A list with some of the elements:

- W:

  Whitening matrix based on the SVD.

- Phi:

  Factor loadings in the original space.

- Psi:

  Standardized loadings.

## See also

[`PCA`](https://yiming-s.github.io/eegwhiten/reference/PCA.md),
[`ZCA`](https://yiming-s.github.io/eegwhiten/reference/ZCA.md),
[`whiten_fit`](https://yiming-s.github.io/eegwhiten/reference/whiten_fit.md)

## Examples

``` r
S <- cov(matrix(rnorm(200 * 6), 200, 6))
res <- SVD(S, n_comp = 3)
dim(res$W)
#> [1] 3 6
```
