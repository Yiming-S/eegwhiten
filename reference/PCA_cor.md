# PCA-based whitening on the correlation matrix

Perform whitening using PCA on the correlation matrix implied by
`Sigma`. This is useful when variables (EEG channels or features) have
very different scales and one wants to whiten in a scale-invariant way.

## Usage

``` r
PCA_cor(
  Sigma,
  n_comp = NULL,
  sign_ref = NULL,
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

  Logical; if `TRUE`, use relaxed checks and faster iterative
  eigensolver settings when available.

## Value

A list with some of the elements:

- W:

  Whitening matrix based on the correlation structure.

- Phi:

  Factor loadings in the original space.

- Psi:

  Standardized loadings.

## See also

[`PCA`](https://yiming-s.github.io/eegwhiten/reference/PCA.md),
[`ZCA_cor`](https://yiming-s.github.io/eegwhiten/reference/ZCA_cor.md),
[`whiten_fit`](https://yiming-s.github.io/eegwhiten/reference/whiten_fit.md)

## Examples

``` r
S <- cov(matrix(rnorm(200 * 6), 200, 6))
res <- PCA_cor(S, n_comp = 3)
dim(res$W)
#> [1] 3 6
```
