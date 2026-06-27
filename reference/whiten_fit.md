# Fit a whitening matrix from a covariance matrix

Compute a whitening matrix `W` from a symmetric positive-definite
covariance matrix using one of the implemented methods.

## Usage

``` r
whiten_fit(
  Sigma,
  method = c("SVD", "ZCA", "ZCA-cor", "PCA", "PCA-cor", "Cholesky"),
  n_comp = NULL,
  var_threshold = NULL,
  lambda = 0,
  eig_method = c("auto", "base", "rspectra"),
  fast = FALSE,
  sign_reference = NULL
)
```

## Arguments

- Sigma:

  Symmetric positive-definite covariance matrix.

- method:

  Whitening method. Underscore aliases `"ZCA_cor"` and `"PCA_cor"` are
  also accepted.

- n_comp:

  Integer (optional); Number of components to keep (for
  PCA/PCA-cor/SVD).

- var_threshold:

  Numeric in (0, 1\]; cumulative explained variance threshold used to
  choose `n_comp` automatically for `"PCA"`, `"PCA-cor"`, and `"SVD"`.

- lambda:

  Numeric in \[0, 1\]; optional shrinkage applied directly to `Sigma`
  toward a scaled identity target before whitening.

- eig_method:

  Eigen solver backend; one of `"auto"`, `"base"`, or `"rspectra"`.

- fast:

  Logical; if `TRUE`, allow faster approximate settings for iterative
  eigensolvers.

- sign_reference:

  Optional reference vectors used for stable component sign orientation
  in `"PCA"`, `"PCA-cor"`, and `"SVD"`.

## Value

Whitening matrix `W`.

## See also

[`whiten_model`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md)
for the full fit/transform workflow,
[`whiten_matrix`](https://yiming-s.github.io/eegwhiten/reference/whiten_matrix.md)
for one-shot whitening.

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(200 * 8), 200, 8)
S <- cov(X)
W <- whiten_fit(S, method = "PCA", n_comp = 4, lambda = 0.1)
dim(W)
#> [1] 4 8
```
