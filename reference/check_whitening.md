# Diagnostic check for whitening quality

Diagnostic check for whitening quality

## Usage

``` r
check_whitening(Z)
```

## Arguments

- Z:

  Whitened data matrix.

## Value

A list with elements:

- diag_range:

  Range (min, max) of diagonal entries of `cov(Z)`.

- diag_mean:

  Mean diagonal value.

- offdiag_frob:

  Frobenius norm of off-diagonal entries.

- cov_dev_frob:

  Frobenius norm of `cov(Z) - I` (total deviation from identity,
  including the diagonal).

- logdet:

  Log-determinant of `cov(Z)` (0 for perfectly white data).

- whiteness:

  Bounded score in `[0, 1]`, equal to 1 when `cov(Z) = I` and decreasing
  as the data departs from white.

- dim:

  Number of components.

- cov_matrix:

  The sample covariance of `Z`.

## See also

[`whiten_model`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md),
[`check_condition`](https://yiming-s.github.io/eegwhiten/reference/check_condition.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(200 * 8), 200, 8)
m <- whiten_model(X, method = "ZCA", lambda = 0.1)
Z <- predict(m, X)
d <- check_whitening(Z)
d$diag_mean
#> [1] 0.9966629
```
