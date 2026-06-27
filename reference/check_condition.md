# Check Condition Number of Data Matrix

Computes the condition number of the sample covariance matrix of `X`,
defined as the ratio of the largest to smallest positive eigenvalue.
High values indicate ill-conditioning.

## Usage

``` r
check_condition(X)
```

## Arguments

- X:

  Numeric matrix with at least 2 rows and 2 columns.

## Value

The condition number (a single positive numeric value, or `Inf` if no
positive eigenvalues exist).

## See also

[`whiten_model`](https://yiming-s.github.io/eegwhiten/reference/whiten_model.md),
[`check_whitening`](https://yiming-s.github.io/eegwhiten/reference/check_whitening.md)

## Examples

``` r
X <- matrix(rnorm(100 * 8), 100, 8)
check_condition(X)
#> [1] 2.275407
```
