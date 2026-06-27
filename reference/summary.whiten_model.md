# Summarize a whitening model

Summarize a whitening model

## Usage

``` r
# S3 method for class 'whiten_model'
summary(object, data = NULL, Z = NULL, ...)
```

## Arguments

- object:

  A `whiten_model` object.

- data:

  Optional numeric matrix in original feature space.

- Z:

  Optional whitened data matrix. If missing and `data` is provided,
  `predict(object, data)` will be used.

- ...:

  Unused.

## Value

An object of class `"summary.whiten_model"`.

## Examples

``` r
X <- matrix(rnorm(200 * 6), 200, 6)
m <- whiten_model(X, method = "PCA", lambda = 0.1)
s <- summary(m, data = X)
print(s)
#> <summary.whiten_model>
#>   method      : PCA 
#>   dimensions  : 6 -> 6 
#>   n_comp      : 6 
#>   var_threshold: NULL 
#>   explained_var: 1 
#>   lambda      : 0.1 
#>   shrink_target: identity 
#>   eig_method  : auto 
#>   fast_mode   : FALSE 
#>   cov_estimator: empirical 
#>   weighted_fit: FALSE 
#>   n_obs       : 200 
#>   updates     : 0 
#>   whitening   : diag_mean=0.9976, offdiag_frob=2.423e-15
```
