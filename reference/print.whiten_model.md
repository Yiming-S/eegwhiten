# Print a whitening model

Print a whitening model

## Usage

``` r
# S3 method for class 'whiten_model'
print(x, ...)
```

## Arguments

- x:

  A `whiten_model` object.

- ...:

  Unused.

## Value

The input object invisibly.

## Examples

``` r
m <- whiten_model(matrix(rnorm(200 * 6), 200, 6), method = "PCA",
                  lambda = 0.1)
print(m)
#> <whiten_model>
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
```
