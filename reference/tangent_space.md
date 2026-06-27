# Map SPD covariance matrices to the tangent space

Projects a set of symmetric positive-definite (SPD) covariance matrices
to the tangent space at a reference point and vectorizes them. This
linearizes the curved SPD manifold so that ordinary Euclidean
classifiers and regressors can be applied to the resulting feature
vectors – the standard Riemannian pipeline for EEG/BCI.

## Usage

``` r
tangent_space(
  covs,
  ref = NULL,
  mean_method = c("riemann", "logeuclid", "euclid"),
  tol = 1e-06,
  max_iter = 50,
  eps = 1e-12
)
```

## Arguments

- covs:

  A list of `[p x p]` SPD matrices, or a 3D array `[n, p, p]`.

- ref:

  Optional reference SPD matrix `[p x p]`. If `NULL`, the reference is
  the mean covariance computed with `mean_method`.

- mean_method:

  Method for the reference mean when `ref = NULL`; one of `"riemann"`
  (affine-invariant), `"logeuclid"`, or `"euclid"`.

- tol:

  Tolerance for `"riemann"` mean convergence.

- max_iter:

  Maximum iterations for `"riemann"` mean.

- eps:

  Numeric stability floor for eigenvalues.

## Value

A numeric matrix `[n x p(p+1)/2]` of tangent-space features. The
reference and metadata are attached as attributes `"reference"`, `"p"`,
and `"mean_method"`.

## Details

For each covariance `C`, the map computes
`L = logm(ref^{-1/2} C ref^{-1/2})` and returns the upper triangle of
`L` with off-diagonal entries scaled by `sqrt(2)`, so that the Euclidean
norm of the feature vector equals the Riemannian distance from the
reference.

## References

Barachant, A., Bonnet, S., Congedo, M., and Jutten, C. (2012).
Multiclass brain-computer interface classification by Riemannian
geometry. IEEE Transactions on Biomedical Engineering.

## See also

[`untangent_space`](https://yiming-s.github.io/eegwhiten/reference/untangent_space.md),
[`epoch_covariances`](https://yiming-s.github.io/eegwhiten/reference/epoch_covariances.md)

Other riemannian:
[`epoch_covariances()`](https://yiming-s.github.io/eegwhiten/reference/epoch_covariances.md),
[`riemann_distance()`](https://yiming-s.github.io/eegwhiten/reference/riemann_distance.md),
[`spd_mean()`](https://yiming-s.github.io/eegwhiten/reference/spd_mean.md),
[`untangent_space()`](https://yiming-s.github.io/eegwhiten/reference/untangent_space.md)

## Examples

``` r
set.seed(1)
X <- array(rnorm(30 * 6 * 64), dim = c(30, 6, 64))
covs <- epoch_covariances(X, lambda = 0.05)
V <- tangent_space(covs, mean_method = "logeuclid")
dim(V)
#> [1] 30 21
```
