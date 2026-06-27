# Euclidean Alignment for EEG Matrices

Aligns a list of raw EEG matrices or covariance matrices to a global
reference covariance. This avoids per-matrix independent whitening,
which can distort between-session/class covariance structure.

## Usage

``` r
euclidean_alignment(
  X_list,
  input = c("auto", "raw", "cov"),
  mean_method = c("logeuclid", "riemann", "euclid"),
  center = TRUE,
  tol = 1e-06,
  max_iter = 50,
  eps = 1e-12
)
```

## Arguments

- X_list:

  List of numeric matrices. Each element is either: raw data
  `[n_samples x n_channels]` or covariance `[n_channels x n_channels]`.

- input:

  Input interpretation: `"auto"`, `"raw"`, or `"cov"`.

- mean_method:

  Reference covariance mean: `"riemann"`, `"logeuclid"`, or `"euclid"`.

- center:

  Logical; whether to center each raw matrix before covariance
  computation when `input = "raw"`.

- tol:

  Tolerance for `"riemann"` mean convergence.

- max_iter:

  Maximum iterations for `"riemann"` mean.

- eps:

  Numeric stability floor for eigenvalues.

## Value

A list with:

- aligned:

  Aligned matrices.

- W:

  Global alignment matrix `C_ref^{-1/2}`.

- reference_cov:

  Reference covariance matrix `C_ref`.

- covariances:

  Covariances used to estimate `C_ref`.

- input_type:

  Resolved input type (`"raw"` or `"cov"`).

- mean_method:

  Mean method used.

- model:

  Alignment model object (class `"ea_model"`).

## References

He, H., and Wu, D. (2019). Transfer learning for brain-computer
interfaces: A Euclidean space data alignment approach. IEEE Transactions
on Biomedical Engineering.

## See also

[`predict.ea_model`](https://yiming-s.github.io/eegwhiten/reference/predict.ea_model.md),
[`inverse_ea`](https://yiming-s.github.io/eegwhiten/reference/inverse_ea.md),
[`whiten_batch`](https://yiming-s.github.io/eegwhiten/reference/whiten_batch.md)
with `mode = "ea"`.

Other alignment:
[`recenter()`](https://yiming-s.github.io/eegwhiten/reference/recenter.md),
[`whiten_relative()`](https://yiming-s.github.io/eegwhiten/reference/whiten_relative.md)

## Examples

``` r
set.seed(1)
X1 <- matrix(rnorm(200 * 6), 200, 6)
X2 <- matrix(rnorm(180 * 6), 180, 6)
ea <- euclidean_alignment(list(X1, X2), input = "raw",
                          mean_method = "logeuclid")
str(ea$aligned)
#> List of 2
#>  $ : num [1:200, 1:6] -0.607 0.237 -0.946 1.63 0.305 ...
#>  $ : num [1:180, 1:6] -1.781 1.915 -1.87 -2.209 0.685 ...
```
