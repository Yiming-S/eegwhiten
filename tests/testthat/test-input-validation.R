############################################################
# Tests for input validation and error handling
############################################################

# --- whiten_model input validation ---

test_that("whiten_model rejects non-matrix input", {
  expect_error(whiten_model(data.frame(a = 1:5, b = 6:10)), "must be a matrix")
})

test_that("whiten_model rejects non-numeric matrix", {
  expect_error(whiten_model(matrix("a", 3, 3)), "must be a numeric")
})

test_that("whiten_model rejects matrix with non-finite values", {
  X <- matrix(c(1, NA, 3, 4), 2, 2)
  expect_error(whiten_model(X), "finite")
})

test_that("whiten_model rejects single-row matrix", {
  expect_error(whiten_model(matrix(1:6, 1, 6)), "at least 2 rows")
})

test_that("whiten_model rejects single-column matrix", {
  expect_error(whiten_model(matrix(1:6, 6, 1)), "more than one column")
})

test_that("whiten_model rejects unknown method", {
  X <- matrix(rnorm(20), 5, 4)
  expect_error(whiten_model(X, method = "INVALID"), "Unknown method")
})

test_that("whiten_model rejects bad center argument", {
  X <- matrix(rnorm(20), 5, 4)
  expect_error(whiten_model(X, center = "yes"), "must be TRUE or FALSE")
})

test_that("whiten_model rejects bad lambda", {
  X <- matrix(rnorm(20), 5, 4)
  expect_error(whiten_model(X, lambda = -1), "between 0 and 1")
  expect_error(whiten_model(X, lambda = 2), "between 0 and 1")
  expect_error(whiten_model(X, lambda = "bad"), "numeric")
})

test_that("whiten_model rejects bad n_comp", {
  X <- matrix(rnorm(100), 20, 5)
  expect_error(whiten_model(X, method = "PCA", n_comp = 0), "must be between 1")
  expect_error(whiten_model(X, method = "PCA", n_comp = 10), "must be between 1")
  expect_error(whiten_model(X, method = "PCA", n_comp = 2.5), "must be an integer")
  expect_error(whiten_model(X, method = "PCA", n_comp = "a"), "single finite numeric")
})

test_that("whiten_model rejects n_comp + var_threshold together", {
  X <- matrix(rnorm(100), 20, 5)
  expect_error(whiten_model(X, method = "PCA", n_comp = 2, var_threshold = 0.9),
               "only one of")
})

test_that("whiten_model warns n_comp for non-reduction methods", {
  X <- matrix(rnorm(100), 20, 5)
  expect_warning(whiten_model(X, method = "ZCA", n_comp = 3, lambda = 0),
                 "n_comp is ignored")
})

test_that("whiten_model warns var_threshold for non-reduction methods", {
  X <- matrix(rnorm(100), 20, 5)
  expect_warning(whiten_model(X, method = "Cholesky", var_threshold = 0.9, lambda = 0),
                 "var_threshold is ignored")
})

test_that("whiten_model rejects bad var_threshold", {
  X <- matrix(rnorm(100), 20, 5)
  expect_error(whiten_model(X, method = "PCA", var_threshold = 0, lambda = 0),
               "must be in \\(0, 1\\]")
  expect_error(whiten_model(X, method = "PCA", var_threshold = 1.5, lambda = 0),
               "must be in \\(0, 1\\]")
})

test_that("whiten_model rejects lambda=auto with sample_weight", {
  X <- matrix(rnorm(100), 20, 5)
  w <- rep(1, 20)
  expect_error(whiten_model(X, lambda = "auto", sample_weight = w),
               "not supported with sample_weight")
})

test_that("whiten_model rejects lambda=auto with non-empirical cov", {
  X <- matrix(rnorm(100), 20, 5)
  expect_error(whiten_model(X, lambda = "auto", cov_estimator = "tyler"),
               "only supported with cov_estimator")
})

test_that("whiten_model rejects bad fast argument", {
  X <- matrix(rnorm(100), 20, 5)
  expect_error(whiten_model(X, fast = "yes"), "must be TRUE or FALSE")
})

# --- whiten_fit input validation ---

test_that("whiten_fit rejects non-PD matrix", {
  S <- matrix(c(1, 2, 2, 1), 2, 2)
  expect_error(whiten_fit(S), "positive-definite")
})

test_that("whiten_fit rejects non-square matrix", {
  S <- matrix(1:6, 2, 3)
  expect_error(whiten_fit(S), "square")
})

# --- sample_weight validation ---

test_that("whiten_model rejects wrong-length sample_weight", {
  X <- matrix(rnorm(100), 20, 5)
  expect_error(whiten_model(X, lambda = 0, sample_weight = rep(1, 10)),
               "length")
})

test_that("whiten_model rejects negative sample_weight", {
  X <- matrix(rnorm(100), 20, 5)
  expect_error(whiten_model(X, lambda = 0, sample_weight = rep(-1, 20)),
               "non-negative")
})

# --- sign_reference validation ---

test_that("whiten_model rejects bad sign_reference", {
  X <- matrix(rnorm(100), 20, 5)
  expect_error(whiten_model(X, method = "PCA", lambda = 0,
                            sign_reference = matrix(1, 3, 3)),
               "nrow.*ncol")
})

test_that("whiten_model warns sign_reference for non-reduction methods", {
  X <- matrix(rnorm(100), 20, 5)
  expect_warning(
    whiten_model(X, method = "ZCA", lambda = 0,
                 sign_reference = diag(5)),
    "sign_reference is ignored"
  )
})

# --- whiten_batch validation ---

test_that("whiten_batch rejects non-list input", {
  expect_error(whiten_batch(matrix(1, 2, 2)), "must be a list")
})

test_that("whiten_batch rejects empty list", {
  expect_error(whiten_batch(list()), "must not be empty")
})

# --- unwhiten validation ---

test_that("unwhiten rejects dimension mismatch", {
  X <- matrix(rnorm(100), 20, 5)
  m <- whiten_model(X, method = "PCA", n_comp = 3, lambda = 0.1)
  Z_bad <- matrix(rnorm(40), 20, 2)
  expect_error(unwhiten(m, Z_bad), "Mismatch")
})

test_that("unwhiten_fast rejects dimension mismatch", {
  X <- matrix(rnorm(100), 20, 5)
  m <- whiten_model(X, method = "PCA", n_comp = 3, lambda = 0.1)
  Z_bad <- matrix(rnorm(40), 20, 2)
  expect_error(unwhiten_fast(m, Z_bad), "Mismatch")
})

# --- check_whitening validation ---

test_that("check_whitening rejects non-numeric matrix", {
  expect_error(check_whitening(matrix("a", 2, 2)), "numeric")
})

# --- euclidean_alignment validation ---

test_that("euclidean_alignment rejects non-list", {
  expect_error(euclidean_alignment(matrix(1, 2, 2)), "non-empty list")
})

test_that("euclidean_alignment rejects mismatched columns", {
  X1 <- matrix(rnorm(20), 10, 2)
  X2 <- matrix(rnorm(30), 10, 3)
  expect_error(euclidean_alignment(list(X1, X2), input = "raw"),
               "same number of columns")
})

# --- whiten_model_update validation ---

test_that("whiten_model_update rejects non-whiten_model", {
  expect_error(whiten_model_update(list(), matrix(1, 2, 2)),
               "must be a 'whiten_model'")
})

test_that("whiten_model_update rejects dimension mismatch", {
  X <- matrix(rnorm(100), 20, 5)
  m <- whiten_model(X, method = "ZCA", lambda = 0)
  X_new <- matrix(rnorm(60), 20, 3)
  expect_error(whiten_model_update(m, X_new), "Mismatch")
})

# --- auto_tune_whitening validation ---

test_that("auto_tune_whitening rejects too few rows", {
  X <- matrix(rnorm(12), 3, 4)
  expect_error(auto_tune_whitening(X), "at least 4 rows")
})

test_that("auto_tune_whitening rejects accuracy without y", {
  X <- matrix(rnorm(100), 20, 5)
  expect_error(auto_tune_whitening(X, scoring = "accuracy"),
               "requires non-NULL y")
})
