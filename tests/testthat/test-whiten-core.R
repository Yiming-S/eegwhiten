############################################################
# Tests for core whitening API: whiten_fit, whiten_model,
# predict, unwhiten, check_whitening, check_condition
############################################################

test_that("whiten_fit returns correct dimensions for all methods", {
  set.seed(10)
  S <- cov(matrix(rnorm(200 * 6), 200, 6))

  for (m in c("PCA", "PCA-cor", "ZCA", "ZCA-cor", "SVD", "Cholesky")) {
    W <- whiten_fit(S, method = m, lambda = 0)
    expect_true(is.matrix(W), info = m)
    expect_equal(ncol(W), 6, info = m)
    expect_true(all(is.finite(W)), info = m)
  }
})

test_that("whiten_fit with n_comp reduces rows", {
  set.seed(11)
  S <- cov(matrix(rnorm(200 * 8), 200, 8))

  for (m in c("PCA", "PCA-cor", "SVD")) {
    W <- whiten_fit(S, method = m, n_comp = 3, lambda = 0)
    expect_equal(nrow(W), 3, info = m)
    expect_equal(ncol(W), 8, info = m)
  }
})

test_that("whiten_fit with var_threshold selects n_comp automatically", {
  set.seed(12)
  S <- cov(matrix(rnorm(200 * 6), 200, 6))

  W <- whiten_fit(S, method = "PCA", var_threshold = 0.5, lambda = 0)
  expect_true(nrow(W) >= 1)
  expect_true(nrow(W) <= 6)
})

test_that("whiten_fit n_comp and var_threshold cannot both be set", {
  S <- cov(matrix(rnorm(200 * 4), 200, 4))
  expect_error(whiten_fit(S, method = "PCA", n_comp = 2, var_threshold = 0.9))
})

test_that("whiten_fit warns for ill-conditioned Sigma with lambda = 0", {
  d <- 4
  S <- diag(c(1e10, 1, 1, 1e-2))
  expect_warning(whiten_fit(S, method = "PCA", lambda = 0),
                 "ill-conditioned")
})

test_that("whiten_fit with lambda > 0 applies shrinkage", {
  set.seed(13)
  S <- cov(matrix(rnorm(200 * 5), 200, 5))
  W0 <- whiten_fit(S, method = "ZCA", lambda = 0)
  W1 <- whiten_fit(S, method = "ZCA", lambda = 0.5)
  expect_false(isTRUE(all.equal(W0, W1)))
})

test_that("whiten_model basic fit and predict roundtrip", {
  set.seed(20)
  X <- matrix(rnorm(300 * 8), 300, 8)
  colnames(X) <- paste0("Ch", 1:8)

  m <- whiten_model(X, method = "ZCA", lambda = 0)
  expect_s3_class(m, "whiten_model")
  expect_equal(m$dim_in, 8)
  expect_equal(m$dim_out, 8)
  expect_equal(m$method, "ZCA")

  Z <- predict(m, X)
  expect_equal(dim(Z), dim(X))
  expect_true(all(is.finite(Z)))

  # Covariance should be close to identity
  S <- cov(Z)
  expect_lt(max(abs(diag(S) - 1)), 0.15)
})

test_that("whiten_model with center = FALSE does not center", {
  set.seed(21)
  X <- matrix(rnorm(200 * 5) + 10, 200, 5)

  m <- whiten_model(X, center = FALSE, method = "PCA", lambda = 0)
  expect_true(all(m$center == 0))

  # predict should not center either
  Z <- predict(m, X)
  expect_true(all(is.finite(Z)))
})

test_that("whiten_model n_comp reduces output dimension", {
  set.seed(22)
  X <- matrix(rnorm(200 * 10), 200, 10)

  m <- whiten_model(X, method = "PCA", n_comp = 4, lambda = 0.1)
  expect_equal(m$dim_out, 4)

  Z <- predict(m, X)
  expect_equal(ncol(Z), 4)
})

test_that("whiten_model assigns column names when missing", {
  set.seed(23)
  X <- matrix(rnorm(100 * 5), 100, 5)
  m <- whiten_model(X, method = "SVD", lambda = 0)
  expect_true(!is.null(names(m$center)))
})

test_that("predict errors on dimension mismatch", {
  set.seed(24)
  X <- matrix(rnorm(100 * 6), 100, 6)
  m <- whiten_model(X, method = "ZCA", lambda = 0)

  X_bad <- matrix(rnorm(100 * 5), 100, 5)
  expect_error(predict(m, X_bad), "Mismatch")
})

test_that("predict errors on non-matrix inputs", {
  set.seed(25)
  X <- matrix(rnorm(100 * 4), 100, 4)
  m <- whiten_model(X, method = "ZCA", lambda = 0)

  expect_error(predict(m, data.frame(X)), "must be a matrix")
  expect_error(predict(m, matrix("a", 2, 4)), "must be a numeric")
})

test_that("unwhiten recovers original data for full-rank method", {
  set.seed(26)
  X <- matrix(rnorm(200 * 6), 200, 6)
  m <- whiten_model(X, method = "ZCA", lambda = 0)
  Z <- predict(m, X)
  X_rec <- unwhiten(Z, m)
  expect_lt(mean((X_rec - X)^2), 1e-20)
})

test_that("unwhiten errors on non-whiten_model", {
  Z <- matrix(1:4, 2, 2)
  expect_error(unwhiten(Z, list()), "must be a 'whiten_model'")
})

test_that("check_whitening returns correct structure", {
  set.seed(27)
  X <- matrix(rnorm(200 * 5), 200, 5)
  m <- whiten_model(X, method = "PCA", lambda = 0.1)
  Z <- predict(m, X)
  d <- check_whitening(Z)

  expect_true(is.list(d))
  expect_true("diag_range" %in% names(d))
  expect_true("diag_mean" %in% names(d))
  expect_true("offdiag_frob" %in% names(d))
  expect_true("dim" %in% names(d))
  expect_true("cov_matrix" %in% names(d))
  expect_equal(d$dim, 5)
})

test_that("check_whitening errors on bad input", {
  expect_error(check_whitening("not a matrix"), "must be a matrix")
  expect_error(check_whitening(matrix(NA, 2, 2)), "numeric")
  expect_error(check_whitening(matrix(c(1, NA, 3, 4), 2, 2)), "finite")
})

test_that("check_condition returns positive numeric", {
  set.seed(28)
  X <- matrix(rnorm(100 * 5), 100, 5)
  cn <- check_condition(X)
  expect_true(is.numeric(cn))
  expect_true(cn > 0)
  expect_true(is.finite(cn))
})

test_that("check_condition errors on bad input", {
  expect_error(check_condition("not a matrix"), "must be a matrix")
  expect_error(check_condition(matrix("a", 2, 2)), "must be a numeric")
  expect_error(check_condition(matrix(1, 1, 3)), "at least 2 rows")
  expect_error(check_condition(matrix(1, 3, 1)), "at least 2 columns")
})
