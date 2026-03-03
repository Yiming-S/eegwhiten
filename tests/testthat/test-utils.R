############################################################
# Tests for internal utility functions
############################################################

test_that(".normalize_whiten_method handles case and aliases", {
  nm <- eegwhiten:::.normalize_whiten_method
  expect_equal(nm("PCA"), "PCA")
  expect_equal(nm("pca"), "PCA")
  expect_equal(nm("PCA-cor"), "PCA-cor")
  expect_equal(nm("PCA_cor"), "PCA-cor")
  expect_equal(nm("pca_cor"), "PCA-cor")
  expect_equal(nm("pca-cor"), "PCA-cor")
  expect_equal(nm("ZCA"), "ZCA")
  expect_equal(nm("zca"), "ZCA")
  expect_equal(nm("ZCA-cor"), "ZCA-cor")
  expect_equal(nm("ZCA_cor"), "ZCA-cor")
  expect_equal(nm("SVD"), "SVD")
  expect_equal(nm("svd"), "SVD")
  expect_equal(nm("Cholesky"), "Cholesky")
  expect_equal(nm("cholesky"), "Cholesky")
  expect_error(nm("INVALID"), "Unknown method")
})

test_that(".is_reduction_method identifies correctly", {
  ir <- eegwhiten:::.is_reduction_method
  expect_true(ir("PCA"))
  expect_true(ir("PCA-cor"))
  expect_true(ir("SVD"))
  expect_false(ir("ZCA"))
  expect_false(ir("ZCA-cor"))
  expect_false(ir("Cholesky"))
})

test_that(".check_unit_interval validates correctly", {
  cui <- eegwhiten:::.check_unit_interval
  expect_silent(cui(0, "x"))
  expect_silent(cui(0.5, "x"))
  expect_silent(cui(1, "x"))
  expect_error(cui(-0.1, "x"), "between 0 and 1")
  expect_error(cui(1.1, "x"), "between 0 and 1")
  expect_error(cui("bad", "x"), "single finite numeric")
  expect_error(cui(NA, "x"), "single finite numeric")

  # auto
  expect_silent(cui("auto", "x", allow_char_auto = TRUE))
  expect_error(cui("bad", "x", allow_char_auto = TRUE), "numeric")
})

test_that(".check_var_threshold validates correctly", {
  cvt <- eegwhiten:::.check_var_threshold
  expect_silent(cvt(NULL))
  expect_silent(cvt(0.5))
  expect_silent(cvt(1))
  expect_error(cvt(0), "must be in \\(0, 1\\]")
  expect_error(cvt(1.5), "must be in \\(0, 1\\]")
  expect_error(cvt("bad"), "single finite numeric")
})

test_that(".validate_sample_weight validates correctly", {
  vsw <- eegwhiten:::.validate_sample_weight
  expect_null(vsw(NULL, 5))
  expect_equal(vsw(c(1, 2, 3), 3), c(1, 2, 3))
  expect_error(vsw(c(1, 2), 3), "length")
  expect_error(vsw(c(-1, 2, 3), 3), "non-negative")
  expect_error(vsw(c(0, 0, 0), 3), "positive")
})

test_that(".pinv works for simple cases", {
  pinv <- eegwhiten:::.pinv
  A <- diag(c(2, 3, 4))
  Ai <- pinv(A)
  expect_lt(max(abs(Ai - diag(c(0.5, 1/3, 0.25)))), 1e-10)

  # Zero matrix
  A0 <- matrix(0, 3, 3)
  Ai0 <- pinv(A0)
  expect_equal(Ai0, matrix(0, 3, 3))
})

test_that(".symmetrize makes matrix symmetric", {
  sym <- eegwhiten:::.symmetrize
  A <- matrix(c(1, 2, 3, 4), 2, 2)
  B <- sym(A)
  expect_lt(max(abs(B - t(B))), 1e-15)
  expect_equal(B[1, 2], 2.5)
})

test_that(".check_symmetric_pd detects non-PD matrices", {
  cspd <- eegwhiten:::.check_symmetric_pd
  # Not PD
  S_bad <- matrix(c(1, 2, 2, 1), 2, 2)
  expect_error(cspd(S_bad, require_pd = TRUE), "positive-definite")

  # Good SPD
  S_good <- diag(c(1, 2, 3))
  expect_silent(cspd(S_good, require_pd = TRUE))
})

test_that(".n_comp_from_threshold works correctly", {
  nct <- eegwhiten:::.n_comp_from_threshold
  vals <- c(10, 5, 3, 1, 0.5)
  # Total = 19.5
  # cumsum/total = 0.513, 0.769, 0.923, 0.974, 1.0

  expect_equal(nct(vals, 0.5), 1L)
  expect_equal(nct(vals, 0.75), 2L)
  expect_equal(nct(vals, 0.9), 3L)
  expect_equal(nct(vals, 1.0), 5L)
})

test_that(".explained_ratio computes correctly", {
  er <- eegwhiten:::.explained_ratio
  vals <- c(10, 5, 3, 1, 0.5)
  expect_equal(er(vals, 2, total = 19.5), 15 / 19.5, tolerance = 1e-10)
  expect_true(is.na(er(vals, NULL)))
})

test_that(".fix_component_sign is deterministic", {
  fcs <- eegwhiten:::.fix_component_sign
  U <- matrix(c(-0.5, 0.3, 0.8,
                0.7, -0.6, 0.1), nrow = 3, ncol = 2)

  U_fixed <- fcs(U)
  # The max-abs element in each column determines sign
  # Col 1: max abs is 0.8 at row 3 → positive, no flip
  # Col 2: max abs is 0.7 at row 1 → positive, no flip
  expect_true(U_fixed[which.max(abs(U[, 1])), 1] >= 0)
  expect_true(U_fixed[which.max(abs(U[, 2])), 2] >= 0)
})

test_that(".sqrtm_spd and .invsqrtm_spd are inverses", {
  sqrtm <- eegwhiten:::.sqrtm_spd
  invsqrtm <- eegwhiten:::.invsqrtm_spd
  S <- diag(c(4, 9, 16))
  Sh <- sqrtm(S)
  Sih <- invsqrtm(S)

  expect_lt(max(abs(Sh %*% Sih - diag(3))), 1e-10)
  expect_lt(max(abs(Sh %*% Sh - S)), 1e-10)
})

test_that(".logm_spd and .expm_sym are inverses", {
  logm <- eegwhiten:::.logm_spd
  expm <- eegwhiten:::.expm_sym
  S <- diag(c(2, 3, 5))
  expect_lt(max(abs(expm(logm(S)) - S)), 1e-10)
})

test_that(".lambda_lw returns value in [0,1]", {
  set.seed(50)
  X <- matrix(rnorm(100 * 20), 100, 20)
  Xc <- sweep(X, 2, colMeans(X), "-")
  lw <- eegwhiten:::.lambda_lw(Xc)
  expect_gte(lw, 0)
  expect_lte(lw, 1)
})

test_that(".lambda_oas returns value in [0,1]", {
  set.seed(51)
  X <- matrix(rnorm(100 * 20), 100, 20)
  Xc <- sweep(X, 2, colMeans(X), "-")
  oas <- eegwhiten:::.lambda_oas(Xc)
  expect_gte(oas, 0)
  expect_lte(oas, 1)
})

test_that(".weighted_moments computes correctly", {
  wm <- eegwhiten:::.weighted_moments
  X <- matrix(1:6, 3, 2)

  # Uniform weights
  res1 <- wm(X)
  expect_equal(res1$sum_w, 3)
  expect_equal(res1$sum_x, colSums(X))
  expect_equal(res1$sum_xx, crossprod(X))

  # Custom weights
  w <- c(1, 2, 3)
  res2 <- wm(X, sample_weight = w)
  expect_equal(res2$sum_w, 6)
  expect_equal(res2$sum_x, colSums(X * w))
})

test_that(".cov_from_moments matches cov()", {
  set.seed(52)
  X <- matrix(rnorm(100 * 4), 100, 4)
  m <- eegwhiten:::.weighted_moments(X)
  cfit <- eegwhiten:::.cov_from_moments(m$sum_w, m$sum_w2, m$sum_x, m$sum_xx)
  expect_equal(cfit$cov, cov(X), tolerance = 1e-10)
  expect_equal(cfit$mean, colMeans(X), tolerance = 1e-10)
})
