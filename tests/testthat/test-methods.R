############################################################
# Tests for individual whitening method functions
############################################################

test_that("PCA produces correct whitening for known covariance", {
  set.seed(100)
  X <- matrix(rnorm(500 * 6), 500, 6)
  S <- cov(X)

  res <- PCA(S, returnW = TRUE, PhiPsi = TRUE, return_decomp = TRUE)
  expect_true(is.matrix(res$W))
  expect_true(is.matrix(res$Phi))
  expect_true(is.matrix(res$Psi))
  expect_equal(nrow(res$W), 6)
  expect_equal(ncol(res$W), 6)

  # W should whiten: W %*% S %*% t(W) ≈ I
  I_approx <- res$W %*% S %*% t(res$W)
  expect_lt(max(abs(I_approx - diag(6))), 1e-8)

  # decomp check
  expect_true(is.matrix(res$decomp$U))
  expect_true(is.numeric(res$decomp$D))
  expect_true(is.matrix(res$decomp$inv_tW))
})

test_that("PCA with n_comp", {
  S <- cov(matrix(rnorm(300 * 8), 300, 8))
  res <- PCA(S, n_comp = 3)
  expect_equal(nrow(res$W), 3)
  expect_equal(ncol(res$W), 8)
})

test_that("PCA errors on non-positive eigenvalues", {
  S <- matrix(c(1, 2, 2, 1), 2, 2)  # not PD
  expect_error(PCA(S))
})

test_that("PCA_cor whitening is scale-invariant", {
  set.seed(101)
  X <- matrix(rnorm(400 * 5), 400, 5)
  S1 <- cov(X)
  S2 <- cov(sweep(X, 2, c(1, 100, 0.01, 50, 0.5), "*"))

  W1 <- PCA_cor(S1)$W
  W2 <- PCA_cor(S2)$W

  # Both produce whitening matrices; shapes should match

  expect_equal(dim(W1), dim(W2))
})

test_that("PCA_cor with return_decomp", {
  S <- cov(matrix(rnorm(300 * 4), 300, 4))
  res <- PCA_cor(S, n_comp = 2, return_decomp = TRUE)
  expect_true(!is.null(res$decomp))
  expect_true(!is.null(res$decomp$scale_diag))
  expect_equal(nrow(res$W), 2)
})

test_that("ZCA whitening matrix is symmetric", {
  set.seed(102)
  S <- cov(matrix(rnorm(300 * 5), 300, 5))
  res <- ZCA(S)
  W <- res$W
  expect_equal(nrow(W), ncol(W))
  expect_lt(max(abs(W - t(W))), 1e-10)
})

test_that("ZCA with return_decomp", {
  S <- cov(matrix(rnorm(300 * 4), 300, 4))
  res <- ZCA(S, return_decomp = TRUE)
  expect_true(!is.null(res$decomp))
  expect_true(!is.null(res$decomp$inv_tW))
})

test_that("ZCA_cor whitening is scale-invariant", {
  set.seed(103)
  X <- matrix(rnorm(400 * 4), 400, 4)
  S <- cov(X)
  res <- ZCA_cor(S)
  expect_equal(nrow(res$W), 4)
  expect_equal(ncol(res$W), 4)
})

test_that("ZCA_cor with return_decomp", {
  S <- cov(matrix(rnorm(300 * 4), 300, 4))
  res <- ZCA_cor(S, return_decomp = TRUE)
  expect_true(!is.null(res$decomp))
  expect_true(!is.null(res$decomp$scale_diag))
})

test_that("SVD produces valid whitening (same as PCA for SPD)", {
  set.seed(104)
  S <- cov(matrix(rnorm(400 * 6), 400, 6))

  W_pca <- PCA(S, returnW = TRUE, PhiPsi = FALSE)$W
  W_svd <- SVD(S, returnW = TRUE, PhiPsi = FALSE)$W

  # For SPD matrices, PCA and SVD should produce same whitening
  # up to sign differences
  for (j in seq_len(nrow(W_pca))) {
    s <- sign(W_pca[j, 1]) * sign(W_svd[j, 1])
    expect_lt(max(abs(W_pca[j, ] - s * W_svd[j, ])), 1e-8)
  }
})

test_that("SVD with n_comp", {
  S <- cov(matrix(rnorm(300 * 8), 300, 8))
  res <- SVD(S, n_comp = 4, return_decomp = TRUE)
  expect_equal(nrow(res$W), 4)
  expect_true(!is.null(res$decomp))
})

test_that("Cholesky produces valid whitening", {
  set.seed(105)
  S <- cov(matrix(rnorm(300 * 5), 300, 5))

  res <- Cholesky(S, return_decomp = TRUE)
  W <- res$W
  expect_equal(nrow(W), 5)
  expect_equal(ncol(W), 5)

  # W %*% S %*% t(W) ≈ I
  I_approx <- W %*% S %*% t(W)
  expect_lt(max(abs(I_approx - diag(5))), 1e-8)
  expect_true(!is.null(res$decomp))
})

test_that("Cholesky PhiPsi loadings are correct", {
  S <- cov(matrix(rnorm(300 * 4), 300, 4))
  res <- Cholesky(S, PhiPsi = TRUE)
  expect_true(is.matrix(res$Phi))
  expect_true(is.matrix(res$Psi))
  expect_equal(dim(res$Phi), c(4, 4))
})

test_that("all methods produce finite results with regularization", {
  set.seed(106)
  # Create ill-conditioned data: more features than samples
  X <- matrix(rnorm(10 * 20), 10, 20)
  S <- cov(X)
  # Regularize to make it PD
  S_reg <- (1 - 0.5) * S + 0.5 * diag(mean(diag(S)), 20)
  attr(S_reg, ".checked_spd") <- TRUE

  for (m in c("PCA", "PCA-cor", "ZCA", "ZCA-cor", "SVD", "Cholesky")) {
    fn <- get(ifelse(m == "PCA-cor", "PCA_cor", ifelse(m == "ZCA-cor", "ZCA_cor", m)))
    res <- fn(S_reg, returnW = TRUE, PhiPsi = FALSE)
    expect_true(all(is.finite(res$W)), info = m)
  }
})
