test_that("unwhiten works for full and reduced whitening models", {
  set.seed(123)
  X <- matrix(rnorm(300 * 10), 300, 10)
  colnames(X) <- paste0("Ch", seq_len(ncol(X)))

  m_full <- whiten_model(X, method = "ZCA", lambda = 0)
  Z_full <- predict(m_full, X)
  X_full_rec <- unwhiten(Z_full, m_full)

  expect_equal(dim(X_full_rec), dim(X))
  expect_lt(mean((X_full_rec - X) ^ 2), 1e-20)

  m_red <- whiten_model(X, method = "PCA", n_comp = 4, lambda = 0.1)
  Z_red <- predict(m_red, X)
  X_red_rec <- unwhiten(Z_red, m_red)

  expect_equal(dim(X_red_rec), dim(X))
  expect_true(all(is.finite(X_red_rec)))

  center_mat <- matrix(colMeans(X), nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
  baseline_mse <- mean((X - center_mat) ^ 2)
  expect_lt(mean((X_red_rec - X) ^ 2), baseline_mse)
})

test_that("unwhiten_fast matches unwhiten", {
  set.seed(77)
  X <- matrix(rnorm(250 * 9), 250, 9)

  m <- whiten_model(X, method = "PCA", n_comp = 5, lambda = 0.1)
  Z <- predict(m, X)

  x_fast <- unwhiten_fast(Z, m)
  x_std <- unwhiten(Z, m)

  expect_equal(x_fast, x_std, tolerance = 1e-10)
})

test_that("whiten_batch supports shared model mode", {
  set.seed(456)
  X1 <- matrix(rnorm(150 * 8), 150, 8)
  X2 <- matrix(rnorm(100 * 8), 100, 8)

  shared <- whiten_batch(
    list(X1, X2),
    method = "PCA",
    n_comp = 3,
    lambda = 0.05,
    shared_model = TRUE
  )

  expect_true(isTRUE(attr(shared, "shared_model")))
  expect_equal(dim(shared[[1]]$W), dim(shared[[2]]$W))
  expect_equal(shared[[1]]$W, shared[[2]]$W)
  expect_equal(ncol(shared[[1]]$Z), 3)
  expect_equal(ncol(shared[[2]]$Z), 3)

  independent <- suppressWarnings(whiten_batch(
    list(X1, X2),
    method = "PCA",
    n_comp = 3,
    lambda = 0.05,
    shared_model = FALSE
  ))

  expect_false(isTRUE(attr(independent, "shared_model")))
  expect_gt(max(abs(independent[[1]]$W - independent[[2]]$W)), 1e-6)
})

test_that("whiten_batch supports EA mode", {
  set.seed(999)
  X1 <- matrix(rnorm(200 * 6), 200, 6)
  X2 <- matrix(rnorm(180 * 6), 180, 6)
  X3 <- matrix(rnorm(220 * 6), 220, 6)

  out <- whiten_batch(
    list(X1, X2, X3),
    mode = "ea",
    ea_mean = "logeuclid"
  )

  expect_equal(length(out), 3)
  expect_equal(attr(out, "alignment"), "ea")
  expect_true(isTRUE(attr(out, "shared_model")))

  Z_all <- do.call(rbind, lapply(out, function(x) x$Z))
  S_all <- cov(Z_all)
  expect_lt(norm(S_all - diag(ncol(S_all)), type = "F"), 0.25)
})

test_that("euclidean_alignment works for covariance input", {
  set.seed(1001)
  X1 <- matrix(rnorm(180 * 5), 180, 5)
  X2 <- matrix(rnorm(160 * 5), 160, 5)
  C1 <- cov(X1)
  C2 <- cov(X2)

  ea <- euclidean_alignment(list(C1, C2), input = "cov", mean_method = "riemann")

  expect_equal(length(ea$aligned), 2)
  expect_equal(dim(ea$W), c(5, 5))
  expect_equal(ea$input_type, "cov")

  C_bar <- (ea$aligned[[1]] + ea$aligned[[2]]) / 2
  expect_lt(norm(C_bar - diag(5), type = "F"), 0.25)
})

test_that("var_threshold automatically selects component count", {
  set.seed(11)
  X <- matrix(rnorm(250 * 7), 250, 7)

  S <- cov(X)
  eig <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  expected_k <- which(cumsum(eig) / sum(eig) >= 0.85)[1]

  m <- whiten_model(X, method = "PCA", var_threshold = 0.85, lambda = 0)
  Z <- predict(m, X)

  expect_equal(m$n_comp, expected_k)
  expect_gte(m$explained_var, 0.85)
  expect_equal(ncol(Z), expected_k)
})

test_that("lambda auto works with OAS and LW", {
  set.seed(99)
  X <- matrix(rnorm(40 * 60), 40, 60)

  m_oas <- whiten_model(X, method = "PCA", lambda = "auto", lambda_method = "oas")
  m_lw <- whiten_model(X, method = "PCA", lambda = "auto", lambda_method = "lw")

  expect_gte(m_oas$lambda, 0)
  expect_lte(m_oas$lambda, 1)
  expect_gte(m_lw$lambda, 0)
  expect_lte(m_lw$lambda, 1)

  expect_equal(m_oas$lambda_input, "auto")
  expect_equal(m_oas$lambda_method, "oas")
  expect_equal(m_lw$lambda_input, "auto")
  expect_equal(m_lw$lambda_method, "lw")

  Z <- predict(m_oas, X)
  diag_res <- check_whitening(Z)
  expect_true(is.finite(diag_res$offdiag_frob))
})

test_that("sign_reference is honored for PCA-like methods", {
  set.seed(1234)
  X <- matrix(rnorm(280 * 7), 280, 7)
  ref <- matrix(rnorm(7 * 4), 7, 4)

  m <- whiten_model(X, method = "PCA", n_comp = 4, sign_reference = ref)
  corr <- colSums(m$sign_basis[, 1:4, drop = FALSE] * ref)

  expect_true(all(corr >= -1e-8))
})

test_that("tensor wrappers work end-to-end", {
  set.seed(2026)
  X_tensor <- array(rnorm(18 * 6 * 40), dim = c(18, 6, 40))

  m <- whiten_model_tensor(X_tensor, method = "ZCA", lambda = 0)
  Z_tensor <- predict_tensor(m, X_tensor)
  X_rec <- unwhiten_tensor(Z_tensor, m)

  expect_equal(dim(Z_tensor), c(18, m$dim_out, 40))
  expect_equal(dim(X_rec), dim(X_tensor))
  expect_lt(mean((X_rec - X_tensor) ^ 2), 1e-20)
})

test_that("print summary and plot methods are available", {
  set.seed(202)
  X <- matrix(rnorm(180 * 6), 180, 6)
  m <- whiten_model(X, method = "PCA", var_threshold = 0.9, lambda = "auto")

  expect_no_error(print(m))

  s <- summary(m, data = X)
  expect_s3_class(s, "summary.whiten_model")
  expect_no_error(print(s))

  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  dev <- plot(m, data = X)
  grDevices::dev.off()

  expect_true(file.exists(tmp))
  expect_true(is.matrix(dev))
  expect_equal(nrow(dev), m$dim_out)
  expect_equal(ncol(dev), m$dim_out)
})

test_that("method aliases with underscores are supported", {
  set.seed(303)
  X <- matrix(rnorm(220 * 6), 220, 6)

  m1 <- whiten_model(X, method = "PCA-cor", lambda = 0)
  m2 <- whiten_model(X, method = "PCA_cor", lambda = 0)
  expect_equal(m2$method, "PCA-cor")
  expect_equal(m1$W, m2$W, tolerance = 1e-10)

  z1 <- whiten_model(X, method = "ZCA-cor", lambda = 0)
  z2 <- whiten_model(X, method = "ZCA_cor", lambda = 0)
  expect_equal(z2$method, "ZCA-cor")
  expect_equal(z1$W, z2$W, tolerance = 1e-10)
})

test_that("whiten_fit supports correlation-based methods", {
  set.seed(404)
  X <- matrix(rnorm(300 * 7), 300, 7)
  S <- cov(X)

  W_pca_cor <- whiten_fit(S, method = "PCA-cor", var_threshold = 0.9)
  W_zca_cor <- whiten_fit(S, method = "ZCA-cor")
  W_pca_cor_alias <- whiten_fit(S, method = "PCA_cor", var_threshold = 0.9)

  expect_true(is.matrix(W_pca_cor))
  expect_true(is.matrix(W_zca_cor))
  expect_equal(ncol(W_pca_cor), ncol(S))
  expect_equal(dim(W_zca_cor), dim(S))
  expect_equal(W_pca_cor, W_pca_cor_alias, tolerance = 1e-10)
})

test_that("whiten_batch parallel mode matches sequential mode", {
  skip_if(.Platform$OS.type != "unix")
  set.seed(505)
  X1 <- matrix(rnorm(180 * 6), 180, 6)
  X2 <- matrix(rnorm(160 * 6), 160, 6)
  X3 <- matrix(rnorm(150 * 6), 150, 6)
  data_list <- list(X1, X2, X3)

  seq_out <- suppressWarnings(whiten_batch(
    data_list,
    method = "PCA",
    n_comp = 4,
    lambda = 0.1,
    parallel = FALSE
  ))
  par_out <- suppressWarnings(whiten_batch(
    data_list,
    method = "PCA",
    n_comp = 4,
    lambda = 0.1,
    parallel = TRUE,
    n_cores = 2
  ))

  expect_equal(length(seq_out), length(par_out))
  for (i in seq_along(seq_out)) {
    expect_equal(seq_out[[i]]$Z, par_out[[i]]$Z, tolerance = 1e-10)
    expect_equal(seq_out[[i]]$W, par_out[[i]]$W, tolerance = 1e-10)
  }
})
