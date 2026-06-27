############################################################
# Tests for the whitening extensions:
# precomp/perf parity, shrink_target, na_action, recenter,
# patterns, diagnostics, forgetting update, analytic tuning,
# tangent space, relative whitening.
############################################################

make_structured <- function(n, p, seed = 1) {
  set.seed(seed)
  A <- matrix(rnorm(p * p), p, p)
  Sig <- crossprod(A) / p + diag(p)
  X <- matrix(rnorm(n * p), n, p) %*% chol(Sig)
  colnames(X) <- paste0("Ch", seq_len(p))
  X
}

# ---- Single eigendecomposition (precomp) preserves whitening ----

test_that("precomp threading yields correct whitening for cov-based methods", {
  X <- make_structured(300, 10, seed = 2)
  # lambda = 0 whitens the training data to exactly the identity, which checks
  # that the reused single eigendecomposition produces the correct W.
  for (m in c("ZCA", "PCA", "SVD")) {
    model <- whiten_model(X, method = m, lambda = 0)
    d <- check_whitening(predict(model, X))
    expect_equal(d$diag_mean, 1, tolerance = 1e-8)
    expect_lt(d$offdiag_frob, 1e-8)
  }
})

test_that("public ZCA matches analytic whiten primitive", {
  X <- make_structured(200, 6, seed = 3)
  S <- stats::cov(X)
  e <- eigen(S, symmetric = TRUE)
  W_method <- ZCA(S, returnW = TRUE, PhiPsi = FALSE)$W
  W_precomp <- ZCA(S, returnW = TRUE, PhiPsi = FALSE,
                   precomp = list(values = e$values, vectors = e$vectors))$W
  expect_equal(unname(W_method), unname(W_precomp), tolerance = 1e-10)
})

# ---- shrink_target ----

test_that("diagonal shrink_target preserves variances, identity does not", {
  X <- make_structured(250, 8, seed = 4)
  S_raw <- stats::cov(X)
  m_id <- whiten_model(X, method = "ZCA", lambda = 0.5, shrink_target = "identity")
  m_dg <- whiten_model(X, method = "ZCA", lambda = 0.5, shrink_target = "diagonal")
  expect_equal(diag(m_dg$cov), diag(S_raw), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(diag(m_id$cov), diag(S_raw))))
  expect_identical(m_dg$shrink_target, "diagonal")
})

test_that("lambda='auto' is accepted with the diagonal target", {
  X <- make_structured(100, 5, seed = 5)
  m <- whiten_model(X, method = "PCA", lambda = "auto", shrink_target = "diagonal")
  expect_true(m$lambda >= 0 && m$lambda <= 1)
  expect_identical(m$shrink_target, "diagonal")
})

# ---- na_action ----

test_that("na_action='omit' drops non-finite rows and keeps weights aligned", {
  X <- make_structured(120, 6, seed = 6)
  X[c(2, 50, 99), 3] <- NA
  expect_error(whiten_model(X, method = "ZCA", lambda = 0.1), "finite")
  expect_warning(
    m <- whiten_model(X, method = "ZCA", lambda = 0.1, na_action = "omit"),
    "Dropped 3"
  )
  expect_equal(m$n_obs, 117)

  w <- runif(120) + 0.1
  expect_warning(
    mw <- whiten_model(X, method = "ZCA", lambda = 0.1, na_action = "omit",
                       sample_weight = w),
    "Dropped 3"
  )
  expect_equal(mw$n_obs, 117)

  # Weights must stay aligned with the surviving rows: equals an explicit
  # reference fit on the kept rows with the kept weights.
  bad <- c(2, 50, 99)
  keep <- setdiff(seq_len(120), bad)
  ref <- whiten_model(X[keep, , drop = FALSE], method = "ZCA", lambda = 0.1,
                      sample_weight = w[keep])
  expect_equal(mw$W, ref$W, tolerance = 1e-10)
  expect_equal(mw$center, ref$center, tolerance = 1e-10)

  # A weight sized to the survivors (not the original rows) is rejected.
  expect_error(
    whiten_model(X, method = "ZCA", lambda = 0.1, na_action = "omit",
                 sample_weight = runif(117) + 0.1),
    "nrow"
  )
})

test_that("na_action='omit' works in whiten_matrix and whiten_batch", {
  X <- make_structured(100, 5, seed = 21)
  X[c(4, 40), 2] <- NA
  # whiten_matrix must not error and returns Z for the surviving rows
  expect_warning(res <- whiten_matrix(X, method = "ZCA", lambda = 0.1,
                                      na_action = "omit"), "Dropped 2")
  expect_equal(nrow(res$Z), 98)

  Xl <- list(X, make_structured(80, 5, seed = 22))
  expect_warning(out <- whiten_batch(Xl, mode = "shared_model", method = "ZCA",
                                     lambda = 0.1, na_action = "omit"))
  expect_equal(nrow(out[[1]]$Z), 98)
  expect_equal(nrow(out[[2]]$Z), 80)
})

test_that("positional sample_weight call still targets sample_weight (signature compat)", {
  X <- make_structured(80, 5, seed = 23)
  w <- runif(80) + 0.1
  # 8th positional arg must still be sample_weight (not a new option)
  m_pos <- whiten_model(X, TRUE, "ZCA", NULL, NULL, 0.1, "oas", w)
  m_named <- whiten_model(X, method = "ZCA", lambda = 0.1, sample_weight = w)
  expect_equal(m_pos$W, m_named$W, tolerance = 1e-12)
  expect_true(isTRUE(m_pos$sample_weighted))
})

# ---- predict recenter ----

test_that("predict(recenter=TRUE) centers on the new data's own mean", {
  X <- make_structured(200, 6, seed = 7)
  m <- whiten_model(X, method = "ZCA", lambda = 0.1)
  Xt <- sweep(make_structured(80, 6, seed = 8), 2, rep(5, 6), "+")
  Z_default <- predict(m, Xt)
  Z_recenter <- predict(m, Xt, self_center = TRUE)
  expect_lt(mean(abs(colMeans(Z_recenter))), mean(abs(colMeans(Z_default))))
  expect_lt(mean(abs(colMeans(Z_recenter))), 1e-8)
  expect_error(predict(m, Xt, self_center = NA), "self_center")
})

# ---- whitening_patterns ----

test_that("whitening_patterns returns filters and patterns with correct shapes", {
  X <- make_structured(200, 8, seed = 9)
  m <- whiten_model(X, method = "PCA", n_comp = 4, lambda = 0.1)
  fp <- whitening_patterns(m)
  expect_equal(dim(fp$filters), c(4, 8))
  expect_equal(dim(fp$patterns), c(8, 4))
  # filters %*% patterns should be identity on the component space
  expect_equal(fp$filters %*% fp$patterns, diag(4), tolerance = 1e-8, ignore_attr = TRUE)
})

# ---- check_whitening diagnostics ----

test_that("check_whitening reports bounded whiteness and identity-zero logdet", {
  X <- make_structured(400, 6, seed = 10)
  m <- whiten_model(X, method = "ZCA", lambda = 0)
  d <- check_whitening(predict(m, X))
  expect_true(all(c("cov_dev_frob", "logdet", "whiteness") %in% names(d)))
  expect_gte(d$whiteness, 0)
  expect_lte(d$whiteness, 1)
  expect_gt(d$whiteness, 0.95)
  expect_equal(d$logdet, 0, tolerance = 1e-6)
})

# ---- forgetting-factor update ----

test_that("decay=1 reproduces the cumulative update", {
  set.seed(11)
  X1 <- matrix(rnorm(200 * 6), 200, 6)
  X2 <- matrix(rnorm(120 * 6), 120, 6)
  m <- whiten_model(X1, method = "ZCA", lambda = 0.1)
  m_a <- whiten_model_update(m, X2)
  m_b <- whiten_model_update(m, X2, decay = 1)
  expect_equal(m_a$W, m_b$W, tolerance = 1e-12)
})

test_that("decay<1 tracks a shifted regime better than cumulative", {
  set.seed(12)
  X0 <- matrix(rnorm(300 * 5), 300, 5)
  m0 <- whiten_model(X0, method = "ZCA", lambda = 0.02)
  Xd <- matrix(rnorm(300 * 5, sd = 3), 300, 5)
  m_cum <- whiten_model_update(m0, Xd, decay = 1)
  m_forg <- whiten_model_update(m0, Xd, decay = 0.3)
  dm_cum <- check_whitening(predict(m_cum, Xd))$diag_mean
  dm_forg <- check_whitening(predict(m_forg, Xd))$diag_mean
  expect_lt(abs(dm_forg - 1), abs(dm_cum - 1))
  expect_error(whiten_model_update(m0, Xd, decay = 0), "decay")
  expect_error(whiten_model_update(m0, Xd, decay = 1.5), "decay")
})

# ---- analytic tuning parity ----

test_that("fast_tune analytic path matches full-refit tuning on structured data", {
  X <- make_structured(200, 8, seed = 13)
  args <- list(
    X = X, methods = c("PCA", "ZCA", "SVD"),
    n_comp_grid = c(4, 6), var_threshold_grid = c(0.9, 0.95),
    lambda_grid = list("auto", 0, 0.05, 0.2), cv_folds = 5, seed = 1, top_n = 500
  )
  rf <- do.call(auto_tune_whitening, c(args, list(fast_tune = TRUE)))$ranking
  rs <- do.call(auto_tune_whitening, c(args, list(fast_tune = FALSE)))$ranking
  key <- function(r) paste(r$method, r$n_comp, r$var_threshold, r$lambda, r$cov_estimator)
  rf$k <- key(rf); rs$k <- key(rs)
  mrg <- merge(rf[, c("k", "mean_score")], rs[, c("k", "mean_score")],
               by = "k", suffixes = c(".f", ".s"))
  expect_gt(nrow(mrg), 0)
  expect_lt(max(abs(mrg$mean_score.f - mrg$mean_score.s)), 1e-8)
})

test_that("auto_tune accepts var_threshold_grid (regression test)", {
  X <- make_structured(120, 6, seed = 14)
  expect_no_error(
    auto_tune_whitening(X, methods = c("PCA"), var_threshold_grid = c(0.8, 0.9),
                        lambda_grid = list(0.1), cv_folds = 3, seed = 1)
  )
})

# ---- tangent space ----

test_that("epoch_covariances returns one SPD matrix per trial", {
  set.seed(15)
  Xt <- array(rnorm(20 * 6 * 64), dim = c(20, 6, 64))
  covs <- epoch_covariances(Xt, lambda = 0.05)
  expect_length(covs, 20)
  expect_equal(dim(covs[[1]]), c(6, 6))
  expect_true(all(eigen(covs[[1]], only.values = TRUE)$values > 0))
})

test_that("tangent_space round-trips and is isometric", {
  set.seed(16)
  Xt <- array(rnorm(25 * 6 * 80), dim = c(25, 6, 80))
  covs <- epoch_covariances(Xt, lambda = 0.05)
  V <- tangent_space(covs, mean_method = "logeuclid")
  expect_equal(dim(V), c(25, 6 * 7 / 2))
  covs_rec <- untangent_space(V)
  expect_equal(covs_rec[[1]], covs[[1]], tolerance = 1e-8)
  # Euclidean norm in tangent space equals Riemannian (affine) distance to ref
  ref <- attr(V, "reference")
  ri <- eegwhiten:::.invsqrtm_spd(ref)
  L1 <- eegwhiten:::.logm_spd(ri %*% covs[[1]] %*% ri)
  expect_equal(sqrt(sum(V[1, ]^2)), norm(L1, "F"), tolerance = 1e-8)
})

test_that("tangent_space accepts an explicit reference and a 3D array", {
  set.seed(17)
  Xt <- array(rnorm(10 * 5 * 64), dim = c(10, 5, 64))
  covs <- epoch_covariances(Xt)
  arr <- array(0, dim = c(10, 5, 5))
  for (i in seq_along(covs)) arr[i, , ] <- covs[[i]]
  V_list <- tangent_space(covs, ref = covs[[1]])
  V_arr <- tangent_space(arr, ref = covs[[1]])
  expect_equal(V_list, V_arr, tolerance = 1e-10, ignore_attr = TRUE)
})

# ---- recenter ----

test_that("recenter maps a matrix's own covariance to the identity", {
  X <- make_structured(200, 6, seed = 18)
  rc <- recenter(X)
  expect_equal(cov(rc$Z), diag(6), tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("whiten_batch mode='recenter' aligns each domain to the identity", {
  set.seed(19)
  Xl <- list(
    matrix(rnorm(150 * 5), 150, 5),
    matrix(rnorm(120 * 5), 120, 5) %*% diag(1:5)
  )
  out <- whiten_batch(Xl, mode = "recenter")
  expect_length(out, 2)
  for (o in out) {
    expect_equal(cov(o$Z), diag(5), tolerance = 1e-8, ignore_attr = TRUE)
    expect_identical(o$method, "recenter")
  }
  expect_identical(attr(out, "alignment"), "recenter")
})

# ---- relative / GED whitening ----

test_that("recenter guards rank-deficient covariance and tiny samples", {
  set.seed(24)
  # rank-deficient: fewer rows than channels
  Xrd <- matrix(rnorm(4 * 8), 4, 8)
  expect_error(recenter(Xrd), "positive-definite")
  # shrinkage rescues it
  expect_silent(rc <- recenter(Xrd, lambda = 0.3))
  expect_equal(ncol(rc$Z), 8)
  # single row without a reference is a clear error
  expect_error(recenter(matrix(rnorm(6), 1, 6)), "at least 2 rows")
  # with an explicit reference a single row is fine
  ref <- diag(6)
  expect_silent(recenter(matrix(rnorm(6), 1, 6), reference = ref))
})

test_that("tangent_space handles a single channel (p = 1) as list and array", {
  set.seed(25)
  covs <- lapply(1:7, function(i) matrix(runif(1, 0.5, 2), 1, 1))
  V <- tangent_space(covs, mean_method = "euclid")
  expect_equal(dim(V), c(7, 1))
  arr <- array(0, dim = c(7, 1, 1))
  for (i in 1:7) arr[i, , ] <- covs[[i]]
  V_arr <- tangent_space(arr, ref = attr(V, "reference"))
  V_ref <- tangent_space(covs, ref = attr(V, "reference"))
  expect_equal(V_arr, V_ref, tolerance = 1e-10, ignore_attr = TRUE)
  rec <- untangent_space(V)
  expect_equal(as.numeric(rec[[1]]), as.numeric(covs[[1]]), tolerance = 1e-8)
})

test_that("epoch_covariances warns when rank-deficient without shrinkage", {
  set.seed(26)
  Xt <- array(rnorm(5 * 8 * 6), dim = c(5, 8, 6))  # n_time(6) < n_channels(8)
  expect_warning(epoch_covariances(Xt), "rank-deficient")
  expect_silent(epoch_covariances(Xt, lambda = 0.1))
})

test_that("fast_tune stays consistent with full refit even with lambda = 1 in the grid", {
  X <- make_structured(160, 6, seed = 27)
  args <- list(X = X, methods = c("PCA", "ZCA"), n_comp_grid = c(3),
               lambda_grid = list(0, 0.1, 1), cv_folds = 4, seed = 1, top_n = 500)
  rf <- do.call(auto_tune_whitening, c(args, list(fast_tune = TRUE)))$ranking
  rs <- do.call(auto_tune_whitening, c(args, list(fast_tune = FALSE)))$ranking
  key <- function(r) paste(r$method, r$n_comp, r$var_threshold, r$lambda, r$cov_estimator)
  rf$k <- key(rf); rs$k <- key(rs)
  mrg <- merge(rf[, c("k", "mean_score")], rs[, c("k", "mean_score")],
               by = "k", suffixes = c(".f", ".s"))
  expect_lt(max(abs(mrg$mean_score.f - mrg$mean_score.s)), 1e-8)
})

test_that("whiten_relative whitens the reference and diagonalizes the target", {
  set.seed(20)
  p <- 6
  Sa <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p)
  Sb <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p)
  gr <- whiten_relative(Sb, Sa)
  expect_equal(gr$W %*% Sa %*% t(gr$W), diag(p), tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(gr$W %*% Sb %*% t(gr$W), diag(gr$values), tolerance = 1e-8, ignore_attr = TRUE)
  expect_false(is.unsorted(rev(gr$values)))
  expect_equal(dim(whiten_relative(Sb, Sa, n_comp = 3)$W), c(3, p))
})

# ---- T1 correctness: Tyler scale ----

test_that("Tyler covariance is scale-restored: Gaussian whitening ~ identity", {
  set.seed(40)
  X <- matrix(rnorm(500 * 6), 500, 6)
  m <- whiten_model(X, method = "ZCA", lambda = 0, cov_estimator = "tyler")
  dm <- check_whitening(predict(m, X))$diag_mean
  expect_equal(dm, 1, tolerance = 0.15)        # was ~9.9 before the scale fix
  # heavy-tailed: scale is O(1), not grossly off
  Xt <- matrix(rt(500 * 6, df = 3), 500, 6)
  mt <- whiten_model(Xt, method = "ZCA", lambda = 0, cov_estimator = "tyler")
  expect_lt(check_whitening(predict(mt, Xt))$diag_mean, 3)
})

# ---- T1 correctness: scale-invariant positive-definiteness ----

test_that("PD check is scale-invariant (tiny-magnitude well-conditioned data accepted)", {
  # A well-conditioned (cond = 2) matrix at tiny magnitude must NOT be rejected.
  S_tiny <- diag(c(1e-10, 2e-10, 1.5e-10))
  expect_silent(W <- whiten_fit(S_tiny, method = "ZCA"))
  expect_equal(dim(W), c(3, 3))

  # EEG-in-volts magnitude whitens to identity just like the unit-scale data.
  set.seed(41)
  X <- matrix(rnorm(300 * 6), 300, 6)
  m_unit <- whiten_model(X, method = "ZCA", lambda = 0)
  m_volt <- whiten_model(X * 1e-6, method = "ZCA", lambda = 0)
  expect_equal(check_whitening(predict(m_unit, X))$diag_mean, 1, tolerance = 1e-8)
  expect_equal(check_whitening(predict(m_volt, X * 1e-6))$diag_mean, 1, tolerance = 1e-8)

  # Genuinely singular / non-PD still errors.
  expect_error(whiten_fit(matrix(1, 4, 4), method = "ZCA"), "positive-definite")

  # Correlation-based methods tolerate extreme heterogeneous channel scales.
  Xs <- matrix(rnorm(300 * 4), 300, 4) %*% diag(c(1e-6, 1, 1e3, 1e6))
  expect_silent(whiten_model(Xs, method = "ZCA-cor", lambda = 0))
})

# ---- T1 test strengthening: rspectra backend, supervised tuning ----

test_that("eig_method='rspectra' agrees with 'base'", {
  skip_if_not_installed("RSpectra")
  X <- make_structured(300, 12, seed = 42)
  m_base <- whiten_model(X, method = "PCA", n_comp = 4, lambda = 0.05, eig_method = "base")
  m_rs   <- whiten_model(X, method = "PCA", n_comp = 4, lambda = 0.05, eig_method = "rspectra")
  # cov(Z) is invariant to component sign, so compare whitening quality
  expect_equal(check_whitening(predict(m_base, X))$cov_matrix,
               check_whitening(predict(m_rs, X))$cov_matrix, tolerance = 1e-6)
})

test_that("supervised tuning (scoring='accuracy') runs and selects a valid model", {
  set.seed(43)
  n <- 120
  y <- factor(rep(c("a", "b"), each = n / 2))
  X <- rbind(
    matrix(rnorm(n / 2 * 6, mean = 0), n / 2, 6),
    matrix(rnorm(n / 2 * 6, mean = 1.2), n / 2, 6)
  )
  colnames(X) <- paste0("Ch", 1:6)
  tuned <- auto_tune_whitening(X, y = y, methods = c("PCA", "ZCA"),
                               n_comp_grid = c(3, 5), lambda_grid = list(0, 0.1),
                               cv_folds = 4, seed = 1, scoring = "accuracy")
  expect_s3_class(tuned$best_model, "whiten_model")
  expect_true(all(tuned$ranking$mean_score >= 0 & tuned$ranking$mean_score <= 1))
  # custom score_fn hook works
  tuned2 <- auto_tune_whitening(X, y = y, methods = "ZCA", lambda_grid = list(0.1),
                                cv_folds = 3, seed = 1, scoring = "accuracy",
                                score_fn = function(Ztr, ytr, Zva, yva) {
                                  mean(yva == ytr[1])  # trivial baseline scorer
                                })
  expect_true(is.finite(tuned2$ranking$mean_score[1]))
})

# ---- T2 API consistency: S3 classes, methods, deprecation shims ----

test_that("recenter/whiten_relative/whiten_tune/ea_model are classed with methods", {
  set.seed(50)
  X <- make_structured(200, 6, seed = 50)
  rc <- recenter(X)
  expect_s3_class(rc, "recenter")
  expect_output(print(rc), "<recenter>")
  Xnew <- make_structured(40, 6, seed = 51)
  expect_equal(dim(predict(rc, Xnew)), c(40, 6))

  gr <- whiten_relative(cov(X), cov(make_structured(200, 6, seed = 52)))
  expect_s3_class(gr, "relative_whiten")
  expect_output(print(gr), "relative_whiten")
  expect_equal(dim(predict(gr, Xnew)), c(40, 6))

  tuned <- auto_tune_whitening(X, methods = c("PCA", "ZCA"), n_comp_grid = c(3),
                               lambda_grid = list(0, 0.1), cv_folds = 3, seed = 1)
  expect_s3_class(tuned, "whiten_tune")
  expect_output(print(tuned), "<whiten_tune>")

  ea <- euclidean_alignment(list(X, Xnew), input = "raw")
  expect_output(print(ea$model), "<ea_model>")
  expect_output(summary(ea$model), "ref cond")
})

test_that("deprecated argument order and names still work with a warning", {
  set.seed(53)
  X <- make_structured(150, 6, seed = 53)
  m <- whiten_model(X, method = "ZCA", lambda = 0)
  Z <- predict(m, X)
  # old data-first order
  expect_warning(rec_old <- unwhiten(Z, m), "deprecated")
  rec_new <- unwhiten(m, Z)
  expect_equal(rec_old, rec_new, tolerance = 1e-10)
  expect_warning(unwhiten_fast(Z, m), "deprecated")

  # old shrinkage name
  Xt <- array(rnorm(10 * 6 * 80), dim = c(10, 6, 80))
  expect_warning(epoch_covariances(Xt, shrinkage = 0.1), "deprecated")

  # old predict recenter arg
  expect_warning(predict(m, X, recenter = TRUE), "deprecated")

  # inverse_ea old order
  ea <- euclidean_alignment(list(X, make_structured(120, 6, seed = 54)), input = "raw")
  Za <- predict(ea$model, X)
  expect_warning(inverse_ea(Za, ea$model), "deprecated")
})

# ---- cleanup: predict.whiten_tune, cor-method scale-invariance ----

test_that("predict.whiten_tune delegates to the best model", {
  X <- make_structured(120, 6, seed = 60)
  tuned <- auto_tune_whitening(X, methods = c("PCA", "ZCA"), n_comp_grid = c(3),
                               lambda_grid = list(0, 0.1), cv_folds = 3, seed = 1)
  expect_equal(predict(tuned, X), predict(tuned$best_model, X), tolerance = 1e-12)
})

test_that("cor-method PD check is scale-invariant in whiten_fit and standalone", {
  # Heterogeneous channel scales: covariance is ill-conditioned but the
  # correlation matrix is well-conditioned, so cor-methods must succeed.
  set.seed(61)
  Xs <- matrix(rnorm(300 * 4), 300, 4) %*% diag(c(1e-6, 1, 1e3, 1e6))
  S <- cov(Xs)
  expect_silent(whiten_fit(S, method = "ZCA-cor"))
  expect_silent(whiten_fit(S, method = "PCA-cor", n_comp = 3))
  expect_silent(ZCA_cor(S, returnW = TRUE, PhiPsi = FALSE))
  expect_silent(PCA_cor(S, returnW = TRUE, PhiPsi = FALSE))
})

# ---- cheap high-leverage additions: diagonal auto, epoch auto, riemann primitives ----

test_that("lambda='auto' works with the diagonal target (Schaefer-Strimmer)", {
  set.seed(70)
  X <- matrix(rnorm(150 * 6), 150, 6) %*% chol(toeplitz(0.5^(0:5)))
  X <- X %*% diag(c(1e-2, 1, 10, 1, 5, 1))   # heterogeneous scales
  m <- whiten_model(X, method = "ZCA", lambda = "auto", shrink_target = "diagonal")
  expect_true(m$lambda >= 0 && m$lambda <= 1)
  # diagonal target preserves the per-channel variances
  expect_equal(diag(m$cov), diag(cov(scale(X, scale = FALSE))),
               tolerance = 1e-8, ignore_attr = TRUE)
  # the chosen intensity matches the standalone estimator
  Xc <- scale(X, scale = FALSE)
  expect_equal(m$lambda, eegwhiten:::.lambda_ss(Xc), tolerance = 1e-10)
})

test_that("epoch_covariances(lambda='auto') yields SPD per-epoch covariances", {
  set.seed(71)
  Xt <- array(rnorm(12 * 8 * 20), dim = c(12, 8, 20))  # n_time ~ n_channels
  covs <- epoch_covariances(Xt, lambda = "auto")
  expect_length(covs, 12)
  expect_true(all(vapply(covs, function(C) min(eigen(C, only.values = TRUE)$values) > 0, logical(1))))
})

test_that("riemann_distance has metric properties", {
  set.seed(72)
  A <- cov(matrix(rnorm(300 * 5), 300, 5))
  B <- cov(matrix(rnorm(300 * 5), 300, 5))
  C <- cov(matrix(rnorm(300 * 5), 300, 5))
  expect_equal(riemann_distance(A, A), 0, tolerance = 1e-8)
  expect_gt(riemann_distance(A, B), 0)
  expect_equal(riemann_distance(A, B), riemann_distance(B, A), tolerance = 1e-10)
  # triangle inequality
  expect_lte(riemann_distance(A, C), riemann_distance(A, B) + riemann_distance(B, C) + 1e-8)
  # affine invariance: d(A,B) == d(GAG', GBG')
  G <- matrix(rnorm(25), 5, 5)
  expect_equal(riemann_distance(A, B),
               riemann_distance(G %*% A %*% t(G), G %*% B %*% t(G)), tolerance = 1e-6)
  # euclid metric equals Frobenius distance
  expect_equal(riemann_distance(A, B, metric = "euclid"), norm(A - B, "F"), tolerance = 1e-10)
})

test_that("spd_mean: euclid is arithmetic; riemann is the Karcher centroid", {
  set.seed(73)
  covs <- lapply(1:6, function(i) cov(matrix(rnorm(200 * 4), 200, 4)))
  expect_equal(spd_mean(covs, "euclid"), Reduce(`+`, covs) / 6, tolerance = 1e-10)
  M <- spd_mean(covs, "riemann")
  expect_true(min(eigen(M, only.values = TRUE)$values) > 0)
  # Karcher mean: sum of tangent log-maps at M is ~ 0
  Mi <- eegwhiten:::.invsqrtm_spd(M)
  S <- Reduce(`+`, lapply(covs, function(C) eegwhiten:::.logm_spd(Mi %*% C %*% Mi)))
  expect_lt(norm(S, "F"), 1e-4)
})

test_that("barycenter_whitener maps the barycenter to identity and inverts", {
  set.seed(74)
  Xt <- array(rnorm(25 * 5 * 200), dim = c(25, 5, 200))
  covs <- epoch_covariances(Xt, lambda = 0.02)
  bw <- barycenter_whitener(covs, metric = "logeuclid")
  expect_s3_class(bw, "ea_model")
  # W M W' == I, i.e. the barycenter is whitened to identity
  expect_equal(bw$W %*% bw$reference_cov %*% t(bw$W), diag(5), tolerance = 1e-8, ignore_attr = TRUE)
  Xnew <- matrix(rnorm(80 * 5), 80, 5)
  Z <- predict(bw, Xnew)
  expect_equal(dim(Z), c(80, 5))
  expect_equal(inverse_ea(bw, Z), Xnew, tolerance = 1e-8, ignore_attr = TRUE)
})

# ---- bug-hunt regressions ----

test_that("lambda_method label reflects the estimator actually used", {
  X <- make_structured(120, 5, seed = 80)
  m_diag <- whiten_model(X, method = "ZCA", lambda = "auto", shrink_target = "diagonal",
                         lambda_method = "oas")
  expect_identical(m_diag$lambda_method, "ss")               # not the ignored 'oas'
  expect_equal(m_diag$lambda, eegwhiten:::.lambda_ss(scale(X, scale = FALSE)),
               tolerance = 1e-10)
  m_id <- whiten_model(X, method = "ZCA", lambda = "auto", lambda_method = "lw")
  expect_identical(m_id$lambda_method, "lw")                 # identity target keeps the choice
})

test_that("whiten_batch(mode='ea') reports the centering vector actually used", {
  set.seed(81)
  X1 <- matrix(rnorm(200 * 5, mean = 3), 200, 5)
  X2 <- matrix(rnorm(150 * 5, mean = -2), 150, 5)
  o <- whiten_batch(list(X1, X2), mode = "ea", center = TRUE)
  expect_equal(o[[1]]$center, colMeans(X1), tolerance = 1e-10)
  expect_equal(o[[2]]$center, colMeans(X2), tolerance = 1e-10)
  # reconstruction using the reported (Z, W, center) recovers the data
  rec <- o[[1]]$Z %*% solve(o[[1]]$W) +
    matrix(o[[1]]$center, nrow(o[[1]]$Z), 5, byrow = TRUE)
  expect_lt(max(abs(rec - X1)), 1e-8)
  # center = FALSE still reports zeros; covariance input reports NULL
  expect_true(all(whiten_batch(list(X1), mode = "ea", center = FALSE)[[1]]$center == 0))
  expect_null(whiten_batch(list(cov(X1)), mode = "ea", ea_input = "cov")[[1]]$center)
})
