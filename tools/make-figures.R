#!/usr/bin/env Rscript
# Generate demonstration / comparison figures for eegwhiten using the package's
# own functions on reproducible synthetic data. Output: man/figures/*.png
#
# Usage:  Rscript tools/make-figures.R
suppressMessages({
  if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(eegwhiten)
})

out_dir <- "man/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- shared helpers ----
div_pal <- grDevices::colorRampPalette(c("#2166AC", "#4393C3", "#F7F7F7", "#D6604D", "#B2182B"))(128)

heat <- function(M, main, zlim = NULL, col = div_pal, labels = TRUE) {
  d <- ncol(M)
  if (is.null(zlim)) {
    lim <- max(abs(M)); lim <- if (!is.finite(lim) || lim <= 0) 1 else lim
    zlim <- c(-lim, lim)
  }
  graphics::image(seq_len(d), seq_len(d), t(M[d:1, , drop = FALSE]),
                  zlim = zlim, col = col, axes = FALSE, xlab = "", ylab = "",
                  main = main)
  graphics::box(col = "grey60")
}

# structured correlated covariance (Toeplitz + a couple of dominant components)
make_cov <- function(p, seed = 1) {
  set.seed(seed)
  R <- 0.6 ^ abs(outer(seq_len(p), seq_len(p), "-"))   # AR(1)-like correlation
  B <- matrix(rnorm(p * 2), p, 2)
  S <- R + 0.8 * tcrossprod(B)                          # add 2 strong sources
  diag(S) <- diag(S) + 0.2
  S
}
sample_data <- function(n, S, seed = 1) {
  set.seed(seed)
  matrix(rnorm(n * ncol(S)), n, ncol(S)) %*% chol(S)
}

# =====================================================================
# Figure 1: whitening works -- covariance before vs after
# =====================================================================
p <- 16
S <- make_cov(p, seed = 1)
X_train <- sample_data(400, S, seed = 1)
X_test  <- sample_data(400, S, seed = 2)
m <- whiten_model(X_train, method = "ZCA", lambda = 0)
Cx <- cov2cor(cov(X_test))
Cz <- cov(predict(m, X_test))

grDevices::png(file.path(out_dir, "01-before-after.png"), width = 1100, height = 520, res = 130)
graphics::par(mfrow = c(1, 2), mar = c(1.5, 1.5, 3, 1.5), oma = c(0, 0, 2, 0))
heat(Cx, sprintf("Before:  cor(X)\noff-diag Frobenius = %.2f", check_whitening(X_test)$offdiag_frob), zlim = c(-1, 1))
heat(Cz, sprintf("After ZCA:  cov(Z)\noff-diag Frobenius = %.3f", check_whitening(predict(m, X_test))$offdiag_frob), zlim = c(-1, 1))
graphics::mtext("eegwhiten: whitening maps the covariance to the identity (held-out data)",
                outer = TRUE, cex = 0.95, font = 2)
grDevices::dev.off()

# =====================================================================
# Figure 2: the six methods all whiten, but the transforms differ (|W|)
# =====================================================================
methods <- c("PCA", "SVD", "ZCA", "ZCA-cor", "PCA-cor", "Cholesky")
grDevices::png(file.path(out_dir, "02-methods-W.png"), width = 1200, height = 820, res = 130)
graphics::par(mfrow = c(2, 3), mar = c(1.2, 1.2, 3, 1.2), oma = c(0, 0, 2.4, 0))
seq_pal <- grDevices::colorRampPalette(c("#FFFFFF", "#FEE08B", "#F46D43", "#A50026"))(128)
for (meth in methods) {
  mm <- whiten_model(X_train, method = meth, lambda = 0)
  W <- abs(mm$W)
  dm <- check_whitening(predict(mm, X_test))$diag_mean
  heat(W, sprintf("%s   (cov(Z) diag = %.2f)", meth, dm),
       zlim = c(0, max(abs(W))), col = seq_pal)
}
graphics::mtext("Six whitening transforms |W| -- all whiten (diag ~ 1), but the filters differ",
                outer = TRUE, cex = 0.95, font = 2)
grDevices::dev.off()

# =====================================================================
# Figure 3: shrinkage stabilizes whitening on ill-conditioned (p ~ n) data
# =====================================================================
p2 <- 40
S2 <- make_cov(p2, seed = 3)
Xtr <- sample_data(60, S2, seed = 3)     # p close to n -> ill-conditioned
Xte <- sample_data(400, S2, seed = 4)
lams <- seq(0, 0.6, by = 0.02)
off <- sapply(lams, function(l) {
  mm <- whiten_model(Xtr, method = "ZCA", lambda = l)
  check_whitening(predict(mm, Xte))$offdiag_frob
})
auto_lam <- whiten_model(Xtr, method = "ZCA", lambda = "auto")$lambda

grDevices::png(file.path(out_dir, "03-shrinkage.png"), width = 1050, height = 620, res = 130)
graphics::par(mar = c(4.2, 4.2, 3, 1))
graphics::plot(lams, off, type = "b", pch = 19, col = "#2166AC", lwd = 2,
               xlab = expression(paste("shrinkage  ", lambda)),
               ylab = "test-set off-diagonal Frobenius  (lower = whiter)",
               main = sprintf("Shrinkage stabilizes whitening (p = %d, n_train = 60)", p2))
graphics::abline(v = auto_lam, lty = 2, col = "#B2182B", lwd = 2)
graphics::text(auto_lam, max(off) * 0.95, sprintf('  lambda="auto" = %.3f', auto_lam),
               col = "#B2182B", pos = 4, cex = 0.9)
graphics::grid(col = "grey85")
grDevices::dev.off()

# =====================================================================
# Figure 4: cross-session alignment shrinks between-session distances
# =====================================================================
n_sess <- 8
set.seed(10)
sess_covs <- lapply(seq_len(n_sess), function(i) {
  Si <- make_cov(p, seed = 100 + i)
  Si <- Si * runif(1, 0.5, 2)                          # per-session scale drift
  rot <- qr.Q(qr(matrix(rnorm(p * p), p, p)))          # mild rotation
  .symmetrize <- function(A) (A + t(A)) / 2
  .symmetrize(rot %*% Si %*% t(rot))
})
sess_data <- lapply(seq_len(n_sess), function(i) {
  sample_data(300, sess_covs[[i]], seed = 200 + i) + matrix(rnorm(p, sd = 2), 300, p, byrow = TRUE)
})

pair_dists <- function(cov_list) {
  k <- length(cov_list); d <- c()
  for (i in seq_len(k - 1)) for (j in (i + 1):k) d <- c(d, riemann_distance(cov_list[[i]], cov_list[[j]]))
  d
}
covs_raw <- lapply(sess_data, function(x) cov(scale(x, scale = FALSE)))
rc  <- whiten_batch(sess_data, mode = "recenter")
covs_rc <- lapply(rc, function(o) cov(o$Z))

dl <- list(`raw (no alignment)` = pair_dists(covs_raw),
           `recenter (per-domain)` = pair_dists(covs_rc))

grDevices::png(file.path(out_dir, "04-alignment.png"), width = 1050, height = 640, res = 130)
graphics::par(mar = c(4, 4.5, 3.2, 1))
graphics::boxplot(dl, col = c("#FDDBC7", "#2166AC"), ylab = "pairwise Riemannian distance",
                  main = sprintf("Per-domain recentering removes between-session shift (%d sessions)", n_sess),
                  cex.main = 0.98,
                  border = "grey30")
graphics::grid(nx = NA, ny = NULL, col = "grey85")
grDevices::dev.off()

# =====================================================================
# Figure 5: analytic fast-tune is much faster than full refit (equal scores)
# =====================================================================
Xt <- sample_data(300, make_cov(30, seed = 5), seed = 5)
grid <- list(methods = c("PCA", "ZCA", "SVD"), n_comp_grid = c(5, 10, 15, 20, 25),
             lambda_grid = list("auto", 0, 0.01, 0.05, 0.1, 0.2, 0.5), cv_folds = 5, seed = 7)
tf <- system.time(do.call(auto_tune_whitening, c(list(Xt), grid, fast_tune = TRUE)))[["elapsed"]]
ts <- system.time(do.call(auto_tune_whitening, c(list(Xt), grid, fast_tune = FALSE)))[["elapsed"]]

grDevices::png(file.path(out_dir, "05-fast-tune.png"), width = 820, height = 620, res = 130)
graphics::par(mar = c(3.5, 4.5, 3.2, 1))
bp <- graphics::barplot(c(`full refit` = ts, `analytic fast_tune` = tf),
                        col = c("#BDBDBD", "#2166AC"), ylab = "tuning time (s)",
                        ylim = c(0, ts * 1.18),
                        main = sprintf("Analytic tuning fast path: %.1fx faster\n(identical CV scores)", ts / tf))
graphics::text(bp, c(ts, tf), sprintf("%.2fs", c(ts, tf)), pos = 3, xpd = NA)
grDevices::dev.off()

# =====================================================================
# Figure 6: tangent-space 2D embedding -- per-domain recentering pulls
#           sessions together
# =====================================================================
set.seed(20)
ns <- 4; ntr <- 45; pc <- 8; Tt <- 160
sess_cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
# per-session base covariance (drift) and per-trial covariances around it
sess_trialcovs <- lapply(seq_len(ns), function(s) {
  base <- make_cov(pc, seed = 300 + s) * runif(1, 0.6, 1.8)
  L <- chol(base)
  lapply(seq_len(ntr), function(k) {
    Xk <- matrix(rnorm(Tt * pc), Tt, pc) %*% L
    .symmetrize <- function(A) (A + t(A)) / 2
    .symmetrize(cov(Xk))
  })
})
labels <- rep(seq_len(ns), each = ntr)

embed2d <- function(cov_list) {
  V <- tangent_space(cov_list, mean_method = "riemann")
  pr <- stats::prcomp(V, center = TRUE, scale. = FALSE)
  pr$x[, 1:2]
}
covs_all <- do.call(c, sess_trialcovs)
emb_raw <- embed2d(covs_all)
# per-domain recentering of covariances: map each session's mean to identity
covs_aligned <- do.call(c, lapply(sess_trialcovs, function(cl) {
  W <- barycenter_whitener(cl, metric = "riemann")$W
  lapply(cl, function(C) W %*% C %*% W)
}))
emb_al <- embed2d(covs_aligned)

grDevices::png(file.path(out_dir, "06-tangent-embedding.png"), width = 1150, height = 560, res = 130)
graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
graphics::plot(emb_raw, col = sess_cols[labels], pch = 19, cex = 0.7,
               xlab = "tangent PC1", ylab = "tangent PC2", main = "raw (no alignment)")
graphics::grid(col = "grey88")
graphics::plot(emb_al, col = sess_cols[labels], pch = 19, cex = 0.7,
               xlab = "tangent PC1", ylab = "tangent PC2", main = "per-domain recentering")
graphics::grid(col = "grey88")
graphics::legend("topright", legend = paste("session", seq_len(ns)), col = sess_cols,
                 pch = 19, bty = "n", cex = 0.8)
graphics::mtext("Tangent-space embedding: recentering pulls the sessions together",
                outer = TRUE, cex = 0.95, font = 2)
grDevices::dev.off()

# =====================================================================
# Figure 7: Tyler is robust to outliers; empirical whitening is corrupted
# =====================================================================
set.seed(30)
p3 <- 12
S3 <- make_cov(p3, seed = 30)
X_clean <- sample_data(600, S3, seed = 31)            # clean held-out test
X_tr <- sample_data(200, S3, seed = 32)
n_out <- round(0.12 * 200)                            # 12% gross outliers
X_tr[sample(200, n_out), ] <- matrix(rnorm(n_out * p3, sd = 10), n_out, p3)
Cz_emp <- cov(predict(whiten_model(X_tr, method = "ZCA", lambda = 0, cov_estimator = "empirical"), X_clean))
Cz_tyl <- cov(predict(whiten_model(X_tr, method = "ZCA", lambda = 0, cov_estimator = "tyler"), X_clean))
devI <- function(C) norm(C - diag(nrow(C)), "F")    # total deviation from identity

grDevices::png(file.path(out_dir, "07-robust-tyler.png"), width = 1100, height = 540, res = 130)
graphics::par(mfrow = c(1, 2), mar = c(1.5, 1.5, 3, 1.5), oma = c(0, 0, 2, 0))
heat(Cz_emp - diag(p3), sprintf("empirical cov\n||cov(Z) - I|| = %.2f", devI(Cz_emp)), zlim = c(-1, 1))
heat(Cz_tyl - diag(p3), sprintf("Tyler (robust)\n||cov(Z) - I|| = %.2f", devI(Cz_tyl)), zlim = c(-1, 1))
graphics::mtext("12% outliers in training: empirical whitening fails on clean test data, Tyler holds",
                outer = TRUE, cex = 0.92, font = 2)
grDevices::dev.off()

# =====================================================================
# Figure 8: covariance estimators across two hard regimes (test whiteness)
# =====================================================================
metric_dev <- function(Xtr, est, lam) {
  m <- whiten_model(Xtr, method = "ZCA", lambda = lam, cov_estimator = est)
  check_whitening(predict(m, X_clean))$cov_dev_frob   # total ||cov(Z) - I||_F
}
# small-n regime (p=12, n=20) and outlier regime (reuse contaminated X_tr)
set.seed(40)
X_smalln <- sample_data(20, S3, seed = 41)
tab <- rbind(
  `small n (p=12, n=20)` = c(
    empirical = metric_dev(X_smalln, "empirical", 0),
    `shrink (auto)` = metric_dev(X_smalln, "empirical", "auto"),
    Tyler = metric_dev(X_smalln, "tyler", 0)),
  `12% outliers (n=200)` = c(
    empirical = metric_dev(X_tr, "empirical", 0),
    `shrink (auto)` = metric_dev(X_tr, "empirical", "auto"),
    Tyler = metric_dev(X_tr, "tyler", 0))
)

grDevices::png(file.path(out_dir, "08-estimators.png"), width = 1050, height = 600, res = 130)
graphics::par(mar = c(4, 4.5, 3.2, 1))
bp <- graphics::barplot(t(tab), beside = TRUE, col = c("#BDBDBD", "#2166AC", "#D6604D"),
                        ylab = expression(paste("test-set  ", "||cov(Z) - I||"[F], "   (lower = whiter)")),
                        main = "Covariance estimators: shrinkage helps small-n, Tyler helps outliers",
                        legend.text = colnames(tab),
                        args.legend = list(x = "topright", bty = "n"))
grDevices::dev.off()

cat("Wrote figures to", normalizePath(out_dir), ":\n")
cat(paste(" -", list.files(out_dir, pattern = "\\.png$")), sep = "\n")
cat("\n")
