
############################################################
# FILENAME: utils_whiten.R
# DESCRIPTION: Internal utilities
# UPDATES: Added check_condition
############################################################

# Set dimnames and method attribute on a matrix
.set_matrix_attr <- function(x, row_names, col_names, method) {
  if (!is.null(row_names)) rownames(x) <- row_names
  if (!is.null(col_names)) colnames(x) <- col_names
  attr(x, "method") <- method
  x
}

# Validate a probability-like scalar in [0, 1]
.check_unit_interval <- function(x, name, allow_char_auto = FALSE) {
  if (allow_char_auto && is.character(x)) {
    if (length(x) != 1L || !(x %in% "auto")) {
      stop(sprintf("%s must be numeric in [0, 1] or the string 'auto'.", name))
    }
    return(invisible(x))
  }

  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop(sprintf("%s must be a single finite numeric value in [0, 1].", name))
  }
  if (x < 0 || x > 1) {
    stop(sprintf("%s must be between 0 and 1.", name))
  }
  invisible(x)
}

# Validate cumulative explained-variance threshold in (0, 1]
.check_var_threshold <- function(var_threshold) {
  if (is.null(var_threshold)) return(invisible(NULL))
  if (!is.numeric(var_threshold) || length(var_threshold) != 1L || !is.finite(var_threshold)) {
    stop("var_threshold must be a single finite numeric value in (0, 1].")
  }
  if (var_threshold <= 0 || var_threshold > 1) {
    stop("var_threshold must be in (0, 1].")
  }
  invisible(var_threshold)
}

# Canonicalize user-facing whitening method names
.normalize_whiten_method <- function(method) {
  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    stop("method must be a single non-missing character string.")
  }

  aliases <- c(
    "SVD" = "SVD",
    "svd" = "SVD",
    "ZCA" = "ZCA",
    "zca" = "ZCA",
    "ZCA-cor" = "ZCA-cor",
    "ZCA_cor" = "ZCA-cor",
    "zca-cor" = "ZCA-cor",
    "zca_cor" = "ZCA-cor",
    "PCA" = "PCA",
    "pca" = "PCA",
    "PCA-cor" = "PCA-cor",
    "PCA_cor" = "PCA-cor",
    "pca-cor" = "PCA-cor",
    "pca_cor" = "PCA-cor",
    "Cholesky" = "Cholesky",
    "cholesky" = "Cholesky"
  )

  key <- trimws(method)
  if (!key %in% names(aliases)) {
    stop("Unknown method. Use one of: SVD, ZCA, ZCA-cor, PCA, PCA-cor, Cholesky.")
  }
  aliases[[key]]
}

# Map canonical method names to concrete function names
.method_to_function <- function(method) {
  method <- .normalize_whiten_method(method)
  switch(
    method,
    "PCA-cor" = "PCA_cor",
    "ZCA-cor" = "ZCA_cor",
    method
  )
}

# Methods that support dimensionality reduction
.is_reduction_method <- function(method) {
  .normalize_whiten_method(method) %in% c("PCA", "PCA-cor", "SVD")
}

# Unified SPD eigendecomposition helper with optional top-k solver.
# precomp may carry a previously computed full decomposition (values, vectors)
# of the same matrix; when it holds at least k components it is reused directly,
# which avoids recomputing the eigendecomposition inside whitening methods.
.eigen_spd <- function(Sigma,
                       k = NULL,
                       method = c("auto", "base", "rspectra"),
                       fast = FALSE,
                       precomp = NULL) {
  method <- match.arg(method)
  if (!is.matrix(Sigma) || !is.numeric(Sigma)) {
    stop("Sigma must be a numeric matrix.")
  }
  d <- ncol(Sigma)
  if (nrow(Sigma) != d) {
    stop("Sigma must be square.")
  }
  if (is.null(k)) {
    k <- d
  } else {
    if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 1 || k > d || k != as.integer(k)) {
      stop("k must be an integer between 1 and ncol(Sigma).")
    }
    k <- as.integer(k)
  }

  if (!is.null(precomp) && !is.null(precomp$values) && !is.null(precomp$vectors) &&
      length(precomp$values) >= k && ncol(precomp$vectors) >= k &&
      nrow(precomp$vectors) == d) {
    return(list(
      values = precomp$values[seq_len(k)],
      vectors = precomp$vectors[, seq_len(k), drop = FALSE],
      engine = if (!is.null(precomp$engine)) precomp$engine else "precomp"
    ))
  }

  rs_available <- requireNamespace("RSpectra", quietly = TRUE)
  use_rs <- FALSE
  if (method == "rspectra") {
    if (!rs_available) {
      stop("eig_method = 'rspectra' requires package 'RSpectra'.")
    }
    use_rs <- k < d
  } else if (method == "auto") {
    use_rs <- rs_available && k < d
  }

  if (use_rs) {
    opts <- if (fast) list(tol = 1e-4, maxitr = 1000) else list()
    rs <- RSpectra::eigs_sym(
      .symmetrize(Sigma),
      k = k,
      which = "LM",
      opts = opts
    )
    ord <- order(Re(rs$values), decreasing = TRUE)
    values <- Re(rs$values[ord])
    vectors <- Re(rs$vectors[, ord, drop = FALSE])
    return(list(values = values, vectors = vectors, engine = "RSpectra"))
  }

  eig <- eigen(Sigma, symmetric = TRUE)
  values <- eig$values
  vectors <- eig$vectors
  if (k < d) {
    values <- values[seq_len(k)]
    vectors <- vectors[, seq_len(k), drop = FALSE]
  }
  list(values = values, vectors = vectors, engine = "eigen")
}

# Validate optional sample weights
.validate_sample_weight <- function(sample_weight, n, allow_null = TRUE) {
  if (is.null(sample_weight)) {
    if (allow_null) return(NULL)
    stop("sample_weight must not be NULL.")
  }
  if (!is.numeric(sample_weight) || length(sample_weight) != n || any(!is.finite(sample_weight))) {
    stop("sample_weight must be a finite numeric vector with length nrow(X).")
  }
  if (any(sample_weight < 0)) {
    stop("sample_weight must be non-negative.")
  }
  if (sum(sample_weight) <= .Machine$double.eps) {
    stop("sum(sample_weight) must be positive.")
  }
  as.numeric(sample_weight)
}

# Weighted first/second moments
.weighted_moments <- function(X, sample_weight = NULL) {
  n <- nrow(X)
  w <- .validate_sample_weight(sample_weight, n, allow_null = TRUE)
  if (is.null(w)) {
    w <- rep(1, n)
  }
  sw <- sum(w)
  sw2 <- sum(w ^ 2)
  xw <- X * w
  list(
    sum_w = sw,
    sum_w2 = sw2,
    sum_x = colSums(xw),
    sum_xx = crossprod(X, xw)
  )
}

# Covariance from moments with optional centering
.cov_from_moments <- function(sum_w, sum_w2, sum_x, sum_xx, center = TRUE) {
  d <- length(sum_x)
  if (center) {
    mu <- sum_x / sum_w
    scatter <- sum_xx - sum_w * tcrossprod(mu)
  } else {
    mu <- rep(0, d)
    scatter <- sum_xx
  }
  denom <- sum_w - (sum_w2 / sum_w)
  if (!is.finite(denom) || denom <= .Machine$double.eps) {
    stop("Insufficient effective sample size for covariance estimation.")
  }
  list(
    cov = .symmetrize(scatter / denom),
    mean = mu,
    denom = denom
  )
}

# Empirical covariance with optional sample weights
.cov_empirical <- function(X, center = TRUE, sample_weight = NULL) {
  m <- .weighted_moments(X, sample_weight = sample_weight)
  cfit <- .cov_from_moments(
    sum_w = m$sum_w,
    sum_w2 = m$sum_w2,
    sum_x = m$sum_x,
    sum_xx = m$sum_xx,
    center = center
  )
  list(
    cov = cfit$cov,
    mean = cfit$mean,
    sum_w = m$sum_w,
    sum_w2 = m$sum_w2,
    sum_x = m$sum_x,
    sum_xx = m$sum_xx,
    denom = cfit$denom
  )
}

# Robust MCD covariance estimator
.cov_mcd <- function(X, center = TRUE, sample_weight = NULL) {
  if (!center) {
    stop("cov_estimator='mcd' requires center = TRUE.")
  }
  if (!is.null(sample_weight)) {
    stop("sample_weight is not supported for cov_estimator='mcd'.")
  }
  if (!requireNamespace("robustbase", quietly = TRUE)) {
    stop("cov_estimator='mcd' requires package 'robustbase'.")
  }
  fit <- robustbase::covMcd(X)
  list(
    cov = .symmetrize(fit$cov),
    mean = fit$center,
    sum_w = NA_real_,
    sum_w2 = NA_real_,
    sum_x = NULL,
    sum_xx = NULL,
    denom = NA_real_
  )
}

# Robust Tyler covariance estimator
.cov_tyler <- function(X,
                       center = TRUE,
                       sample_weight = NULL,
                       tol = 1e-6,
                       max_iter = 200,
                       eps = 1e-6) {
  if (!center) {
    stop("cov_estimator='tyler' requires center = TRUE.")
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tyler_tol must be a positive finite scalar.")
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 1) {
    stop("tyler_max_iter must be a positive integer.")
  }
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("tyler_eps must be a positive finite scalar.")
  }
  max_iter <- as.integer(max_iter)

  n <- nrow(X)
  d <- ncol(X)
  w <- .validate_sample_weight(sample_weight, n, allow_null = TRUE)
  if (is.null(w)) {
    w <- rep(1, n)
  }
  sw <- sum(w)
  sw2 <- sum(w ^ 2)
  mu <- colSums(X * w) / sw
  Xc <- sweep(X, 2, mu, "-")

  S <- crossprod(Xc) / max(n - 1L, 1L)
  S <- .symmetrize(S + eps * diag(d))
  trS <- sum(diag(S))
  if (!is.finite(trS) || trS <= eps) {
    S <- diag(1, d)
  } else {
    S <- S * (d / trS)
  }

  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    Sinv <- tryCatch(
      chol2inv(chol(S)),
      error = function(e) chol2inv(chol(S + eps * diag(d)))
    )
    q <- rowSums((Xc %*% Sinv) * Xc)
    q <- pmax(q, eps)
    Xw <- Xc * sqrt(w / q)
    S_new <- (d / sw) * crossprod(Xw)
    S_new <- .symmetrize(S_new + eps * diag(d))
    tr_new <- sum(diag(S_new))
    if (!is.finite(tr_new) || tr_new <= eps) {
      S_new <- diag(1, d)
    } else {
      S_new <- S_new * (d / tr_new)
    }
    rel <- norm(S_new - S, type = "F") / pmax(norm(S, type = "F"), eps)
    S <- S_new
    if (rel < tol) {
      converged <- TRUE
      break
    }
  }
  if (!converged) {
    warning("Tyler covariance estimator did not converge within max_iter.")
  }

  # Tyler's iteration estimates only the SHAPE (normalized to trace = d); it is
  # scale-free, so whitening with it would not produce unit-variance data.
  # Restore the scale with a robust, consistency-correcting factor: under an
  # elliptical model x^T S_shape^{-1} x / sigma^2 ~ chi^2_d, so
  # sigma^2 = median_i(x_i^T S_shape^{-1} x_i) / qchisq(0.5, d).
  Sinv <- tryCatch(
    chol2inv(chol(S)),
    error = function(e) chol2inv(chol(S + eps * diag(d)))
  )
  q_scale <- rowSums((Xc %*% Sinv) * Xc)
  q_scale <- q_scale[is.finite(q_scale) & q_scale > 0]
  if (length(q_scale) > 0L) {
    sigma2 <- stats::median(q_scale) / stats::qchisq(0.5, df = d)
    if (is.finite(sigma2) && sigma2 > 0) {
      S <- S * sigma2
    }
  }

  list(
    cov = S,
    mean = mu,
    sum_w = sw,
    sum_w2 = sw2,
    sum_x = colSums(X * w),
    sum_xx = crossprod(X, X * w),
    denom = sw - (sw2 / sw)
  )
}

# Unified covariance estimator dispatch
.estimate_covariance <- function(X,
                                 center = TRUE,
                                 sample_weight = NULL,
                                 cov_estimator = c("empirical", "mcd", "tyler"),
                                 tyler_tol = 1e-6,
                                 tyler_max_iter = 200,
                                 tyler_eps = 1e-6) {
  cov_estimator <- match.arg(cov_estimator)
  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    stop("center must be TRUE or FALSE.")
  }
  if (cov_estimator == "empirical") {
    return(.cov_empirical(X, center = center, sample_weight = sample_weight))
  }
  if (cov_estimator == "mcd") {
    return(.cov_mcd(X, center = center, sample_weight = sample_weight))
  }
  .cov_tyler(
    X,
    center = center,
    sample_weight = sample_weight,
    tol = tyler_tol,
    max_iter = tyler_max_iter,
    eps = tyler_eps
  )
}

# Default unsupervised whitening score
.default_whitening_score <- function(Z) {
  d <- check_whitening(Z)
  -(abs(d$diag_mean - 1) + d$offdiag_frob)
}

# Default supervised score: nearest-centroid accuracy
.default_supervised_score <- function(Z_train, y_train, Z_valid, y_valid) {
  y_train <- as.factor(y_train)
  y_valid <- as.factor(y_valid)
  cls <- levels(y_train)
  centers <- do.call(
    rbind,
    lapply(cls, function(cl) {
      colMeans(Z_train[y_train == cl, , drop = FALSE])
    })
  )
  rownames(centers) <- cls
  dmat <- as.matrix(stats::dist(rbind(centers, Z_valid)))
  nc <- length(cls)
  dv <- dmat[seq_len(nc), (nc + 1):ncol(dmat), drop = FALSE]
  pred <- cls[max.col(-t(dv), ties.method = "first")]
  mean(pred == as.character(y_valid))
}

# Build fold indices for cross-validation
.make_cv_folds <- function(n, k = 5L, y = NULL, seed = NULL) {
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 2) {
    stop("cv_folds must be an integer >= 2.")
  }
  k <- as.integer(k)
  if (k > n) stop("cv_folds cannot exceed nrow(X).")
  if (!is.null(seed)) set.seed(seed)

  if (is.null(y)) {
    idx <- sample.int(n)
    return(split(idx, rep(seq_len(k), length.out = n)))
  }

  if (length(y) != n) {
    stop("y must have length nrow(X).")
  }

  y <- as.factor(y)
  folds <- vector("list", k)
  for (i in seq_len(k)) folds[[i]] <- integer(0)

  for (cl in levels(y)) {
    idx <- which(y == cl)
    idx <- sample(idx)
    bin <- split(idx, rep(seq_len(k), length.out = length(idx)))
    for (i in seq_len(k)) {
      folds[[i]] <- c(folds[[i]], bin[[i]])
    }
  }
  lapply(folds, sort)
}

# Return component spectrum used for dimensionality decisions
.component_spectrum_info <- function(Sigma,
                                     method,
                                     k = NULL,
                                     eig_method = c("auto", "base", "rspectra"),
                                     fast = FALSE) {
  method <- .normalize_whiten_method(method)
  eig_method <- match.arg(eig_method)
  if (method == "PCA-cor") {
    R <- stats::cov2cor(Sigma)
    dec <- .eigen_spd(R, k = k, method = eig_method, fast = fast)
    total <- ncol(R)
  } else if (method %in% c("PCA", "SVD")) {
    dec <- .eigen_spd(Sigma, k = k, method = eig_method, fast = fast)
    total <- sum(diag(Sigma))
  } else {
    stop("Component spectrum is only defined for PCA, PCA-cor, and SVD.")
  }
  list(
    values = pmax(dec$values, 0),
    total = max(total, .Machine$double.eps),
    engine = dec$engine
  )
}

# Pick minimum k such that cumulative explained variance reaches threshold
.n_comp_from_threshold <- function(values, var_threshold) {
  if (length(values) == 0L) stop("values must not be empty.")
  denom <- sum(values)
  if (denom <= 0) {
    return(length(values))
  }
  csum <- cumsum(values) / denom
  k <- which(csum >= var_threshold)[1L]
  if (is.na(k)) length(values) else as.integer(k)
}

# Explained variance ratio for first k components
.explained_ratio <- function(values, k, total = NULL) {
  if (is.null(k)) return(NA_real_)
  if (is.null(total) && k >= length(values)) return(1)
  denom <- if (is.null(total)) sum(values) else total
  if (denom <= 0) return(NA_real_)
  k_use <- min(k, length(values))
  sum(values[seq_len(k_use)]) / denom
}

# Deterministic sign alignment for eigen/singular vectors
.fix_component_sign <- function(U, sign_ref = NULL) {
  if (!is.matrix(U)) stop("U must be a matrix.")
  k <- ncol(U)
  if (k == 0L) return(U)

  if (!is.null(sign_ref)) {
    if (!is.matrix(sign_ref)) stop("sign_ref must be a matrix.")
    if (nrow(sign_ref) != nrow(U)) {
      stop("sign_ref must have the same number of rows as U.")
    }
  }

  s <- rep(1, k)
  for (j in seq_len(k)) {
    if (!is.null(sign_ref) && j <= ncol(sign_ref)) {
      dotv <- sum(U[, j] * sign_ref[, j])
      if (is.finite(dotv) && abs(dotv) > .Machine$double.eps) {
        s[j] <- if (dotv >= 0) 1 else -1
        next
      }
    }

    idx <- which.max(abs(U[, j]))
    anchor <- U[idx, j]
    s[j] <- if (anchor >= 0) 1 else -1
  }
  sweep(U, 2, s, "*")
}

# Apply covariance shrinkage toward a target.
#   "identity": (1 - lambda) * S + lambda * (tr(S)/d) * I   (default)
#   "diagonal": (1 - lambda) * S + lambda * diag(diag(S))   (preserves variances,
#               shrinks correlations toward zero)
# The identity target preserves the eigenvectors of S (only the eigenvalues are
# shifted affinely), which the analytic shrinkage path in tuning relies on.
.shrink_cov <- function(S, lambda, target = c("identity", "diagonal")) {
  target <- match.arg(target)
  if (lambda <= 0) return(S)
  d <- ncol(S)
  if (target == "identity") {
    target_val <- mean(diag(S))
    S_t <- (1 - lambda) * S + lambda * diag(target_val, d)
  } else {
    S_t <- (1 - lambda) * S
    diag(S_t) <- diag(S)
  }
  .symmetrize(S_t)
}

# Shrunk eigenvalues for the identity target, derived analytically from the
# unshrunk eigenvalues. Adding lambda * c * I (c = mean eigenvalue = tr(S)/d)
# leaves the eigenvectors unchanged and maps each eigenvalue d_i to
# (1 - lambda) * d_i + lambda * c. Used by the fast tuning path.
.shrink_eigvals_identity <- function(d_vals, lambda) {
  if (lambda <= 0) return(d_vals)
  cc <- mean(d_vals)
  (1 - lambda) * d_vals + lambda * cc
}

# Build whitening outputs for covariance-based methods (PCA, SVD, ZCA) directly
# from an eigendecomposition S = U diag(d_vals) U^T of the (already shrunk)
# covariance. Returns W, its inverse map inv_tW, and the retained U / D.
# This mirrors the math in PCA()/SVD()/ZCA() and is shared by the fast tuning path.
.whiten_from_cov_eig <- function(method, U, d_vals, n_comp = NULL,
                                 sign_ref = NULL, col_names = NULL) {
  method <- .normalize_whiten_method(method)
  if (!method %in% c("PCA", "SVD", "ZCA")) {
    stop("'.whiten_from_cov_eig' supports only PCA, SVD, and ZCA.")
  }
  d <- nrow(U)

  if (method == "ZCA") {
    inv_sqrt <- 1 / sqrt(d_vals)
    W <- tcrossprod(sweep(U, 2L, inv_sqrt, "*"), U)
    inv_tW <- tcrossprod(sweep(U, 2L, sqrt(d_vals), "*"), U)
    row_names <- paste0("L", seq_len(d))
    U_out <- U
    D_out <- d_vals
  } else {
    k <- if (is.null(n_comp)) length(d_vals) else as.integer(n_comp)
    Uk <- U[, seq_len(k), drop = FALSE]
    dk <- d_vals[seq_len(k)]
    Uk <- .fix_component_sign(Uk, sign_ref = sign_ref)
    W <- sweep(t(Uk), 1L, 1 / sqrt(dk), "*")
    inv_tW <- sweep(t(Uk), 1L, sqrt(dk), "*")
    row_names <- paste0(if (method == "SVD") "SV" else "PC", seq_len(k))
    U_out <- Uk
    D_out <- dk
  }

  W <- .set_matrix_attr(W, row_names = row_names, col_names = col_names, method = method)
  list(W = W, inv_tW = inv_tW, U = U_out, D = D_out)
}

# Subtract a per-column vector from every row without sweep() overhead.
# rep(mu, each = n) lays mu out in column-major order to match X.
.center_cols <- function(X, mu) {
  n <- nrow(X)
  X - rep(mu, each = n)
}

# Drop rows containing any non-finite value; report how many were removed.
.drop_nonfinite_rows <- function(X, na_action = c("error", "omit"),
                                 what = "X") {
  na_action <- match.arg(na_action)
  bad <- !is.finite(X)
  if (!any(bad)) return(X)
  if (na_action == "error") {
    stop(sprintf("%s must contain only finite values.", what))
  }
  keep <- rowSums(bad) == 0
  n_drop <- sum(!keep)
  if (n_drop > 0L) {
    warning(sprintf("Dropped %d row(s) of %s containing non-finite values.", n_drop, what))
  }
  X[keep, , drop = FALSE]
}

# Indices of the upper triangle (including diagonal) in column-major order.
.upper_tri_index <- function(p) {
  which(upper.tri(matrix(0, p, p), diag = TRUE))
}

# Vectorize a symmetric tangent matrix L into a vector of length p(p+1)/2.
# Off-diagonal entries are scaled by sqrt(2) so that the Euclidean norm of the
# vector equals the Frobenius norm of L (isometry of the tangent inner product).
.spd_tangent_vec <- function(L) {
  p <- ncol(L)
  w <- L
  off <- sqrt(2)
  w[upper.tri(w)] <- w[upper.tri(w)] * off
  w[.upper_tri_index(p)]
}

# Inverse of .spd_tangent_vec: rebuild the symmetric matrix from its vector.
.spd_tangent_unvec <- function(v, p) {
  L <- matrix(0, p, p)
  idx <- .upper_tri_index(p)
  L[idx] <- v
  off <- sqrt(2)
  ut <- upper.tri(L)
  L[ut] <- L[ut] / off
  L <- L + t(L)
  diag(L) <- diag(L) / 2
  L
}

# Relative tolerance for the scale-invariant positive-definiteness test: a
# symmetric matrix is treated as PD when its smallest eigenvalue exceeds
# .PD_REL_TOL times its largest (i.e. condition number below ~1 / .PD_REL_TOL).
.PD_REL_TOL <- 1e-14

# Basic check: Sigma must be symmetric; optionally positive-definite
.check_symmetric_pd <- function(Sigma, tol = 1e-8, require_pd = TRUE, return_eigen = TRUE) {
  if (!is.matrix(Sigma)) stop("Sigma must be a matrix.")
  if (nrow(Sigma) != ncol(Sigma)) stop("Sigma must be a square matrix.")
  if (!is.numeric(Sigma)) stop("Sigma must be numeric.")
  if (any(!is.finite(Sigma))) stop("Sigma must contain only finite values.")
  if (!isTRUE(all.equal(Sigma, t(Sigma), tolerance = tol))) {
    stop("Sigma must be symmetric.")
  }

  if (!return_eigen) {
    if (require_pd) {
      chol_ok <- tryCatch(
        {
          chol(Sigma)
          TRUE
        },
        error = function(e) FALSE
      )
      if (!chol_ok) {
        eig <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
        min_eig <- min(eig)
        stop(sprintf("Sigma must be positive-definite. Min eigenvalue = %g. Try using lambda > 0 in whiten_model().", min_eig))
      }
    }
    return(invisible(NULL))
  }

  eig <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  if (require_pd) {
    # Scale-invariant positive-definiteness test: compare the smallest
    # eigenvalue to the largest (a numerical-rank tolerance) rather than to an
    # absolute floor, so a well-conditioned matrix at any magnitude is accepted.
    max_eig <- max(eig)
    if (!is.finite(max_eig) || max_eig <= 0 || min(eig) <= max_eig * .PD_REL_TOL) {
      stop(sprintf("Sigma must be positive-definite (min/max eigenvalue = %g / %g). Try using lambda > 0 in whiten_model().", min(eig), max_eig))
    }
  }
  invisible(eig)
}

# Moore-Penrose pseudo-inverse using SVD
.pinv <- function(A, tol = NULL) {
  if (!is.matrix(A)) stop("A must be a matrix.")
  if (!is.numeric(A)) stop("A must be numeric.")

  s <- svd(A)
  if (is.null(tol)) {
    tol <- sqrt(.Machine$double.eps) * max(dim(A)) * max(s$d)
  }
  nz <- s$d > tol

  if (!any(nz)) {
    return(matrix(0, nrow = ncol(A), ncol = nrow(A)))
  }

  s$v[, nz, drop = FALSE] %*%
    diag(1 / s$d[nz], nrow = sum(nz), ncol = sum(nz)) %*%
    t(s$u[, nz, drop = FALSE])
}

# Ledoit-Wolf style shrinkage intensity toward scaled identity
.lambda_lw <- function(Xc) {
  if (!is.matrix(Xc)) stop("Xc must be a matrix.")
  n <- nrow(Xc)
  p <- ncol(Xc)
  if (n <= 1L || p <= 1L) return(0)

  S <- crossprod(Xc) / n
  mu <- sum(diag(S)) / p
  target <- diag(mu, p)
  gamma_hat <- sum((S - target) ^ 2)
  if (gamma_hat <= .Machine$double.eps) return(1)

  x2_sum <- rowSums(Xc ^ 2)
  term1 <- sum(x2_sum ^ 2)
  XS <- Xc %*% S
  term2 <- 2 * sum(XS * Xc)
  term3 <- n * sum(S ^ 2)

  beta_hat <- (term1 - term2 + term3) / (n ^ 2)
  lambda <- beta_hat / gamma_hat
  min(max(lambda, 0), 1)
}

# Schaefer-Strimmer (2005) analytic shrinkage intensity toward the diagonal
# target (shrink correlations toward zero, preserve variances). Computed on the
# correlation scale, so it is scale-invariant -- the right companion to
# shrink_target = "diagonal". Xc is centered (not standardized) data.
.lambda_ss <- function(Xc) {
  n <- nrow(Xc)
  p <- ncol(Xc)
  if (n <= 2L || p <= 1L) return(0)
  sds <- sqrt(colSums(Xc ^ 2) / (n - 1L))
  if (any(sds <= .Machine$double.eps)) return(1)
  W <- sweep(Xc, 2L, sds, "/")              # standardized variables
  R <- crossprod(W) / (n - 1L)              # sample correlation matrix
  wbar <- R * ((n - 1L) / n)
  sum_m2 <- crossprod(W ^ 2)                # sum_k (w_ik w_jk)^2
  var_r <- (n / (n - 1L) ^ 3) * (sum_m2 - n * wbar ^ 2)
  off <- row(R) != col(R)
  num <- sum(var_r[off])
  den <- sum((R ^ 2)[off])
  if (!is.finite(den) || den <= .Machine$double.eps) return(1)
  min(max(num / den, 0), 1)
}

# OAS shrinkage intensity toward scaled identity
.lambda_oas <- function(Xc) {
  if (!is.matrix(Xc)) stop("Xc must be a matrix.")
  n <- nrow(Xc)
  p <- ncol(Xc)
  if (n <= 1L || p <= 1L) return(0)

  S <- crossprod(Xc) / n
  mu <- sum(diag(S)) / p
  alpha <- mean(S ^ 2)

  num <- alpha + mu ^ 2
  den <- (n + 1) * (alpha - (mu ^ 2 / p))
  if (!is.finite(den) || den <= .Machine$double.eps) return(1)

  lambda <- num / den
  min(max(lambda, 0), 1)
}

# Ensure matrix is numerically symmetric
.symmetrize <- function(A) {
  (A + t(A)) / 2
}

# Resolve (model, data) arguments, accepting the deprecated (data, model) order.
# The model is identified by its S3 class; a deprecation warning is emitted when
# the old order is detected.
.order_model_first <- function(a, b, model_class, fn) {
  if (inherits(a, model_class)) {
    return(list(model = a, data = b))
  }
  if (inherits(b, model_class)) {
    warning(sprintf("%s(data, model) is deprecated; use %s(model, data).", fn, fn),
            call. = FALSE)
    return(list(model = b, data = a))
  }
  stop(sprintf("%s() requires a '%s' object.", fn, model_class))
}

# Matrix sqrt for SPD matrix
.sqrtm_spd <- function(S, eps = 1e-12) {
  eig <- eigen(.symmetrize(S), symmetric = TRUE)
  vals <- pmax(eig$values, eps)
  eig$vectors %*% diag(sqrt(vals), nrow = length(vals), ncol = length(vals)) %*% t(eig$vectors)
}

# Matrix inverse sqrt for SPD matrix
.invsqrtm_spd <- function(S, eps = 1e-12) {
  eig <- eigen(.symmetrize(S), symmetric = TRUE)
  vals <- pmax(eig$values, eps)
  eig$vectors %*% diag(1 / sqrt(vals), nrow = length(vals), ncol = length(vals)) %*% t(eig$vectors)
}

# Matrix log for SPD matrix
.logm_spd <- function(S, eps = 1e-12) {
  eig <- eigen(.symmetrize(S), symmetric = TRUE)
  vals <- pmax(eig$values, eps)
  eig$vectors %*% diag(log(vals), nrow = length(vals), ncol = length(vals)) %*% t(eig$vectors)
}

# Matrix exponential for symmetric matrix
.expm_sym <- function(S) {
  eig <- eigen(.symmetrize(S), symmetric = TRUE)
  eig$vectors %*% diag(exp(eig$values), nrow = length(eig$values), ncol = length(eig$values)) %*% t(eig$vectors)
}

# Affine-invariant Riemannian mean of SPD matrices
.riemann_mean_spd <- function(C_list, tol = 1e-6, max_iter = 50, eps = 1e-12) {
  n <- length(C_list)
  if (n == 0L) stop("C_list must not be empty.")

  G <- Reduce("+", C_list) / n
  G <- .symmetrize(G)

  for (iter in seq_len(max_iter)) {
    # Share eigendecomposition for sqrtm and invsqrtm
    eig_G <- eigen(.symmetrize(G), symmetric = TRUE)
    vals_G <- pmax(eig_G$values, eps)
    G_half <- eig_G$vectors %*% diag(sqrt(vals_G), nrow = length(vals_G)) %*% t(eig_G$vectors)
    G_inv_half <- eig_G$vectors %*% diag(1 / sqrt(vals_G), nrow = length(vals_G)) %*% t(eig_G$vectors)

    delta <- matrix(0, nrow(G), ncol(G))
    for (C in C_list) {
      tmp <- G_inv_half %*% C %*% G_inv_half
      delta <- delta + .logm_spd(tmp, eps = eps)
    }
    delta <- delta / n
    crit <- norm(delta, type = "F")
    G <- G_half %*% .expm_sym(delta) %*% G_half
    G <- .symmetrize(G)

    if (crit < tol) break
  }
  G
}

#' Check Condition Number of Data Matrix
#'
#' Computes the condition number of the sample covariance matrix of
#' \code{X}, defined as the ratio of the largest to smallest positive
#' eigenvalue. High values indicate ill-conditioning.
#'
#' @param X Numeric matrix with at least 2 rows and 2 columns.
#'
#' @return The condition number (a single positive numeric value, or
#'   \code{Inf} if no positive eigenvalues exist).
#'
#' @seealso \code{\link{whiten_model}}, \code{\link{check_whitening}}
#'
#' @examples
#' X <- matrix(rnorm(100 * 8), 100, 8)
#' check_condition(X)
#'
#' @export
check_condition <- function(X) {
  if (!is.matrix(X)) stop("X must be a matrix.")
  if (!is.numeric(X)) stop("X must be a numeric matrix.")
  if (any(!is.finite(X))) stop("X must contain only finite values.")
  if (nrow(X) < 2L) stop("X must have at least 2 rows.")
  if (ncol(X) < 2L) stop("X must have at least 2 columns.")
  ev <- eigen(cov(X), symmetric = TRUE, only.values = TRUE)$values
  ev <- ev[ev > 0]
  if (length(ev) == 0L) return(Inf)
  max(ev) / min(ev)
}
