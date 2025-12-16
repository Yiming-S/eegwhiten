############################################################
# FILENAME: whiten_core.R
# DESCRIPTION: Core whitening API for EEG data
# UPDATES: Added Regularization (lambda) and Dim Reduction (n_comp)
############################################################

#' Fit a whitening matrix from a covariance matrix
#'
#' Compute a whitening matrix \code{W} from a symmetric
#' positive-definite covariance matrix using one of the
#' implemented methods.
#'
#' @param Sigma Symmetric positive-definite covariance matrix.
#' @param method Whitening method.
#' @param n_comp Integer (optional); Number of components to keep (for PCA/SVD).
#'
#' @return Whitening matrix \code{W}.
#'
#' @export
whiten_fit <- function(Sigma,
                       method = c("SVD", "ZCA", "ZCA-cor",
                                  "PCA", "PCA-cor", "Cholesky"),
                       n_comp = NULL) {
  method <- match.arg(method)
  
  # All implemented methods expect symmetric PD covariance
  .check_symmetric_pd(Sigma, require_pd = TRUE)
  
  # Prepare arguments
  args <- list(Sigma = Sigma, returnW = TRUE, PhiPsi = FALSE)
  
  # Only pass n_comp to methods that support it to avoid "unused argument" errors
  if (method %in% c("PCA", "PCA-cor", "SVD") && !is.null(n_comp)) {
    args$n_comp <- n_comp
  }
  
  # Call the specific function (e.g., PCA(Sigma, ...))
  res <- do.call(method, args)
  
  if (is.null(res$W)) stop("Whitening matrix W is NULL.")
  res$W
}

#' Fit a whitening model for EEG data
#'
#' The input matrix \code{X} is assumed to have trials/epochs in rows
#' and channels or features in columns. This function estimates
#' a centering vector and a whitening matrix, optionally using
#' regularization (shrinkage) and dimensionality reduction.
#'
#' @param X Numeric matrix [n_trials x n_channels].
#' @param center Logical; whether to subtract column means before
#'   computing the covariance.
#' @param method Whitening method; one of \code{"SVD"}, \code{"ZCA"},
#'   \code{"ZCA-cor"}, \code{"PCA"}, \code{"PCA-cor"}, \code{"Cholesky"}.
#' @param n_comp Integer (optional); Number of components to keep.
#'   Only applicable for "PCA", "PCA-cor", and "SVD".
#' @param lambda Numeric (0 to 1); Regularization parameter for covariance shrinkage.
#'   0 = No regularization (Empirical covariance).
#'   1 = Full shrinkage to target (scaled Identity).
#'   Values between 0 and 1 mix the two. Useful for high-dimensional/low-sample EEG data.
#'
#' @return An object of class \code{"whiten_model"}.
#'
#' @export
whiten_model <- function(X, center = TRUE,
                         method = c("SVD", "ZCA", "ZCA-cor",
                                    "PCA", "PCA-cor", "Cholesky"),
                         n_comp = NULL,
                         lambda = 0) {
  method <- match.arg(method)
  
  if (!is.matrix(X)) stop("X must be a matrix.")
  n <- nrow(X)
  d <- ncol(X)
  
  if (d <= 1L) {
    stop("X must have more than one column.")
  }
  
  # Validation for n_comp
  if (!is.null(n_comp)) {
    if (n_comp < 1 || n_comp > d) stop("n_comp must be between 1 and ncol(X).")
    if (method %in% c("ZCA", "ZCA-cor", "Cholesky") && n_comp != d) {
      warning("Dimensionality reduction (n_comp) is ignored for ZCA/Cholesky methods. Use PCA if you need dimension reduction.")
      n_comp <- NULL 
    }
  }
  
  # Ensure column names exist
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("V", seq_len(d))
  }
  
  # 1. Centering
  mu <- if (center) colMeans(X) else rep(0, d)
  names(mu) <- colnames(X)
  
  # Centered data
  Xc <- sweep(X, 2, mu, "-")
  
  # 2. Covariance Estimation with Regularization
  # Empirical covariance of centered data
  S <- crossprod(Xc) / (n - 1L)
  
  # Apply Shrinkage if lambda > 0
  if (lambda > 0) {
    if (lambda > 1) stop("lambda must be <= 1.")
    # Target: average eigenvalue * Identity
    target_val <- mean(diag(S))
    target_mat <- diag(target_val, d)
    
    # Shrinkage formula: (1 - lambda) * S + lambda * Target
    S <- (1 - lambda) * S + lambda * target_mat
  }
  
  # 3. Fit Whitening Matrix
  # We use the internal whiten_fit helper or call directly
  # Using whiten_fit logic here for clarity
  
  # Prepare arguments
  fit_args <- list(Sigma = S, returnW = TRUE, PhiPsi = FALSE)
  if (method %in% c("PCA", "PCA-cor", "SVD") && !is.null(n_comp)) {
    fit_args$n_comp <- n_comp
  }
  
  res <- do.call(method, fit_args)
  W   <- res$W
  
  structure(
    list(
      W      = W,
      center = mu,
      method = method,
      n_comp = n_comp,
      lambda = lambda,
      dim_in = d,
      dim_out = ncol(W) # Note: W is usually (k x d) in definition, but here we store as needed
    ),
    class = "whiten_model"
  )
}

#' Apply a whitening model to new EEG data
#'
#' @param object A \code{whiten_model} object.
#' @param newdata Numeric matrix with the same number of columns
#'   as the data used to fit the model.
#' @param ... Unused.
#'
#' @return Whitened data matrix (possibly with reduced dimensions).
#'
#' @export
#' @method predict whiten_model
predict.whiten_model <- function(object, newdata, ...) {
  if (!is.matrix(newdata)) stop("newdata must be a matrix.")
  
  d <- length(object$center)
  if (ncol(newdata) != d) {
    stop(sprintf("Mismatch: Model expects %d columns, but newdata has %d.", d, ncol(newdata)))
  }
  
  # Center using stored mean
  Xc <- sweep(newdata, 2, object$center, "-")
  
  # Z = centered X times whitening matrix
  # W is typically [k x d]. tcrossprod(Xc, W) -> [n x d] * [d x k] = [n x k]
  Z <- tcrossprod(Xc, object$W)
  
  # Update column names based on output dimension
  k <- ncol(Z)
  colnames(Z) <- paste0("PC", seq_len(k)) # Generic name, or "L"
  attr(Z, "method") <- object$method
  
  Z
}

#' Inverse transform from whitened space back to original EEG space
#'
#' @param Z Whitened data matrix.
#' @param model A \code{whiten_model} object.
#'
#' @return Approximate reconstruction of the original data matrix.
#'
#' @export
unwhiten <- function(Z, model) {
  if (!inherits(model, "whiten_model")) {
    stop("model must be a 'whiten_model' object.")
  }
  if (!is.matrix(Z)) stop("Z must be a matrix.")
  
  # Z = Xc %*% t(W)
  # If W is square (d x d) and invertible: Xc = Z %*% solve(t(W))
  # If W is rectangular (k x d) (PCA reduction): Xc = Z %*% W (pseudo-inverse concept for PCA)
  # Note: For PCA, W = D^{-1/2} U^T.  t(W) = U D^{-1/2}.
  # Inverse direction: Z -> Xc. 
  # Xc_hat = Z %*% (t(W))^-1 is hard for rectangular.
  # Better approach for PCA: Xc ~ Z %*% sqrt(Lambda) %*% U^T 
  # But we only stored W. 
  # For PCA whitening, the pseudo-inverse of W is involved.
  # This simple 'solve' only works for square matrices (ZCA, full PCA).
  
  if (model$dim_in != model$dim_out) {
    # Rectangular case (Dimensionality Reduction)
    # We use the Moore-Penrose generalized inverse of t(W)
    # t(W) is [d x k]. solve(t(W)) won't work.
    # We use ginv from MASS or simple lstsq. 
    # Here we assume user has MASS or we implement a simple solver.
    # For now, we error out or warn.
    # A simple reconstruction for PCA whiten is:
    # W = L^-1/2 U^T.  W^+ = U L^1/2.
    # Z = X U L^-1/2. -> X_rec = Z L^1/2 U^T.
    # We can compute pinv(t(W)).
    
    # Let's try to find a solution Xc such that Xc %*% t(W) = Z
    # This is Xc = Z %*% pinv(t(W))
    
    # Minimal dependency implementation of pinv using SVD
    tW <- t(model$W)
    s <- svd(tW)
    # Filter small singular values
    tol <- sqrt(.Machine$double.eps) * max(dim(tW)) * max(s$d)
    nz <- s$d > tol
    if (!any(nz)) {
      inv_tW <- matrix(0, nrow=ncol(tW), ncol=nrow(tW))
    } else {
      inv_tW <- s$u[, nz] %*% diag(1/s$d[nz], sum(nz)) %*% t(s$v[, nz])
    }
    
    Xc <- Z %*% t(inv_tW) # because pinv(A) approx A^-1
    
  } else {
    # Square case
    Xc <- Z %*% solve(t(model$W))
  }
  
  X  <- sweep(Xc, 2, model$center, "+")
  
  colnames(X) <- names(model$center)
  X
}

#' Diagnostic check for whitening quality
#'
#' @param Z Whitened data matrix.
#'
#' @return List of diagnostics.
#' @export
check_whitening <- function(Z) {
  if (!is.matrix(Z)) stop("Z must be a matrix.")
  
  S <- cov(Z)
  d <- ncol(Z)
  
  diag_vals  <- diag(S)
  diag_range <- range(diag_vals)
  
  off_diag <- S
  diag(off_diag) <- 0
  off_norm <- norm(off_diag, type = "F")
  
  list(
    diag_range   = diag_range,
    diag_mean    = mean(diag_vals),
    offdiag_frob = off_norm,
    dim          = d,
    cov_matrix   = S
  )
}