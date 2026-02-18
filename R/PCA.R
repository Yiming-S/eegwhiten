
############################################################
# FILENAME: PCA.R
# DESCRIPTION: PCA-based whitening
# UPDATES: Added n_comp support
############################################################

#' PCA-based whitening on the covariance matrix
#'
#' Perform whitening based on the eigen-decomposition.
#' Supports dimensionality reduction via \code{n_comp}.
#'
#' @param Sigma Symmetric positive-definite covariance matrix.
#' @param n_comp Integer; number of components to keep.
#' @param sign_ref Optional reference vectors used to stabilize component
#'   signs across runs. Must have compatible dimensions with eigenvectors.
#' @param returnW Logical; if TRUE, return the whitening matrix \code{W}.
#' @param PhiPsi Logical; if TRUE, return factor loadings.
#' @param return_decomp Logical; if TRUE, return decomposition terms used for
#'   fast inverse transformation.
#'
#' @export
PCA <- function(Sigma,
                n_comp = NULL,
                sign_ref = NULL,
                returnW = TRUE,
                PhiPsi = TRUE,
                return_decomp = FALSE) {
  .check_symmetric_pd(Sigma, require_pd = TRUE)
  
  v <- diag(Sigma)
  if (any(v <= 0)) {
    stop("Diagonal elements of Sigma must be positive.")
  }
  
  eSigma <- eigen(Sigma, symmetric = TRUE)
  U      <- eSigma$vectors
  lambda <- eSigma$values
  
  # --- Dimensionality Reduction Logic ---
  d <- length(lambda)
  if (!is.null(n_comp)) {
    if (n_comp < 1 || n_comp > d) stop("n_comp out of range.")
    k <- n_comp
  } else {
    k <- d
  }
  
  # Truncate
  U      <- U[, 1:k, drop = FALSE]
  lambda <- lambda[1:k]
  
  # Check positive eigenvalues after truncation (usually handled by check_pd but safe to check)
  if (any(lambda <= 0)) stop("Non-positive eigenvalues detected in PCA.")
  
  # Fix sign ambiguity with deterministic / reference-driven alignment
  U <- .fix_component_sign(U, sign_ref = sign_ref)
  
  result <- list()
  
  if (returnW) {
    # W: [k x d]. 
    # Transforms X (n x d) -> Z (n x k) via Z = X %*% t(W)
    # W = D^-1/2 * U^T
    W <- tcrossprod(diag(1 / sqrt(lambda), nrow = k, ncol = k), U)
    
    result$W <- .set_matrix_attr(
      W,
      row_names = paste0("PC", seq_len(k)),
      col_names = colnames(Sigma),
      method    = "PCA"
    )
  }
  
  if (PhiPsi) {
    # Phi: Loadings in original space [d x k]
    Phi <- U %*% diag(sqrt(lambda), nrow = k, ncol = k)
    
    # Psi: Standardized loadings [d x k]
    # Psi_ij = Phi_ij / sigma_i
    Psi <- diag(1 / sqrt(v)) %*% Phi
    
    row_names <- colnames(Sigma)
    col_names <- paste0("PC", seq_len(k))
    
    result$Phi <- .set_matrix_attr(Phi, row_names, col_names, "PCA")
    result$Psi <- .set_matrix_attr(Psi, row_names, col_names, "PCA")
  }

  if (return_decomp) {
    result$decomp <- list(
      U = U,
      D = lambda,
      inv_tW = diag(sqrt(lambda), nrow = k, ncol = k) %*% t(U)
    )
  }
  
  result
}
