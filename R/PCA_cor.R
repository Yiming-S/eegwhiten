############################################################
# FILENAME: PCA_cor.R
# DESCRIPTION: PCA on correlation matrix
# UPDATES: Added n_comp support and stability checks
############################################################

#' PCA-based whitening on the correlation matrix
#'
#' Perform whitening using PCA on the correlation matrix implied
#' by \code{Sigma}. This is useful when variables (EEG channels or
#' features) have very different scales and one wants to whiten in
#' a scale-invariant way.
#'
#' @param Sigma Symmetric positive-definite covariance matrix.
#' @param n_comp Integer; number of components to keep. If NULL, keeps all.
#' @param sign_ref Optional reference vectors used to stabilize component
#'   signs across runs.
#' @param returnW Logical; if TRUE, return the whitening matrix \code{W}.
#' @param PhiPsi Logical; if TRUE, return factor loadings \code{Phi}
#'   and standardized loadings \code{Psi}.
#' @param return_decomp Logical; if TRUE, return decomposition terms used for
#'   fast inverse transformation.
#' @param eig_method Eigen solver backend; one of \code{"auto"},
#'   \code{"base"}, or \code{"rspectra"}.
#' @param fast Logical; if \code{TRUE}, use relaxed checks and faster iterative
#'   eigensolver settings when available.
#'
#' @return A list with some of the elements:
#'   \item{W}{Whitening matrix based on the correlation structure.}
#'   \item{Phi}{Factor loadings in the original space.}
#'   \item{Psi}{Standardized loadings.}
#'
#' @export
PCA_cor <- function(Sigma,
                    n_comp = NULL,
                    sign_ref = NULL,
                    returnW = TRUE,
                    PhiPsi = TRUE,
                    return_decomp = FALSE,
                    eig_method = c("auto", "base", "rspectra"),
                    fast = FALSE) {
  if (!isTRUE(attr(Sigma, ".checked_spd"))) {
    .check_symmetric_pd(Sigma, require_pd = TRUE)
  }
  eig_method <- match.arg(eig_method)
  if (!is.logical(fast) || length(fast) != 1L || is.na(fast)) {
    stop("fast must be TRUE or FALSE.")
  }
  
  v <- diag(Sigma)
  if (any(v <= 0)) {
    stop("Diagonal elements of Sigma must be positive.")
  }
  
  R <- stats::cov2cor(Sigma)
  
  # Diagonal of R should be numerically 1
  tol_diag <- if (fast) 1e-5 else sqrt(.Machine$double.eps)
  if (any(abs(diag(R) - 1) > tol_diag)) {
    stop("Diagonal elements of the correlation matrix must be approximately 1.")
  }
  
  # --- Dimensionality Reduction and Stability Check ---
  d <- ncol(R)
  if (!is.null(n_comp)) {
    if (n_comp < 1 || n_comp > d) stop("n_comp out of range.")
    k <- n_comp
  } else {
    k <- d
  }
  dec <- .eigen_spd(R, k = k, method = eig_method, fast = fast)
  G <- dec$vectors
  theta <- dec$values
  
  # Crucial stability check for PCA whitening
  # If eigenvalues of R are <= 0, 1/sqrt(theta) will fail or produce NaNs
  if (any(theta <= .Machine$double.eps)) {
    stop("PCA_cor: Non-positive eigenvalues detected in correlation matrix. Try reducing n_comp or increasing regularization lambda in whiten_model().")
  }
  
  # Fix sign ambiguity with deterministic / reference-driven alignment
  G <- .fix_component_sign(G, sign_ref = sign_ref)
  
  result <- list()
  
  if (returnW) {
    # W for PCA-cor with reduction
    # W = Theta^-1/2 * G^T * V^-1/2
    # Dimensions: (k x k) * (k x d) * (d x d) = (k x d)
    # Using diag(..., nrow=k, ncol=k) ensures correct behavior even if k=1
    
    W <- sweep(t(G), 1L, 1 / sqrt(theta), "*")
    W <- sweep(W, 2L, 1 / sqrt(v), "*")
    
    result$W <- .set_matrix_attr(
      W,
      row_names = paste0("PC", seq_len(k)),
      col_names = colnames(Sigma),
      method    = "PCA-cor"
    )
  }
  
  if (PhiPsi) {
    # Psi (Standardized Loadings): G * Theta^1/2
    Psi <- sweep(G, 2L, sqrt(theta), "*")
    
    # Phi (Loadings): V^1/2 * Psi
    Phi <- sweep(Psi, 1L, sqrt(v), "*")
    
    row_names <- colnames(Sigma)
    col_names <- paste0("PC", seq_len(k))
    
    result$Phi <- .set_matrix_attr(Phi, row_names, col_names, "PCA-cor")
    result$Psi <- .set_matrix_attr(Psi, row_names, col_names, "PCA-cor")
  }

  if (return_decomp) {
    inv_tW <- sweep(t(G), 1L, sqrt(theta), "*")
    inv_tW <- sweep(inv_tW, 2L, sqrt(v), "*")
    result$decomp <- list(
      U = G,
      D = theta,
      scale_diag = v,
      inv_tW = inv_tW
    )
  }
  
  result
}
