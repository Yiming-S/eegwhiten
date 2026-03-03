############################################################
# Whitening via SVD of the covariance matrix
############################################################

#' SVD-based whitening on the covariance matrix
#'
#' Perform whitening using the singular value decomposition (SVD)
#' of a symmetric positive-definite covariance matrix. For SPD
#' matrices, SVD and eigen-decomposition coincide up to ordering,
#' but this function may be convenient in some workflows.
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
#' @param fast Logical; if \code{TRUE}, allow faster approximate settings
#'   for iterative eigensolvers.
#'
#' @return A list with some of the elements:
#'   \item{W}{Whitening matrix based on the SVD.}
#'   \item{Phi}{Factor loadings in the original space.}
#'   \item{Psi}{Standardized loadings.}
#'
#' @export
SVD <- function(Sigma,
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
  if (any(v <= 0)) stop("Diagonal elements of Sigma must be positive.")
  
  # --- Dimensionality Reduction ---
  d_full <- ncol(Sigma)
  if (!is.null(n_comp)) {
    if (n_comp < 1 || n_comp > d_full) stop("n_comp out of range.")
    k <- n_comp
  } else {
    k <- d_full
  }
  dec <- .eigen_spd(Sigma, k = k, method = eig_method, fast = fast)
  U <- dec$vectors
  D <- dec$values
  
  if (any(D <= 0)) stop("Sigma must be positive-definite.")
  
  U <- .fix_component_sign(U, sign_ref = sign_ref)
  
  result <- list()
  
  if (returnW) {
    # W: [k x d]
    # W = D^-1/2 * U^T
    W <- sweep(t(U), 1L, 1 / sqrt(D), "*")
    
    result$W <- .set_matrix_attr(
      W,
      row_names = paste0("SV", seq_len(k)),
      col_names = colnames(Sigma),
      method    = "SVD"
    )
  }
  
  if (PhiPsi) {
    # Phi: U * D^1/2
    Phi <- sweep(U, 2L, sqrt(D), "*")
    
    # Psi: V^-1/2 * Phi
    Psi <- sweep(Phi, 1L, 1 / sqrt(v), "*")
    
    row_names <- colnames(Sigma)
    col_names <- paste0("SV", seq_len(k))
    
    result$Phi <- .set_matrix_attr(Phi, row_names, col_names, "SVD")
    result$Psi <- .set_matrix_attr(Psi, row_names, col_names, "SVD")
  }

  if (return_decomp) {
    result$decomp <- list(
      U = U,
      D = D,
      inv_tW = sweep(t(U), 1L, sqrt(D), "*")
    )
  }
  
  result
}
