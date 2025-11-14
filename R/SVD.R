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
#' @param returnW Logical; if TRUE, return the whitening matrix \code{W}.
#' @param PhiPsi Logical; if TRUE, return factor loadings \code{Phi}
#'   and standardized loadings \code{Psi}.
#'
#' @return A list with some of the elements:
#'   \item{W}{Whitening matrix based on the SVD.}
#'   \item{Phi}{Factor loadings in the original space.}
#'   \item{Psi}{Standardized loadings.}
#'
#' @export
SVD <- function(Sigma, returnW = TRUE, PhiPsi = TRUE) {
  .check_symmetric_pd(Sigma, require_pd = TRUE)
  
  v <- diag(Sigma)
  if (any(v <= 0)) {
    stop("Diagonal elements of Sigma must be positive.")
  }
  
  svdSigma <- svd(Sigma)
  U <- svdSigma$u
  D <- svdSigma$d
  
  if (any(D <= 0)) {
    stop("Sigma must be positive-definite for SVD-based whitening.")
  }
  
  result <- list()
  
  if (returnW) {
    # For symmetric SPD Sigma, U D U^T = Sigma
    # Whitening: W = U D^{-1/2} U^T
    W <- U %*% diag(1 / sqrt(D)) %*% t(U)
    result$W <- .set_matrix_attr(
      W,
      row_names = paste0("L", seq_len(ncol(Sigma))),
      col_names = colnames(Sigma),
      method    = "SVD"
    )
  }
  
  if (PhiPsi) {
    Phi <- U %*% diag(sqrt(D))
    Psi <- diag(1 / sqrt(v)) %*% Phi
    
    row_names <- colnames(Sigma)
    col_names <- paste0("L", seq_len(ncol(Sigma)))
    
    result$Phi <- .set_matrix_attr(Phi, row_names, col_names, "SVD")
    result$Psi <- .set_matrix_attr(Psi, row_names, col_names, "SVD")
  }
  
  result
}