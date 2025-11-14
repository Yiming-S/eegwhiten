############################################################
# PCA-based whitening on the correlation matrix
############################################################

#' PCA-based whitening on the correlation matrix
#'
#' Perform whitening using PCA on the correlation matrix implied
#' by \code{Sigma}. This is useful when variables (EEG channels or
#' features) have very different scales and one wants to whiten in
#' a scale-invariant way.
#'
#' @param Sigma Symmetric positive-definite covariance matrix.
#' @param returnW Logical; if TRUE, return the whitening matrix \code{W}.
#' @param PhiPsi Logical; if TRUE, return factor loadings \code{Phi}
#'   and standardized loadings \code{Psi}.
#'
#' @return A list with some of the elements:
#'   \item{W}{Whitening matrix based on the correlation structure.}
#'   \item{Phi}{Factor loadings in the original space.}
#'   \item{Psi}{Standardized loadings.}
#'
#' @export
PCA_cor <- function(Sigma, returnW = TRUE, PhiPsi = TRUE) {
  .check_symmetric_pd(Sigma, require_pd = TRUE)
  
  v <- diag(Sigma)
  if (any(v <= 0)) {
    stop("Diagonal elements of Sigma must be positive.")
  }
  
  R <- cov2cor(Sigma)
  
  # Diagonal of R should be numerically 1
  if (any(abs(diag(R) - 1) > sqrt(.Machine$double.eps))) {
    stop("Diagonal elements of the correlation matrix must be approximately 1.")
  }
  
  eR    <- eigen(R, symmetric = TRUE)
  G     <- eR$vectors
  theta <- eR$values
  
  # Fix sign ambiguity
  G <- sweep(G, 2, sign(diag(G)), "*")
  
  result <- list()
  
  if (returnW) {
    # W: components x original features
    W <- diag(1 / sqrt(theta)) %*% t(G) %*% diag(1 / sqrt(v))
    result$W <- .set_matrix_attr(
      W,
      row_names = paste0("L", seq_len(ncol(Sigma))),
      col_names = colnames(Sigma),
      method    = "PCA-cor"
    )
  }
  
  if (PhiPsi) {
    Psi <- G %*% diag(sqrt(theta))
    Phi <- diag(sqrt(v)) %*% Psi
    
    row_names <- colnames(Sigma)
    col_names <- paste0("L", seq_len(ncol(Sigma)))
    
    result$Phi <- .set_matrix_attr(Phi, row_names, col_names, "PCA-cor")
    result$Psi <- .set_matrix_attr(Psi, row_names, col_names, "PCA-cor")
  }
  
  result
}