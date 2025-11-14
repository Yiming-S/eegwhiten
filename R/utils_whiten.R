
############################################################
# Internal utilities for whitening
############################################################

# Set dimnames and method attribute on a matrix
.set_matrix_attr <- function(x, row_names, col_names, method) {
  if (!is.null(row_names)) rownames(x) <- row_names
  if (!is.null(col_names)) colnames(x) <- col_names
  attr(x, "method") <- method
  x
}

# Basic check: Sigma must be symmetric; optionally positive-definite
.check_symmetric_pd <- function(Sigma, tol = 1e-8, require_pd = TRUE) {
  if (!is.matrix(Sigma)) stop("Sigma must be a matrix.")
  if (!isTRUE(all.equal(Sigma, t(Sigma)))) {
    stop("Sigma must be symmetric.")
  }
  eig <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  if (require_pd && any(eig <= tol)) {
    stop("Sigma must be positive-definite (eigenvalues <= tol).")
  }
  invisible(eig)
}