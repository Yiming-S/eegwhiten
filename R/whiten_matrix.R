############################################################
# Convenience wrappers for single-matrix and batch whitening
############################################################

#' One-shot whitening for a single EEG matrix
#'
#' Convenience wrapper around \code{whiten_model()} that fits
#' a whitening model on \code{X} and immediately applies it to
#' \code{X}, returning both the whitened data and the underlying
#' model parameters.
#'
#' @param X Numeric matrix [n_trials x n_channels].
#' @param center Logical; whether to subtract the column means.
#' @param method Whitening method; one of \code{"SVD"}, \code{"ZCA"},
#'   \code{"ZCA-cor"}, \code{"PCA"}, \code{"PCA-cor"}, \code{"Cholesky"}.
#' @param n_comp Integer (optional); number of components to keep for
#'   \code{"PCA"}, \code{"PCA-cor"}, and \code{"SVD"}.
#' @param var_threshold Numeric in (0, 1]; cumulative explained-variance
#'   threshold for automatic \code{n_comp} selection.
#' @param lambda Numeric in [0, 1] or \code{"auto"}; covariance shrinkage
#'   strength.
#' @param lambda_method Method used when \code{lambda = "auto"}.
#'   One of \code{"oas"} or \code{"lw"}.
#' @param sign_reference Optional reference vectors used for stable component
#'   sign orientation in \code{"PCA"}, \code{"PCA-cor"}, and \code{"SVD"}.
#'
#' @return A list with components:
#'   \item{Z}{Whitened data matrix.}
#'   \item{W}{Whitening matrix.}
#'   \item{center}{Centering vector used.}
#'   \item{method}{Method name.}
#'   \item{n_comp}{Number of retained components (if applicable).}
#'   \item{var_threshold}{Variance threshold used for automatic selection.}
#'   \item{explained_var}{Cumulative explained variance of retained components.}
#'   \item{lambda}{Shrinkage value used during fitting.}
#'   \item{lambda_method}{Automatic shrinkage estimator, if used.}
#'   \item{model}{The fitted \code{whiten_model} object.}
#'
#' @examples
#' set.seed(42)
#' X <- matrix(rnorm(200 * 16), 200, 16)
#' colnames(X) <- paste0("Ch", seq_len(16))
#' res <- whiten_matrix(X, center = TRUE, method = "PCA")
#' str(res$Z)
#'
#' @export
whiten_matrix <- function(X, center = TRUE,
                          method = c("SVD", "ZCA", "ZCA-cor",
                                     "PCA", "PCA-cor", "Cholesky"),
                          n_comp = NULL,
                          var_threshold = NULL,
                          lambda = "auto",
                          lambda_method = c("oas", "lw"),
                          sign_reference = NULL) {
  method <- match.arg(method)
  lambda_method <- match.arg(lambda_method)
  
  model <- whiten_model(
    X,
    center = center,
    method = method,
    n_comp = n_comp,
    var_threshold = var_threshold,
    lambda = lambda,
    lambda_method = lambda_method,
    sign_reference = sign_reference
  )
  Z     <- predict(model, X)
  
  list(
    Z      = Z,
    W      = model$W,
    center = model$center,
    method = model$method,
    n_comp = model$n_comp,
    var_threshold = model$var_threshold,
    explained_var = model$explained_var,
    lambda = model$lambda,
    lambda_method = model$lambda_method,
    model  = model
  )
}

#' Batch whitening for a list of EEG matrices
#'
#' Each matrix in \code{X_list} is whitened independently using
#' \code{whiten_matrix}. For proper train/test usage, fit a model
#' on training data with \code{whiten_model()} and call
#' \code{predict()} on test data instead of using this function.
#'
#' @param X_list List of numeric matrices, each [n_trials x n_channels].
#' @param center Logical; whether to subtract the column means in each matrix.
#' @param method Whitening method; one of \code{"SVD"}, \code{"ZCA"},
#'   \code{"ZCA-cor"}, \code{"PCA"}, \code{"PCA-cor"}, \code{"Cholesky"}.
#' @param n_comp Integer (optional); number of components to keep for
#'   \code{"PCA"}, \code{"PCA-cor"}, and \code{"SVD"}.
#' @param var_threshold Numeric in (0, 1]; cumulative explained-variance
#'   threshold for automatic \code{n_comp} selection.
#' @param lambda Numeric in [0, 1] or \code{"auto"}; covariance shrinkage
#'   strength.
#' @param lambda_method Method used when \code{lambda = "auto"}.
#'   One of \code{"oas"} or \code{"lw"}.
#' @param sign_reference Optional reference vectors used for stable component
#'   sign orientation in \code{"PCA"}, \code{"PCA-cor"}, and \code{"SVD"}.
#' @param shared_model Logical; if \code{TRUE}, fit one model on the first
#'   matrix and apply it to all matrices in \code{X_list}. If \code{FALSE},
#'   each matrix is whitened independently.
#' @param mode Batch whitening strategy:
#'   \code{"independent"} (fit each matrix separately),
#'   \code{"shared_model"} (fit once on first matrix), or
#'   \code{"ea"} (Euclidean Alignment with global reference covariance).
#' @param ea_mean Mean type used in \code{"ea"} mode; one of
#'   \code{"riemann"}, \code{"logeuclid"}, or \code{"euclid"}.
#' @param ea_input Input interpretation in \code{"ea"} mode; one of
#'   \code{"auto"}, \code{"raw"}, or \code{"cov"}.
#' @param ea_tol Tolerance for Riemannian mean convergence in \code{"ea"} mode.
#' @param ea_max_iter Maximum iterations for Riemannian mean in \code{"ea"} mode.
#'
#' @return A list of results as returned by \code{whiten_matrix()}.
#'
#' @export
whiten_batch <- function(X_list, center = TRUE,
                         method = c("SVD", "ZCA", "ZCA-cor",
                                    "PCA", "PCA-cor", "Cholesky"),
                         n_comp = NULL,
                         var_threshold = NULL,
                         lambda = "auto",
                         lambda_method = c("oas", "lw"),
                         sign_reference = NULL,
                         shared_model = FALSE,
                         mode = c("independent", "shared_model", "ea"),
                         ea_mean = c("riemann", "logeuclid", "euclid"),
                         ea_input = c("auto", "raw", "cov"),
                         ea_tol = 1e-6,
                         ea_max_iter = 50) {
  method <- match.arg(method)
  lambda_method <- match.arg(lambda_method)
  mode <- match.arg(mode)
  ea_mean <- match.arg(ea_mean)
  ea_input <- match.arg(ea_input)
  if (!is.list(X_list)) stop("X_list must be a list of matrices.")
  if (length(X_list) == 0L) stop("X_list must not be empty.")
  if (!is.logical(shared_model) || length(shared_model) != 1L || is.na(shared_model)) {
    stop("shared_model must be TRUE or FALSE.")
  }
  if (shared_model && identical(mode, "independent")) {
    mode <- "shared_model"
  }
  if (mode == "independent" && length(X_list) > 1L) {
    warning("Independent per-matrix whitening may remove between-session covariance structure. Consider mode = 'ea' or mode = 'shared_model'.")
  }
  
  if (mode == "ea") {
    ea <- euclidean_alignment(
      X_list,
      input = ea_input,
      mean_method = ea_mean,
      center = center,
      tol = ea_tol,
      max_iter = ea_max_iter
    )
    
    out <- lapply(ea$aligned, function(Z) {
      list(
        Z = Z,
        W = ea$W,
        center = if (ea$input_type == "raw") rep(0, ncol(Z)) else NULL,
        method = "EA",
        n_comp = NULL,
        var_threshold = NULL,
        explained_var = NA_real_,
        lambda = NA_real_,
        lambda_method = NA_character_,
        model = ea$model
      )
    })
    
    attr(out, "shared_model") <- TRUE
    attr(out, "model") <- ea$model
    attr(out, "alignment") <- "ea"
    return(out)
  }
  
  if (mode == "shared_model") {
    model <- whiten_model(
      X_list[[1]],
      center = center,
      method = method,
      n_comp = n_comp,
      var_threshold = var_threshold,
      lambda = lambda,
      lambda_method = lambda_method,
      sign_reference = sign_reference
    )
    
    out <- lapply(X_list, function(X) {
      Z <- predict(model, X)
      list(
        Z      = Z,
        W      = model$W,
        center = model$center,
        method = model$method,
        n_comp = model$n_comp,
        var_threshold = model$var_threshold,
        explained_var = model$explained_var,
        lambda = model$lambda,
        lambda_method = model$lambda_method,
        model  = model
      )
    })
    
    attr(out, "shared_model") <- TRUE
    attr(out, "model") <- model
    return(out)
  }
  
  lapply(
    X_list,
    function(X) {
      whiten_matrix(
        X,
        center = center,
        method = method,
        n_comp = n_comp,
        var_threshold = var_threshold,
        lambda = lambda,
        lambda_method = lambda_method,
        sign_reference = sign_reference
      )
    }
  )
}
