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
#'   Underscore aliases \code{"ZCA_cor"} and \code{"PCA_cor"} are also accepted.
#' @param n_comp Integer (optional); number of components to keep for
#'   \code{"PCA"}, \code{"PCA-cor"}, and \code{"SVD"}.
#' @param var_threshold Numeric in (0, 1]; cumulative explained-variance
#'   threshold for automatic \code{n_comp} selection.
#' @param lambda Numeric in [0, 1] or \code{"auto"}; covariance shrinkage
#'   strength.
#' @param lambda_method Method used when \code{lambda = "auto"}.
#'   One of \code{"oas"} or \code{"lw"}.
#' @param sample_weight Optional non-negative sample weights with
#'   length \code{nrow(X)}.
#' @param cov_estimator Covariance estimator; one of \code{"empirical"},
#'   \code{"mcd"}, or \code{"tyler"}.
#' @param tyler_tol Convergence tolerance for \code{cov_estimator = "tyler"}.
#' @param tyler_max_iter Maximum iterations for Tyler covariance estimation.
#' @param tyler_eps Numerical stabilization floor for Tyler estimation.
#' @param eig_method Eigen solver backend; one of \code{"auto"},
#'   \code{"base"}, or \code{"rspectra"}.
#' @param fast Logical; if \code{TRUE}, allow faster approximate settings
#'   for iterative eigensolvers.
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
#' @seealso \code{\link{whiten_model}} for the fit/transform workflow,
#'   \code{\link{whiten_batch}} for batch whitening.
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
                          sample_weight = NULL,
                          cov_estimator = c("empirical", "mcd", "tyler"),
                          tyler_tol = 1e-6,
                          tyler_max_iter = 200,
                          tyler_eps = 1e-6,
                          eig_method = c("auto", "base", "rspectra"),
                          fast = FALSE,
                          sign_reference = NULL) {
  if (length(method) > 1L) method <- method[[1L]]
  method <- .normalize_whiten_method(method)
  lambda_method <- match.arg(lambda_method)
  cov_estimator <- match.arg(cov_estimator)
  eig_method <- match.arg(eig_method)
  
  model <- whiten_model(
    X,
    center = center,
    method = method,
    n_comp = n_comp,
    var_threshold = var_threshold,
    lambda = lambda,
    lambda_method = lambda_method,
    sample_weight = sample_weight,
    cov_estimator = cov_estimator,
    tyler_tol = tyler_tol,
    tyler_max_iter = tyler_max_iter,
    tyler_eps = tyler_eps,
    eig_method = eig_method,
    fast = fast,
    sign_reference = sign_reference
  )
  Z     <- stats::predict(model, X)
  
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
#'   Underscore aliases \code{"ZCA_cor"} and \code{"PCA_cor"} are also accepted.
#' @param n_comp Integer (optional); number of components to keep for
#'   \code{"PCA"}, \code{"PCA-cor"}, and \code{"SVD"}.
#' @param var_threshold Numeric in (0, 1]; cumulative explained-variance
#'   threshold for automatic \code{n_comp} selection.
#' @param lambda Numeric in [0, 1] or \code{"auto"}; covariance shrinkage
#'   strength.
#' @param lambda_method Method used when \code{lambda = "auto"}.
#'   One of \code{"oas"} or \code{"lw"}.
#' @param sample_weight Optional sample weights.
#'   For a single matrix, provide one numeric vector.
#'   For multiple matrices in \code{mode = "independent"}, provide a list
#'   with one vector per matrix.
#' @param cov_estimator Covariance estimator; one of \code{"empirical"},
#'   \code{"mcd"}, or \code{"tyler"}.
#' @param tyler_tol Convergence tolerance for \code{cov_estimator = "tyler"}.
#' @param tyler_max_iter Maximum iterations for Tyler covariance estimation.
#' @param tyler_eps Numerical stabilization floor for Tyler estimation.
#' @param eig_method Eigen solver backend; one of \code{"auto"},
#'   \code{"base"}, or \code{"rspectra"}.
#' @param fast Logical; if \code{TRUE}, allow faster approximate settings
#'   for iterative eigensolvers.
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
#' @param parallel Logical; if \code{TRUE}, use multicore processing for
#'   per-matrix operations on Unix-like systems.
#' @param n_cores Integer; number of worker processes when
#'   \code{parallel = TRUE}.
#'
#' @return A list of results as returned by \code{whiten_matrix()}.
#'
#' @seealso \code{\link{whiten_matrix}}, \code{\link{whiten_model}},
#'   \code{\link{euclidean_alignment}} for direct EA access.
#'
#' @examples
#' set.seed(1)
#' X_list <- list(
#'   matrix(rnorm(200 * 8), 200, 8),
#'   matrix(rnorm(180 * 8), 180, 8)
#' )
#' out <- whiten_batch(X_list, mode = "shared_model",
#'                     method = "PCA", n_comp = 4, lambda = 0.1)
#' length(out)
#'
#' @export
whiten_batch <- function(X_list, center = TRUE,
                         method = c("SVD", "ZCA", "ZCA-cor",
                                    "PCA", "PCA-cor", "Cholesky"),
                         n_comp = NULL,
                         var_threshold = NULL,
                         lambda = "auto",
                         lambda_method = c("oas", "lw"),
                         sample_weight = NULL,
                         cov_estimator = c("empirical", "mcd", "tyler"),
                         tyler_tol = 1e-6,
                         tyler_max_iter = 200,
                         tyler_eps = 1e-6,
                         eig_method = c("auto", "base", "rspectra"),
                         fast = FALSE,
                         sign_reference = NULL,
                         shared_model = FALSE,
                         mode = c("independent", "shared_model", "ea"),
                         ea_mean = c("logeuclid", "riemann", "euclid"),
                         ea_input = c("auto", "raw", "cov"),
                         ea_tol = 1e-6,
                         ea_max_iter = 50,
                         parallel = FALSE,
                         n_cores = 1L) {
  if (length(method) > 1L) method <- method[[1L]]
  method <- .normalize_whiten_method(method)
  lambda_method <- match.arg(lambda_method)
  cov_estimator <- match.arg(cov_estimator)
  eig_method <- match.arg(eig_method)
  mode <- match.arg(mode)
  ea_mean <- match.arg(ea_mean)
  ea_input <- match.arg(ea_input)
  if (!is.list(X_list)) stop("X_list must be a list of matrices.")
  if (length(X_list) == 0L) stop("X_list must not be empty.")
  if (!is.logical(shared_model) || length(shared_model) != 1L || is.na(shared_model)) {
    stop("shared_model must be TRUE or FALSE.")
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("parallel must be TRUE or FALSE.")
  }
  if (!is.numeric(n_cores) || length(n_cores) != 1L || !is.finite(n_cores) || n_cores < 1) {
    stop("n_cores must be a positive integer.")
  }
  n_cores <- as.integer(n_cores)
  if (parallel && .Platform$OS.type != "unix") {
    warning("parallel = TRUE is only supported on Unix-like systems; falling back to sequential processing.")
    parallel <- FALSE
  }

  apply_list <- function(data, fun) {
    if (!parallel || n_cores <= 1L || length(data) <= 1L) {
      return(lapply(data, fun))
    }
    parallel::mclapply(data, fun, mc.cores = n_cores)
  }

  resolve_weights <- function(w_arg, mats, mode) {
    n_list <- length(mats)
    out <- vector("list", n_list)
    if (is.null(w_arg)) return(out)

    if (is.list(w_arg)) {
      if (length(w_arg) != n_list) {
        stop("When sample_weight is a list, it must have the same length as X_list.")
      }
      for (i in seq_len(n_list)) {
        out[[i]] <- .validate_sample_weight(w_arg[[i]], nrow(mats[[i]]), allow_null = TRUE)
      }
      return(out)
    }

    if (!is.numeric(w_arg) || any(!is.finite(w_arg))) {
      stop("sample_weight must be NULL, a numeric vector, or a list of numeric vectors.")
    }

    if (n_list == 1L || identical(mode, "shared_model")) {
      out[[1L]] <- .validate_sample_weight(w_arg, nrow(mats[[1L]]), allow_null = FALSE)
      return(out)
    }

    stop("For multiple matrices in independent mode, sample_weight must be a list with one weight vector per matrix.")
  }

  if (shared_model && identical(mode, "independent")) {
    warning("'shared_model' parameter is deprecated; use mode = 'shared_model' instead.", call. = FALSE)
    mode <- "shared_model"
  }
  if (mode == "independent" && length(X_list) > 1L) {
    warning("Independent per-matrix whitening may remove between-session covariance structure. Consider mode = 'ea' or mode = 'shared_model'.")
  }
  weight_list <- resolve_weights(sample_weight, X_list, mode)
  
  if (mode == "ea") {
    if (any(vapply(weight_list, Negate(is.null), logical(1)))) {
      warning("sample_weight is ignored when mode = 'ea'.")
    }
    ea <- euclidean_alignment(
      X_list,
      input = ea_input,
      mean_method = ea_mean,
      center = center,
      tol = ea_tol,
      max_iter = ea_max_iter
    )
    
    out <- apply_list(ea$aligned, function(Z) {
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
      sample_weight = weight_list[[1L]],
      cov_estimator = cov_estimator,
      tyler_tol = tyler_tol,
      tyler_max_iter = tyler_max_iter,
      tyler_eps = tyler_eps,
      eig_method = eig_method,
      fast = fast,
      sign_reference = sign_reference
    )
    
    out <- apply_list(X_list, function(X) {
      Z <- stats::predict(model, X)
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
  
  apply_list(
    seq_along(X_list),
    function(i) {
      whiten_matrix(
        X_list[[i]],
        center = center,
        method = method,
        n_comp = n_comp,
        var_threshold = var_threshold,
        lambda = lambda,
        lambda_method = lambda_method,
        sample_weight = weight_list[[i]],
        cov_estimator = cov_estimator,
        tyler_tol = tyler_tol,
        tyler_max_iter = tyler_max_iter,
        tyler_eps = tyler_eps,
        eig_method = eig_method,
        fast = fast,
        sign_reference = sign_reference
      )
    }
  )
}
