#' Plug-in variance for CC (complete-case) estimator of (beta, sigma)
#'
#' @description
#' Computes the plug-in variance–covariance matrix for the complete-case estimator
#' \code{\link{get_CC}}. Supports four modes matching that estimator:
#' \itemize{
#'   \item \strong{CASE A}: estimate \eqn{(\beta,\sigma)} with \code{z_data}
#'   \item \strong{CASE B}: estimate \eqn{(\beta,\sigma)} without \code{z_data}
#'   \item \strong{CASE C}: estimate \eqn{\beta} with \code{z_data} and fixed \code{sigma}
#'   \item \strong{CASE D}: estimate \eqn{\beta} without \code{z_data} and fixed \code{sigma}
#' }
#'
#' @param beta Numeric vector of estimated \eqn{\beta}.
#' @param y_data Numeric vector.
#' @param w_data Numeric vector of \eqn{W=\min(X,C)}.
#' @param delta_data Logical or \{0,1\} vector: \eqn{\Delta=1\{X \le C\}}.
#' @param z_data Optional numeric vector (default \code{NULL}); if provided, must match length of \code{y_data}.
#' @param sigma Numeric scalar \eqn{\hat\sigma}. Ignored when \code{sigma_fixed} is provided.
#' @param sigma_fixed Optional positive scalar. If supplied, \eqn{\sigma} is treated as known
#'   and the returned variance is for \eqn{\beta} only.
#'
#' @returns
#' A variance–covariance matrix.
#' \itemize{
#'   \item CASE A/B (\code{sigma_fixed = NULL}): matrix of size \code{length(beta)+1}.
#'   \item CASE C/D (\code{sigma_fixed} provided): matrix of size \code{length(beta)}.
#' }
#'
#' @importFrom MASS ginv
#' @importFrom stats var
#' @export
variance_CC <- function(beta,
                        y_data,
                        w_data,
                        delta_data,
                        z_data = NULL,
                        sigma,
                        sigma_fixed = NULL) {

  ## ---- Basic checks --------------------------------------------------------
  stopifnot(is.numeric(beta), is.numeric(y_data), is.numeric(w_data))
  n <- length(y_data)
  if (!(length(w_data) == n && length(delta_data) == n)) {
    stop("`y_data`, `w_data`, and `delta_data` must have equal length.", call. = FALSE)
  }
  if (!is.logical(delta_data)) {
    if (!all(delta_data %in% c(0, 1))) stop("`delta_data` must be logical or {0,1}.", call. = FALSE)
    delta_data <- as.logical(delta_data)
  }
  if (!is.null(z_data) && (!is.numeric(z_data) || length(z_data) != n)) {
    stop("If provided, `z_data` must be numeric and same length as `y_data`.", call. = FALSE)
  }
  if (is.null(sigma_fixed)) {
    if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
      stop("`sigma` must be a positive numeric scalar.", call. = FALSE)
    }
  } else {
    if (!is.numeric(sigma_fixed) || length(sigma_fixed) != 1L || !is.finite(sigma_fixed) || sigma_fixed <= 0) {
      stop("`sigma_fixed` must be a positive numeric scalar.", call. = FALSE)
    }
  }

  ## ---- Full-data scores evaluated at X = W (complete cases only) ----------
  # With z: beta = (β0, β1, β2, β3); Y = β0 + β1 X + β2 Z + β3 ZX + ε, ε ~ N(0, σ^2)
  S_beta_f <- function(beta, y, x, z, sigma) {
    res <- y - (beta[1] + beta[2]*x + beta[3]*z + beta[4]*z*x)
    c(res * c(1, x, z, z*x) / sigma^2,
      res^2 / (2*sigma^4) - 1/(2*sigma^2))
  }
  S_beta_f <- Vectorize(S_beta_f, vectorize.args = "x")

  # No z: beta = (β0, β1); Y = β0 + β1 X + ε, ε ~ N(0, σ^2)
  S_beta_f_nz <- function(beta, y, x, sigma) {
    res <- y - (beta[1] + beta[2]*x)
    c(res * c(1, x) / sigma^2,
      res^2 / (2*sigma^4) - 1/(2*sigma^2))
  }
  S_beta_f_nz <- Vectorize(S_beta_f_nz, vectorize.args = "x")

  ## ---- Case dispatch -------------------------------------------------------
  z_is_null      <- is.null(z_data)
  sigma_is_fixed <- !is.null(sigma_fixed)
  idx_cc <- which(delta_data)
  n_cc   <- length(idx_cc)
  if (n_cc < 2L) stop("Not enough complete cases (need at least 2) to estimate variance.", call. = FALSE)

  if (!sigma_is_fixed && !z_is_null) {
    ## ===== CASE A: (beta, sigma) WITH z ====================================
    len.beta <- length(beta)
    Smat <- matrix(NA_real_, nrow = n_cc, ncol = len.beta + 1L)
    k <- 1L
    for (i in idx_cc) {
      Smat[k, ] <- S_beta_f(beta, y = y_data[i], x = w_data[i], z = z_data[i], sigma = sigma)
      k <- k + 1L
    }
    B <- stats::var(Smat)
    return(MASS::ginv(B) / n_cc)

  } else if (!sigma_is_fixed && z_is_null) {
    ## ===== CASE B: (beta, sigma) WITHOUT z =================================
    len.beta <- length(beta) # expected 2
    Smat <- matrix(NA_real_, nrow = n_cc, ncol = len.beta + 1L)
    k <- 1L
    for (i in idx_cc) {
      Smat[k, ] <- S_beta_f_nz(beta, y = y_data[i], x = w_data[i], sigma = sigma)
      k <- k + 1L
    }
    B <- stats::var(Smat)
    return(MASS::ginv(B) / n_cc)

  } else if (sigma_is_fixed && !z_is_null) {
    ## ===== CASE C: beta ONLY WITH z, sigma fixed ===========================
    len.beta  <- length(beta)
    sigma_use <- as.numeric(sigma_fixed)
    # build full (beta,sigma) score then take beta–beta block
    Smat <- matrix(NA_real_, nrow = n_cc, ncol = len.beta + 1L)
    k <- 1L
    for (i in idx_cc) {
      Smat[k, ] <- S_beta_f(beta, y = y_data[i], x = w_data[i], z = z_data[i], sigma = sigma_use)
      k <- k + 1L
    }
    B_full <- stats::var(Smat)
    V_full <- MASS::ginv(B_full) / n_cc
    return(V_full[seq_len(len.beta), seq_len(len.beta), drop = FALSE])

  } else {
    ## ===== CASE D: beta ONLY WITHOUT z, sigma fixed ========================
    len.beta  <- length(beta) # expected 2
    sigma_use <- as.numeric(sigma_fixed)
    Smat <- matrix(NA_real_, nrow = n_cc, ncol = len.beta + 1L)
    k <- 1L
    for (i in idx_cc) {
      Smat[k, ] <- S_beta_f_nz(beta, y = y_data[i], x = w_data[i], sigma = sigma_use)
      k <- k + 1L
    }
    B_full <- stats::var(Smat)
    V_full <- MASS::ginv(B_full) / n_cc
    return(V_full[seq_len(len.beta), seq_len(len.beta), drop = FALSE])
  }
}
