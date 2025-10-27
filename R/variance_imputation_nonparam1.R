#' Plug-in variance for imputation estimator of (beta, sigma) with nonparametrically estimated \eqn{f_{X\mid Z}}
#'
#' @description
#' Computes the plug-in variance–covariance matrix for the imputation estimator
#' \code{\link{get_imputation_nonparam1}}. Supports four modes matching that estimator:
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
#' @param surv_c_vals Numeric vector, length \code{length(y_data)}: conditional survival \eqn{S_{C\mid Y,Z}(W_i\mid Y_i,Z_i)} (or \eqn{S_{C\mid Y}(W_i\mid Y_i)} when \code{z_data} is \code{NULL}), evaluated at each observation.
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
#' @examples
#' library(truncnorm)
#' set.seed(1)
#' n = 500; beta = c(0, 3); sigma = 4
#' # Data generation
#' # Y|X ~ N(beta_0 + beta_1 X, 1)
#' x_data = rtruncnorm(n, a = -1, b = 1, mean = 0, sd = 1)
#' y_data = cbind(1,x_data) %*% beta + sigma * rnorm(n)
#' c_data = rtruncnorm(n, a = -1, b = 1, mean = 0, sd = 1)
#' z_data = rbinom(n, size = 1, prob = 0.6)
#' # Generate W and Delta from X and C
#' w_data = pmin(c_data,x_data)
#' delta_data = as.numeric(x_data<=c_data)
#' # Compute survival function
#' kern_gaussian = function(cond, cond_data, h = sd(cond_data)){ #Gaussian kernel
#'   exp(- 0.5 * apply(as.matrix(cond_data), 1, function(x) sum((cond - x)^2)) / h^2) / (sqrt(2*pi)*h)
#' }
#' surv_c_vals = sapply(1:n, function(i){
#'   conditional_Kaplan_Meier(w_data[i], w_data, 1-delta_data,
#'                            z_c = y_data[i], z_c_data = y_data,
#'                            z_d = z_data[i], z_d_data = z_data,
#'                            kern = kern_gaussian)})
#' surv_c_vals_nz = sapply(1:n, function(i){
#'   conditional_Kaplan_Meier(w_data[i], w_data, 1-delta_data,
#'                            z_c = y_data[i], z_c_data = y_data,
#'                            kern = kern_gaussian)})
#' # CASE A: true = (0, 3, 0, 0, 4)
#' variance_imputation_nonparam1(beta = c(0.5, 2.5, 0.1, 0.1),
#'                               y_data,
#'                               w_data,
#'                               delta_data,
#'                               z_data,
#'                               surv_c_vals = surv_c_vals,
#'                               sigma = 3.5)
#' # CASE B: true = (0, 3, 4)
#' variance_imputation_nonparam1(beta = c(0.5, 2.5),
#'                               y_data,
#'                               w_data,
#'                               delta_data,
#'                               surv_c_vals = surv_c_vals_nz,
#'                               sigma = 3.5)
#' # CASE C: true = (0, 3, 0, 0)
#' variance_imputation_nonparam1(beta = c(0.5, 2.5, 0.1, 0.1),
#'                               y_data,
#'                               w_data,
#'                               delta_data,
#'                               z_data,
#'                               surv_c_vals = surv_c_vals,
#'                               sigma_fixed = 4)
#' # CASE D: true = (0, 3)
#' variance_imputation_nonparam1(beta = c(0.5, 2.5),
#'                               y_data,
#'                               w_data,
#'                               delta_data,
#'                               surv_c_vals = surv_c_vals_nz,
#'                               sigma_fixed = 4)
#' @importFrom MASS ginv
#' @importFrom stats dnorm var
#' @export
variance_imputation_nonparam1 <- function(beta,
                                          y_data,
                                          w_data,
                                          delta_data,
                                          z_data = NULL,
                                          surv_c_vals,
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
  if (!is.numeric(surv_c_vals) || length(surv_c_vals) != n) {
    stop("`surv_c_vals` must be a numeric vector of length equal to `y_data`.", call. = FALSE)
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

  ## ---- Full-data scores (WITH / WITHOUT z) --------------------------------
  S_beta_f <- function(beta, y, x, z, sigma) {
    res <- y - (beta[1] + beta[2]*x + beta[3]*z + beta[4]*z*x)
    c(res * c(1, x, z, z*x) / sigma^2,
      res^2 / (2*sigma^4) - 1/(2*sigma^2))
  }
  S_beta_f <- Vectorize(S_beta_f, vectorize.args = "x")

  S_beta_f_nz <- function(beta, y, x, sigma) {
    res <- y - (beta[1] + beta[2]*x)
    c(res * c(1, x) / sigma^2,
      res^2 / (2*sigma^4) - 1/(2*sigma^2))
  }
  S_beta_f_nz <- Vectorize(S_beta_f_nz, vectorize.args = "x")

  ## ---- Kernel-imputation scores (WITH / WITHOUT z) ------------------------
  S_beta_kern1_imp <- function(beta, y, w, delta, z,
                               y_data, w_data, delta_data, z_data,
                               surv_c_vals,
                               sigma) {
    if (!delta) {
      idx_list <- which((delta_data == 1) & (w_data > w) & (z_data == z))
      if (length(idx_list) != 0L) {
        dens  <- stats::dnorm(y,
                              mean = cbind(1, w_data[idx_list], z, z*w_data[idx_list]) %*% beta,
                              sd   = sigma)
        wt    <- dens / surv_c_vals[idx_list]
        denom <- sum(wt)
        if (denom == 0) return(S_beta_f(beta, y, w, z, sigma))
        w_imp <- as.numeric(sum(w_data[idx_list] * wt) / denom)
        return(S_beta_f(beta, y, w_imp, z, sigma))
      } else {
        return(S_beta_f(beta, y, w, z, sigma))
      }
    } else {
      return(S_beta_f(beta, y, w, z, sigma))
    }
  }

  S_beta_kern1_imp_nz <- function(beta, y, w, delta,
                                  y_data, w_data, delta_data,
                                  surv_c_vals,
                                  sigma) {
    if (!delta) {
      idx_list <- which((delta_data == 1) & (w_data > w))
      if (length(idx_list) != 0L) {
        dens  <- stats::dnorm(y,
                              mean = cbind(1, w_data[idx_list]) %*% beta,
                              sd   = sigma)
        wt    <- dens / surv_c_vals[idx_list]
        denom <- sum(wt)
        if (denom == 0) return(S_beta_f_nz(beta, y, w, sigma))
        w_imp <- as.numeric(sum(w_data[idx_list] * wt) / denom)
        return(S_beta_f_nz(beta, y, w_imp, sigma))
      } else {
        return(S_beta_f_nz(beta, y, w, sigma))
      }
    } else {
      return(S_beta_f_nz(beta, y, w, sigma))
    }
  }

  ## ---- Case dispatch (kernel imputation) ----------------------------------
  z_is_null      <- is.null(z_data)
  sigma_is_fixed <- !is.null(sigma_fixed)

  if (!sigma_is_fixed && !z_is_null) {
    ## ===== CASE A: (beta, sigma) WITH z — kernel-imputation score =====
    len.beta <- length(beta)

    S_imp <- matrix(NA_real_, nrow = n, ncol = len.beta + 1L)
    S_ref <- matrix(NA_real_, nrow = n, ncol = len.beta + 1L)
    for (i in seq_len(n)) {
      S_imp[i, ] <- S_beta_kern1_imp(beta, y_data[i], w_data[i], delta_data[i], z_data[i],
                                     y_data, w_data, delta_data, z_data,
                                     surv_c_vals, sigma)
      S_ref[i, ] <- S_beta_f(        beta, y_data[i], w_data[i],                 z_data[i], sigma)
    }
    Ahat <- -stats::var(cbind(S_imp, S_ref))[seq_len(len.beta+1), -(seq_len(len.beta+1))]
    Bhat <-  stats::var(S_imp)
    return(MASS::ginv(Ahat) %*% Bhat %*% t(MASS::ginv(Ahat)) / n)

  } else if (!sigma_is_fixed && z_is_null) {
    ## ===== CASE B: (beta, sigma) WITHOUT z — kernel-imputation score =====
    len.beta <- length(beta)  # expected 2

    S_imp <- matrix(NA_real_, nrow = n, ncol = len.beta + 1L)
    S_ref <- matrix(NA_real_, nrow = n, ncol = len.beta + 1L)
    for (i in seq_len(n)) {
      S_imp[i, ] <- S_beta_kern1_imp_nz(beta, y_data[i], w_data[i], delta_data[i],
                                        y_data, w_data, delta_data,
                                        surv_c_vals, sigma)
      S_ref[i, ] <- S_beta_f_nz(          beta, y_data[i], w_data[i],                        sigma)
    }
    Ahat <- -stats::var(cbind(S_imp, S_ref))[seq_len(len.beta+1), -(seq_len(len.beta+1))]
    Bhat <-  stats::var(S_imp)
    return(MASS::ginv(Ahat) %*% Bhat %*% t(MASS::ginv(Ahat)) / n)

  } else if (sigma_is_fixed && !z_is_null) {
    ## ===== CASE C: beta ONLY WITH z, sigma fixed — kernel-imputation score =====
    len.beta  <- length(beta)
    sigma_use <- as.numeric(sigma_fixed)

    S_imp_full <- matrix(NA_real_, nrow = n, ncol = len.beta + 1L)
    S_ref_full <- matrix(NA_real_, nrow = n, ncol = len.beta + 1L)
    for (i in seq_len(n)) {
      s_imp <- S_beta_kern1_imp(beta, y_data[i], w_data[i], delta_data[i], z_data[i],
                                y_data, w_data, delta_data, z_data,
                                surv_c_vals, sigma_use)
      s_ref <- S_beta_f(        beta, y_data[i], w_data[i],                 z_data[i], sigma_use)
      S_imp_full[i, ] <- s_imp
      S_ref_full[i, ] <- s_ref
    }
    A_full <- -stats::var(cbind(S_imp_full, S_ref_full))[seq_len(len.beta+1), -(seq_len(len.beta+1))]
    B_full <-  stats::var(S_imp_full)
    V_full <- MASS::ginv(A_full) %*% B_full %*% t(MASS::ginv(A_full)) / n
    return(V_full[seq_len(len.beta), seq_len(len.beta), drop = FALSE])

  } else {
    ## ===== CASE D: beta ONLY WITHOUT z, sigma fixed — kernel-imputation score =====
    len.beta  <- length(beta) # expected 2
    sigma_use <- as.numeric(sigma_fixed)

    S_imp_full <- matrix(NA_real_, nrow = n, ncol = len.beta + 1L)
    S_ref_full <- matrix(NA_real_, nrow = n, ncol = len.beta + 1L)
    for (i in seq_len(n)) {
      s_imp <- S_beta_kern1_imp_nz(beta, y_data[i], w_data[i], delta_data[i],
                                   y_data, w_data, delta_data,
                                   surv_c_vals, sigma_use)
      s_ref <- S_beta_f_nz(          beta, y_data[i], w_data[i],                      sigma_use)
      S_imp_full[i, ] <- s_imp
      S_ref_full[i, ] <- s_ref
    }
    A_full <- -stats::var(cbind(S_imp_full, S_ref_full))[seq_len(len.beta+1), -(seq_len(len.beta+1))]
    B_full <-  stats::var(S_imp_full)
    V_full <- MASS::ginv(A_full) %*% B_full %*% t(MASS::ginv(A_full)) / n
    return(V_full[seq_len(len.beta), seq_len(len.beta), drop = FALSE])
  }
}
