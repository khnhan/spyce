#' MLE (maximum likelihood estimator) for (beta, sigma) with nonparametrically estimated \eqn{f_{X\mid Z}}
#'
#' @description
#' Solves \eqn{\sum_{i=1}^n S_\theta(Y_i,W_i,\Delta_i,Z_i;\theta)=0} via \code{nleqslv}.
#'
#' Supported modes:
#' - (A) estimate \eqn{(\beta,\sigma)} with \code{z_data}
#' - (B) estimate \eqn{(\beta,\sigma)} without \code{z_data}
#' - (C) estimate \eqn{\beta} with \code{z_data} and fixed \code{sigma}
#' - (D) estimate \eqn{\beta} without \code{z_data} and fixed \code{sigma}
#'
#' @section Model and nuisance estimation:
#' - With \code{z_data} present, then \eqn{Y\mid X,Z \sim \mathcal{N}\!\big(\beta_0+\beta_1 X+\beta_2 Z+\beta_3 ZX,\ \sigma^2\big)} and \eqn{X\mid Z} is estimated nonparametrically.
#' - Without \code{z_data}, then \eqn{Y\mid X \sim \mathcal{N}\!\big(\beta_0+\beta_1 X,\ \sigma^2\big)} and \eqn{X} is estimated nonparametrically.
#'
#' @param y_data Numeric vector.
#' @param w_data Numeric vector of \eqn{W=\min(X,C)}.
#' @param delta_data Logical or \{0,1\} vector: \eqn{\Delta=1\{X \le C\}}.
#' @param z_data Optional numeric vector (default \code{NULL}); if provided, must match length of \code{y_data}.
#' @param surv_c_vals Numeric vector, length \code{length(y_data)}: conditional survival function estimates \eqn{S_{C\mid Y,Z}(W_i\mid Y_i,Z_i)} (or \eqn{S_{C\mid Y}(W_i\mid Y_i)} when \code{z_data} is \code{NULL}), evaluated at each observation.
#' @param init Numeric vector of free parameters:
#'   - Case A: \code{c(beta0,beta1,beta2,beta3, sigma)}
#'   - Case B: \code{c(beta0,beta1, sigma)}
#'   - Case C: \code{c(beta0,beta1,beta2,beta3)}
#'   - Case D: \code{c(beta0,beta1)}
#' @param sigma_fixed Optional positive scalar; if supplied, \eqn{\sigma} is treated as known.
#' @param nleqslv_args Named list of extra args passed to \code{nleqslv::nleqslv()}.
#'
#' @returns list with \code{$beta_MLE_nonparam1}, \code{$sigma_MLE_nonparam1}, \code{$root}.
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
#' get_MLE_nonparam1(y_data,
#'                   w_data,
#'                   delta_data,
#'                   z_data,
#'                   surv_c_vals = surv_c_vals,
#'                   init = c(0.5, 2.5, 0.1, 0.1, 3.5),
#'                   nleqslv_args = list())
#' # CASE B: true = (0, 3, 4)
#' get_MLE_nonparam1(y_data,
#'                   w_data,
#'                   delta_data,
#'                   surv_c_vals = surv_c_vals_nz,
#'                   init = c(0.5, 2.5, 3.5),
#'                   nleqslv_args = list())
#' # CASE C: true = (0, 3, 0, 0)
#' get_MLE_nonparam1(y_data,
#'                   w_data,
#'                   delta_data,
#'                   z_data,
#'                   surv_c_vals = surv_c_vals,
#'                   init = c(0.5, 2.5, 0.1, 0.1),
#'                   sigma_fixed = 4,
#'                   nleqslv_args = list())
#' # CASE D: true = (0, 3)
#' get_MLE_nonparam1(y_data,
#'                   w_data,
#'                   delta_data,
#'                   surv_c_vals = surv_c_vals_nz,
#'                   init = c(0.5, 2.5),
#'                   sigma_fixed = 4,
#'                   nleqslv_args = list())
#' @importFrom nleqslv nleqslv
#' @importFrom stats dnorm
#' @export
get_MLE_nonparam1 <- function(y_data,
                              w_data,
                              delta_data,
                              z_data = NULL,
                              surv_c_vals,
                              init,
                              sigma_fixed = NULL,
                              nleqslv_args = list()) {

  ## ---- Basic checks ----
  stopifnot(is.numeric(y_data), is.numeric(w_data))
  n <- length(y_data)
  if (!(length(w_data) == n && length(delta_data) == n))
    stop("`y_data`, `w_data`, and `delta_data` must have equal length.", call. = FALSE)
  if (!is.logical(delta_data)) {
    if (!all(delta_data %in% c(0, 1))) stop("`delta_data` must be logical or {0,1}.", call. = FALSE)
    delta_data <- as.logical(delta_data)
  }
  if (!is.null(z_data) && (!is.numeric(z_data) || length(z_data) != n))
    stop("If provided, `z_data` must be numeric and same length as `y_data`.", call. = FALSE)
  if (!is.numeric(surv_c_vals) || length(surv_c_vals) != n)
    stop("`surv_c_vals` must be a numeric vector of length equal to `y_data`.", call. = FALSE)
  if (!is.numeric(init) || length(init) < 1L)
    stop("`init` must be the numeric vector of FREE parameters (length >= 1).", call. = FALSE)
  if (!is.null(sigma_fixed)) {
    if (!is.numeric(sigma_fixed) || length(sigma_fixed) != 1L || !is.finite(sigma_fixed) || sigma_fixed <= 0)
      stop("`sigma_fixed` must be a positive numeric scalar.", call. = FALSE)
  }

  ## ---- Score kernel helpers (nonparametric f_{X|Z}) -----------------------

  S_beta_f <- function(beta, y, x, z, sigma){
    res <- y - (beta[1] + beta[2]*x + beta[3]*z + beta[4]*z*x)
    c(res * c(1, x, z, z*x) / sigma^2,
      res^2 / (2*sigma^4) - 1/(2*sigma^2))
  }
  S_beta_f <- Vectorize(S_beta_f, vectorize.args = "x")

  S_beta_kern1 <- function(beta, y, w, delta, z,
                           y_data, w_data, delta_data, z_data,
                           surv_c_vals,
                           sigma) {
    if (!delta) {
      idx_list <- which((delta_data == 1) & (w_data > w) & (z_data == z))
      if (length(idx_list) != 0L) {
        dens  <- stats::dnorm(y,
                              mean = cbind(1, w_data[idx_list], z, z*w_data[idx_list]) %*% beta,
                              sd = sigma)
        wt    <- dens / surv_c_vals[idx_list]
        num   <- sapply(idx_list, function(idx) {
          S_beta_f(beta, y, w_data[idx], z, sigma)
        }) %*% wt
        denom <- sum(wt)
        if (denom == 0) return(S_beta_f(beta, y, w, z, sigma))
        return(as.numeric(num / denom))
      } else {
        return(S_beta_f(beta, y, w, z, sigma))
      }
    } else {
      return(S_beta_f(beta, y, w, z, sigma))
    }
  }

  S_beta_f_nz <- function(beta, y, x, sigma) {
    res <- y - (beta[1] + beta[2]*x)
    c(res * c(1, x) / sigma^2,
      res^2 / (2*sigma^4) - 1/(2*sigma^2))
  }
  S_beta_f_nz <- Vectorize(S_beta_f_nz, vectorize.args = "x")

  S_beta_kern1_nz <- function(beta, y, w, delta,
                              y_data, w_data, delta_data,
                              surv_c_vals,
                              sigma) {
    if (!delta) {
      idx_list <- which((delta_data == 1) & (w_data > w))
      if (length(idx_list) != 0L) {
        dens  <- stats::dnorm(y,
                              mean = cbind(1, w_data[idx_list]) %*% beta,
                              sd = sigma)
        wt    <- dens / surv_c_vals[idx_list]
        num   <- sapply(idx_list, function(idx) {
          S_beta_f_nz(beta, y, w_data[idx], sigma)
        }) %*% wt
        denom <- sum(wt)
        if (denom == 0) return(S_beta_f_nz(beta, y, w, sigma))
        return(as.numeric(num / denom))
      } else {
        return(S_beta_f_nz(beta, y, w, sigma))
      }
    } else {
      return(S_beta_f_nz(beta, y, w, sigma))
    }
  }

  ## ---- Estimating function: A, B, C, D (kernel versions) ------------------
  pe_gauss_MLE_nonparam1 <- function(theta, y_data, w_data, delta_data, z_data,
                                     surv_c_vals, sigma_fixed) {
    z_is_null      <- is.null(z_data)
    sigma_is_fixed <- !is.null(sigma_fixed)

    if (!sigma_is_fixed && !z_is_null) {
      # ===== CASE A: estimate (beta, sigma), WITH z → kernel score =====
      beta  <- theta[-length(theta)]
      sigma <- theta[length(theta)]
      val <- rep(0, length(theta))
      for (i in seq_along(y_data)) {
        val <- val + S_beta_kern1(beta, y_data[i], w_data[i], delta_data[i], z_data[i],
                                  y_data, w_data, delta_data, z_data,
                                  surv_c_vals, sigma)
      }
      return(val)

    } else if (!sigma_is_fixed && z_is_null) {
      # ===== CASE B: estimate (beta0, beta1, sigma), WITHOUT z → kernel score =====
      if (length(theta) != 3L)
        stop("Case B expects init of length 3: c(beta0, beta1, sigma).", call. = FALSE)
      beta  <- theta[1:2]
      sigma <- theta[3]
      val <- rep(0, 3L)
      for (i in seq_along(y_data)) {
        val <- val + S_beta_kern1_nz(beta, y_data[i], w_data[i], delta_data[i],
                                     y_data, w_data, delta_data,
                                     surv_c_vals, sigma)
      }
      return(val)

    } else if (sigma_is_fixed && !z_is_null) {
      # ===== CASE C: estimate beta only, WITH z (sigma fixed) → kernel score =====
      beta  <- theta
      sigma <- as.numeric(sigma_fixed)
      val <- rep(0, length(theta))
      for (i in seq_along(y_data)) {
        s_i <- S_beta_kern1(beta, y_data[i], w_data[i], delta_data[i], z_data[i],
                            y_data, w_data, delta_data, z_data,
                            surv_c_vals, sigma)
        val <- val + s_i[seq_along(theta)]  # beta-score components only
      }
      return(val)

    } else {
      # ===== CASE D: estimate beta only, WITHOUT z (sigma fixed) → kernel score =====
      if (length(theta) != 2L)
        stop("Case D expects init of length 2: c(beta0, beta1).", call. = FALSE)
      beta  <- theta
      sigma <- as.numeric(sigma_fixed)
      val <- rep(0, length(theta))
      for (i in seq_along(y_data)) {
        s_i <- S_beta_kern1_nz(beta, y_data[i], w_data[i], delta_data[i],
                               y_data, w_data, delta_data,
                               surv_c_vals, sigma)
        val <- val + s_i[seq_along(theta)]  # beta-score components only
      }
      return(val)
    }
  }

  ## ---- Root solve ----
  nleqslv_call <- c(list(
    x  = init,
    fn = pe_gauss_MLE_nonparam1,
    y_data = y_data,
    w_data = w_data,
    delta_data = delta_data,
    z_data = z_data,           # may be NULL
    surv_c_vals = surv_c_vals, # required
    sigma_fixed = sigma_fixed  # may be NULL
  ), nleqslv_args)

  root <- try(do.call(nleqslv::nleqslv, nleqslv_call), silent = TRUE)
  if (inherits(root, "try-error")) {
    stop("nleqslv failed: ", conditionMessage(attr(root, "condition")), call. = FALSE)
  }
  if (!is.null(root$termcd) && root$termcd > 2) {
    warning("nleqslv did not report strong convergence; inspect `root`.", call. = FALSE)
  }

  theta_hat <- root$x
  if (is.null(sigma_fixed)) {
    sigma_MLE_nonparam1 <- as.numeric(tail(theta_hat, 1L))
    beta_MLE_nonparam1  <- as.numeric(head(theta_hat, -1L))
  } else {
    sigma_MLE_nonparam1 <- as.numeric(sigma_fixed)
    beta_MLE_nonparam1  <- as.numeric(theta_hat)
  }

  out <- list(
    beta_MLE_nonparam1  = beta_MLE_nonparam1,
    sigma_MLE_nonparam1 = sigma_MLE_nonparam1,
    root = root
  )
  class(out) <- c("mle_nonparam1", class(out))
  out
}

