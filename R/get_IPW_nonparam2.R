#' IPW (inverse probability weighting) estimator for (beta, sigma) with nonparametrically estimated \eqn{f_{C\mid Y,Z}}
#'
#' @description
#' Solves \eqn{\sum_{i=1}^n \frac{\Delta_i}{S_{C\mid Y,Z}(W_i\mid Y_i,Z_i)}S_\theta(Y_i,W_i,\Delta_i,Z_i;\theta)=0} via \code{nleqslv}.
#'
#' Supported modes:
#' - (A) estimate \eqn{(\beta,\sigma)} with \code{z_data}
#' - (B) estimate \eqn{(\beta,\sigma)} without \code{z_data}
#' - (C) estimate \eqn{\beta} with \code{z_data} and fixed \code{sigma}
#' - (D) estimate \eqn{\beta} without \code{z_data} and fixed \code{sigma}
#'
#' @section Model and nuisance estimation:
#' - With \code{z_data} present, \eqn{Y\mid X,Z \sim \mathcal{N}\!\big(\beta_0+\beta_1 X+\beta_2 Z+\beta_3 ZX,\ \sigma^2\big)} and the conditional survival \eqn{S_{C\mid Y,Z}(W\mid Y,Z)} is supplied nonparametrically via \code{surv_c_vals}.
#' - Without \code{z_data}, \eqn{Y\mid X \sim \mathcal{N}\!\big(\beta_0+\beta_1 X,\ \sigma^2\big)} and \code{surv_c_vals} supplies \eqn{S_{C\mid Y}(W\mid Y)}.
#'
#' @param y_data Numeric vector.
#' @param w_data Numeric vector of \eqn{W=\min(X,C)}.
#' @param delta_data Logical or \{0,1\} vector: \eqn{\Delta=1\{X \le C\}}.
#' @param z_data Optional numeric vector (default \code{NULL}); if provided, must match length of \code{y_data}.
#' @param surv_c_vals Numeric vector, length \code{length(y_data)}: conditional survival \eqn{S_{C\mid Y,Z}(W_i\mid Y_i,Z_i)}
#'   (or \eqn{S_{C\mid Y}(W_i\mid Y_i)} when \code{z_data} is \code{NULL}), evaluated at each observation.
#' @param init Numeric vector of free parameters:
#'   - Case A: \code{c(beta0,beta1,beta2,beta3, sigma)}
#'   - Case B: \code{c(beta0,beta1, sigma)}
#'   - Case C: \code{c(beta0,beta1,beta2,beta3)}
#'   - Case D: \code{c(beta0,beta1)}
#' @param sigma_fixed Optional positive scalar; if supplied, \eqn{\sigma} is treated as known.
#' @param nleqslv_args Named list of extra args passed to \code{nleqslv::nleqslv()}.
#'
#' @returns list with \code{$beta_IPW_nonparam2}, \code{$sigma_IPW_nonparam2}, \code{$root}.
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
#' get_IPW_nonparam2(y_data,
#'                   w_data,
#'                   delta_data,
#'                   z_data,
#'                   surv_c_vals = surv_c_vals,
#'                   init = c(0.5, 2.5, 0.1, 0.1, 3.5),
#'                   nleqslv_args = list())
#' # CASE B: true = (0, 3, 4)
#' get_IPW_nonparam2(y_data,
#'                   w_data,
#'                   delta_data,
#'                   surv_c_vals = surv_c_vals_nz,
#'                   init = c(0.5, 2.5, 3.5),
#'                   nleqslv_args = list())
#' # CASE C: true = (0, 3, 0, 0)
#' get_IPW_nonparam2(y_data,
#'                   w_data,
#'                   delta_data,
#'                   z_data,
#'                   surv_c_vals = surv_c_vals,
#'                   init = c(0.5, 2.5, 0.1, 0.1),
#'                   sigma_fixed = 4,
#'                   nleqslv_args = list())
#' # CASE D: true = (0, 3)
#' get_IPW_nonparam2(y_data,
#'                   w_data,
#'                   delta_data,
#'                   surv_c_vals = surv_c_vals_nz,
#'                   init = c(0.5, 2.5),
#'                   sigma_fixed = 4,
#'                   nleqslv_args = list())
#' @importFrom nleqslv nleqslv
#' @importFrom stats dnorm
#' @export
get_IPW_nonparam2 <- function(y_data,
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
    stop("`init` must be the numeric vector of FREE parameters.", call. = FALSE)
  if (!is.null(sigma_fixed)) {
    if (!is.numeric(sigma_fixed) || length(sigma_fixed) != 1L || !is.finite(sigma_fixed) || sigma_fixed <= 0)
      stop("`sigma_fixed` must be a positive numeric scalar.", call. = FALSE)
  }

  ## ---- Full-data scores evaluated at X = W (for complete cases) -----------
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

  ## ---- IPW estimating function (kernel survival weights) ------------------
  pe_gauss_IPW_kernel2 <- function(theta, y_data, w_data, delta_data, z_data,
                                   surv_c_vals, sigma_fixed) {
    sigma_is_fixed <- !is.null(sigma_fixed)

    if (!sigma_is_fixed) {
      # ===== CASE A: estimate (beta, sigma) WITH z =====
      beta  <- theta[-length(theta)]
      sigma <- theta[length(theta)]
      val <- rep(0, length(theta))
      for (i in which(delta_data)) {
        w <- max(surv_c_vals[i], .Machine$double.eps)
        val <- val + S_beta_f(beta, y = y_data[i], x = w_data[i], z = z_data[i], sigma = sigma) / w
      }
      return(val)

    } else {
      # ===== CASE C: estimate beta ONLY WITH z (sigma fixed) =====
      beta  <- theta
      sigma <- as.numeric(sigma_fixed)
      val <- rep(0, length(beta))
      for (i in which(delta_data)) {
        w <- max(surv_c_vals[i], .Machine$double.eps)
        s_full <- S_beta_f(beta, y = y_data[i], x = w_data[i], z = z_data[i], sigma = sigma) / w
        val    <- val + s_full[seq_along(beta)]  # beta components only
      }
      return(val)
    }
  }

  pe_gauss_IPW_kernel2_nz <- function(theta, y_data, w_data, delta_data,
                                      surv_c_vals, sigma_fixed) {
    sigma_is_fixed <- !is.null(sigma_fixed)

    if (!sigma_is_fixed) {
      # ===== CASE B: estimate (beta, sigma) WITHOUT z =====
      if (length(theta) != 3L)
        stop("Case B expects init of length 3: c(beta0, beta1, sigma).", call. = FALSE)
      beta  <- theta[1:2]
      sigma <- theta[3]
      val <- rep(0, 3L)
      for (i in which(delta_data)) {
        w <- max(surv_c_vals[i], .Machine$double.eps)
        val <- val + S_beta_f_nz(beta, y = y_data[i], x = w_data[i], sigma = sigma) / w
      }
      return(val)

    } else {
      # ===== CASE D: estimate beta ONLY WITHOUT z (sigma fixed) =====
      if (length(theta) != 2L)
        stop("Case D expects init of length 2: c(beta0, beta1).", call. = FALSE)
      beta  <- theta
      sigma <- as.numeric(sigma_fixed)
      val <- rep(0, length(beta))
      for (i in which(delta_data)) {
        w <- max(surv_c_vals[i], .Machine$double.eps)
        s_full <- S_beta_f_nz(beta, y = y_data[i], x = w_data[i], sigma = sigma) / w
        val    <- val + s_full[seq_along(beta)]  # beta components only
      }
      return(val)
    }
  }

  ## ---- Root solve ----
  nleqslv_call <- if (is.null(z_data)) {
    c(list(
      x  = init,
      fn = pe_gauss_IPW_kernel2_nz,
      y_data = y_data,
      w_data = w_data,
      delta_data = delta_data,
      surv_c_vals = surv_c_vals,
      sigma_fixed = sigma_fixed
    ), nleqslv_args)
  } else {
    c(list(
      x  = init,
      fn = pe_gauss_IPW_kernel2,
      y_data = y_data,
      w_data = w_data,
      delta_data = delta_data,
      z_data = z_data,
      surv_c_vals = surv_c_vals,
      sigma_fixed = sigma_fixed
    ), nleqslv_args)
  }

  root <- try(do.call(nleqslv::nleqslv, nleqslv_call), silent = TRUE)
  if (inherits(root, "try-error")) {
    stop("nleqslv failed: ", conditionMessage(attr(root, "condition")), call. = FALSE)
  }
  if (!is.null(root$termcd) && root$termcd > 2) {
    warning("nleqslv did not report strong convergence; inspect `root`.", call. = FALSE)
  }

  theta_hat <- root$x
  if (is.null(sigma_fixed)) {
    sigma_IPW_nonparam2 <- as.numeric(tail(theta_hat, 1L))
    beta_IPW_nonparam2  <- as.numeric(head(theta_hat, -1L))
  } else {
    sigma_IPW_nonparam2 <- as.numeric(sigma_fixed)
    beta_IPW_nonparam2  <- as.numeric(theta_hat)
  }

  out <- list(
    beta_IPW_nonparam2  = beta_IPW_nonparam2,
    sigma_IPW_nonparam2 = sigma_IPW_nonparam2,
    root = root
  )
  class(out) <- c("ipw_nonparam2", class(out))
  out
}



