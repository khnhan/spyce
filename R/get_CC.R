#' CC (complete-case) estimator for (beta, sigma)
#'
#' @description
#' Solves the complete-case estimating equations
#' \deqn{\sum_{i=1}^n \Delta_i\, S_\theta^F\!\big(Y_i, W_i, Z_i;\theta\big)=0}
#' via \code{nleqslv}, where \eqn{\Delta_i=1\{X_i \le C_i\}\) and \(W_i=\min(X_i,C_i)}.
#'
#' Supported modes:
#' - (A) estimate \eqn{(\beta,\sigma)} with \code{z_data}
#' - (B) estimate \eqn{(\beta,\sigma)} without \code{z_data}
#' - (C) estimate \eqn{\beta} with \code{z_data} and fixed \code{sigma}
#' - (D) estimate \eqn{\beta} without \code{z_data} and fixed \code{sigma}
#'
#' @section Model and nuisance estimation:
#' - With \code{z_data} present, then \eqn{Y\mid X,Z \sim \mathcal{N}\!\big(\beta_0+\beta_1 X+\beta_2 Z+\beta_3 ZX,\ \sigma^2\big)}.
#' - Without \code{z_data}, then \eqn{Y\mid X \sim \mathcal{N}\!\big(\beta_0+\beta_1 X,\ \sigma^2\big)}.
#' This complete-case estimator does **not** estimate \eqn{\eta_1} or \eqn{\eta_2}; it uses only observations with \eqn{\Delta=1} and plugs \eqn{X=W}.
#'
#' @param y_data Numeric vector.
#' @param w_data Numeric vector of \eqn{W=\min(X,C)}.
#' @param delta_data Logical or \{0,1\} vector: \eqn{\Delta=1\{X \le C\}}.
#' @param z_data Optional numeric vector (default \code{NULL}); if provided, must match length of \code{y_data}.
#' @param init Numeric vector of free parameters:
#'   \itemize{
#'     \item Case A: \code{c(beta0,beta1,beta2,beta3, sigma)}
#'     \item Case B: \code{c(beta0,beta1, sigma)}
#'     \item Case C: \code{c(beta0,beta1,beta2,beta3)}
#'     \item Case D: \code{c(beta0,beta1)}
#'   }
#' @param sigma_fixed Optional positive scalar; if supplied, \eqn{\sigma} is treated as known.
#' @param nleqslv_args Named list passed to \code{nleqslv::nleqslv()}.
#'
#' @returns list with \code{$beta_CC}, \code{$sigma_CC}, \code{$root}.
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
#' # CASE A: true = (0, 3, 0, 0, 4)
#' get_CC(y_data,
#'        w_data,
#'        delta_data,
#'        z_data,
#'        init = c(0.5, 2.5, 0.1, 0.1, 3.5),
#'        nleqslv_args = list())
#' # CASE B: true = (0, 3, 4)
#' get_CC(y_data,
#'        w_data,
#'        delta_data,
#'        init = c(0.5, 2.5, 3.5),
#'        nleqslv_args = list())
#' # CASE C: true = (0, 3, 0, 0)
#' get_CC(y_data,
#'        w_data,
#'        delta_data,
#'        z_data,
#'        init = c(0.5, 2.5, 0.1, 0.1),
#'        sigma_fixed = 4,
#'        nleqslv_args = list())
#' # CASE D: true = (0, 3)
#' get_CC(y_data,
#'        w_data,
#'        delta_data,
#'        init = c(0.5, 2.5),
#'        sigma_fixed = 4,
#'        nleqslv_args = list())
#'
#' @importFrom nleqslv nleqslv
#' @export
get_CC <- function(y_data,
                   w_data,
                   delta_data,
                   z_data = NULL,
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
  if (!is.numeric(init) || length(init) < 1L)
    stop("`init` must be the numeric vector of FREE parameters.", call. = FALSE)
  if (!is.null(sigma_fixed)) {
    if (!is.numeric(sigma_fixed) || length(sigma_fixed) != 1L || !is.finite(sigma_fixed) || sigma_fixed <= 0)
      stop("`sigma_fixed` must be a positive numeric scalar.", call. = FALSE)
  }

  ## ---- Full-data scores (evaluate at X = W; no nuisance) ------------------

  # With z: beta = (β0, β1, β2, β3); model: Y = β0 + β1 X + β2 Z + β3 ZX + ε, ε ~ N(0, σ^2)
  S_beta_f <- function(beta, y, x, z, sigma) {
    res <- y - (beta[1] + beta[2]*x + beta[3]*z + beta[4]*z*x)
    c(res * c(1, x, z, z*x) / sigma^2,
      res^2 / (2*sigma^4) - 1/(2*sigma^2))
  }
  S_beta_f <- Vectorize(S_beta_f, vectorize.args = "x")

  # No z: beta = (β0, β1); model: Y = β0 + β1 X + ε, ε ~ N(0, σ^2)
  S_beta_f_nz <- function(beta, y, x, sigma) {
    res <- y - (beta[1] + beta[2]*x)
    c(res * c(1, x) / sigma^2,
      res^2 / (2*sigma^4) - 1/(2*sigma^2))
  }
  S_beta_f_nz <- Vectorize(S_beta_f_nz, vectorize.args = "x")

  ## ---- Estimating function for CC -----------------------------------------
  pe_gauss_CC <- function(theta, y_data, w_data, delta_data, z_data, sigma_fixed) {
    z_is_null      <- is.null(z_data)
    sigma_is_fixed <- !is.null(sigma_fixed)

    if (!sigma_is_fixed && !z_is_null) {
      # CASE A: estimate (beta, sigma) WITH z
      beta  <- theta[-length(theta)]
      sigma <- theta[length(theta)]
      val <- rep(0, length(theta))
      for (i in which(delta_data)) {
        val <- val + S_beta_f(beta, y = y_data[i], x = w_data[i], z = z_data[i], sigma = sigma)
      }
      return(val)

    } else if (!sigma_is_fixed && z_is_null) {
      # CASE B: estimate (beta0, beta1, sigma) WITHOUT z
      if (length(theta) != 3L)
        stop("Case B expects init of length 3: c(beta0, beta1, sigma).", call. = FALSE)
      beta  <- theta[1:2]
      sigma <- theta[3]
      val <- rep(0, 3L)
      for (i in which(delta_data)) {
        val <- val + S_beta_f_nz(beta, y = y_data[i], x = w_data[i], sigma = sigma)
      }
      return(val)

    } else if (sigma_is_fixed && !z_is_null) {
      # CASE C: estimate beta ONLY WITH z (sigma fixed)
      beta  <- theta
      sigma <- as.numeric(sigma_fixed)
      val <- rep(0, length(beta))
      for (i in which(delta_data)) {
        s_full <- S_beta_f(beta, y = y_data[i], x = w_data[i], z = z_data[i], sigma = sigma)
        val    <- val + s_full[seq_along(beta)]  # beta components only
      }
      return(val)

    } else {
      # CASE D: estimate beta ONLY WITHOUT z (sigma fixed)
      if (length(theta) != 2L)
        stop("Case D expects init of length 2: c(beta0, beta1).", call. = FALSE)
      beta  <- theta
      sigma <- as.numeric(sigma_fixed)
      val <- rep(0, length(beta))
      for (i in which(delta_data)) {
        s_full <- S_beta_f_nz(beta, y = y_data[i], x = w_data[i], sigma = sigma)
        val    <- val + s_full[seq_along(beta)]  # beta components only
      }
      return(val)
    }
  }

  ## ---- Root solve ----------------------------------------------------------
  nleqslv_call <- c(list(
    x  = init,
    fn = pe_gauss_CC,
    y_data = y_data,
    w_data = w_data,
    delta_data = delta_data,
    z_data = z_data,
    sigma_fixed = sigma_fixed
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
    sigma_CC <- as.numeric(tail(theta_hat, 1L))
    beta_CC  <- as.numeric(head(theta_hat, -1L))
  } else {
    sigma_CC <- as.numeric(sigma_fixed)
    beta_CC  <- as.numeric(theta_hat)
  }

  out <- list(beta_CC = beta_CC, sigma_CC = sigma_CC, root = root)
  class(out) <- c("cc_param1", class(out))
  out
}
