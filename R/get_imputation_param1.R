#' Imputation estimator for (beta, sigma) with parametrically estimated \eqn{f_{X\mid Z}}
#'
#' @description
#' Solves \eqn{\sum_{i=1}^n S_\theta(Y_i,\tilde W_i,\Delta_i,Z_i;\theta)=0}
#' via \code{nleqslv}, where censored observations (\eqn{\Delta=0}) are handled by
#' imputing a single \eqn{X}-value inside the score using \eqn{E[X \mid X>W, Y, Z]}
#' under a parametric working model for \eqn{X\mid Z}.
#'
#' Supported modes:
#' - (A) estimate \eqn{(\beta,\sigma)} with \code{z_data}
#' - (B) estimate \eqn{(\beta,\sigma)} without \code{z_data}
#' - (C) estimate \eqn{\beta} with \code{z_data} and fixed \code{sigma}
#' - (D) estimate \eqn{\beta} without \code{z_data} and fixed \code{sigma}
#'
#' @section Model and nuisance estimation:
#' - With \code{z_data} present, then \eqn{Y\mid X,Z \sim \mathcal{N}\!\big(\beta_0+\beta_1 X+\beta_2 Z+\beta_3 ZX,\ \sigma^2\big)} and \eqn{X\mid Z \sim \mathrm{TN}\!\big(\alpha_{10}+\alpha_{11} Z,\ \tau_1^2;\ w_{\min}, w_{\max}\big)}.
#' - Without \code{z_data}, then \eqn{Y\mid X \sim \mathcal{N}\!\big(\beta_0+\beta_1 X,\ \sigma^2\big)} and \eqn{X \sim \mathrm{TN}\!\big(\alpha_1,\ \tau_1^2;\ w_{\min}, w_{\max}\big)}.
#' In all cases, the truncated-normal nuisance parameters \eqn{(\alpha_1,\tau_1)} for \eqn{X\mid Z} are estimated by maximum likelihood on \eqn{[w_{\min}, w_{\max}]}
#' at each iterate of \eqn{(\beta,\sigma)} (or \eqn{\beta} when \code{sigma_fixed} is supplied).
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
#' @param w_min,w_max Numeric bounds for the truncation interval used inside the nuisance likelihood.
#' @param nleqslv_args Named list of extra args passed to \code{nleqslv::nleqslv()}.
#'
#' @returns list with \code{$beta_imp_param1}, \code{$sigma_imp_param1}, \code{$root}.
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
#' # Set w_min and w_max
#' w_min = min(w_data) - sd(w_data)
#' w_max = max(w_data) + sd(w_data)
#' # CASE A: true = (0, 3, 0, 0, 4)
#' get_imputation_param1(y_data,
#'                       w_data,
#'                       delta_data,
#'                       z_data,
#'                       init = c(0.5, 2.5, 0.1, 0.1, 3.5),
#'                       w_min = w_min, w_max = w_max,
#'                       nleqslv_args = list())
#' # CASE B: true = (0, 3, 4)
#' get_imputation_param1(y_data,
#'                       w_data,
#'                       delta_data,
#'                       init = c(0.5, 2.5, 3.5),
#'                       w_min = w_min, w_max = w_max,
#'                       nleqslv_args = list())
#' # CASE C: true = (0, 3, 0, 0)
#' get_imputation_param1(y_data,
#'                       w_data,
#'                       delta_data,
#'                       z_data,
#'                       init = c(0.5, 2.5, 0.1, 0.1),
#'                       w_min = w_min, w_max = w_max,
#'                       sigma_fixed = 4,
#'                       nleqslv_args = list())
#' # CASE D: true = (0, 3)
#' get_imputation_param1(y_data,
#'                       w_data,
#'                       delta_data,
#'                       init = c(0.5, 2.5),
#'                       w_min = w_min, w_max = w_max,
#'                       sigma_fixed = 4,
#'                       nleqslv_args = list())
#' @importFrom nleqslv nleqslv
#' @importFrom stats dnorm pnorm optim
#' @importFrom truncnorm dtruncnorm ptruncnorm
#' @export
get_imputation_param1 <- function(y_data,
                                  w_data,
                                  delta_data,
                                  z_data = NULL,
                                  init,
                                  sigma_fixed = NULL,
                                  w_min, w_max,
                                  nleqslv_args = list()) {

  ## ---- Basic checks ----
  if (!is.numeric(w_min) || !is.numeric(w_max) || length(w_min)!=1L || length(w_max)!=1L || !(w_min < w_max))
    stop("`w_min` and `w_max` must be numeric scalars with w_min < w_max.", call. = FALSE)

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
    stop("`init` must be the numeric vector of FREE parameters (length >= 1).", call. = FALSE)
  if (!is.null(sigma_fixed)) {
    if (!is.numeric(sigma_fixed) || length(sigma_fixed) != 1L || !is.finite(sigma_fixed) || sigma_fixed <= 0)
      stop("`sigma_fixed` must be a positive numeric scalar.", call. = FALSE)
  }

  ## ---- Full-data score (used inside imputation) ---------------------------

  # With z: beta = (β0, β1, β2, β3); model: Y = β0 + β1 X + β2 Z + β3 ZX + ε, ε ~ N(0, σ^2)
  S_beta_f <- function(beta, y, x, z, sigma) {
    res <- y - (beta[1] + beta[2]*x + beta[3]*z + beta[4]*z*x)
    c(res * c(1, x, z, z*x) / sigma^2,
      res^2 / (2*sigma^4) - 1/(2*sigma^2))
  }
  S_beta_f <- Vectorize(S_beta_f, vectorize.args = "x")

  # No z: beta = (β0, β1)
  S_beta_f_nz <- function(beta, y, x, sigma) {
    res <- y - (beta[1] + beta[2]*x)
    c(res * c(1, x) / sigma^2,
      res^2 / (2*sigma^4) - 1/(2*sigma^2))
  }
  S_beta_f_nz <- Vectorize(S_beta_f_nz, vectorize.args = "x")

  ## ---- Imputation scores S_beta_imp (with z) and S_beta_imp_nz (no z) ----

  # With z: replace censored X by w_imp = E[X | X > w, Y, Z] under N(eta_x, v_x) truncated to [w, w_max]
  S_beta_imp <- function(beta, y, w, delta, z, alpha1_star,
                         sigma, tau1,
                         w_min, w_max) {
    v_x   <- 1/((beta[2] + beta[4]*z)^2 / sigma^2 + 1/tau1^2)
    eta_x <- v_x * ((beta[2] + beta[4]*z) * (y - (beta[1] + beta[3]*z)) / sigma^2 +
                      (alpha1_star[1] + alpha1_star[2]*z) / tau1^2)
    if (!delta) {
      if (w > w_max) {
        S_beta_f(beta, y, w, z, sigma) # any finite value; no mass above w_max
      } else {
        x_grid <- seq(w, w_max, length.out = 20)
        d      <- stats::dnorm((x_grid - eta_x) / sqrt(v_x))
        w_imp  <- sum(x_grid * d) / sum(d)  # E[X | X>w, Y, Z] under the working normal
        S_beta_f(beta, y, w_imp, z, sigma)
      }
    } else {
      S_beta_f(beta, y, w, z, sigma)
    }
  }

  # No z: analogous imputation; alpha1 is scalar
  S_beta_imp_nz <- function(beta, y, w, delta, alpha1_star,
                            sigma, tau1,
                            w_min, w_max) {
    v_x   <- 1 / ((beta[2]^2) / sigma^2 + 1 / tau1^2)
    eta_x <- v_x * (beta[2] * (y - beta[1]) / sigma^2 + alpha1_star / tau1^2)
    if (!delta) {
      if (w > w_max) {
        S_beta_f_nz(beta, y, w, sigma)
      } else {
        x_grid <- seq(w, w_max, length.out = 20)
        d      <- stats::dnorm((x_grid - eta_x) / sqrt(v_x))
        w_imp  <- sum(x_grid * d) / sum(d)
        S_beta_f_nz(beta, y, w_imp, sigma)
      }
    } else {
      S_beta_f_nz(beta, y, w, sigma)
    }
  }

  ## ---- Nuisance MLEs for X|Z (alpha1, tau1) -------------------------------

  # With z: alpha1 = (α10, α11), tau1 > 0
  find_alpha1_MLE <- function(beta, y_data, w_data, delta_data, z_data,
                              sigma, w_min, w_max) {
    log_llhd <- function(alpha1_star, tau1_star, y, w, delta, z) {
      v_x   <- 1/((beta[2] + beta[4]*z)^2 / sigma^2 + 1/tau1_star^2)
      eta_x <- v_x * ((beta[2] + beta[4]*z) * (y - (beta[1] + beta[3]*z)) / sigma^2 +
                        (alpha1_star[1] + alpha1_star[2]*z) / tau1_star^2)
      if (delta) {
        log(truncnorm::dtruncnorm(w, a = w_min, b = w_max, mean = eta_x, sd = sqrt(v_x)))
      } else {
        log(stats::pnorm(w_max, mean = eta_x, sd = sqrt(v_x)) -
              stats::pnorm(w,     mean = eta_x, sd = sqrt(v_x))) +
          stats::dnorm(
            y,
            mean = sum(c(1,
                         alpha1_star[1] + alpha1_star[2]*z,
                         z,
                         z*(alpha1_star[1] + alpha1_star[2]*z)) * beta),
            sd = sqrt(sigma^2 + (beta[2] + z*beta[4])^2 * tau1_star^2),
            log = TRUE
          ) -
          log(stats::pnorm(w_max, mean = alpha1_star[1] + alpha1_star[2]*z, sd = tau1_star) -
                stats::pnorm(w_min, mean = alpha1_star[1] + alpha1_star[2]*z, sd = tau1_star))
      }
    }
    obj <- function(par) {
      a1 <- par[1:2]; t1 <- par[3]
      -sum(mapply(log_llhd,
                  MoreArgs = list(alpha1_star = a1, tau1_star = t1),
                  y = y_data, w = w_data, delta = delta_data, z = z_data))
    }
    init <- c(mean(w_data[delta_data == 1]), 0.1, stats::sd(w_data[delta_data == 1]))
    stats::optim(init, obj)$par
  }

  # No z: alpha1 scalar, tau1 > 0
  find_alpha1_MLE_nz <- function(beta, y_data, w_data, delta_data,
                                 sigma, w_min, w_max) {
    log_llhd <- function(alpha1_star, tau1_star, y, w, delta) {
      v_x   <- 1 / ((beta[2]^2) / sigma^2 + 1 / tau1_star^2)
      eta_x <- v_x * (beta[2] * (y - beta[1]) / sigma^2 + alpha1_star / tau1_star^2)
      if (delta) {
        log(truncnorm::dtruncnorm(w, a = w_min, b = w_max, mean = eta_x, sd = sqrt(v_x)))
      } else {
        log(stats::pnorm(w_max, mean = eta_x, sd = sqrt(v_x)) -
              stats::pnorm(w,     mean = eta_x, sd = sqrt(v_x))) +
          stats::dnorm(
            y,
            mean = beta[1] + beta[2]*alpha1_star,
            sd   = sqrt(sigma^2 + (beta[2]^2) * tau1_star^2),
            log  = TRUE
          ) -
          log(stats::pnorm(w_max, mean = alpha1_star, sd = tau1_star) -
                stats::pnorm(w_min, mean = alpha1_star, sd = tau1_star))
      }
    }
    obj <- function(par) {
      a1 <- par[1]; t1 <- par[2]
      -sum(mapply(log_llhd,
                  MoreArgs = list(alpha1_star = a1, tau1_star = t1),
                  y = y_data, w = w_data, delta = delta_data))
    }
    init <- c(mean(w_data[delta_data == 1]), stats::sd(w_data[delta_data == 1]))
    stats::optim(init, obj)$par
  }

  ## ---- Estimating function: A, B, C, D (imputation) -----------------------
  pe_gauss_imp <- function(theta, y_data, w_data, delta_data, z_data, sigma_fixed, w_min, w_max) {
    z_is_null      <- is.null(z_data)
    sigma_is_fixed <- !is.null(sigma_fixed)

    if (!sigma_is_fixed && !z_is_null) {
      # ===== CASE A: estimate (beta, sigma), WITH z =====
      beta  <- theta[-length(theta)]
      sigma <- theta[length(theta)]
      a1t1  <- find_alpha1_MLE(beta, y_data, w_data, delta_data, z_data, sigma, w_min, w_max)
      alpha1_star <- a1t1[1:2]; tau1 <- a1t1[3]

      val <- rep(0, length(theta))
      for (i in seq_along(y_data)) {
        val <- val + S_beta_imp(beta, y_data[i], w_data[i], delta_data[i], z_data[i],
                                alpha1_star, sigma, tau1, w_min, w_max)
      }
      return(val)

    } else if (!sigma_is_fixed && z_is_null) {
      # ===== CASE B: estimate (beta0, beta1, sigma), WITHOUT z =====
      if (length(theta) != 3L)
        stop("Case B expects init of length 3: c(beta0, beta1, sigma).", call. = FALSE)
      beta  <- theta[1:2]
      sigma <- theta[3]
      a1t1  <- find_alpha1_MLE_nz(beta, y_data, w_data, delta_data, sigma, w_min, w_max)
      alpha1_star <- a1t1[1]; tau1 <- a1t1[2]

      val <- rep(0, 3L)
      for (i in seq_along(y_data)) {
        val <- val + S_beta_imp_nz(beta, y_data[i], w_data[i], delta_data[i],
                                   alpha1_star, sigma, tau1, w_min, w_max)
      }
      return(val)

    } else if (sigma_is_fixed && !z_is_null) {
      # ===== CASE C: estimate beta only, WITH z (sigma fixed) =====
      beta  <- theta
      sigma <- as.numeric(sigma_fixed)
      a1t1  <- find_alpha1_MLE(beta, y_data, w_data, delta_data, z_data, sigma, w_min, w_max)
      alpha1_star <- a1t1[1:2]; tau1 <- a1t1[3]

      val <- rep(0, length(theta))
      for (i in seq_along(y_data)) {
        s_i <- S_beta_imp(beta, y_data[i], w_data[i], delta_data[i], z_data[i],
                          alpha1_star, sigma, tau1, w_min, w_max)
        val <- val + s_i[seq_along(theta)] # beta components only
      }
      return(val)

    } else {
      # ===== CASE D: estimate beta only, WITHOUT z (sigma fixed) =====
      if (length(theta) != 2L)
        stop("Case D expects init of length 2: c(beta0, beta1).", call. = FALSE)
      beta  <- theta
      sigma <- as.numeric(sigma_fixed)
      a1t1  <- find_alpha1_MLE_nz(beta, y_data, w_data, delta_data, sigma, w_min, w_max)
      alpha1_star <- a1t1[1]; tau1 <- a1t1[2]

      val <- rep(0, length(theta))
      for (i in seq_along(y_data)) {
        s_i <- S_beta_imp_nz(beta, y_data[i], w_data[i], delta_data[i],
                             alpha1_star, sigma, tau1, w_min, w_max)
        val <- val + s_i[seq_along(theta)] # beta components only
      }
      return(val)
    }
  }

  ## ---- Root solve ----
  nleqslv_call <- c(list(
    x  = init,
    fn = pe_gauss_imp,
    y_data = y_data,
    w_data = w_data,
    delta_data = delta_data,
    z_data = z_data,           # may be NULL
    sigma_fixed = sigma_fixed, # may be NULL
    w_min = w_min,
    w_max = w_max
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
    sigma_imp_param1 <- as.numeric(tail(theta_hat, 1L))
    beta_imp_param1  <- as.numeric(head(theta_hat, -1L))
  } else {
    sigma_imp_param1 <- as.numeric(sigma_fixed)
    beta_imp_param1  <- as.numeric(theta_hat)
  }

  out <- list(
    beta_imp_param1  = beta_imp_param1,
    sigma_imp_param1 = sigma_imp_param1,
    root = root
  )
  class(out) <- c("imp_param1", class(out))
  out
}
