#' Plug-in variance for MLE (maximum likelihood estimator) of (beta, sigma) with parametrically estimated \eqn{f_{X\mid Z}}
#'
#' @description
#' Computes the plug-in variance–covariance matrix for the MLE obtained by
#' \code{\link{get_MLE_param1}}. Supports four modes matching that estimator:
#' \itemize{
#'   \item \strong{CASE A}: estimate \eqn{(\beta,\sigma)} with \code{z_data}
#'   \item \strong{CASE B}: estimate \eqn{(\beta,\sigma)} without \code{z_data}
#'   \item \strong{CASE C}: estimate \eqn{\beta} with \code{z_data} and fixed \code{sigma}
#'   \item \strong{CASE D}: estimate \eqn{\beta} without \code{z_data} and fixed \code{sigma}
#' }
#' It is recommended to use the same \code{w_min} and \code{w_max} as in \code{\link{get_MLE_param1}}.
#'
#' @param beta Numeric vector of estimated \eqn{\beta}.
#' @param y_data Numeric vector.
#' @param w_data Numeric vector of \eqn{W=\min(X,C)}.
#' @param delta_data Logical or \{0,1\} vector: \eqn{\Delta=1\{X \le C\}}.
#' @param z_data Optional numeric vector (default \code{NULL}); if provided, must match length of \code{y_data}.
#' @param sigma Numeric scalar \eqn{\hat\sigma}. Ignored when \code{sigma_fixed} is provided.
#' @param sigma_fixed Optional positive scalar. If supplied, \eqn{\sigma} is treated as known
#'   and the returned variance is for \eqn{\beta} only.
#' @param w_min,w_max Numeric bounds used internally (same role as in \code{\link{get_MLE_param1}}).
#'
#' @returns
#' A variance–covariance matrix.
#' \itemize{
#'   \item CASE A/B (\code{sigma_fixed = NULL}): matrix of size \code{length(beta)+1}.
#'   \item CASE C/D (\code{sigma_fixed} provided): matrix of size \code{length(beta)}.
#' }
#'
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
#' variance_MLE_param1(beta = c(0.5, 2.5, 0.1, 0.1),
#'                     y_data,
#'                     w_data,
#'                     delta_data,
#'                     z_data,
#'                     sigma = 3.5,
#'                     w_min = w_min, w_max = w_max)
#' # CASE B: true = (0, 3, 4)
#' variance_MLE_param1(beta = c(0.5, 2.5),
#'                     y_data,
#'                     w_data,
#'                     delta_data,
#'                     sigma = 3.5,
#'                     w_min = w_min, w_max = w_max)
#' # CASE C: true = (0, 3, 0, 0)
#' variance_MLE_param1(beta = c(0.5, 2.5, 0.1, 0.1),
#'                     y_data,
#'                     w_data,
#'                     delta_data,
#'                     z_data,
#'                     sigma_fixed = 4,
#'                     w_min = w_min, w_max = w_max)
#' # CASE D: true = (0, 3)
#' variance_MLE_param1(beta = c(0.5, 2.5),
#'                     y_data,
#'                     w_data,
#'                     delta_data,
#'                     sigma_fixed = 4,
#'                     w_min = w_min, w_max = w_max)
#' @importFrom MASS ginv
#' @importFrom stats dnorm pnorm optim var
#' @importFrom truncnorm dtruncnorm ptruncnorm
#' @export
variance_MLE_param1 <- function(beta,
                                y_data,
                                w_data,
                                delta_data,
                                z_data = NULL,
                                sigma,
                                sigma_fixed = NULL,
                                w_min, w_max) {

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
  if (!is.null(sigma_fixed)) {
    if (!is.numeric(sigma_fixed) || length(sigma_fixed) != 1L || !is.finite(sigma_fixed) || sigma_fixed <= 0) {
      stop("`sigma_fixed` must be a positive numeric scalar.", call. = FALSE)
    }
  }

  ## ---- Local helpers (WITH z) ---------------------------------------------
  S_beta_f <- function(beta, y, x, z, sigma) {
    res <- y - (beta[1] + beta[2]*x + beta[3]*z + beta[4]*z*x)
    c(res * c(1, x, z, z*x) / sigma^2,
      res^2 / (2 * sigma^4) - 1 / (2 * sigma^2))
  }
  S_beta_f <- Vectorize(S_beta_f, vectorize.args = "x")

  S_beta <- function(beta, y, w, delta, z, alpha1_star,
                     sigma, tau1,
                     w_min, w_max) {
    v_x   <- 1 / ((beta[2] + beta[4]*z)^2 / sigma^2 + 1 / tau1^2)
    eta_x <- v_x * ((beta[2] + beta[4]*z) * (y - (beta[1] + beta[3]*z)) / sigma^2 +
                      (alpha1_star[1] + alpha1_star[2]*z) / tau1^2)

    if (!delta) {
      if (w > w_max) {
        S_beta_f(beta, y, w, z, sigma)
      } else {
        x_norm <- (seq(w, w_max, length.out = 20) - eta_x) / sqrt(v_x)
        d      <- stats::dnorm(x_norm)
        num    <- S_beta_f(beta, y, x_norm * sqrt(v_x) + eta_x, z, sigma) %*% d
        denom  <- sum(d)
        if (denom == 0) return(S_beta_f(beta, y, w, z, sigma))
        as.numeric(num / denom)
      }
    } else {
      S_beta_f(beta, y, w, z, sigma)
    }
  }

  find_alpha1_MLE <- function(beta, y_data, w_data, delta_data, z_data,
                              sigma, w_min, w_max) {
    log_llhd <- function(alpha1_star, tau1_star, y, w, delta, z) {
      v_x   <- 1 / ((beta[2] + beta[4]*z)^2 / sigma^2 + 1 / tau1_star^2)
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

    log_llhd_sum <- function(alpha1tau1) {
      alpha1 <- alpha1tau1[1:2]
      tau1   <- alpha1tau1[3]
      llhd   <- 0
      for (i in seq_along(y_data)) {
        llhd <- llhd + log_llhd(alpha1, tau1, y_data[i], w_data[i], delta_data[i], z_data[i])
      }
      -llhd
    }

    init <- c(mean(w_data[delta_data == 1]), 0.1, stats::sd(w_data[delta_data == 1]))
    stats::optim(init, log_llhd_sum)$par
  }

  find_alpha2_MLE <- function(y_data, w_data, delta_data, z_data,
                              w_min, w_max) {
    log_llhd <- function(alpha2_star, tau2_star, y, w, delta, z) {
      mu <- sum(c(1, y, z) * alpha2_star)
      if (delta) {
        log(1 - truncnorm::ptruncnorm(w, a = w_min, b = w_max, mean = mu, sd = tau2_star))
      } else {
        log(truncnorm::dtruncnorm(w, a = w_min, b = w_max, mean = mu, sd = tau2_star))
      }
    }

    log_llhd_sum <- function(alpha2tau2) {
      alpha2 <- alpha2tau2[1:3]
      tau2   <- alpha2tau2[4]
      llhd   <- 0
      for (i in seq_along(y_data)) {
        llhd <- llhd + log_llhd(alpha2, tau2, y_data[i], w_data[i], delta_data[i], z_data[i])
      }
      -llhd
    }

    init <- c(mean(w_data[delta_data == 0]), 0.1, 0.1, stats::sd(w_data[delta_data == 0]))
    stats::optim(init, log_llhd_sum)$par
  }

  ## ---- Local helpers (NO z) -----------------------------------------------
  S_beta_f_nz <- function(beta, y, x, sigma) {
    res <- y - (beta[1] + beta[2] * x)
    c(res * c(1, x) / sigma^2,
      res^2 / (2 * sigma^4) - 1 / (2 * sigma^2))
  }
  S_beta_f_nz <- Vectorize(S_beta_f_nz, vectorize.args = "x")

  S_beta_nz <- function(beta, y, w, delta, alpha1_star,
                        sigma, tau1,
                        w_min, w_max) {
    v_x   <- 1 / ((beta[2]^2) / sigma^2 + 1 / tau1^2)
    eta_x <- v_x * (beta[2] * (y - beta[1]) / sigma^2 + alpha1_star / tau1^2)

    if (!delta) {
      if (w > w_max) {
        S_beta_f_nz(beta, y, w, sigma)
      } else {
        x_norm <- (seq(w, w_max, length.out = 20) - eta_x) / sqrt(v_x)
        d      <- stats::dnorm(x_norm)
        num    <- S_beta_f_nz(beta, y, x_norm * sqrt(v_x) + eta_x, sigma) %*% d
        denom  <- sum(d)
        if (denom == 0) return(S_beta_f_nz(beta, y, w, sigma))
        as.numeric(num / denom)
      }
    } else {
      S_beta_f_nz(beta, y, w, sigma)
    }
  }

  find_alpha1_MLE_nz <- function(beta, y_data, w_data, delta_data,
                                 sigma, w_min, w_max) {
    log_llhd <- function(alpha1_star, tau1_star, y, w, delta) {
      v_x   <- 1 / ((beta[2]^2) / sigma^2 + 1 / tau1_star^2)
      eta_x <- v_x * (beta[2] * (y - beta[1]) / sigma^2 + alpha1_star / tau1_star^2)

      if (delta) {
        log(truncnorm::dtruncnorm(w, a = w_min, b = w_max,
                                  mean = eta_x, sd = sqrt(v_x)))
      } else {
        log(stats::pnorm(w_max, mean = eta_x, sd = sqrt(v_x)) -
              stats::pnorm(w,     mean = eta_x, sd = sqrt(v_x))) +
          stats::dnorm(
            y,
            mean = beta[1] + beta[2] * alpha1_star,
            sd   = sqrt(sigma^2 + (beta[2]^2) * tau1_star^2),
            log  = TRUE
          ) -
          log(stats::pnorm(w_max, mean = alpha1_star, sd = tau1_star) -
                stats::pnorm(w_min, mean = alpha1_star, sd = tau1_star))
      }
    }

    log_llhd_sum <- function(par) {
      alpha1 <- par[1]
      tau1   <- par[2]
      llhd   <- 0
      for (i in seq_along(y_data)) {
        llhd <- llhd + log_llhd(alpha1, tau1, y_data[i], w_data[i], delta_data[i])
      }
      -llhd
    }

    init <- c(mean(w_data[delta_data == 1]), stats::sd(w_data[delta_data == 1]))
    stats::optim(init, log_llhd_sum)$par
  }

  ## ---- Case dispatch -------------------------------------------------------
  z_is_null      <- is.null(z_data)
  sigma_is_fixed <- !is.null(sigma_fixed)

  if (!sigma_is_fixed && !z_is_null) {
    ## ===== CASE A: (beta, sigma) WITH z =====
    len.beta <- length(beta)

    get_Sigma <- function(beta, y_data, w_data, delta_data, z_data,
                          alpha1_star, alpha2_star,
                          sigma, tau1, tau2,
                          w_min, w_max) {
      val <- matrix(NA_real_, nrow = n, ncol = len.beta + 1)
      for (i in seq_len(n)) {
        val[i, ] <- S_beta(beta,
                           y = y_data[i],
                           w = w_data[i],
                           delta = delta_data[i],
                           z = z_data[i],
                           alpha1_star = alpha1_star,
                           sigma = sigma, tau1 = tau1,
                           w_min = w_min, w_max = w_max)
      }
      stats::var(val)
    }

    alpha1tau1 <- find_alpha1_MLE(beta, y_data, w_data, delta_data, z_data,
                                  sigma, w_min, w_max)
    alpha1_MLE <- alpha1tau1[1:2]
    tau1       <- alpha1tau1[3]

    alpha2tau2 <- find_alpha2_MLE(y_data, w_data, delta_data, z_data,
                                  w_min, w_max)
    alpha2_MLE <- alpha2tau2[1:3]
    tau2       <- alpha2tau2[4]

    Sigma <- get_Sigma(beta, y_data, w_data, delta_data, z_data,
                       alpha1_star = alpha1_MLE, alpha2_star = alpha2_MLE,
                       sigma = sigma, tau1 = tau1, tau2 = tau2,
                       w_min = w_min, w_max = w_max)

    return(MASS::ginv(Sigma) / n)

  } else if (!sigma_is_fixed && z_is_null) {
    ## ===== CASE B: (beta, sigma) WITHOUT z  — length(beta) should be 2 =====
    len.beta <- length(beta)  # expected 2

    alpha1tau1 <- find_alpha1_MLE_nz(beta, y_data, w_data, delta_data,
                                     sigma, w_min, w_max)
    alpha1_MLE <- alpha1tau1[1]
    tau1       <- alpha1tau1[2]

    get_Sigma_nz <- function(beta, y_data, w_data, delta_data,
                             alpha1_star,
                             sigma, tau1,
                             w_min, w_max) {
      val <- matrix(NA_real_, nrow = n, ncol = len.beta + 1)  # (β0,β1,σ)
      for (i in seq_len(n)) {
        val[i, ] <- S_beta_nz(beta,
                              y = y_data[i],
                              w = w_data[i],
                              delta = delta_data[i],
                              alpha1_star = alpha1_star,
                              sigma = sigma, tau1 = tau1,
                              w_min = w_min, w_max = w_max)
      }
      stats::var(val)
    }

    Sigma <- get_Sigma_nz(beta, y_data, w_data, delta_data,
                          alpha1_star = alpha1_MLE,
                          sigma = sigma, tau1 = tau1,
                          w_min = w_min, w_max = w_max)

    return(MASS::ginv(Sigma) / n)

  } else if (sigma_is_fixed && !z_is_null) {
    ## ===== CASE C: beta ONLY WITH z, sigma fixed =====
    len.beta  <- length(beta)
    sigma_use <- as.numeric(sigma_fixed)

    get_Sigma <- function(beta, y_data, w_data, delta_data, z_data,
                          alpha1_star, alpha2_star,
                          sigma, tau1, tau2,
                          w_min, w_max) {
      val <- matrix(NA_real_, nrow = n, ncol = len.beta + 1)
      for (i in seq_len(n)) {
        val[i, ] <- S_beta(beta,
                           y = y_data[i],
                           w = w_data[i],
                           delta = delta_data[i],
                           z = z_data[i],
                           alpha1_star = alpha1_star,
                           sigma = sigma, tau1 = tau1,
                           w_min = w_min, w_max = w_max)
      }
      stats::var(val)
    }

    alpha1tau1 <- find_alpha1_MLE(beta, y_data, w_data, delta_data, z_data,
                                  sigma_use, w_min, w_max)
    alpha1_MLE <- alpha1tau1[1:2]
    tau1       <- alpha1tau1[3]

    alpha2tau2 <- find_alpha2_MLE(y_data, w_data, delta_data, z_data,
                                  w_min, w_max)
    alpha2_MLE <- alpha2tau2[1:3]
    tau2       <- alpha2tau2[4]

    Sigma_full <- get_Sigma(beta, y_data, w_data, delta_data, z_data,
                            alpha1_star = alpha1_MLE, alpha2_star = alpha2_MLE,
                            sigma = sigma_use, tau1 = tau1, tau2 = tau2,
                            w_min = w_min, w_max = w_max)

    V_full <- MASS::ginv(Sigma_full) / n
    return(V_full[seq_len(len.beta), seq_len(len.beta), drop = FALSE])

  } else {
    ## ===== CASE D: beta ONLY WITHOUT z, sigma fixed =====
    len.beta  <- length(beta)           # expected 2
    sigma_use <- as.numeric(sigma_fixed)

    alpha1tau1 <- find_alpha1_MLE_nz(beta, y_data, w_data, delta_data,
                                     sigma_use, w_min, w_max)
    alpha1_MLE <- alpha1tau1[1]
    tau1       <- alpha1tau1[2]

    get_Sigma_nz <- function(beta, y_data, w_data, delta_data,
                             alpha1_star,
                             sigma, tau1,
                             w_min, w_max) {
      val <- matrix(NA_real_, nrow = n, ncol = len.beta + 1)  # (β0,β1,σ)
      for (i in seq_len(n)) {
        val[i, ] <- S_beta_nz(beta,
                              y = y_data[i],
                              w = w_data[i],
                              delta = delta_data[i],
                              alpha1_star = alpha1_star,
                              sigma = sigma, tau1 = tau1,
                              w_min = w_min, w_max = w_max)
      }
      stats::var(val)
    }

    Sigma_full <- get_Sigma_nz(beta, y_data, w_data, delta_data,
                               alpha1_star = alpha1_MLE,
                               sigma = sigma_use, tau1 = tau1,
                               w_min = w_min, w_max = w_max)

    V_full <- MASS::ginv(Sigma_full) / n
    return(V_full[seq_len(len.beta), seq_len(len.beta), drop = FALSE])
  }
}
