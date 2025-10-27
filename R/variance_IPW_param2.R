#' Plug-in variance for IPW (inverse probability weighting) estimator of (beta, sigma) with parametrically estimated \eqn{f_{C\mid Y,Z}}
#'
#' @description
#' Computes the plug-in variance–covariance matrix for the IPW estimator
#' \code{\link{get_IPW_param2}}. Supports four modes matching that estimator:
#' \itemize{
#'   \item \strong{CASE A}: estimate \eqn{(\beta,\sigma)} with \code{z_data}
#'   \item \strong{CASE B}: estimate \eqn{(\beta,\sigma)} without \code{z_data}
#'   \item \strong{CASE C}: estimate \eqn{\beta} with \code{z_data} and fixed \code{sigma}
#'   \item \strong{CASE D}: estimate \eqn{\beta} without \code{z_data} and fixed \code{sigma}
#' }
#' It is recommended to use the same \code{w_min} and \code{w_max} as in \code{\link{get_IPW_param2}}.
#'
#' @param beta Numeric vector of estimated \eqn{\beta}.
#' @param y_data Numeric vector.
#' @param w_data Numeric vector of \eqn{W=\min(X,C)}.
#' @param delta_data Logical or \{0,1\} vector: \eqn{\Delta=1\{X \le C\}}.
#' @param z_data Optional numeric vector (default \code{NULL}); if provided, must match length of \code{y_data}.
#' @param sigma Numeric scalar \eqn{\hat\sigma}. Ignored when \code{sigma_fixed} is provided.
#' @param sigma_fixed Optional positive scalar. If supplied, \eqn{\sigma} is treated as known
#'   and the returned variance is for \eqn{\beta} only.
#' @param w_min,w_max Numeric bounds used internally (same role as in \code{\link{get_IPW_param2}}).
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
#' variance_IPW_param2(beta = c(0.5, 2.5, 0.1, 0.1),
#'                     y_data,
#'                     w_data,
#'                     delta_data,
#'                     z_data,
#'                     sigma = 3.5,
#'                     w_min = w_min, w_max = w_max)
#' # CASE B: true = (0, 3, 4)
#' variance_IPW_param2(beta = c(0.5, 2.5),
#'                     y_data,
#'                     w_data,
#'                     delta_data,
#'                     sigma = 3.5,
#'                     w_min = w_min, w_max = w_max)
#' # CASE C: true = (0, 3, 0, 0)
#' variance_IPW_param2(beta = c(0.5, 2.5, 0.1, 0.1),
#'                     y_data,
#'                     w_data,
#'                     delta_data,
#'                     z_data,
#'                     sigma_fixed = 4,
#'                     w_min = w_min, w_max = w_max)
#' # CASE D: true = (0, 3)
#' variance_IPW_param2(beta = c(0.5, 2.5),
#'                     y_data,
#'                     w_data,
#'                     delta_data,
#'                     sigma_fixed = 4,
#'                     w_min = w_min, w_max = w_max)
#'
#' @importFrom MASS ginv
#' @importFrom stats dnorm pnorm optim var
#' @importFrom truncnorm dtruncnorm ptruncnorm
#' @export
variance_IPW_param2 <- function(beta,
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
  if (is.null(sigma_fixed)) {
    if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
      stop("`sigma` must be a positive numeric scalar.", call. = FALSE)
    }
  } else {
    if (!is.numeric(sigma_fixed) || length(sigma_fixed) != 1L || !is.finite(sigma_fixed) || sigma_fixed <= 0) {
      stop("`sigma_fixed` must be a positive numeric scalar.", call. = FALSE)
    }
  }

  ## ---- Shared helper: full-data score with z ------------------------------
  # (used by CASE A; extend similarly for other cases as needed)
  S_beta_f <- function(beta, y, x, z, sigma) {
    res <- y - (beta[1] + beta[2]*x + beta[3]*z + beta[4]*z*x)
    c(res * c(1, x, z, z*x) / sigma^2,
      res^2 / (2*sigma^4) - 1/(2*sigma^2))
  }

  ## ---- Shared helper: MLE for C | (Y,Z) -----------------------------------
  find_alpha2_MLE <- function(y_data, w_data, delta_data, z_data,
                              w_min = 0, w_max = 12) {
    log_llhd <- function(alpha2_star, tau2_star, y, w, delta, z) {
      if (delta == 1) { # log P(C > w | y,z)
        return(log(1 - truncnorm::ptruncnorm(
          w, a = w_min, b = w_max,
          mean = sum(c(1, y, z) * alpha2_star),
          sd = tau2_star
        )))
      } else {          # log f_C(w | y,z)
        return(log(truncnorm::dtruncnorm(
          w, a = w_min, b = w_max,
          mean = sum(c(1, y, z) * alpha2_star),
          sd = tau2_star
        )))
      }
    }
    log_llhd_sum <- function(alpha2tau2) {
      alpha2 <- alpha2tau2[1:3]
      tau2   <- alpha2tau2[4]
      llhd <- 0
      for (i in seq_along(y_data)) {
        llhd <- llhd + log_llhd(alpha2, tau2, y_data[i], w_data[i], delta_data[i], z_data[i])
      }
      -llhd
    }
    init <- c(mean(w_data[delta_data == 0]), 0.1, 0.1, stats::sd(w_data[delta_data == 0]))
    stats::optim(init, log_llhd_sum)$par
  }

  ## ---- Case dispatch -------------------------------------------------------
  z_is_null      <- is.null(z_data)
  sigma_is_fixed <- !is.null(sigma_fixed)

  if (!sigma_is_fixed && !z_is_null) {
    ## ===== CASE A: (beta, sigma) WITH z =====
    len.beta <- length(beta)

    # Fit censoring model and build survival function with 1/n floor
    a2t2 <- find_alpha2_MLE(y_data, w_data, delta_data, z_data, w_min, w_max)
    alpha2_MLE <- a2t2[1:3]
    tau2       <- a2t2[4]
    surv_c3 <- function(t, y, z, alpha2_star, tau2) {
      max(1 - truncnorm::ptruncnorm(
        t, a = w_min, b = w_max,
        mean = sum(c(1, y, z) * alpha2_star),
        sd = tau2
      ), 1/n)
    }

    # Build IPW score matrix Sbeta_IPW and reference Sbeta (unweighted)
    Sbeta_IPW <- matrix(0, nrow = n, ncol = len.beta + 1L)
    Sbeta_ref <- matrix(0, nrow = n, ncol = len.beta + 1L)
    for (i in seq_len(n)) {
      if (delta_data[i]) {
        wt <- surv_c3(w_data[i], y_data[i], z_data[i], alpha2_MLE, tau2)
        Sbeta_IPW[i, ] <- S_beta_f(beta, y_data[i], w_data[i], z_data[i], sigma) / wt
        Sbeta_ref[i, ] <- S_beta_f(beta, y_data[i], w_data[i], z_data[i], sigma)
      }
    }

    Ahat  <- -stats::cov(Sbeta_IPW, Sbeta_ref)
    Ainv  <- MASS::ginv(Ahat)
    Bhat  <-  stats::cov(Sbeta_IPW)
    return(Ainv %*% Bhat %*% t(Ainv) / n)
  } else if (!sigma_is_fixed && z_is_null) {
    ## ===== CASE B: (beta, sigma) WITHOUT z =====
    len.beta <- length(beta)  # typically 2

    # No-Z score
    S_beta_f_nz <- function(beta, y, x, sigma) {
      res <- y - (beta[1] + beta[2]*x)
      c(res * c(1, x) / sigma^2,
        res^2 / (2*sigma^4) - 1/(2*sigma^2))
    }

    # MLE for C|Y (no Z)
    find_alpha2_MLE_nz <- function(y_data, w_data, delta_data,
                                   w_min = 0, w_max = 12) {
      log_llhd <- function(alpha2_star, tau2_star, y, w, delta) {
        mu <- alpha2_star[1] + alpha2_star[2]*y
        if (delta) {
          log(1 - truncnorm::ptruncnorm(w, a = w_min, b = w_max, mean = mu, sd = tau2_star))
        } else {
          log(truncnorm::dtruncnorm(w, a = w_min, b = w_max, mean = mu, sd = tau2_star))
        }
      }
      obj <- function(par) {
        a2 <- par[1:2]; t2 <- par[3]
        -sum(mapply(log_llhd,
                    MoreArgs = list(alpha2_star = a2, tau2_star = t2),
                    y = y_data, w = w_data, delta = delta_data))
      }
      init <- c(mean(w_data[delta_data == 0]), 0.1, stats::sd(w_data[delta_data == 0]))
      stats::optim(init, obj)$par
    }

    a2t2 <- find_alpha2_MLE_nz(y_data, w_data, delta_data, w_min, w_max)
    alpha2_MLE <- a2t2[1:2]
    tau2       <- a2t2[3]

    surv_c2 <- function(t, y, alpha2_star, tau2) {
      mu <- alpha2_star[1] + alpha2_star[2]*y
      max(1 - truncnorm::ptruncnorm(t, a = w_min, b = w_max, mean = mu, sd = tau2), 1/n)
    }

    Sbeta_IPW <- matrix(0, nrow = n, ncol = len.beta + 1L)
    Sbeta_ref <- matrix(0, nrow = n, ncol = len.beta + 1L)
    for (i in seq_len(n)) {
      if (delta_data[i]) {
        wt <- surv_c2(w_data[i], y_data[i], alpha2_MLE, tau2)
        Sbeta_IPW[i, ] <- S_beta_f_nz(beta, y = y_data[i], x = w_data[i], sigma = sigma) / wt
        Sbeta_ref[i, ] <- S_beta_f_nz(beta, y = y_data[i], x = w_data[i], sigma = sigma)
      }
    }

    Ahat <- -stats::cov(Sbeta_IPW, Sbeta_ref)
    Ainv <- MASS::ginv(Ahat)
    Bhat <-  stats::cov(Sbeta_IPW)
    return(Ainv %*% Bhat %*% t(Ainv) / n)
  } else if (sigma_is_fixed && !z_is_null) {
    ## ===== CASE C: beta ONLY WITH z, sigma fixed =====
    len.beta  <- length(beta)
    sigma_use <- as.numeric(sigma_fixed)

    # censoring model (with Z) reused from CASE A helper in your file:
    a2t2 <- find_alpha2_MLE(y_data, w_data, delta_data, z_data, w_min, w_max)
    alpha2_MLE <- a2t2[1:3]
    tau2       <- a2t2[4]
    surv_c3 <- function(t, y, z, alpha2_star, tau2) {
      max(1 - truncnorm::ptruncnorm(
        t, a = w_min, b = w_max,
        mean = sum(c(1, y, z) * alpha2_star),
        sd = tau2
      ), 1/n)
    }

    # build full (β,σ) score matrices and take β–β block at the end
    Sbeta_IPW_full <- matrix(0, nrow = n, ncol = len.beta + 1L)
    Sbeta_ref_full <- matrix(0, nrow = n, ncol = len.beta + 1L)
    for (i in seq_len(n)) {
      if (delta_data[i]) {
        wt <- surv_c3(w_data[i], y_data[i], z_data[i], alpha2_MLE, tau2)
        s_full <- S_beta_f(beta, y = y_data[i], x = w_data[i], z = z_data[i], sigma = sigma_use)
        Sbeta_IPW_full[i, ] <- s_full / wt
        Sbeta_ref_full[i, ] <- s_full
      }
    }

    A_full <- -stats::cov(Sbeta_IPW_full, Sbeta_ref_full)
    Ainv   <- MASS::ginv(A_full)
    B_full <-  stats::cov(Sbeta_IPW_full)
    V_full <- Ainv %*% B_full %*% t(Ainv) / n
    return(V_full[seq_len(len.beta), seq_len(len.beta), drop = FALSE])
  } else {
    ## ===== CASE D: beta ONLY WITHOUT z, sigma fixed =====
    len.beta  <- length(beta) # typically 2
    sigma_use <- as.numeric(sigma_fixed)

    # No-Z score
    S_beta_f_nz <- function(beta, y, x, sigma) {
      res <- y - (beta[1] + beta[2]*x)
      c(res * c(1, x) / sigma^2,
        res^2 / (2*sigma^4) - 1/(2*sigma^2))
    }

    # MLE for C|Y (no Z)
    find_alpha2_MLE_nz <- function(y_data, w_data, delta_data,
                                   w_min = 0, w_max = 12) {
      log_llhd <- function(alpha2_star, tau2_star, y, w, delta) {
        mu <- alpha2_star[1] + alpha2_star[2]*y
        if (delta) {
          log(1 - truncnorm::ptruncnorm(w, a = w_min, b = w_max, mean = mu, sd = tau2_star))
        } else {
          log(truncnorm::dtruncnorm(w, a = w_min, b = w_max, mean = mu, sd = tau2_star))
        }
      }
      obj <- function(par) {
        a2 <- par[1:2]; t2 <- par[3]
        -sum(mapply(log_llhd,
                    MoreArgs = list(alpha2_star = a2, tau2_star = t2),
                    y = y_data, w = w_data, delta = delta_data))
      }
      init <- c(mean(w_data[delta_data == 0]), 0.1, stats::sd(w_data[delta_data == 0]))
      stats::optim(init, obj)$par
    }

    a2t2 <- find_alpha2_MLE_nz(y_data, w_data, delta_data, w_min, w_max)
    alpha2_MLE <- a2t2[1:2]
    tau2       <- a2t2[3]

    surv_c2 <- function(t, y, alpha2_star, tau2) {
      mu <- alpha2_star[1] + alpha2_star[2]*y
      max(1 - truncnorm::ptruncnorm(t, a = w_min, b = w_max, mean = mu, sd = tau2), 1/n)
    }

    Sbeta_IPW_full <- matrix(0, nrow = n, ncol = len.beta + 1L)
    Sbeta_ref_full <- matrix(0, nrow = n, ncol = len.beta + 1L)
    for (i in seq_len(n)) {
      if (delta_data[i]) {
        wt <- surv_c2(w_data[i], y_data[i], alpha2_MLE, tau2)
        s_full <- S_beta_f_nz(beta, y = y_data[i], x = w_data[i], sigma = sigma_use)
        Sbeta_IPW_full[i, ] <- s_full / wt
        Sbeta_ref_full[i, ] <- s_full
      }
    }

    A_full <- -stats::cov(Sbeta_IPW_full, Sbeta_ref_full)
    Ainv   <- MASS::ginv(A_full)
    B_full <-  stats::cov(Sbeta_IPW_full)
    V_full <- Ainv %*% B_full %*% t(Ainv) / n
    return(V_full[seq_len(len.beta), seq_len(len.beta), drop = FALSE])
  }
}
