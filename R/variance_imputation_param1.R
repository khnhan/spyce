#' Plug-in variance for imputation estimator of (beta, sigma) with parametrically estimated \eqn{f_{X\mid Z}}
#'
#' @description
#' Computes the plug-in variance–covariance matrix for the imputation estimator
#' \code{\link{get_imputation_param1}}. Supports four modes matching that estimator:
#' \itemize{
#'   \item \strong{CASE A}: estimate \eqn{(\beta,\sigma)} with \code{z_data}
#'   \item \strong{CASE B}: estimate \eqn{(\beta,\sigma)} without \code{z_data}
#'   \item \strong{CASE C}: estimate \eqn{\beta} with \code{z_data} and fixed \code{sigma}
#'   \item \strong{CASE D}: estimate \eqn{\beta} without \code{z_data} and fixed \code{sigma}
#' }
#' It is recommended to use the same \code{w_min} and \code{w_max} as in \code{\link{get_imputation_param1}}.
#'
#' @param beta Numeric vector of estimated \eqn{\beta}.
#' @param y_data Numeric vector.
#' @param w_data Numeric vector of \eqn{W=\min(X,C)}.
#' @param delta_data Logical or \{0,1\} vector: \eqn{\Delta=1\{X \le C\}}.
#' @param z_data Optional numeric vector (default \code{NULL}); if provided, must match length of \code{y_data}.
#' @param sigma Numeric scalar \eqn{\hat\sigma}. Ignored when \code{sigma_fixed} is provided.
#' @param sigma_fixed Optional positive scalar. If supplied, \eqn{\sigma} is treated as known
#'   and the returned variance is for \eqn{\beta} only.
#' @param w_min,w_max Numeric bounds used internally (same role as in \code{\link{get_imputation_param1}}).
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
#' variance_imputation_param1(beta = c(0.5, 2.5, 0.1, 0.1),
#'                            y_data,
#'                            w_data,
#'                            delta_data,
#'                            z_data,
#'                            sigma = 3.5,
#'                            w_min = w_min, w_max = w_max)
#' # CASE B: true = (0, 3, 4)
#' variance_imputation_param1(beta = c(0.5, 2.5),
#'                            y_data,
#'                            w_data,
#'                            delta_data,
#'                            sigma = 3.5,
#'                            w_min = w_min, w_max = w_max)
#' # CASE C: true = (0, 3, 0, 0)
#' variance_imputation_param1(beta = c(0.5, 2.5, 0.1, 0.1),
#'                            y_data,
#'                            w_data,
#'                            delta_data,
#'                            z_data,
#'                            sigma_fixed = 4,
#'                            w_min = w_min, w_max = w_max)
#' # CASE D: true = (0, 3)
#' variance_imputation_param1(beta = c(0.5, 2.5),
#'                            y_data,
#'                            w_data,
#'                            delta_data,
#'                            sigma_fixed = 4,
#'                            w_min = w_min, w_max = w_max)
#' @importFrom MASS ginv
#' @importFrom stats dnorm pnorm optim var
#' @importFrom truncnorm dtruncnorm ptruncnorm
#' @export
variance_imputation_param1 <- function(beta,
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

  ## ---- Scores: full-data and imputation (WITH z) --------------------------
  S_beta_f <- function(beta, y, x, z, sigma) {
    res <- y - (beta[1] + beta[2]*x + beta[3]*z + beta[4]*z*x)
    c(res * c(1, x, z, z*x) / sigma^2,
      res^2 / (2*sigma^4) - 1/(2*sigma^2))
  }
  S_beta_f <- Vectorize(S_beta_f, vectorize.args = "x")

  S_beta_imp <- function(beta, y, w, delta, z, alpha1_star,
                         sigma, tau1,
                         w_min, w_max) {
    v_x   <- 1/((beta[2] + beta[4]*z)^2 / sigma^2 + 1/tau1^2)
    eta_x <- v_x * ((beta[2] + beta[4]*z) * (y - (beta[1] + beta[3]*z)) / sigma^2 +
                      (alpha1_star[1] + alpha1_star[2]*z) / tau1^2)
    if (!delta) {
      if (w > w_max) {
        S_beta_f(beta, y, w, z, sigma)
      } else {
        x_grid <- seq(w, w_max, length.out = 20)
        d      <- stats::dnorm((x_grid - eta_x) / sqrt(v_x))
        w_imp  <- sum(x_grid * d) / sum(d)
        S_beta_f(beta, y, w_imp, z, sigma)
      }
    } else {
      S_beta_f(beta, y, w, z, sigma)
    }
  }

  ## ---- Scores: full-data and imputation (NO z) ----------------------------
  S_beta_f_nz <- function(beta, y, x, sigma) {
    res <- y - (beta[1] + beta[2]*x)
    c(res * c(1, x) / sigma^2,
      res^2 / (2*sigma^4) - 1/(2*sigma^2))
  }
  S_beta_f_nz <- Vectorize(S_beta_f_nz, vectorize.args = "x")

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

  ## ---- Nuisance: X|Z MLEs (same as estimator) -----------------------------
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
    obj <- function(par) {
      a1 <- par[1:2]; t1 <- par[3]
      -sum(mapply(log_llhd,
                  MoreArgs = list(alpha1_star = a1, tau1_star = t1),
                  y = y_data, w = w_data, delta = delta_data, z = z_data))
    }
    init <- c(mean(w_data[delta_data == 1]), 0.1, stats::sd(w_data[delta_data == 1]))
    stats::optim(init, obj)$par
  }

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

  ## ---- Case dispatch -------------------------------------------------------
  z_is_null      <- is.null(z_data)
  sigma_is_fixed <- !is.null(sigma_fixed)

  if (!sigma_is_fixed && !z_is_null) {
    ## ===== CASE A: (beta, sigma) WITH z =====
    len.beta <- length(beta)
    a1t1 <- find_alpha1_MLE(beta, y_data, w_data, delta_data, z_data, sigma, w_min, w_max)
    a1   <- a1t1[1:2]; t1 <- a1t1[3]

    S_imp <- matrix(NA_real_, nrow = n, ncol = len.beta + 1L)
    S_ref <- matrix(NA_real_, nrow = n, ncol = len.beta + 1L)
    for (i in seq_len(n)) {
      S_imp[i, ] <- S_beta_imp(beta, y_data[i], w_data[i], delta_data[i], z_data[i], a1, sigma, t1, w_min, w_max)
      S_ref[i, ] <- S_beta_f(  beta, y_data[i], w_data[i],                         z_data[i],       sigma)
    }
    Ahat <- -stats::var(cbind(S_imp, S_ref)) [seq_len(len.beta+1), -(seq_len(len.beta+1))]
    Bhat <-  stats::var(S_imp)
    return(MASS::ginv(Ahat) %*% Bhat %*% t(MASS::ginv(Ahat)) / n)

  } else if (!sigma_is_fixed && z_is_null) {
    ## ===== CASE B: (beta, sigma) WITHOUT z =====
    len.beta <- length(beta)  # expected 2
    a1t1 <- find_alpha1_MLE_nz(beta, y_data, w_data, delta_data, sigma, w_min, w_max)
    a1   <- a1t1[1]; t1 <- a1t1[2]

    S_imp <- matrix(NA_real_, nrow = n, ncol = len.beta + 1L)
    S_ref <- matrix(NA_real_, nrow = n, ncol = len.beta + 1L)
    for (i in seq_len(n)) {
      S_imp[i, ] <- S_beta_imp_nz(beta, y_data[i], w_data[i], delta_data[i], a1, sigma, t1, w_min, w_max)
      S_ref[i, ] <- S_beta_f_nz(  beta, y_data[i], w_data[i],                        sigma)
    }
    Ahat <- -stats::var(cbind(S_imp, S_ref)) [seq_len(len.beta+1), -(seq_len(len.beta+1))]
    Bhat <-  stats::var(S_imp)
    return(MASS::ginv(Ahat) %*% Bhat %*% t(MASS::ginv(Ahat)) / n)

  } else if (sigma_is_fixed && !z_is_null) {
    ## ===== CASE C: beta ONLY WITH z, sigma fixed =====
    len.beta  <- length(beta)
    sigma_use <- as.numeric(sigma_fixed)
    a1t1 <- find_alpha1_MLE(beta, y_data, w_data, delta_data, z_data, sigma_use, w_min, w_max)
    a1   <- a1t1[1:2]; t1 <- a1t1[3]

    # build full (β,σ) scores then take β–β block
    S_imp_full <- matrix(NA_real_, nrow = n, ncol = len.beta + 1L)
    S_ref_full <- matrix(NA_real_, nrow = n, ncol = len.beta + 1L)
    for (i in seq_len(n)) {
      s_imp <- S_beta_imp(beta, y_data[i], w_data[i], delta_data[i], z_data[i], a1, sigma_use, t1, w_min, w_max)
      s_ref <- S_beta_f(  beta, y_data[i], w_data[i],                         z_data[i],        sigma_use)
      S_imp_full[i, ] <- s_imp
      S_ref_full[i, ] <- s_ref
    }
    A_full <- -stats::var(cbind(S_imp_full, S_ref_full)) [seq_len(len.beta+1), -(seq_len(len.beta+1))]
    B_full <-  stats::var(S_imp_full)
    V_full <- MASS::ginv(A_full) %*% B_full %*% t(MASS::ginv(A_full)) / n
    return(V_full[seq_len(len.beta), seq_len(len.beta), drop = FALSE])

  } else {
    ## ===== CASE D: beta ONLY WITHOUT z, sigma fixed =====
    len.beta  <- length(beta) # expected 2
    sigma_use <- as.numeric(sigma_fixed)
    a1t1 <- find_alpha1_MLE_nz(beta, y_data, w_data, delta_data, sigma_use, w_min, w_max)
    a1   <- a1t1[1]; t1 <- a1t1[2]

    S_imp_full <- matrix(NA_real_, nrow = n, ncol = len.beta + 1L)
    S_ref_full <- matrix(NA_real_, nrow = n, ncol = len.beta + 1L)
    for (i in seq_len(n)) {
      s_imp <- S_beta_imp_nz(beta, y_data[i], w_data[i], delta_data[i], a1, sigma_use, t1, w_min, w_max)
      s_ref <- S_beta_f_nz(  beta, y_data[i], w_data[i],                      sigma_use)
      S_imp_full[i, ] <- s_imp
      S_ref_full[i, ] <- s_ref
    }
    A_full <- -stats::var(cbind(S_imp_full, S_ref_full)) [seq_len(len.beta+1), -(seq_len(len.beta+1))]
    B_full <-  stats::var(S_imp_full)
    V_full <- MASS::ginv(A_full) %*% B_full %*% t(MASS::ginv(A_full)) / n
    return(V_full[seq_len(len.beta), seq_len(len.beta), drop = FALSE])
  }
}
