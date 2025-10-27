#' Plug-in variance for SPYCE of (beta, sigma) with nonparametrically estimated \eqn{f_{X\mid Z}} and \eqn{f_{C\mid Y,Z}}
#'
#' @description
#' Computes the plug-in variance–covariance matrix for the nonparametric SPYCE estimator
#' \code{\link{get_SPYCE_nonparam12}}. Supports four modes matching that estimator:
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
#' @param z_data Optional numeric vector (binary \{0,1\}); if provided, must match length of \code{y_data}.
#' @param surv_c_vals Numeric vector (length \code{length(y_data)}): \eqn{S_{C\mid Y,Z}(W_i\mid Y_i,Z_i)} (or \eqn{S_{C\mid Y}} if \code{z_data} is \code{NULL}).
#' @param surv_x_vals Numeric vector (length \code{length(y_data)}): \eqn{S_{X\mid Y,Z}(W_i\mid Y_i,Z_i)} (or \eqn{S_{X\mid Y}} if \code{z_data} is \code{NULL}).
#' @param sigma Numeric scalar \eqn{\hat\sigma}. Ignored when \code{sigma_fixed} is provided.
#' @param sigma_fixed Optional positive scalar; if supplied, \eqn{\sigma} is treated as known and the returned variance is for \eqn{\beta} only.
#' @param tt Number of Gaussian–quadrature grid points; defaults to 20.
#' @param m Number of grid points on \code{[min(w_data), max(w_data)]} to solve the integral equation; defaults to 20.
#' @param h3 Positive bandwidth used for kernel weighting for \eqn{Y} variable.
#'
#' @returns
#' A variance–covariance matrix:
#' \itemize{
#'   \item CASE A/B: \code{(length(beta)+1) x (length(beta)+1)}
#'   \item CASE C/D: \code{length(beta) x length(beta)}
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
#' surv_x_vals = sapply(1:n, function(i){
#'   conditional_Kaplan_Meier(w_data[i], w_data, delta_data,
#'                            z_c = y_data[i], z_c_data = y_data,
#'                            z_d = z_data[i], z_d_data = z_data,
#'                            kern = kern_gaussian)})
#' surv_c_vals_nz = sapply(1:n, function(i){
#'   conditional_Kaplan_Meier(w_data[i], w_data, 1-delta_data,
#'                            z_c = y_data[i], z_c_data = y_data,
#'                            kern = kern_gaussian)})
#' surv_x_vals_nz = sapply(1:n, function(i){
#'   conditional_Kaplan_Meier(w_data[i], w_data, delta_data,
#'                            z_c = y_data[i], z_c_data = y_data,
#'                            kern = kern_gaussian)})
#' # CASE A: true = (0, 3, 0, 0, 4)
#' variance_SPYCE_nonparam12(beta = c(0.5, 2.5, 0.1, 0.1),
#'                           y_data,
#'                           w_data,
#'                           delta_data,
#'                           z_data,
#'                           surv_c_vals = surv_c_vals,
#'                           surv_x_vals = surv_x_vals,
#'                           sigma = 3.5,
#'                           h3 = 3*sd(y_data))
#' # CASE B: true = (0, 3, 4)
#' variance_SPYCE_nonparam12(beta = c(0.5, 2.5),
#'                           y_data,
#'                           w_data,
#'                           delta_data,
#'                           surv_c_vals = surv_c_vals_nz,
#'                           surv_x_vals = surv_x_vals_nz,
#'                           sigma = 3.5,
#'                           h3 = 3*sd(y_data))
#' # CASE C: true = (0, 3, 0, 0)
#' variance_SPYCE_nonparam12(beta = c(0.5, 2.5, 0.1, 0.1),
#'                           y_data,
#'                           w_data,
#'                           delta_data,
#'                           z_data,
#'                           surv_c_vals = surv_c_vals,
#'                           surv_x_vals = surv_x_vals,
#'                           sigma_fixed = 4,
#'                           h3 = 3*sd(y_data))
#' # CASE D: true = (0, 3)
#' variance_SPYCE_nonparam12(beta = c(0.5, 2.5),
#'                           y_data,
#'                           w_data,
#'                           delta_data,
#'                           surv_c_vals = surv_c_vals_nz,
#'                           surv_x_vals = surv_x_vals_nz,
#'                           sigma_fixed = 4,
#'                           h3 = 3*sd(y_data))
#' @importFrom MASS ginv
#' @importFrom stats var dnorm approx
#' @export
variance_SPYCE_nonparam12 <- function(beta,
                                      y_data,
                                      w_data,
                                      delta_data,
                                      z_data = NULL,
                                      surv_c_vals,
                                      surv_x_vals,
                                      sigma,
                                      sigma_fixed = NULL,
                                      tt = 20, m = 20,
                                      h3) {

  ## ---- Basic checks ----
  stopifnot(is.numeric(beta), is.numeric(y_data), is.numeric(w_data))
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
    stop("`surv_c_vals` must be numeric and the same length as `y_data`.", call. = FALSE)
  if (!is.numeric(surv_x_vals) || length(surv_x_vals) != n)
    stop("`surv_x_vals` must be numeric and the same length as `y_data`.", call. = FALSE)
  if (is.null(sigma_fixed)) {
    if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma <= 0)
      stop("`sigma` must be a positive numeric scalar.", call. = FALSE)
  } else {
    if (!is.numeric(sigma_fixed) || length(sigma_fixed) != 1L || !is.finite(sigma_fixed) || sigma_fixed <= 0)
      stop("`sigma_fixed` must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(h3) || length(h3) != 1L || h3 <= 0)
    stop("`h3` must be a positive scalar bandwidth.", call. = FALSE)

  z_is_null      <- is.null(z_data)
  sigma_is_fixed <- !is.null(sigma_fixed)

  ## ---- Shared helpers (match get_SPYCE_nonparam12) ------------------------
  S_beta_f <- function(beta, y, x, z, sigma){
    res <- y - (beta[1] + beta[2]*x + beta[3]*z + beta[4]*z*x)
    c(res * c(1, x, z, z*x) / sigma^2,
      res^2 / (2*sigma^4) - 1/(2*sigma^2))
  }
  S_beta_f <- Vectorize(S_beta_f, vectorize.args = "x")

  S_beta_f_nz <- function(beta, y, x, sigma){
    res <- y - (beta[1] + beta[2]*x)
    c(res * c(1, x) / sigma^2,
      res^2 / (2*sigma^4) - 1/(2*sigma^2))
  }
  S_beta_f_nz <- Vectorize(S_beta_f_nz, vectorize.args = "x")

  kern <- function(cond, cond_data, h){
    exp(-0.5 * ((cond - cond_data)^2) / h^2) / (sqrt(2*pi)*h)
  }
  gauss <- function(tt, len = 3){
    grid <- seq(-len, len, length.out = tt)
    d    <- stats::dnorm(grid)
    list(x = grid, w = d/sum(d))
  }
  find_weight <- function(x, x_grid){
    m <- length(x_grid)
    w <- rep(0, m)
    if (x <= x_grid[1]) {
      w[1] <- 1
    } else if (x >= x_grid[m]) {
      w[m] <- 1
    } else {
      idx <- which.max(x < x_grid)
      w[idx]   <- (x - x_grid[idx-1]) / (x_grid[idx] - x_grid[idx-1])
      w[idx-1] <- (x_grid[idx] - x)   / (x_grid[idx] - x_grid[idx-1])
    }
    w
  }

  ## ---- With-Z blocks ------------------------------------------------------
  S_beta_kern1 <- function(beta, y, w, delta, z,
                           y_data, w_data, delta_data, z_data,
                           surv_c_vals,
                           sigma){
    if (!delta) {
      idx <- which((delta_data == 1) & (w_data > w) & (z_data == z))
      if (length(idx) != 0L) {
        dens  <- stats::dnorm(y,
                              mean = cbind(1, w_data[idx], z, z*w_data[idx]) %*% beta,
                              sd   = sigma)
        wt    <- dens / surv_c_vals[idx]
        denom <- sum(wt)
        if (denom == 0) return(S_beta_f(beta, y, w, z, sigma))
        num   <- sapply(idx, function(ii) S_beta_f(beta, y, w_data[ii], z, sigma)) %*% wt
        return(as.numeric(num / denom))
      } else {
        return(S_beta_f(beta, y, w, z, sigma))
      }
    } else {
      return(S_beta_f(beta, y, w, z, sigma))
    }
  }

  b_xz_gauss_kern12 <- function(beta, x_a, z,
                                y_data, w_data, delta_data, z_data,
                                surv_c_vals, surv_x_vals,
                                h3,
                                sigma,
                                tt = 20){
    n <- length(y_data)
    sbeta_kern1_vals <- matrix(0, nrow = n, ncol = length(beta)+1L)
    for (idx in which((delta_data == 0) & (z_data == z))) {
      sbeta_kern1_vals[idx, ] <- S_beta_kern1(beta, y_data[idx], w_data[idx], 0, z_data[idx],
                                              y_data, w_data, delta_data, z_data,
                                              surv_c_vals, sigma)
    }
    cc <- gauss(tt)

    temp2 <- function(x, y){
      idx0 <- which((delta_data == 0) & (z_data == z))
      if (length(idx0) != 0L) {
        kv    <- kern(y, y_data[idx0], h3)
        denom <- sum(kv / surv_x_vals[idx0])
        if (denom == 0) return(rep(0, length(beta)+1L))
        idx1 <- which((delta_data == 0) & (x <= w_data) & (z_data == z))
        num1 <- if (length(idx1) != 0L) {
          S_beta_f(beta, y, x, z, sigma) *
            sum(kv[y_data[idx0] %in% y_data[idx1]] / surv_x_vals[idx0][y_data[idx0] %in% y_data[idx1]])
        } else rep(0, length(beta)+1L)
        idx2 <- which((delta_data == 0) & (x > w_data) & (z_data == z))
        num2 <- if (length(idx2) != 0L) {
          sapply(idx2, function(ii) sbeta_kern1_vals[ii, ]) %*% (kern(y, y_data[idx2], h3) / surv_x_vals[idx2])
        } else rep(0, length(beta)+1L)
        as.numeric((num1 + num2) / denom)
      } else {
        S_beta_f(beta, y, x, z, sigma)
      }
    }

    temp3 <- function(x){
      yc <- sigma * cc$x + sum(c(1, x, z, z*x) * beta)
      sapply(yc, function(yy) temp2(x, yy)) %*% cc$w
    }
    temp3 <- Vectorize(temp3, vectorize.args = "x")
    temp3(x_a)
  }

  L_xz_gauss_kern12 <- function(beta, x_a, z,
                                y_data, w_data, delta_data, z_data,
                                surv_c_vals, surv_x_vals, h3,
                                sigma,
                                tt = 20){
    cc <- gauss(tt)

    temp <- function(x_a, c, y){
      idx1 <- which((delta_data == 1) & (z_data == z) & (w_data > c))
      if (length(idx1) != 0L) {
        dens  <- stats::dnorm(y,
                              mean = cbind(1, w_data[idx1], z, z*w_data[idx1]) %*% beta,
                              sd   = sigma)
        Wmat  <- sapply(idx1, function(ii) find_weight(w_data[ii], x_a))
        wt    <- dens / surv_c_vals[idx1]
        denom <- sum(wt)
        if (denom == 0) return(rep(0, length(x_a)))
        as.numeric(Wmat %*% wt / denom)
      } else {
        find_weight(c, x_a)
      }
    }

    temp2 <- function(x_a, x, y){
      idx0 <- which((delta_data == 0) & (z_data == z))
      if (length(idx0) != 0L) {
        kv <- kern(y, y_data[idx0], h3)
        denom <- sum(kv / surv_x_vals[idx0])
        if (denom == 0) return(rep(0, length(x_a)))
        idx2 <- which(x > w_data[idx0])
        if (length(idx2) != 0L) {
          num2 <- sapply(idx2, function(k) temp(x_a, w_data[idx0][k], y)) %*%
            (kv[idx2] / surv_x_vals[idx0][idx2])
        } else num2 <- rep(0, length(x_a))
        as.numeric(num2 / denom)
      } else {
        rep(0, length(x_a))
      }
    }

    temp3 <- function(x_a, x){
      yc <- sigma * cc$x + sum(c(1, x, z, z*x) * beta)
      sapply(yc, function(yy) temp2(x_a, x, yy)) %*% cc$w
    }
    temp3 <- Vectorize(temp3, vectorize.args = "x")

    temp4 <- function(x, y){
      idx0 <- which((delta_data == 0) & (z_data == z))
      if (length(idx0) != 0L) {
        kv <- kern(y, y_data[idx0], h3)
        denom <- sum(kv / surv_x_vals[idx0])
        if (denom == 0) return(0)
        idx1 <- which(x <= w_data[idx0])
        if (length(idx1) != 0L) {
          num1 <- sum(kv[idx1] / surv_x_vals[idx0][idx1])
          return(num1 / denom)
        } else return(0)
      } else return(0)
    }
    temp4 <- Vectorize(temp4, vectorize.args = "y")

    temp5 <- function(x){
      yc <- sigma * cc$x + sum(c(1, x, z, z*x) * beta)
      sum(temp4(x, yc) * cc$w)
    }
    temp5 <- Vectorize(temp5, vectorize.args = "x")

    diag(temp5(x_a)) + t(temp3(x_a, x_a))
  }

  a_gauss_kern12 <- function(beta, x_a, z,
                             y_data, w_data, delta_data, z_data,
                             surv_c_vals, surv_x_vals, h3,
                             sigma,
                             tt = 20){
    MASS::ginv(L_xz_gauss_kern12(beta, x_a, z,
                                 y_data, w_data, delta_data, z_data,
                                 surv_c_vals, surv_x_vals, h3,
                                 sigma, tt)) %*%
      t(b_xz_gauss_kern12(beta, x_a, z,
                          y_data, w_data, delta_data, z_data,
                          surv_c_vals, surv_x_vals, h3,
                          sigma, tt))
  }

  S_eff_gauss_kern12 <- function(beta, y, w, delta, z,
                                 x_a, a0,
                                 y_data, w_data, delta_data, z_data,
                                 surv_c_vals, surv_x_vals,
                                 sigma,
                                 tt = 20){
    if (!delta) {
      temp <- function(x_a, c, y){
        idx1 <- which((delta_data == 1) & (z_data == z) & (w_data > c))
        if (length(idx1) != 0L) {
          dens <- stats::dnorm(y,
                               mean = cbind(1, w_data[idx1], z, z*w_data[idx1]) %*% beta,
                               sd   = sigma)
          Wmat <- sapply(idx1, function(ii) find_weight(w_data[ii], x_a))
          wt   <- dens / surv_c_vals[idx1]
          denom <- sum(wt)
          if (denom == 0) return(rep(0, length(x_a)))
          as.numeric(Wmat %*% wt / denom)
        } else {
          find_weight(c, x_a)
        }
      }
      sbeta <- S_beta_kern1(beta, y, w, 0, z,
                            y_data, w_data, delta_data, z_data,
                            surv_c_vals, sigma)
      return(sbeta - t(a0) %*% temp(x_a, w, y))
    } else {
      a0_w <- numeric(length(beta) + 1L)
      for (j in seq_len(length(beta) + 1L)) {
        a0_w[j] <- stats::approx(x_a, a0[, j], w, rule = 2)$y
      }
      sbeta <- S_beta_f(beta, y, w, z, sigma)
      return(sbeta - a0_w)
    }
  }

  ## ---- No-Z blocks --------------------------------------------------------
  S_beta_kern1_nz <- function(beta, y, w, delta,
                              y_data, w_data, delta_data,
                              surv_c_vals,
                              sigma){
    if (!delta) {
      idx <- which((delta_data == 1) & (w_data > w))
      if (length(idx) != 0L) {
        dens  <- stats::dnorm(y,
                              mean = cbind(1, w_data[idx]) %*% beta,
                              sd   = sigma)
        wt    <- dens / surv_c_vals[idx]
        denom <- sum(wt)
        if (denom == 0) return(S_beta_f_nz(beta, y, w, sigma))
        num   <- sapply(idx, function(ii) S_beta_f_nz(beta, y, w_data[ii], sigma)) %*% wt
        return(as.numeric(num / denom))
      } else {
        return(S_beta_f_nz(beta, y, w, sigma))
      }
    } else {
      return(S_beta_f_nz(beta, y, w, sigma))
    }
  }

  b_x_gauss_kern12_nz <- function(beta, x_a,
                                  y_data, w_data, delta_data,
                                  surv_c_vals, surv_x_vals,
                                  h3,
                                  sigma,
                                  tt = 20){
    n <- length(y_data)
    sbeta_kern1_vals <- matrix(0, nrow = n, ncol = length(beta)+1L)
    for (idx in which(delta_data == 0)) {
      sbeta_kern1_vals[idx, ] <- S_beta_kern1_nz(beta, y_data[idx], w_data[idx], 0,
                                                 y_data, w_data, delta_data,
                                                 surv_c_vals, sigma)
    }
    cc <- gauss(tt)

    temp2 <- function(x, y){
      idx0 <- which(delta_data == 0)
      if (length(idx0) != 0L) {
        kv    <- kern(y, y_data[idx0], h3)
        denom <- sum(kv / surv_x_vals[idx0])
        if (denom == 0) return(rep(0, length(beta)+1L))
        idx1 <- which((delta_data == 0) & (x <= w_data))
        num1 <- if (length(idx1) != 0L) {
          S_beta_f_nz(beta, y, x, sigma) *
            sum(kv[y_data[idx0] %in% y_data[idx1]] / surv_x_vals[idx0][y_data[idx0] %in% y_data[idx1]])
        } else rep(0, length(beta)+1L)
        idx2 <- which((delta_data == 0) & (x > w_data))
        num2 <- if (length(idx2) != 0L) {
          sapply(idx2, function(ii) sbeta_kern1_vals[ii, ]) %*% (kern(y, y_data[idx2], h3) / surv_x_vals[idx2])
        } else rep(0, length(beta)+1L)
        as.numeric((num1 + num2) / denom)
      } else {
        S_beta_f_nz(beta, y, x, sigma)
      }
    }

    temp3 <- function(x){
      yc <- sigma * cc$x + sum(c(1, x) * beta)
      sapply(yc, function(yy) temp2(x, yy)) %*% cc$w
    }
    temp3 <- Vectorize(temp3, vectorize.args = "x")
    temp3(x_a)
  }

  L_x_gauss_kern12_nz <- function(beta, x_a,
                                  y_data, w_data, delta_data,
                                  surv_c_vals, surv_x_vals, h3,
                                  sigma,
                                  tt = 20){
    cc <- gauss(tt)

    temp <- function(x_a, c, y){
      idx1 <- which((delta_data == 1) & (w_data > c))
      if (length(idx1) != 0L) {
        dens <- stats::dnorm(y,
                             mean = cbind(1, w_data[idx1]) %*% beta,
                             sd   = sigma)
        Wmat <- sapply(idx1, function(ii) find_weight(w_data[ii], x_a))
        wt   <- dens / surv_c_vals[idx1]
        denom <- sum(wt)
        if (denom == 0) return(rep(0, length(x_a)))
        as.numeric(Wmat %*% wt / denom)
      } else {
        find_weight(c, x_a)
      }
    }

    temp2 <- function(x_a, x, y){
      idx0 <- which(delta_data == 0)
      if (length(idx0) != 0L) {
        kv <- kern(y, y_data[idx0], h3)
        denom <- sum(kv / surv_x_vals[idx0])
        if (denom == 0) return(rep(0, length(x_a)))
        idx2 <- which(x > w_data[idx0])
        if (length(idx2) != 0L) {
          num2 <- sapply(idx2, function(k) temp(x_a, w_data[idx0][k], y)) %*%
            (kv[idx2] / surv_x_vals[idx0][idx2])
        } else num2 <- rep(0, length(x_a))
        as.numeric(num2 / denom)
      } else {
        rep(0, length(x_a))
      }
    }

    temp3 <- function(x_a, x){
      yc <- sigma * cc$x + sum(c(1, x) * beta)
      sapply(yc, function(yy) temp2(x_a, x, yy)) %*% cc$w
    }
    temp3 <- Vectorize(temp3, vectorize.args = "x")

    temp4 <- function(x, y){
      idx0 <- which(delta_data == 0)
      if (length(idx0) != 0L) {
        kv <- kern(y, y_data[idx0], h3)
        denom <- sum(kv / surv_x_vals[idx0])
        if (denom == 0) return(0)
        idx1 <- which(x <= w_data[idx0])
        if (length(idx1) != 0L) {
          num1 <- sum(kv[idx1] / surv_x_vals[idx0][idx1])
          return(num1 / denom)
        } else return(0)
      } else return(0)
    }
    temp4 <- Vectorize(temp4, vectorize.args = "y")

    temp5 <- function(x){
      yc <- sigma * cc$x + sum(c(1, x) * beta)
      sum(temp4(x, yc) * cc$w)
    }
    temp5 <- Vectorize(temp5, vectorize.args = "x")

    diag(temp5(x_a)) + t(temp3(x_a, x_a))
  }

  a_gauss_kern12_nz <- function(beta, x_a,
                                y_data, w_data, delta_data,
                                surv_c_vals, surv_x_vals, h3,
                                sigma,
                                tt = 20){
    MASS::ginv(L_x_gauss_kern12_nz(beta, x_a,
                                   y_data, w_data, delta_data,
                                   surv_c_vals, surv_x_vals, h3,
                                   sigma, tt)) %*%
      t(b_x_gauss_kern12_nz(beta, x_a,
                            y_data, w_data, delta_data,
                            surv_c_vals, surv_x_vals, h3,
                            sigma, tt))
  }

  S_eff_gauss_kern12_nz <- function(beta, y, w, delta,
                                    x_a, a0,
                                    y_data, w_data, delta_data,
                                    surv_c_vals, surv_x_vals,
                                    sigma,
                                    tt = 20){
    if (!delta) {
      temp <- function(x_a, c, y){
        idx1 <- which((delta_data == 1) & (w_data > c))
        if (length(idx1) != 0L) {
          dens <- stats::dnorm(y,
                               mean = cbind(1, w_data[idx1]) %*% beta,
                               sd   = sigma)
          Wmat <- sapply(idx1, function(ii) find_weight(w_data[ii], x_a))
          wt   <- dens / surv_c_vals[idx1]
          denom <- sum(wt)
          if (denom == 0) return(rep(0, length(x_a)))
          as.numeric(Wmat %*% wt / denom)
        } else {
          find_weight(c, x_a)
        }
      }
      sbeta <- S_beta_kern1_nz(beta, y, w, 0,
                               y_data, w_data, delta_data,
                               surv_c_vals, sigma)
      return(sbeta - t(a0) %*% temp(x_a, w, y))
    } else {
      a0_w <- numeric(length(beta) + 1L)
      for (j in seq_len(length(beta) + 1L)) {
        a0_w[j] <- stats::approx(x_a, a0[, j], w, rule = 2)$y
      }
      sbeta <- S_beta_f_nz(beta, y, w, sigma)
      return(sbeta - a0_w)
    }
  }

  ## ---- Build grid for integral equation -----------------------------------
  x_a <- seq(min(w_data, na.rm = TRUE), max(w_data, na.rm = TRUE), length.out = m)

  ## ---- Case dispatch: compute var of S_eff and invert ---------------------
  if (!sigma_is_fixed && !z_is_null) {
    # CASE A: with Z, free sigma
    len.beta <- length(beta)
    sigma_use <- sigma

    a_cache <- list(
      `0` = a_gauss_kern12(beta, x_a, 0, y_data, w_data, delta_data, z_data,
                           surv_c_vals, surv_x_vals, h3, sigma_use, tt),
      `1` = a_gauss_kern12(beta, x_a, 1, y_data, w_data, delta_data, z_data,
                           surv_c_vals, surv_x_vals, h3, sigma_use, tt)
    )

    S_eff <- matrix(0, nrow = n, ncol = len.beta + 1L)
    for (i in seq_len(n)) {
      a0 <- a_cache[[as.character(z_data[i])]]
      S_eff[i, ] <- S_eff_gauss_kern12(beta, y_data[i], w_data[i], delta_data[i], z_data[i],
                                       x_a, a0,
                                       y_data, w_data, delta_data, z_data,
                                       surv_c_vals, surv_x_vals,
                                       sigma_use, tt)
    }
    Sigma <- stats::var(S_eff)
    return(MASS::ginv(Sigma) / n)

  } else if (!sigma_is_fixed && z_is_null) {
    # CASE B: no Z, free sigma
    len.beta <- length(beta)
    sigma_use <- sigma

    a0 <- a_gauss_kern12_nz(beta, x_a,
                            y_data, w_data, delta_data,
                            surv_c_vals, surv_x_vals, h3,
                            sigma_use, tt)

    S_eff <- matrix(0, nrow = n, ncol = len.beta + 1L)
    for (i in seq_len(n)) {
      S_eff[i, ] <- S_eff_gauss_kern12_nz(beta, y_data[i], w_data[i], delta_data[i],
                                          x_a, a0,
                                          y_data, w_data, delta_data,
                                          surv_c_vals, surv_x_vals,
                                          sigma_use, tt)
    }
    Sigma <- stats::var(S_eff)
    return(MASS::ginv(Sigma) / n)

  } else if (sigma_is_fixed && !z_is_null) {
    # CASE C: with Z, sigma fixed -> take β–β block
    len.beta  <- length(beta)
    sigma_use <- as.numeric(sigma_fixed)

    a_cache <- list(
      `0` = a_gauss_kern12(beta, x_a, 0, y_data, w_data, delta_data, z_data,
                           surv_c_vals, surv_x_vals, h3, sigma_use, tt),
      `1` = a_gauss_kern12(beta, x_a, 1, y_data, w_data, delta_data, z_data,
                           surv_c_vals, surv_x_vals, h3, sigma_use, tt)
    )

    S_eff_full <- matrix(0, nrow = n, ncol = len.beta + 1L)
    for (i in seq_len(n)) {
      a0 <- a_cache[[as.character(z_data[i])]]
      S_eff_full[i, ] <- S_eff_gauss_kern12(beta, y_data[i], w_data[i], delta_data[i], z_data[i],
                                            x_a, a0,
                                            y_data, w_data, delta_data, z_data,
                                            surv_c_vals, surv_x_vals,
                                            sigma_use, tt)
    }
    V_full <- MASS::ginv(stats::var(S_eff_full)) / n
    return(V_full[seq_len(len.beta), seq_len(len.beta), drop = FALSE])

  } else {
    # CASE D: no Z, sigma fixed -> β–β block
    len.beta  <- length(beta)
    sigma_use <- as.numeric(sigma_fixed)

    a0 <- a_gauss_kern12_nz(beta, x_a,
                            y_data, w_data, delta_data,
                            surv_c_vals, surv_x_vals, h3,
                            sigma_use, tt)

    S_eff_full <- matrix(0, nrow = n, ncol = len.beta + 1L)
    for (i in seq_len(n)) {
      S_eff_full[i, ] <- S_eff_gauss_kern12_nz(beta, y_data[i], w_data[i], delta_data[i],
                                               x_a, a0,
                                               y_data, w_data, delta_data,
                                               surv_c_vals, surv_x_vals,
                                               sigma_use, tt)
    }
    V_full <- MASS::ginv(stats::var(S_eff_full)) / n
    return(V_full[seq_len(len.beta), seq_len(len.beta), drop = FALSE])
  }
}
