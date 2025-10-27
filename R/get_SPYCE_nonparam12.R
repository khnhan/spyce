#' SPYCE for (beta, sigma) with nonparametrically estimated \eqn{f_{X\mid Z}} and \eqn{f_{C\mid Y,Z}}
#'
#' @description
#' Solves the efficient estimating equation
#' \deqn{\sum_{i=1}^n S_{\mathrm{eff}}(Y_i,W_i,\Delta_i,Z_i;\theta)=0}
#' via \code{nleqslv}, where the efficient score uses kernel-based components.
#' This nonparametric version expects survival weights for both \eqn{C\mid Y,Z} (or \eqn{C\mid Y})
#' and \eqn{X\mid Y,Z} (or \eqn{X\mid Y}) supplied as \code{surv_c_vals} and \code{surv_x_vals}.
#'
#' Supported modes:
#' - (A) estimate \eqn{(\beta,\sigma)} with \code{z_data}
#' - (B) estimate \eqn{(\beta,\sigma)} without \code{z_data}
#' - (C) estimate \eqn{\beta} with \code{z_data} and fixed \code{sigma}
#' - (D) estimate \eqn{\beta} without \code{z_data} and fixed \code{sigma}
#'
#' @param y_data Numeric vector.
#' @param w_data Numeric vector of \eqn{W=\min(X,C)}.
#' @param delta_data Logical or \{0,1\} vector: \eqn{\Delta=1\{X \le C\}}.
#' @param z_data Optional numeric vector (binary \{0,1\}); if provided, must match length of \code{y_data}.
#' @param surv_c_vals Numeric vector, length \code{length(y_data)}: conditional survival
#'   \eqn{S_{C\mid Y,Z}(W_i\mid Y_i,Z_i)} (or \eqn{S_{C\mid Y}(W_i\mid Y_i)}) at each observation.
#' @param surv_x_vals Numeric vector, length \code{length(y_data)}: conditional survival
#'   \eqn{S_{X\mid Y,Z}(W_i\mid Y_i,Z_i)} (or \eqn{S_{X\mid Y}(W_i\mid Y_i)}) at each observation.
#' @param init Numeric vector of free parameters:
#'   - Case A: \code{c(beta0,beta1,beta2,beta3, sigma)}
#'   - Case B: \code{c(beta0,beta1, sigma)}
#'   - Case C: \code{c(beta0,beta1,beta2,beta3)}
#'   - Case D: \code{c(beta0,beta1)}
#' @param sigma_fixed Optional positive scalar; if supplied, \eqn{\sigma} is treated as known.
#' @param tt Number of grid points where the integral is computed; defaults to 20.
#' @param m Number of grid points over \code{[min(w_data), max(w_data)]} to solve the Fredholm equation; default 20.
#' @param h3 Positive bandwidth used for kernel weighting for \eqn{Y} variable.
#' @param nleqslv_args Named list of extra args passed to \code{nleqslv::nleqslv()}.
#'
#' @returns list with \code{$beta_SPYCE_nonparam12}, \code{$sigma_SPYCE_nonparam12}, \code{$root}.
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
#' get_SPYCE_nonparam12(y_data,
#'                      w_data,
#'                      delta_data,
#'                      z_data,
#'                      surv_c_vals = surv_c_vals,
#'                      surv_x_vals = surv_x_vals,
#'                      init = c(0.5, 2.5, 0.1, 0.1, 3.5),
#'                      h3 = 3*sd(y_data),
#'                      nleqslv_args = list())
#' # CASE B: true = (0, 3, 4)
#' get_SPYCE_nonparam12(y_data,
#'                      w_data,
#'                      delta_data,
#'                      surv_c_vals = surv_c_vals_nz,
#'                      surv_x_vals = surv_x_vals_nz,
#'                      init = c(0.5, 2.5, 3.5),
#'                      h3 = 3*sd(y_data),
#'                      nleqslv_args = list())
#' # CASE C: true = (0, 3, 0, 0)
#' get_SPYCE_nonparam12(y_data,
#'                      w_data,
#'                      delta_data,
#'                      z_data,
#'                      surv_c_vals = surv_c_vals,
#'                      surv_x_vals = surv_x_vals,
#'                      init = c(0.5, 2.5, 0.1, 0.1),
#'                      sigma_fixed = 4,
#'                      h3 = 3*sd(y_data),
#'                      nleqslv_args = list())
#' # CASE D: true = (0, 3)
#' get_SPYCE_nonparam12(y_data,
#'                      w_data,
#'                      delta_data,
#'                      surv_c_vals = surv_c_vals_nz,
#'                      surv_x_vals = surv_x_vals_nz,
#'                      init = c(0.5, 2.5),
#'                      sigma_fixed = 4,
#'                      h3 = 3*sd(y_data),
#'                      nleqslv_args = list())
#' @importFrom nleqslv nleqslv
#' @importFrom MASS ginv
#' @importFrom stats dnorm approx
#' @export
get_SPYCE_nonparam12 <- function(y_data,
                                 w_data,
                                 delta_data,
                                 z_data = NULL,
                                 surv_c_vals,
                                 surv_x_vals,
                                 init,
                                 sigma_fixed = NULL,
                                 tt = 20, m = 20,
                                 h3,
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
    stop("`surv_c_vals` must be numeric and the same length as `y_data`.", call. = FALSE)
  if (!is.numeric(surv_x_vals) || length(surv_x_vals) != n)
    stop("`surv_x_vals` must be numeric and the same length as `y_data`.", call. = FALSE)
  if (!is.null(sigma_fixed)) {
    if (!is.numeric(sigma_fixed) || length(sigma_fixed) != 1L || !is.finite(sigma_fixed) || sigma_fixed <= 0)
      stop("`sigma_fixed` must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(h3) || length(h3) != 1L || h3 <= 0)
    stop("`h3` must be a positive scalar bandwidth.", call. = FALSE)

  z_is_null      <- is.null(z_data)
  sigma_is_fixed <- !is.null(sigma_fixed)

  ## ---- Core helpers (shared) ----------------------------------------------

  # Full-data scores
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

  # Kernel over Y (Gaussian)
  kern <- function(cond, cond_data, h){
    exp(-0.5 * ((cond - cond_data)^2) / h^2) / (sqrt(2*pi)*h)
  }

  gauss <- function(tt, len = 3){
    grid <- seq(-len, len, length.out = tt)
    d    <- stats::dnorm(grid)
    list(x = grid, w = d / sum(d))
  }

  # Linear weights to project a(C,·) onto grid x_a
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

  ## ---- With-Z building blocks ---------------------------------------------

  # Imputation score using kernel weights (with Z)
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
    Lmat <- L_xz_gauss_kern12(beta, x_a, z,
                              y_data, w_data, delta_data, z_data,
                              surv_c_vals, surv_x_vals, h3,
                              sigma, tt)
    bmat <- b_xz_gauss_kern12(beta, x_a, z,
                              y_data, w_data, delta_data, z_data,
                              surv_c_vals, surv_x_vals, h3,
                              sigma, tt)
    MASS::ginv(Lmat) %*% t(bmat)
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

  ## ---- No-Z building blocks -----------------------------------------------

  # Imputation score using kernel weights (NO Z)
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
    Lmat <- L_x_gauss_kern12_nz(beta, x_a,
                                y_data, w_data, delta_data,
                                surv_c_vals, surv_x_vals, h3,
                                sigma, tt)
    bmat <- b_x_gauss_kern12_nz(beta, x_a,
                                y_data, w_data, delta_data,
                                surv_c_vals, surv_x_vals, h3,
                                sigma, tt)
    MASS::ginv(Lmat) %*% t(bmat)
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

  ## ---- Estimating functions for the four cases ----------------------------

  # A: with Z, free sigma
  pe_gauss_kern12_A <- function(betasigma, y_data, w_data, delta_data, z_data,
                                x_a, h3, surv_c_vals, surv_x_vals,
                                tt = 20){
    beta  <- betasigma[1:4]
    sigma <- betasigma[5]
    a_cache <- list(
      `0` = a_gauss_kern12(beta, x_a, 0, y_data, w_data, delta_data, z_data,
                           surv_c_vals, surv_x_vals, h3, sigma, tt),
      `1` = a_gauss_kern12(beta, x_a, 1, y_data, w_data, delta_data, z_data,
                           surv_c_vals, surv_x_vals, h3, sigma, tt)
    )
    val <- rep(0, 5)
    for (i in seq_len(length(y_data))) {
      a0  <- a_cache[[as.character(z_data[i])]]
      val <- val + S_eff_gauss_kern12(beta, y_data[i], w_data[i], delta_data[i], z_data[i],
                                      x_a, a0, y_data, w_data, delta_data, z_data,
                                      surv_c_vals, surv_x_vals, sigma, tt)
    }
    val
  }

  # B: no Z, free sigma
  pe_gauss_kern12_B <- function(theta, y_data, w_data, delta_data,
                                x_a, h3, surv_c_vals, surv_x_vals,
                                tt = 20){
    if (length(theta) != 3L)
      stop("Case B expects init = c(beta0,beta1,sigma).", call. = FALSE)
    beta  <- theta[1:2]
    sigma <- theta[3]
    a0 <- a_gauss_kern12_nz(beta, x_a, y_data, w_data, delta_data,
                            surv_c_vals, surv_x_vals, h3, sigma, tt)
    val <- rep(0, 3L)
    for (i in seq_len(length(y_data))) {
      val <- val + S_eff_gauss_kern12_nz(beta, y_data[i], w_data[i], delta_data[i],
                                         x_a, a0, y_data, w_data, delta_data,
                                         surv_c_vals, surv_x_vals, sigma, tt)
    }
    val
  }

  # C: with Z, sigma fixed -> return beta components only
  pe_gauss_kern12_C <- function(beta, y_data, w_data, delta_data, z_data,
                                x_a, h3, surv_c_vals, surv_x_vals,
                                sigma_fixed, tt = 20){
    sigma <- as.numeric(sigma_fixed)
    a_cache <- list(
      `0` = a_gauss_kern12(beta, x_a, 0, y_data, w_data, delta_data, z_data,
                           surv_c_vals, surv_x_vals, h3, sigma, tt),
      `1` = a_gauss_kern12(beta, x_a, 1, y_data, w_data, delta_data, z_data,
                           surv_c_vals, surv_x_vals, h3, sigma, tt)
    )
    val_full <- rep(0, length(beta) + 1L)  # build full then drop last
    for (i in seq_len(length(y_data))) {
      a0 <- a_cache[[as.character(z_data[i])]]
      s  <- S_eff_gauss_kern12(beta, y_data[i], w_data[i], delta_data[i], z_data[i],
                               x_a, a0, y_data, w_data, delta_data, z_data,
                               surv_c_vals, surv_x_vals, sigma, tt)
      val_full <- val_full + s
    }
    val_full[seq_len(length(beta))]
  }

  # D: no Z, sigma fixed -> return beta components only
  pe_gauss_kern12_D <- function(beta, y_data, w_data, delta_data,
                                x_a, h3, surv_c_vals, surv_x_vals,
                                sigma_fixed, tt = 20){
    sigma <- as.numeric(sigma_fixed)
    a0 <- a_gauss_kern12_nz(beta, x_a, y_data, w_data, delta_data,
                            surv_c_vals, surv_x_vals, h3, sigma, tt)
    val_full <- rep(0, length(beta) + 1L)
    for (i in seq_len(length(y_data))) {
      s <- S_eff_gauss_kern12_nz(beta, y_data[i], w_data[i], delta_data[i],
                                 x_a, a0, y_data, w_data, delta_data,
                                 surv_c_vals, surv_x_vals, sigma, tt)
      val_full <- val_full + s
    }
    val_full[seq_len(length(beta))]
  }

  ## ---- Root solve (dispatch by case) --------------------------------------
  x_a <- seq(min(w_data, na.rm = TRUE), max(w_data, na.rm = TRUE), length.out = m)

  if (!sigma_is_fixed && !z_is_null) {
    # CASE A
    if (!is.numeric(init) || length(init) != 5L)
      stop("Case A expects init = c(beta0,beta1,beta2,beta3,sigma).", call. = FALSE)
    nleqslv_call <- c(list(
      x  = init,
      fn = pe_gauss_kern12_A,
      y_data = y_data,
      w_data = w_data,
      delta_data = delta_data,
      z_data = z_data,
      x_a = x_a,
      h3 = h3,
      surv_c_vals = surv_c_vals,
      surv_x_vals = surv_x_vals,
      tt = tt
    ), nleqslv_args)

  } else if (!sigma_is_fixed && z_is_null) {
    # CASE B
    if (!is.numeric(init) || length(init) != 3L)
      stop("Case B expects init = c(beta0,beta1,sigma).", call. = FALSE)
    nleqslv_call <- c(list(
      x  = init,
      fn = pe_gauss_kern12_B,
      y_data = y_data,
      w_data = w_data,
      delta_data = delta_data,
      x_a = x_a,
      h3 = h3,
      surv_c_vals = surv_c_vals,
      surv_x_vals = surv_x_vals,
      tt = tt
    ), nleqslv_args)

  } else if (sigma_is_fixed && !z_is_null) {
    # CASE C
    if (!is.numeric(init) || length(init) != 4L)
      stop("Case C expects init = c(beta0,beta1,beta2,beta3).", call. = FALSE)
    nleqslv_call <- c(list(
      x  = init,
      fn = pe_gauss_kern12_C,
      y_data = y_data,
      w_data = w_data,
      delta_data = delta_data,
      z_data = z_data,
      x_a = x_a,
      h3 = h3,
      surv_c_vals = surv_c_vals,
      surv_x_vals = surv_x_vals,
      sigma_fixed = sigma_fixed,
      tt = tt
    ), nleqslv_args)

  } else {
    # CASE D
    if (!is.numeric(init) || length(init) != 2L)
      stop("Case D expects init = c(beta0,beta1).", call. = FALSE)
    nleqslv_call <- c(list(
      x  = init,
      fn = pe_gauss_kern12_D,
      y_data = y_data,
      w_data = w_data,
      delta_data = delta_data,
      x_a = x_a,
      h3 = h3,
      surv_c_vals = surv_c_vals,
      surv_x_vals = surv_x_vals,
      sigma_fixed = sigma_fixed,
      tt = tt
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
  if (!sigma_is_fixed && !z_is_null) {
    # A
    sigma_SPYCE_nonparam12 <- as.numeric(tail(theta_hat, 1L))
    beta_SPYCE_nonparam12  <- as.numeric(head(theta_hat, -1L))
  } else if (!sigma_is_fixed && z_is_null) {
    # B
    sigma_SPYCE_nonparam12 <- as.numeric(tail(theta_hat, 1L))
    beta_SPYCE_nonparam12  <- as.numeric(head(theta_hat, -1L))
  } else if (sigma_is_fixed && !z_is_null) {
    # C
    sigma_SPYCE_nonparam12 <- as.numeric(sigma_fixed)
    beta_SPYCE_nonparam12  <- as.numeric(theta_hat)
  } else {
    # D
    sigma_SPYCE_nonparam12 <- as.numeric(sigma_fixed)
    beta_SPYCE_nonparam12  <- as.numeric(theta_hat)
  }

  out <- list(
    beta_SPYCE_nonparam12  = beta_SPYCE_nonparam12,
    sigma_SPYCE_nonparam12 = sigma_SPYCE_nonparam12,
    root = root
  )
  class(out) <- c("spyce_nonparam12", class(out))
  out
}
