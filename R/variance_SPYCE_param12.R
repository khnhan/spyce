#' Plug-in variance for SPYCE of (beta, sigma) with parametrically estimated \eqn{f_{X\mid Z}} and \eqn{f_{C\mid Y,Z}}
#'
#' @description
#' Computes the plug-in variance–covariance matrix for the SPYCE estimator
#' \code{\link{get_SPYCE_param12}}. Supports four modes matching that estimator:
#' \itemize{
#'   \item \strong{CASE A}: estimate \eqn{(\beta,\sigma)} with \code{z_data}
#'   \item \strong{CASE B}: estimate \eqn{(\beta,\sigma)} without \code{z_data}
#'   \item \strong{CASE C}: estimate \eqn{\beta} with \code{z_data} and fixed \code{sigma}
#'   \item \strong{CASE D}: estimate \eqn{\beta} without \code{z_data} and fixed \code{sigma}
#' }
#' It is recommended to use the same \code{w_min}, \code{w_max}, \code{tt}, and \code{m} as in \code{\link{get_SPYCE_param12}}.
#'
#' @param beta Numeric vector of estimated \eqn{\beta}.
#' @param y_data Numeric vector.
#' @param w_data Numeric vector of \eqn{W=\min(X,C)}.
#' @param delta_data Logical or \{0,1\} vector: \eqn{\Delta=1\{X \le C\}}.
#' @param z_data Optional numeric vector (default \code{NULL}); if provided, must match length of \code{y_data}.
#' @param sigma Numeric scalar \eqn{\hat\sigma}. Ignored when \code{sigma_fixed} is provided.
#' @param sigma_fixed Optional positive scalar; if supplied, \eqn{\sigma} is treated as known and the returned variance is for \eqn{\beta} only.
#' @param w_min,w_max Numeric bounds used internally (same role as in \code{\link{get_SPYCE_param12}}).
#' @param tt Number of grid points where the integral is computed; defaults to 20.
#' @param m Number of grid points between \code{w_min} and \code{w_max} for which the integral equation is solved; defaults to 20.
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
#' # Set w_min and w_max
#' w_min = min(w_data) - sd(w_data)
#' w_max = max(w_data) + sd(w_data)
#' # CASE A: true = (0, 3, 0, 0, 4)
#' variance_SPYCE_param12(beta = c(0.5, 2.5, 0.1, 0.1),
#'                        y_data,
#'                        w_data,
#'                        delta_data,
#'                        z_data,
#'                        sigma = 3.5,
#'                        w_min = w_min, w_max = w_max)
#' # CASE B: true = (0, 3, 4)
#' variance_SPYCE_param12(beta = c(0.5, 2.5),
#'                        y_data,
#'                        w_data,
#'                        delta_data,
#'                        sigma = 3.5,
#'                        w_min = w_min, w_max = w_max)
#' # CASE C: true = (0, 3, 0, 0)
#' variance_SPYCE_param12(beta = c(0.5, 2.5, 0.1, 0.1),
#'                        y_data,
#'                        w_data,
#'                        delta_data,
#'                        z_data,
#'                        sigma_fixed = 4,
#'                        w_min = w_min, w_max = w_max)
#' # CASE D: true = (0, 3)
#' variance_SPYCE_param12(beta = c(0.5, 2.5),
#'                        y_data,
#'                        w_data,
#'                        delta_data,
#'                        sigma_fixed = 4,
#'                        w_min = w_min, w_max = w_max)
#' @importFrom MASS ginv
#' @importFrom stats var cov dnorm pnorm optim approx
#' @importFrom truncnorm dtruncnorm ptruncnorm
#' @export
variance_SPYCE_param12 <- function(beta,
                                   y_data,
                                   w_data,
                                   delta_data,
                                   z_data = NULL,
                                   sigma,
                                   sigma_fixed = NULL,
                                   w_min, w_max,
                                   tt = 20, m = 20) {

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
  if (is.null(sigma_fixed)) {
    if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma <= 0)
      stop("`sigma` must be a positive numeric scalar.", call. = FALSE)
  } else {
    if (!is.numeric(sigma_fixed) || length(sigma_fixed) != 1L || !is.finite(sigma_fixed) || sigma_fixed <= 0)
      stop("`sigma_fixed` must be a positive numeric scalar.", call. = FALSE)
  }

  ## ---- Core building blocks (with Z) ----
  S_beta_f <- function(beta, y, x, z, sigma) {
    res <- y - (beta[1] + beta[2]*x + beta[3]*z + beta[4]*z*x)
    c(res * c(1, x, z, z*x) / sigma^2,
      res^2 / (2*sigma^4) - 1/(2*sigma^2))
  }
  S_beta_f <- Vectorize(S_beta_f, vectorize.args = "x")

  S_beta <- function(beta, y, w, delta, z, alpha1_star,
                     sigma, tau1, w_min, w_max) {
    v_x   <- 1/((beta[2] + beta[4]*z)^2 / sigma^2 + 1/tau1^2)
    eta_x <- v_x * ((beta[2] + beta[4]*z) * (y - (beta[1] + beta[3]*z)) / sigma^2 +
                      (alpha1_star[1] + alpha1_star[2]*z) / tau1^2)
    if (!delta) {
      if (w > w_max) {
        S_beta_f(beta, y, w, z, sigma)
      } else {
        x_norm <- (seq(w, w_max, length.out = 20) - eta_x)/sqrt(v_x)
        d      <- stats::dnorm(x_norm)
        num    <- S_beta_f(beta, y, x_norm*sqrt(v_x)+eta_x, z, sigma) %*% d
        denom  <- sum(d)
        if (denom == 0) return(S_beta_f(beta, y, w, z, sigma))
        as.numeric(num/denom)
      }
    } else {
      S_beta_f(beta, y, w, z, sigma)
    }
  }

  gauss <- function(tt, len = 3) {
    grid <- seq(-len, len, length.out = tt)
    d    <- stats::dnorm(grid)
    list(x = grid, w = d / sum(d))
  }

  b_xz_gauss_param12 <- function(beta, x_a, z,
                                 alpha1_star, alpha2_star,
                                 sigma, tau1, tau2,
                                 tt = 20, w_min = 0, w_max = 12) {
    cc  <- gauss(tt)
    len <- 20
    temp2 <- function(x, y) {
      v_x   <- 1/((beta[2] + beta[4]*z)^2 / sigma^2 + 1/tau1^2)
      eta_x <- v_x * ((beta[2] + beta[4]*z) * (y - (beta[1] + beta[3]*z)) / sigma^2 +
                        (alpha1_star[1] + alpha1_star[2]*z) / tau1^2)
      if (x < w_min) {
        S_beta_f(beta, y, x, z, sigma)
      } else {
        c_grid <- seq(w_min, min(x, w_max), length.out = len)
        dens_c <- truncnorm::dtruncnorm(c_grid, a = w_min, b = w_max,
                                        mean = sum(c(1, y, z) * alpha2_star), sd = tau2)
        Sbeta  <- sapply(c_grid, function(w)
          S_beta(beta, y, w, 0, z, alpha1_star, sigma, tau1, w_min, w_max))
        Sbeta[, c(1, len)] <- Sbeta[, c(1, len)]/2
        by <- c_grid[2] - c_grid[1]
        S_beta_f(beta, y, x, z, sigma) *
          (1 - truncnorm::ptruncnorm(x, a = w_min, b = w_max,
                                     mean = sum(c(1, y, z) * alpha2_star), sd = tau2)) +
          Sbeta %*% dens_c * by
      }
    }
    temp3 <- function(x) {
      sapply(sigma * cc$x + sum(c(1, x, z, z*x) * beta),
             function(y_norm) temp2(x, y_norm)) %*% cc$w
    }
    temp3 <- Vectorize(temp3, vectorize.args = "x")
    temp3(x_a)
  }

  L_xz_gauss_param12 <- function(beta, x_a, z,
                                 alpha1_star, alpha2_star,
                                 sigma, tau1, tau2,
                                 tt = 20, w_min = 0, w_max = 12) {
    cc <- gauss(tt)
    temp <- function(x_a, c, y) {
      v_x   <- 1/((beta[2] + beta[4]*z)^2 / sigma^2 + 1/tau1^2)
      eta_x <- v_x * ((beta[2] + beta[4]*z) * (y - (beta[1] + beta[3]*z)) / sigma^2 +
                        (alpha1_star[1] + alpha1_star[2]*z) / tau1^2)
      p   <- truncnorm::dtruncnorm(x_a, a = w_min, b = w_max, mean = eta_x, sd = sqrt(v_x))
      num <- (x_a > c) * p
      denom <- sum(num)
      ifelse(is.nan(num/denom), 0, num/denom)
    }
    temp2 <- function(x_a, x, y) {
      c_grid <- seq(w_min, w_max, length = 20)
      dens   <- truncnorm::dtruncnorm(c_grid, a = w_min, b = w_max,
                                      mean = sum(c(1, y, z) * alpha2_star), sd = tau2)
      dens <- dens/sum(dens)
      sapply(c_grid, function(c_norm) (x > c_norm) * temp(x_a, c_norm, y)) %*% dens
    }
    temp3 <- function(x_a, x) {
      sapply(sigma * cc$x + sum(beta * c(1, x, z, z*x)), function(y_norm) {
        temp2(x_a, x, y_norm)
      }) %*% cc$w
    }
    temp3 <- Vectorize(temp3, vectorize.args = "x")
    temp4 <- function(x, y) {
      c_grid <- seq(w_min, w_max, length = 20)
      dens   <- truncnorm::dtruncnorm(c_grid, a = w_min, b = w_max,
                                      mean = sum(c(1, y, z) * alpha2_star), sd = tau2)
      dens <- dens/sum(dens)
      sum((x <= c_grid) * dens)
    }
    temp4 <- Vectorize(temp4, vectorize.args = "y")
    temp5 <- function(x) {
      sum(temp4(x, sigma * cc$x + sum(beta * c(1, x, z, z*x))) * cc$w)
    }
    temp5 <- Vectorize(temp5, vectorize.args = "x")
    diag(temp5(x_a)) + t(temp3(x_a, x_a))
  }

  a_gauss_param12 <- function(beta, x_a, z,
                              alpha1_star, alpha2_star,
                              sigma, tau1, tau2,
                              tt = 20, w_min = 0, w_max = 12) {
    MASS::ginv(L_xz_gauss_param12(beta, x_a, z, alpha1_star, alpha2_star,
                                  sigma, tau1, tau2, tt, w_min, w_max)) %*%
      t(b_xz_gauss_param12(beta, x_a, z, alpha1_star, alpha2_star,
                           sigma, tau1, tau2, tt, w_min, w_max))
  }

  S_eff_gauss_param12 <- function(beta, y, w, delta, z,
                                  x_a, a0,
                                  alpha1_star,
                                  sigma, tau1, tau2,
                                  tt = 20, w_min = 0, w_max = 12) {
    len_beta <- length(beta) + 1
    if (!delta) {
      temp <- function(x_a, c, y) {
        v_x   <- 1/((beta[2] + beta[4]*z)^2 / sigma^2 + 1/tau1^2)
        eta_x <- v_x * ((beta[2] + beta[4]*z) * (y - (beta[1] + beta[3]*z)) / sigma^2 +
                          (alpha1_star[1] + alpha1_star[2]*z) / tau1^2)
        p   <- truncnorm::dtruncnorm(x_a, a = w_min, b = w_max, mean = eta_x, sd = sqrt(v_x))
        num <- (x_a > c) * p
        denom <- sum(num)
        ifelse(is.nan(num/denom), 0, num/denom)
      }
      sbeta <- S_beta(beta, y, w, 0, z, alpha1_star, sigma, tau1, w_min, w_max)
      sbeta - t(a0) %*% temp(x_a, w, y)
    } else {
      a0_w <- numeric(len_beta)
      for (j in 1:len_beta) {
        a0_w[j] <- stats::approx(x_a, a0[, j], w, rule = 2)$y
      }
      sbeta <- S_beta_f(beta, y, w, z, sigma)
      sbeta - a0_w
    }
  }

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

  find_alpha2_MLE <- function(y_data, w_data, delta_data, z_data,
                              w_min = 0, w_max = 12) {
    log_llhd <- function(alpha2_star, tau2_star, y, w, delta, z) {
      if (delta == 1) {
        log(1 - truncnorm::ptruncnorm(w, a = w_min, b = w_max,
                                      mean = sum(c(1, y, z) * alpha2_star), sd = tau2_star))
      } else {
        log(truncnorm::dtruncnorm(w, a = w_min, b = w_max,
                                  mean = sum(c(1, y, z) * alpha2_star), sd = tau2_star))
      }
    }
    obj <- function(par) {
      a2 <- par[1:3]; t2 <- par[4]
      -sum(mapply(log_llhd,
                  MoreArgs = list(alpha2_star = a2, tau2_star = t2),
                  y = y_data, w = w_data, delta = delta_data, z = z_data))
    }
    init <- c(mean(w_data[delta_data == 0]), 0.1, 0.1, stats::sd(w_data[delta_data == 0]))
    stats::optim(init, obj)$par
  }

  find_alpha1_MLE_nz <- function(beta,y_data,w_data,delta_data,sigma,w_min,w_max){
    logL <- function(alpha1_star,tau1_star,y,w,delta){
      v_x <- 1/((beta[2]^2)/sigma^2+1/tau1_star^2)
      eta_x <- v_x*(beta[2]*(y-beta[1])/sigma^2+alpha1_star/tau1_star^2)
      if (delta)
        log(truncnorm::dtruncnorm(w,a=w_min,b=w_max,mean=eta_x,sd=sqrt(v_x)))
      else
        log(pnorm(w_max,mean=eta_x,sd=sqrt(v_x))-pnorm(w,mean=eta_x,sd=sqrt(v_x)))+
        dnorm(y,mean=beta[1]+beta[2]*alpha1_star,
              sd=sqrt(sigma^2+(beta[2]^2)*tau1_star^2),log=TRUE)
    }
    obj <- function(par){
      a1<-par[1];t1<-par[2]
      -sum(mapply(logL,MoreArgs=list(alpha1_star=a1,tau1_star=t1),
                  y=y_data,w=w_data,delta=delta_data))
    }
    stats::optim(c(mean(w_data[delta_data==1]),sd(w_data[delta_data==1])),obj)$par
  }

  find_alpha2_MLE_nz <- function(y_data,w_data,delta_data,w_min,w_max){
    logL <- function(alpha2_star,tau2_star,y,w,delta){
      if (delta==1)
        log(1 - truncnorm::ptruncnorm(w,a=w_min,b=w_max,
                                      mean=alpha2_star[1]+alpha2_star[2]*y,sd=tau2_star))
      else
        log(truncnorm::dtruncnorm(w,a=w_min,b=w_max,
                                  mean=alpha2_star[1]+alpha2_star[2]*y,sd=tau2_star))
    }
    obj <- function(par){
      a2<-par[1:2];t2<-par[3]
      -sum(mapply(logL,MoreArgs=list(alpha2_star=a2,tau2_star=t2),
                  y=y_data,w=w_data,delta=delta_data))
    }
    stats::optim(c(mean(w_data[delta_data==0]),0.1,sd(w_data[delta_data==0])),obj)$par
  }

  S_beta_f_nz <- function(beta, y, x, sigma) {
    res <- y - (beta[1] + beta[2]*x)
    c(res * c(1, x) / sigma^2,
      res^2 / (2*sigma^4) - 1/(2*sigma^2))
  }
  S_beta_f_nz <- Vectorize(S_beta_f_nz, vectorize.args = "x")

  S_beta_nz <- function(beta, y, w, delta, alpha1_star,
                        sigma, tau1, w_min, w_max) {
    v_x   <- 1/( (beta[2]^2)/sigma^2 + 1/tau1^2 )
    eta_x <- v_x * ( beta[2]*(y - beta[1])/sigma^2 + alpha1_star/tau1^2 )
    if (!delta) {
      if (w > w_max) return(S_beta_f_nz(beta, y, w, sigma))
      x_norm <- (seq(w, w_max, length.out = 20) - eta_x)/sqrt(v_x)
      d      <- dnorm(x_norm)
      num    <- S_beta_f_nz(beta, y, x_norm*sqrt(v_x)+eta_x, sigma) %*% d
      as.numeric(num/sum(d))
    } else S_beta_f_nz(beta, y, w, sigma)
  }

  b_xz_gauss_param12_nz <- function(beta, x_a,
                                    alpha1_star, alpha2_star,
                                    sigma, tau1, tau2,
                                    tt = 20, w_min = 0, w_max = 12) {
    cc  <- gauss(tt)
    len <- 20
    temp2 <- function(x, y) {
      v_x   <- 1/( (beta[2]^2)/sigma^2 + 1/tau1^2 )
      eta_x <- v_x * ( beta[2]*(y - beta[1])/sigma^2 + alpha1_star/tau1^2 )
      if (x < w_min) {
        S_beta_f_nz(beta, y, x, sigma)
      } else {
        c_grid <- seq(w_min, min(x, w_max), length.out = len)
        dens_c <- truncnorm::dtruncnorm(c_grid, a = w_min, b = w_max,
                                        mean = alpha2_star[1] + alpha2_star[2]*y,
                                        sd = tau2)
        Sbeta  <- sapply(c_grid, function(w)
          S_beta_nz(beta, y, w, 0, alpha1_star, sigma, tau1, w_min, w_max))
        Sbeta[, c(1, len)] <- Sbeta[, c(1, len)]/2
        by <- c_grid[2] - c_grid[1]
        S_beta_f_nz(beta, y, x, sigma) *
          (1 - truncnorm::ptruncnorm(x, a = w_min, b = w_max,
                                     mean = alpha2_star[1] + alpha2_star[2]*y, sd = tau2)) +
          Sbeta %*% dens_c * by
      }
    }
    temp3 <- function(x) {
      sapply(sigma * cc$x + (beta[1] + beta[2]*x),
             function(y_norm) temp2(x, y_norm)) %*% cc$w
    }
    temp3 <- Vectorize(temp3, vectorize.args = "x")
    temp3(x_a)
  }

  L_xz_gauss_param12_nz <- function(beta, x_a,
                                    alpha1_star, alpha2_star,
                                    sigma, tau1, tau2,
                                    tt = 20, w_min = 0, w_max = 12) {
    cc <- gauss(tt)
    temp <- function(x_a, c, y) {
      v_x   <- 1/( (beta[2]^2)/sigma^2 + 1/tau1^2 )
      eta_x <- v_x * ( beta[2]*(y - beta[1])/sigma^2 + alpha1_star/tau1^2 )
      p   <- truncnorm::dtruncnorm(x_a, a = w_min, b = w_max,
                                   mean = eta_x, sd = sqrt(v_x))
      num <- (x_a > c) * p
      denom <- sum(num)
      ifelse(is.nan(num/denom), 0, num/denom)
    }
    temp2 <- function(x_a, x, y) {
      c_grid <- seq(w_min, w_max, length = 20)
      dens   <- truncnorm::dtruncnorm(c_grid, a = w_min, b = w_max,
                                      mean = alpha2_star[1] + alpha2_star[2]*y,
                                      sd = tau2)
      dens <- dens/sum(dens)
      sapply(c_grid, function(c_norm) (x > c_norm) * temp(x_a, c_norm, y)) %*% dens
    }
    temp3 <- function(x_a, x) {
      sapply(sigma * cc$x + (beta[1] + beta[2]*x),
             function(y_norm) temp2(x_a, x, y_norm)) %*% cc$w
    }
    temp3 <- Vectorize(temp3, vectorize.args = "x")
    temp4 <- function(x, y) {
      c_grid <- seq(w_min, w_max, length = 20)
      dens   <- truncnorm::dtruncnorm(c_grid, a = w_min, b = w_max,
                                      mean = alpha2_star[1] + alpha2_star[2]*y, sd = tau2)
      dens <- dens/sum(dens)
      sum((x <= c_grid) * dens)
    }
    temp4 <- Vectorize(temp4, vectorize.args = "y")
    temp5 <- function(x) {
      sum(temp4(x, sigma * cc$x + (beta[1] + beta[2]*x)) * cc$w)
    }
    temp5 <- Vectorize(temp5, vectorize.args = "x")
    diag(temp5(x_a)) + t(temp3(x_a, x_a))
  }

  a_gauss_param12_nz <- function(beta, x_a,
                                 alpha1_star, alpha2_star,
                                 sigma, tau1, tau2,
                                 tt = 20, w_min = 0, w_max = 12) {
    MASS::ginv(L_xz_gauss_param12_nz(beta, x_a,
                                     alpha1_star, alpha2_star,
                                     sigma, tau1, tau2,
                                     tt, w_min, w_max)) %*%
      t(b_xz_gauss_param12_nz(beta, x_a,
                              alpha1_star, alpha2_star,
                              sigma, tau1, tau2,
                              tt, w_min, w_max))
  }

  S_eff_gauss_param12_nz <- function(beta, y, w, delta,
                                     x_a, a0,
                                     alpha1_star,
                                     sigma, tau1, tau2,
                                     tt = 20, w_min = 0, w_max = 12) {
    len_beta <- length(beta) + 1
    if (!delta) {
      temp <- function(x_a, c, y) {
        v_x   <- 1/( (beta[2]^2)/sigma^2 + 1/tau1^2 )
        eta_x <- v_x * ( beta[2]*(y - beta[1])/sigma^2 + alpha1_star/tau1^2 )
        p   <- truncnorm::dtruncnorm(x_a, a = w_min, b = w_max,
                                     mean = eta_x, sd = sqrt(v_x))
        num <- (x_a > c) * p
        denom <- sum(num)
        ifelse(is.nan(num/denom), 0, num/denom)
      }
      sbeta <- S_beta_nz(beta, y, w, 0, alpha1_star, sigma, tau1, w_min, w_max)
      sbeta - t(a0) %*% temp(x_a, w, y)
    } else {
      a0_w <- numeric(len_beta)
      for (j in seq_len(len_beta))
        a0_w[j] <- approx(x_a, a0[, j], w, rule = 2)$y
      sbeta <- S_beta_f_nz(beta, y, w, sigma)
      sbeta - a0_w
    }
  }

  ## ---- Case dispatch ----
  z_is_null      <- is.null(z_data)
  sigma_is_fixed <- !is.null(sigma_fixed)

  x_a <- seq(w_min, w_max, length.out = m)

  if (!sigma_is_fixed && !z_is_null) {
    ## ===== CASE A: (beta, sigma) WITH z =====
    len.beta <- length(beta)

    a1t1 <- find_alpha1_MLE(beta, y_data, w_data, delta_data, z_data, sigma, w_min, w_max)
    a1   <- a1t1[1:2]; t1 <- a1t1[3]
    a2t2 <- find_alpha2_MLE(y_data, w_data, delta_data, z_data, w_min, w_max)
    a2   <- a2t2[1:3]; t2 <- a2t2[4]

    # cache a(x,z) by unique z
    z_levels <- sort(unique(z_data))
    a_cache  <- setNames(vector("list", length(z_levels)), as.character(z_levels))
    for (zk in z_levels) {
      a_cache[[as.character(zk)]] <-
        a_gauss_param12(beta, x_a, zk, a1, a2, sigma, t1, t2, tt, w_min, w_max)
    }

    S_eff <- matrix(NA_real_, nrow = n, ncol = len.beta + 1L)
    for (i in seq_len(n)) {
      a0 <- a_cache[[as.character(z_data[i])]]
      S_eff[i, ] <- S_eff_gauss_param12(beta, y_data[i], w_data[i], delta_data[i], z_data[i],
                                        x_a, a0, a1, sigma, t1, t2, tt, w_min, w_max)
    }
    Sigma <- stats::var(S_eff)
    return(MASS::ginv(Sigma) / n)

  } else if (!sigma_is_fixed && z_is_null) {
    ## ===== CASE B: (beta, sigma) WITHOUT z =====
    len.beta  <- length(beta)
    a1t1 <- find_alpha1_MLE_nz(beta, y_data, w_data, delta_data, sigma, w_min, w_max)
    a1   <- a1t1[1]; t1 <- a1t1[2]
    a2t2 <- find_alpha2_MLE_nz(y_data, w_data, delta_data, w_min, w_max)
    a2   <- a2t2[1:2]; t2 <- a2t2[3]

    a0 <- a_gauss_param12_nz(beta, x_a, a1, a2, sigma, t1, t2, tt, w_min, w_max)

    val <- rep(0, 3)
    for (i in seq_along(y_data)) {
      val <- val + S_eff_gauss_param12_nz(beta, y_data[i], w_data[i], delta_data[i],
                                          x_a, a0, a1, sigma, t1, t2, tt, w_min, w_max)
    }

    S_eff <- matrix(NA_real_, nrow = n, ncol = len.beta + 1L)
    for (i in seq_len(n)) {
      S_eff[i, ] <- S_eff_gauss_param12_nz(beta, y_data[i], w_data[i], delta_data[i],
                                           x_a, a0, a1, sigma, t1, t2, tt, w_min, w_max)
    }
    Sigma <- stats::var(S_eff)
    return(MASS::ginv(Sigma) / n)

  } else if (sigma_is_fixed && !z_is_null) {
    ## ===== CASE C: beta ONLY WITH z, sigma fixed =====
    len.beta  <- length(beta)
    sigma_use <- as.numeric(sigma_fixed)

    a1t1 <- find_alpha1_MLE(beta, y_data, w_data, delta_data, z_data, sigma_use, w_min, w_max)
    a1   <- a1t1[1:2]; t1 <- a1t1[3]
    a2t2 <- find_alpha2_MLE(y_data, w_data, delta_data, z_data, w_min, w_max)
    a2   <- a2t2[1:3]; t2 <- a2t2[4]

    z_levels <- sort(unique(z_data))
    a_cache  <- setNames(vector("list", length(z_levels)), as.character(z_levels))
    for (zk in z_levels) {
      a_cache[[as.character(zk)]] <-
        a_gauss_param12(beta, x_a, zk, a1, a2, sigma_use, t1, t2, tt, w_min, w_max)
    }

    # build full (β,σ) scores then extract β–β block
    S_eff_full <- matrix(NA_real_, nrow = n, ncol = len.beta + 1L)
    for (i in seq_len(n)) {
      a0 <- a_cache[[as.character(z_data[i])]]
      S_eff_full[i, ] <- S_eff_gauss_param12(beta, y_data[i], w_data[i], delta_data[i], z_data[i],
                                             x_a, a0, a1, sigma_use, t1, t2, tt, w_min, w_max)
    }
    Sigma_full <- stats::var(S_eff_full)
    V_full <- MASS::ginv(Sigma_full) / n
    return(V_full[seq_len(len.beta), seq_len(len.beta), drop = FALSE])

  } else {
    ## ===== CASE D: beta ONLY WITHOUT z, sigma fixed =====
    len.beta  <- length(beta)
    sigma_use <- as.numeric(sigma_fixed)

    a1t1 <- find_alpha1_MLE_nz(beta, y_data, w_data, delta_data, sigma_use, w_min, w_max)
    a1   <- a1t1[1]; t1 <- a1t1[2]
    a2t2 <- find_alpha2_MLE_nz(y_data, w_data, delta_data, w_min, w_max)
    a2   <- a2t2[1:2]; t2 <- a2t2[3]

    a0 <- a_gauss_param12_nz(beta, x_a, a1, a2, sigma_use, t1, t2, tt, w_min, w_max)

    S_eff_full <- matrix(NA_real_, nrow = n, ncol = len.beta + 1L)
    for (i in seq_len(n)) {
      S_eff_full[i, ] <- S_eff_gauss_param12_nz(beta, y = y_data[i], w = w_data[i],
                                                delta = delta_data[i],
                                                x_a, a0, a1, sigma_use, t1, t2, tt, w_min, w_max)
    }
    Sigma_full <- stats::var(S_eff_full)
    V_full <- MASS::ginv(Sigma_full) / n
    return(V_full[seq_len(len.beta), seq_len(len.beta), drop = FALSE])
  }
}
