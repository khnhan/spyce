#' Compute the conditional Kaplan-Meier function to estimate \eqn{S_{X|Z_c,Z_d}(t|z_c,z_d)}
#'
#' The function `conditional_Kaplan_Meier` finds the conditional Kaplan-Meier estimator of conditional survival function \eqn{S_{X|Z_c,Z_d}(t|z_c, z_d)}, where \eqn{Z_c} is a continuous random variable, and \eqn{Z_d} is a discrete random variable. If `z_c` and `z_c_data` or `z_d` and `z_d_data` are not specified,
#'
#' To find \eqn{S_{C|Z_c,Z_d}(t|z_c,z_d)}, apply `delta_data = 1-delta_data` beforehand to switch the role of \eqn{X} and \eqn{C}.
#'
#' @param t value of \eqn{t} at which \eqn{S_{X|Z_c,Z_d}(t|z_c,z_d)} is computed.
#' @param w_data vector of \eqn{W=min(X,C)} with length \eqn{n}.
#' @param delta_data vector of \eqn{\Delta=I(X\le C)} with length \eqn{n}.
#' @param z_c value of a continuous random variable \eqn{Z_c} at which \eqn{S_{X|Z_c,Z_d}(t|z_c,z_d)} is computed. Default is `NULL`.
#' @param z_c_data vector (or matrix) of \eqn{Z_c} with length (or row length) \eqn{n}. Default is `NULL`.
#' @param z_d value of a discrete random variable \eqn{Z_d} at which \eqn{S_{X|Z_c,Z_d}(t|z_c,z_d)} is computed. Default is `NULL`.
#' @param z_d_data vector (or matrix) of \eqn{Z_d} with length (or row length) \eqn{n}. \eqn{Z_d} is a discrete random variable. Default is `NULL`.
#' @param kern kernel function applied to \eqn{Z_c}. In `conditional_Kaplan_Meier` function, `kern(z_c, z_c_data)` should return a length \eqn{n} vector. See `kern_gaussian` in the example below.
#' @return value of \eqn{S_{X|Z_c,Z_d}(t|z_c,z_d)}.
#' @examples
#' set.seed(11)
#' x_data = rexp(500)
#' c_data = rexp(500, rate = 2)
#' w_data = pmin(x_data, c_data)
#' delta_data = (x_data <= c_data)
#' z_c_data = rnorm(500)
#' z_d_data = rbinom(500, size = 1, prob = 0.3)
#' kern_gaussian = function(cond, cond_data, h = sd(cond_data)){ #Gaussian kernel
#'   exp(- 0.5 * apply(as.matrix(cond_data), 1, function(x) sum((cond - x)^2)) / h^2) / (sqrt(2*pi)*h)
#' }
#' kern_gaussian(z_c = 1, z_c_data = z_c_data)
#'
#' conditional_Kaplan_Meier(t = 1, w_data, delta_data,
#'                          kern = kern_gaussian) #compare with exp(-1) = 0.368
#'
#' conditional_Kaplan_Meier(t = 1, w_data, delta_data,
#'                          z_c = 1, z_c_data = z_c_data,
#'                          kern = kern_gaussian) #compare with exp(-1) = 0.368
#' conditional_Kaplan_Meier(t = 1, w_data, delta_data,
#'                          z_d = 0, z_d_data = z_d_data,
#'                          kern = kern_gaussian) #compare with exp(-1) = 0.368
#' conditional_Kaplan_Meier(t = 1, w_data, delta_data,
#'                          z_c = 1, z_c_data = z_c_data, z_d = 0, z_d_data = z_d_data,
#'                          kern = kern_gaussian) #compare with exp(-1) = 0.368
#'
#' @export
conditional_Kaplan_Meier = function(t, w_data, delta_data,
                                    z_c = NULL, z_c_data = NULL, z_d = NULL, z_d_data = NULL,
                                    kern){
  n = length(w_data)
  rows_match_na <- function(X, pattern) {
    # If X is a plain vector, treat it as one-column matrix
    if (is.null(pattern)){
      return(rep(TRUE, n))
    }
    if (is.atomic(X) && is.null(dim(X))) {
      X <- matrix(X, ncol = 1)
    }
    stopifnot(ncol(X) == length(pattern))
    M <- matrix(rep(pattern, each = nrow(X)), nrow = nrow(X))
    eq <- (X == M) | (is.na(X) & is.na(M))
    (rowSums(eq) == ncol(X))
  }

  idx_j = which((w_data <= t) & (delta_data == 1)) #index of j in z_d_data==z_d
  if(length(idx_j)!=0){
    kernel_vals = ifelse(rep(is.null(z_c), n), rep(1, n), kern(z_c, z_c_data)) * rows_match_na(z_d_data, z_d) #used in denominator
    denom = sapply((w_data)[idx_j], function(w) {
      sum(kernel_vals * (w_data >= w))})
    log_vals = log(1 - kernel_vals[idx_j] / denom)
    log_vals = ifelse(is.nan(log_vals), 0, log_vals)
    exp(sum(log_vals))
  } else{
    1
  }
}
