# ============================================================================
# depCenR: Internal helper functions
# ============================================================================

#' Direct Clayton copula conditional inverse (vectorised)
#'
#' Given \eqn{u_1, u_2 \sim \text{Uniform}(0,1)} independent, returns a matrix
#' \code{cbind(u, v)} where the pair \eqn{(u, v)} has Clayton(\eqn{\theta})
#' dependence.
#'
#' The conditional CDF of Clayton(\eqn{\theta}) is:
#' \deqn{h(v \mid u) = u^{-(\theta+1)} \bigl(u^{-\theta} + v^{-\theta} - 1\bigr)^{-(1+\theta)/\theta}}
#' Inverting for \eqn{v}:
#' \deqn{v = \bigl(u^{-\theta} (u_2^{-\theta/(1+\theta)} - 1) + 1\bigr)^{-1/\theta}}
#'
#' This is pure vectorised arithmetic with no copula-package dispatch overhead,
#' making it suitable for repeated calls inside an optimisation loop.
#'
#' @param u1 Numeric vector of uniform(0,1) draws.
#' @param u2 Numeric vector of uniform(0,1) draws (same length as \code{u1}).
#' @param theta Non-negative Clayton dependence parameter.
#'
#' @return A two-column matrix \code{cbind(u, v)}.
#'
#' @keywords internal
#' @noRd
clayton_inv_cond <- function(u1, u2, theta) {
  if (theta < 1e-10) return(cbind(u1, u2))
  v <- (u1^(-theta) * (u2^(-theta / (1 + theta)) - 1) + 1)^(-1 / theta)
  cbind(u1, v)
}

#' Fast sample skewness
#'
#' Computes the (biased) sample skewness without the overhead of
#' \code{moments::skewness()}.
#'
#' @param x Numeric vector.
#'
#' @return Scalar skewness estimate.
#'
#' @keywords internal
#' @noRd
skew_fast <- function(x) {
  n <- length(x)
  if (n < 3L) return(0)
  m <- mean(x)
  s <- sd(x)
  if (!is.finite(s) || s <= 0) return(0)
  mean(((x - m) / s)^3)
}

#' Inverse generalised-gamma survival transform
#'
#' Maps uniform(0,1) quantiles to generalised-gamma random variates using
#' the AFT (log-linear) parameterisation of \pkg{flexsurv}: location
#' \eqn{\mu}, scale \eqn{\sigma}, and shape \eqn{Q}.
#'
#' @param u Numeric vector of uniform(0,1) quantiles.
#' @param mu  Location parameter (log-scale).
#' @param sigma Scale parameter (> 0).
#' @param q Shape parameter (Q != 0).
#'
#' @return Numeric vector of survival times.
#'
#' @keywords internal
#' @noRd
inv_gg_surv <- function(u, mu, sigma, q) {
  exp(sigma / q * (2 * log(q) + log(qgamma(1 - u, shape = 1 / q^2))) + mu)
}

#' Construct the Bagging-MME objective function
#'
#' Returns a closure that evaluates the weighted sum-of-squared-errors between
#' target summary statistics and those produced by a candidate parameter vector.
#'
#' @param n_sim  Number of Monte Carlo draws per evaluation.
#' @param targets Named list with elements \code{pct_T}, \code{mu_T},
#'   \code{mu_C}, \code{sd_T}, \code{sd_C}, \code{q_T}, \code{q_C}.
#' @param weights Numeric vector of length 7 (loss weights; will be normalised).
#'
#' @return A function \code{f(params)} suitable for \code{GenSA::GenSA()}.
#'
#' @keywords internal
#' @noRd
make_objective <- function(n_sim, targets, weights) {
  w <- weights / sum(weights)
  function(params) {
    mu_T    <- max(0, params[1])
    mu_C    <- max(0, params[2])
    sigma_T <- max(0, params[3])
    sigma_C <- max(0, params[4])
    q_T     <- params[5]
    q_C     <- params[6]
    rho     <- max(0, params[7])

    # Correlated uniforms via direct Clayton inverse
    uv <- clayton_inv_cond(runif(n_sim), runif(n_sim), rho)

    # Inverse generalised-gamma survival transform
    X1 <- inv_gg_surv(uv[, 1], mu_T, sigma_T, q_T)
    X2 <- inv_gg_surv(uv[, 2], mu_C, sigma_C, q_C)

    delta <- as.integer(X1 <= X2)

    sq_errs <- c(
      (targets$pct_T - mean(delta))^2,
      (targets$mu_T  - mean(X1[delta == 1]))^2,
      (targets$mu_C  - mean(X2[delta == 0]))^2,
      (targets$sd_T  - sd(X1[delta == 1]))^2,
      (targets$sd_C  - sd(X2[delta == 0]))^2,
      (targets$q_T   - skew_fast(X1[delta == 1]))^2,
      (targets$q_C   - skew_fast(X2[delta == 0]))^2
    )
    sum(w * sq_errs)
  }
}
