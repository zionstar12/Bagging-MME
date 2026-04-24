# ============================================================================
# depCenR: depTest — Test for dependence between T and C
# ============================================================================

#' Test for the presence of dependence between event and censoring times
#'
#' Performs a hypothesis test of \eqn{H_0: \tau = 0} (independence between
#' event time \eqn{T} and censoring time \eqn{C}) versus
#' \eqn{H_1: \tau > 0} (positive dependence), using the Bagging-MME
#' framework.
#'
#' @inheritParams depEst
#' @param alpha Significance level for the test (default 0.05).
#'
#' @return A list of class \code{"depTest"} containing:
#'   \describe{
#'     \item{tau_hat}{Bagging-MME point estimate of Kendall's \eqn{\tau}.}
#'     \item{ci_95}{Bootstrap 95\% confidence interval for \eqn{\tau}.}
#'     \item{reject_H0}{Logical; \code{TRUE} if the CI excludes zero.}
#'     \item{alpha}{Significance level used.}
#'     \item{depEst_result}{The full \code{"depEst"} object (for inspection).}
#'     \item{call}{The matched call.}
#'   }
#'
#' @details
#' The test rejects \eqn{H_0} when the lower bound of the bootstrap
#' \eqn{(1 - \alpha)} confidence interval for \eqn{\tau} exceeds zero.
#' This is a simple inversion of the Bagging-MME confidence interval and
#' controls type I error at level \eqn{\alpha} for the one-sided alternative.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n <- 200
#' obs_time  <- rexp(n, rate = 0.1)
#' event_ind <- rbinom(n, 1, 0.6)
#'
#' test <- depTest(time = obs_time, event = event_ind,
#'                 B = 49, n_cores = 2)
#' print(test)
#' }
depTest <- function(time,
                    event,
                    B = 199,
                    n_sim = 1000,
                    n_cores = max(1L, parallel::detectCores() - 1L),
                    n_ranges = 4,
                    rho_breaks = c(-0.1818182, 0.3529412, 1.333333,
                                   3.714286, 18),
                    marg_bounds = NULL,
                    loss_weights = c(4, 2, 2, 1, 1, 0.5, 0.5),
                    maxit = 100,
                    temperature = 700,
                    alpha = 0.05,
                    seed = NULL) {

  cl_call <- match.call()

  # Run the full Bagging-MME estimation
  est <- depEst(
    time = time, event = event, B = B, n_sim = n_sim,
    n_cores = n_cores, n_ranges = n_ranges,
    rho_breaks = rho_breaks, marg_bounds = marg_bounds,
    loss_weights = loss_weights, maxit = maxit,
    temperature = temperature, seed = seed
  )

  # CI-based test: reject H0 if lower bound of (1 - alpha) CI > 0
  ci <- quantile(est$tau_boot, probs = c(alpha / 2, 1 - alpha / 2))
  reject <- ci[1] > 0

  structure(
    list(
      tau_hat       = est$tau_hat,
      ci_95         = ci,
      reject_H0     = reject,
      alpha         = alpha,
      depEst_result = est,
      call          = cl_call
    ),
    class = "depTest"
  )
}


#' @export
print.depTest <- function(x, ...) {
  cat("Bagging-MME Dependence Test\n")
  cat("---------------------------\n")
  cat(sprintf("H0: tau = 0  vs  H1: tau > 0\n"))
  cat(sprintf("Kendall's tau (est.):  %.4f\n", x$tau_hat))
  cat(sprintf("%.0f%% bootstrap CI:     [%.4f, %.4f]\n",
              (1 - x$alpha) * 100, x$ci_95[1], x$ci_95[2]))
  cat(sprintf("Decision:              %s H0 at alpha = %.2f\n",
              ifelse(x$reject_H0, "REJECT", "Fail to reject"), x$alpha))
  invisible(x)
}
