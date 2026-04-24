# ============================================================================
# depCenR: depEst — Bagging-MME dependence estimation
# ============================================================================

#' Estimate dependence between event and censoring times via Bagging-MME
#'
#' Applies the Bootstrap-Aggregated Method-of-Moments Estimation (Bagging-MME)
#' algorithm to estimate Kendall's \eqn{\tau} between a survival time \eqn{T}
#' and a censoring time \eqn{C} from right-censored data assuming bivariate
#' generalised-gamma marginals linked by a Clayton copula.
#'
#' @param time Numeric vector of observed times \eqn{\tilde{T} = \min(T, C)}.
#' @param event Integer or logical vector of event indicators (1 = event,
#'   0 = censored).
#' @param B Number of bootstrap resamples (default 99). The original sample is
#'   always included as the \eqn{(B+1)}-th replicate.
#' @param n_sim Number of Monte Carlo draws per GenSA evaluation
#'   (default 1000).
#' @param n_cores Number of parallel workers (default:
#'   \code{parallel::detectCores() - 1}).
#' @param n_ranges Number of Clayton-parameter sub-ranges to search over
#'   (default 4).
#' @param rho_breaks Numeric vector of length \code{n_ranges + 1} giving the
#'   boundaries of the Clayton-parameter sub-ranges.
#' @param marg_bounds A 6 × 2 matrix of lower/upper bounds for the six
#'   marginal parameters (\eqn{\mu_T, \mu_C, \sigma_T, \sigma_C, Q_T, Q_C}).
#'   If \code{NULL} (default), bounds are estimated from the data using
#'   univariate generalised-gamma fits.
#' @param loss_weights Numeric vector of length 7 giving the relative weights
#'   for each summary statistic in the objective function (event proportion,
#'   conditional means, SDs, and skewness).
#' @param maxit Maximum GenSA iterations per sub-range (default 100).
#' @param temperature GenSA temperature parameter (default 700).
#' @param seed Optional random seed for reproducibility.
#'
#' @return A list of class \code{"depEst"} containing:
#'   \describe{
#'     \item{tau_hat}{Point estimate of Kendall's \eqn{\tau} (from the
#'       original sample).}
#'     \item{tau_boot}{Numeric vector of length \eqn{B+1} with all bootstrap
#'       \eqn{\tau} estimates.}
#'     \item{ci_95}{Two-element vector with the bootstrap 2.5\% and 97.5\%
#'       percentiles.}
#'     \item{rho_hat}{Estimated Clayton parameter (original sample).}
#'     \item{best_range}{Index of the majority-vote winning sub-range.}
#'     \item{range_counts}{Integer vector of sub-range selection frequencies.}
#'     \item{call}{The matched call.}
#'   }
#'
#' @details
#' The algorithm proceeds as follows:
#' \enumerate{
#'   \item Draw \eqn{B} bootstrap samples; append the original data.
#'   \item For each replicate, compute seven target summary statistics:
#'     event proportion, conditional means, SDs, and skewness of the event and
#'     censored subgroups.
#'   \item For each of \code{n_ranges} pre-defined Clayton-parameter
#'     sub-ranges, run \code{GenSA::GenSA()} to minimise the weighted
#'     sum-of-squared-errors between simulated and target statistics.
#'   \item Select the best sub-range for each replicate (lowest objective
#'     value) and determine the overall winning sub-range by majority vote.
#'   \item Extract the Clayton parameter from the winning sub-range for all
#'     replicates; convert to Kendall's \eqn{\tau = \rho / (\rho + 2)}.
#'   \item Report the point estimate, bootstrap distribution, and 95\% CI.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulated example
#' set.seed(42)
#' n <- 200
#' T_event   <- rexp(n, rate = 0.1)
#' T_censor  <- rexp(n, rate = 0.15)
#' obs_time  <- pmin(T_event, T_censor)
#' event_ind <- as.integer(T_event <= T_censor)
#'
#' result <- depEst(time = obs_time, event = event_ind,
#'                  B = 49, n_cores = 2)
#' print(result)
#' }
depEst <- function(time,
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
                   seed = NULL) {

  cl_call <- match.call()

  # ---- Input validation ----
  stopifnot(
    is.numeric(time), length(time) > 10,
    is.numeric(event) || is.logical(event),
    length(time) == length(event),
    B >= 1, n_sim >= 100, n_ranges >= 1,
    length(rho_breaks) == n_ranges + 1
  )
  event <- as.integer(event)
  n <- length(time)

  if (!is.null(seed)) set.seed(seed)

  # ---- Reconstruct marginal observations ----
  T1 <- time                            # proxy for event time
  T2 <- time                            # proxy for censoring time
  d1 <- event
  d2 <- 1L - event

  # ---- Marginal parameter bounds (placeholder: user-supplied or fitted) ----
  if (is.null(marg_bounds)) {
    message("Note: marg_bounds not supplied; using default placeholder bounds. ",
            "For production use, supply data-driven bounds from univariate ",
            "generalised-gamma fits (e.g. via flexsurv::flexsurvreg).")
    marg_bounds <- matrix(
      c(4.119, 3.918, 0.276, 0.180, -0.121, -0.181,
        4.359, 4.662, 0.574, 0.359,  1.472,  0.869),
      nrow = 6, ncol = 2,
      dimnames = list(
        c("mu1", "mu2", "sigma1", "sigma2", "q1", "q2"),
        c("lower", "upper")
      )
    )
  }

  rho_lo <- rho_breaks[1:n_ranges]
  rho_hi <- rho_breaks[2:(n_ranges + 1)]

  # ---- Bootstrap samples ----
  boot_list <- replicate(B, {
    idx <- sample.int(n, size = n, replace = TRUE)
    data.frame(T1 = T1[idx], T2 = T2[idx], d1 = d1[idx], d2 = d2[idx])
  }, simplify = FALSE)
  boot_list[[B + 1]] <- data.frame(T1 = T1, T2 = T2, d1 = d1, d2 = d2)

  # ---- Pre-compute target statistics for each bootstrap sample ----
  targets_list <- lapply(boot_list, function(df) {
    list(
      pct_T = sum(df$d1) / nrow(df),
      mu_T  = mean(df$T1[df$d1 == 1]),
      mu_C  = mean(df$T2[df$d2 == 1]),
      sd_T  = sd(df$T1[df$d1 == 1]),
      sd_C  = sd(df$T2[df$d2 == 1]),
      q_T   = skew_fast(df$T1[df$d1 == 1]),
      q_C   = skew_fast(df$T2[df$d2 == 1])
    )
  })

  # ---- Parallel computation ----
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  num_params <- 7
  pkgs <- c("GenSA")

  simul_res <- foreach::foreach(
    num = seq_len(B + 1),
    .combine = "rbind",
    .packages = pkgs
  ) %dopar% {

    # Re-define helpers inside each worker
    clayton_inv_cond_w <- function(u1, u2, theta) {
      if (theta < 1e-10) return(cbind(u1, u2))
      v <- (u1^(-theta) * (u2^(-theta / (1 + theta)) - 1) + 1)^(-1 / theta)
      cbind(u1, v)
    }
    skew_fast_w <- function(x) {
      nn <- length(x)
      if (nn < 3L) return(0)
      m <- mean(x); s <- sd(x)
      if (!is.finite(s) || s <= 0) return(0)
      mean(((x - m) / s)^3)
    }
    inv_gg_surv_w <- function(u, mu, sigma, q) {
      exp(sigma / q * (2 * log(q) + log(qgamma(1 - u, shape = 1 / q^2))) + mu)
    }

    tgt <- targets_list[[num]]
    w   <- loss_weights / sum(loss_weights)

    # Objective function
    obj_fn <- function(params) {
      mu_T    <- max(0, params[1])
      mu_C    <- max(0, params[2])
      sigma_T <- max(0, params[3])
      sigma_C <- max(0, params[4])
      q_T     <- params[5]
      q_C     <- params[6]
      rho     <- max(0, params[7])

      uv <- clayton_inv_cond_w(runif(n_sim), runif(n_sim), rho)
      X1 <- inv_gg_surv_w(uv[, 1], mu_T, sigma_T, q_T)
      X2 <- inv_gg_surv_w(uv[, 2], mu_C, sigma_C, q_C)
      delta <- as.integer(X1 <= X2)

      sq_errs <- c(
        (tgt$pct_T - mean(delta))^2,
        (tgt$mu_T  - mean(X1[delta == 1]))^2,
        (tgt$mu_C  - mean(X2[delta == 0]))^2,
        (tgt$sd_T  - sd(X1[delta == 1]))^2,
        (tgt$sd_C  - sd(X2[delta == 0]))^2,
        (tgt$q_T   - skew_fast_w(X1[delta == 1]))^2,
        (tgt$q_C   - skew_fast_w(X2[delta == 0]))^2
      )
      sum(w * sq_errs)
    }

    # Marginal bounds
    marg_lo <- marg_bounds[, 1]
    marg_hi <- marg_bounds[, 2]

    ctrl <- list(threshold.stop = 1e-04, maxit = maxit,
                 temperature = temperature, verbose = FALSE)

    range_results <- vector("list", n_ranges)
    range_scores  <- numeric(n_ranges)

    for (r in seq_len(n_ranges)) {
      lower_r <- c(marg_lo, rho_lo[r])
      upper_r <- c(marg_hi, rho_hi[r])
      par_r   <- runif(num_params, lower_r, upper_r)

      s <- GenSA::GenSA(par = par_r, fn = obj_fn,
                        lower = lower_r, upper = upper_r,
                        control = ctrl)

      range_results[[r]] <- c(s$value, s$par)
      range_scores[r]    <- s$value
    }

    c(range_results, list(which.min(range_scores)))
  }

  # ---- Post-processing: majority-vote sub-range selection ----
  genSA_scores <- numeric(B + 1)
  genSA_params <- vector("list", B + 1)
  rg_count     <- integer(n_ranges)
  best_idx     <- integer(B + 1)

  for (i in seq_len(B + 1)) {
    best_idx[i] <- simul_res[[i, n_ranges + 1]]
    rg_count[best_idx[i]] <- rg_count[best_idx[i]] + 1
    genSA_scores[i] <- simul_res[[i, best_idx[i]]][1]
    genSA_params[[i]] <- simul_res[[i, best_idx[i]]][2:(num_params + 1)]
  }
  rho_range <- which.max(rg_count)

  # ---- Extract rho estimates from the winning sub-range ----
  rho_est <- numeric(B + 1)
  for (i in seq_len(B + 1)) {
    rho_est[i] <- simul_res[[i, rho_range]][num_params + 1]
  }

  tau_boot <- rho_est / (rho_est + 2)
  tau_hat  <- tau_boot[B + 1]
  ci_95    <- quantile(tau_boot, probs = c(0.025, 0.975))

  # ---- Return ----
  structure(
    list(
      tau_hat      = tau_hat,
      tau_boot     = tau_boot,
      ci_95        = ci_95,
      rho_hat      = rho_est[B + 1],
      best_range   = rho_range,
      range_counts = rg_count,
      call         = cl_call
    ),
    class = "depEst"
  )
}


#' @export
print.depEst <- function(x, ...) {
  cat("Bagging-MME Dependence Estimation\n")
  cat("----------------------------------\n")
  cat(sprintf("Kendall's tau (point est.): %.4f\n", x$tau_hat))
  cat(sprintf("95%% bootstrap CI:          [%.4f, %.4f]\n",
              x$ci_95[1], x$ci_95[2]))
  cat(sprintf("Clayton rho (point est.):   %.4f\n", x$rho_hat))
  cat(sprintf("Winning sub-range:          %d (votes: %s)\n",
              x$best_range, paste(x$range_counts, collapse = ", ")))
  invisible(x)
}
