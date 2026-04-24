###############################################################################
# Bagging-MME Simulation Code
# depCenR package — Standalone reproduction script
#
# Purpose:
#   Demonstrates the Bagging-MME algorithm on simulated data of 
#   bivariate generalized-gamma marginals with Clayton copula dependence, 
#   estimating Kendall's tau and bootstrap 95% confidence intervals.
#
# Key optimization:
#   Replaces copula::cCopula() with a direct Clayton conditional inverse 
#   for the inner GenSA loop (thousands of evaluations per run).
#   The conditional CDF of Clayton(theta) is:
#     h(v|u) = u^{-(theta+1)} * (u^{-theta} + v^{-theta} - 1)^{-(1+theta)/theta}
#   Inverting for v given (u1, u2) ~ iid Uniform(0,1):
#     v = ( u^{-theta} * (u2^{-theta/(1+theta)} - 1) + 1 )^{-1/theta}
#
# Author:  Hyun-Soo Zhang
# Date:    2026
# License: GPL-3
###############################################################################

# ---- Environment setup ----
rm(list = ls())

library(parallel)
library(doParallel)
library(copula)       # used only for data generation (not inside the optimizer)
library(GenSA)
library(compound.Cox)
library(flexsurv)
library(moments)

# ---- Cluster configuration ----
n_workers <- min(20, parallel::detectCores() - 1L)
cl <- makeCluster(n_workers)
registerDoParallel(cl)
cat("Registered", getDoParWorkers(), "parallel workers\n")


# ============================================================================
# Section 1: Internal helper functions
# ============================================================================

#' Direct Clayton conditional inverse (vectorized)
#' Replaces copula::cCopula() with pure arithmetic & no dispatch overhead
clayton_inv_cond <- function(u1, u2, theta) {
  if (theta < 1e-10) return(cbind(u1, u2))
  v <- (u1^(-theta) * (u2^(-theta / (1 + theta)) - 1) + 1)^(-1 / theta)
  cbind(u1, v)
}

#' Fast sample skewness (avoids moments::skewness overhead)
skew_fast <- function(x) {
  n <- length(x)
  if (n < 3L) return(0)
  m <- mean(x)
  s <- sd(x)
  if (!is.finite(s) || s <= 0) return(0)
  mean(((x - m) / s)^3)
}

#' Inverse generalised-gamma survival transform (AFT-type parameterization)
inv_gg_surv <- function(u, mu, sigma, q) {
  exp(sigma / q * (2 * log(q) + log(qgamma(1 - u, shape = 1 / q^2))) + mu)
}


# ============================================================================
# Section 2: Objective function and loss weights
# ============================================================================

# Loss weights for the seven summary statistics:
#   (event proportion, mu_T, mu_C, sd_T, sd_C, skew_T, skew_C)
# Higher weight = statistic contributes more to the objective.
LOSS_WEIGHTS <- c(4, 2, 2, 1, 1, 0.5, 0.5)
LOSS_WEIGHTS <- LOSS_WEIGHTS / sum(LOSS_WEIGHTS)

obj.fcn <- function(params) {
  mu.T    <- max(0, params[1])
  mu.C    <- max(0, params[2])
  sigma.T <- max(0, params[3])
  sigma.C <- max(0, params[4])
  q.T     <- params[5]
  q.C     <- params[6]
  rho     <- max(0, params[7])

  uv <- clayton_inv_cond(runif(n, 0, 1), runif(n), rho)

  X1 <- inv_gg_surv(uv[, 1], mu.T, sigma.T, q.T)
  X2 <- inv_gg_surv(uv[, 2], mu.C, sigma.C, q.C)

  delta <- as.integer(X1 <= X2)

  sq_errs <- c(
    (true_pct.T  - mean(delta))^2,
    (true_mu.Tbar - mean(X1[delta == 1]))^2,
    (true_mu.Cbar - mean(X2[delta == 0]))^2,
    (true_sd.Tbar - sd(X1[delta == 1]))^2,
    (true_sd.Cbar - sd(X2[delta == 0]))^2,
    (true_q.Tbar  - skew_fast(X1[delta == 1]))^2,
    (true_q.Cbar  - skew_fast(X2[delta == 0]))^2
  )
  sum(w * sq_errs)
}


# ============================================================================
# Section 3: True parameters and data generation
# ============================================================================

# Clayton dependence scenarios (indexed 1–4)
clayton_cors <- c(1e-09, 0.8571429, 2, 8)
cor_strength <- 4
cor_param    <- clayton_cors[cor_strength]
num_params   <- 7

n <- 1000
true_params <- c(4.23, 4.44, 0.51, 0.29, 0.31, 0.52, cor_param)
true_theta  <- true_params[7]

# Generate bivariate generalized-gamma sample
set.seed(20260320)
u1 <- runif(n)
uv <- cCopula(cbind(u1, runif(n)),
              copula = claytonCopula(true_theta, dim = 2),
              inverse = TRUE)
u <- uv[, 1]; v <- uv[, 2]

T1 <- inv_gg_surv(u, true_params[1], true_params[3], true_params[5])
T2 <- inv_gg_surv(v, true_params[2], true_params[4], true_params[6])

d1 <- as.integer(T1 <= T2)
d2 <- 1L - d1
df <- data.frame(T = pmin(T1, T2), T1 = T1, T2 = T2, d1 = d1, d2 = d2)


# ============================================================================
# Section 4: Marginal parameter bounds: NEEDS SEPARATE ESTIMATION!!
# ============================================================================

marg_par_ranges <- matrix(
  c(4.119, 3.918, 0.276, 0.180, -0.121, -0.181,
    4.359, 4.662, 0.574, 0.359,  1.472,  0.869),
  nrow = 6, ncol = 2,
  dimnames = list(
    c("mu1", "mu2", "sigma1", "sigma2", "q1", "q2"),
    c("25%", "75%")
  )
)

mu1.low    <- marg_par_ranges[1, 1]; mu1.up    <- marg_par_ranges[1, 2]
mu2.low    <- marg_par_ranges[2, 1]; mu2.up    <- marg_par_ranges[2, 2]
sigma1.low <- marg_par_ranges[3, 1]; sigma1.up <- marg_par_ranges[3, 2]
sigma2.low <- marg_par_ranges[4, 1]; sigma2.up <- marg_par_ranges[4, 2]
q1.low     <- marg_par_ranges[5, 1]; q1.up     <- marg_par_ranges[5, 2]
q2.low     <- marg_par_ranges[6, 1]; q2.up     <- marg_par_ranges[6, 2]


# ============================================================================
# Section 5: Clayton-parameter sub-ranges and bootstrap
# ============================================================================

rho_lo   <- c(-0.1818182, 0.3529412, 1.333333, 3.714286)
rho_hi   <- c( 0.3529412, 1.333333,  3.714286, 18)
n_ranges <- 4

# Bootstrap samples (B = 99 resamples + 1 original)
b <- 199
boot_list <- replicate(b, {
  data.frame(df[sample(seq_len(nrow(df)), size = n, replace = TRUE), ])
})
boot_list[[length(boot_list) + 1]] <- df

# Pre-compute target summary statistics for each bootstrap sample
truPct.T <- numeric(length(boot_list))
truMu.T  <- numeric(length(boot_list)); truMu.C <- numeric(length(boot_list))
truSD.T  <- numeric(length(boot_list)); truSD.C <- numeric(length(boot_list))
truQ.T   <- numeric(length(boot_list)); truQ.C  <- numeric(length(boot_list))

for (i in seq_along(boot_list)) {
  bdf <- boot_list[[i]]
  truPct.T[i] <- sum(bdf$d1) / n
  truMu.T[i]  <- mean(bdf$T1[bdf$d1 == 1])
  truMu.C[i]  <- mean(bdf$T2[bdf$d2 == 1])
  truSD.T[i]  <- sd(bdf$T1[bdf$d1 == 1])
  truSD.C[i]  <- sd(bdf$T2[bdf$d2 == 1])
  truQ.T[i]   <- skew_fast(bdf$T1[bdf$d1 == 1])
  truQ.C[i]   <- skew_fast(bdf$T2[bdf$d2 == 1])
}


# ============================================================================
# Section 6: Main parallel loop (Bagging-MME)
# ============================================================================

pkgs <- c("copula", "GenSA", "compound.Cox", "flexsurv", "moments")

elapsed_time <- system.time({
  simul_res <- foreach(
    num = seq_len(b + 1),
    .combine = "rbind",
    .packages = pkgs
  ) %dopar% {

    # ---- Re-define helpers inside each worker ----
    clayton_inv_cond <- function(u1, u2, theta) {
      if (theta < 1e-10) return(cbind(u1, u2))
      v <- (u1^(-theta) * (u2^(-theta / (1 + theta)) - 1) + 1)^(-1 / theta)
      cbind(u1, v)
    }
    skew_fast <- function(x) {
      n <- length(x)
      if (n < 3L) return(0)
      m <- mean(x); s <- sd(x)
      if (!is.finite(s) || s <= 0) return(0)
      mean(((x - m) / s)^3)
    }
    inv_gg_surv <- function(u, mu, sigma, q) {
      exp(sigma / q * (2 * log(q) + log(qgamma(1 - u, shape = 1 / q^2))) + mu)
    }

    w <- LOSS_WEIGHTS

    # Target statistics for this bootstrap sample
    true_pct.T   <- truPct.T[num]
    true_mu.Tbar <- truMu.T[num];  true_mu.Cbar <- truMu.C[num]
    true_sd.Tbar <- truSD.T[num];  true_sd.Cbar <- truSD.C[num]
    true_q.Tbar  <- truQ.T[num];   true_q.Cbar  <- truQ.C[num]

    # Marginal bounds (shared across all sub-ranges)
    marg_lo <- c(mu1.low, mu2.low, sigma1.low, sigma2.low, q1.low, q2.low)
    marg_hi <- c(mu1.up,  mu2.up,  sigma1.up,  sigma2.up,  q1.up,  q2.up)

    # GenSA control parameters
    maxit <- 100
    temp  <- 700
    ctrl  <- list(threshold.stop = 1e-04, maxit = maxit,
                  temperature = temp, verbose = FALSE)

    range_results <- vector("list", n_ranges)
    range_scores  <- numeric(n_ranges)

    for (r in seq_len(n_ranges)) {
      lower_r <- c(marg_lo, rho_lo[r])
      upper_r <- c(marg_hi, rho_hi[r])
      par_r   <- runif(num_params, lower_r, upper_r)

      s <- GenSA(par = par_r, fn = obj.fcn,
                 lower = lower_r, upper = upper_r,
                 control = ctrl)

      range_results[[r]] <- c(s$value, s$par)
      range_scores[r]    <- s$value
    }

    c(range_results, list(which.min(range_scores)))
  }
})[3]

cat(sprintf("Elapsed time: %.2f minutes\n", elapsed_time / 60))


# ============================================================================
# Section 7: Post-processing and evaluation
# ============================================================================

# ---- Majority-vote sub-range selection ----
genSA_scores <- numeric(b + 1)
genSA_params <- replicate(b + 1, rep(0, num_params), simplify = FALSE)
rg_count     <- integer(n_ranges)
best_index   <- integer(b + 1)

for (i in seq_len(b + 1)) {
  best_index[i] <- simul_res[[i, n_ranges + 1]]
  rg_count[best_index[i]] <- rg_count[best_index[i]] + 1
  genSA_scores[i]  <- simul_res[[i, best_index[i]]][1]
  genSA_params[[i]] <- simul_res[[i, best_index[i]]][2:(num_params + 1)]
}
rho_range <- which.max(rg_count)

# ---- Extract rho estimates from winning sub-range ----
rho_est <- numeric(b + 1)
for (i in seq_len(b + 1)) {
  rho_est[i] <- simul_res[[i, rho_range]][num_params + 1]
}

# ---- Convert Clayton rho to Kendall's tau ----
tau_hat <- rho_est / (rho_est + 2)
true_tau <- cor_param / (cor_param + 2)

# ---- Results ----
cat("\n=== Bagging-MME Results ===\n")
cat(sprintf("True Kendall's tau:    %.4f\n", true_tau))
cat(sprintf("Estimated tau (orig.): %.4f\n", tau_hat[b + 1]))
cat(sprintf("Bootstrap mean tau:    %.4f (SD = %.4f)\n",
            mean(tau_hat), sd(tau_hat)))
cat(sprintf("MAE:                   %.4f\n", mean(abs(true_tau - tau_hat))))

ci_95 <- quantile(tau_hat, probs = c(0.025, 0.975))
cat(sprintf("95%% bootstrap CI:      [%.4f, %.4f]\n", ci_95[1], ci_95[2]))
cat(sprintf("Winning sub-range:     %d (counts: %s)\n",
            rho_range, paste(rg_count, collapse = ", ")))

# ---- Cleanup ----
parallel::stopCluster(cl)
