# depCenR <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/zionstar12/Bagging-MME/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zionstar12/Bagging-MME/actions/workflows/R-CMD-check.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

## Overview

**depCenR** implements the **Bagging-MME** (Bootstrap-Aggregated Method-of-Moments Estimation) framework for estimating correlation and marginal distributions from **dependently censored survival data**. It provides tools for:
  
- Test whether dependence exists between a survival time *T* and a censoring time *C* (i.e., Kendall's τ = 0 vs. τ > 0).
- Estimate the strength of that dependence (Kendall's τ with bootstrap 95% CIs).
- Augmented estimation for small samples via diffusion-model–based data augmentation.
> *(Under development)*

The methodology is described in:

> **Method-of-moments bootstrap-aggregating for the estimation of correlation and marginal distributions in dependently censored survival data**
> *(viewable at arxiv.org)*

## Installation

The development version of **depCenR** will be installable from GitHub as:

```r
# install.packages("devtools")
devtools::install_github("zionstar12/Bagging-MME")
```

## Main Functions

| Function | Purpose |
|---|---|
| `depTest()` | Test for the presence of dependence between *T* and *C* (H₀: τ = 0 vs. H₁: τ > 0) via the Bagging-MME procedure. |
| `depEst()` | Estimate Kendall's τ with bootstrap 95% confidence intervals using the full Bagging-MME algorithm. |
| `depEst.aug()` | Augmented estimation for small samples: leverages diffusion-model–generated synthetic data to improve estimation of low-to-moderate correlations (τ ≈ 0.3–0.5). *(Under development)* |

## Quick Start

```r
library(depCenR)

# --- Example: Estimate dependence ---
# Suppose 'obs_time' is the observed (minimum) time,
# and 'event' is the event indicator (1 = event, 0 = censored).
result <- depEst(time = obs_time, event = event, B = 200, n_cores = 10)
print(result)

# --- Example: Test for dependence ---
test_result <- depTest(time = obs_time, event = event, B = 200, n_cores = 10)
print(test_result)

# --- Example: Small-sample augmented estimation ---
aug_result <- depEst.aug(time = obs_time, event = event, B = 200,
                         n_aug = 500, n_cores = 10)
print(aug_result)
```

## Methodology

The Bagging-MME approach works by:

1. **Bootstrap resampling** — Drawing *B* bootstrap samples from the observed data.
2. **Computing target summary statistics** — For each bootstrap sample, calculating the observed event proportion, conditional means, SDs, and skewness of the event and censored subgroups.
3. **Stochastic global optimization** — Using Generalized Simulated Annealing (GenSA) to find the bivariate generalized-gamma parameters (including a parametric copula dependence parameter) that best reproduce the target statistics, across partitioned dependence-parameter sub-ranges.
4. **Aggregation** — Selecting the best-fitting sub-range by majority vote and aggregating the bootstrap estimates of Kendall's τ.

## Conference Presentation

A PDF of the oral presentation introducing the Bagging-MME proposal, delivered at the International Biometric Conference (IBC) 2024, is available in [`inst/presentations/`](inst/presentations/).

## Dependencies

depCenR relies on the following R packages:

- [GenSA](https://CRAN.R-project.org/package=GenSA) — Generalized Simulated Annealing
- [copula](https://CRAN.R-project.org/package=copula) — Multivariate dependence via copulas
- [flexsurv](https://CRAN.R-project.org/package=flexsurv) — Flexible parametric survival models
- [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/parallel-package.html) / [doParallel](https://CRAN.R-project.org/package=doParallel) — Parallel computation
- [compound.Cox](https://CRAN.R-project.org/package=compound.Cox) — Univariate survival fits

## License

GPL-3. See [LICENSE](LICENSE.md) for details.

## Citation

If you use **depCenR** in published work, please cite:

```
@article{depCenR,
  title   = {Method-of-moments bootstrap-aggregating for the estimation of
             correlation and marginal distributions in dependently censored
             survival data},
  author  = {[Hyun-Soo Zhang, Inkyung Jung, and Chung Mo Nam]},
  journal = {[https://arxiv.org/abs/2604.04032]},
  year    = {2026},
  note    = {R package version 0.1.0}
}
```
