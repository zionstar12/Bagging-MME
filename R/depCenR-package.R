#' depCenR: Bagging-MME for Dependently Censored Survival Data
#'
#' The \pkg{depCenR} package implements the Bootstrap-Aggregated
#' Method-of-Moments Estimation (Bagging-MME) framework for estimating
#' correlation (Kendall's tau) and marginal distributions from dependently
#' censored survival data.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{depTest}}}{Test for the presence of dependence
#'     (H0: tau = 0 vs H1: tau > 0).}
#'   \item{\code{\link{depEst}}}{Estimate Kendall's tau with bootstrap 95\%
#'     confidence intervals.}
#'   \item{\code{\link{depEst.aug}}}{Augmented estimation for small samples
#'     using diffusion-model data augmentation.}
#' }
#'
#' @docType package
#' @name depCenR-package
#' @aliases depCenR
"_PACKAGE"
