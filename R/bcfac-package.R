## usethis namespace: start
#' @importFrom rlang .data
#' @importFrom Rcpp sourceCpp
#' @useDynLib bcfac, .registration = TRUE
## usethis namespace: end
NULL

#' bcfac: Bayesian Additive Regression Trees for Confounder Selection
#'
#' Fit Bayesian Regression Additive Trees (BART) models to
#' select relevant confounders among a large set of potential confounders
#' and to estimate average treatment effect. For more information, see
#' Kim et al. (2022).
#'
#' Functions in `bcfac` serve one of three purposes.
#' \enumerate{
#'   \item Functions for fitting: `sbart()`, `mbart()` and `bcf()`.
#'   \item Functions for summary: `summary()`, `plot()` and `gelman_rubin()`.
#'   \item Utility function for OpenMP: `count_omp_thread()`.
#' }
#'
#' @references
#' Kim, C., Tec, M., & Zigler, C. M. (2022).
#' Bayesian Nonparametric Adjustment of Confounding.
#' *arXiv preprint arXiv:2203.11798*.
#' \doi{10.48550/arXiv.2203.11798}
#'
#' @docType package
#' @name bcfac-package
NULL
