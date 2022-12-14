% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcfac-package.R
\docType{package}
\name{bcfac-package}
\alias{bcfac-package}
\title{bcfac: Bayesian Additive Regression Trees for Confounder Selection}
\description{
Fit Bayesian Regression Additive Trees (BART) models to
select relevant confounders among a large set of potential confounders
and to estimate average treatment effect. For more information, see
Kim et al. (2022).
}
\details{
Functions in \code{bcfac} serve one of three purposes.
\enumerate{
\item Functions for fitting: \code{sbart()}, \code{mbart()} and \code{bcf()}.
\item Functions for summary: \code{summary()}, \code{plot()} and \code{gelman_rubin()}.
\item Utility function for OpenMP: \code{count_omp_thread()}.
}
}
\references{
Kim, C., Tec, M., & Zigler, C. M. (2022).
Bayesian Nonparametric Adjustment of Confounding.
\emph{arXiv preprint arXiv:2203.11798}.
\doi{10.48550/arXiv.2203.11798}
}
