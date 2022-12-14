% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gelman_rubin.R
\name{gelman_rubin}
\alias{gelman_rubin}
\title{Gelman-Rubin diagnostic for \code{bcfac} objects.}
\usage{
gelman_rubin(x)
}
\arguments{
\item{x}{A \code{bcfac} object.}
}
\value{
Gelman-Rubin diagnostic value.
}
\description{
\code{gelman_rubin()} computes Gelman-Rubin diagnostic for \code{bcfac} objects.
}
\examples{
data(ihdp, package = "bcfac")
x <- mbart(
  Y               = ihdp$y_factual,
  trt             = ihdp$treatment,
  X               = ihdp[, 6:30],
  num_tree        = 10,
  num_chain       = 2,
  num_post_sample = 20,
  num_burn_in     = 10,
  verbose         = FALSE
)

gelman_rubin(x)

}
