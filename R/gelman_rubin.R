#' Gelman-Rubin diagnostic for \code{bcfac} objects.
#'
#' `gelman_rubin()` computes Gelman-Rubin diagnostic for `bcfac` objects.
#'
#' @param x A `bcfac` object.
#'
#' @return Gelman-Rubin diagnostic value.
#' @examples
#' data(ihdp, package = "bcfac")
#' x <- mbart(
#'   Y               = ihdp$y_factual,
#'   trt             = ihdp$treatment,
#'   X               = ihdp[, 6:30],
#'   num_tree        = 10,
#'   num_chain       = 2,
#'   num_post_sample = 20,
#'   num_burn_in     = 10,
#'   verbose         = FALSE
#' )
#'
#' gelman_rubin(x)
#'
#' @export
gelman_rubin <- function(x) {
  UseMethod("gelman_rubin")
}

#' @exportS3Method
gelman_rubin.bcfac <- function(x) {
  # stop if x or y is not bcfac object
  num_chain <- x$params$num_chain
  if (num_chain == 1) {
    message("`num_chain` must be greater than 1 for Gelman-Rubin diagnostics.")
    return (NA)
  }

  # placeholder for sequence result
  seq_mean <- rep(0, 2 * num_chain)
  seq_var  <- rep(0, 2 * num_chain)

  # compute mean and variance of each sequence
  for (i in seq_len(num_chain)) {
    estimand        <- x$chains[[i]]$ATE
    num_post_sample <- x$params$num_post_sample
    seq_length      <- num_post_sample %/% 2
    first_half      <- 1:seq_length
    second_half     <- (seq_length + 1):num_post_sample

    seq_mean[2 * i - 1] <- mean(estimand[first_half])
    seq_mean[2 * i]     <- mean(estimand[second_half])
    seq_var[2 * i - 1]  <- stats::var(estimand[first_half])
    seq_var[2 * i]      <- stats::var(estimand[second_half])
  }

  # compute between-sequence variance and within-sequence variance
  between_var <- num_chain * stats::var(seq_mean)
  within_var  <- mean(seq_var)

  # compute marginal posterior variance
  mar_post_var <- within_var * (seq_length - 1) / seq_length +
    between_var / seq_length

  # return
  sqrt(mar_post_var / within_var)
}
