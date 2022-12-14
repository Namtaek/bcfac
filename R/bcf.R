#' @rdname bart
#' @usage NULL
#' @export
bcf <- function(
  Y, trt, X,
  trt_treated     = 1,
  trt_control     = 0,
  num_tree        = 100,
  num_chain       = 4,
  num_burn_in     = 100,
  num_thin        = 10,
  num_post_sample = 200,
  step_prob       = c(0.28, 0.28, 0.44),
  alpha           = 0.95,
  beta            = 2,
  nu              = 3,
  alpha2          = 0.25,
  beta2           = 3,
  nu2             = 10,
  num_tree_mod    = 25,
  q               = 0.95,
  dir_alpha       = 5,
  boot_size       = NULL,
  parallel        = NULL,
  verbose         = TRUE
) {

  # ---- check input ----
  check_input2(
    Y, trt, X, trt_treated, trt_control,
    num_tree, num_chain,
    num_burn_in, num_thin, num_post_sample,
    step_prob, alpha, beta, nu,
    alpha2, beta2, nu2, num_tree_mod,
    q, dir_alpha, verbose
  )



  # ---- data preprocessing ----
  n <- nrow(X)
  p <- ncol(X)

  # check for factor variable then change it to dummy variables
  if (sum(vapply(X, is.factor, TRUE))) {
    X <- fct_to_dummy(X)
    p <- ncol(X)
  }

  # convert to numeric vector and matrix
  if (!is.numeric(Y))
    Y <- as.numeric(Y)
  if (!is.numeric(trt))
    trt <- as.numeric(trt)
  if (!is.matrix(X))
    X <- as.matrix(X)

  # shift and rescale to [-0.5, 0.5]
  Y_mean <- mean(Y)
  Y      <- Y - Y_mean
  Y_max  <- max(Y)
  Y_min  <- min(Y)
  Y      <- (Y - Y_min) / (Y_max - Y_min) - 0.5

  # initialize bootstrap sample size and parallel
  if (is.null(boot_size))
    boot_size <- 2 * n
  if (is.null(parallel))
    parallel <- ifelse(n < 15e4, TRUE, FALSE)

  # assign variable names if there are no name
  if (is.null(colnames(X)))
    colnames(X) <- paste0("X", seq_len(p))


  # ---- bcf specific preprocessing step ----
  # check whether it is binary treatment
  is_binary_trt <- isTRUE(
    all.equal(sort(unique(trt)), sort(c(trt_treated, trt_control)))
  )


  # ---- calculate lambda before MCMC iterations ----
  sigma2_exp <- ifelse(is_binary_trt, 1, stats::var(Y))
  sigma2_out <- stats::var(Y)


  if (is_binary_trt) {
    lambda_exp <- 0 # arbitrary value
  } else {
    f <- function(lambda) {
      invgamma::qinvgamma(
        q, nu / 2, rate = lambda * nu / 2, lower.tail = TRUE, log.p = FALSE
      ) - sqrt(sigma2_exp)
    }
    lambda_exp <- rootSolve::uniroot.all(f, c(0.1^5, 10))
  }

  f <- function(lambda) {
    invgamma::qinvgamma(
      q, nu / 2, rate = lambda * nu / 2, lower.tail = TRUE, log.p = FALSE
    ) - sqrt(sigma2_out)
  }
  lambda_out <- rootSolve::uniroot.all(f, c(0.1^5, 10))

  f <- function(lambda) {
    invgamma::qinvgamma(
      q, nu2 / 2, rate = lambda * nu2 / 2, lower.tail = TRUE, log.p = FALSE
    ) - sqrt(sigma2_out)
  }
  lambda_mod <- rootSolve::uniroot.all(f, c(0.1^5, 10))


  # ---- run MCMC and save result of each chain ----
  chains         <- list()
  num_chain_iter <- num_burn_in + (num_thin + 1) * num_post_sample
  if (verbose) {
    cat(
      "\n",
      "Fitting ", num_chain, " chains with ", num_chain_iter, " iters each...",
      "\n\n",
      sep = ""
    )
  }

  for (chain_idx in seq_len(num_chain)) {
    # placeholder for MCMC samples
    # Y1 and Y0 are potential outcomes. ATE = Y1 - Y0
    Y1  <- vector(mode = "numeric", length = num_post_sample)
    Y0  <- vector(mode = "numeric", length = num_post_sample)

    # placeholder for inclusion probabilities and initialize var_prob
    var_count <- matrix(0, nrow = num_post_sample, ncol = p)
    var_prob  <- MCMCpack::rdirichlet(1, rep(dir_alpha, p))

    # placeholder for parameters
    sigma2_out_hist    <- vector(mode = "numeric", length = num_chain_iter + 1)
    dir_alpha_hist     <- vector(mode = "numeric", length = num_chain_iter + 1)
    sigma2_out_hist[1] <- stats::var(Y)
    dir_alpha_hist[1]  <- dir_alpha

    # call Rcpp implementation
    fit_bcf_hetero(
      Y1, Y0, var_count, var_prob,
      sigma2_exp, sigma2_out_hist, dir_alpha_hist,
      Y, trt, X, trt_treated, trt_control,
      chain_idx, num_chain,
      num_chain_iter, num_burn_in, num_thin, num_post_sample,
      num_tree, step_prob, alpha, beta, nu,
      alpha2, beta2, nu2, num_tree_mod,
      lambda_exp, lambda_out, lambda_mod,
      boot_size, is_binary_trt, parallel, verbose
    )

    # post-processing for each MCMC chain
    var_prob <- colMeans(ifelse(var_count > 1, 1, 0))

    # rescale result
    Y1  <- (Y1  + 0.5) * (Y_max - Y_min) + Y_min + Y_mean
    Y0  <- (Y0  + 0.5) * (Y_max - Y_min) + Y_min + Y_mean
    ATE <-  Y1 - Y0

    chains[[chain_idx]] <- list(
      ATE = ATE, Y1 = Y1, Y0 = Y0,
      var_count  = var_count,
      var_prob   = var_prob,
      sigma2_out = sigma2_out_hist,
      dir_alpha  = dir_alpha_hist
    )
  }


  # ---- post processing ----
  names(chains) <- paste0("chain", seq_len(num_chain))

  # merge result
  ATE <- Y1 <- Y0 <- NULL
  var_prob <- vector(mode = "numeric", length = p)
  for (chain_idx in seq_len(num_chain)) {
    ATE <- c(ATE, chains[[chain_idx]]$ATE)
    Y1  <- c(Y1,  chains[[chain_idx]]$Y1)
    Y0  <- c(Y0,  chains[[chain_idx]]$Y0)
    var_prob <- var_prob + chains[[chain_idx]]$var_prob
  }
  var_prob        <- var_prob / num_chain
  names(var_prob) <- c(colnames(X))

  cat("\n")

  # return as bcfac object
  structure(
    list(
      ATE = ATE, Y1 = Y1, Y0 = Y0,
      var_prob = var_prob,
      chains   = chains,
      model    = "bcf",
      label    = c(colnames(X)),
      params   = list(
        trt_treated     = trt_treated,
        trt_control     = trt_control,
        num_tree        = num_tree,
        num_chain_iter  = num_chain_iter,
        num_chain       = num_chain,
        num_burn_in     = num_burn_in,
        num_thin        = num_thin,
        num_post_sample = num_post_sample,
        step_prob       = step_prob,
        alpha           = alpha,
        beta            = beta,
        nu              = nu,
        #alpha2          = alpha2,
        #beta2           = beta2,
        #nu2             = nu2,
        #num_tree_mod    = num_tree_mod,
        q               = q
      )
    ),
    class = "bcfac"
  )
}

bcf
