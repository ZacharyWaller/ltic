#' Non-parametric estimate for left-truncated interval-censored (LTIC) data
#'
#' @description
#' Calculate the non-parametric maximum-likelihood estimate of left-truncated
#' interval-censored (LTIC) data
#'
#' @param left Left end-point of observation interval. Last time of negative
#' diagnoses. -\code{Inf} for let-censored data.
#' @param right Right end-point of observation interval. First time of positive
#' diagnoses. \code{Inf} for right-censored data.
#' @param trunc Truncation (entry) time. May be 0 or \code{NULL} for no
#' truncation.
#' @param tol Tolerance used for convergence of maximisation algorithm.
#' @param max_it Maximum number of iterations.
#' @param open_L Open interval left end-point. See details.
#' @param open_T Open interval from truncation time. See details.
#' @param method Maximization method. See details.
#'
#' @details
#' If \code{trunc} is supplied then it must be supplied for all observations.
#' The default observation intervals are open on the left and closed on the
#' right (L, R] as suggested by Ng (2002). Truncation times have been a cause of
#' continued confusion in respect to interval-censored data, see Yu (2023) for
#' an in-depth discussion. Here truncation times are treated as closed.
#'
#' The maximization methods available via the \code{method} argument include a
#' novel product-limit EM algorithm (\code{prodlim}) and combined EM-ICM
#' algorithm based on it (\code{both}). Efficient implementations of a number
#' of existing algorithms are also available, including EM algorithms
#' (\code{turnbull, shen, yu, breslow}), modified ICM algorithms
#' (\code{icm, surv}), EM-ICM algorithms (\code{turn_comb, yu_comb,
#' bres_comb}) and a quasi-Newton algorithm (\code{optim}). The product-limit
#' EM-ICM is the default as it is the fastest and most accurate across most
#' situations. See Waller et al. (2025) for details.
#'
#' @author Zachary Waller
#'
#' @references
#' Ng, M. P. (2002). A modification of Peto's nonparametric estimation of survival curves for interval-censored data. Biometrics, 58(2), 439-442.
#'
#' Yu, Q. (2023). The generalised MLE with truncated interval-censored data. Journal of Nonparametric Statistics, 35(2), 266-282.
#'
#' Waller, Z., Marshall, A. H., Kee, F., & Lamrock, F. (2025). A fast and stable algorithm for calculating the non-parametric maximum likelihood estimate of left-truncated and interval-censored data.
#'
#' @export
#'
#' @examples
#' trunc <- c(0, 2, 4)
#' left <- c(3, 5, 7)
#' right <- c(Inf, 10, 10)
#' ltic_np(left, right, trunc)
#'
ltic_np <- function(left, right, trunc = NULL, tol = 1e-7,
  max_it = 1e5, open_L = TRUE, open_T = FALSE,
  method = c("both", "turnbull", "shen",
             "yu", "breslow", "prodlim", "bres_comb", "binomial", "optim",
             "turn_comb", "yu_comb", "icm", "surv", "new")
) {

  method <- match.arg(method)

  # truncation times
  if (is.null(trunc)) {
    trunc <- rep(0, length(left))
  } else if (length(trunc) != length(left)) {
    stop("Error: length(trunc) != length(left) - if supplied, there must be as
         many truncation times as left interval times.")
  }

  # calculate inner intervals
  intervals <- inner_intervals(
    left, right, trunc,
    open_L = open_L, open_T = open_T
  )

  alpha <- indicator_matrix(intervals$II, intervals$Oi)
  beta <- indicator_matrix(intervals$II, intervals$Ti)

  if (method %in% c("prodlim", "both", "breslow", "new", "bres_comb")) {
    remove_rcens <- TRUE
  } else {
    remove_rcens <- FALSE
  }

  r_star <- intervals$Oi$right

  if (remove_rcens) {
    r_star[r_star == Inf] <- intervals$Oi$left[r_star == Inf]
  }
  gamma_int <- data.frame(
    left = intervals$Ti$left,
    right = r_star,
    open_left = FALSE,
    open_right = !open_L
  )
  gamma <- indicator_matrix(intervals$II, gamma_int)
  deriv_1_0 <- colSums(gamma)

  # Left truncation ------------------------------------------------------------
  # number of participants who have entered. One for each inner interval
  y_0 <- colSums(beta)

  # Right censoring ------------------------------------------------------------
  # number of participants at each time point
  r_cens <- intervals$Oi$right == Inf
  if (remove_rcens & any(r_cens)) {

    right_cens <- alpha[r_cens, , drop = FALSE]

    if (class(right_cens)[1] == "integer") {
      right_cens <- t(matrix(right_cens))
    }

    leaving <- colSums(right_cens)
    y_0 <- y_0 - leaving
  }

  # Event matrix ---------------------------------------------------------------
  # n potential events in interval
  if (remove_rcens) {
    event <- alpha[!r_cens, , drop = FALSE]
  } else {
    event <- alpha
  }
  if (class(event)[1] == "integer") {
    event <- t(matrix(event))
  }

  # Initialise values ----------------------------------------------------------
  l <- apply(event, 1, function(x) min(which(x == 1)) - 1)
  r <- apply(event, 1, function(x) max(which(x == 1)))
  if (remove_rcens & any(r_cens)) {
    t <- apply(
      beta[!r_cens, , drop = FALSE],
      1, function(x) min(which(x == 1)) - 1)
  } else {
    t <- apply(beta, 1, function(x) min(which(x == 1)) - 1)
  }

  init <- rep(1 / ncol(alpha), ncol(alpha))

  # Checks
  if (length(init) != length(y_0)) {
    message("y_0")
    stop()
  }

  if (length(init) != length(deriv_1_0)) {
    message("deriv_1_0")
    stop()
  }

  # Call maximisation algorithm ------------------------------------------------
  type <- "surv"
  if (method == "both") {
    ## PL-ICM ----
    l_full <- apply(alpha, 1, function(x) min(which(x == 1)) - 1)
    r_full <- apply(alpha, 1, function(x) max(which(x == 1)))
    t_full <- apply(beta, 1, function(x) min(which(x == 1)) - 1)
    cat("\n both \n")
    t0 <- Sys.time()
    res <- ltic_r(init, l, r, t, y_0, l_full, r_full, t_full, tol, max_it)
    time <- Sys.time() - t0
    type <- "expo"

  } else if (method == "surv") {
    ## ICM (surv) ----
    l_full <- apply(alpha, 1, function(x) min(which(x == 1)) - 1)
    r_full <- apply(alpha, 1, function(x) max(which(x == 1)))
    t_full <- apply(beta, 1, function(x) min(which(x == 1)) - 1)
    cat("\n Go \n")
    t0 <- Sys.time()
    res <- ltic_s_r(init, l, r, t, y_0, l_full, r_full, t_full, tol, max_it)
    time <- Sys.time() - t0

  } else if (method == "turn_comb") {
    ## Turnbull-ICM ----
    t0 <- Sys.time()
    res <- ltic_turn_r(init, l, r, t, tol, max_it)
    time <- Sys.time() - t0

  } else if (method == "turnbull") {
    ## Turnbull ----
    t0 <- Sys.time()
    res <- turnbull_r(init, l, r, t, tol, max_it)
    time <- Sys.time() - t0

  } else if (method == "shen") {
    ## Shen ----
    t0 <- Sys.time()
    res <- shen_r(init, l, r, t, tol, max_it)
    time <- Sys.time() - t0

  } else if (method == "yu") {
    ## Yu ----
    t0 <- Sys.time()
    res <- yu_r(init, l, r, t, tol, max_it)
    time <- Sys.time() - t0

  } else if (method == "yu_comb") {
    ## Yu-ICM ----
    t0 <- Sys.time()
    res <- ltic_yu_r(init, l, r, t, tol, max_it)
    time <- Sys.time() - t0

  } else if (method == "breslow") {
    ## Breslow ----
    l_full <- apply(alpha, 1, function(x) min(which(x == 1)) - 1)
    r_full <- apply(alpha, 1, function(x) max(which(x == 1)))
    t_full <- apply(beta, 1, function(x) min(which(x == 1)) - 1)

    t0 <- Sys.time()
    res <- breslow_r(init, l, r, t, deriv_1_0, l_full, r_full, t_full, tol, max_it)
    time <- Sys.time() - t0
    type <- "expo"

  } else if (method == "prodlim") {
    ## Prodlim ----
    l_full <- apply(alpha, 1, function(x) min(which(x == 1)) - 1)
    r_full <- apply(alpha, 1, function(x) max(which(x == 1)))
    t_full <- apply(beta, 1, function(x) min(which(x == 1)) - 1)

    t0 <- Sys.time()
    res <- prodlim_r(init, l, r, t, y_0, l_full, r_full, t_full, tol, max_it)
    time <- Sys.time() - t0
    type <- "prodlim"

  } else if (method == "new") {
    ## Experimental method using Laurent expansion ----
    l_full <- apply(alpha, 1, function(x) min(which(x == 1)) - 1)
    r_full <- apply(alpha, 1, function(x) max(which(x == 1)))
    t_full <- apply(beta, 1, function(x) min(which(x == 1)) - 1)

    t0 <- Sys.time()
    res <- test_r(init, l, r, t, y_0, l_full, r_full, t_full, tol, max_it)
    time <- Sys.time() - t0

  } else if (method == "bres_comb") {
    ## Breslow-ICM ----
    l_full <- apply(alpha, 1, function(x) min(which(x == 1)) - 1)
    r_full <- apply(alpha, 1, function(x) max(which(x == 1)))
    t_full <- apply(beta, 1, function(x) min(which(x == 1)) - 1)

    t0 <- Sys.time()
    res <- bres_comb_r(init, l, r, t, deriv_1_0, l_full, r_full, t_full, tol, max_it)
    time <- Sys.time() - t0
    type <- "expo"

  } else if (method == "binomial") {
    ## Binomial ----
    l_full <- apply(alpha, 1, function(x) min(which(x == 1)) - 1)
    r_full <- apply(alpha, 1, function(x) max(which(x == 1)))
    t_full <- apply(beta, 1, function(x) min(which(x == 1)) - 1)
    t0 <- Sys.time()
    res <- binomial_r(init, l, r, t, deriv_1_0, l_full, r_full, t_full, tol, max_it)
    time <- Sys.time() - t0

  } else if (method == "optim") {
    ## optim() ----
    t0 <- Sys.time()
    res <- optim_method(init, alpha, beta, tol)
    time <- Sys.time() - t0
    type <- "surv"

  } else if (method == "icm") {
    ## ICM (hazards) ----
    t0 <- Sys.time()
    res <- icm_r(init, l, r, t, tol, max_it)
    time <- Sys.time() - t0
    type <- "expo"
  }

  # Return results -------------------------------------------------------------
  output <- list(
    res = res,
    l = l,
    r = r,
    t = t,
    time = time,
    intervals = intervals,
    method = method,
    type = type
  )

  class(output) <- "ltic"

  output
}
