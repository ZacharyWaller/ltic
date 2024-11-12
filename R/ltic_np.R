#' Calculate maximum likelihood of survival function
#'
#' @param left
#' @param right
#' @param trunc
#' @param tol
#' @param init
#' @param max_iterations
#' @param remove_rcens
#' @param open_L
#' @param open_T
#' @param constr
#'
#' @return
#' @export
#'
#' @examples
ltic_np <- function(left, right, trunc = NULL, tol = 1e-7, init = NULL,
  max_it = 1e5, remove_rcens = TRUE, open_L = TRUE, open_T = FALSE,
  constr = NULL,
  method = c("both", "turnbull", "shen",
             "yu", "breslow", "prodlim", "bres_comb", "binomial", "optim",
             "turn_comb", "yu_comb", "icm")
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
  remove <- colSums(alpha) == 0
  alpha <- alpha[, !remove, drop = FALSE]
  beta <- indicator_matrix(intervals$II, intervals$Ti)
  beta <- beta[, !remove, drop = FALSE]

  r_star <- intervals$Oi$right
  r_star[r_star == Inf] <- intervals$Oi$left[r_star == Inf]
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
  r_cens <- alpha[, ncol(alpha)] == 1
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


  cat("\n Starting algorithm...\n")
  # Call maximisation algorithm ------------------------------------------------
  if (method == "both") {
    l_full <- apply(alpha, 1, function(x) min(which(x == 1)) - 1)
    r_full <- apply(alpha, 1, function(x) max(which(x == 1)))
    t_full <- apply(beta, 1, function(x) min(which(x == 1)) - 1)
    cat("\n Go \n")
    t0 <- Sys.time()
    res <- ltic_r(init, l, r, t, y_0, l_full, r_full, t_full, tol, max_it)
    time <- Sys.time() - t0

  } else if (method == "turn_comb") {
    t0 <- Sys.time()
    res <- ltic_turn_r(init, l, r, t, tol, max_it)
    time <- Sys.time() - t0

  } else if (method == "turnbull") {
    t0 <- Sys.time()
    res <- turnbull_r(init, l, r, t, tol, max_it)
    time <- Sys.time() - t0

  } else if (method == "shen") {
    t0 <- Sys.time()
    res <- shen_r(init, l, r, t, tol, max_it)
    time <- Sys.time() - t0

  } else if (method == "yu") {
    t0 <- Sys.time()
    res <- yu_r(init, l, r, t, tol, max_it)
    time <- Sys.time() - t0

  } else if (method == "yu_comb") {
    t0 <- Sys.time()
    res <- ltic_yu_r(init, l, r, t, tol, max_it)
    time <- Sys.time() - t0

  } else if (method == "breslow") {
    l_full <- apply(alpha, 1, function(x) min(which(x == 1)) - 1)
    r_full <- apply(alpha, 1, function(x) max(which(x == 1)))
    t_full <- apply(beta, 1, function(x) min(which(x == 1)) - 1)
    t0 <- Sys.time()
    res <- breslow_r(init, l, r, t, deriv_1_0, l_full, r_full, t_full, tol, max_it)
    time <- Sys.time() - t0

  } else if (method == "prodlim") {
    l_full <- apply(alpha, 1, function(x) min(which(x == 1)) - 1)
    r_full <- apply(alpha, 1, function(x) max(which(x == 1)))
    t_full <- apply(beta, 1, function(x) min(which(x == 1)) - 1)
    cat("\n Go \n")
    t0 <- Sys.time()
    res <- prodlim_r(init, l, r, t, y_0, l_full, r_full, t_full, tol, max_it)
    time <- Sys.time() - t0

  } else if (method == "bres_comb") {
    l_full <- apply(alpha, 1, function(x) min(which(x == 1)) - 1)
    r_full <- apply(alpha, 1, function(x) max(which(x == 1)))
    t_full <- apply(beta, 1, function(x) min(which(x == 1)) - 1)
    t0 <- Sys.time()
    res <- bres_comb_r(init, l, r, t, deriv_1_0, l_full, r_full, t_full, tol, max_it)
    time <- Sys.time() - t0

  } else if (method == "binomial") {
    l_full <- apply(alpha, 1, function(x) min(which(x == 1)) - 1)
    r_full <- apply(alpha, 1, function(x) max(which(x == 1)))
    t_full <- apply(beta, 1, function(x) min(which(x == 1)) - 1)
    t0 <- Sys.time()
    res <- binomial_r(init, l, r, t, deriv_1_0, l_full, r_full, t_full)
    time <- Sys.time() - t0

  } else if (method == "optim") {
    t0 <- Sys.time()
    res <- optim_method(init, alpha, beta, tol)
    time <- Sys.time() - t0
  } else if (method == "icm") {
    t0 <- Sys.time()
    res <- icm_r(init, l, r, t, tol, max_it)
    time <- Sys.time() - t0
  }
  cat("done")
  list(
    res = res,
    l = l,
    r = r,
    t = t,
    deriv_1_0 = deriv_1_0,
    y_0 = y_0,
    time = time
  )
}
