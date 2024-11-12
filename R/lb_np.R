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
lb_np <- function(left, right, tol = 1e-7, init = NULL,
  max_it = 1e5, open_L = TRUE, open_T = FALSE,
  constr = NULL
) {

  # calculate inner intervals
  intervals <- inner_intervals(
    left, right,
    open_L = open_L
  )

  max_t <- max(intervals$II$left)
  delta_t <- c(intervals$II$left, Inf) - c(0, intervals$II$left)
  delta_t <- delta_t[-length(delta_t)] / max_t

  alpha <- indicator_matrix(intervals$II, intervals$Oi)
  remove <- colSums(alpha) == 0
  alpha <- alpha[, !remove, drop = FALSE]
  beta <- indicator_matrix(intervals$II, intervals$Ti)
  beta <- beta[, !remove, drop = FALSE]

  y_0 <- colSums(beta)
  # Event matrix ---------------------------------------------------------------
  # n potential events in interval
  event <- alpha
  if (class(event)[1] == "integer") {
    event <- t(matrix(event))
  }

  # Initialise values ----------------------------------------------------------
  l <- apply(event, 1, function(x) min(which(x == 1)) - 1)
  r <- apply(event, 1, function(x) max(which(x == 1)))
  t <- apply(beta, 1, function(x) min(which(x == 1)) - 1)

  init <- rep(1 / ncol(alpha), ncol(alpha))


  cat("\n Starting algorithm...\n")
  # Call maximisation algorithm ------------------------------------------------
  l_full <- apply(alpha, 1, function(x) min(which(x == 1)) - 1)
  r_full <- apply(alpha, 1, function(x) max(which(x == 1)))
  t_full <- apply(beta, 1, function(x) min(which(x == 1)) - 1)
  cat("\n Go \n")
  t0 <- Sys.time()
  res <- length_bias_r(init, l, r, t, y_0, delta_t, tol, max_it)
  time <- Sys.time() - t0


  cat("done")
  list(
    res = res,
    l = l,
    r = r,
    t = t,
    y_0 = y_0,
    time = time
  )
}
