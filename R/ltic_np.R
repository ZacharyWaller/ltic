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
                              max_iterations = 1e5, remove_rcens = TRUE,
                              open_L = TRUE, open_T = FALSE, constr = NULL,
                              method = c("comb", "em", "newton", "both")) {


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

  gamma_int <- data.frame(
      left = 0,
      right = intervals$Oi$left,
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
  if (remove_rcens & any(right == Inf)) {
    right_cens <- alpha[right == Inf, , drop = FALSE]

    if (class(right_cens)[1] == "integer") {
      right_cens <- t(matrix(right_cens))
    }

    leaving <- colSums(right_cens)
    y_0 <- y_0 - leaving
  }

  # Event matrix ---------------------------------------------------------------
  # n potential events in interval
  if (remove_rcens) {
    event <- alpha[right != Inf, , drop = FALSE]
  } else {
    event <- alpha
  }
  if (class(event)[1] == "integer") {
    event <- t(matrix(event))
  }

  # Initialise values ----------------------------------------------------------
  # if (is.null(init)) {
  #   # good first guess for events
  #   n_events <- apply(event / apply(event, 1, sum), 2, sum)
  #   n_cum <- cumsum(n_events)
  #   y_risk <- y_0 - c(0, n_cum[-length(n_cum)])
  #   lambda <- n_events / y_risk
  #   lambda[y_risk == 0] <- 1
  # }

  l <- apply(event, 1, function(x) min(which(x == 1)) - 1)
  r <- apply(event, 1, function(x) max(which(x == 1)))
  t <- apply(beta, 1, function(x) min(which(x == 1)) - 1)

  lambda <- rep(1 / ncol(alpha), ncol(alpha))
  lambda[is.nan(lambda)] <- 9999

  cat("\n Starting algorithm...")
  # Call maximisation algorithm ------------------------------------------------
  if (method == "newton") {
    res <- newton_algorithm(lambda, l, r, deriv_1_0)
  } else if (method == "em") {
    res <- em_algorithm(lambda, l, r, y_0)
  } else if (method == "comb") {
    res <- combined_algorithm(lambda, l, r, deriv_1_0, y_0)
  } else if (method == "both") {
    res <- ltic_r(lambda, l, r, t, y_0)
  }

  list(
    res = res,
    l = l,
    r = r,
    t = t,
    deriv_1_0 = deriv_1_0,
    y_0 = y_0
  )
}
