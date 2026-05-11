# Plot results -----------------------------------------------------------------
#' Plot estimates from self-consistency algorithm
#'
#' @param x ltic object returned from \code{ltic_np}
#' @param cond Index of inner-interval to condition survival on
#' @param plot_type Should a new plot be made (\code{"new"}) or should the
#' plot go over an existing one (\code{"over"})?
#' @param ends Choose between step-changes at right (\code{"r"}) or left
#' (\code{"l"}) end-points or inner-intervals.
#' @param ... Other arguments supplied to plot function.
#'
#' @export
#' @examples
#' est <- ltic_np(mhcps$Ui, mhcps$Vi, mhcps$Ti)
#' plot(est, cond = 1)
#'
plot.ltic <- function(
  x, cond = NULL, plot_type = c("new", "over"), ends = c("r", "l"), ...
) {

  est <- x

  plot_type <- match.arg(plot_type)
  ends <- match.arg(ends)

  int <- est$intervals
  ii_left <- int$II$left
  ii_right <- int$II$right
  surv_type <- est$type

  max_x <- max(ii_right[ii_right != Inf])
  min_x <- min(int$Ti$left)

  if (!is.null(cond)) {
    if (surv_type == "expo") {
      len <- length(est$res$lambda)
      adj <- rep(est$res$lambda[cond + 1], len - cond - 1)
      cond_lambda <- c(est$res$lambda[1:(cond + 1)], adj)
      steps <- exp(-est$res$lambda + cond_lambda)
    } else if (surv_type == "prodlim") {
      len <- length(est$res$lambda)
      adj <- rep(0, cond)
      cond_lambda <- c(adj, est$res$lambda[(cond + 1):len])
      steps <- c(1, cumprod(1 - cond_lambda))
    } else if (surv_type == "surv") {
      steps_0 <- 1 - cumsum(est$res$s)
      len <- length(steps_0)
      adj <- rep(1, cond)
      fact <- 1 - sum(est$res$s[1:cond])
      steps <- c(1, c(adj, steps_0[(cond + 1):len] / fact))
    }

  } else {
    if (surv_type == "expo") {
      steps <- exp(-est$res$lambda)
    } else if (surv_type == "prodlim") {
      steps <- c(1, cumprod(1 - est$res$lambda))
    } else if (surv_type == "surv") {
      steps <- c(1, 1 - cumsum(est$res$s))
    } else {
      steps <- est$res$surv
    }
  }

  # Plot -----------------------------------------------------------------------
  if (ends %in% c("r", "b")) {

    r_time_points <- sort(c(0, rep(ii_right, each = 2), Inf))
    x <- r_time_points
    ord <- order(x)
    y <- rep(steps, each = 2)
    x <- x[ord]
    y <- y[ord]

    if (plot_type == "over") {
      lines(
        x = x, y = y, lwd = 2, ...
      )
    } else {
      plot(
        x = x, y = y, ylim = c(0, 1), xlim = c(min_x, max_x), type = "l",
        lwd = 2, ...
      )
    }
  }

  if (ends %in% c("l", "b")) {

    l_time_points <- sort(c(0, rep(ii_left, each = 2), max(ii_right)))

    x <- l_time_points
    ord <- order(x)
    y <- rep(steps, each = 2)
    x <- x[ord]
    y <- y[ord]

    lines(
      x = x, y = y, lwd = 2, ...
    )
  }
}
