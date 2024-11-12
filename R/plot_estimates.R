# Plot results -----------------------------------------------------------------
plot_estimate <- function(estimate, plot_over = FALSE, plot_type = c("gaps", "filled"), xlim = NULL, cond = NULL, ...) {
  #' Plot estimates from self-consistency algorithm
  #'
  #' @param turnbulls Data frame containing Turnbull intervals
  #' @param s Estimates of step changes from EM aglorithm
  #' @param plot_type Either \code{"filled"} for Turnbull intervals filled with
  #' line or \code{"gaps"} for Turnbull intervals left blank.
  #' @param right Only necessary when \code{plot_type} is \code{"filled"} to
  #' get a right upper limit for the plot. Note: Turnbull intervals are where we
  #' don't know the value
  #'
  #' @export
  #' @examples
  #'
  #'

  plot_type <- match.arg(plot_type)

  intervals <- estimate$intervals

  # Turnbull intervals left as gaps
  if (plot_type == "gaps") {
    surv <- estimate$surv

    left <- intervals$II$left
    right <- intervals$II$right
    end <- max(left, right[right != Inf])
    start_points <- c(0, right)
    end_points <- c(left, end)

    if (is.null(xlim)) {
      xlim <- c(0, end)
    }


    if (plot_over) {
      segments(x0 = start_points, y0 = surv, x1 = end_points, y1 = surv, ...)
    } else {
      plot(
        x = 0, y = 0,
        xlim = xlim, ylim = c(0,1), ...
      )
      segments(x0 = start_points, y0 = surv, x1 = end_points, y1 = surv, ...)
    }

  } else {

    left <- intervals$II$left
    right <- intervals$II$right
    # Turnbull intervals filled with diagonal line
    time_points <- sort(c(0, left, right, Inf))
    steps <- estimate$surv

    x <- c(time_points, left)
    ord <- order(x)

    y <- rep(steps, each = 2)
    y <- c(y, steps[-1])

    x <- x[ord]
    y <- y[ord]

    if (plot_over) {
      lines(
        x = x, y = y,
        ylim = c(0,1), type = "l",
        xlab = "Time", ylab = "Survival",
        ...
      )
    } else {
      plot(
        x = x, y = y,
        ylim = c(0,1), type = "l",
        xlab = "Time", ylab = "Survival",
        ...
      )
    }
  }



}
