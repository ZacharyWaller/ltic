# Plot results -----------------------------------------------------------------
plot_estimate <- function(estimate, cond = NULL, ...) {
  #' Plot estimates from self-consistency algorithm
  #'
  #' @param estimate ltic object returned from \code{ltic_np}
  #' @param cond index of inner-interval to condition survival on
  #' @param ... Other arguments supplied to plot function.
  #'
  #' @export
  #' @examples
  #'

  int <- estimate$intervals
  ii_left <- int$II$left
  ii_right <- int$II$right

  max_x <- max(ii_right[ii_right != Inf])
  min_x <- min(int$Ti$left)


  if (!is.null(cond)) {
    len <- length(estimate$res$lambda)
    adj <- rep(estimate$res$lambda[cond], len - cond)
    cond_lambda <- c(estimate$res$lambda[1:cond], adj)
    steps <- exp(-estimate$res$lambda + cond_lambda)
  } else {
    steps <- exp(-estimate$res$lambda)
  }
  
  l_time_points <- sort(c(0, rep(ii_left, each = 2), max(ii_right)))

  x <- l_time_points
  ord <- order(x)
  y <- rep(steps, each = 2)
  x <- x[ord]
  y <- y[ord]

  plot(
    x = x, y = y, ylim = c(0, 1), xlim = c(min_x, max_x), type = "l", 
    xlab = "Time", ylab = "Survival", 
    col = "red", lwd = 2, ...
  )

  r_time_points <- sort(c(0, rep(ii_right, each = 2), Inf))
  x <- r_time_points
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]

  lines(
    x = x, y = y, col = "red", lwd = 2
  )

}
