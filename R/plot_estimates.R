# Plot results -----------------------------------------------------------------
plot_estimate <- function(estimate, surv_type, cond = NULL, plot_type = "new", ...) {
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
    if (surv_type == "exp") {
      len <- length(estimate$res$lambda)
      adj <- rep(estimate$res$lambda[cond], len - cond)
      cond_lambda <- c(estimate$res$lambda[1:cond], adj)
      steps <- exp(-estimate$res$lambda + cond_lambda)
    } else if (surv_type == "prod") {
      len <- length(estimate$res$lambda)
      adj <- rep(0, cond)
      cond_lambda <- c(adj, estimate$res$lambda[cond:len])
      steps <- c(1, cumprod(1 - cond_lambda))
    } else if (surv_type == "step") {
      len <- length(estimate$res$s)
      adj <- rep(1, cond)
      fact <- 1 - sum(estimate$res$s[1:cond])
      cond_s <- c(adj, 1 - cumsum(estimate$res$s[cond:len] / fact))
      steps <- c(1, cond_s)
    }

  } else {
    if (surv_type == "exp") {
      steps <- exp(-estimate$res$lambda)
    } else if (surv_type == "prod") {
      steps <- c(1, cumprod(1 - estimate$res$lambda))
    } else if (surv_type == "step") {
      steps <- c(1, 1 - cumsum(estimate$res$s))
    } else {
      steps <- estimate$res$surv
    }
  }
  
  l_time_points <- sort(c(0, rep(ii_left, each = 2), max(ii_right)))

  x <- l_time_points
  ord <- order(x)
  y <- rep(steps, each = 2)
  x <- x[ord]
  y <- y[ord]

  if (plot_type == "over") {
    lines(
      x = x, y = y, col = "red", lwd = 2
    )
  } else {
    plot(
      x = x, y = y, ylim = c(0, 1), xlim = c(min_x, max_x), type = "l",
      xlab = "Time", ylab = "Survival",
      col = "red", lwd = 2, ...
    )
  }


  r_time_points <- sort(c(0, rep(ii_right, each = 2), Inf))
  x <- r_time_points
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]

  lines(
    x = x, y = y, col = "red", lwd = 2
  )

}
