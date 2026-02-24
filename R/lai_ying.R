lai_ying <- function(x, c = 1, alpha = 0.25) {

  if (x$method != "PL-ICM") {
    stop("The Lai-Ying method is currently only supported when using the for 
          estimates taken from the PL-ICM algorithm.")
  }

  n <- nrow(x$intervals$Oi)
  low_lim <- c * n^alpha
  events <- x$res$n
  risk_set <- c(x$R_0 - c(0, cumsum(x$res$n[-length(x$res$n)])))

  remove <- risk_set < low_lim
  events[remove] <- 0
  h <- events / risk_set
  surv <- c(1, cumprod(1 - h))

  list(
    surv = surv,
    events = events,
    risk_set = risk_set,
    n = n,
    ly_remove = remove
  )
}

# Lai-Ying plots ---------------------------------------------------------------
lai_ying_plot <- function(x, c = 1, alpha = 0.25, ...) {

  est <- lai_ying(x, c, alpha)

  par(mfrow = c(2, 1))

  # survival plot ----
  ii_right <- x$intervals$II$right
  r_time_points <- sort(c(0, rep(ii_right, each = 2), Inf))
  x_plot <- r_time_points
  ord <- order(x_plot)

  steps <- est$surv
  y_plot <- rep(steps, each = 2)
  y_plot <- y_plot[ord]
  plot(
    x_plot, y_plot, type = "l", col = "red", lwd = 1.5,
    ylim = c(0, 1),
    xlab = "Age", ylab = "Survival",
    ...
  )

  # risk-set plot ----
  # events
  steps <- c(0, est$events)
  y_plot <- rep(steps, each = 2)
  y_plot <- y_plot[ord]
  max_y <- max(est$events, est$risk_set)
  plot(
    x_plot, y_plot, type = "l", col = "red", lwd = 1.5,
    ylim = c(0, max_y),
    xlab = "Age", ylab = "Expected Number",
    ...
  )

  # risk set
  steps <- c(0, est$risk_set)
  y_plot <- rep(steps, each = 2)
  x_plot <- x_plot[ord]
  y_plot <- y_plot[ord]
  lines(x_plot, y_plot, type = "l", lwd = 1.5, lty = 3)

  # highlight low risk sets
  ii_left <- x$intervals$II$left[est$ly_remove]
  ii_right <- x$intervals$II$right[est$ly_remove]
  for (j in seq_along(ii_right)) {
    x_1 <- ii_left[j]
    x_2 <- ii_right[j]
    rect(
      x_1, 0, x_2, max_y,
      col = rgb(0.75, 0.75, 0.75, alpha = 0.5), border = NA
    )
  }

  #par(mar = c(0, 0, 0, 0), new = TRUE)
  legend(
    "topright", legend = c("Events", "Risk-set", "Small risk-set"),
    inset = 0.02,  xpd = TRUE,
    col = c("red", "black", rgb(0.75, 0.75, 0.75, alpha = 0.5)),
    fill = c(0, 0, rgb(0.75, 0.75, 0.75, alpha = 0.5)),
    border = c(NA, NA, NA),
    lty = c(1, 2, 0), cex = 0.8, box.lty = 0, bg = "white"
  )

}


lai_ying_plot(fem_comb, xlim = c(65, 100))
  
  