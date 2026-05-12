library(dplyr)
library(readr)

# MHCPS is included in this package
data <- mhcps

# Split into female and male
data_fem <- data %>% filter(Zi == 0)
data_mal <- data %>% filter(Zi == 1)

max_it <- 1E6
tol <- 1e-8

methods <- c(
  "Product-Limit", "PL-ICM", "Breslow", "Breslow-ICM", "Turnbull",
  "Turnbull-ICM", "Yu", "Yu-ICM", "Shen", "ICM", "Quasi-Newton"
)

# Female -----------------------------------------------------------------------
left <- data_fem$Ui
right <- data_fem$Vi
trunc <- data_fem$Ti

## Run estimates ----
fem_results <- lapply(
  methods,
  function(method) {
    ltic_np(
      left, right, trunc, method = method, tol = tol, max_it = max_it
    )
  }
)
names(fem_results) <- methods

saveRDS(fem_results, "outputs/mhcps_female.RDS")

# Male -------------------------------------------------------------------------
## Run estimates ----
left <- data_mal$Ui
right <- data_mal$Vi
trunc <- data_mal$Ti

## Run estimates ----
mal_results <- lapply(
  methods,
  function(method) {
    ltic_np(
      left, right, trunc, method = method, tol = tol, max_it = max_it
    )
  }
)
names(mal_results) <- methods

saveRDS(mal_results, "outputs/mhcps_male.RDS")

# Unconditional Plots ----------------------------------------------------------
pdf(
  file = "outputs/mhcps_plots.pdf",
  width = 8, height = 4
)

par(mfrow = c(1, 2))
plot(
  fem_results[["Product-Limit"]], end = "r",
  xlab = "Age", ylab = "Survival probability", main = "Females"
)
plot(fem_results[["Quasi-Newton"]], col = "red", lty = 2, plot_type = "over")
plot(fem_results[["Turnbull"]], col = "blue", lty = 3, plot_type = "over")

plot(
  mal_results[["Product-Limit"]], end = "r",
  xlab = "Age", ylab = "Survival probability", main = "Males"
)
plot(mal_results[["Quasi-Newton"]], col = "red", lty = 2, plot_type = "over")
plot(mal_results[["Turnbull"]], col = "blue", lty = 3, plot_type = "over")

dev.off()

# Conditional likelihood -------------------------------------------------------
pdf(
  file = "outputs/mhcps_plots_cond.pdf",
  width = 8, height = 4
)

par(mfrow = c(1, 2))
plot(
  fem_results[["Product-Limit"]], end = "r", cond = 1,
  xlab = "Age", ylab = "Survival probability", main = "Females"
)
plot(
  fem_results[["Quasi-Newton"]], cond = 1, col = "red", lty = 2,
  plot_type = "over"
)
plot(
  fem_results[["Turnbull"]], cond = 1, col = "blue", lty = 3,
  plot_type = "over"
)

plot(
  mal_results[["Product-Limit"]], cond = 1, end = "r",
  xlab = "Age", ylab = "Survival probability", main = "Males"
)
plot(
  mal_results[["Quasi-Newton"]], cond = 1, col = "red", lty = 2,
  plot_type = "over"
)
plot(
  mal_results[["Turnbull"]], cond = 1, col = "blue", lty = 3,
  plot_type = "over"
)

dev.off()

# Conditional vs Trubnull ------------------------------------------------------
pdf(
  file = "outputs/mhcps_plots_cond_turnbull.pdf",
  width = 8, height = 4
)

par(mfrow = c(1, 2))
plot(
  fem_results[["Product-Limit"]], end = "r", cond = 1,
  xlab = "Age", ylab = "Survival probability", main = "Females"
)
plot(fem_results[["Turnbull"]], col = "red", lty = 2, plot_type = "over")

plot(
  mal_results[["Product-Limit"]], cond = 1, end = "r",
  xlab = "Age", ylab = "Survival probability", main = "Males"
)
plot(mal_results[["Turnbull"]], col = "red", lty = 2, plot_type = "over")

legend(
  "topright", legend = c("Conditional", "Turnbull"), inset = 0.02,  xpd = TRUE,
  col = c("black", "red"), lty = 1:2, cex = 0.8, box.lty = 0
)

dev.off()

# Males and Females together ---------------------------------------------------
pdf(
  file = "outputs/mhcps_plots_male_female.pdf",
  width = 4, height = 4
)

par(mar = c(4.1, 4.1, 2.1, 2.1))
plot(
  fem_results[["Product-Limit"]], cond = 1,
  xlab = "Age", ylab = "Survival", main = "Males and Females"
)
plot(
  mal_results[["Product-Limit"]], cond = 1,
  plot_type = "over", col = "red"
)

legend(
  "topright", legend = c("Females", "Males"), inset = 0.02,  xpd = TRUE,
  col = c("black", "red"), lty = 1, cex = 0.8, box.lty = 0
)

dev.off()


# Risk sets --------------------------------------------------------------------
pdf(
  file = "outputs/mhcps_risk_set.pdf",
  width = 8, height = 4
)

ests <- list(
  "Females" = results_fem$comb,
  "Males"   = results_mal$comb
)
ests <- list(
  "Females" = fem_comb,
  "Males"   = mal_comb
)
max_y <- 0

par(mfrow = c(1, 2), mar = c(4.1, 4.1, 2.1, 2.1))
for (i in seq_along(ests)) {

  est <- ests[[i]]
  max_y <- max(max_y, ceiling(max(est$y_0)))

  ii_right <- est$intervals$II$right
  r_time_points <- sort(c(0, rep(ii_right, each = 2), Inf))
  x <- r_time_points
  ord <- order(x)

  steps <- c(0, est$res$n)
  y <- rep(steps, each = 2)
  y <- y[ord]
  plot(
    x, y, type = "l", col = "red", lwd = 1.5,
    xlim = c(64, 100), ylim = c(0, max_y),
    xlab = "Age", ylab = "Expected Number", main = names(ests)[[i]]
  )

  steps <- c(0, est$y_0 - c(0, cumsum(est$res$n[-length(est$res$n)])))

  y <- rep(steps, each = 2)
  x <- x[ord]
  y <- y[ord]

  lines(x, y, type = "l", lwd = 1.5, lty = 3)

  # Lai-Yin stuff
  n <- nrow(est$intervals$Oi)
  low_lim <- n^0.25
  risk_set <- c(est$y_0 - c(0, cumsum(est$res$n[-length(est$res$n)])))
  lai_yin <- risk_set < low_lim
  ii_left <- est$intervals$II$left[lai_yin]
  ii_right <- est$intervals$II$right[lai_yin]
  for (j in 1:length(ii_right)) {
    x_1 <- ii_left[j]
    x_2 <- ii_right[j]
    rect(
      x_1, 0, x_2, max_y,
      col = rgb(0.75, 0.75, 0.75, alpha = 0.5), border = NA)
  }
}

par(mar = c(0, 0, 0, 0), new = TRUE)
legend(
  "topright", legend = c("Events", "Risk-set", "Small risk-set"),
  inset = 0.02,  xpd = TRUE,
  col = c("red", "black", rgb(0.75, 0.75, 0.75, alpha = 0.5)),
  fill = c(0, 0, rgb(0.75, 0.75, 0.75, alpha = 0.5)),
  border = c(NA, NA, NA),
  lty = c(1, 2, 0), cex = 0.8, box.lty = 0, bg = "white"
)

dev.off()

# Lai-Ying Estimator -----------------------------------------------------------
pdf(
  file = "outputs/mhcps_lai_ying.pdf",
  width = 8, height = 4
)

ests <- list(
  "Females" = fem_results$`PL-ICM`,
  "Males"   = mal_results$`PL-ICM`
)

ly_fem <- lai_ying(fem_results[["PL-ICM"]])
ly_mal <- lai_ying(mal_results[["PL-ICM"]])

plot(ly_fem, xlim = c(65, 100))
plot(ly_mal, xlim = c(65, 100))


par(mfrow = c(1, 2), mar = c(4.1, 4.1, 2.1, 2.1))
for (i in seq_along(ests)) {

  est <- ests[[i]]

  events <- est$res$n
  risk_set <- c(est$R_0 - c(0, cumsum(est$res$n[-length(est$res$n)])))

  n <- nrow(est$intervals$Oi)
  low_lim <- n^0.25
  lai_yin <- risk_set < low_lim
  events[lai_yin] <- 0
  h <- events / risk_set
  surv <- c(1, cumprod(1 - h))

  ii_right <- est$intervals$II$right
  r_time_points <- sort(c(0, rep(ii_right, each = 2), Inf))
  x <- r_time_points
  ord <- order(x)

  steps <- surv
  y <- rep(steps, each = 2)
  y <- y[ord]
  plot(
    x, y, type = "l", col = "red", lwd = 1.5,
    xlim = c(64, 100), ylim = c(0, 1),
    xlab = "Age", ylab = "Survival", main = names(ests)[[i]]
  )

  cond_events <- est$res$n
  cond_events[1] <- 0
  h <- cond_events / risk_set
  surv <- c(1, cumprod(1 - h))

  y <- rep(surv, each = 2)
  x <- x[ord]
  y <- y[ord]

  lines(x, y, type = "l", lwd = 1.5, lty = 2)
}

par(mar = c(0, 0, 0, 0), new = TRUE)
legend(
  "topright", legend = c("Conditional", "Lai-Ying"), inset = 0.02,  xpd = TRUE,
  col = c("black", "red"), lty = 2:1, cex = 0.8, box.lty = 0
)

dev.off()

# Turnbull augmented risk set --------------------------------------------------
# This gives an idea of the size of the augmented data in Turnbull's aproach.
# The augmented data includes the number of observed values and the ghosts used 
# to make truncation work.

## Females ----
int <- fem_results[["Turnbull"]]$intervals
alpha <- indicator_matrix(int$II, int$Oi)
beta <- indicator_matrix(int$II, int$Ti)

s <- fem_results[["Turnbull"]]$res$s

mu_ij <- t(t(alpha) * s) / colSums(t(alpha) * s)
nu_ij <- t(t(1 - beta) * s) / colSums(t(beta) * s)

# Number of observed values
sum(mu_ij)
# equal to the number of females in the data set (up to some tolerance)
abs(sum(mu_ij) - filter(mhcps, Zi == 0) %>% nrow()) < 1E-9

fem_aug <- scales::scientific(sum(mu_ij + nu_ij))
fem_aug

## Males ----
int <- mal_results[["Turnbull"]]$intervals
alpha <- indicator_matrix(int$II, int$Oi)
beta <- indicator_matrix(int$II, int$Ti)

s <- mal_results[["Turnbull"]]$res$s

mu_ij <- t(t(alpha) * s) / colSums(t(alpha) * s)
nu_ij <- t(t(1 - beta) * s) / colSums(t(beta) * s)

# Number of observed values
sum(mu_ij)
# equal to the number of females in the data set (up to some tolerance)
abs(sum(mu_ij) - filter(mhcps, Zi == 1) %>% nrow()) < 1E-9

mal_aug <- scales::scientific(sum(mu_ij + nu_ij))

## Both ----
fem_aug
mal_aug
