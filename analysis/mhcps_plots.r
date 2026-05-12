
# Load results -----------------------------------------------------------------
results_fem <- readRDS("outputs/results_fem.RDS")
results_mal <- readRDS("outputs/results_mal.RDS")

# Make plots -------------------------------------------------------------------
devtools::load_all()
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(xtable)

data <- results_mal
sim_name <- "mhcps_mal"

# Extract and collate times ----------------------------------------------------
t_results <- data.frame(
  prod = data$prod$time,
  icm  = data$icm$time,
  comb = data$comb$time,
  tu_c = data$tu_c$time,
  yu_c = data$yu_c$time,
  br_c = data$br_c$time,
  bres = data$bres$time,
  opti = data$opti$time,
  yu   = data$yu$time,
  shen = data$shen$time,
  turn = data$turn$time
)

# Extract and collate numits ---------------------------------------------------
numit_results <- data.frame(
  prod = data$prod$res$it,
  icm  = data$icm$res$it,
  comb = 11 * data$comb$res$it,
  tu_c = 11 * data$tu_c$res$it,
  yu_c = 11 * data$yu_c$res$it,
  br_c = 11 * data$br_c$res$it,
  bres = data$bres$res$it,
  opti = data$opti$res$numit[[1]],
  yu   = data$yu$res$it,
  shen = data$shen$res$it,
  turn = data$turn$res$it
)

# Extract and collate likelihoods ----------------------------------------------
like_results <- data.frame(
  prod = data$prod$res$llike,
  icm  = data$icm$res$llike,
  comb = data$comb$res$llike,
  tu_c = data$tu_c$res$llike,
  yu_c = data$yu_c$res$llike,
  br_c = data$br_c$res$llike,
  bres = data$bres$res$llike,
  opti = data$opti$res$like,
  yu   = data$yu$res$llike,
  shen = data$shen$res$llike,
  turn = data$turn$res$llike
)

# Plot helpers -----------------------------------------------------------------
labels <- c(
  "prod" = "Product limit",
  "icm"  = "ICM",
  "comb" = "Product limit + ICM",
  "tu_c" = "Turnbull + ICM",
  "yu_c" = "Yu  + ICM",
  "br_c" = "Breslow + ICM",
  "bres" = "Breslow style",
  "opti" = "Quasi-Newton",
  "yu"   = "Yu",
  "shen" = "Shen",
  "turn" = "Turnbull"
)

# Plots ------------------------------------------------------------------------
## Times ----
plot_time <- bind_rows(t_results) %>%
  mutate(run = row_number()) %>%
  pivot_longer(c(prod, icm, comb, tu_c, yu_c, br_c, bres, opti, yu, shen, turn)) %>%
  mutate(value = as.numeric(value)) %>%
  ggplot(aes(name, value)) +
  scale_y_log10(
    name = "CPU Time (s)"
  ) +
  scale_x_discrete(
    name = "Method", labels = labels
  ) +
  geom_boxplot() +
  theme_minimal() +
  coord_flip()

## Iterations ----
plot_numit <- bind_rows(numit_results) %>%
  mutate(run = row_number()) %>%
  pivot_longer(c(prod, icm, comb, tu_c, yu_c, br_c, bres, opti, yu, shen, turn)) %>%
  ggplot(aes(name, value)) +
  scale_y_log10(
    name = "Iterations"
  ) +
  scale_x_discrete(
    name = "Method", labels = labels
  ) +
  geom_boxplot() +
  theme_minimal() +
  coord_flip()

# Tables -----------------------------------------------------------------------
## Likeilhood Ranks ----
table_like_ranks <- bind_rows(like_results) %>%
  rowwise() %>%
  mutate(
    ranks = list(rank(-c(prod, icm, comb, tu_c, yu_c, br_c, bres, opti, yu, shen, turn)))
  ) %>%
  select(ranks) %>%
  unnest_wider(ranks, names_sep = "_") %>%
  summarise(
    across(everything(),
      function(x) formatC(mean(x), format = "f", digits = 1, big.mark = " ")
    )
  ) %>%
  setNames(labels)

## Likelihood ----
like_mean <- bind_rows(like_results) %>%
  summarise(
    across(everything(),
      function(x) formatC(mean(x, na.rm = TRUE), format = "f", digits = 2, big.mark = " ")
    )
  ) %>%
  setNames(labels)

like_sd <- bind_rows(like_results) %>%
  summarise(
    across(everything(),
      function(x) formatC(sd(x, na.rm = TRUE), format = "f", digits = 2, big.mark = " ")
    )
  )

## Iterations ----
numit_mean <- bind_rows(numit_results) %>%
  summarise(
    across(everything(),
      function(x) formatC(mean(x), format = "f", digits = 0, big.mark = " ")
    )
  )
numit_sd <- bind_rows(numit_results) %>%
  summarise(
    across(everything(),
      function(x) formatC(sd(x, na.rm = TRUE), format = "f", digits = 0, big.mark = " ")
    )
  )

## % Converged
conv_per <- bind_rows(numit_results) %>%
  mutate(across(everything(), function(x) x < 1e5 & !is.nan(x))) %>%
  summarise(across(everything(), function(x) 100 * mean(x)))

## Times ----
time_mean <- bind_rows(t_results) %>%
  summarise(
    across(everything(),
      function(x) formatC(mean(x), format = "f", digits = 3, big.mark = " ")
    )
  )
time_sd <- bind_rows(t_results) %>%
  summarise(
    across(everything(),
      function(x) formatC(sd(x), format = "f", digits = 3, big.mark = " ")
    )
  )

res_table <- data.frame(
  "log-likelihood" = t(like_mean),
  "Ierations" = t(numit_mean),
  "CPU time (s)" = t(time_mean),
  check.names = FALSE
) %>%
  mutate(sim = sim_name) %>%
  rownames_to_column(var = "Method")

# Save -------------------------------------------------------------------------
saveRDS(res_table, paste0("outputs/", sim_name, ".RDS"))

# Collate tables together ------------------------------------------------------
table_1 <- readRDS("outputs/mhcps_mal.RDS")
table_2 <- readRDS("outputs/mhcps_fem.RDS")

table_1 %>%
  arrange(
    desc(sim), Method
  ) %>%
  select(-sim) %>%
  xtable() %>%
  print(include.rownames = FALSE)

table_2 %>%
  arrange(
    desc(sim), Method
  ) %>%
  select(-sim) %>%
  xtable() %>%
  print(include.rownames = FALSE)

# Unconitional -----------------------------------------------------------------
pdf(
  file = "outputs/mhcps_plots.pdf",
  width = 8, height = 4
)
par(mfrow = c(1, 2))
par(mar = c(4.1, 4.1, 2.1, 2.1))

plot.ltic(results_fem$prod, xlab = "Age", ylab = "Survival", main = "Females")
plot.ltic(results_fem$opti, col = "red", lty = 2, plot_type = "over")
plot.ltic(results_fem$turn, col = "blue", lty = 3, plot_type = "over")

plot_estimate(results_mal$prod, surv_type = "prod", xlab = "Age", ylab = "Survival", main = "Males")
plot_estimate(results_mal$opti, surv_type = "step", col = "red", lty = 2, plot_type = "over")
plot_estimate(results_mal$turn, surv_type = "step", col = "blue", lty = 3, plot_type = "over")

legend(
  "topright", legend = c("PL-ICM", "Quasi-Newton", "Turnbull"), inset = 0.02,  xpd = TRUE,
  col = c("black", "red", "blue"), lty = 1:3, cex = 0.8, box.lty = 0
)

dev.off()

# Conditional ------------------------------------------------------------------
pdf(
  file = "outputs/mhcps_plots_cond.pdf",
  width = 8, height = 4
)
par(mfrow = c(1, 2))
par(mar = c(4.1, 4.1, 2.1, 2.1))

plot_estimate(results_fem$prod, surv_type = "prod", cond = 1, xlab = "Age", ylab = "Survival", main = "Females")
plot_estimate(results_mal$prod, surv_type = "prod", cond = 1, xlab = "Age", ylab = "Survival", main = "Males")

dev.off()

# Conditional vs Turnbull ------------------------------------------------------
pdf(
  file = "outputs/mhcps_plots_cond_turnbull.pdf",
  width = 8, height = 4
)

par(mfrow = c(1, 2))
par(mar = c(4.1, 4.1, 2.1, 2.1))

plot_estimate(results_fem$prod, surv_type = "prod", cond = 1, xlab = "Age", ylab = "Survival", main = "Females")
plot_estimate(results_fem$turn, surv_type = "step", col = "red", lty = 2, plot_type = "over")

plot_estimate(results_mal$prod, surv_type = "prod", cond = 1, xlab = "Age", ylab = "Survival", main = "Males")
plot_estimate(results_mal$turn, surv_type = "step", col = "red", lty = 2, plot_type = "over")

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

#par(mfrow = c(1, 2), mar = c(4.1, 4.1, 2.1, 2.1))
par(mar = c(4.1, 4.1, 2.1, 2.1))
plot_estimate(
  results_fem$prod, surv_type = "prod", cond = 1,
  xlab = "Age", ylab = "Survival", main = "Males and Females"
)
plot_estimate(
  results_mal$prod, surv_type = "prod", cond = 1,
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

  events <- est$res$n
  risk_set <- c(est$y_0 - c(0, cumsum(est$res$n[-length(est$res$n)])))

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

# Turnbull augmented risk set ----
## Females ----
data <- results_fem
int <- data$turn$intervals
alpha <- indicator_matrix(int$II, int$Oi)
beta <- indicator_matrix(int$II, int$Ti)

s <- data$turn$res$s

mu_ij <- t(t(alpha) * s) / colSums(t(alpha) * s)
nu_ij <- t(t(1 - beta) * s) / colSums(t(beta) * s)

fem_aug <- scales::scientific(sum(mu_ij + nu_ij))
1 - sum(nu_ij) / sum(mu_ij + nu_ij)

## Males ----
data <- results_mal
int <- data$turn$intervals
alpha <- indicator_matrix(int$II, int$Oi)
beta <- indicator_matrix(int$II, int$Ti)

s <- data$turn$res$s

mu_ij <- t(t(alpha) * s) / colSums(t(alpha) * s)
nu_ij <- t(t(1 - beta) * s) / colSums(t(beta) * s)

mal_aug <- scales::scientific(sum(mu_ij + nu_ij))
1 - sum(nu_ij) / sum(mu_ij + nu_ij)

## Both ----
fem_aug
mal_aug

