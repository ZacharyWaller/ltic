# Load results -----------------------------------------------------------------
fem_results <- readRDS("outputs/mhcps_female.RDS")
mal_results <- readRDS("outputs/mhcps_male.RDS")

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(xtable)

# Make tables ------------------------------------------------------------------
methods <- c(
  "Breslow", "Product-Limit", "Shen", "Turnbull", "Yu",
  "Breslow-ICM", "PL-ICM", "Turnbull-ICM", "Yu-ICM", "ICM", "Quasi-Newton"
)

table_female <- lapply(
  methods,
  function(method) {
    est <- fem_results[[method]]
    data.frame(
      method = method,
      time = as.numeric(est$time, units = "secs"),
      it = est$res$it,
      llike = est$res$llike
    )
  }
)

table_male <- lapply(
  methods,
  function(method) {
    est <- mal_results[[method]]
    data.frame(
      method = method,
      time = as.numeric(est$time, units = "secs"),
      it = est$res$it,
      llike = est$res$llike
    )
  }
)

## for paper ----
paper_table_female <- bind_rows(table_female) %>%
  filter(
    !is.nan(llike)
  ) %>%
  mutate(
    time = formatC(time, format = "f", digits = 3),
    llike = formatC(llike, format = "f", digits = 2),
    it = formatC(it, big.mark = ","),
  ) %>%
  select(
    "Method" = method, "log-likelihood" = llike, "Iterations" = it,
    "CPU time (s)" = time
  )

paper_table_male <- bind_rows(table_male) %>%
  filter(
    !is.nan(llike)
  ) %>%
  mutate(
    time = formatC(time, format = "f", digits = 3),
    llike = formatC(llike, format = "f", digits = 2),
    it = formatC(it, big.mark = ","),
  ) %>%
  select(
    "Method" = method, "log-likelihood" = llike, "Iterations" = it,
    "CPU time (s)" = time
  )

paper_table_female
paper_table_male

## save ----
saveRDS(paper_table_female, paste0("outputs/mhcps_female_table.RDS"))
saveRDS(paper_table_male, paste0("outputs/mhcps_male_table.RDS"))

# Plots ------------------------------------------------------------------------
## Unconditional Plots ----
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

## Conditional likelihood ----
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

## Conditional vs Trubnull ----
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

## Males and Females together ----
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

## Lai-Ying Estimator ----
ests <- list(
  "Females" = fem_results$`PL-ICM`,
  "Males"   = mal_results$`PL-ICM`
)

ly_fem <- lai_ying(fem_results[["PL-ICM"]])
ly_mal <- lai_ying(mal_results[["PL-ICM"]])

plot(ly_fem, xlim = c(65, 100))
plot(ly_mal, xlim = c(65, 100))

dev.off()

## Turnbull augmented risk set ----
# This gives an idea of the size of the augmented data in Turnbull's aproach.
# The augmented data includes the number of observed values and the ghosts used 
# to make truncation work.

### Females ----
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

### Males ----
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

### Both ----
fem_aug
mal_aug
