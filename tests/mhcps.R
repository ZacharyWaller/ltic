devtools::load_all()
library(dplyr)
library(readr)

data <- read_delim("../icprodlim/data/mhcps.txt", delim = " ")
data$Vi[data$Vi == 9999] <- Inf
data$Ui[data$Ui == 7755] <- 77.55
data$Ui[data$Ui == 873] <- 87.3
data$Ui[data$Ui == 1.55] <- 71.55
data$Ui[data$Ui == 3.30] <- 83.3
data$Ti[data$Ti == 689] <- 68.9

data_fem <- data %>% filter(Zi == 0)
data_mal <- data %>% filter(Zi == 1)

max_it <- 1e6
tol <- 1e-8

# Female -----------------------------------------------------------------------
left <- data_fem$Ui
right <- data_fem$Vi
trunc <- data_fem$Ti

## Run estimates ----
fem_turn_comb <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "turn_comb",
                    tol = tol, max_it = 500)
fem_prod <- ltic_np(left, right, trunc, remove_rcens = TRUE, method = "prodlim",
                tol = tol, max_it = max_it)
fem_bres <- ltic_np(left, right, trunc, remove_rcens = TRUE, method = "breslow",
              tol = tol, max_it = max_it)
fem_bres_comb <- ltic_np(left, right, trunc, remove_rcens = TRUE, method = "bres_comb",
                    tol = tol, max_it = max_it)
fem_comb <- ltic_np(left, right, trunc, remove_rcens = TRUE, method = "both",
                tol = tol, max_it = max_it)
fem_icm <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "icm",
              tol = tol, max_it = max_it)
fem_yu_comb <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "yu_comb",
                  tol = tol, max_it = 100)
fem_shen <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "shen",
                tol = tol, max_it = max_it)
fem_opti <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "optim",
                tol = tol, max_it = max_it)
fem_yu <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "yu",
              tol = tol, max_it = max_it)
fem_turn <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "turnbull",
                tol = tol, max_it = max_it)

results_fem <- list(
  prod = c(fem_prod),
  icm = c(fem_icm),
  comb = c(fem_comb),
  tu_c = c(fem_turn_comb),
  yu_c = c(fem_yu_comb),
  bres = c(fem_bres),
  opti = c(fem_opti),
  yu = c(fem_yu),
  shen = c(fem_shen),
  turn = c(fem_turn)
)

saveRDS(results_fem, "outputs/results_fem.RDS")

# Male -------------------------------------------------------------------------
## Run estimates ----
left <- data_mal$Ui
right <- data_mal$Vi
trunc <- data_mal$Ti

## Run estimates ----
mal_turn_comb <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "turn_comb",
                tol = tol, max_it = 100)
mal_prod <- ltic_np(left, right, trunc, remove_rcens = TRUE, method = "prodlim",
                tol = tol, max_it = max_it)
mal_bres <- ltic_np(left, right, trunc, remove_rcens = TRUE, method = "breslow",
                tol = tol, max_it = max_it)
mal_comb <- ltic_np(left, right, trunc, remove_rcens = TRUE, method = "both",
                tol = tol, max_it = max_it)
mal_icm <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "icm",
                tol = tol, max_it = max_it)
mal_yu_comb <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "yu_comb",
                tol = tol, max_it = 500)
mal_shen <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "shen",
                tol = tol, max_it = max_it)
mal_opti <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "optim",
                tol = tol, max_it = max_it)
mal_yu <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "yu",
              tol = tol, max_it = max_it)
mal_turn <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "turnbull",
                tol = tol, max_it = max_it)

results_mal <- list(
  prod = c(mal_prod),
  icm = c(mal_icm),
  comb = c(mal_comb),
  tu_c = c(mal_turn_comb),
  yu_c = c(mal_yu_comb),
  bres = c(mal_bres),
  opti = c(mal_opti),
  yu = c(mal_yu),
  shen = c(mal_shen),
  turn = c(mal_turn)
)

saveRDS(results_mal, "outputs/results_mal.RDS")

# Unconditional Plots ----------------------------------------------------------
pdf(
  file = "outputs/mhcps_plots.pdf",
  width = 8, height = 4
)

par(mfrow = c(1, 2))
plot_estimate(
  fem_prod, surv_type = "prod", end = "r",
  xlab = "Age", ylab = "Survival probability", main = "Females"
)
plot_estimate(fem_opti, surv_type = "step", col = "red", lty = 2, plot_type = "over")
plot_estimate(fem_turn, surv_type = "step", col = "blue", lty = 3, plot_type = "over")

plot_estimate(
  mal_prod, surv_type = "prod", end = "r",
  xlab = "Age", ylab = "Survival probability", main = "Males"
)
plot_estimate(mal_opti, surv_type = "step", col = "red", lty = 2, plot_type = "over")
plot_estimate(mal_turn, surv_type = "step", col = "blue", lty = 3, plot_type = "over")

dev.off()

# Conditional likelihood -------------------------------------------------------
pdf(
  file = "outputs/mhcps_plots_cond.pdf",
  width = 8, height = 4
)

par(mfrow = c(1, 2))
plot_estimate(
  fem_prod, surv_type = "prod", end = "r", cond = 1,
  xlab = "Age", ylab = "Survival probability", main = "Females"
)
plot_estimate(fem_opti, cond = 1, surv_type = "step", col = "red", lty = 2, plot_type = "over")
plot_estimate(fem_turn, cond = 1, surv_type = "step", col = "blue", lty = 3, plot_type = "over")

plot_estimate(
  mal_prod, cond = 1, surv_type = "prod", end = "r",
  xlab = "Age", ylab = "Survival probability", main = "Males"
)
plot_estimate(mal_opti, cond = 1, surv_type = "step", col = "red", lty = 2, plot_type = "over")
plot_estimate(mal_turn, cond = 1, surv_type = "step", col = "blue", lty = 3, plot_type = "over")

dev.off()
