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
turn_comb <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "turn_comb",
                     tol = tol, max_it = 100)
prod <- ltic_np(left, right, trunc, remove_rcens = TRUE, method = "prodlim",
                tol = tol, max_it = max_it)
bres <- ltic_np(left, right, trunc, remove_rcens = TRUE, method = "breslow",
                tol = tol, max_it = max_it)
comb <- ltic_np(left, right, trunc, remove_rcens = TRUE, method = "both",
                tol = tol, max_it = max_it)
icm <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "icm",
               tol = tol, max_it = max_it)
yu_comb <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "yu_comb",
                   tol = tol, max_it = 500)
shen <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "shen",
                tol = tol, max_it = max_it)
opti <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "optim",
                tol = tol, max_it = max_it)
yu <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "yu",
              tol = tol, max_it = max_it)
turn <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "turnbull",
                tol = tol, max_it = max_it)
vardi <- lb_np(left, right, tol = tol, max_it = max_it)


results_fem <- list(
  prod = c(prod),
  icm = c(icm),
  comb = c(comb),
  tu_c = c(turn_comb),
  yu_c = c(yu_comb),
  bres = c(bres),
  opti = c(opti),
  yu = c(yu),
  shen = c(shen),
  turn = c(turn)
)

saveRDS(results_fem, "outputs/results_fem.RDS")

# Male -------------------------------------------------------------------------
## Run estimates ----
left <- data_mal$Ui
right <- data_mal$Vi
trunc <- data_mal$Ti

## Run estimates ----
turn_comb <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "turn_comb",
                tol = tol, max_it = 100)
prod <- ltic_np(left, right, trunc, remove_rcens = TRUE, method = "prodlim",
                tol = tol, max_it = max_it)
bres <- ltic_np(left, right, trunc, remove_rcens = TRUE, method = "breslow",
                tol = tol, max_it = max_it)
comb <- ltic_np(left, right, trunc, remove_rcens = TRUE, method = "both",
                tol = tol, max_it = max_it)
icm <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "icm",
                tol = tol, max_it = max_it)
yu_comb <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "yu_comb",
                tol = tol, max_it = 500)
shen <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "shen",
                tol = tol, max_it = max_it)
opti <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "optim",
                tol = tol, max_it = max_it)
yu <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "yu",
              tol = tol, max_it = max_it)
turn <- ltic_np(left, right, trunc, remove_rcens = FALSE, method = "turnbull",
                tol = tol, max_it = max_it)

results_mal <- list(
  prod = c(prod),
  icm = c(icm),
  comb = c(comb),
  tu_c = c(turn_comb),
  yu_c = c(yu_comb),
  bres = c(bres),
  opti = c(opti),
  yu = c(yu),
  shen = c(shen),
  turn = c(turn)
)

saveRDS(results_mal, "outputs/results_mal.RDS")

# Unconditional Plots ----------------------------------------------------------
pdf(
  file = "outputs/mhcps_plots.pdf",
  width = 8, height = 4
)

par(mfrow = c(1, 2))
plot_estimate(fem_opt, lwd = 2, xlim = c(65, 100), xlab = "Age", ylab = "Survival probability", main = "Females")
plot_estimate(fem_prod, lwd = 2, col = "red", plot_over = TRUE)
plot_estimate(fem_turn, lwd = 2, col = "blue", plot_over = TRUE)

plot_estimate(mal_opt, lwd = 2, xlim = c(65, 100), xlab = "Age", ylab = "Survival probability", main = "Males")
plot_estimate(mal_prod, lwd = 2, col = "red", plot_over = TRUE)
plot_estimate(mal_turn, lwd = 2, col = "blue", plot_over = TRUE)

dev.off()

# Conditional likelihood -------------------------------------------------------
data_fem_cond <- data %>% filter(Ti > 65.1, Zi == 0)
data_mal_cond <- data %>% filter(Ti > 65.1, Zi == 1)

fem_prod_cond <- sc_prodlim_proper(data_fem_cond$Ui, data_fem_cond$Vi, data_fem_cond$Ti, tol = 1e-9)
fem_prod_cond$numit
plot_estimate(fem_prod_cond, lwd = 3)
fem_prod_cond$like

fem_turn_cond <- sc_turnbull_method(data_fem_cond$Ui, data_fem_cond$Vi, data_fem_cond$Ti, tol = 1e-6)
fem_turn_cond$numit
plot_estimate(fem_turn_cond, lwd = 3)
fem_turn_cond$like

t0 <- Sys.time()
mal_prod_cond <- sc_prodlim_proper(data_mal_cond$Ui, data_mal_cond$Vi, data_mal_cond$Ti, tol = 1e-9)
Sys.time() - t0
mal_prod_cond$numit
plot_estimate(mal_prod_cond, lwd = 3)
mal_prod_cond$like

## Optimization using generic algorithm ----
ints <- inner_intervals(data_fem_cond$Ui, data_fem_cond$Vi, data_fem_cond$Ti)
alpha <- indicator_matrix(ints$II, ints$Oi)
beta <- indicator_matrix(ints$II, ints$Ti)
n <- ncol(alpha)
tol <- 1e-10
t0 <- Sys.time()
fem_opt <- optim(
  par = rep(1/n, n),
  fn = optim_likelihood,
  gr = like_gradient,
  lower = tol, upper = 1 - tol,
  alpha = alpha, beta = beta
)
Sys.time() - t0
-fem_opt$value
fem_prod$like

opt_list <- list(intervals = ints, surv = 1 - c(0, cumsum(fem_opt$par / sum(fem_opt$par))))
plot_estimate(opt_list, lwd = 2)
plot_estimate(fem_prod_cond, lwd = 2, col = "red", plot_over = TRUE)

fem_prod_cond$surv - opt_list$surv


fem_prod_cond$like - -fem_opt$value

data %>%
  filter(Zi == 0) %>%
  mutate(
    Vi = if_else(Vi == Inf, 9999, Vi)
  ) %>%
  select(-Zi) %>%
  write_delim("../npmle_c/mhcps_male.ssv", col_names = FALSE)


# Bootstraps -------------------------------------------------------------------
n_rep <- 100

t0 <- Sys.time()
results <- replicate(
  n_rep,
  {
    boot_data <- dplyr::sample_frac(data_fem_cond, replace = TRUE)
    boot_est <- sc_turnbull_method(boot_data$Ui, boot_data$Vi, boot_data$Ti)

    boot_est
  },
  simplify = FALSE
)
Sys.time() - t0


boot_results <- lapply(
  results,
  function(result) {
    surv <- result$surv
    intervals <- result$intervals$II

    list(
      x_values = sort(c(0, intervals$left, intervals$right, Inf)),
      surv = rep(surv, each = 2)
    )
  }
)

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

source("analysis/point_wise_average.R")
source("analysis/linear_interpolate.R")


times <- lapply(boot_results, function(x) x$x_values)
ests_surv <- lapply(boot_results, function(x) x$surv)
par_function <- function(times) rep(NA, length(times))
# Total survival
plot_surv <- point_wise_mean(
  times, ests_surv, par_function, plot_y_max = 1,
  plot_title = "Survival with functional independence", y_title = "Survival"
) + xlim(60, 100) + labs(x = "Age")
plot_surv
ggsave("../icnpmle_analysis/images/mhcps.pdf", plot_surv, "pdf", width = 10, height = 4)

# Digammma Haz -----------------------------------------------------------------
digamma_haz <- function(y, n) {
  digamma(y + 1) - digamma(y - n + 1)
}

ine_haz <- function(y, n) {
  n / y
}

plot(product_survival(fem_prod$lambda), type = "l")
lines(exp(-c(0, cumsum(digamma_haz(fem_prod$y, fem_prod$n)))), col = "red")
lines(exp(-c(0, cumsum(ine_haz(fem_prod$y, fem_prod$n)))), col = "blue")


npmle_haz <- -log(1 - fem_prod$lambda)
x <- which(!is.nan(npmle_haz) & !is.infinite(npmle_haz))
npmle_haz[which(is.nan(npmle_haz) | is.infinite(npmle_haz))] <- 1
lines(cumsum(npmle_haz), col = "red")




pan_data <- read_delim(
  "../npmle_c/mhcps_male.ssvNPMLERes0.000010",
  delim = " ", col_names = c("L", "R", "S")
)

plot_estimate(fem_turn, xlim = c(60, 100))
plot(c(pan_data$L, 99), 1 - c(0, cumsum(pan_data$S)))



pan_data$L
mal_prod$intervals$II$left

# Conditional plots ------------------------------------------------------------
pdf(
  file = "outputs/mhcps_cond_plots.pdf",
  width = 8, height = 4
)

par(mfrow = c(1, 3))
fem_cond <- list(
  intervals = list(II = fem_prod$intervals$II[-1, ]),
  surv = product_survival(fem_prod$lambda[-1])
)
plot_estimate(fem_cond, lwd = 2, xlim = c(65, 100), xlab = "Age", ylab = "Survival probability", main = "Females", col = "red")

mal_cond <- list(
  intervals = list(II = mal_prod$intervals$II[-1, ]),
  surv = product_survival(mal_prod$lambda[-1])
)
plot_estimate(mal_cond, lwd = 2, xlim = c(65, 100), xlab = "Age", ylab = "Survival probability", main = "Males", col = "red")

dev.off()

plot_estimate(fem_opt_cond, lwd = 2, xlim = c(65, 100), xlab = "Age", ylab = "Survival probability", main = "Females", col = "red")



fem_opt_cond <- list(
  intervals = list(II = fem_opt$intervals$II[-1, ]),
  surv = fem_opt$surv[-c(1)] / fem_opt$surv[2]
)


plot_estimate(fem_turn, lwd = 2, xlim = c(65, 100), xlab = "Age", ylab = "Survival probability", main = "Males", col = "red")


plot(cumprod(c(1, 1 - fem_prod$lambda)))
lines()


prodlim_rem <- ltic_np(data_fem$Ui, data_fem$Vi, data_fem$Ti, remove_rcens = TRUE,
                  method = "prodlim")

prodlim <- ltic_np(data_fem$Ui, data_fem$Vi, data_fem$Ti, remove_rcens = FALSE,
                  method = "prodlim")

prodlim$res$llike
prodlim_rem$res$llike

plot(exp(-diff(c(0, prodlim$res$lambda))))
points(exp(-diff(c(0, prodlim_rem$res$lambda))), col = "red")
