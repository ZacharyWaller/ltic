
# Load results -----------------------------------------------------------------
results_fem <- readRDS("outputs/results_fem.RDS")
results_mal <- readRDS("outputs/results_mal.RDS")

# Make plots -------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(xtable)

data <- results_fem
sim_name <- "mhcps_fem"

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
  theme_minimal()

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
  theme_minimal()

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
plot_estimate(prod_fem, surv_type = "prod", xlab = "Age", ylab = "Survival", main = "Females")
plot_estimate(opti_fem, surv_type = "step", col = "red", lty = 2, plot_type = "over")
plot_estimate(turn_fem, surv_type = "step", col = "blue", lty = 3, plot_type = "over")

plot_estimate(prod_mal, surv_type = "prod", xlab = "Age", ylab = "Survival", main = "Males")
plot_estimate(opti_mal, surv_type = "step", col = "red", lty = 2, plot_type = "over")
plot_estimate(turn_mal, surv_type = "step", col = "blue", lty = 3, plot_type = "over")

dev.off()

plot_estimate(prod, surv_type = "prod", cond = 1)
plot_estimate(opti, surv_type = "step", cond = 1, col = "red", lty = 2, plot_type = "over")
plot_estimate(turn, surv_type = "step", cond = 0, col = "blue", lty = 3, plot_type = "over")