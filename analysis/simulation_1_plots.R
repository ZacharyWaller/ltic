
# Load results -----------------------------------------------------------------
results_100 <- readRDS("outputs/results_0.RDS")
results_075 <- readRDS("outputs/results_075.RDS")
results_050 <- readRDS("outputs/results_050.RDS")
results_025 <- readRDS("outputs/results_025.RDS")
results <- readRDS("outputs/results.RDS")
results_mal <- readRDS("outputs/results_mal.RDS")

# Make plots -------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(xtable)

data <- results_100
sim_name <- "sim_100"
 
# Extract and collate times ----------------------------------------------------
t_results <- lapply(
  data,
  function(x) {
    data.frame(
      prod = x$prod$time,
      icm  = x$icm$time,
      comb = x$comb$time,
      tu_c = x$tu_c$time,
      yu_c = x$yu_c$time,
      br_c = x$br_c$time,
      bres = x$bres$time,
      opti = x$opti$time,
      yu   = x$yu$time,
      shen = x$shen$time,
      turn = x$turn$time
    )
  }
)

# Extract and collate numits ---------------------------------------------------
numit_results <- lapply(
  data,
  function(x) {
    data.frame(
      prod = x$prod$res$it,
      icm  = x$icm$res$it,
      comb = 11 * x$comb$res$it,
      tu_c = 11 * x$tu_c$res$it,
      yu_c = 11 * x$yu_c$res$it,
      br_c = 11 * x$br_c$res$it,
      bres = x$bres$res$it,
      opti = x$opti$res$numit[[1]],
      yu   = x$yu$res$it,
      shen = x$shen$res$it,
      turn = x$turn$res$it
    )
  }
)

# Extract and collate likelihoods ----------------------------------------------
like_results <- lapply(
  data,
  function(x) {
    data.frame(
      prod = x$prod$res$llike,
      icm  = x$icm$res$llike,
      comb = x$comb$res$llike,
      tu_c = x$tu_c$res$llike,
      yu_c = x$yu_c$res$llike,
      br_c = x$br_c$res$llike,
      bres = x$bres$res$llike,
      opti = x$opti$res$like,
      yu   = x$yu$res$llike,
      shen = x$shen$res$llike,
      turn = x$turn$res$llike
    )
  }
)
# Plot helpers -----------------------------------------------------------------
labels <- c(
  "prod" = "Product limit",
  "icm"  = "ICM",
  "comb" = "Product limit + ICM",
  "tu_c" = "Turnbull + ICM",
  "yu_c" = "Yu  + ICM",
  "br_c" = "Breslow + ICM",
  "bres" = "Breslow",
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
  #pivot_longer(c(prod, comb, tu_c, yu_c, br_c, bres, yu, shen, turn)) %>%
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
  #pivot_longer(c(prod, comb, tu_c, yu_c, br_c, bres, yu, shen, turn)) %>%
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
    #ranks = list(rank(-c(prod, comb, tu_c, yu_c, br_c, bres, yu, shen, turn)))
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
  )

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
  "Av. Rank" = t(table_like_ranks),
  "log-likelihood" = paste0(like_mean, "(", like_sd, ")"),
  "Iterations" = paste0(numit_mean, "(", numit_sd, ")"),
  "CPU time (s)" = paste0(time_mean, "(", time_sd, ")"),
  "% converged" = paste0(conv_per, "%"),
  check.names = FALSE
) %>%
  mutate(sim = sim_name) %>%
  rownames_to_column(var = "Method")

# Save -------------------------------------------------------------------------
saveRDS(res_table, paste0("outputs/", sim_name, ".RDS"))

# Collate tables together ------------------------------------------------------
table_1 <- readRDS("outputs/sim_100.RDS")
table_2 <- readRDS("outputs/sim_075.RDS")
table_3 <- readRDS("outputs/sim_050.RDS")

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

table_3 %>%
  arrange(
    desc(sim), Method
  ) %>%
  select(-sim) %>%
  xtable() %>%
  print(include.rownames = FALSE)


res_table %>%
  arrange(
    desc(sim), Method
  ) %>%
  select(-sim) %>%
  xtable() %>%
  print(include.rownames = FALSE)

res_table %>%
  arrange(
    desc(sim), Method
  ) %>%
  select(-sim) %>%
  xtable() %>%
  print(include.rownames = FALSE)

