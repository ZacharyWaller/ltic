
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
  "bres" = "Breslow",
  "br_c" = "Breslow-ICM",
  "icm"  = "ICM",
  "prod" = "PL",
  "comb" = "PL-ICM",
  "opti" = "Quasi-Newton",
  "shen" = "Shen",
  "turn" = "Turnbull",
  "tu_c" = "Turnbull-ICM",
  "yu"   = "Yu",
  "yu_c" = "Yu-ICM"
)

levels <- names(labels)

# Plots ------------------------------------------------------------------------
## Times ----
plot_time <- bind_rows(t_results) %>%
  mutate(run = row_number()) %>%
  pivot_longer(c(bres, br_c, icm, prod, comb, opti, shen, turn, tu_c, yu, yu_c)) %>%
  mutate(value = as.numeric(value)) %>%
  ggplot(aes(name, value)) +
  scale_y_log10(
    name = "CPU Time (s)", limits = c(0.003, 100)
  ) +
  scale_x_discrete(
    name = "Method", labels = labels, limits = rev(levels)
  ) +
  geom_boxplot() +
  theme_minimal() +
  coord_flip()

saveRDS(plot_time, paste0("outputs/", sim_name, "_plot.RDS"))

## Iterations ----
plot_numit <- bind_rows(numit_results) %>%
  mutate(run = row_number()) %>%
  pivot_longer(c(bres, br_c, icm, prod, comb, opti, shen, turn, tu_c, yu, yu_c)) %>%
  ggplot(aes(name, value)) +
  scale_y_log10(
    name = "Iterations"
  ) +
  scale_x_discrete(
    name = "Method", labels = labels, limits = rev(levels)
  ) +
  geom_boxplot() +
  theme_minimal() +
  coord_flip()

# Tables -----------------------------------------------------------------------
## Likeilhood Ranks ----
like_all <- bind_rows(like_results) %>%
  mutate_all(
    function(x) if_else(is.nan(x), -Inf, x)
  )

table_like_ranks <- like_all %>%
  rowwise() %>%
  mutate(
    ranks = list(rank(-c(bres, br_c, icm, prod, comb, opti, shen, turn, tu_c, yu, yu_c)))
  ) %>%
  select(ranks) %>%
  unnest_wider(ranks, names_sep = "_") %>%
  summarise(
    across(everything(),
      function(x) formatC(median(x), format = "f", digits = 1, big.mark = " ")
    )
  ) %>%
  setNames(labels)

## Likelihood ----
like_mean <- like_all %>%
  summarise(
    across(everything(),
      function(x) formatC(median(x), format = "f", digits = 2, big.mark = " ")
    )
  )

like_sd <- like_all %>%
  mutate_all(
    function(x) if_else(is.nan(x), -Inf, x)
  ) %>%
  summarise(
    across(everything(),
      function(x) formatC(IQR(x), format = "f", digits = 2, big.mark = " ")
    )
  )

## Iterations ----
numit_mean <- bind_rows(numit_results) %>%
  summarise(
    across(everything(),
      function(x) formatC(median(x), format = "f", digits = 0, big.mark = " ")
    )
  )
numit_sd <- bind_rows(numit_results) %>%
  summarise(
    across(everything(),
      function(x) formatC(IQR(x), format = "f", digits = 0, big.mark = " ") 
    )
  )

## % Converged
finite_like <- is.finite(as.matrix(like_all))
identical(conv_per, test2)
not_max_it <- bind_rows(numit_results) %>%
  mutate(across(everything(), function(x) x < 1e5)) %>% 
  as.matrix()

conv_per <- (not_max_it * finite_like) %>% colSums()

## Times ----
time_mean <- bind_rows(t_results) %>%
  summarise(
    across(everything(),
      function(x) formatC(median(x), format = "f", digits = 3, big.mark = " ")
    )
  )
time_sd <- bind_rows(t_results) %>%
  summarise(
    across(everything(),
      function(x) formatC(IQR(x), format = "f", digits = 3, big.mark = " ")
    )
  )

order <- levels

res_table <- data.frame(
  #"Av. Rank" = t(table_like_ranks),
  "Method" = labels[order],
  "log-likelihood" = paste0(like_mean[order], "(", like_sd[order], ")"),
  "Iterations" = paste0(numit_mean[order], "(", numit_sd[order], ")"),
  "CPU time (s)" = paste0(time_mean[order], "(", time_sd[order], ")"),
  "% converged" = paste0(conv_per[order], "%"),
  check.names = FALSE
) %>%
  mutate(sim = sim_name)

# Save -------------------------------------------------------------------------
saveRDS(res_table, paste0("outputs/", sim_name, ".RDS"))

# Collate tables together ------------------------------------------------------
files <- list.files("outputs/sim", full.names = TRUE)
methods <- c(
  "prod", "icm", "comb", "tu_c", "yu_c", "br_c", "bres", "opti",
  "yu", "shen", "turn"
)

for (file in files) {
  cat(file, "\n")
  data <- readRDS(file)
  name <- gsub("(.*/)|(\\..*)", "", file)

  # Extract all relevant information from each run
  table <- lapply(
    data,
    function(run) {
      lapply(
        methods,
        function(method) {
          est <- run[[method]]
          data.frame(
            method = method,
            time = as.numeric(est$time, units = "secs"),
            it = est$res$it[[1]],
            llike = est$res$llike
          )
        }
      ) %>%
        bind_rows()
    }
  )

  # combine into a nice table
  paper_table <- bind_rows(table) %>%
    mutate(
      converged = if_else(!is.nan(llike), 1, 0)
    ) %>%
    group_by(method) %>%
    summarise(
      time_med = mean(time), time_iqr = IQR(time),
      it_med = median(it), it_iqr = IQR(it),
      llike_med = median(llike, na.rm = TRUE),
      llike_iqr = IQR(llike, na.rm = TRUE),
      converged = mean(converged)
    )

  saveRDS(paper_table, paste0("outputs/sim/", name, "_table.RDS"))
}

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

# Plots ----
plot_075 <- readRDS("outputs/sim_075_plot.RDS")
plot_050 <- readRDS("outputs/sim_050_plot.RDS")

library(cowplot)
time_plot <- plot_grid(
  plot_075, plot_050,
  labels = c("a) 25% truncation", "b) 50% truncation"),
  label_size = 10, hjust = -0.1
)

ggsave(
  "outputs/time.pdf", time_plot, width = 10, height = 5
)


