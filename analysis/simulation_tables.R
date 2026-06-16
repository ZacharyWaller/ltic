# Make plots -------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(xtable)

# Collate tables together ------------------------------------------------------
files <- list.files("outputs/res", full.names = TRUE)

methods <- c(
  "Breslow", "Product-Limit", "Shen", "Turnbull", "Yu",
  "Breslow-ICM", "PL-ICM", "Turnbull-ICM", "Yu-ICM", "ICM", "Quasi-Newton"
)


for (file in files) {
  cat(file, "\n")
  data <- readRDS(file)
  name <- gsub("(.*/)|(\\..*)", "", file)

  data <- lapply(
    data,
    function(x) {
      names(x) <- methods
      x
    }
  )

  # Extract all relevant information from each run
  table <- lapply(
    data,
    function(run) {
      lapply(
        methods,
        function(method) {
          cat(method, "\\n")
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
      converged = if_else(!is.nan(llike), 1, 0),
      llike = if_else(is.nan(llike), -Inf, llike)
    ) %>%
    group_by(method) %>%
    summarise(
      time_med = median(time * converged),
      time_iqr = IQR(time * converged),
      it_med = median(it * converged),
      it_iqr = IQR(it * converged),
      llike_med = median(llike),
      llike_iqr = IQR(llike),
      converged = mean(converged)
    )

  saveRDS(paper_table, paste0("outputs/tables/", name, "_table.RDS"))
}

## for paper ----
files <- list.files("outputs/tables", full.names = TRUE)

paper_table <- readRDS(files[[3]]) %>%
  mutate(
    time_med = formatC(time_med, format = "f", digits = 3),
    time_iqr = formatC(time_iqr, format = "f", digits = 3),
    llike_med = formatC(llike_med, format = "f", digits = 2),
    llike_iqr = formatC(llike_iqr, format = "f", digits = 2),
    it_med = formatC(it_med, big.mark = ",", format = "f", digits = 0),
    it_iqr = formatC(it_iqr, big.mark = ",", format = "f", digits = 0),
    converged = paste0(100 * converged, "%"),
    # Combine
    time = paste0(time_med, "(", time_iqr, ")"),
    llike = paste0(llike_med, "(", llike_iqr, ")"),
    it = paste0(it_med, "(", it_iqr, ")"),
    method = factor(method, levels = methods)
  ) %>%
  arrange(method) %>%
  select(
    "Method" = method, "log-likelihood" = llike, "Iterations" = it,
    "CPU time (s)" = time, "% converged" = converged
  ) %>%
  xtable() %>%
  print(include.rownames = FALSE)

