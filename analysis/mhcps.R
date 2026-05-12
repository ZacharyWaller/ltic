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
