library(ltic)

## Simulation ------------------------------------------------------------------
# Vary amount of truncation with interval censored data

# Make directory for results
if (!dir.exists("outputs/sim")) {
  dir.create("outputs/sim")
}


# Simulate data ----------------------------------------------------------------
simulate_data <- function(n, runif_max, seed) {
  set.seed(seed)
  len <- 0.5
  #trunc <- sample(seq(0, runif_max, by = runif_max / 52), n, replace = TRUE)
  trunc <- runif(n, 0, runif_max)
  failure_time <- rexp(n, 1)

  #first_obs <- sample(seq(0, 2, 2 / 52), n, replace = TRUE)#
  first_obs <- runif(n, 0, 2)
  second_obs <- first_obs + len

  # ignore failure times lower than first observation time
  keep <- trunc <= failure_time
  failure <- failure_time[keep]
  y_1 <- first_obs[keep]
  y_2 <- second_obs[keep]

  obs_list <- lapply(
    seq_along(failure),
    function(x) {
      t(between_obs(failure[[x]], c(y_1[[x]], y_2[[x]])))
    }
  )
  data.frame(Reduce(rbind, obs_list), trunc = trunc[keep])
}

# Simulation run ---------------------------------------------------------------
sim_run <- function(n, trunc_max, seed, methods) {
  data <- simulate_data(n = n, trunc_max, seed)

  left <- data$lower
  right <- data$upper
  trunc <- data$trunc

  runs <- lapply(
    methods,
    function(method) {
      ltic_np(left, right, trunc, max_it = max_it, tol = tol, method = method)
    }
  )
  names(runs) <- methods

  runs
}

# Initial values ---------------------------------------------------------------
# set seeds
run_date <- "2025-10-01 11:08:59 CET"
seed <- as.integer(as.POSIXct(run_date))
set.seed(seed)
seeds <- sample(seed,  1201)
# algorithm variables
tol <- 1e-5
max_it <- 1e5
methods <- c(
  "Breslow", "Product-Limit", "Shen", "Turnbull", "Yu",
  "Breslow-ICM", "PL-ICM", "Turnbull-ICM", "Yu-ICM", "ICM", "Quasi-Newton"
)

# data sim variables
trunc_max <- c(0, 0.6, 1.6) # maximum truncation time
n_factor <- c(1, 0.75, 0.5) # rough prop. of observations left after truncation
n <- c(500, 1000, 2000) # sample sizes
# other stuff
rep <- 100
i <- 1 # keeping track of seeds
j <- 1 # keep track of initial sample size
k <- 1 # keeping track of truncation
# Run --------------------------------------------------------------------------
for (j in 1:length(n)) {
  for (k in 1:length(trunc_max)) {
    # number of observations before truncation
    data_n <- n[j] / n_factor[k]

    ## Run actual simulations rep times ----
    results <- replicate(
      n = rep,
      expr = {
        cat("\n", i, j, k)
        i <<- i + 1
        sim_run(data_n, trunc_max[k], seeds[i], methods)
      },
      simplify = FALSE
    )

    # save
    file_name <- paste0("outputs/sim/res_", n[j], "_", 100 * n_factor[k], ".RDS")
    saveRDS(results, file_name)
  }
}

# Exact % left from truncation ----
# To end up with 75% left we need trunc_max = 0.6058601
# uniroot(
#   f = function(x, p) p + exp(-x) / x - 1 / x,
#   c(0.1, 4), p = c(0.75)
# )$root

# # To end up with 50% left we need trunc_max = 1.593615
# uniroot(
#   f = function(x, p) p + exp(-x) / x - 1 / x,
#   c(0.1, 4), p = c(0.5)
# )$root
