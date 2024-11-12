# Simulation 1 -----------------------------------------------------------------
# Vary amount of truncation with interval censored data, leaving a remaining 
# sample size of around 1000
devtools::load_all()

# Simulate data ----------------------------------------------------------------
simulate_data <- function(n, runif_max, seed) {
  set.seed(seed)
  len <- 0.5
  trunc <- runif(n, 0, runif_max)
  failure_time <- rexp(n, 1)

  delta <- min(diff(sort(failure_time)))
  first_obs <- failure_time - delta / 3.
  second_obs <- failure_time + delta / 3.

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

seeds <- 1:500
tol <- 1e-5
max_it <- 1e5
i <- 1
# Simulate data and run NPMLE algorithms ---------------------------------------
## 75% left --------------------------------------------------------------------
n <- 1000
data <- simulate_data(n, 0.6, sample(1000, 1))

left <- data$lower
right <- data$upper
trunc <- data$trunc

cond <- ltic_np(left, right, trunc, tol = tol, max_it = max_it)
vardi <- lb_np(left, right, tol = tol, max_it =  max_it)
i <- i + 1
