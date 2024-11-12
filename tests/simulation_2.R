# Example from Pan and Chappell 1998
devtools::load_all()

# Simulate data ----------------------------------------------------------------
simulate_data <- function(n) {
  len <- 0.5
  first_obs <- runif(n, 0, 4)
  failure_time <- rgamma(n, 2, 1)

  # ignore failure times lower than first observation time
  keep <- first_obs < failure_time
  failure <- failure_time[keep]
  t_0 <- first_obs[keep]

  before_obs_1 <- failure <= t_0 + len
  after_obs_1 <- failure > t_0 + len

  left <- c(t_0[before_obs_1], t_0[after_obs_1] + len)
  right <- c(t_0[before_obs_1] + len, t_0[after_obs_1] + Inf)
  trunc <- c(t_0[before_obs_1], t_0[after_obs_1])

  data.frame(left = left, right = right, trunc = trunc)
}

tol <- 1e-5
max_it <- 1e5
## Run simulation --------------------------------------------------------------
results <- replicate(
  n = 100,
  expr = {
    data <- simulate_data(n = 5000)
    left <- data$left
    right <- data$right
    trunc <- data$trunc

    prod <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "prodlim")
    comb <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "both")
    turn_comb <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "turn_comb", remove_rcens = FALSE)
    yu_comb <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "yu_comb", remove_rcens = FALSE)
    bres <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "breslow")
    opti <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "optim", remove_rcens = FALSE)
    yu <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "yu", remove_rcens = FALSE)
    shen <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "shen", remove_rcens = FALSE)
    turn <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "turnbull", remove_rcens = FALSE)

    list(
      prod = c(prod),
      comb = c(comb),
      tu_c = c(turn_comb),
      yu_c = c(yu_comb),
      bres = c(bres),
      opti = c(opti),
      yu = c(yu),
      shen = c(shen),
      turn = c(turn)
    )

  },
  simplify = FALSE
)
saveRDS(results, "outputs/results.RDS")
rm(results)
