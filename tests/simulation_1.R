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

seeds <- 1:500
tol <- 1e-5
max_it <- 1e5
i <- 1
# Simulate data and run NPMLE algorithms ---------------------------------------
## 100% left -------------------------------------------------------------------
results_0 <- replicate(
  n = 100,
  expr = {
    data <- simulate_data(n = 1000, 0.0, seeds[i])
    i <<- i + 1

    left <- data$lower
    right <- data$upper
    trunc <- data$trunc
    
    prod <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "prodlim")
    icm <- ltic_np(left, right, trunc, tol = tol, max_it = max_it, method = "icm", remove_rcens = FALSE)
    comb <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "both")
    turn_comb <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "turn_comb", remove_rcens = FALSE)
    yu_comb <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "yu_comb", remove_rcens = FALSE)
    bres_comb <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "bres_comb")
    bres <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "breslow")
    opti <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "optim", remove_rcens = FALSE)
    yu <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "yu", remove_rcens = FALSE)
    shen <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "shen", remove_rcens = FALSE)
    turn <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "turnbull", remove_rcens = FALSE)

    list(
      prod = c(prod),
      icm = c(icm),
      comb = c(comb),
      tu_c = c(turn_comb),
      yu_c = c(yu_comb),
      br_c = c(bres_comb),
      bres = c(bres),
      opti = c(opti),
      yu = c(yu),
      shen = c(shen),
      turn = c(turn)
    )
  },
  simplify = FALSE
)
saveRDS(results_0, "outputs/results_0.RDS")
rm(results_0)

i <- 101
## 75% left --------------------------------------------------------------------
results_075 <- replicate(
  n = 100,
  expr = {
    data <- simulate_data(n = 1333, 0.6, seeds[i])
    i <<- i + 1

    left <- data$lower
    right <- data$upper
    trunc <- data$trunc

    prod <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "prodlim")
    icm <- ltic_np(left, right, trunc, tol = tol, max_it = max_it, method = "icm", remove_rcens = FALSE)
    comb <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "both")
    turn_comb <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "turn_comb", remove_rcens = FALSE)
    yu_comb <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "yu_comb", remove_rcens = FALSE)
    bres_comb <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "bres_comb")
    bres <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "breslow")
    opti <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "optim", remove_rcens = FALSE)
    yu <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "yu", remove_rcens = FALSE)
    shen <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "shen", remove_rcens = FALSE)
    turn <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "turnbull", remove_rcens = FALSE)

    list(
      prod = c(prod),
      icm = c(icm),
      comb = c(comb),
      tu_c = c(turn_comb),
      yu_c = c(yu_comb),
      br_c = c(bres_comb),
      bres = c(bres),
      opti = c(opti),
      yu = c(yu),
      shen = c(shen),
      turn = c(turn)
    )
  },
  simplify = FALSE
)
saveRDS(results_075, "outputs/results_075.RDS")
rm(results_075)

i <- 201
## 50% left --------------------------------------------------------------------
results_050 <- replicate(
  n = 100,
  expr = {
    data <- simulate_data(n = 2000, 1.6, seeds[i])
    i <<- i + 1

    left <- data$lower
    right <- data$upper
    trunc <- data$trunc

    prod <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "prodlim")
    icm <- ltic_np(left, right, trunc, tol = tol, max_it = max_it, method = "icm", remove_rcens = FALSE)
    comb <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "both")
    turn_comb <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "turn_comb", remove_rcens = FALSE)
    yu_comb <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "yu_comb", remove_rcens = FALSE)
    bres_comb <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "bres_comb")
    bres <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "breslow")
    opti <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "optim", remove_rcens = FALSE)
    yu <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "yu", remove_rcens = FALSE)
    shen <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "shen", remove_rcens = FALSE)
    turn <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "turnbull", remove_rcens = FALSE)

    list(
      prod = c(prod),
      icm = c(icm),
      comb = c(comb),
      tu_c = c(turn_comb),
      yu_c = c(yu_comb),
      br_c = c(bres_comb),
      bres = c(bres),
      opti = c(opti),
      yu = c(yu),
      shen = c(shen),
      turn = c(turn)
    )

  },
  simplify = FALSE
)
saveRDS(results_050, "outputs/results_050.RDS")
rm(results_050)

## 25% left --------------------------------------------------------------------
results_025 <- replicate(
  n = 100,
  expr = {
    data <- simulate_data(n = 4000, 3.9, seeds[i])
    i <<- i + 1

    left <- data$lower
    right <- data$upper
    trunc <- data$trunc

    prod <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "prodlim", remove_rcens = TRUE)
    icm <- ltic_np(left, right, trunc, tol = tol, max_it = max_it, method = "icm", remove_rcens = FALSE)
    comb <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it / 10, method = "both", remove_rcens = TRUE)
    turn_comb <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it / 10, method = "turn_comb", remove_rcens = FALSE)
    yu_comb <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it / 10, method = "yu_comb", remove_rcens = FALSE)
    bres_comb <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it / 10, method = "bres_comb")
    bres <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "breslow")
    opti <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "optim", remove_rcens = FALSE)
    yu <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "yu", remove_rcens = FALSE)
    shen <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "shen", remove_rcens = FALSE)
    turn <- ltic_np(left, right, trunc, tol = tol, max_it =  max_it, method = "turnbull", remove_rcens = FALSE)

    list(
      prod = c(prod),
      comb = c(comb),
      icm  = c(icm),
      tu_c = c(turn_comb),
      yu_c = c(yu_comb),
      br_c = c(bres_comb),
      bres = c(bres),
      opti = c(opti),
      yu = c(yu),
      shen = c(shen),
      turn = c(turn)
    )

  },
  simplify = FALSE
)
saveRDS(results_025, "outputs/results_025.RDS")
rm(results_025)
