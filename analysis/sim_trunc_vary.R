devtools::load_all()
library(microbenchmark)

# Simulate data ----------------------------------------------------------------
simulate_data <- function(n, runif_max) {
  len <- 0.5
  trunc <- runif(n, 0, runif_max)
  failure_time <- rexp(n, 1)

  first_obs <- runif(n, 0, 2)
  second_obs <- first_obs + len

  # ignore failure times lower than first observation time
  keep <- trunc <= failure_time
  failure <- failure_time[keep]
  Y_1 <- first_obs[keep]
  Y_2 <- second_obs[keep]

  obs_list <- lapply(
    seq_along(failure),
    function(x) {
      t(between_obs(failure[[x]], c(Y_1[[x]], Y_2[[x]])))
    }
  )

  data.frame(Reduce(rbind, obs_list), trunc = trunc[keep])

}

tol <- 1e-5
max_it <- 1e5
# Simulate data and run NPMLE algorithms ---------------------------------------
results_0 <- replicate(
  n = 100,
  expr = {
    data <- simulate_data(n = 2000, 0)
    microbenchmark(
      both <- ltic_np(data$lower, data$upper, data$trunc, method = "both"),
      both_rcens <- ltic_np(data$lower, data$upper, data$trunc, method = "both", remove_rcens = FALSE),
      bres <- ltic_np(data$lower, data$upper, data$trunc, method = "breslow"),
      prod <- ltic_np(data$lower, data$upper, data$trunc, method = "prodlim"),
      bres_comb <- ltic_np(data$lower, data$upper, data$trunc, method = "bres_comb"),
      times = 2
    )
    t_prodlim <- system.time({est_prodlim <- sc_prodlim_method(data$lower, data$upper, data$trunc, tol = tol, max_iterations = max_it)})
    t_robust <- system.time({est_robust <- sc_prodlim_proper(data$lower, data$upper, data$trunc, tol = tol, max_iterations = max_it)})
    t_turnbull <- system.time({est_turnbull <- sc_turnbull_method(data$lower, data$upper, data$trunc, tol = tol, max_iterations = max_it)})
    t_yu <- system.time({est_yu <- sc_yu_method(data$lower, data$upper, data$trunc, tol = tol, max_iterations = max_it)})
    t_shen <- system.time({est_shen <- sc_shen_method(data$lower, data$upper, data$trunc, tol = tol, max_iterations = max_it)})
    t_optim <- system.time({est_optim <- optim_method(data$lower, data$upper, data$trunc)})

    list(
      prodlim = c(t_problim = t_prodlim, est_prodlim),
      robust = c(t_robust = t_robust, est_robust),
      turnbull = c(t_turnbull = t_turnbull, est_turnbull),
      yu = c(t_yu = t_yu, est_yu),
      shen = c(t_shen = t_shen, est_shen),
      optim = c(t_optim = t_optim, est_optim)
    )

  },
  simplify = FALSE
)

results_075 <- replicate(
  n = 100,
  expr = {
    data <- simulate_data(n = 1333, 0.6)
    t_prodlim <- system.time({est_prodlim <- sc_prodlim_method(data$lower, data$upper, data$trunc, tol = tol, max_iterations = max_it)})
    t_robust <- system.time({est_robust <- sc_prodlim_proper(data$lower, data$upper, data$trunc, tol = tol, max_iterations = max_it)})
    t_turnbull <- system.time({est_turnbull <- sc_turnbull_method(data$lower, data$upper, data$trunc, tol = tol, max_iterations = max_it)})
    #t_yu <- system.time({est_yu <- sc_yu_method(data$lower, data$upper, data$trunc, tol = tol, max_iterations = max_it)})
    t_shen <- system.time({est_shen <- sc_shen_method(data$lower, data$upper, data$trunc, tol = tol, max_iterations = max_it)})
    t_optim <- system.time({est_optim <- optim_method(data$lower, data$upper, data$trunc)})

    list(
      prodlim = c(t_problim = t_prodlim, est_prodlim),
      robust = c(t_robust = t_robust, est_robust),
      turnbull = c(t_turnbull = t_turnbull, est_turnbull),
      yu = c(t_yu = t_yu, est_yu),
      shen = c(t_shen = t_shen, est_shen),
      optim = c(t_optim = t_optim, est_optim)
    )

  },
  simplify = FALSE
)

results_050 <- replicate(
  n = 100,
  expr = {
    data <- simulate_data(n = 2000, 1.6)
    t_prodlim <- system.time({est_prodlim <- sc_prodlim_method(data$lower, data$upper, data$trunc, tol = tol, max_iterations = max_it)})
    t_robust <- system.time({est_robust <- sc_prodlim_proper(data$lower, data$upper, data$trunc, tol = tol, max_iterations = max_it)})
    t_turnbull <- system.time({est_turnbull <- sc_turnbull_method(data$lower, data$upper, data$trunc, tol = tol, max_iterations = max_it)})
    #t_yu <- system.time({est_yu <- sc_yu_method(data$lower, data$upper, data$trunc, tol = tol, max_iterations = max_it)})
    t_shen <- system.time({est_shen <- sc_shen_method(data$lower, data$upper, data$trunc, tol = tol, max_iterations = max_it)})
    t_optim <- system.time({est_optim <- optim_method(data$lower, data$upper, data$trunc)})

    list(
      prodlim = c(t_problim = t_prodlim, est_prodlim),
      robust = c(t_robust = t_robust, est_robust),
      turnbull = c(t_turnbull = t_turnbull, est_turnbull),
      yu = c(t_yu = t_yu, est_yu),
      shen = c(t_shen = t_shen, est_shen),
      optim = c(t_optim = t_optim, est_optim)
    )

  },
  simplify = FALSE
)


results_025 <- replicate(
  n = 100,
  expr = {
    data <- simulate_data(n = 4000, 3.9)
    t_prodlim <- system.time({est_prodlim <- sc_prodlim_method(data$lower, data$upper, data$trunc, tol = tol, max_iterations = max_it)})
    t_robust <- system.time({est_robust <- sc_prodlim_proper(data$lower, data$upper, data$trunc, tol = tol, max_iterations = max_it)})
    t_turnbull <- system.time({est_turnbull <- sc_turnbull_method(data$lower, data$upper, data$trunc, tol = tol, max_iterations = max_it)})
    #t_yu <- system.time({est_yu <- sc_yu_method(data$lower, data$upper, data$trunc, tol = tol, max_iterations = max_it)})
    t_shen <- system.time({est_shen <- sc_shen_method(data$lower, data$upper, data$trunc, tol = tol, max_iterations = max_it)})
    t_optim <- system.time({est_optim <- optim_method(data$lower, data$upper, data$trunc)})

    list(
      prodlim = c(t_problim = t_prodlim, est_prodlim),
      robust = c(t_robust = t_robust, est_robust),
      turnbull = c(t_turnbull = t_turnbull, est_turnbull),
      yu = c(t_yu = t_yu, est_yu),
      shen = c(t_shen = t_shen, est_shen),
      optim = c(t_optim = t_optim, est_optim)
    )

  },
  simplify = FALSE
)



# Extract and collate times ----------------------------------------------------
t_results <- lapply(
  results_075,
  function(x) {
    data.frame(
      prodlim = x$prodlim$t_problim.elapsed,
      robust = x$robust$t_robust.elapsed,
      turnbull = x$turnbull$t_turnbull.elapsed,
      yu = x$yu$t_yu.elapsed,
      shen = x$shen$t_shen.elapsed,
      optim = x$optim$t_optim.elapsed
    )
  }
)

# Extract and collate numits ---------------------------------------------------
numit_results <- lapply(
  results_075,
  function(x) {
    data.frame(
      prodlim = x$prodlim$numit,
      robust = x$robust$numit,
      turnbull = x$turnbull$numit,
      yu = x$yu$numit,
      shen = x$shen$numit,
      optim = x$optim$numit[[1]]
    )
  }
)

# results_025# Make plots -------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)

## Times ----
bind_rows(t_results) %>%
  mutate(run = row_number()) %>%
  pivot_longer(c(prodlim, robust, turnbull, shen, optim)) %>%
  ggplot(aes(name, value)) +
  geom_boxplot()

## numits ----
bind_rows(numit_results) %>%
  mutate(run = row_number()) %>%
  pivot_longer(c(prodlim, robust, turnbull, shen, optim)) %>%
  ggplot(aes(name, value)) +
  geom_boxplot()


numit_mean <- bind_rows(numit_results) %>% summarise(across(everything(), mean))
numit_sd <- bind_rows(numit_results) %>% summarise(across(everything(), sd))
paste0(numit_mean, "(", numit_sd, ") &", collapse = " ")

bind_rows(numit_results) %>% summarise(across(everything(), median))
bind_rows(numit_results) %>% summarise(across(everything(), IQR))

time_mean <- bind_rows(t_results) %>% summarise(across(everything(), mean))
time_sd <- bind_rows(t_results) %>% summarise(across(everything(), sd))
paste0(time_mean, "(", time_sd, ") &", collapse = " ")

bind_rows(t_results) %>% summarise(across(everything(), median))
bind_rows(t_results) %>% summarise(across(everything(), IQR))

# Save results -----------------------------------------------------------------
saveRDS(results_0, "outputs/sim_trunc_0.rds")
saveRDS(results_050, "outputs/sim_trunc_050.rds")
saveRDS(results_075, "outputs/sim_trunc_075.rds")



like_results <- lapply(
  results_075,
  function(x) {
    data.frame(
      prodlim = x$prodlim$like,
      robust = x$robust$like,
      turnbull = x$turnbull$like,
      yu = x$yu$like,
      shen = x$shen$like,
      optim = x$optim$like
    )
  }
)

bind_rows(like_results) %>%
  rowwise() %>%
  mutate(
    best = which.max(c(prodlim, robust, turnbull, shen, optim)),
  ) %>% View()

bind_rows(like_results) %>%
  mutate(run = row_number()) %>%
  rowwise() %>%
  mutate(range = diff(range(c(robust, optim)))) %>%
  ungroup() %>%
  pivot_longer(c(robust, optim)) %>%
  ggplot(aes(run, value, colour = name)) +
  geom_line()

like_mean <- bind_rows(like_results) %>% summarise(across(everything(), mean))
like_sd <- bind_rows(like_results) %>% summarise(across(everything(), sd))
paste0(like_mean, "(", like_sd, ") &", collapse = " ")


seed <- sample(1:1000, 1)
set.seed(seed)
m <- 100

tol <- 1e-5
max_it <- 1e6

results <- replicate(
  1000,
  {
    t_1 <- Sys.time()
    turnbull_est <- sc_turnbull_method(data$left, data$right, data$trunc, tol = tol, max_iterations = max_it)
    Sys.time() - t_1

    t_1 <- Sys.time()
    my_est <- sc_prodlim_method(data$left, data$right, data$trunc, tol = tol, max_iterations = max_it)
    Sys.time() - t_1

    plot_estimate(my_est)
    plot_estimate(turnbull_est)

    t_1 <- Sys.time()
    yu_est <- sc_yu_method(data$left, data$right, data$trunc, tol = tol, max_iterations = max_it)
    Sys.time() - t_1

    plot_estimate(yu_est)

    t_1 <- Sys.time()
    test_est <- sc_prodlim_proper(data$left, data$right, data$trunc, tol = tol, max_iterations = max_it)
    Sys.time() - t_1
  }
)
data <- simulate_data(3000)



test_est$like
turnbull_est$like
yu_est$like
my_est$like

yu_est$numit
my_est$numit
turnbull_est$numit

plot(my_est$surv)
lines(yu_est$surv, col = "red")
lines(turnbull_est$surv, col = "blue")

# likelihood
likelihood <- function(steps, alpha, beta) {

  steps <- steps / sum(steps)
  top <- colSums(t(alpha) * steps)
  bottom <- colSums(t(beta) * steps)

  like <- sum(log(top) - log(bottom))
  #like <- prod(top / bottom)

  - like
}

intervals <- inner_intervals(data$left, data$right, data$trunc)
alpha <- indicator_matrix(intervals$II, intervals$Oi)
beta <- indicator_matrix(intervals$II, intervals$Ti)

-likelihood(yu_est$steps, alpha, beta)
-likelihood(turnbull_est$steps, alpha, beta)

steps_my <- my_est$steps
steps_my[steps_my == 0] <- 1e-8
-likelihood(steps_my, alpha, beta)


t1 <- Sys.time()
opt <- optim_method(data$left, data$right, data$trunc)
Sys.time() - t1


plot_estimate(opt)
plot_estimate(test_est, plot_over = TRUE, col = "red")

opt$like
test_est$like
test_est$lambda

plot(product_survival(test_est$lambda[-1]), type = "l", lwd = 2)
lines(opt$surv)
