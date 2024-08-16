# Some speed tests using simulated data ----------------------------------------
library(dplyr)
devtools::document()
devtools::load_all()

n <- 10000
failure_times <- rweibull(n, 1, 1)

intervals <- sapply(
  failure_times,
  function(failure_time) {
    observations <- c(0, runif(1, 0, 1))
    t(between_obs(failure_time, observations))
  },
  simplify = FALSE
)

transitions <- as.data.frame(do.call(rbind, intervals))

left <- transitions$lower
right <- transitions$upper

t0 <- Sys.time()
ltic_est <- ltic_np(left, right, method = "em", remove_rcens = TRUE)
Sys.time() - t0

t0 <- Sys.time()
icen_est <- icenReg::ic_np(
  formula = cbind(left, right)
)
Sys.time() - t0

icen_est

icen_est$llk
ltic_est$like

1 - cumsum(icen_est$p_hat)
exp(-cumsum(ltic_est$lambda_1))

ltic_est$it
icen_est$iterations

