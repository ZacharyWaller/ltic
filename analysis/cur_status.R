# Examples from Jewell
n_sim <- 100
results <- replicate(
  n_sim,
  {
    obs_times <- seq(0.5, 5.0, by = 0.5)
    n <- 20000

    lambda <- 0.3

    fail <- rexp(n, lambda)
    obs_i <- sample(obs_times, n, replace = TRUE)

    time_1 <- fail
    time_2 <- fail
    time_1[fail < obs_i] <- -Inf
    time_2[fail > obs_i] <- Inf

    prodlim <- ltic_np(time_1, time_2, remove_rcens = FALSE,
                  method = "prodlim")
    bres <- ltic_np(time_1, time_2, remove_rcens = TRUE,
                  method = "breslow")


    turn <- ltic_np(time_1, time_2, remove_rcens = FALSE,
                  method = "turnbull")
    both <- ltic_np(time_1, time_2, remove_rcens = TRUE,
                  method = "both")
    bres_c <- ltic_np(time_1, time_2, remove_rcens = TRUE,
                  method = "bres_comb")

    f_1cont <- est$event_1
    f_2cont <- est$event_2

    f_1 <- cumsum(est$steps * colSums(est$event * as.logical(type_1[!r_cens])) / colSums(est$event))
    f_2 <- cumsum(est$steps * colSums(est$event * as.logical(type_2[!r_cens])) / colSums(est$event))
    f <- 1 - est$surv[-1]
    t <- est$intervals$II$right
    numit <- est$numit

    list(f_1 = f_1, f_2 = f_2, f = f, t = t, numit = numit, f_1cont = f_1cont, f_2cont = f_2cont)
  },
  simplify = FALSE
)

t <- unlist(sapply(results, function(x) x$t))
f_1 <- unlist(sapply(results, function(x) x$f_1))
f_2 <- unlist(sapply(results, function(x) x$f_2))
f_1_cont <- unlist(sapply(results, function(x) x$f_1cont))
f_2_cont <- unlist(sapply(results, function(x) x$f_2cont))
f <- unlist(sapply(results, function(x) x$f))

library(dplyr)
data <- data.frame(
  t = t, f_1 = f_1, f_2 = f_2, f = f, f_1_cont = f_1_cont, f_2_cont = f_2_cont
) %>%
  group_by(t) %>%
  summarise(
    m_f = mean(f),
    m_f_1 = mean(f_1),
    m_f_1_cont = mean(f_1_cont),
    m_f_2_cont = mean(f_2_cont),
    m_f_2 = mean(f_2),
    sd_f = sd(f),
    sd_f_1 = sd(f_1),
    sd_f_2 = sd(f_2),
    sd_f_1_cont = sd(f_1_cont),
    sd_f_2_cont = sd(f_2_cont)
  )
data
