library(Rcpp)

#sourceCpp("src/newton.cpp")

lambda <- c(0.1, 0.1, 0.1)
cum_lambda <- cumsum(lambda)
l <- c(1, 1, 2)
r <- c(2, 3, 3)
q <- c(1, 2, 3)

calc_deriv(cum_lambda, l, r, q)
