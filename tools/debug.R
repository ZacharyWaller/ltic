devtools::clean_dll()
devtools::load_all()

#sourceCpp("src/newton.cpp")

lambda <- c(0.1, 0.1, 0.1)
cum_lambda <- cumsum(lambda)
l <- c(1, 1, 2)
r <- c(2, 3, 3)
q <- c(1, 2, 3)

calc_deriv_R <- function(cum_lambda, l, r, q) {
    .Call('src/newton.o', cum_lambda, l, r, q)
}
calc_deriv_R(cum_lambda, l, r, q)
