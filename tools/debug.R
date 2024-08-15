devtools::load_all()

# Debug example 1 --------------------------------------------------------------
lambda <- c(0.1, 0.1, 0.1)
l <- c(0, 0, 1, 2)
r <- c(1, 3, 2, 3)
q <- c(1, 2, 3)
R0 <- c(4, 4, 4, 4)

newt <- newton_algorithm(lambda, l, r, q)
exp(-cumsum(newt$lambda_0))
newt

em <- em_algorithm(lambda, l, r, q, R0)
cumprod(1 - em$h)
em

comb <- combined_algorithm(lambda, l, r, q, R0)
comb
exp(-cumsum(comb$lambda_0))

which.max(c(comb$like, em$like, newt$like))
which.min(c(comb$it, em$it, newt$it))

comb$it
em$it
newt$it

# Debug example 2 --------------------------------------------------------------
l <- c(0, 0, 1, 2)
r <- c(1, 2, 3, 3)
q <- c(1, 2, 3)
lambda <- c(0.1, 0.1, 0.1)
R0 <- c(4, 4, 4, 4)

newt <- newton_algorithm(lambda, l, r, q)
exp(-cumsum(newt$lambda_0))
newt

em <- em_algorithm(lambda, l, r, q, R0)
cumprod(1 - em$h)
em

comb <- combined_algorithm(lambda, l, r, q, R0)
comb
exp(-cumsum(comb$lambda_0))

which.max(c(comb$like, em$like, newt$like))
which.min(c(comb$it, em$it, newt$it))

comb$it
em$it
newt$it

