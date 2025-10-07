#' NPMLE for left-truncated interval-censored data using quasi-Newton method
#'
#' @description Wrapper for the \code{optim()} function to calculate the NPMLE
#'
#' @param init Initial values for the probability masses in each inner-interval
#' @param alpha Indicator matrix of events in each inner-interval
#' @param beta Indicator matrix of truncation times for each inner-interval
#' @param tol Tolerance
#'
#' @details Using a method as roughly described by Hudgens (2005) (the exact
#' details aren't given there)
#'
#' @author Zachary Waller
#'
#' @references
#' Hudgens, M. G. (2005). On nonparametric maximum likelihood estimation with interval censoring and left truncation. Journal of the Royal Statistical Society Series B: Statistical Methodology, 67(4), 573-587.
#'
optim_method <- function(init, alpha, beta, tol = 1e-9) {

  n <- ncol(alpha)
  opt <- optim(
    par = apply(alpha/apply(alpha, 1, sum), 2, mean),
    fn = optim_likelihood,
    gr = like_gradient,
    lower = rep(tol, n), upper = rep(1 - tol, n),
    control = list(ndeps = rep(0, n), factr = 0, pgtol = 0),
    method = "L-BFGS-B",
    alpha = alpha, beta = beta
  )

  list(
    surv = 1 - c(0, cumsum(opt$par / sum(opt$par))),
    s    = opt$par / sum(opt$par),
    like = -opt$value,
    numit = opt$counts
  )

}


like_gradient <- function(steps, alpha, beta) {

  top <- alpha / colSums(t(alpha) * steps)
  bottom <- beta / colSums(t(beta) * steps)

  -colSums(top - bottom)

}

optim_likelihood <- function(steps, alpha, beta) {

  steps <- steps / sum(steps)
  top <- colSums(t(alpha) * steps)
  bottom <- colSums(t(beta) * steps)

  like <- sum(log(top) - log(bottom))

  -like
}
