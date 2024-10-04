#' @param left
#' @param right
#' @param trunc
#' @param tol
#'
#' @return
#' @export
#'
#' @examples
optim_method <- function(init, left, right, trunc, tol = 1e-9) {


  n_int <- length(init)
  n_obs <- length(left)
  opt <- optim(
    par = steps,
    fn = optim_likelihood,
    gr = like_gradient,
    lower = rep(tol, n), upper = rep(1 - tol, n),
    control = list(ndeps = rep(0, n), factr = 0, pgtol = 0),
    method = "L-BFGS-B",
    alpha = alpha, beta = beta
  )

  list(
    surv = 1 - c(0, cumsum(opt$par / sum(opt$par))),
    like = -opt$value,
    numit = opt$counts
  )

}

#' @param steps
#' @param alpha
#' @param beta
#'
#' @return
#' @export
#'
#' @examples
like_gradient <- function(steps, alpha, beta) {

  top <- alpha / colSums(t(alpha) * steps)
  bottom <- beta / colSums(t(beta) * steps)

  -colSums(top - bottom)

}

#' @param steps
#' @param alpha
#' @param beta
#'
#' @return
#' @export
#'
#' @examples
optim_likelihood <- function(steps, alpha, beta) {

  steps <- steps / sum(steps)
  top <- colSums(t(alpha) * steps)
  bottom <- colSums(t(beta) * steps)

  like <- sum(log(top) - log(bottom))

  -like
}