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
    par = init,
    fn = calc_like_r,
    gr = calc_derivs_r,
    lower = rep(tol, n_int), upper = rep(1 - tol, n_int),
    control = list(ndeps = rep(0, n_int), factr = 0, pgtol = 0),
    method = "L-BFGS-B",
    left = left, right = right, trun = trunc, n_obs = n_obs, n_int = n_int
  )

  list(
    surv = 1 - cumsum(opt$par / sum(opt$par)),
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