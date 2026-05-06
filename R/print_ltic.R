#' Print an NPMLE estimate
#'
#' @param x Estimate from \code{ltic_np()}
#' @param ... Other printing options
#'
#' @export
#'
#' @details
#' Gives method used for calculating the NPMLE, details about the
#' inner-intervals, log-likelihood and number of iterations. For left-truncated
#' data it will also warn about inner-intervals that may cause a large drop in
#' survival estimate. See Hudgens (2005) and Waller et al. (2025) for
#' details.
#'
#' @references
#' Hudgens, Michael G. "On nonparametric maximum likelihood estimation with interval censoring and left truncation." Journal of the Royal Statistical Society Series B: Statistical Methodology 67.4 (2005): 573-587.
#' Waller, Zachary, et al. "A fast and stable algorithm for calculating the non-parametric maximum likelihood estimate ofleft-truncated and interval-censored data." (2025).
#'
#' @examples
#' est <- ltic_np(mhcps$Ui, mhcps$Vi, mhcps$Ti)
#' print(est)
#'
print.ltic <- function(x, ...) {

  llike <- x$res$llike
  it <- x$res$it

  cat(
    "Non-parametric model using", x$method, "algorithm",
    "\nNumber of inner-intervals:", length(x$intervals$II$left),
    "\nLog-likelihood:", llike,
    "\nIterations:", it
  )

  if (x$method %in% c("Product-Limit", "PL-ICM", "Breslow-ICM", "Breslow")) {
    n_int <- length(x$res$lambda)

    if (x$type == "expo") {
      diff_lambda <- diff(x$res$lambda)
      condition <- which(diff_lambda[-(n_int - 1)] >= 10.)
    } else {
      condition <- which(x$res$lambda[-n_int] >= 1 - 1E-9)
    }

    if (length(condition) > 0) {
      cat(
        "\nInner-interval(s)",
        paste0(condition, collapse = ", "),
        "have a very high interval-hazard"
      )
    }
  }

}
