#' @export
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