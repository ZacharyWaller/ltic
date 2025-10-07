#' @export
print.ltic <- function(x, ...) {

  llike <- x$res$llike
  it <- x$res$it


  cat(
    "Non-parametric model",
    "\nNumber of inner-intervals:", length(x$intervals$II$left),
    "\nLog-likelihood:", llike,
    "\nIterations:", it
  )

  if (x$type == "prodlim") {
    n_int <- length(x$res$lambda)
    condition <- which(x$res$lambda[-n_int] >= 1 - 1E-9)

    if (length(condition) > 0) {
      cat(
        "\nInner-interval(s)",
        paste0(condition, collapse = ", "),
        "have a very high interval-hazard"
      )
    }
  }

}