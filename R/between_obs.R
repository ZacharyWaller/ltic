# between_obs() ----------------------------------------------------------------
between_obs <- function(event_time, assess) {
  #' Check which assess a point lies between
  #'
  #' @param event_time A vector of event times
  #' @param assess A vector of assessment times to check. Not necessarily
  #' in order
  #'
  #' @return A matrix of lower and upper bounds from \code{assess} that
  #' each \code{event_time} lies between.
  #'
  #' @details This treats the intervals between observation times as open on
  #' the bottom and closed on top. If \code{event_time} lies exactly on an
  #' observation, that observation will be the upper.
  #' If \code{event_time} doesn't lie between any of the assess, +/-
  #' \code{Inf}
  #' will be returned along with the closest observation.
  #'
  #' @export
  #' @examples
  #' between_obs(0.5, c(0, 1, 2, 3, 4)) # returns 0, 1
  #' between_obs(5, c(1, 2, 3, 4, 5, 6)) # returns 5, 5
  #' between_obs(0, c(1, 2, 3, 4, 5, 6)) # returns -Inf, 1
  #' between_obs(9, c(1, 2, 3, 4, 5, 6)) # returns 6, Inf

  between_inner <- function(event_time, assess) {

    bigger <- assess >= event_time
    smaller <- assess < event_time

    if (all(!bigger)) {
      upper <- Inf
    } else {
      upper <- min(assess[bigger])
    }

    if (all(!smaller)) {
      lower <- -Inf
    } else {
      lower <- max(assess[smaller])
    }

    c("lower" = lower, "upper" = upper)

  }


  sapply(event_time, function(x) between_inner(x, assess), simplify = "vector")

}

