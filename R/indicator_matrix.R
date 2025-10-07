# Indicator matrix -------------------------------------------------------
indicator_matrix <- function(
  inner_int, interval, delta = NULL, type = c("integer", "logical")
) {
  #' Check which inner-intervals are subsets of other intervals. These other
  #' intervals could be observation intervals or truncation intervals, for
  #' example.
  #'
  #' @return Matrix of \code{M} x \code{N} where \code{M} is the number of
  #' censored intervals and \code{N} is the number of inner intervals.
  #' 1 where inner interval is subset of censored interval and 0 otherwise.
  #'
  #' @examples
  #' left <- c(0, 1, 2, 3, 4)
  #' right <- c(left + 1.5)
  #' turnbulls <- inner_intervals(left, right)
  #' indicator_matrix(turnbulls, left, right)

  type <- match.arg(type)

  points <- inner_int
  left <- interval$left
  right <- interval$right

  # Get minium difference between time points
  if (is.null(delta)) {
    delta <- min(
      diff(sort(unique(c(inner_int$left, inner_int$right, left, right)))),
      na.rm = TRUE
    ) / 2
  }

  if (is.infinite(delta)) {
    delta <- 1
  }

  inner_int_right <- points$right
  inner_int_left <- points$left

  # correct for open/closed intervals
  inner_int_left[points$open_left] <- inner_int_left[points$open_left] + delta
  inner_int_right[points$open_right] <- inner_int_right[points$open_right] - delta
  inner_int_right[inner_int_right == Inf] <- max(inner_int_right[inner_int_right != Inf]) * 10

  left[interval$open_left] <- left[interval$open_left] + delta
  right[interval$open_right] <- right[interval$open_right] - delta

  # use tolerances because some values can change by a small amount when adding/
  # subtracting delta
  bigger_than_left <- sapply(inner_int_left, function(x) x - left >= 0)
  smaller_than_right <- sapply(inner_int_right, function(x) x - right <= 0)

  if (type == "integer") {
    matrix <- bigger_than_left * smaller_than_right
  } else {
    matrix <- bigger_than_left & smaller_than_right
  }

  # turn into matrix if one dimensional (i.e. has only one participant)
  if (class(matrix)[1] %in% c("integer", "logical")) {
    matrix <- t(matrix(matrix))
  }

  matrix
}
