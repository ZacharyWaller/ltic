# Indicator matrix -------------------------------------------------------
indicator_matrix <- function(II, interval, delta = NULL, type = c("integer", "logical")) {
  #' Check which inner intervals are subsets of other intervals. These other
  #' intervals could be observation intervals or truncation intervals, for example.
  #'
  #' @return Matrix of \code{M} x \code{N} where \code{M} is the number of
  #' censored intervals and \code{N} is the number of inner intervals.
  #' 1 where inner interval is subset of censored interval and 0 otherwise.
  #'
  #' @export
  #' @examples
  #' left <- c(0, 1, 2, 3, 4)
  #' right <- c(left + 1.5)
  #' turnbulls <- inner_intervals(left, right)
  #' indicator_matrix(turnbulls, left, right)

  type <- match.arg(type)

  points <- II
  left <- interval$left
  right <- interval$right


  if (is.null(delta)) {
    delta <- min(diff(sort(unique(c(II$left, II$right, left, right)))), na.rm = TRUE) / 2
  }

  if (is.infinite(delta)) {
    delta <- 1
  }

  II_right <- points$right
  II_left <- points$left

  # correct for open/closed intervals
  II_left[points$open_left] <- II_left[points$open_left] + delta
  II_right[points$open_right] <- II_right[points$open_right] - delta
  II_right[II_right == Inf] <- max(II_right[II_right != Inf])*10

  left[interval$open_left] <- left[interval$open_left] + delta
  right[interval$open_right] <- right[interval$open_right] - delta

  # use tolerances because some values can change by a small amount when adding/
  # subtracting delta
  bigger_than_left <- sapply(II_left, function(x) x - left >= -1e-10)
  smaller_than_right <- sapply(II_right, function(x) x - right <= 1e-10)

  if (type == "integer") {
    matrix <- bigger_than_left * smaller_than_right
  } else {
    matrix <- bigger_than_left & smaller_than_right
  }

  # turn into matrix if one dimensional (i.e. has only one participant)
  if (class(matrix)[1] %in% c("integer", "logical")) {
    matrix <- t(matrix(matrix))
  }

  # make labels ----------------------------------------------------------------
  # left_label <- ifelse(points$open_left, "(", "[")
  # right_label <- ifelse(points$open_right, ")", "]")
  #
  # oi_left_label <- ifelse(interval$open_left, "(", "[")
  # oi_right_label <- ifelse(interval$open_right, ")", "]")

  # dimnames(matrix) <- list(
  #   "between" = paste0(oi_left_label, interval$left, ", ", interval$right, oi_right_label),
  #   "interval" = paste0(left_label, points$left, ", ", points$right, right_label)
  # )

  matrix
}
