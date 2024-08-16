# Find the Inner intervals -----------------------------------------------------
inner_intervals <- function(left, right, trunc = NULL, open_L = TRUE, open_T = FALSE) {
  #' Find inner-intervals from observation times
  #'
  #' @param left Left end point of interval
  #' @param right Right end point of interval
  #' @param trunc Truncation time
  #'
  #' @return Data-frame containing left and right inner-intervals and open/closed
  #' status of each end-point.
  #' @export
  #' @references Qiqing Yu "The generalised MLE with truncated interval-censored data 2022"
  #'
  #' @examples
  #' trunc <- c(0,   0, 0,   1,    1,   1)
  #' left <- c(-Inf, 1, 2,   -Inf, 1,   2)
  #' right <- c(1,   2, Inf, 1,    2, Inf)
  #' inner_intervals(left, right, trunc)
  #' turnbull_est(left, right, trunc)

  if (length(left) != length(right)) {
    stop("left and right should be the same length.")
  }

  if (is.null(trunc)) {
    trunc <- rep(0, length(left))
  }

  delta <- min(diff(sort(unique(c(left, right, trunc)))), na.rm = TRUE) / 2

  # renaming in line with Yu (2022) notation
  L_i <- left
  R_i <- right
  T_i <- trunc

  # left made from max of [T_i, and (L_i
  W_i <- pmax(L_i, T_i)

  # TRUE for open, FALSE for closed
  open_W_i <- rep(open_L, length(W_i))
  open_W_i[W_i != L_i] <- open_T

  T_min <- rep(min(T_i), length(T_i))
  open_T_min <- rep(open_T, length(T_min))

  all_left <- c(W_i, T_min)
  open_left <- c(open_W_i, open_T_min)

  all_right <- c(R_i, T_i)
  open_right <- c(
    rep(FALSE, length(R_i)),
    rep(!open_T, length(T_i))
  )

  # keep for end
  Oi <- data.frame(
    left = W_i,
    right = R_i,
    open_left = open_W_i,
    open_right = rep(FALSE, length(R_i))
  )

  # add small number to account for open intervals
  all_right[open_right] <- all_right[open_right] - delta
  all_left[open_left] <- all_left[open_left] + delta

  # looking for inner intervals - left followed by right
  all <- c(all_left, all_right)
  all_open <- c(open_left, open_right)
  left_right <- c(rep(0, length(all_left)), rep(1, length(all_right)))

  ord <- order(all)

  left_then_right <- diff(left_right[ord])
  II_flag <- which(left_then_right == 1)
  II_left <- all[ord][II_flag]
  II_right <- all[ord][II_flag + 1]

  # keep track of open/closed ends of inner intervals
  II_left_open <- all_open[ord][II_flag]
  II_right_open <- all_open[ord][II_flag + 1]

  # remove small number for display
  II_left[II_left_open] <- II_left[II_left_open] - delta
  II_right[II_right_open] <- II_right[II_right_open] + delta

  # make into data frame
  list(
    II = data.frame(
      left = II_left,
      right = II_right,
      open_left = II_left_open,
      open_right = II_right_open
    ),
    Oi = Oi,
    Ti = data.frame(
      left = T_i,
      right = Inf,
      open_left = open_T,
      open_right = !open_T
    ),
    delta = delta
  )
}
