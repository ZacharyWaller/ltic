#' Massachusetts Health Care Panel Study (MHCPS) data
#'
#' @description
#' Dataset on the loss of functional independence in older adults.
#'
#' @format ## `mhcps`
#' A data frame with 1,030 rows and 4 columns:
#' \describe{
#'   \item{Ti}{Age at entry into study (truncation time)}
#'   \item{Ui}{Last age at assessment when still functionally independent (left
#'   end-point)}
#'   \item{Vi}{First age at assessment when not functionally indepdendent (right
#'   end-point)}
#'   \item{Zi}{Sex, 0 for female, 1 for male.}
#' }
#' @references
#' Pan, Wei, and Rick Chappell. "Estimation in the Cox proportional hazards model with left‐truncated and interval‐censored data." Biometrics 58.1 (2002): 64-70.
#'
#' @source <https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fj.0006-341X.2002.00064.x&file=BIOM_64_sm_010423.txt>
"mhcps"
