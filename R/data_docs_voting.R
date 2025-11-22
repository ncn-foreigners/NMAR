#' Aggregated Exit Poll Data for Gangdong-Gap (2012)
#'
#' This dataset contains the aggregated exit poll results for the
#' Gangdong-Gap district in Seoul from the 2012 nineteenth South Korean
#' legislative election. The data is transcribed directly from Table 9 of
#' Riddles, Kim, and Im (2016).
#'
#' In the paper's application, `Gender` is used as the nonresponse instrumental
#' variable and `Age_group` is the primary auxiliary variable .
#' @keywords dataset
#' @format A data frame with 8 rows and 7 variables:
#' \describe{
#'   \item{Gender}{Factor. The gender of the voter ("Male", "Female").}
#'   \item{Age_group}{Character. The age group of the voter.}
#'   \item{Voted_A}{Numeric. Count of respondents voting for Party A.}
#'   \item{Voted_B}{Numeric. Count of respondents voting for Party B.}
#'   \item{Other}{Numeric. Count of respondents voting for another party.}
#'   \item{Refusal}{Numeric. Count of sampled individuals who refused
#'                  to respond (this is the nonresponse count).}
#'   \item{Total}{Numeric. Total individuals sampled in the group
#'                (Responders + Refusals).}
#' }
#' @source Riddles, M. K., Kim, J. K., & Im, J. (2016).
#'         A Propensity-Score-Adjustment Method for Nonignorable Nonresponse.
#'         *Journal of Survey Statistics and Methodology*, 4(1), 1–31.
#'         (Data from Table 9, p. 20).
"voting"
