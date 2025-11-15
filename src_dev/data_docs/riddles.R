#' Riddles Simulation, Case 1: Linear Mean
#'
#' A simulated dataset of 500 observations based on
#' Simulation Study I (Model 1, Case 1) of Riddles, Kim, and Im (2016) [cite: 363-365].
#' The data features a nonignorable nonresponse mechanism where the
#' response probability depends on the study variable `y`.
#'
#' The mean structure is linear: `y = -1 + x + e`.
#' @keywords dataset
#' @format A data frame with 500 rows and 4 variables:
#' \describe{
#'   \item{x}{Numeric. The auxiliary variable, `x ~ N(0, 0.5)`.}
#'   \item{y}{Numeric. The study variable with nonignorable nonresponse.
#'              `y` contains `NA`s for nonrespondents.}
#'   \item{y_true}{Numeric. The complete, true value of `y` before
#'                 missingness was introduced.}
#'   \item{delta}{Integer. The response indicator (1=responded, 0=nonresponse).}
#' }
#' @source Riddles, M. K., Kim, J. K., & Im, J. (2016).
#'         A Propensity-Score-Adjustment Method for Nonignorable Nonresponse.
#'         *Journal of Survey Statistics and Methodology*, 4(1), 1–31.
"riddles_case1"


#' Riddles Simulation, Case 2: Exponential Mean
#'
#' A simulated dataset of 500 observations based on
#' Simulation Study I (Model 1, Case 2) of Riddles, Kim, and Im (2016) [cite: 363-365].
#' The data features a nonignorable nonresponse mechanism where the
#' response probability depends on the study variable `y`.
#'
#' The mean structure is exponential: `y = -2 + 0.5 * exp(x) + e`.
#' @keywords dataset
#' @format A data frame with 500 rows and 4 variables:
#' \describe{
#'   \item{x}{Numeric. The auxiliary variable, `x ~ N(0, 0.5)`.}
#'   \item{y}{Numeric. The study variable with nonignorable nonresponse.
#'              `y` contains `NA`s for nonrespondents.}
#'   \item{y_true}{Numeric. The complete, true value of `y` before
#'                 missingness was introduced.}
#'   \item{delta}{Integer. The response indicator (1=responded, 0=nonresponse).}
#' }
#' @source Riddles, M. K., Kim, J. K., & Im, J. (2016).
#'         A Propensity-Score-Adjustment Method for Nonignorable Nonresponse.
#'         *Journal of Survey Statistics and Methodology*, 4(1), 1–31.
"riddles_case2"


#' Riddles Simulation, Case 3: Sine Wave Mean
#'
#' A simulated dataset of 500 observations based on
#' Simulation Study I (Model 1, Case 3) of Riddles, Kim, and Im (2016) [cite: 363-365].
#' The data features a nonignorable nonresponse mechanism where the
#' response probability depends on the study variable `y`.
#'
#' The mean structure is a sine wave: `y = -1 + sin(2 * x) + e`.
#' @keywords dataset
#' @format A data frame with 500 rows and 4 variables:
#' \describe{
#'   \item{x}{Numeric. The auxiliary variable, `x ~ N(0, 0.5)`.}
#'   \item{y}{Numeric. The study variable with nonignorable nonresponse.
#'              `y` contains `NA`s for nonrespondents.}
#'   \item{y_true}{Numeric. The complete, true value of `y` before
#'                 missingness was introduced.}
#'   \item{delta}{Integer. The response indicator (1=responded, 0=nonresponse).}
#' }
#' @source Riddles, M. K., Kim, J. K., & Im, J. (2016).
#'         A Propensity-Score-Adjustment Method for Nonignorable Nonresponse.
#'         *Journal of Survey Statistics and Methodology*, 4(1), 1–31.
"riddles_case3"


#' Riddles Simulation, Case 4: Cubic Mean
#'
#' A simulated dataset of 500 observations based on
#' Simulation Study I (Model 1, Case 4) of Riddles, Kim, and Im (2016) [cite: 363-365].
#' The data features a nonignorable nonresponse mechanism where the
#' response probability depends on the study variable `y`.
#'
#' The mean structure is cubic: `y = -1 + 0.4 * x^3 + e`.
#' @keywords dataset
#' @format A data frame with 500 rows and 4 variables:
#' \describe{
#'   \item{x}{Numeric. The auxiliary variable, `x ~ N(0, 0.5)`.}
#'   \item{y}{Numeric. The study variable with nonignorable nonresponse.
#'              `y` contains `NA`s for nonrespondents.}
#'   \item{y_true}{Numeric. The complete, true value of `y` before
#'                 missingness was introduced.}
#'   \item{delta}{Integer. The response indicator (1=responded, 0=nonresponse).}
#' }
#' @source Riddles, M. K., Kim, J. K., & Im, J. (2016).
#'         A Propensity-Score-Adjustment Method for Nonignorable Nonresponse.
#'         *Journal of Survey Statistics and Methodology*, 4(1), 1–31.
"riddles_case4"
