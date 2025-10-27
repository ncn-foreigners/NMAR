#' Internal stats imports for R CMD check
#'
#' Centralized imports for base/statistics functions used across the package
#' to satisfy R CMD check visibility. No runtime code.
#'
#' @name nmar_internal_stats_imports
#' @keywords internal
#' @noRd
#' @importFrom stats AIC Gamma confint cov dexp dlnorm gaussian getCall glm lm median model.matrix plogis pnorm rnorm update var resid
#' @importFrom utils modifyList
NULL
