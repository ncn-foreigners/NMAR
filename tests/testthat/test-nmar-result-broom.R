skip_if_not_installed("generics")

make_rich_result <- function() {
  res <- NMAR:::new_nmar_result(
    estimate = 0.5,
    estimate_name = "Y_mean",
    se = 0.1,
    converged = TRUE,
    model = list(
      coefficients = c(beta0 = 1.2, beta1 = -0.3),
      vcov = matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2,
                    dimnames = list(c("beta0", "beta1"), c("beta0", "beta1")))
    ),
    weights_info = list(values = c(0.2, 0.3, 0.5), trimmed_fraction = 0.25),
    sample = list(n_total = 100L, n_respondents = 60L, is_survey = FALSE, design = NULL),
    inference = list(variance_method = "delta", df = NA_real_, message = NA_character_),
    diagnostics = list(),
    meta = list(engine_name = "empirical_likelihood"),
    class = "nmar_result_el"
  )
  NMAR:::validate_nmar_result(res, "nmar_result_el")
}

test_that("vcov exposes all covariance components", {
  res <- make_rich_result()
# vcov now returns covariance of model coefficients (matches coef())
  beta <- coef(res)
  expect_true(is.numeric(beta) && length(beta) == 2)
  vc <- vcov(res)
  expect_true(is.matrix(vc))
  expect_equal(dim(vc), c(length(beta), length(beta)))
  expect_equal(rownames(vc), names(beta))
  expect_equal(colnames(vc), names(beta))
})

test_that("tidy/glance follow broom conventions and use canonical slots", {
  res <- make_rich_result()

  tidy_df <- generics::tidy(res, conf.level = 0.9)
  expect_true(all(c(
    "term", "estimate", "std.error", "statistic",
    "p.value", "conf.low", "conf.high", "component"
  ) %in% names(tidy_df)))
  expect_equal(tidy_df$component[1], "estimand")
  expect_true(all(tidy_df$component %in% c("estimand", "response")))
  resp_rows <- subset(tidy_df, component == "response")
  expect_equal(resp_rows$term, c("beta0", "beta1"))
  expect_true(all(is.finite(resp_rows$std.error)))
  expect_true(all(!is.na(resp_rows$conf.low)))

  glance_df <- generics::glance(res)
  expect_equal(glance_df$trimmed_fraction, 0.25)
  expect_equal(glance_df$estimate, 0.5)
})
