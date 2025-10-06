test_that("variance_method = 'none' skips variance for IID and survey", {
  cases <- list(
    list(kind = "iid", builder = function() {
      b <- make_iid_nmar(n = 50, alpha = 0.4, seed = 1)
      list(data = data.frame(Y = b$Y_miss, X = b$X), formula = Y ~ X, is_survey = FALSE)
    }),
    list(kind = "survey", builder = function() {
      skip_if_not_installed("survey")
      set.seed(2)
      b <- make_iid_nmar(n = 60, alpha = 0.4, seed = 2)
      df <- data.frame(Y = b$Y_miss, X = b$X, w = runif(nrow(b), 0.5, 2))
      des <- survey::svydesign(ids = ~1, weights = ~w, data = df)
      list(data = des, formula = Y ~ X, is_survey = TRUE)
    })
  )

  for (cs in cases) {
    built <- cs$builder()
    if (cs$kind == "iid") {
      fit <- NMAR:::el.data.frame(built$data, built$formula,
                                  standardize = TRUE,
                                  auxiliary_means = NULL,
                                  variance_method = "none",
                                  bootstrap_reps = 10)
    } else {
      fit <- suppressWarnings(NMAR:::el.survey.design(built$data, built$formula,
                                  standardize = TRUE,
                                  auxiliary_means = NULL,
                                  variance_method = "none",
                                  bootstrap_reps = 10))
    }
    expect_true(is.na(fit$se))
    expect_true(all(is.na(fit$model$vcov)))
    expect_match(fit$diagnostics$vcov_message, "Variance skipped")
  }
})
