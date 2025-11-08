test_that("respondents-only data.frame with auxiliaries requires auxiliary_means", {
  df <- data.frame(Y_miss = rnorm(50), X = rnorm(50))
# No NAs in Y_miss => respondents-only
  expect_error(
    nmar(
      Y_miss ~ X,
      data = df,
      engine = el_engine(
        auxiliary_means = NULL,
        n_total = nrow(df),
        variance_method = "none"
      )
    ),
    "Respondents-only data detected.*auxiliary_means",
    fixed = FALSE
  )
# Provide auxiliary means -> should run
  expect_s3_class(
    nmar(
      Y_miss ~ X,
      data = df,
      engine = el_engine(
        auxiliary_means = c(X = mean(df$X)),
        n_total = nrow(df),
        variance_method = "none"
      )
    ),
    "nmar_result_el"
  )
})

test_that("respondents-only survey design with auxiliaries requires auxiliary_means", {
  skip_if_not_installed("survey")
  set.seed(1)
  df <- data.frame(Y_miss = rnorm(60), X = rnorm(60))
  des <- survey::svydesign(ids = ~1, weights = ~1, data = df)
  expect_error(
    nmar(
      Y_miss ~ X,
      data = des,
      engine = el_engine(
        auxiliary_means = NULL,
        n_total = sum(weights(des)),
        variance_method = "none"
      )
    ),
    "Respondents-only survey design.*auxiliary_means",
    fixed = FALSE
  )
  expect_s3_class(
    nmar(
      Y_miss ~ X,
      data = des,
      engine = el_engine(
        auxiliary_means = c(X = mean(df$X)),
        n_total = sum(weights(des)),
        variance_method = "none"
      )
    ),
    "nmar_result_el"
  )
})
