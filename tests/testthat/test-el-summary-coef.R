test_that("coef(summary()) returns response-model table for EL (IID)", {
  df <- make_iid_nmar(n = 200, alpha = 0.5, seed = 6601)
  fit <- nmar(
    formula = Y_miss ~ X,
    data = df,
    engine = make_engine(auxiliary_means = c(X = 0), variance_method = "delta", standardize = FALSE)
  )
  expect_true(fit$converged)
  sm <- summary(fit)
  tb <- coef(sm)
  expect_true(is.null(tb) || is.data.frame(tb))
  if (is.data.frame(tb)) {
    nm <- names(tb)
    expect_true(all(c("Estimate", "Std. Error") %in% nm))
# IID path uses z-statistics
    expect_true(any(nm == "z value"))
    expect_true(any(nm == "Pr(>|z|)"))
    expect_true(nrow(tb) >= 1)
  }
  ci <- confint(sm)
  if (!is.null(ci)) {
    expect_true(is.matrix(ci))
    expect_true(nrow(ci) >= 1)
    expect_true(all(grepl("%$", colnames(ci))))
  }
})

test_that("coef(summary()) returns t-based labels for EL (survey)", {
  skip_if_not_installed("survey")
  set.seed(6602)
  N <- 120
  x <- stats::rchisq(N, df = 2)
  eps <- stats::rnorm(N)
  y <- x + eps * sqrt(x) / 5
  pr <- stats::plogis(0.2 * y - (-1.5))
  r <- stats::rbinom(N, 1, pr)
  df <- data.frame(y_miss = ifelse(r == 1, y, NA_real_), x = x)
  design <- survey::svydesign(ids = ~1, weights = ~1, data = df)
  fit <- nmar(
    formula = y_miss ~ x,
    data = design,
    engine = make_engine(auxiliary_means = c(x = mean(df$x)), variance_method = "delta", standardize = FALSE)
  )
  expect_true(fit$converged)
  sm <- summary(fit)
  tb <- coef(sm)
  if (is.data.frame(tb)) {
    nm <- names(tb)
    expect_true(any(nm == "t value"))
    expect_true(any(nm == "Pr(>|t|)"))
  }
})
