test_that("tidy/glance produce expected shapes", {
  set.seed(42)
  N <- 200
  X <- rnorm(N)
  Y <- 2 + 0.5 * X + rnorm(N)
  p <- plogis(-1 + 0.4 * scale(Y)[, 1])
  R <- runif(N) < p
  df <- data.frame(Y_miss = Y, X = X)
  df[!R, "Y_miss"] <- NA_real_
  eng <- el_engine(auxiliary_means = c(X = 0), variance_method = "delta")
  fml <- Y_miss ~ X
  fit <- nmar(formula = fml, data = df, engine = eng)

  td <- tidy(fit)
  gl <- glance(fit)
  expect_true(is.data.frame(td) && nrow(td) >= 1)
  expect_true(all(c("estimate", "std.error") %in% names(td)))
  expect_true(is.data.frame(gl) && nrow(gl) == 1)
  expect_true(all(c("estimate", "std.error", "converged") %in% names(gl)))
})

test_that("plot/autoplot run (skip ggplot2 if missing)", {
  set.seed(43)
  N <- 150
  X <- rnorm(N)
  Y <- 2 + 0.5 * X + rnorm(N)
  p <- plogis(-1 + 0.4 * scale(Y)[, 1])
  R <- runif(N) < p
  df <- data.frame(Y_miss = Y, X = X)
  df[!R, "Y_miss"] <- NA_real_
  eng <- el_engine(auxiliary_means = c(X = 0), variance_method = "delta")
  fml <- Y_miss ~ X
  fit <- nmar(formula = fml, data = df, engine = eng)

  # Open a temporary PDF device to avoid creating Rplots.pdf in the test dir
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(file = tmp)
  on.exit(
    {
      grDevices::dev.off()
      unlink(tmp)
    },
    add = TRUE
  )

  expect_error(suppressWarnings(plot(fit, which = "weights")), NA)
  expect_error(suppressWarnings(plot(fit, which = "fitted")), NA)

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    g1 <- autoplot(fit, type = "weights")
    g2 <- autoplot(fit, type = "fitted")
    expect_true(inherits(g1, "ggplot"))
    expect_true(inherits(g2, "ggplot"))
  } else {
    skip("ggplot2 not installed")
  }
})
