test_that("family score_eta equals d/deta log p for logit and probit", {
  set.seed(1001)
  eta <- c(-10, -5, -2, 0, 2, 5, 10)

  # Logit: score should be 1 - p
  flogit <- nmar:::logit_family()
  p_logit <- flogit$linkinv(eta)
  s_logit <- flogit$score_eta(eta, 1)
  expect_equal(s_logit, 1 - p_logit, tolerance = 1e-12)

  # Probit: score should be phi / Phi (compute stably via logs)
  fprobit <- nmar:::probit_family()
  p_probit <- fprobit$linkinv(eta)
  s_probit <- fprobit$score_eta(eta, 1)
  expected <- exp(dnorm(eta, log = TRUE) - pnorm(eta, log.p = TRUE))
  expect_equal(s_probit, expected, tolerance = 1e-6)
})
