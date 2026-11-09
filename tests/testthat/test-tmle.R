test_that("alpha = 1 without censoring is close to the sample mean", {
  sim <- simulate_lii(n = 400, K = 2, seed = 4, censoring = FALSE)
  d <- lii_data(sim$data, "id", "visit", "R", "A", "Y", covariates = "L", baseline = "W")
  fit <- lii_tmle(d, alphas = c(0.5, 1, 1.5), learner = learner_glm, instrument_markov = TRUE,
                  verbose = FALSE)
  e <- fit$estimates
  expect_lt(abs(e$estimate[e$alpha == 1] - mean(d$Y)), 0.5 * e$se[e$alpha == 1])
  expect_equal(e$diff[e$alpha == 1], 0)
  expect_true(all(abs(colMeans(fit$ic)) < 1e-6))
})

test_that("estimates are close to the truth in a larger simulation (glm)", {
  sim <- simulate_lii(n = 3000, K = 2, seed = 5)
  d <- lii_data(sim$data, "id", "visit", "R", "A", "Y", covariates = "L", baseline = "W")
  fit <- lii_tmle(d, alphas = c(0.5, 1.5), learner = learner_glm, verbose = FALSE)
  tr <- sim$truth(fit$estimates$alpha, n_mc = 1e5)
  expect_true(all(abs(fit$estimates$estimate - tr) < 4 * fit$estimates$se))
})

test_that("the implied intervention leaves pi0 untouched and is affine in alpha", {
  sim <- simulate_lii(n = 300, K = 1, seed = 6)
  d <- lii_data(sim$data, "id", "visit", "R", "A", "Y", covariates = "L", baseline = "W")
  nu <- lii_nuisance(d, learner = learner_glm, verbose = FALSE)
  im <- lii_implied(d, nu, c(0, 0.5, 1))
  v0 <- im[im$visit == 0, ]
  o <- d$obs[, 1]
  expect_equal(v0$implied_p_treat[v0$alpha == 0], mean(nu$pi0[o, 1]))
  expect_equal(diff(diff(v0$implied_p_treat)), 0, tolerance = 1e-12)
})

test_that("treatment-count outcome and calibration run", {
  sim <- simulate_lii(n = 400, K = 2, seed = 7)
  d <- lii_data(sim$data, "id", "visit", "R", "A", "Y", covariates = "L", baseline = "W")
  dt <- lii_treatment_outcome(d)
  ft <- lii_tmle(dt, alphas = c(0.5, 1.5), learner = learner_glm, verbose = FALSE)
  expect_true(all(diff(ft$estimates$estimate) > 0))
  ca <- calibrate_alpha(ft, 1.1)
  expect_true(ca$alpha > 1 && ca$alpha < 1.5)
})

test_that("bootstrap returns finite SEs", {
  sim <- simulate_lii(n = 200, K = 1, seed = 8)
  d <- lii_data(sim$data, "id", "visit", "R", "A", "Y", covariates = "L", baseline = "W")
  b <- lii_bootstrap(d, alphas = c(0.5, 1), learner = learner_glm, B = 20, verbose = FALSE)
  expect_true(all(is.finite(b$estimates$boot_se)))
  expect_equal(b$estimates$diff_boot_se[b$estimates$alpha == 1], 0)
})
