test_that("lii_data builds the panel and handles censoring", {
  sim <- simulate_lii(n = 200, K = 2, seed = 1)
  d <- lii_data(sim$data, "id", "visit", "R", "A", "Y", covariates = "L", baseline = "W")
  expect_s3_class(d, "lii_data")
  expect_equal(d$n, 200); expect_equal(d$K, 2)
  nv <- table(sim$data$id)
  expect_equal(unname(rowSums(d$obs)), as.integer(pmin(nv, 3)))
  expect_true(all(d$y_obs <= (rowSums(d$obs) == 3)))
  expect_true(all(is.na(d$R[!d$obs])))
})

test_that("lii_data rejects non-binary instruments and NA covariates", {
  sim <- simulate_lii(n = 50, K = 1, seed = 2)
  bad <- sim$data; bad$R <- bad$R * 0.5
  expect_error(lii_data(bad, "id", "visit", "R", "A", "Y", covariates = "L"), "0/1")
  bad <- sim$data; bad$L[3] <- NA
  expect_error(lii_data(bad, "id", "visit", "R", "A", "Y", covariates = "L"), "NA")
})

test_that("continuous outcomes are rescaled and reported on their scale", {
  sim <- simulate_lii(n = 300, K = 1, seed = 3, censoring = FALSE)
  sim$data$Y <- 10 + 5 * sim$data$Y + ifelse(is.na(sim$data$Y), NA, stats::runif(nrow(sim$data)))
  d <- lii_data(sim$data, "id", "visit", "R", "A", "Y", covariates = "L", baseline = "W")
  expect_false(d$y_binary)
  fit <- lii_tmle(d, alphas = 1, learner = learner_glm, instrument_markov = TRUE, verbose = FALSE)
  expect_true(fit$estimates$estimate > 10 && fit$estimates$estimate < 16)
  expect_lt(abs(fit$estimates$estimate - mean(d$Y)), 0.5 * fit$estimates$se)
})
