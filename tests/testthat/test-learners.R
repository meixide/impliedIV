test_that("learner_hapc calls cv.hapc with binomial, norm = '1', and handles soft labels", {
  skip_if_not_installed("hapc"); skip_if_not_installed("glmnet")
  set.seed(1); n <- 300; X <- matrix(runif(n * 2), n, dimnames = list(NULL, c("a", "b")))
  p <- plogis(sin(pi * X[, 1]) - 1 + X[, 2])
  lii_reset_log()
  for (Y in list(rbinom(n, 1, p), p)) {
    q <- learner_hapc(X, Y, X, label = "t")
    expect_length(q, n); expect_true(all(q > 0 & q < 1))
    expect_gt(cor(q, p), 0.7)
  }
  lg <- lii_fit_log()
  expect_equal(lg$learner, c("hapc", "hapc"))
  expect_true(all(lg$lambda_position %in% c("interior", "lower_edge", "upper_edge", "null_model")))
  ## the call really uses norm = "1": its argument list is fixed in the source
  src <- paste(deparse(body(learner_hapc)), collapse = "")
  expect_match(src, 'family = "binomial", norm = "1"', fixed = TRUE)
})
