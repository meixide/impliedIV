test_that("preference_instrument matches a brute-force computation", {
  set.seed(1)
  o <- data.frame(doc = sample(1:3, 300, TRUE), pat = sample(1:40, 300, TRUE),
                  t = sample(1:200, 300, TRUE), abx = rbinom(300, 1, 0.4))
  p <- preference_instrument(o, "doc", "pat", "t", "abx", window = 30, min_orders = 1)
  brute <- function(i, from, to) {
    s <- o$doc == o$doc[i] & o$pat != o$pat[i] & o$t < o$t[i] - to & o$t >= o$t[i] - from
    c(if (sum(s)) mean(o$abx[s]) else NA, sum(s))
  }
  b <- t(sapply(seq_len(nrow(o)), brute, from = 30, to = 0))
  expect_equal(p$n_orders, b[, 2])
  expect_equal(p$pref, b[, 1])
  bp <- t(sapply(seq_len(nrow(o)), brute, from = 60, to = 30))
  expect_equal(p$pref_prev, bp[, 1])
})

test_that("lag_sufficiency_test runs and does not reject when the assumption holds", {
  sim <- simulate_lii(n = 1500, K = 3, seed = 9, censoring = FALSE)
  lt <- lag_sufficiency_test(sim$data, "id", "visit", "R", "A", covariates = "L", lags = 2)
  expect_equal(lt$df, 2)
  expect_true(lt$p_value > 0.001)
  expect_true(lt$current_term$estimate > 1)
})

test_that("binarize_instrument splits at the median", {
  expect_equal(binarize_instrument(c(0.1, 0.2, 0.3, 0.4)), c(0L, 0L, 1L, 1L))
})
