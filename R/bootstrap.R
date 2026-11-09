#' Nonparametric bootstrap for the fixed parameter
#'
#' Resamples subjects with replacement and reruns the whole pipeline,
#' refitting every nuisance, including \eqn{\hat\lambda}, in each resample.
#' This targets \eqn{\psi(\alpha)} defined with the *true* instrument
#' intensity \eqn{\lambda}, whose influence curve has an extra term (the
#' intervention itself depends on the data distribution); the
#' influence-curve standard errors of [lii_tmle()] omit it. The bootstrap is
#' expected to be valid with HAL-type learners (hapc); it is slow: `B` full
#' fits.
#' @param x An object from [lii_data()].
#' @param alphas,learner,learner_args,lags,instrument_markov,gbound Passed to
#'   [lii_tmle()].
#' @param B Number of resamples.
#' @param seed Seed.
#' @param verbose Print progress every 10 resamples.
#' @return A list: `estimates` (alpha, point estimate on the full data,
#'   bootstrap SE, 2.5% and 97.5% percentiles, same for the contrast with
#'   alpha = 1) and `draws` (B x alphas matrix).
#' @export
lii_bootstrap <- function(x, alphas = c(0.5, 0.75, 1, 1.25, 1.5), learner = learner_hapc,
                          learner_args = list(), lags = 1L, instrument_markov = FALSE,
                          gbound = 0.01, B = 200L, seed = 1L, verbose = TRUE) {
  full <- lii_tmle(x, alphas, learner, learner_args, lags, instrument_markov,
                   gbound = gbound, verbose = FALSE)
  al <- full$estimates$alpha
  set.seed(seed)
  draws <- t(vapply(seq_len(B), function(b) {
    if (verbose && b %% 10 == 0) message("bootstrap ", b, "/", B)
    xb <- .subset_lii(x, sample.int(x$n, replace = TRUE))
    tryCatch(lii_tmle(xb, al, learner, learner_args, lags, instrument_markov,
                      gbound = gbound, verbose = FALSE)$estimates$estimate,
             error = function(e) rep(NA_real_, length(al)))
  }, numeric(length(al))))
  if (length(al) == 1) draws <- matrix(draws, ncol = 1)
  colnames(draws) <- paste0("alpha_", al)
  ref <- which(al == 1)
  dd <- draws - draws[, ref]
  qs <- function(M, p) apply(M, 2, stats::quantile, p, na.rm = TRUE)
  list(estimates = data.frame(
    alpha = al, estimate = full$estimates$estimate, boot_se = apply(draws, 2, stats::sd, na.rm = TRUE),
    lower = qs(draws, 0.025), upper = qs(draws, 0.975),
    diff = full$estimates$diff, diff_boot_se = apply(dd, 2, stats::sd, na.rm = TRUE),
    diff_lower = qs(dd, 0.025), diff_upper = qs(dd, 0.975),
    failed = colSums(is.na(draws)), row.names = NULL),
    draws = draws)
}
