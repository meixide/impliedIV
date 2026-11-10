#' Comparator that ignores the instrument (for simulations and sensitivity)
#'
#' Estimates the mean outcome under the stochastic intervention that sets the
#' probability of *treatment* at each visit to the implied
#' \eqn{g_j(\alpha)(F_j) = \pi_0 + \alpha\lambda_j(\pi_1 - \pi_0)}, as if the
#' treatment were unconfounded given the covariate-treatment history
#' \eqn{F_j} (a standard longitudinal TMLE for a stochastic intervention on
#' the treatment, with the instrument left out of the analysis). This is what
#' an analyst without the instrument would compute for the same change in the
#' treatment mechanism. It is biased when unmeasured factors drive both
#' treatment and outcome; the difference with [lii_tmle()] is a diagnostic
#' of such confounding. At \eqn{\alpha = 1} it estimates the observed mean.
#' @param x An object from [lii_data()].
#' @param alphas Scaling factors (1 is added).
#' @param nuisance Output of [lii_nuisance()] for `x` (it provides the implied
#'   treatment probabilities).
#' @param learner,learner_args,gbound,min_arm,verbose As in [lii_tmle()].
#' @return A `data.frame`: `alpha`, `estimate`, `se`, `lower`, `upper`, and the
#'   contrast with alpha = 1 (`diff`, `diff_se`, `diff_lower`, `diff_upper`).
#' @export
lii_naive <- function(x, alphas, nuisance, learner = learner_hapc, learner_args = list(),
                      gbound = 0.01, min_arm = 30L, verbose = TRUE) {
  stopifnot(inherits(x, "lii_data"), inherits(nuisance, "lii_nuisance"))
  alphas <- sort(unique(c(alphas, 1)))
  xn <- x
  xn$R <- x$A                      # the intervention node is now the treatment
  xn$A <- ifelse(is.na(x$A), NA, 0)
  nun <- lii_nuisance(xn, learner, learner_args, nuisance$lags, instrument_markov = FALSE,
                      verbose = FALSE, pooled = isTRUE(nuisance$pooled), parts = c("g_R", "g_C"))
  h1 <- lapply(alphas, function(a)
    nuisance$pi0 + pmin(a * nuisance$lambda, 1) * (nuisance$pi1 - nuisance$pi0))
  res <- .tmle_core(xn, nun, NULL, paste0("naive_alpha", alphas), learner, learner_args,
                    gbound, min_arm, verbose, h1_list = h1)
  ic <- sapply(res, `[[`, "ic"); est <- vapply(res, `[[`, 0, "psi")
  ref <- which(alphas == 1); z <- stats::qnorm(0.975)
  se <- apply(ic, 2, stats::sd) / sqrt(x$n); dse <- apply(ic - ic[, ref], 2, stats::sd) / sqrt(x$n)
  data.frame(alpha = alphas, estimate = est, se = se, lower = est - z * se, upper = est + z * se,
             diff = est - est[ref], diff_se = dse, diff_lower = est - est[ref] - z * dse,
             diff_upper = est - est[ref] + z * dse)
}
