#' Implied intervention on the treatment
#'
#' Under alpha-scaling of the instrument, the implied probability of
#' treatment at visit `j` given the covariate-treatment history is
#' \deqn{g_j(\alpha)(F_j) = \pi_0(F_j) + \alpha\lambda_j(F_j)\,[\pi_1(F_j) - \pi_0(F_j)],}
#' with \eqn{\pi_r(F_j) = P(A(j) = 1 \mid R(j) = r, F_j)}. The baseline
#' \eqn{\pi_0} (treatment without the instrument) is untouched: alpha moves
#' only the part of the treatment probability that responds to the
#' instrument. (Reading \eqn{\pi_0} as "always-takers" and
#' \eqn{\pi_1 - \pi_0} as "compliers" requires monotonicity, an extra
#' assumption; the formula itself does not.)
#'
#' This function averages \eqn{g_j(\alpha)} over the *observed* histories of
#' the subjects present at each visit. It describes how the treatment
#' mechanism changes; it is not the probability of treatment in the
#' intervened world, where histories also change. For the latter (for
#' example to calibrate alpha to "+10% treatment"), estimate the mean number
#' of treatments with [lii_treatment_outcome()] and [lii_tmle()], and use
#' [calibrate_alpha()].
#' @param x An object from [lii_data()].
#' @param nuisance An object from [lii_nuisance()].
#' @param alphas Scaling factors.
#' @return A `data.frame`: `alpha`, `visit`, `n`, `natural_p_treat` (mean of
#'   `A(j)`), `implied_p_treat`, `mean_lambda`, `mean_first_stage`
#'   (\eqn{\pi_1 - \pi_0}).
#' @export
lii_implied <- function(x, nuisance, alphas) {
  do.call(rbind, lapply(alphas, function(a) do.call(rbind, lapply(0:x$K, function(j) {
    o <- x$obs[, j + 1]
    lam <- nuisance$lambda[o, j + 1]; p0 <- nuisance$pi0[o, j + 1]; p1 <- nuisance$pi1[o, j + 1]
    data.frame(alpha = a, visit = j, n = sum(o), natural_p_treat = mean(x$A[o, j + 1]),
               implied_p_treat = mean(p0 + pmin(a * lam, 1) * (p1 - p0)),
               mean_lambda = mean(lam), mean_first_stage = mean(p1 - p0))
  }))))
}

#' Replace the outcome by the mean number of treatments
#'
#' Returns a copy of `x` whose outcome is the proportion of visits
#' `0..K` with `A(j) = 1` (for subjects observed through visit `K`).
#' [lii_tmle()] on it estimates how the expected amount of treatment changes
#' with alpha *in the intervened world*, which is the natural x-axis for
#' interpreting alpha (see [calibrate_alpha()]).
#' @param x An object from [lii_data()].
#' @return An object of class `lii_data`.
#' @export
lii_treatment_outcome <- function(x) {
  stopifnot(inherits(x, "lii_data"))
  full <- rowSums(x$obs) == x$K + 1
  x$Y <- ifelse(full, rowMeans(x$A), NA_real_)
  x$y_obs <- full
  x$y_range <- c(0, 1); x$y_binary <- FALSE
  x$names$outcome <- paste0("mean_of_", x$names$treatment)
  x
}

#' Choose alpha by the change in the amount of treatment
#'
#' Given a fit of [lii_tmle()] on [lii_treatment_outcome()], finds by linear
#' interpolation over the alpha grid the alpha whose expected amount of
#' treatment is `ratio` times the natural one (e.g. `ratio = 1.1` for 10% more
#' treatment). The result is a data-adaptive choice: inference for the
#' outcome at that alpha from [lii_tmle()] treats it as fixed.
#' @param fit_treatment Output of [lii_tmle()] on [lii_treatment_outcome()].
#' @param ratio Target ratio(s) of treatment relative to alpha = 1.
#' @return A `data.frame` with `ratio` and `alpha` (`NA` if outside the grid).
#' @export
calibrate_alpha <- function(fit_treatment, ratio) {
  e <- fit_treatment$estimates
  rel <- e$estimate / e$estimate[e$alpha == 1]
  if (any(diff(rel) <= 0)) warning("The treatment curve is not increasing in alpha on this grid.")
  data.frame(ratio = ratio,
             alpha = stats::approx(rel, e$alpha, xout = ratio, ties = mean)$y)
}
