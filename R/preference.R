#' Build a physician-preference instrument from an order log
#'
#' For every encounter, the proportion of the provider's *other patients'*
#' orders in the window `[t - window, t)` that received the treatment of
#' interest:
#' * only orders strictly **before** the encounter time count (no same-day
#'   or same-encounter leakage);
#' * the **current patient's own orders are excluded** (leave-one-patient-out):
#'   otherwise the instrument partly contains the patient's own treatment
#'   decision and is no longer exogenous; with few patients per provider the
#'   instrument would almost equal the treatment;
#' * `n_orders` counts the orders used; set the preference to `NA` below
#'   `min_orders`;
#' * with `trend = TRUE`, also the change from the previous window
#'   `[t - 2 window, t - window)`: a (recent level, trend) pair can make
#'   "only the current instrument acts on the current treatment" more
#'   plausible than a long window (check with [lag_sufficiency_test()]).
#'
#' The result is continuous. [lii_tmle()] needs a binary instrument
#' (alpha-scaling of a probability); see [binarize_instrument()] and the
#' package README for what that changes in the interpretation.
#' @param orders `data.frame` with one row per order or encounter.
#' @param provider,patient,time,treated Column names: provider id, patient id,
#'   time (numeric or Date), and 0/1 indicator of the treatment of interest.
#' @param window Window length, in the units of `time` (days for a `Date`).
#' @param min_orders Minimum number of orders in the window.
#' @param trend Also return the change from the previous window.
#' @return `orders` with columns `pref`, `n_orders` and, if `trend`,
#'   `pref_prev` and `pref_trend`.
#' @examples
#' o <- data.frame(doc = c(1, 1, 1, 1, 2, 2), pat = c(1, 2, 3, 1, 4, 5),
#'                 t = c(1, 2, 3, 4, 1, 2), abx = c(1, 0, 1, 1, 0, 1))
#' preference_instrument(o, "doc", "pat", "t", "abx", window = 10, min_orders = 1)
#' @export
preference_instrument <- function(orders, provider, patient, time, treated, window,
                                  min_orders = 10L, trend = TRUE) {
  orders <- as.data.frame(orders)
  tt <- as.numeric(orders[[time]]); yy <- as.numeric(orders[[treated]])
  if (anyNA(tt) || anyNA(yy)) stop("`time` and `treated` must not contain NA.")
  win_sum <- function(groups, from, to) {
    ## sum and count of yy over rows of the same group with time in [t - from, t - to)
    S <- N <- numeric(length(tt))
    for (ix in split(seq_along(tt), groups, drop = TRUE)) {
      o <- ix[order(tt[ix])]; ts <- tt[o]; cs <- c(0, cumsum(yy[o]))
      hi <- findInterval(tt[ix] - to, ts, left.open = TRUE)   # times < t - to
      lo <- findInterval(tt[ix] - from, ts, left.open = TRUE) # times < t - from
      S[ix] <- cs[hi + 1] - cs[lo + 1]; N[ix] <- hi - lo
    }
    list(S = S, N = N)
  }
  grp <- orders[[provider]]; gp <- interaction(orders[[provider]], orders[[patient]], drop = TRUE)
  pref_at <- function(from, to) {
    a <- win_sum(grp, from, to); b <- win_sum(gp, from, to)
    N <- a$N - b$N
    list(p = ifelse(N >= min_orders, (a$S - b$S) / pmax(N, 1), NA_real_), N = N)
  }
  cur <- pref_at(window, 0)
  orders$pref <- cur$p; orders$n_orders <- cur$N
  if (trend) {
    prev <- pref_at(2 * window, window)
    orders$pref_prev <- prev$p
    orders$pref_trend <- orders$pref - orders$pref_prev
  }
  orders
}

#' Binarise a continuous instrument
#'
#' Turns a continuous instrument (e.g. a prescribing preference) into the
#' indicator `z > cut`, so that [lii_tmle()] can scale its probability. The
#' estimand then concerns being seen by a "high-preference" provider, and the
#' dichotomisation can create a direct path to the outcome through the
#' discarded variation (a violation of the exclusion restriction): report it
#' as a sensitivity analysis and compare several cut-offs.
#' @param z Numeric vector.
#' @param cut Threshold, or `"median"`.
#' @return Integer 0/1 vector (`NA` where `z` is `NA`).
#' @export
binarize_instrument <- function(z, cut = "median") {
  if (identical(cut, "median")) cut <- stats::median(z, na.rm = TRUE)
  as.integer(z > cut)
}
