#' Diagnostics of a fit
#'
#' Read these before the estimates.
#' * `instrument`: distribution of the fitted instrument probabilities
#'   (`lambda`, `g_R`) and the largest alpha that never needs capping,
#'   `1 / max(lambda)`.
#' * `first_stage`: distribution of \eqn{\pi_1 - \pi_0}, the change in the
#'   probability of treatment when the instrument fires. Values near 0
#'   everywhere mean a weak instrument: alpha then barely moves the
#'   treatment and the curve in alpha will be flat whatever the effect of the
#'   treatment.
#' * `weights`: for each alpha, the smallest effective sample size of the
#'   cumulative weights over visits and the largest weight.
#' * `censoring`: distribution of the probability of remaining observed.
#' * `fits`: the learner log ([lii_fit_log()]).
#' @param fit An object from [lii_tmle()].
#' @return A list of `data.frame`s.
#' @export
lii_diagnostics <- function(fit) {
  nu <- fit$nuisance
  q <- function(v) { v <- v[!is.na(v)]
    data.frame(min = min(v), q05 = stats::quantile(v, 0.05), median = stats::median(v),
               mean = mean(v), q95 = stats::quantile(v, 0.95), max = max(v), row.names = NULL) }
  fs <- nu$pi1 - nu$pi0
  list(
    instrument = rbind(cbind(quantity = "lambda", q(nu$lambda)),
                       cbind(quantity = "g_R", q(nu$g_R))),
    alpha_max_uncapped = 1 / max(nu$lambda, na.rm = TRUE),
    first_stage = cbind(q(fs), share_below_0.05 = mean(fs[!is.na(fs)] < 0.05)),
    weights = fit$estimates[, c("alpha", "min_n_eff", "max_weight", "p_capped")],
    censoring = q(nu$g_C),
    fits = lii_fit_log())
}

#' Test that only the current instrument acts on the current treatment
#'
#' The identification of the implied intervention uses
#' \eqn{A(j) \perp (R(j-1), \dots, R(j-p)) \mid R(j), F_j}: once the current
#' instrument and the covariate-treatment history are known, past
#' instruments carry no further information about the treatment. This
#' function fits the pooled logistic regression of `A(j)` on the current
#' instrument, `lags` past values of it, the current covariates, the
#' previous treatment and the visit index, and tests that the coefficients of
#' the past instruments are jointly zero (Wald test with a variance robust to
#' clustering within subject).
#'
#' It works on the raw long data, so the instrument may be continuous (for
#' example a physician's prescribing preference). With a windowed preference,
#' consecutive windows overlap, and a rejection says that the older window
#' still pulls on today's decision beyond the current one: consider a shorter
#' window or adding a trend term to the current instrument (see
#' [preference_instrument()]). Neither AUROC stability across windows nor an
#' autoregressive model for the instrument tests this assumption. A
#' non-rejection is not proof; the test is only as good as the covariates.
#' @param data Long `data.frame`, one row per subject and visit.
#' @param id,time,instrument,treatment Column names.
#' @param covariates Covariate columns (current visit) to adjust for.
#' @param lags Number of past instrument values tested.
#' @return A list with the coefficient table of the lag terms, the Wald
#'   statistic, degrees of freedom, p-value and the number of rows used.
#' @export
lag_sufficiency_test <- function(data, id, time, instrument, treatment, covariates = NULL, lags = 2L) {
  data <- as.data.frame(data)
  data <- data[order(data[[id]], data[[time]]), , drop = FALSE]
  g <- data[[id]]
  lagv <- function(v, l) stats::ave(v, g, FUN = function(z) c(rep(NA, l), utils::head(z, -l)))
  X <- data.frame(Z_now = data[[instrument]])
  for (l in seq_len(lags)) X[[paste0("Z_lag", l)]] <- lagv(data[[instrument]], l)
  X$A_lag1 <- lagv(data[[treatment]], 1)
  X$visit <- stats::ave(seq_len(nrow(data)), g, FUN = seq_along)
  if (length(covariates)) X <- cbind(X, data[, covariates, drop = FALSE])
  ok <- stats::complete.cases(X) & !is.na(data[[treatment]])
  if (sum(ok) < 50) stop("Too few rows with ", lags, " past instrument values.")
  Xm <- stats::model.matrix(~ ., data = X[ok, , drop = FALSE])
  y <- data[[treatment]][ok]; cl <- g[ok]
  f <- stats::glm.fit(Xm, y, family = stats::binomial())
  b <- f$coefficients; keep <- !is.na(b)
  Xm <- Xm[, keep, drop = FALSE]; b <- b[keep]
  mu <- as.numeric(stats::plogis(Xm %*% b))
  bread <- solve(crossprod(Xm * sqrt(mu * (1 - mu))))
  U <- rowsum(Xm * (y - mu), cl)
  V <- bread %*% crossprod(U) %*% bread
  lag_names <- intersect(paste0("Z_lag", seq_len(lags)), names(b))
  bl <- b[lag_names]; Vl <- V[lag_names, lag_names, drop = FALSE]
  wald <- as.numeric(t(bl) %*% solve(Vl, bl))
  list(lag_terms = data.frame(term = lag_names, estimate = bl, robust_se = sqrt(diag(Vl)),
                              row.names = NULL),
       current_term = data.frame(term = "Z_now", estimate = b["Z_now"],
                                 robust_se = sqrt(V["Z_now", "Z_now"]), row.names = NULL),
       wald = wald, df = length(bl), p_value = stats::pchisq(wald, length(bl), lower.tail = FALSE),
       n_rows = sum(ok))
}
