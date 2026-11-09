#' Fit the nuisance models that do not depend on alpha
#'
#' Pooled over visits (the visit index is a regressor), on the rows of
#' subjects observed at each visit:
#' * `g_R`: \eqn{P(R(j) = 1 \mid} full past, including past instruments\eqn{)},
#'   the propensity of the instrument used in the inverse weights;
#' * `lambda`: \eqn{\lambda_j(F_j) = P(R(j) = 1 \mid F_j)}, the instrument
#'   intensity given the covariate-treatment history only (no past
#'   instruments). The intervention is \eqn{h^*_j(1 \mid F_j) =
#'   \alpha\lambda_j(F_j)}. With `instrument_markov = TRUE` the natural
#'   instrument law is assumed to depend on the past only through \eqn{F_j},
#'   so one model serves both purposes;
#' * `g_C`: probability of remaining observed at the next visit (or of
#'   having the outcome observed, after visit `K`), given the full past and
#'   `R(j)`, `A(j)`. Skipped when nobody is censored;
#' * `pi0`, `pi1`: \eqn{\pi_r(F_j) = P(A(j) = 1 \mid R(j) = r, F_j)}, fitted
#'   separately in each instrument arm. They give the implied intervention on
#'   the treatment, \eqn{g_j(\alpha) = \pi_0 + \alpha\lambda_j(\pi_1 - \pi_0)},
#'   and the first-stage strength \eqn{\pi_1 - \pi_0}.
#'
#' Usually called by [lii_tmle()]; call it directly to reuse the fits across
#' several calls (argument `nuisance` of [lii_tmle()]).
#' @param x An object from [lii_data()].
#' @param learner,learner_args Learner function and extra arguments (see
#'   [learners]).
#' @param lags Number of lags of `L`, `A` and `R` in the history summaries
#'   (counts of past `A` and `R` are always included).
#' @param instrument_markov Logical; see above.
#' @param verbose Print progress.
#' @return An object of class `lii_nuisance`: matrices (subjects x visits) of
#'   predictions, `NA` where a subject is not observed.
#' @export
lii_nuisance <- function(x, learner = learner_hapc, learner_args = list(), lags = 1L,
                         instrument_markov = FALSE, verbose = TRUE) {
  stopifnot(inherits(x, "lii_data"))
  n <- x$n; K <- x$K
  empty <- function() matrix(NA_real_, n, K + 1)
  put <- function(M, P, v) { M[cbind(P$i, P$j + 1)] <- v; M }
  say <- function(...) if (verbose) message(...)

  ## instrument models
  Pf <- .pooled(x, lags, full = FALSE)
  yR <- x$R[cbind(Pf$i, Pf$j + 1)]
  say("Fitting lambda = P(R(j)=1 | F_j) on ", nrow(Pf$X), " visit rows")
  lam <- put(empty(), Pf, .fit_predict(learner, learner_args, Pf$X, yR, list(Pf$X), "lambda")[[1]])
  if (instrument_markov) {
    gR <- lam
  } else {
    Pa <- .pooled(x, lags, full = TRUE)
    say("Fitting g_R = P(R(j)=1 | full past)")
    gR <- put(empty(), Pa, .fit_predict(learner, learner_args, Pa$X, yR, list(Pa$X), "g_R")[[1]])
  }

  ## censoring
  nxt <- cbind(x$obs[, -1, drop = FALSE], x$y_obs)
  if (all(nxt[x$obs])) {
    gC <- ifelse(x$obs, 1, NA_real_)
  } else {
    RA <- function(j) cbind(R_now = x$R[, j + 1], A_now = x$A[, j + 1])
    Pc <- .pooled(x, lags, full = TRUE, extra = RA)
    yC <- as.numeric(nxt[cbind(Pc$i, Pc$j + 1)])
    say("Fitting g_C = P(observed at next visit | past, R(j), A(j))")
    gC <- put(empty(), Pc, .fit_predict(learner, learner_args, Pc$X, yC, list(Pc$X), "g_C")[[1]])
  }

  ## first stage, by instrument arm
  yA <- x$A[cbind(Pf$i, Pf$j + 1)]
  pis <- lapply(0:1, function(r) {
    s <- yR == r
    if (sum(s) < 20) {
      Xr <- cbind(Pf$X, R_now = yR)
      nd <- cbind(Pf$X, R_now = r)
      p <- .fit_predict(learner, learner_args, Xr, yA, list(nd), paste0("pi", r, "_joint"))[[1]]
    } else {
      p <- .fit_predict(learner, learner_args, Pf$X[s, , drop = FALSE], yA[s], list(Pf$X), paste0("pi", r))[[1]]
    }
    put(empty(), Pf, p)
  })
  say("Fitted pi_0, pi_1 = P(A(j)=1 | R(j)=r, F_j)")
  structure(list(g_R = gR, lambda = lam, g_C = gC, pi0 = pis[[1]], pi1 = pis[[2]],
                 lags = lags, instrument_markov = instrument_markov),
            class = "lii_nuisance")
}
