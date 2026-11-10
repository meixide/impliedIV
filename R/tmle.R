#' Mean outcome under alpha-scaling of the instrument (longitudinal TMLE)
#'
#' Estimates \eqn{\psi(\alpha) = E[Y^{h^*_\alpha}]}: the mean outcome, had
#' nobody been censored, under the stochastic intervention that sets the
#' probability of the instrument at every visit to
#' \eqn{h^*_j(1 \mid F_j) = \alpha\,\lambda_j(F_j)}, where
#' \eqn{\lambda_j(F_j) = P(R(j) = 1 \mid F_j)} is its natural probability
#' given the covariate-treatment history. \eqn{\alpha = 1} is the natural
#' world when the instrument's law depends on the past only through
#' \eqn{F_j} (otherwise it already removes the dependence on past
#' instruments; contrasts are always with \eqn{\alpha = 1});
#' \eqn{\alpha = 1.2} means "the instrument fires 20% more often than
#' it does now, for every type of patient". The intervention never asks for
#' the instrument where it never occurs, so there are no positivity problems
#' for \eqn{\alpha \le 1}; for \eqn{\alpha > 1} the probability is capped at 1
#' where \eqn{\alpha\lambda_j > 1} (the share of capped visits is reported).
#'
#' **What is identified.** Only the instrument is intervened on; the treatment
#' and covariates follow their natural (possibly confounded) mechanisms. The
#' estimand needs sequential randomisation of the instrument given the
#' observed past and positivity of `g_R`; it does *not* need the treatment to
#' be unconfounded. Under the exclusion restriction (the instrument affects
#' the outcome only through the treatment) and the assumption that the
#' instrument acts only on the current treatment, \eqn{\psi(\alpha)} is also
#' the effect of the *implied* stochastic intervention on the treatment,
#' reported in `implied` (see [lii_implied()]).
#'
#' **Estimator.** Sequential regression (iterated conditional expectations)
#' with a logistic targeting step at every visit, on subjects still
#' observed. The outcome regressions are fitted separately in each instrument
#' arm. The intervention uses the estimated \eqn{\hat\lambda}; standard
#' errors come from the influence curve and are valid for the
#' *data-adaptive* parameter that treats \eqn{\alpha\hat\lambda} as given
#' (Hubbard, Kherad-Pajouh and van der Laan, 2016).
#'
#' **Fixed parameter.** With `fixed = TRUE` (default) the output also contains
#' the one-step estimator and standard errors for the parameter defined with
#' the *true* \eqn{\lambda} (columns ending in `_fix`). The efficient
#' influence curve of that parameter has one more term, because the
#' intervention depends on the data distribution:
#' \eqn{\sum_j \alpha\, b_j(F_j)\{R(j) - \lambda_j(F_j)\}} with
#' \eqn{b_j(F_j) = E[W_{j-1}\{Q_j(1, H_j) - Q_j(0, H_j)\} \mid F_j]}, where
#' \eqn{W_{j-1}} is the cumulative weight before visit \eqn{j} and
#' \eqn{Q_j(r, H_j)} the sequential regression at instrument value \eqn{r}.
#' It costs one more regression per visit and per alpha. [lii_bootstrap()]
#' is the resampling alternative.
#'
#' @param x An object from [lii_data()].
#' @param alphas Scaling factors. 1 is always added (it is the reference of
#'   the contrasts).
#' @param learner,learner_args Learner for all regressions and its extra
#'   arguments (see [learners]). Default `learner_hapc`.
#' @param lags,instrument_markov Passed to [lii_nuisance()].
#' @param nuisance Optional output of [lii_nuisance()] to reuse.
#' @param gbound Lower bound for `g_R`, `1 - g_R` and `g_C` in the weights.
#' @param min_arm Smallest number of rows to fit an outcome regression within
#'   one instrument arm; below it, one regression with the instrument as a
#'   regressor is used.
#' @param verbose Print progress.
#' @param fixed Also compute the estimator and standard errors for the
#'   parameter defined with the true instrument intensity (see Details).
#' @param pooled Passed to [lii_nuisance()] when `nuisance` is `NULL`.
#' @return An object of class `lii_fit` with
#'   * `estimates`: one row per alpha: `estimate`, `se`, 95% `lower`/`upper`;
#'     the contrast with alpha = 1 (`diff`, `diff_se`, `diff_lower`,
#'     `diff_upper`); `p_capped` (share of visits where
#'     \eqn{\alpha\hat\lambda > 1}); `min_n_eff` and `max_weight` of the
#'     cumulative weights over visits;
#'   * `implied`: the implied treatment probabilities by visit (see
#'     [lii_implied()]);
#'     with `fixed = TRUE`, the same quantities for the fixed parameter in
#'     `estimate_fix`, `se_fix`, `lower_fix`, `upper_fix`, `diff_fix`,
#'     `diff_se_fix`, `diff_lower_fix`, `diff_upper_fix`;
#'   * `ic`: influence-curve matrix (subjects x alphas), on the outcome scale;
#'   * `ic_lambda`: the additional influence-curve term of the fixed parameter;
#'   * `nuisance`, `settings`.
#' @examples
#' sim <- simulate_lii(n = 300, K = 2, seed = 2)
#' d <- lii_data(sim$data, "id", "visit", "R", "A", "Y", covariates = "L", baseline = "W")
#' fit <- lii_tmle(d, alphas = c(0.5, 1.5), learner = learner_glm, verbose = FALSE)
#' fit
#' @export
lii_tmle <- function(x, alphas = c(0.5, 0.75, 1, 1.25, 1.5), learner = learner_hapc,
                     learner_args = list(), lags = 1L, instrument_markov = FALSE,
                     nuisance = NULL, gbound = 0.01, min_arm = 30L, verbose = TRUE,
                     fixed = TRUE, pooled = TRUE) {
  stopifnot(inherits(x, "lii_data"))
  if (any(alphas < 0)) stop("alphas must be >= 0.")
  alphas <- sort(unique(c(alphas, 1)))
  if (is.null(nuisance)) nuisance <- lii_nuisance(x, learner, learner_args, lags,
                                                   instrument_markov, verbose, pooled = pooled)
  n <- x$n; K <- x$K
  res <- .tmle_core(x, nuisance, lapply(alphas, rep, x$K + 1), paste0("alpha", alphas),
                    learner, learner_args, gbound, min_arm, verbose, fixed = fixed)

  ic <- sapply(res, `[[`, "ic"); if (!is.matrix(ic)) ic <- matrix(ic, ncol = 1)
  colnames(ic) <- paste0("alpha_", alphas)
  est <- vapply(res, `[[`, 0, "psi")
  se <- apply(ic, 2, stats::sd) / sqrt(n)
  ref <- which(alphas == 1)
  dse <- apply(ic - ic[, ref], 2, stats::sd) / sqrt(n)
  z <- stats::qnorm(0.975)
  estimates <- data.frame(
    alpha = alphas, estimate = est, se = se, lower = est - z * se, upper = est + z * se,
    diff = est - est[ref], diff_se = dse, diff_lower = est - est[ref] - z * dse,
    diff_upper = est - est[ref] + z * dse,
    p_capped = vapply(res, `[[`, 0, "p_capped"),
    min_n_eff = vapply(res, `[[`, 0, "min_n_eff"),
    max_weight = vapply(res, `[[`, 0, "max_weight"))
  ic_lambda <- NULL
  if (fixed) {
    ic_lambda <- sapply(res, `[[`, "ic_lambda"); if (!is.matrix(ic_lambda)) ic_lambda <- matrix(ic_lambda, ncol = 1)
    colnames(ic_lambda) <- colnames(ic)
    full <- ic + ic_lambda
    ef <- est + colMeans(ic_lambda)
    sf <- apply(full, 2, stats::sd) / sqrt(n)
    dsf <- apply(full - full[, ref], 2, stats::sd) / sqrt(n)
    estimates <- cbind(estimates, estimate_fix = ef, se_fix = sf, lower_fix = ef - z * sf,
                       upper_fix = ef + z * sf, diff_fix = ef - ef[ref], diff_se_fix = dsf,
                       diff_lower_fix = ef - ef[ref] - z * dsf, diff_upper_fix = ef - ef[ref] + z * dsf)
  }
  structure(list(estimates = estimates, implied = lii_implied(x, nuisance, alphas),
                 ic = ic, ic_lambda = ic_lambda, nuisance = nuisance,
                 settings = list(n = n, K = K, alphas = alphas, gbound = gbound,
                                 lags = nuisance$lags, learner = deparse(substitute(learner)),
                                 outcome = x$names$outcome)),
            class = "lii_fit")
}

## LTMLE for a list of per-visit scaling vectors `avecs` (each of length K + 1)
## `h1_list` (optional): list of n x (K + 1) matrices of intervention probabilities, used instead
## of alpha * lambda (then `avecs` is ignored). `fixed = TRUE` also returns the component of the
## efficient influence curve that accounts for the estimation of lambda,
## sum_j alpha_j b_j(F_j) (R(j) - lambda_j(F_j)), b_j(F_j) = E[W_{j-1} {Q_j(1, H_j) - Q_j(0, H_j)} | F_j].
.tmle_core <- function(x, nuisance, avecs, labels, learner, learner_args, gbound, min_arm, verbose,
                       fixed = FALSE, h1_list = NULL) {
  n <- x$n; K <- x$K
  lo <- x$y_range[1]; sc <- diff(x$y_range)
  Ys <- (x$Y - lo) / sc
  gR <- pmin(pmax(nuisance$g_R, gbound), 1 - gbound)
  gC <- pmax(nuisance$g_C, gbound)
  lam <- pmin(pmax(nuisance$lambda, 0), 1)
  nxt <- cbind(x$obs[, -1, drop = FALSE], x$y_obs)
  Hist <- lapply(0:K, function(j) .history(x, j, nuisance$lags, full = TRUE))
  cache_K <- NULL

  if (!is.null(h1_list)) { avecs <- h1_list; fixed <- FALSE }
  res <- lapply(seq_along(avecs), function(k) {
    a <- avecs[[k]]
    if (verbose) message("intervention ", labels[k])
    h1 <- if (is.null(h1_list)) pmin(sweep(lam, 2, a, `*`), 1) else pmin(pmax(a, 0), 1)
    Dl <- rep(0, n)
    ratio <- ifelse(x$R == 1, h1 / gR, (1 - h1) / (1 - gR)) / gC
    ratio[!x$obs] <- NA
    Hcum <- t(apply(ratio, 1, function(r) cumprod(ifelse(is.na(r), 0, r))))
    if (K == 0) Hcum <- matrix(Hcum, ncol = 1)
    Qnext <- Ys
    D <- rep(0, n); neff <- maxw <- numeric(K + 1)
    for (j in K:0) {
      o <- x$obs[, j + 1]; f <- nxt[, j + 1] & o
      Rj <- x$R[, j + 1]
      if (j == K && !is.null(cache_K)) {
        Qb <- cache_K
      } else {
        Qb <- .q_arms(learner, learner_args, Hist[[j + 1]], Rj, Qnext, o, f, min_arm,
                      paste0("Q", j, if (j < K) paste0("_", labels[k])))
        if (j == K) cache_K <<- Qb
      }
      QbR <- ifelse(Rj == 1, Qb[, 2], Qb[, 1])
      w <- Hcum[f, j + 1]
      eps <- .fluctuate(Qnext[f], stats::qlogis(QbR[f]), w)
      Qs <- stats::plogis(stats::qlogis(Qb) + eps)
      QsR <- ifelse(Rj == 1, Qs[, 2], Qs[, 1])
      D[f] <- D[f] + w * (Qnext[f] - QsR[f])
      neff[j + 1] <- if (sum(w) > 0) sum(w)^2 / sum(w^2) else 0
      maxw[j + 1] <- max(w, 0)
      if (fixed) {
        Wprev <- if (j == 0) rep(1, n) else Hcum[, j]
        Fh <- .history(x, j, nuisance$lags, full = FALSE)[o, , drop = FALSE]
        b <- .fit_real(learner, learner_args, Fh, Wprev[o] * (Qs[o, 2] - Qs[o, 1]),
                       paste0("b", j, "_", labels[k]))
        open <- a[j + 1] * nuisance$lambda[o, j + 1] < 1
        Dl[o] <- Dl[o] + a[j + 1] * open * b * (Rj[o] - lam[o, j + 1])
      }
      Qj <- rep(NA_real_, n)
      Qj[o] <- h1[o, j + 1] * Qs[o, 2] + (1 - h1[o, j + 1]) * Qs[o, 1]
      Qnext <- Qj
    }
    psi <- mean(Qnext)
    list(psi = lo + sc * psi, ic = sc * (D + Qnext - psi), ic_lambda = sc * Dl,
         p_capped = if (is.null(h1_list)) mean((sweep(nuisance$lambda, 2, a, `*`) > 1)[x$obs]) else 0,
         min_n_eff = min(neff), max_weight = max(maxw))
  })
  res
}

## regression of a real-valued response, through the [0, 1] learners (response rescaled by its range)
.fit_real <- function(learner, learner_args, X, Z, label) {
  rg <- range(Z)
  if (diff(rg) < 1e-10) return(rep(mean(Z), length(Z)))
  rg[1] + diff(rg) * .fit_predict(learner, learner_args, X, (Z - rg[1]) / diff(rg), list(X), label)[[1]]
}

## outcome regression by instrument arm; returns n x 2 matrix (r = 0, r = 1), NA outside `o`
.q_arms <- function(learner, learner_args, H, Rj, Qnext, o, f, min_arm, label) {
  out <- matrix(NA_real_, length(Rj), 2)
  arms_ok <- all(vapply(0:1, function(r) sum(f & Rj == r) >= min_arm, TRUE))
  if (arms_ok) {
    for (r in 0:1) {
      s <- f & Rj == r
      out[o, r + 1] <- .fit_predict(learner, learner_args, H[s, , drop = FALSE], Qnext[s],
                                    list(H[o, , drop = FALSE]), paste0(label, "_R", r))[[1]]
    }
  } else {
    X <- cbind(H, R_now = Rj)
    nd <- lapply(0:1, function(r) cbind(H[o, , drop = FALSE], R_now = r))
    p <- .fit_predict(learner, learner_args, X[f, , drop = FALSE], Qnext[f], nd, paste0(label, "_joint"))
    out[o, 1] <- p[[1]]; out[o, 2] <- p[[2]]
  }
  pmin(pmax(out, 1e-5), 1 - 1e-5)
}

## intercept-only weighted logistic fluctuation: solves sum w (y - expit(off + eps)) = 0
.fluctuate <- function(y, off, w) {
  if (!length(y) || sum(w) <= 0) return(0)
  s <- function(e) sum(w * (y - stats::plogis(off + e)))
  lo <- s(-30); hi <- s(30)
  if (lo <= 0) return(-30)
  if (hi >= 0) return(30)
  stats::uniroot(s, c(-30, 30), tol = 1e-10)$root
}

#' @export
print.lii_fit <- function(x, digits = 4, ...) {
  cat("Implied interventions: mean of", x$settings$outcome,
      "under alpha-scaling of the instrument\n")
  cat("n =", x$settings$n, " visits 0..", x$settings$K, sep = "")
  cat("  (SEs from the influence curve: data-adaptive parameter)\n\n")
  e <- x$estimates
  print(format(data.frame(alpha = e$alpha, estimate = e$estimate, se = e$se,
                          ci95 = sprintf("[%.*f, %.*f]", digits, e$lower, digits, e$upper),
                          diff_vs_1 = e$diff,
                          diff_ci95 = sprintf("[%.*f, %.*f]", digits, e$diff_lower, digits, e$diff_upper),
                          p_capped = e$p_capped, min_n_eff = round(e$min_n_eff)),
               digits = digits), row.names = FALSE)
  if (!is.null(e$estimate_fix)) {
    cat("\nParameter defined with the true instrument intensity (full efficient influence curve):\n")
    print(format(data.frame(alpha = e$alpha, estimate = e$estimate_fix, se = e$se_fix,
                            ci95 = sprintf("[%.*f, %.*f]", digits, e$lower_fix, digits, e$upper_fix),
                            diff_vs_1 = e$diff_fix,
                            diff_ci95 = sprintf("[%.*f, %.*f]", digits, e$diff_lower_fix, digits, e$diff_upper_fix)),
                 digits = digits), row.names = FALSE)
  }
  if (any(e$p_capped > 0.05))
    cat("\nNote: for some alphas more than 5% of visits have alpha*lambda > 1 (capped at 1).\n")
  invisible(x)
}

#' @export
summary.lii_fit <- function(object, ...) object$estimates

#' Plot the dose-response curve in alpha
#' @param x An object from [lii_tmle()].
#' @param diff Plot the contrast with alpha = 1 instead of the mean.
#' @param ... Passed to [graphics::plot()].
#' @export
plot.lii_fit <- function(x, diff = FALSE, ...) {
  e <- x$estimates
  y <- if (diff) e$diff else e$estimate
  lo <- if (diff) e$diff_lower else e$lower
  hi <- if (diff) e$diff_upper else e$upper
  graphics::plot(e$alpha, y, type = "b", pch = 19, ylim = range(lo, hi),
                 xlab = expression(alpha ~ "(scaling of the instrument probability)"),
                 ylab = if (diff) "Difference with alpha = 1" else paste("Mean of", x$settings$outcome), ...)
  graphics::arrows(e$alpha, lo, e$alpha, hi, angle = 90, code = 3, length = 0.04)
  graphics::abline(v = 1, lty = 3)
  if (diff) graphics::abline(h = 0, lty = 2)
  invisible(x)
}
