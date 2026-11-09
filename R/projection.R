#' Project a desired treatment intervention onto the instrument interventions
#'
#' The question an analyst usually has is about the **treatment**: "what if
#' patients took the drug 20% more often?". Without an unconfounded
#' treatment that question is not identified. What is identified is the effect
#' of changing the **instrument**, and each change of the instrument implies a
#' change of the treatment mechanism,
#' \deqn{g_j(\alpha_j)(F_j) = \pi_0(F_j) + \alpha_j\lambda_j(F_j)[\pi_1(F_j) - \pi_0(F_j)].}
#' This function finds the instrument intervention whose implied treatment
#' mechanism is closest to the desired one, \eqn{g^*_j}, in Kullback–Leibler
#' divergence averaged over the observed histories:
#' \deqn{\alpha^\dagger_j = \arg\min_{\alpha_j \ge 0} E_{F_j}\, \mathrm{KL}\bigl(g^*_j(\cdot \mid F_j),\, g_j(\alpha_j)(\cdot \mid F_j)\bigr),}
#' one problem per visit (`per_visit = TRUE`), or the sum over visits with a
#' common \eqn{\alpha} (`per_visit = FALSE`). Then it estimates the causal
#' effect of that projected intervention, \eqn{E[Y^{\alpha^\dagger}]}, and its
#' contrast with the current world (\eqn{\alpha = 1}), with [lii_tmle()]'s
#' estimator.
#'
#' **Reachability.** Scaling the instrument can move the treatment probability
#' only between \eqn{\pi_0} (instrument never fires) and \eqn{\pi_1} (always
#' fires). A desired intervention outside that range is not reachable, and the
#' projection stops at the boundary. Always compare `target_p_treat` with
#' `implied_p_treat` and with the reachable range in `projection` before
#' interpreting the effect as "the effect of the desired intervention": it is
#' the effect of the closest identifiable one.
#'
#' **Inference.** The projection uses the nuisance fits but not the outcome,
#' so the influence-curve CIs are valid for the data-adaptive parameter
#' \eqn{E[Y^{\hat\alpha^\dagger}]}, treating \eqn{\hat\alpha^\dagger} and
#' \eqn{\hat\lambda} as given.
#'
#' @param x An object from [lii_data()].
#' @param target The desired intervention on the treatment, \eqn{g^*_j(1 \mid F_j)}:
#'   * a number `c`: multiply the current probability of treatment,
#'     \eqn{P(A(j) = 1 \mid F_j)}, by `c` at every visit, capped at 1
#'     (e.g. `1.2` = "20% more treatment, proportionally to current practice");
#'   * a vector of length `K + 1`: one such factor per visit;
#'   * a function `function(visit, history)` returning the desired probability
#'     for each row of `history` (a `data.frame` with the history at that
#'     visit: `visit`, baseline covariates, `<covariate>_now`, lags, `A_lag1`,
#'     `A_count`, ...). Example: `function(j, h) as.numeric(h$severity_now > 2)`.
#' @param per_visit One alpha per visit (`TRUE`) or a common alpha (`FALSE`).
#'   Default: `TRUE` when `target` has one factor per visit, `FALSE` otherwise.
#' @param alpha_max Upper end of the search interval for alpha.
#' @param nuisance Optional output of [lii_nuisance()] to reuse.
#' @param learner,learner_args,lags,instrument_markov,gbound,min_arm,verbose
#'   As in [lii_tmle()].
#' @return An object of class `lii_projection`:
#'   * `alpha`: projected scaling per visit;
#'   * `projection`: by visit, mean probability of treatment now
#'     (`natural_p_treat`), desired (`target_p_treat`), under the projection
#'     (`implied_p_treat`), the reachable range (`reachable_min` = mean
#'     \eqn{\pi_0}, `reachable_max` = mean \eqn{\pi_1}) and the mean KL
#'     divergence to the target now and after projection;
#'   * `effect`: estimate, SE and 95% CI of the mean outcome under the projected
#'     intervention and under the current world, and of their difference;
#'   * `ic`, `nuisance`, `settings`.
#' @examples
#' sim <- simulate_lii(n = 500, K = 2, seed = 3)
#' d <- lii_data(sim$data, "id", "visit", "R", "A", "Y", covariates = "L", baseline = "W")
#' pr <- lii_project(d, target = 1.2, learner = learner_glm, verbose = FALSE)
#' pr
#' @export
lii_project <- function(x, target, per_visit = NULL, alpha_max = 20, nuisance = NULL,
                        learner = learner_hapc, learner_args = list(), lags = 1L,
                        instrument_markov = FALSE, gbound = 0.01, min_arm = 30L, verbose = TRUE) {
  stopifnot(inherits(x, "lii_data"))
  K <- x$K
  if (is.null(nuisance)) nuisance <- lii_nuisance(x, learner, learner_args, lags,
                                                   instrument_markov, verbose)
  if (is.null(per_visit)) per_visit <- is.numeric(target) && length(target) == K + 1 && K > 0
  ## desired probabilities g*_j on the observed rows
  tg <- lapply(0:K, function(j) {
    o <- x$obs[, j + 1]
    lam <- nuisance$lambda[o, j + 1]; p0 <- nuisance$pi0[o, j + 1]; p1 <- nuisance$pi1[o, j + 1]
    nat <- p0 + lam * (p1 - p0)
    g <- if (is.function(target)) {
      h <- as.data.frame(.history(x, j, nuisance$lags, full = FALSE)[o, , drop = FALSE])
      as.numeric(target(j, h))
    } else if (is.numeric(target) && length(target) %in% c(1, K + 1)) {
      pmin(1, nat * if (length(target) == 1) target else target[j + 1])
    } else stop("`target` must be a number, a vector of length K + 1 or a function(visit, history).")
    if (length(g) != sum(o) || anyNA(g) || any(g < 0 | g > 1))
      stop("The target at visit ", j, " must give one probability in [0, 1] per subject observed.")
    list(lam = lam, p0 = p0, p1 = p1, nat = nat, g = g)
  })
  kl <- function(gs, g) {
    gs <- pmin(pmax(gs, 1e-6), 1 - 1e-6); g <- pmin(pmax(g, 1e-6), 1 - 1e-6)
    mean(gs * log(gs / g) + (1 - gs) * log((1 - gs) / (1 - g)))
  }
  implied <- function(t, a) t$p0 + pmin(a * t$lam, 1) * (t$p1 - t$p0)
  upper <- function(t) min(alpha_max, 1 / max(min(t$lam), 1e-6))
  if (per_visit) {
    alpha <- vapply(tg, function(t) stats::optimize(function(a) kl(t$g, implied(t, a)),
                                                    c(0, upper(t)))$minimum, 0)
  } else {
    up <- min(vapply(tg, upper, 0))
    a1 <- stats::optimize(function(a) sum(vapply(tg, function(t) kl(t$g, implied(t, a)), 0)),
                          c(0, up))$minimum
    alpha <- rep(a1, K + 1)
  }
  projection <- do.call(rbind, lapply(0:K, function(j) {
    t <- tg[[j + 1]]
    data.frame(visit = j, alpha = alpha[j + 1], n = length(t$g),
               natural_p_treat = mean(t$nat), target_p_treat = mean(t$g),
               implied_p_treat = mean(implied(t, alpha[j + 1])),
               reachable_min = mean(t$p0), reachable_max = mean(t$p1),
               kl_natural = kl(t$g, t$nat), kl_projected = kl(t$g, implied(t, alpha[j + 1])))
  }))
  res <- .tmle_core(x, nuisance, list(alpha, rep(1, K + 1)), c("projected", "natural"),
                    learner, learner_args, gbound, min_arm, verbose)
  ic <- cbind(projected = res[[1]]$ic, natural = res[[2]]$ic)
  est <- c(res[[1]]$psi, res[[2]]$psi)
  se <- apply(ic, 2, stats::sd) / sqrt(x$n)
  dse <- stats::sd(ic[, 1] - ic[, 2]) / sqrt(x$n)
  z <- stats::qnorm(0.975)
  effect <- data.frame(
    intervention = c("projected", "natural (alpha = 1)", "difference"),
    estimate = c(est, est[1] - est[2]), se = c(se, dse),
    lower = c(est, est[1] - est[2]) - z * c(se, dse),
    upper = c(est, est[1] - est[2]) + z * c(se, dse),
    min_n_eff = c(res[[1]]$min_n_eff, res[[2]]$min_n_eff, NA),
    max_weight = c(res[[1]]$max_weight, res[[2]]$max_weight, NA))
  structure(list(alpha = alpha, projection = projection, effect = effect, ic = ic,
                 nuisance = nuisance,
                 settings = list(n = x$n, K = K, per_visit = per_visit,
                                 target = if (is.function(target)) "function" else target,
                                 outcome = x$names$outcome)),
            class = "lii_projection")
}

#' @export
print.lii_projection <- function(x, digits = 3, ...) {
  cat("Projection of the desired treatment intervention onto instrument scalings\n")
  cat("alpha", if (x$settings$per_visit) "per visit:" else "(common):",
      paste(round(x$alpha, digits), collapse = ", "), "\n\n")
  print(format(x$projection, digits = digits), row.names = FALSE)
  far <- abs(x$projection$implied_p_treat - x$projection$target_p_treat) > 0.02
  if (any(far)) cat("\nNote: at visit(s)", paste(x$projection$visit[far], collapse = ", "),
                    "the projection does not reach the target (it is outside what the instrument can move).\n")
  cat("\nMean of", x$settings$outcome, "(SEs from the influence curve: data-adaptive parameter)\n")
  print(format(x$effect[, 1:6], digits = digits), row.names = FALSE)
  invisible(x)
}
