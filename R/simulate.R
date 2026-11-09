#' Simulate data with a hidden confounder and a binary instrument
#'
#' Prescription -> adherence toy model with visits `0..K`:
#' * `W` baseline covariate, `U` unmeasured (affects adherence and outcome);
#' * `L(j)`: severity at visit `j`, depends on `L(j-1)`, `W`, `U` and
#'   previous adherence;
#' * `R(j)`: prescription, \eqn{P(R(j)=1) = \mathrm{expit}(r_0 + 0.6L(j) -
#'   0.5A(j-1))}, so its natural law depends on the past only through
#'   \eqn{F_j} and \eqn{\lambda_j} is known;
#' * `A(j)`: adherence, depends on `R(j)` (strongly), `L(j)`, `A(j-1)` and `U`
#'   (hidden confounding), not on past prescriptions;
#' * censoring after each visit with probability depending on `L(j)` (if
#'   `censoring = TRUE`);
#' * `Y`: binary outcome after visit `K`, depends on the number of visits
#'   with adherence, `L(K)`, `W` and `U`, never on `R` (exclusion).
#'
#' `truth(alphas)` returns Monte Carlo values of \eqn{\psi(\alpha)} (no
#' censoring), for checking estimators.
#' @param n Number of subjects.
#' @param K Last visit index.
#' @param r0 Intercept of the prescription model (controls how often it fires).
#' @param effect Log-odds effect on `Y` of each visit with adherence.
#' @param censoring Simulate censoring.
#' @param seed Seed, or `NULL`.
#' @return A list with `data` (long `data.frame`: `id`, `visit`, `W`, `L`, `R`,
#'   `A`, `Y` on the last row of each subject) and the function `truth`.
#' @export
simulate_lii <- function(n = 1000, K = 2, r0 = -0.5, effect = -0.6, censoring = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  gen <- function(n, alpha = NULL, cens = censoring) {
    W <- stats::rnorm(n); U <- stats::rnorm(n)
    L <- R <- A <- matrix(NA_real_, n, K + 1)
    alive <- rep(TRUE, n); Aprev <- rep(0, n); Lprev <- rep(0, n); nA <- rep(0, n)
    last <- rep(K, n)
    for (j in 0:K) {
      L[, j + 1] <- 0.5 * Lprev + 0.4 * W + 0.4 * U - 0.5 * Aprev + stats::rnorm(n, sd = 0.7)
      lam <- stats::plogis(r0 + 0.6 * L[, j + 1] - 0.5 * Aprev)
      p <- if (is.null(alpha)) lam else pmin(alpha * lam, 1)
      R[, j + 1] <- stats::rbinom(n, 1, p)
      A[, j + 1] <- stats::rbinom(n, 1, stats::plogis(-1.5 + 2.2 * R[, j + 1] + 0.4 * L[, j + 1] +
                                                        0.8 * Aprev + 0.8 * U))
      Aprev <- A[, j + 1]; Lprev <- L[, j + 1]; nA <- nA + A[, j + 1]
      if (cens) {
        stay <- stats::rbinom(n, 1, stats::plogis(2.5 - 0.4 * L[, j + 1])) == 1
        newly <- alive & !stay
        last[newly] <- j
        alive <- alive & stay
      }
    }
    Y <- stats::rbinom(n, 1, stats::plogis(-0.3 + effect * nA + 0.5 * L[, K + 1] + 0.4 * W + 0.9 * U))
    list(W = W, L = L, R = R, A = A, Y = Y, alive = alive, last = last)
  }
  s <- gen(n)
  rows <- lapply(seq_len(n), function(i) {
    v <- 0:s$last[i]
    data.frame(id = i, visit = v, W = s$W[i], L = s$L[i, v + 1], R = s$R[i, v + 1],
               A = s$A[i, v + 1], Y = c(rep(NA, length(v) - 1), if (s$alive[i]) s$Y[i] else NA))
  })
  data <- do.call(rbind, rows)
  truth <- function(alphas, n_mc = 2e5, seed = 1) {
    if (!is.null(seed)) set.seed(seed)
    vapply(alphas, function(a) mean(gen(n_mc, alpha = a, cens = FALSE)$Y), 0)
  }
  list(data = data, truth = truth)
}
