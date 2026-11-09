.lii_log <- new.env()
.lii_log$fits <- list()

#' Learners for the nuisance regressions
#'
#' A learner is a function `learner(X, Y, newdata, label)` that fits
#' \eqn{P(Y = 1 \mid X)} for `Y` in \[0, 1\] (soft labels allowed: the
#' sequential regressions have fractional pseudo-outcomes) and returns the
#' predicted probabilities for the rows of `newdata` (same columns as `X`).
#'
#' * `learner_hapc()` (default): Highly Adaptive Principal Components with
#'   an L1 penalty and cross-validated lambda, via [hapc::cv.hapc()]
#'   (`family = "binomial"`, `norm = "1"`). `norm = "1"` (logistic lasso on the
#'   HAPC basis, through glmnet) is fixed on purpose: the sequential
#'   regressions have soft labels in \[0, 1\], which hapc supports only with
#'   `norm = "1"` or `"2"`. Defaults (`max_degree = 2`, `floor(n^0.7)`
#'   components, 5 folds, log-lambda grid from -8 to -1 with 15 points,
#'   extended up to 3 times when the optimum is on an edge) are those of
#'   the `tmle2t` package. The HAL kernel is n x n, so fits
#'   with more than `max_n` rows use a random subsample of `max_n` rows
#'   (recorded in [lii_fit_log()]). If a fit fails, the main-terms logistic
#'   regression is used instead and the failure is recorded.
#' * `learner_glm()`: main-terms logistic regression (quasibinomial). Fast;
#'   use it for a first pass, for simulations and as a sensitivity analysis.
#'
#' Any function with the same signature can be passed to [lii_tmle()].
#' @param X Numeric matrix of training regressors.
#' @param Y Numeric vector in \[0, 1\].
#' @param newdata Numeric matrix with the columns of `X`.
#' @param label Character tag recorded in [lii_fit_log()].
#' @param max_degree,npc_pow,nfolds,log_lambda_min,log_lambda_max,grid_length,extend_grid
#'   hapc settings: HAL kernel degree, number of principal components
#'   `floor(n^npc_pow)`, CV folds and lambda grid (see [hapc::cv.hapc()]).
#' @param max_n Largest number of rows used in one hapc fit.
#' @return Numeric vector of predictions, one per row of `newdata`.
#' @name learners
NULL

#' @rdname learners
#' @export
learner_glm <- function(X, Y, newdata, label = "") {
  X <- .clean(X); keep <- .nonconstant(X)
  if (stats::sd(Y) < 1e-12 || !any(keep)) return(rep(mean(Y), nrow(newdata)))
  X <- X[, keep, drop = FALSE]; newdata <- .clean(newdata)[, keep, drop = FALSE]
  f <- suppressWarnings(stats::glm.fit(cbind(1, X), Y, family = stats::quasibinomial()))
  b <- f$coefficients; b[is.na(b)] <- 0
  .log_fit(label, "glm", nrow(X), ncol(X), NA_real_, NA_character_, mean(Y), 0)
  as.numeric(stats::plogis(cbind(1, newdata) %*% b))
}

#' @rdname learners
#' @export
learner_hapc <- function(X, Y, newdata, label = "", max_degree = 2L, npc_pow = 0.7,
                         nfolds = 5L, log_lambda_min = -8, log_lambda_max = -1,
                         grid_length = 15L, extend_grid = 3L, max_n = 4000L) {
  if (!requireNamespace("hapc", quietly = TRUE))
    stop("learner_hapc() needs the 'hapc' package: remotes::install_github('meixide/hapc'). ",
         "Use learner = learner_glm for a quick run without it.")
  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("learner_hapc() needs 'glmnet' (hapc's norm = \"1\" path): install.packages('glmnet').")
  X <- .clean(X); newdata <- .clean(newdata); keep <- .nonconstant(X)
  if (stats::sd(Y) < 1e-12 || !any(keep)) return(rep(mean(Y), nrow(newdata)))
  X <- X[, keep, drop = FALSE]; newdata <- newdata[, keep, drop = FALSE]
  sub <- nrow(X) > max_n
  if (sub) { s <- sort(sample.int(nrow(X), max_n)); X <- X[s, , drop = FALSE]; Y <- Y[s] }
  if (stats::sd(Y) < 1e-12) return(rep(mean(Y), nrow(newdata)))
  t0 <- Sys.time()
  cv <- tryCatch(withCallingHandlers(
    hapc::cv.hapc(X, Y, family = "binomial", norm = "1", max_degree = max_degree,
                  npcs = max(1L, floor(nrow(X)^npc_pow)), nfolds = nfolds,
                  log_lambda_min = log_lambda_min, log_lambda_max = log_lambda_max,
                  grid_length = grid_length, extend_grid = extend_grid, predict = newdata),
    warning = function(w) if (grepl("Non-binary labels|edge", conditionMessage(w)))
      invokeRestart("muffleWarning")),
    error = function(e) NULL)
  secs <- round(as.numeric(Sys.time() - t0, units = "secs"), 1)
  if (is.null(cv)) {
    .log_fit(label, "hapc_failed_glm_used", nrow(X), ncol(X), NA_real_, "error", mean(Y), secs, sub)
    return(learner_glm(X, Y, newdata, label = paste0(label, " (fallback)")))
  }
  .log_fit(label, "hapc", nrow(X), ncol(X), cv$best_lambda,
           if (is.null(cv$lambda_position)) NA_character_ else cv$lambda_position, mean(Y), secs, sub)
  as.numeric(cv$predictions)
}

.clean <- function(M) { M <- as.matrix(M); storage.mode(M) <- "double"; dimnames(M) <- list(NULL, colnames(M)); M }
.nonconstant <- function(X) apply(X, 2, function(v) length(unique(v)) > 1)

.log_fit <- function(label, learner, n, p, lambda, position, mean_y, secs, subsampled = FALSE) {
  .lii_log$fits[[length(.lii_log$fits) + 1]] <- data.frame(
    model = label, learner = learner, n = n, p = p, lambda = lambda,
    lambda_position = position, mean_y = round(mean_y, 4), seconds = secs,
    subsampled = subsampled)
}

#' Log of the nuisance fits
#'
#' One row per learner call since the last [lii_reset_log()]: model label,
#' learner used, rows and columns, selected lambda and its position on the
#' grid (`"lower_edge"` / `"upper_edge"` suggest widening the grid;
#' `"null_model"` means the learner kept only the intercept), mean of the
#' response, time and whether the rows were subsampled.
#' @return A `data.frame`.
#' @export
lii_fit_log <- function() {
  if (!length(.lii_log$fits)) return(data.frame())
  do.call(rbind, .lii_log$fits)
}

#' @rdname lii_fit_log
#' @export
lii_reset_log <- function() { .lii_log$fits <- list(); invisible(NULL) }

## fit once, predict several newdata matrices
.fit_predict <- function(learner, learner_args, X, Y, newdata_list, label) {
  nd <- do.call(rbind, newdata_list)
  p <- do.call(learner, c(list(X = X, Y = Y, newdata = nd, label = label), learner_args))
  split(p, rep(seq_along(newdata_list), vapply(newdata_list, nrow, 1L)))
}
