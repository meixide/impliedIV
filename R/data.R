#' Prepare longitudinal data for the implied-intervention estimators
#'
#' Converts a long data set (one row per subject and visit) into the internal
#' structure used by [lii_tmle()]. Visits are indexed `j = 0, ..., K` by their
#' order within each subject (the `time` column is only used for ordering).
#' Within each visit the temporal order of the nodes is
#' `L(j) -> R(j) -> A(j)`: covariates measured at the visit, then the
#' instrument (for example "was the drug prescribed at this visit?"), then the
#' treatment it encourages (for example "did the patient take it?"). The
#' outcome `Y` is measured after the last visit `K`.
#'
#' **Censoring.** A subject with fewer than `K + 1` visits is censored after
#' their last visit; a subject with `K + 1` visits but a missing outcome is
#' censored after visit `K`. Censoring is assumed monotone and explainable by
#' the observed past (it is intervened on: estimands are "had nobody been
#' censored"). Subjects with more than `K + 1` visits are truncated at visit
#' `K`.
#'
#' **Outcome.** Taken from the last row of each subject (the last non-missing
#' value among that subject's rows). Binary or continuous; continuous outcomes
#' are rescaled to \[0, 1\] internally and estimates are reported on the
#' original scale.
#'
#' @param data A `data.frame` in long format.
#' @param id Name of the subject identifier column.
#' @param time Name of the column that orders the visits within a subject.
#' @param instrument Name of the binary (0/1) instrument column `R(j)`.
#' @param treatment Name of the binary (0/1) treatment column `A(j)`.
#' @param outcome Name of the outcome column (see Details).
#' @param covariates Names of the time-varying covariate columns `L(j)`,
#'   measured at the visit before the instrument. Numeric or factor; factors
#'   are expanded into indicator columns.
#' @param baseline Names of baseline covariates (constant within subject;
#'   the value at the first visit is used), or `NULL`.
#' @param K Index of the last visit used (visits `0..K`). Default: the
#'   largest number of visits observed minus one. Choose `K` so that enough
#'   subjects reach visit `K` (see [lii_describe()]).
#' @return An object of class `lii_data`.
#' @examples
#' sim <- simulate_lii(n = 200, K = 2, seed = 1)
#' d <- lii_data(sim$data, id = "id", time = "visit", instrument = "R",
#'               treatment = "A", outcome = "Y", covariates = "L",
#'               baseline = "W")
#' d
#' @export
lii_data <- function(data, id, time, instrument, treatment, outcome,
                     covariates, baseline = NULL, K = NULL) {
  data <- as.data.frame(data)
  need <- c(id, time, instrument, treatment, outcome, covariates, baseline)
  miss <- setdiff(need, names(data))
  if (length(miss)) stop("Columns not found in `data`: ", paste(miss, collapse = ", "))
  if (anyNA(data[[id]]) || anyNA(data[[time]])) stop("`id` and `time` must not contain NA.")
  data <- data[order(data[[id]], data[[time]]), , drop = FALSE]
  if (anyDuplicated(data[, c(id, time)])) stop("Duplicated (id, time) rows.")

  ids <- unique(data[[id]])
  n <- length(ids)
  idx <- match(data[[id]], ids)
  visit <- stats::ave(seq_len(nrow(data)), idx, FUN = seq_along) - 1L
  nvis <- tabulate(idx, nbins = n)
  if (is.null(K)) K <- max(nvis) - 1L
  K <- as.integer(K)
  if (K < 0L) stop("K must be >= 0.")
  keep <- visit <= K
  data <- data[keep, , drop = FALSE]; idx <- idx[keep]; visit <- visit[keep]
  nvis <- pmin(nvis, K + 1L)

  ## covariate design (factors -> indicators), shared columns across visits
  Lmat <- .numeric_design(data[, covariates, drop = FALSE], "covariates")
  Wmat <- NULL
  if (length(baseline)) {
    first <- !duplicated(idx)
    Wmat <- .numeric_design(data[first, baseline, drop = FALSE], "baseline")
    Wmat <- Wmat[order(idx[first]), , drop = FALSE]
  }

  ## instrument and treatment must be 0/1 where the visit is observed
  Rv <- data[[instrument]]; Av <- data[[treatment]]
  for (nm in c("instrument", "treatment")) {
    v <- if (nm == "instrument") Rv else Av
    if (anyNA(v)) stop("The ", nm, " column has NA at observed visits. ",
                       "Drop or impute those visits before calling lii_data().")
    if (!all(v %in% c(0, 1))) stop("The ", nm, " column must be 0/1 ",
                                   "(a continuous instrument can be binarised with binarize_instrument()).")
  }

  obs <- matrix(FALSE, n, K + 1L)
  obs[cbind(idx, visit + 1L)] <- TRUE
  Rm <- Am <- matrix(NA_real_, n, K + 1L)
  Rm[cbind(idx, visit + 1L)] <- Rv
  Am[cbind(idx, visit + 1L)] <- Av
  L <- lapply(0:K, function(j) {
    M <- matrix(NA_real_, n, ncol(Lmat), dimnames = list(NULL, colnames(Lmat)))
    r <- visit == j
    M[idx[r], ] <- Lmat[r, , drop = FALSE]
    M
  })

  ## outcome: last non-missing value per subject, only if observed through K
  yv <- data[[outcome]]
  Y <- rep(NA_real_, n)
  for (i in which(nvis == K + 1L)) {
    yi <- yv[idx == i]; yi <- yi[!is.na(yi)]
    if (length(yi)) Y[i] <- yi[length(yi)]
  }
  if (!is.numeric(Y)) stop("The outcome must be numeric (0/1 or continuous).")
  y_obs <- !is.na(Y)
  if (!any(y_obs)) stop("No subject has an observed outcome after visit K = ", K, ".")
  yr <- range(Y[y_obs])
  binary <- all(Y[y_obs] %in% c(0, 1))
  y_range <- if (binary || (yr[1] >= 0 && yr[2] <= 1)) c(0, 1) else yr
  if (diff(y_range) == 0) stop("The outcome is constant.")

  structure(list(
    n = n, K = K, ids = ids, W = Wmat, L = L, R = Rm, A = Am, obs = obs,
    Y = Y, y_obs = y_obs, y_range = y_range, y_binary = binary,
    names = list(id = id, time = time, instrument = instrument, treatment = treatment,
                 outcome = outcome, covariates = covariates, baseline = baseline)
  ), class = "lii_data")
}

.numeric_design <- function(df, what) {
  for (v in names(df)) if (is.character(df[[v]]) || is.logical(df[[v]])) df[[v]] <- factor(df[[v]])
  if (anyNA(df)) stop("The ", what, " contain NA at observed visits. ",
                      "Impute them (or add missingness indicators) before calling lii_data().")
  M <- stats::model.matrix(~ ., data = df)[, -1, drop = FALSE]
  colnames(M) <- make.names(colnames(M), unique = TRUE)
  M
}

#' @export
print.lii_data <- function(x, ...) {
  cat("lii_data:", x$n, "subjects, visits 0..", x$K, "\n", sep = "")
  cat("  instrument:", x$names$instrument, " treatment:", x$names$treatment,
      " outcome:", x$names$outcome, if (x$y_binary) "(binary)" else "(continuous)", "\n")
  cat("  covariates per visit:", ncol(x$L[[1]]),
      " baseline covariates:", if (is.null(x$W)) 0 else ncol(x$W), "\n")
  cat("  outcome observed for", sum(x$y_obs), "subjects (", round(100 * mean(x$y_obs), 1), "% )\n")
  invisible(x)
}

#' Describe a prepared data set visit by visit
#'
#' Number of subjects still observed, proportion with the instrument and with
#' the treatment, and the crude difference in treatment proportion between
#' visits with and without the instrument (a first look at relevance).
#' @param x An object from [lii_data()].
#' @return A `data.frame` with one row per visit.
#' @export
lii_describe <- function(x) {
  stopifnot(inherits(x, "lii_data"))
  do.call(rbind, lapply(0:x$K, function(j) {
    o <- x$obs[, j + 1]
    r <- x$R[o, j + 1]; a <- x$A[o, j + 1]
    data.frame(visit = j, n_observed = sum(o), p_instrument = mean(r), p_treatment = mean(a),
               p_treat_if_R1 = if (any(r == 1)) mean(a[r == 1]) else NA_real_,
               p_treat_if_R0 = if (any(r == 0)) mean(a[r == 0]) else NA_real_)
  }))
}

## subset subjects (used by the bootstrap)
.subset_lii <- function(x, i) {
  x$n <- length(i); x$ids <- x$ids[i]
  if (!is.null(x$W)) x$W <- x$W[i, , drop = FALSE]
  x$L <- lapply(x$L, function(M) M[i, , drop = FALSE])
  x$R <- x$R[i, , drop = FALSE]; x$A <- x$A[i, , drop = FALSE]
  x$obs <- x$obs[i, , drop = FALSE]; x$Y <- x$Y[i]; x$y_obs <- x$y_obs[i]
  x
}
