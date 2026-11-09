## History summaries used as regressors.
##
## F-history at visit j (main.tex's F_j = (L-bar(j), A-bar(j-1)), WITHOUT instruments):
##   baseline W, visit index j, L(j), lags L(j-1..j-lags), A(j-1..j-lags), count of past A.
## Full history at visit j (everything observed before R(j)): F-history plus
##   R(j-1..j-lags) and the count of past R.
## Lags before visit 0 are padded with zeros (the visit index tells the learner).
## Rows of subjects not observed at visit j are NA.

.lag_cols <- function(M, j, lags, prefix) {
  if (lags < 1) return(NULL)
  out <- lapply(seq_len(lags), function(l) {
    v <- if (j - l >= 0) M[, j - l + 1] else rep(0, nrow(M))
    v[is.na(v)] <- 0
    v
  })
  out <- do.call(cbind, out)
  colnames(out) <- paste0(prefix, "_lag", seq_len(lags))
  out
}

.past_count <- function(M, j) {
  if (j == 0) return(rep(0, nrow(M)))
  rowSums(M[, seq_len(j), drop = FALSE], na.rm = TRUE)
}

.history <- function(x, j, lags = 1L, full = TRUE) {
  n <- x$n
  parts <- list(visit = matrix(j, n, 1, dimnames = list(NULL, "visit")))
  if (!is.null(x$W)) parts$W <- x$W
  Lj <- x$L[[j + 1]]; colnames(Lj) <- paste0(colnames(Lj), "_now")
  parts$L <- Lj
  if (lags >= 1) for (l in seq_len(lags)) {
    Ll <- if (j - l >= 0) x$L[[j - l + 1]] else matrix(0, n, ncol(Lj))
    Ll[is.na(Ll)] <- 0
    colnames(Ll) <- paste0(colnames(x$L[[1]]), "_lag", l)
    parts[[paste0("L", l)]] <- Ll
  }
  parts$A <- .lag_cols(x$A, j, lags, "A")
  parts$nA <- matrix(.past_count(x$A, j), n, 1, dimnames = list(NULL, "A_count"))
  if (full) {
    parts$R <- .lag_cols(x$R, j, lags, "R")
    parts$nR <- matrix(.past_count(x$R, j), n, 1, dimnames = list(NULL, "R_count"))
  }
  H <- do.call(cbind, parts[!vapply(parts, is.null, TRUE)])
  H[!x$obs[, j + 1], ] <- NA
  H
}

## pooled design over visits for the models of R(j), A(j) and censoring
.pooled <- function(x, lags, full, extra = NULL, rows = function(j) x$obs[, j + 1]) {
  out <- lapply(0:x$K, function(j) {
    r <- rows(j)
    H <- .history(x, j, lags, full)[r, , drop = FALSE]
    if (!is.null(extra)) H <- cbind(H, extra(j)[r, , drop = FALSE])
    list(X = H, i = which(r), j = rep(j, sum(r)))
  })
  list(X = do.call(rbind, lapply(out, `[[`, "X")),
       i = unlist(lapply(out, `[[`, "i")), j = unlist(lapply(out, `[[`, "j")))
}
