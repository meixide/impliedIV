## One end-to-end run with learner_hapc, timing and comparison with the truth.
## Rscript inst/validation/run_hapc.R [n]
args <- commandArgs(TRUE); n <- if (length(args)) as.integer(args[1]) else 1000L
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(quiet = TRUE)
} else library(impliedIV)
sim <- simulate_lii(n = n, K = 2, seed = 11)
d <- lii_data(sim$data, "id", "visit", "R", "A", "Y", covariates = "L", baseline = "W")
t0 <- Sys.time()
fit <- lii_tmle(d, alphas = c(0.5, 1.5), learner = learner_hapc)
cat("seconds:", round(as.numeric(Sys.time() - t0, units = "secs")), "\n")
print(fit)
truth <- sim$truth(fit$estimates$alpha, n_mc = 1e6, seed = 99)
print(cbind(fit$estimates[, c("alpha", "estimate", "lower", "upper")], truth = truth))
dg <- lii_diagnostics(fit)
print(dg$first_stage); print(dg$fits)
