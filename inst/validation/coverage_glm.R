## Coverage of the influence-curve CIs with learner_glm on simulate_lii().
## Run from the package root: Rscript inst/validation/coverage_glm.R [reps] [n]
args <- commandArgs(TRUE)
reps <- if (length(args) >= 1) as.integer(args[1]) else 200L
n <- if (length(args) >= 2) as.integer(args[2]) else 1000L
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(quiet = TRUE)
} else library(impliedIV)
alphas <- c(0.5, 0.75, 1, 1.25, 1.5)
truth <- simulate_lii(n = 10, K = 2, seed = 1)$truth(alphas, n_mc = 1e6, seed = 99)
dtruth <- truth - truth[alphas == 1]
res <- lapply(seq_len(reps), function(r) {
  sim <- simulate_lii(n = n, K = 2, seed = 1000 + r)
  d <- lii_data(sim$data, "id", "visit", "R", "A", "Y", covariates = "L", baseline = "W")
  e <- lii_tmle(d, alphas, learner = learner_glm, verbose = FALSE)$estimates
  data.frame(rep = r, alpha = e$alpha, est = e$estimate, se = e$se,
             cover = e$lower <= truth & truth <= e$upper,
             dest = e$diff, dse = e$diff_se,
             dcover = e$diff_lower <= dtruth & dtruth <= e$diff_upper)
})
res <- do.call(rbind, res)
summ <- do.call(rbind, lapply(split(res, res$alpha), function(s) {
  a <- s$alpha[1]
  data.frame(alpha = a, truth = truth[alphas == a], bias = mean(s$est) - truth[alphas == a],
             sd = sd(s$est), mean_se = mean(s$se), coverage = mean(s$cover),
             diff_truth = dtruth[alphas == a], diff_bias = mean(s$dest) - dtruth[alphas == a],
             diff_sd = sd(s$dest), diff_mean_se = mean(s$dse), diff_coverage = mean(s$dcover))
}))
cat("reps =", reps, " n =", n, " K = 2, learner_glm\n")
print(format(summ, digits = 3), row.names = FALSE)
dir.create("inst/validation/results", showWarnings = FALSE)
write.csv(summ, sprintf("inst/validation/results/coverage_glm_n%d_reps%d.csv", n, reps), row.names = FALSE)
