## Coverage of lii_bootstrap() percentile CIs (fixed parameter) with learner_glm.
## Rscript inst/validation/coverage_bootstrap_glm.R [reps] [B] [n]
args <- commandArgs(TRUE)
reps <- if (length(args) >= 1) as.integer(args[1]) else 100L
B <- if (length(args) >= 2) as.integer(args[2]) else 100L
n <- if (length(args) >= 3) as.integer(args[3]) else 1000L
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(quiet = TRUE)
} else library(impliedIV)
alphas <- c(0.5, 1, 1.5)
truth <- simulate_lii(n = 10, K = 2, seed = 1)$truth(alphas, n_mc = 1e6, seed = 99)
dtruth <- truth - truth[alphas == 1]
res <- do.call(rbind, lapply(seq_len(reps), function(r) {
  sim <- simulate_lii(n = n, K = 2, seed = 1000 + r)
  d <- lii_data(sim$data, "id", "visit", "R", "A", "Y", covariates = "L", baseline = "W")
  b <- lii_bootstrap(d, alphas, learner = learner_glm, B = B, seed = r, verbose = FALSE)$estimates
  data.frame(alpha = b$alpha, cover = b$lower <= truth & truth <= b$upper,
             dcover = b$diff_lower <= dtruth & dtruth <= b$diff_upper, boot_se = b$boot_se,
             diff_boot_se = b$diff_boot_se)
}))
summ <- aggregate(cbind(cover, dcover, boot_se, diff_boot_se) ~ alpha, res, mean)
names(summ) <- c("alpha", "coverage", "diff_coverage", "mean_boot_se", "mean_diff_boot_se")
cat("reps =", reps, " B =", B, " n =", n, "\n"); print(summ)
dir.create("inst/validation/results", showWarnings = FALSE)
write.csv(summ, sprintf("inst/validation/results/coverage_bootstrap_glm_n%d_reps%d_B%d.csv", n, reps, B),
          row.names = FALSE)
