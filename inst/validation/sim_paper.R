## Simulation study of the paper. One process = one (scenario, learner, n) cell and a block of
## replicates; results are appended to a CSV, one row per replicate and alpha.
## Rscript inst/validation/sim_paper.R <scenario: fresh|persistent> <learner: glm|hapc> <n> <first_rep> <last_rep> <outdir>
args <- commandArgs(TRUE)
scen <- args[1]; lrn <- args[2]; n <- as.integer(args[3])
reps <- as.integer(args[4]):as.integer(args[5]); outdir <- args[6]
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(quiet = TRUE)
} else library(impliedIV)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out <- file.path(outdir, sprintf("sim_%s_%s_n%d_r%d-%d.csv", scen, lrn, n, min(reps), max(reps)))
alphas <- if (lrn == "glm") c(0.5, 0.75, 1, 1.25, 1.5) else c(0.5, 1, 1.5)
learner <- if (lrn == "glm") learner_glm else learner_hapc
for (r in reps) {
  t0 <- Sys.time()
  sim <- simulate_lii(n = n, K = 2, seed = 5000 + r, censoring = FALSE, confounder = scen)
  d <- lii_data(sim$data, "id", "visit", "R", "A", "Y", covariates = "L", baseline = "W")
  row <- tryCatch({
    nu <- lii_nuisance(d, learner = learner, pooled = FALSE, verbose = FALSE)
    fit <- lii_tmle(d, alphas, learner = learner, nuisance = nu, fixed = TRUE, verbose = FALSE)
    e <- fit$estimates
    nv <- lii_naive(d, alphas, nu, learner = learner, verbose = FALSE)
    imp <- tapply(fit$implied$implied_p_treat, fit$implied$alpha, mean)
    trt <- if (lrn == "glm") lii_tmle(lii_treatment_outcome(d), alphas, learner = learner, nuisance = nu,
                                     fixed = FALSE, verbose = FALSE)$estimates$estimate else NA
    lagp <- lag_sufficiency_test(sim$data, "id", "visit", "R", "A", covariates = c("L", "W"), lags = 1)$p_value
    data.frame(scenario = scen, learner = lrn, n = n, rep = r, alpha = e$alpha,
               est = e$estimate, se = e$se, diff = e$diff, diff_se = e$diff_se,
               est_fix = e$estimate_fix, se_fix = e$se_fix, diff_fix = e$diff_fix, diff_se_fix = e$diff_se_fix,
               naive = nv$estimate, naive_se = nv$se, naive_diff = nv$diff, naive_diff_se = nv$diff_se,
               implied_treat = as.numeric(imp), tmle_treat = trt, lag_p = lagp,
               min_n_eff = e$min_n_eff, p_capped = e$p_capped,
               seconds = round(as.numeric(Sys.time() - t0, units = "secs"), 1))
  }, error = function(err) { message("rep ", r, " failed: ", conditionMessage(err)); NULL })
  if (!is.null(row)) write.table(row, out, sep = ",", row.names = FALSE, col.names = !file.exists(out),
                                 append = file.exists(out))
}
