## Summaries of the simulation study: LaTeX tables and a figure for the manuscript.
## Rscript inst/validation/make_tables.R <results dir> <output dir for tables> <output dir for figures>
args <- commandArgs(TRUE)
res_dir <- args[1]; tab_dir <- args[2]; fig_dir <- args[3]
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(quiet = TRUE)
} else library(impliedIV)
files <- list.files(res_dir, pattern = "^sim_.*\\.csv$", full.names = TRUE)
files <- files[file.size(files) > 0]
D <- do.call(rbind, lapply(files, read.csv))
D <- D[!duplicated(D[, c("scenario", "learner", "n", "rep", "alpha")]), ]
alphas <- c(0.5, 0.75, 1, 1.25, 1.5)
truth_file <- file.path(res_dir, "truth.rds")
if (file.exists(truth_file)) TR <- readRDS(truth_file) else {
  TR <- lapply(c(persistent = "persistent", fresh = "fresh"), function(sc) {
    s <- simulate_lii(n = 10, K = 2, censoring = FALSE, confounder = sc, seed = 1)
    list(psi = s$truth(alphas, n_mc = 4e6, seed = 2026),
         treat = rowMeans(s$truth_treatment(alphas, n_mc = 4e6, seed = 2027)))
  })
  saveRDS(TR, truth_file)
}
z <- qnorm(0.975)
summ <- do.call(rbind, lapply(split(D, list(D$scenario, D$learner, D$n, D$alpha), drop = TRUE), function(s) {
  sc <- s$scenario[1]; a <- s$alpha[1]; k <- match(a, alphas)
  psi <- TR[[sc]]$psi[k]; dpsi <- psi - TR[[sc]]$psi[3]; tau <- TR[[sc]]$treat[k]
  cov <- function(est, se, truth) mean(abs(est - truth) <= z * se)
  data.frame(scenario = sc, learner = s$learner[1], n = s$n[1], alpha = a, reps = nrow(s),
             truth = psi, bias = mean(s$est_fix) - psi, sd = sd(s$est_fix), se = mean(s$se_fix),
             cover = cov(s$est_fix, s$se_fix, psi), cover_da = cov(s$est, s$se, psi),
             dtruth = dpsi, dbias = mean(s$diff_fix) - dpsi, dsd = sd(s$diff_fix), dse = mean(s$diff_se_fix),
             dcover = cov(s$diff_fix, s$diff_se_fix, dpsi), dcover_da = cov(s$diff, s$diff_se, dpsi),
             nbias = mean(s$naive_diff) - dpsi, nsd = sd(s$naive_diff), ncover = cov(s$naive_diff, s$naive_diff_se, dpsi),
             tau = tau, implied = mean(s$implied_treat), tmle_treat = mean(s$tmle_treat),
             lag_reject = mean(s$lag_p < 0.05), min_n_eff = mean(s$min_n_eff), seconds = mean(s$seconds))
}))
summ <- summ[order(summ$scenario, summ$learner, summ$n, summ$alpha), ]
write.csv(summ, file.path(res_dir, "summary.csv"), row.names = FALSE)
f3 <- function(x) formatC(x, format = "f", digits = 3); f2 <- function(x) formatC(x, format = "f", digits = 2)
sgn <- function(x) ifelse(round(x, 3) < 0, paste0("$-$", f3(abs(x))), f3(abs(x)))
lab <- c(persistent = "Persistent", fresh = "Renewed")

## Table 1: contrasts, logistic regression
t1 <- summ[summ$learner == "glm" & summ$alpha %in% c(0.5, 1.5) & summ$n %in% c(500, 1000, 2000), ]
rows <- with(t1, paste(lab[scenario], n, alpha, sgn(dtruth), sgn(dbias), f3(dsd), f3(dse), f2(dcover), f2(dcover_da),
                       sgn(nbias), f2(ncover), sep = " & "))
writeLines(paste0(rows, " \\\\"), file.path(tab_dir, "sim_glm.tex"))
## Table 2: highly adaptive learner
t2 <- summ[summ$learner == "hapc" & summ$alpha %in% c(0.5, 1.5), ]
rows <- with(t2, paste(n, reps, alpha, sgn(dtruth), sgn(dbias), f3(dsd), f3(dse), f2(dcover), f2(dcover_da),
                       sgn(nbias), f2(ncover), sep = " & "))
writeLines(paste0(rows, " \\\\"), file.path(tab_dir, "sim_hapc.tex"))
## Table 3: treatment scale, logistic regression, n = 2000 (or the largest available)
nmax <- max(summ$n[summ$learner == "glm" & summ$n <= 2000])
t3 <- summ[summ$learner == "glm" & summ$n == nmax, ]
rows <- with(t3, paste(lab[scenario], alpha, f3(tau), f3(tmle_treat), f3(implied), sep = " & "))
writeLines(paste0(rows, " \\\\"), file.path(tab_dir, "sim_treat.tex"))
## sentence fragments quoted in the text
frag <- function(txt, file) cat(txt, "%\n", sep = "", file = file.path(tab_dir, file))
frag(paste(f2(t2$dcover), collapse = " and "), "hapc_cov.tex")
pc <- function(x) paste0(round(100 * x), "\\%")
lr <- unique(summ[summ$learner == "glm" & summ$n %in% c(500, 1000, 2000), c("scenario", "n", "lag_reject")])
fr <- lr$lag_reject[lr$scenario == "fresh"]; pe <- lr$lag_reject[lr$scenario == "persistent"]
frag(paste0("between ", pc(min(fr)), " and ", pc(max(fr))), "lag_renewed.tex")
frag(paste0(pc(pe[1]), ", ", pc(pe[2]), " and ", pc(pe[3])), "lag_persistent.tex")
## levels (supplementary numbers quoted in the text)
print(summ[, c("scenario", "learner", "n", "alpha", "reps", "truth", "bias", "sd", "se", "cover", "cover_da",
               "dtruth", "dbias", "dsd", "dse", "dcover", "dcover_da", "nbias", "ncover")], digits = 3, row.names = FALSE)
print(summ[, c("scenario", "learner", "n", "alpha", "tau", "tmle_treat", "implied", "lag_reject", "min_n_eff", "seconds")],
      digits = 3, row.names = FALSE)

## Figure: dose-response curve, persistent scenario, logistic regression
fg <- summ[summ$learner == "glm" & summ$scenario == "persistent" & summ$n == nmax, ]
pdf(file.path(fig_dir, "sim_curve.pdf"), width = 5.2, height = 3.6)
par(mar = c(4, 4, 0.5, 0.5), mgp = c(2.3, 0.8, 0))
est <- fg$dtruth + fg$dbias; nv <- fg$dtruth + fg$nbias
plot(fg$alpha, fg$dtruth, type = "l", lwd = 2, ylim = range(fg$dtruth, est - 2 * fg$dsd, est + 2 * fg$dsd, nv),
     xlab = expression(alpha), ylab = expression(psi(alpha) - psi(1)))
abline(h = 0, lty = 3)
arrows(fg$alpha, est - z * fg$dsd, fg$alpha, est + z * fg$dsd, angle = 90, code = 3, length = 0.03)
points(fg$alpha, est, pch = 19)
points(fg$alpha, nv, pch = 2); lines(fg$alpha, nv, lty = 2)
legend("topright", bty = "n", lty = c(1, NA, 2), pch = c(NA, 19, 2), lwd = c(2, NA, 1),
       legend = c("truth", "proposed estimator", "analysis ignoring the instrument"))
invisible(dev.off())
