## Exact numerical check of the theory in a fully discrete model (K = 1, seven binary
## variables L0, R0, A0, L1, R1, A1, Y; 128 cells):
##  (1) the efficient influence curve D* = D_h + D_lambda of psi(alpha), with the intervention
##      h_j(1 | F_j) = alpha * P(R(j) = 1 | F_j) depending on P: numerical pathwise derivative
##      versus E[D* s];
##  (2) the exact second-order expansion: Psi(P) - Psi(P0) + P0 D*(P) equals the remainder R2;
##  (3) psi(alpha) is a polynomial of degree <= K + 1 = 2 in alpha;
##  (4) without the D_lambda term the derivative does NOT match.
## Rscript inst/validation/eic_check.R
set.seed(1)
G <- expand.grid(L0 = 0:1, R0 = 0:1, A0 = 0:1, L1 = 0:1, R1 = 0:1, A1 = 0:1, Y = 0:1)
cE <- function(p, f, by) {            # E_p[f | by], one value per row
  k <- if (length(by)) interaction(G[, by, drop = FALSE], drop = TRUE) else factor(rep(1, nrow(G)))
  (tapply(p * f, k, sum) / tapply(p, k, sum))[k]
}
set1 <- function(v, val) { H <- G; H[[v]] <- val; H }
at <- function(f, v, val) {           # f evaluated at row with variable v set to val
  key <- do.call(paste, G); f[match(do.call(paste, set1(v, val)), key)]
}
pieces <- function(p, alpha) {
  lam0 <- cE(p, G$R0, "L0"); lam1 <- cE(p, G$R1, c("L0", "A0", "L1"))
  g0 <- lam0; g1 <- cE(p, G$R1, c("L0", "R0", "A0", "L1"))
  h0 <- ifelse(G$R0 == 1, alpha * lam0, 1 - alpha * lam0)
  h1 <- ifelse(G$R1 == 1, alpha * lam1, 1 - alpha * lam1)
  w0 <- h0 / ifelse(G$R0 == 1, g0, 1 - g0); w1 <- h1 / ifelse(G$R1 == 1, g1, 1 - g1)
  Q1 <- cE(p, G$Y, c("L0", "R0", "A0", "L1", "R1"))
  Q1_1 <- at(Q1, "R1", 1); Q1_0 <- at(Q1, "R1", 0)
  m1 <- alpha * lam1 * Q1_1 + (1 - alpha * lam1) * Q1_0
  Q0 <- cE(p, m1, c("L0", "R0")); Q0_1 <- at(Q0, "R0", 1); Q0_0 <- at(Q0, "R0", 0)
  m0 <- alpha * lam0 * Q0_1 + (1 - alpha * lam0) * Q0_0
  psi <- sum(p * m0)
  Dh <- w0 * w1 * (G$Y - Q1) + w0 * (m1 - Q0) + m0 - psi
  b0 <- Q0_1 - Q0_0
  b1 <- cE(p, w0 * (Q1_1 - Q1_0), c("L0", "A0", "L1"))
  Dl <- alpha * b0 * (G$R0 - lam0) + alpha * b1 * (G$R1 - lam1)
  list(psi = psi, Dh = Dh, Dl = Dl, lam0 = lam0, lam1 = lam1, w0 = w0, w1 = w1,
       Q1 = Q1, Q0 = Q0, m1 = m1, b0 = b0, b1 = b1, Q1_1 = Q1_1, Q1_0 = Q1_0)
}
p0 <- rexp(nrow(G)); p0 <- p0 / sum(p0)
alpha <- 1.3
stopifnot(max(alpha * cE(p0, G$R0, "L0"), alpha * cE(p0, G$R1, c("L0", "A0", "L1"))) < 1)
P0 <- pieces(p0, alpha)
cat("mean of D_h and D_lambda under P0:", signif(sum(p0 * P0$Dh), 3), signif(sum(p0 * P0$Dl), 3), "\n")

## (1) pathwise derivative along 5 random scores
res <- t(sapply(1:5, function(i) {
  s <- rnorm(nrow(G)); s <- s - sum(p0 * s)
  path <- function(e) { q <- p0 * (1 + e * s); pieces(q / sum(q), alpha)$psi }
  e <- 1e-5
  c(numerical = (path(e) - path(-e)) / (2 * e), with_D_lambda = sum(p0 * (P0$Dh + P0$Dl) * s),
    without_D_lambda = sum(p0 * P0$Dh * s))
}))
print(signif(res, 6))
cat("max |numerical - E[D* s]|        :", signif(max(abs(res[, 1] - res[, 2])), 3), "\n")
cat("max |numerical - E[D_h s]| (wrong):", signif(max(abs(res[, 1] - res[, 3])), 3), "\n")

## (2) exact remainder: Psi(P) - Psi(P0) + P0 D*(P) = R2(P, P0)
remainder <- function(p) {
  P <- pieces(p, alpha)
  lhs <- P$psi - P0$psi + sum(p0 * (P$Dh + P$Dl))
  ## first part: sum_j E0[(W_j - W0_j)(Qtilde0_j - Q_j)], weights W use h(P) and g(P); W0 uses h(P) and g0
  g0_0 <- cE(p0, G$R0, "L0"); g1_0 <- cE(p0, G$R1, c("L0", "R0", "A0", "L1"))
  h0 <- ifelse(G$R0 == 1, alpha * P$lam0, 1 - alpha * P$lam0)
  h1 <- ifelse(G$R1 == 1, alpha * P$lam1, 1 - alpha * P$lam1)
  w0t <- h0 / ifelse(G$R0 == 1, g0_0, 1 - g0_0); w1t <- h1 / ifelse(G$R1 == 1, g1_0, 1 - g1_0)
  Qt1 <- cE(p0, G$Y, c("L0", "R0", "A0", "L1", "R1"))            # E0[m_2 | R1, H1], m_2 = Y
  Qt0 <- cE(p0, P$m1, c("L0", "R0"))                              # E0[m_1(P) | R0, H0]
  r_dr <- sum(p0 * (P$w0 * P$w1 - w0t * w1t) * (Qt1 - P$Q1)) + sum(p0 * (P$w0 - w0t) * (Qt0 - P$Q0))
  ## second part: alpha * sum_j E0[(bmix_j - b_j)(lam_j - lam0_j)]
  ## bmix_j: true (P0) weights with h(P0) before j, true (P0) regressions with h(P) after j
  h0_00 <- ifelse(G$R0 == 1, alpha * P0$lam0, 1 - alpha * P0$lam0)
  w0_00 <- h0_00 / ifelse(G$R0 == 1, g0_0, 1 - g0_0)
  Qt1_1 <- at(Qt1, "R1", 1); Qt1_0 <- at(Qt1, "R1", 0)
  bmix1 <- cE(p0, w0_00 * (Qt1_1 - Qt1_0), c("L0", "A0", "L1"))
  m1mix <- alpha * P$lam1 * Qt1_1 + (1 - alpha * P$lam1) * Qt1_0  # true Q1, estimated h_1
  Q0mix <- cE(p0, m1mix, c("L0", "R0"))
  bmix0 <- at(Q0mix, "R0", 1) - at(Q0mix, "R0", 0)
  r_l <- alpha * sum(p0 * (bmix0 - P$b0) * (P$lam0 - P0$lam0)) +
         alpha * sum(p0 * (bmix1 - P$b1) * (P$lam1 - P0$lam1))
  c(lhs = lhs, R2 = r_dr + r_l)
}
cat("\nExact expansion (lhs = Psi(P) - Psi(P0) + P0 D*(P); R2 = closed-form remainder):\n")
for (e in c(0.3, 0.1, 0.03)) {
  s <- rnorm(nrow(G)); s <- s - sum(p0 * s); q <- p0 * exp(e * s); q <- q / sum(q)
  r <- remainder(q); cat(sprintf("  eps = %.2f  lhs = % .3e  R2 = % .3e  diff = %.1e\n", e, r[1], r[2], abs(r[1] - r[2])))
}

## (3) polynomial of degree <= 2 in alpha: third differences vanish
a <- seq(0.2, 1.4, by = 0.3); ps <- sapply(a, function(x) pieces(p0, x)$psi)
cat("\npsi(alpha):", signif(ps, 6), "\nthird differences:", signif(diff(ps, differences = 3), 3), "\n")
