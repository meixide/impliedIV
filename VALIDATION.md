# impliedIV — validation status

All checks use `simulate_lii()`: prescription → adherence, K = 2, and an
unmeasured confounder of adherence and outcome. The prescription law depends
on the past only through `F_j`. Censoring after each visit depends on `L(j)`;
in this DGP about 23% of subjects are censored before the outcome. Truth:
Monte Carlo with 10^6 draws of the intervened system, without censoring.
Scripts are in `inst/validation/`; results are in
`inst/validation/results/`.

## 1. Unit tests (`devtools::test()`)

28 expectations, all passing (≈ 4 s).
- Data preparation and censoring.
- Rejection of non-binary instruments and of NA covariates.
- Rescaling of continuous outcomes.
- At α = 1 without censoring, the estimate is within 0.5 SE of the sample mean.
- The influence curve has mean 0.
- Agreement with the truth at n = 3000.
- The implied intervention is affine in α and equals π0 at α = 0.
- The treatment-count outcome and calibration.
- The bootstrap.
- `preference_instrument()` against a brute-force computation.
- `lag_sufficiency_test()`.

`R CMD check --no-manual`: no errors. The only warning came from the
sandbox's lack of network access to CRAN (unstated test dependencies `glue`
and `waldo`, now declared in Suggests).

## 2. Influence-curve CIs, `learner_glm` (`coverage_glm.R 200 1000`)

200 replications, n = 1000.

| α | truth | bias | sd | mean SE | coverage | contrast with α = 1: bias | sd | mean SE | coverage |
|---|---|---|---|---|---|---|---|---|---|
| 0.50 | 0.3172 | −0.0004 | 0.0222 | 0.0211 | 0.925 | 0.0004 | 0.0134 | 0.0128 | 0.920 |
| 0.75 | 0.3028 | −0.0011 | 0.0181 | 0.0174 | 0.935 | −0.0002 | 0.0061 | 0.0058 | 0.920 |
| 1.00 | 0.2882 | −0.0009 | 0.0166 | 0.0160 | 0.950 | — | — | — | — |
| 1.25 | 0.2756 | −0.0020 | 0.0174 | 0.0165 | 0.935 | −0.0011 | 0.0059 | 0.0055 | 0.915 |
| 1.50 | 0.2631 | −0.0024 | 0.0202 | 0.0189 | 0.925 | −0.0015 | 0.0121 | 0.0112 | 0.930 |

- Bias is negligible.
- The influence-curve SE is about 95% of the Monte Carlo sd, so coverage is
  0.92–0.95. The Monte Carlo error of a coverage with 200 replications is
  about ±0.03.
- The shortfall is the one expected: the truth here is the parameter defined
  with the *true* instrument intensity λ, whose influence curve has an extra
  term (λ̂ is estimated). That term is ignored by the influence-curve SE,
  which targets the data-adaptive parameter.

## 3. Bootstrap CIs, `learner_glm` (`coverage_bootstrap_glm.R 100 100 1000`)

PENDING_BOOTSTRAP

## 4. End-to-end run with `learner_hapc` (`run_hapc.R 1000`)

n = 1000, K = 2, α ∈ {0.5, 1, 1.5}; hapc 2.7.1 (`cv.hapc`, binomial,
`norm = "1"`, `max_degree = 2`, `n^0.7` PCs, 5 folds).

| α | estimate | 95% CI | truth |
|---|---|---|---|
| 0.5 | 0.344 | [0.300, 0.389] | 0.317 |
| 1.0 | 0.291 | [0.259, 0.323] | 0.289 |
| 1.5 | 0.253 | [0.215, 0.291] | 0.262 |

- Every lambda was interior to its grid; no fit failed.
- First stage: π1 − π0 had median 0.45.
- **Time: 403 s**, almost all of it in the five pooled nuisance fits (≈ 2 min
  each with 2735 rows). The outcome regressions take about 1 s each.
- hapc cost grows roughly like n³ in the rows of a fit, so the default cap
  `max_n = 4000` rows per fit means ≈ 5–6 min per pooled fit on a laptop.
- This is a single run: coverage with hapc has not been studied.

## Not validated yet

- Coverage with hapc, and of the bootstrap with hapc.
- K > 2, rare instruments, strong censoring, continuous outcomes in
  simulation.
- `binarize_instrument()` (by design only a sensitivity analysis).
- Interpretation of `lag_sufficiency_test()` with overlapping windows: its
  power, and the size under heavy overlap.
