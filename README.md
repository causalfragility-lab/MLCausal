 # MLCausal

**Causal Inference Methods for Multilevel and Clustered Data**

<!-- badges: start -->
[![R-CMD-check](https://github.com/causalfragility-lab/MLCausal/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/causalfragility-lab/MLCausal/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/MLCausal)](https://CRAN.R-project.org/package=MLCausal)
<!-- badges: end -->

`MLCausal` provides a structured workflow for estimating causal effects in
clustered observational data, such as students within schools, patients within
hospitals, or employees within firms.

The package integrates the full causal design pipeline for multilevel settings:
cluster-aware propensity score estimation, inverse-probability weighting,
within-cluster matching, covariate balance diagnostics at both the individual
and cluster levels, overlap assessment, outcome modelling with cluster-robust
standard errors, and sensitivity analysis for unmeasured confounding.

---

## Installation

```r
# From CRAN
install.packages("MLCausal")

# Development version from GitHub
# install.packages("remotes")
remotes::install_github("causalfragility-lab/MLCausal")

```

Required dependencies installed automatically: `sandwich`, `lmtest`,
`ggplot2`, `rlang`.

---

## Workflow

```
simulate_ml_data()
        |
     ml_ps()                  cluster-aware propensity score estimation
      /     \
ml_weight() ml_match()        inverse-probability weighting  OR  matching
      \     /
   balance_ml()               balance diagnostics (individual + cluster levels)
   plot_overlap_ml()          overlap / positivity assessment
        |
estimate_att_ml()             outcome model with cluster-robust SEs
        |
     sens_ml()                tipping-point sensitivity analysis
```

---

## Functions

| Function | Description |
|---|---|
| `simulate_ml_data()` | Generate clustered observational data with a known data-generating process |
| `ml_ps()` | Cluster-aware propensity score estimation: Mundlak (default), fixed-effects, or single-level |
| `ml_weight()` | ATT / ATE inverse-probability weights, stabilised and with optional trimming |
| `ml_match()` | Within-cluster nearest-neighbour matching with a dual-balance composite distance |
| `balance_ml()` | Standardised mean differences at both the individual and cluster-mean level |
| `plot_overlap_ml()` | Propensity score overlap plot, overall or faceted by cluster |
| `estimate_att_ml()` | Weighted linear outcome model with cluster-robust (HC1) standard errors |
| `sens_ml()` | Tipping-point sensitivity analysis for omitted cluster-level confounding |

---

## Quick start

```r
library(MLCausal)

# -- 1. Simulate clustered data ----------------------------------------------
dat <- simulate_ml_data(
  n_clusters   = 30,
  cluster_size = 20,
  seed         = 42
)

# -- 2. Propensity score estimation (Mundlak method) -------------------------
ps <- ml_ps(
  data       = dat,
  treatment  = "z",
  covariates = c("x1", "x2", "x3"),
  cluster    = "school_id",
  method     = "mundlak"
)

# -- 3. Overlap (positivity) check -------------------------------------------
plot_overlap_ml(ps)

# -- 4. Inverse-probability weighting ----------------------------------------
dat_w <- ml_weight(
  ps_fit    = ps,
  estimand  = "ATT",
  stabilize = TRUE,
  trim      = 10
)

# -- 5. Balance diagnostics --------------------------------------------------
bal <- balance_ml(
  data       = dat_w,
  treatment  = "z",
  covariates = c("x1", "x2", "x3"),
  cluster    = "school_id",
  weights    = "weights"
)
print(bal)

# -- 6. Treatment effect estimation ------------------------------------------
est <- estimate_att_ml(
  data       = dat_w,
  outcome    = "y",
  treatment  = "z",
  cluster    = "school_id",
  covariates = c("x1", "x2", "x3"),
  weights    = "weights"
)
print(est)

# -- 7. Sensitivity analysis -------------------------------------------------
sens <- sens_ml(
  estimate = est$estimate,
  se       = est$se
)

# First confounder strength that would nullify significance
sens[sens$crosses_null, ][1, ]
```

---

## Matching path (alternative to weighting)

```r
# Within-cluster matching with dual-balance penalty (lambda = 1)
matched <- ml_match(
  ps_fit  = ps,
  ratio   = 1,
  caliper = 0.5,
  lambda  = 1
)

# Balance on the matched sample
bal_m <- balance_ml(
  data       = matched$data_matched,
  treatment  = "z",
  covariates = c("x1", "x2", "x3"),
  cluster    = "school_id",
  weights    = "match_weight"
)
print(bal_m)

# Effect estimate on the matched sample
est_m <- estimate_att_ml(
  data      = matched$data_matched,
  outcome   = "y",
  treatment = "z",
  cluster   = "school_id",
  weights   = "match_weight"
)
print(est_m)
```

The `lambda` argument controls the dual-balance penalty. `lambda = 0`
recovers standard within-cluster propensity score matching; larger values
increasingly penalise matches that worsen cluster-mean covariate balance.

---

## Citation

If you use `MLCausal` in published research, please cite:

```
Hait, S. (2026). MLCausal: Causal inference methods for multilevel and clustered data. R package.
```

---

## License
MIT © 2026 causalfragility-lab
