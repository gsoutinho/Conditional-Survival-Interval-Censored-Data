# Conditional Survival Estimation with Interval-Censored Event Times

This repository contains the R code used in the simulation studies and real-data application for the paper on nonparametric conditional survival estimation with interval-censored event times.

Two settings are considered: (i) the survival function of an interval-censored response conditional on a continuous covariate, and (ii) the survival function of a terminal event conditional on an interval-censored intermediate event falling within a given time window $(a,b]$. For each setting, the proposed Turnbull-based estimator is compared with a midpoint-imputation benchmark based on the Beran estimator.

## Contents

- `cond_surv_ic.R`: conditional survival estimation given a continuous covariate. The argument `method` selects between the kernel-weighted conditional Turnbull estimator (`"Turnbull"`) and the midpoint-based Beran estimator (`"Beran"`).
- `cond_aZb_surv_ic.R`: survival conditional on an interval-censored intermediate event, $S(t \mid a < Z \le b)$. The argument `method` selects between the Turnbull-based reconstruction estimator (`"Turnbull"`) and the midpoint-based benchmark (`"Beran"`).
- `results.R`: scripts reproducing all results reported in the manuscript, including the simulation studies and the application to the breast cancer data.

## Usage

```r
source("cond_surv_ic.R")
source("cond_aZb_surv_ic.R")
```

See `results.R` for complete examples.

## Related software

These functions will be included in the R package `icSurvML`, which provides a unified framework for the analysis of interval-censored survival data.
