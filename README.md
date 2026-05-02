# Composite Categorical Regression (CCR)

This repository contains R code for the Composite Categorical Regression (CCR) method used in the manuscript. The code supports reproducibility of the proposed method, including simulation studies, an additional bootstrap-based type I error evaluation, and the BRFSS application. Only the CCR method is included; competing penalization methods are not provided in this repository.

## Files

| File | Purpose |
|---|---|
| `CCR_simulation_unadjusted.R` | Runs the unadjusted simulation study for Gaussian and binary outcomes. |
| `CCR_simulation_adjusted.R` | Runs the adjusted simulation study for Gaussian and binary outcomes. |
| `CCR_type1_bootstrap.R` | Runs the bootstrap-based type I error simulation under the null setting. |
| `CCR_application_logistic_estimation.R` | Fits the survey-weighted CCR logistic regression model in the BRFSS application. |
| `CCR_application_logistic_bootstrap.R` | Performs bootstrap resampling for inference in the BRFSS application. |
| `CCR_application_bootstrap_summary.R` | Summarizes the bootstrap results and exports point estimates and confidence intervals. |
| `data_1923_ungrp.RData` | Preprocessed BRFSS analytic dataset used by the application scripts. |

## Software requirements

The code was written in R. The main packages used across the scripts include:

```r
MASS
pROC
openxlsx
dplyr
boot
tidyverse
glmnet
writexl
```

Package usage differs by script; not every package is required for every file.

## Simulation studies

The main simulation studies are contained in:

```r
source("CCR_simulation_unadjusted.R")
source("CCR_simulation_adjusted.R")
```

Both scripts consider Gaussian and binary outcomes under three exposure-correlation settings:

```r
rho_values <- c(0.2, 0.5, 0.8)
```

The additional bootstrap-based type I error simulation can be run separately:

```r
source("CCR_type1_bootstrap.R")
```

Simulation outputs are written to the working directory. These outputs summarize exposure selection, prediction performance, coefficient estimation, zero-weight frequencies, and bootstrap-based type I error evaluation for the CCR method.

## BRFSS application

The BRFSS application is organized into three scripts and should be run in the following order:

```r
source("CCR_application_logistic_estimation.R")
source("CCR_application_logistic_bootstrap.R")
source("CCR_application_bootstrap_summary.R")
```

The application scripts use the preprocessed BRFSS analytic dataset included with the repository. The expected input file is:

```text
data_1923_ungrp.RData
```

The file contains the analytic data object used by the application scripts. Output files, including estimated CCR weights, fitted model objects, bootstrap replicates, and confidence-interval summaries, are saved to the working directory.

## Method summary

CCR constructs a composite exposure index from dummy variables representing categorical exposures. The dummy-specific CCR weights are constrained to be nonnegative and normalized to sum to one. In the BRFSS application, this index is included in a survey-weighted logistic regression model along with adjustment covariates.

## Notes

The simulation scripts can be run without external data. The BRFSS application scripts use the included preprocessed analytic dataset and may take additional time when running the bootstrap procedure.
