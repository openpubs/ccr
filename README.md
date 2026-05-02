# Composite Categorical Regression (CCR)

This repository contains R code for the Composite Categorical Regression (CCR) method used in the manuscript. The code supports reproducibility of the simulation studies and the BRFSS application. Only the CCR method is included; competing penalization methods are not included in this repository.

## Files

| File | Purpose |
|---|---|
| `CCR_simulation_unadjusted.R` | Runs the unadjusted simulation study for Gaussian and binary outcomes. |
| `CCR_simulation_adjusted.R` | Runs the adjusted simulation study for Gaussian and binary outcomes. |
| `CCR_simulation_test.R` | Runs an additional simulation check for the CCR method. |
| `CCR_application_logistic_estimation.R` | Fits the survey-weighted CCR logistic regression model for the BRFSS application. |
| `CCR_application_logistic_bootstrap.R` | Performs bootstrap resampling for inference in the BRFSS application. |
| `CCR_application_bootstrap_summary.R` | Summarizes the bootstrap results and exports point estimates and confidence intervals. |

## Software requirements

The code was written in R. The required packages may vary by script, but the main packages used are:

```r
MASS
pROC
openxlsx
boot
dplyr
tidyverse
glmnet
writexl
```

## Simulation studies

The main simulation scripts can be run directly:

```r
source("CCR_simulation_unadjusted.R")
source("CCR_simulation_adjusted.R")
```

An additional simulation check can be run using:

```r
source("CCR_simulation_test.R")
```

The simulation scripts evaluate the CCR method under Gaussian and binary outcomes and different correlation settings among categorical exposures. Output files are saved to the working directory.

## BRFSS application

The BRFSS application is organized into three scripts and should be run in the following order:

```r
source("CCR_application_logistic_estimation.R")
source("CCR_application_logistic_bootstrap.R")
source("CCR_application_bootstrap_summary.R")
```

The application scripts fit the survey-weighted CCR logistic regression model, perform bootstrap resampling, and summarize the estimated CCR weights and confidence intervals. The required preprocessed BRFSS analytic dataset should be placed in the working directory before running the application scripts.

## Method summary

CCR constructs a composite exposure index from dummy variables representing categorical exposures. The dummy-specific CCR weights are constrained to be nonnegative and normalized to sum to one. In the BRFSS application, the CCR index is included in a survey-weighted logistic regression model with adjustment covariates.

## Notes

The simulation scripts can be run without external data. The BRFSS application scripts require the preprocessed analytic dataset and may take additional time when running the bootstrap procedure.
