# ISTE-TL
Transfer Learning for Individual Survival Treatment Effect

# Overview

This repository implements ISTE-TL (Individualized Survival Treatment Effect – Transfer Learning), a unified framework for estimating individualized treatment effects (ITE) in survival settings using transfer learning.

The method is designed for scenarios where:

-The target population has limited sample size;

-A larger source population is available;

-Treatment effects may be partially heterogeneous across populations;

-Direct model transport may induce negative transfer.

This repository provides:
The full execution pipeline for two transfer learning method: ISTE-TL, and TransCox (for estimating ITE); we also provide benchmark implementations:Target-only Cox, Source-only Cox, and Stratified Cox

# R Packages and Python Environment (Required for TransCox)

The following R packages are required:
survival, dplyr,,tidyr, purrr, gbm, rms, reticulate, TransCox

TransCox relies on a TensorFlow backend.
You must configure a conda environment named:
TransCoxEnvi
Inside this environment, ensure TensorFlow and required Python dependencies are installed.

The script uses:
use_condaenv("TransCoxEnvi", required = TRUE)

# Data Requirements
You must provide two datasets:

-data_t.rds — Target population data (with small sample size)

-data_s.rds — Source population data

Each dataset should contain:

-Survival time (Y)

-Event indicator (status)

-Treatment indicator (A)

-Baseline covariates

Example covariates used in the pipeline:

covs_X <- c("X1","X2","X3","X4","X5","X6","X7","X8")

Modify according to your study design.

# Quick Start
1. Load the Environment and R packages (see step 1-2 in file Main ISTE analysis.R)
2. Place your data files:
```r
data_t <- readRDS("data_t.rds")
data_s <- readRDS("data_s.rds")
data_result <- data_t
```

4. Run the main pipeline
```r
covs_X <- c("X1","X2","X3","X4","X5","X6","X7","X8")
# Example: Run only ISTE-TL
data_result <- run_transfer_analysis(method="ISTE-TL", data_t, data_s, covs_X)
saveRDS(data_result, file = "data_result.rds")
```
###
Please note that different methods can be applied to the dataset. The argument method = "ISTE-TL" can take the following options:

"ISTE-TL", "TransCox", "Target-only", "Source-only", and "Stratified-Cox".

The resulting dataset data_result will contain, for each individual in the target population:

-The estimated median survival time if treated (pred_method1);

-The estimated median survival time if untreated (pred_method0);

-The individualized treatment effect defined as the difference between the two (ite_method = pred_method1 - pred_method0).

These estimated individualized survival treatment effects can then be used for subsequent downstream analyses.

