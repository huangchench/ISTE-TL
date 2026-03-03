# ISTE-TL
Transfer Learning for Individual Survival Treatment Effect

# Overview

This repository implements ISTE-TL (Individualized Survival Treatment Effect – Transfer Learning), a unified framework for estimating individualized treatment effects (ITE) in survival settings using transfer learning.

The method is designed for scenarios where:

- The target population has limited sample size;
- A larger source population is available;
- Treatment effects may be partially heterogeneous across populations;
- Direct model transport may induce negative transfer.

This repository provides:

- A complete execution pipeline for two transfer learning approaches (**ISTE-TL** and **TransCox**) for individualized treatment effect estimation;
- Standard benchmark models, including **Target-only Cox**, **Source-only Cox**, and **Stratified Cox**.

# R Packages and Python Environment (Required for TransCox)

The following R packages are required:
- `survival`  
- `dplyr`  
- `tidyr`  
- `purrr`  
- `gbm`  
- `rms`  
- `reticulate`  
- `TransCox`


**Python Environment (only required if using TransCox):**  
TransCox relies on a TensorFlow backend. You must configure a conda environment named:
```r
TransCoxEnvi
```
Inside this environment, ensure TensorFlow and required Python dependencies are installed.

The script uses:
```r
use_condaenv("TransCoxEnvi", required = TRUE)
```
**Note:** If you do not plan to use the `TransCox` method and only run ISTE-TL or benchmark Cox models, configuring the Python environment is **not required**.

# Data Requirements
You must provide two datasets:

- data_t: Target population data (with small sample size)
- data_s: Source population data

Each dataset should contain:

- Survival time (Y)
- Event indicator (status)
- Treatment indicator (A)
- Baseline covariates

   Example covariates used in the pipeline:
  
   `covs_X <- c("X1","X2","X3","X4","X5","X6","X7","X8")`
  
   Modify according to your study design.

# Quick Start
All analyses are executed through the main script:
`Main ISTE analysis.R`

**1. Load the Environment and R packages**

Follow Steps 1–2 in `Main ISTE analysis.R` to:
- Load required R packages
- Configure the Python environment for TransCox (if used)
- Source the core function file ISTE-TL_funcs.R
  
**2. Place your dataset:**

```r
data_t <- readRDS("data_t.rds")
data_s <- readRDS("data_s.rds")
data_result <- data_t
```
Ensure that both datasets contain:
- Survival time
- Event indicator
- Treatment indicator
- Baseline covariates

**3. Specify Covariates**

Define the baseline covariates used in the analysis:
```r
covs_X <- c("X1","X2","X3","X4","X5","X6","X7","X8")
```
Modify this according to your study.

**4. Run the Transfer Analysis**
```r
# Example: Run only ISTE-TL
data_result <- run_transfer_analysis(method="ISTE-TL", data_t, data_s, covs_X)
```
Please note that different methods can be applied to the dataset. The argument `method = "ISTE-TL"` can take the following options:

`"ISTE-TL"`, `"TransCox"`, `"Target-only"`, `"Source-only"`, and `"Stratified-Cox"`.

**5. Output**

The resulting dataset `data_result` will contain, for each individual in the target population:

- The estimated median survival time if treated (`pred_method1`);
- The estimated median survival time if untreated (`pred_method0`);
- The individualized treatment effect defined as the difference between the two (`ite_method = pred_method1 - pred_method0`).

You may then proceed with downstream analyses based on the estimated individualized survival treatment effects.

