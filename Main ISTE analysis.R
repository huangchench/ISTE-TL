################################################################################
# Project: ISTE Transfer Learning for Individual Treatment Effect (ITE)
# Purpose: Main execution pipeline for ISTE-TL and benchmark methods
# Maintainer: Chen Huang
# Date: 2026-03-02
################################################################################

# --- 1. Load Dependencies ---
library(survival)
library(dplyr)
library(tidyr)
library(purrr)
library(gbm)
library(rms)
library(reticulate)
library(TransCox)

# --- 2. Environment Configuration ---
# Modular path management: Use relative paths for better portability on GitHub
BASE_DIR <- getwd() 
source(file.path(BASE_DIR, "ISTE-TL_funcs.R"))

# Configure Python environment for TransCox (TensorFlow backend)
use_condaenv("TransCoxEnvi", required = TRUE)
py_config()
source_python(system.file("python", "TransCoxFunction.py", package = "TransCox"))

# --- 3. Data Acquisition ---
# Ensure these files exist in your working directory or provide specific paths
data_t <- readRDS("data_t.rds")
data_s <- readRDS("data_s.rds")
# Initialize a data frame to store estimation results and individual treatment effects
data_result <- data_t

# --- 4. Parameter Specification ---
# Define baseline covariates for the specific clinical study
covs_X <- c("X1","X2","X3","X4","X5","X6","X7","X8")

# --- 5. User-defined Method Selection ---
# Options: "ISTE-TL", "TransCox", "Target-only", "Source-only", "Stratified-Cox"
# Example: Run only ISTE-TL
data_result <- run_transfer_analysis(method="ISTE-TL", data_t, data_s, covs_X)

# --- 6. Export ---
saveRDS(data_result, file = "data_result.rds")