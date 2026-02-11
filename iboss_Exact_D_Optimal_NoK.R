# A Column Generation Approach to Exact Experimental Design
# Authors: Selin Ahipasaoglu, Stefano Cipolla, Jacek Gondzio
# arXiv ID: 2507.03210 (July 2025)
#
# R translation of matlab/iboss_Exact_D_Optimal_NoK.m
# Uses the original R IBOSS package (iboss.od with compiled C++ backend)
#
# Usage:
#   From the project root directory:
#     Rscript iboss_Exact_D_Optimal_NoK.R
#   Or within an R session:
#     source("iboss_Exact_D_Optimal_NoK.R")
#
# Prerequisites:
#   install.packages("R.matlab")  # for reading .mat files
#   devtools::install(".")        # install the IBOSS package from source

# --- Load required packages ---
library(R.matlab)

# Load the IBOSS package.
# Option 1: If installed as a package
# library(IBOSS)
# Option 2: If running from the IBOSS project root (development mode)
devtools::load_all(".")

# --- Configuration ---
# The path on which all the data are:
problems_path <- "/scratch/sc9c23/Dataset/Gen1"

# Find all the problems and store their names
files <- sort(list.files(problems_path, pattern = "^Generator_10M_.*\\.mat$",
                         full.names = FALSE))

# N offsets
NN <- c(0)

# Total number of experiments
expN <- length(files) * length(NN)

# --- Initialize results data frame ---
resultsIBOSS <- data.frame(
  Problem = character(expN),
  n       = numeric(expN),
  m       = numeric(expN),
  Time    = numeric(expN),
  N       = numeric(expN),
  Obj     = numeric(expN),
  stringsAsFactors = FALSE
)

# --- Main experiment loop ---
exp_num <- 0

for (k in seq_along(files)) {
  # Load .mat file and extract the data matrix
  mat_data <- readMat(file.path(problems_path, files[k]))
  X <- mat_data$genmatrix

  cat(files[k], "\n")

  # X is (n x m): n = number of features (rows), m = number of observations (columns)
  n <- nrow(X)
  m <- ncol(X)

  for (ll in seq_along(NN)) {
    cat(ll, "\n")

    N <- 2 * n + NN[ll]
    exp_num <- exp_num + 1

    # Run IBOSS: transpose X so that iboss.od receives (observations x features)
    # Y = ones(m,1) is a dummy response (only index selection matters here)
    time_start <- proc.time()
    fit <- iboss.od(t(X), rep(1, m), N)
    time_IBOSS <- (proc.time() - time_start)[["elapsed"]]

    # Select the columns of X corresponding to the chosen observation indices
    X_S <- X[, fit$index, drop = FALSE]

    # D-optimality objective: log(det(X_S %*% t(X_S)))
    obj_IBOSS <- log(det(X_S %*% t(X_S)))

    # Record results
    resultsIBOSS$Problem[exp_num] <- files[k]
    resultsIBOSS$n[exp_num]       <- n
    resultsIBOSS$m[exp_num]       <- m
    resultsIBOSS$Time[exp_num]    <- time_IBOSS
    resultsIBOSS$N[exp_num]       <- N
    resultsIBOSS$Obj[exp_num]     <- obj_IBOSS
  }
}

# --- Write results to CSV ---
write.csv(resultsIBOSS, file.path("Results", "R_IBOSS_10M_NoK.csv"), row.names = FALSE)

cat("Results written to Results/IBOSS_10M_NoK.csv\n")
