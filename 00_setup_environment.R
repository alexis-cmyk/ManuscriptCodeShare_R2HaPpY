# 00_setup_environment.R

# Ensure yaml is available (needed by rmarkdown/knitr and renv's dependency parsing)
if (!requireNamespace("yaml", quietly = TRUE)) {
  install.packages("yaml", repos = "https://cran.r-project.org")
}

# Ensure renv is available globally
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cran.r-project.org")
}


# ---------------------------------------------------------
# Install BiocManager if missing
# ---------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# ---------------------------------------------------------
# Set the correct Bioconductor version for R 4.4.2
# (R 4.4.x → Bioc 3.19)
# ---------------------------------------------------------
# BiocManager::install(version = "3.19")
# 
# # Ensure repos are set correctly (renv-safe)
# options(repos = BiocManager::repositories(version = "3.19"))

# ---------------------------------------------------------
# Helper: install package if missing
# ---------------------------------------------------------
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

# ---------------------------------------------------------
# Install Bioconductor packages
# ---------------------------------------------------------
install_if_missing("limma")
install_if_missing("biomaRt")


# Load package
library(limma)
library(biomaRt)
message("Setup complete. Now restart R in this project and run renv::restore().")

