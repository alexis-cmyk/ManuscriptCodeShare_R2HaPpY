# 00_setup_environment.R

# Ensure yaml is available (needed by rmarkdown/knitr and renv's dependency parsing)
if (!requireNamespace("yaml", quietly = TRUE)) {
  install.packages("yaml", repos = "https://cran.r-project.org")
}

# Ensure renv is available globally
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cran.r-project.org")
}

message("Setup complete. Now restart R in this project and run renv::restore().")