# 00_setup_environment.R

if (!requireNamespace("yaml", quietly = TRUE)) {
  install.packages("yaml", repos = "https://cran.r-project.org")
}

install.packages("renv", repos = "https://cran.r-project.org")
renv::restore()

### One-time setup

# From the project root in R:
#   
# ```r
# source("00_setup_environment.R")
# ```