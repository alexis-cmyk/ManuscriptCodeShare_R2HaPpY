#Force known good repositories for renv boostrap
options(repos = c(CRAN = "https://cran.r-project.org"))

#Ensure renv uses a working download method. Choose the one that works for install.packages() on your computer
options(renv.download.method = "libcurl")   # OR "wininet"


# Use CRAN + Bioconductor 3.22 repos
options(repos = c(
  CRAN       = "https://cran.r-project.org",
  BioCsoft   = "https://bioconductor.org/packages/3.22/bioc",
  BioCann    = "https://bioconductor.org/packages/3.22/data/annotation",
  BioCexp    = "https://bioconductor.org/packages/3.22/data/experiment"
))

# Tell renv: do NOT use its Bioconductor integration (no BiocManager dance).
options(renv.config.bioconductor = FALSE)

# Tell renv not to query BiocManager/BiocVersion
options(
  renv.bioconductor.repos = c(
    BioCsoft      = "https://bioconductor.org/packages/3.22/bioc",
    BioCann       = "https://bioconductor.org/packages/3.22/data/annotation",
    BioCexp       = "https://bioconductor.org/packages/3.22/data/experiment",
    BioCworkflows = "https://bioconductor.org/packages/3.22/workflows"
  )
)


# Ensure 'yaml' is available (needed by rmarkdown/knitr and renv dependency parsing)
if (!requireNamespace("yaml", quietly = TRUE)) {
  message("Installing 'yaml' (needed for this project)...")
  try(
    install.packages("yaml", repos = "https://cran.r-project.org"),
    silent = TRUE
  )
}

# Disable renv Bioconductor integration completely
options(renv.bioconductor.enabled = FALSE)
options(renv.config.bioconductor = FALSE)

if (file.exists("renv/activate.R"))
  source("renv/activate.R")
