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



source("renv/activate.R")
