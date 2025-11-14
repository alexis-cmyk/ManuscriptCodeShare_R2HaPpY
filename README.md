
# Phosphotyrosine Proteomics Using a Scalable SH2 Superbinder Enrichment Strategy

This repository contains all analysis scripts, figure-generation code, and search results used to reproduce the results of the manuscript:

**“Phosphotyrosine Proteomics Using a Scalable SH2 Superbinder Enrichment Strategy.”**

The codebase includes:
- R scripts used to process phosphoproteomics, global phosphoproteome, and proteome datasets  
- Figure-generation pipelines for all Main Figures and Supplementary Figures  
- Shiny-based interactive visualization modules  
- Processed data tables required for reproducing plots  

Large raw CSV files (Comet search results, AScore outputs, redundancy tables, etc.) are not stored in the repository due to GitHub file-size limits, but are provided as release assets (see “Data availability” below).

---

## Interactive Dataset

An interactive, browser-based version of the full R2HaPpY dataset is hosted at:

➡ **https://r2happy.gs.washington.edu**

This site allows exploration of phosphorylation patterns, gene-level filtering, replicate variability, and downstream functional annotation analyses.

---

## Protocols

The full experimental workflow, including SH2 Superbinder expression, bead preparation, enrichment, and LC-MS/MS acquisition is documented on protocols.io:

➡ **[https://www.protocols.io/private/48FD48C3977A11F0A3900A58A9FEAC02]**

---

## Data availability

Raw phosphoproteomics and proteome CSVs for all figures are provided as release assets:

- **[Download R2HaPpY raw data (v1.0)](https://github.com/alexis-cmyk/ManuscriptCodeShare_R2HaPpY/releases/tag/v1.0.0-data)**

Each archive contains the full set of raw search engine outputs (Comet), AScore assignments, peptide-level redundancy tables, and meta-data required to reproduce the manuscript analyses.

Processed and intermediate data tables used in the figure-generation scripts (e.g., scaled intensities, summarized peptide tables, cluster assignments) are included directly in this repository under the `modified_data/` subdirectory.

---

## Reproducing figures

All figures in the manuscript can be regenerated using:

- `R/00_master_run_all_figures.R`
- Individual figure-specific scripts in `Scripts`
- A Shiny application enabling interactive re-analysis and visualization (`shiny/`)

Dependencies are managed using **renv**.  
To reproduce the full environment:


```{r}
renv::restore()
```

