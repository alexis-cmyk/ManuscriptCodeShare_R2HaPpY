---

# Phosphotyrosine Proteomics Using a Scalable SH2 Superbinder Enrichment Strategy

[![R Version](https://img.shields.io/badge/R-4.4.2-blue.svg)]()
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.22-green.svg)]()
[![Reproducibility](https://img.shields.io/badge/Reproducible-renv-orange.svg)]()
[![License](https://img.shields.io/badge/License-MIT-lightgrey.svg)]()

This repository contains all analysis scripts, figure-generation pipelines, and processed datasets supporting:

**“Phosphotyrosine Proteomics Using a Scalable SH2 Superbinder Enrichment Strategy.”**
**[bioRxiv Preprint](https://www.biorxiv.org/content/10.1101/2025.05.14.653984v1)**

<img src="./README_VisualAbstract.png" width="650">

---

# Table of Contents

1. [Repository Overview](#repository-overview)
2. [Installation & Execution](#installation--execution)

   * [1. Download Input Files](#1-download-input-files)
   * [2. Install Required R Version](#2-install-required-r-version)
   * [3. Initialize Environment](#3-initialize-environment)
   * [4. Restore Software Environment](#4-restore-software-environment)
   * [5. Generate All Figures](#5-generate-all-figures)
3. [Interactive Dataset](#interactive-dataset)
4. [Experimental Protocols](#experimental-protocols)
5. [Data Availability](#data-availability)
6. [Contact](#contact)

---

# Repository Overview

* R scripts for phosphopeptide, phosphoproteome, and proteome data analysis
* Code to generate all main and supplementary figures
* Processed data tables required for reproducibility
* Expected output directory structure

Large raw datasets (Comet results, AScore tables, FASTA files) are provided separately as GitHub Release assets.

---

# Installation & Execution

## 1. Download Input Files

Download the Release archive containing:

* **raw_data/**
* **modified_data/**
* **output/**

Place all directories in the project root.

---

## 2. Install Required R Version

This project was developed and tested using:

* **R 4.4.2**
* **Bioconductor 3.22**

Matching these versions ensures compatibility with the **renv** lockfile.

- [Download R version 4.4.2 here](https://https://cran.r-project.org/bin/windows/base/old/4.4.2/)
---

## 3. Initialize Environment

From the project directory:

```r
source("00_setup_environment.R")
```

Restart R after running this script.

---

## 4. Restore Software Environment

```r
renv::restore()
```

This installs all CRAN and Bioconductor packages pinned in the `renv.lock` file.

---

## 5. Generate All Figures

```r
source("00_master_run_all_figures.R")
```

Outputs:

* Figures → **output/** subdirectories
* Intermediate tables → **modified_data/** subdirectories

---

# Interactive Dataset

A browser-based interactive version of the R2HaPpY dataset is available at:

➡ **[https://r2happy.gs.washington.edu](https://r2happy.gs.washington.edu)**

This site supports phosphosite-level browsing, replicate comparison, gene filtering, and functional annotation.

---

# Experimental Protocols

Detailed laboratory procedures—including SH2 Superbinder expression, bead preparation, enrichment, and LC-MS/MS acquisition—are available at:

➡ **[Protocols.io](https://www.protocols.io/private/48FD48C3977A11F0A3900A58A9FEAC02)**

---

# Data Availability

Raw phosphoproteomics and proteome CSVs:

**[Download Raw Data (v1.0.1)](https://github.com/alexis-cmyk/ManuscriptCodeShare_R2HaPpY/releases/tag/v1.0.1-InputData_and_OutputFolders)**

Each archive contains Comet outputs, AScore assignments, peptide redundancy tables, and required metadata.

Processed datasets required by figure scripts are included under **modified_data/**.

---

# Contact

**Alexis Chang**
Villen Laboratory, Genome Sciences
University of Washington, Seattle

Lab website: **[https://villenlab.gs.washington.edu/wordpress/](https://villenlab.gs.washington.edu/wordpress/)**

---


