#Phosphotyrosine Proteomics Using a Scalable SH2 Superbinder Enrichment Strategy


This repository contains all analysis scripts, figure-generation code, and processed data supporting:

**“Phosphotyrosine Proteomics Using a Scalable SH2 Superbinder Enrichment Strategy.”**
**[BioRxiv Preprint](https://www.biorxiv.org/content/10.1101/2025.05.14.653984v1)**

<img src="./README_VisualAbstract.png" width="650">

## 🔗 Quick Navigation

- Installation & Execution

- Interactive Dataset

- Experimental Protocols

- Data Availability

- Contact

## 📦 Repository Overview

- R scripts for phosphopeptide, phosphoproteome, and proteome data processing

- Code for generating all main and supplementary manuscript figures

- Processed data tables required for analysis reproducibility

- Expected output directory structure

- Large raw datasets (Comet outputs, Ascore assignments, FASTAs) are provided as GitHub Release assets.

## ⚙️ Installation & Execution
<details> <summary><strong>Click to expand instructions</strong></summary>
### 1. <u>Download Input Files</u>

- Download the Release archive containing the folders:

  - raw_data/

  - modified_data/

  - output/

- Place all directories in the project root.

### 2. <u>Install Required R Version</u>

- This project was developed using:

  - **R 4.4.2**

  - **Bioconductor 3.22**

  _- Matching these versions ensures compatibility with the renv.lock file._

### 3. <u>Initialize Environment</u>

From the project directory, run:

`source("00_setup_environment.R")`


Restart R after this step.

### 4. <u>Restore Software Environment</u>

`renv::restore()`


- This installs all CRAN and Bioconductor packages pinned in the lockfile.

### 5. <u>Generate All Manuscript Figures</u>
`source("00_master_run_all_figures.R")`


- Figures will be saved into subdirectories of output/

- Intermediate tables appear in modified_data/

</details>
## 🌐 Interactive Dataset
<details> <summary><strong>Launch the interactive web application</strong></summary>

Explore the full R2HaPpY phosphotyrosine dataset, including phosphosite dynamics, gene-level filtering, replicate comparisons, and functional annotation:

➡ **https://r2happy.gs.washington.edu**

</details>
## 🧪 Experimental Protocols
<details> <summary><strong>View experimental workflow</strong></summary>

Full experimental protocols (SH2 Superbinder expression, bead preparation, enrichment, LC-MS/MS acquisition) are available on protocols.io:

➡ Protocols.io Link

</details>
## 📁 Data Availability
<details> <summary><strong>Download raw and processed datasets</strong></summary>

Raw proteome and phosphoproteome CSVs:

**[R2HaPpY raw data (v1.0.1 )](https://github.com/alexis-cmyk/ManuscriptCodeShare_R2HaPpY/releases/tag/v1.0.1-InputData_and_OutputFolders)**

- The zip archive includes:

  - Search-engine outputs (Comet)
  
  - Ascore assignments
  
  - Peptide-level redundancy mappings
  
  - Required metadata
  
  - Directory structure for figure and table outputs
  
**Unzip and place the three folders (raw_data/ , output/ , and modified_data/) into this project directory.

</details>
## 📫 Contact

Alexis Chang
Villen Laboratory, Genome Sciences
University of Washington, Seattle

Lab website: https://villenlab.gs.washington.edu/wordpress/