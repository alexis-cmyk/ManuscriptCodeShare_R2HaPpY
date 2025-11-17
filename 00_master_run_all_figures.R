### This RScript points to each figure-specific script to run sequentially with one line ### -----------
###  To run this code, in console within this Rproject, run: source("00_master_run_all_figures.R")
library(dplyr)

source(knitr::purl("01_Script_MainAndSuppFig1_NHS_vs_Halo_1mgInput_Triplicate.Rmd"))
source(knitr::purl("02_Script_MainFig2_VaryPeptideConcentration_HighpTyrContent.Rmd"))
source(knitr::purl("03_Script_MainFig3to6_SuppFig9to15_ProteomeLevel.Rmd"))
source(knitr::purl("04_Script_MainFig3to6_SuppFig9to15_GlobalPhosphoSiteLevel.Rmd"))
source(knitr::purl("05_Script_MainFig3to6_SuppFig9to15_GlobalPhosphoIsoformLevel.Rmd"))
source(knitr::purl("06_Script_MainFig3to6_SuppFig9to15_pTyrSiteLevel.Rmd"))
source(knitr::purl("07_Script_MainFig3to6_SuppFig9to15_pTyrIsoformLevel.Rmd"))
source(knitr::purl("08_Script_MainFig3to6_SuppFig9to15_Reddy2016_EGFRTMT_sectimecourse.Rmd"))
source(knitr::purl("09_Script_SuppFig2_NHS_vs_Halo_VaryPeptideInput.Rmd"))
source(knitr::purl("10_Script_SuppFig5_VaryBeadAmount.Rmd"))
source(knitr::purl("11_Script_SuppFig6_VaryPeptideConcentration_LowpTyrContent.Rmd"))
source(knitr::purl("12_Script_SuppFig7_EasyPositiveControlSample_Yeast_vSrc.Rmd"))
source(knitr::purl("13_Script_SuppFig8_BeadAge_4months.Rmd"))





