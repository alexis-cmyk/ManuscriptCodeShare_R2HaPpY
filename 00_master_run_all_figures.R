### This RScript points to each figure-specific script to run sequentially with one line ### -----------

###  To run this code, in console within this Rproject, run: source("00_master_run_all_figures.R")

library(rmarkdown)

rmarkdown::render("01_Script_MainAndSuppFig1_NHS_vs_Halo_1mgInput_Triplicate.Rmd")

rmarkdown::render("02_Script_MainFig2_VaryPeptideConcentration_HighpTyrContent")

rmarkdown::render("03_Script_MainFig3to6_SuppFig9to15_ProteomeLevel")

rmarkdown::render("04_Script_MainFig3to6_SuppFig9to15_GlobalPhosphoSiteLevel")

rmarkdown::render("05_Script_MainFig3to6_SuppFig9to15_GlobalPhosphoIsoformLevel")

rmarkdown::render("06_Script_MainFig3to6_SuppFig9to15_pTyrSiteLevel")

rmarkdown::render("07_Script_MainFig3to6_SuppFig9to15_pTyrIsoformLevel")

rmarkdown::render("08_Script_MainFig3to6_SuppFig9to15_Reddy2016_EGFRTMT_sectimecourse")

rmarkdown::render("09_Script_SuppFig2_NHS_vs_Halo_VaryPeptideInput")

rmarkdown::render("10_Script_SuppFig5_VaryBeadAmount")

rmarkdown::render("11_Script_SuppFig6_VaryPeptideConcentration_LowpTyrContent")

rmarkdown::render("12_Script_SuppFig7_EasyPositiveControlSample_Yeast_vSrc")

rmarkdown::render("13_Script_SuppFig8_BeadAge_4months")

