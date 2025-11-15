## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(janitor)
library(readr)
library(gridExtra)
library(ggrepel)
library(cowplot)
library(reshape2)
library(eulerr)
library(protti)
library(iq)
library(RColorBrewer)
#for PCA loadings
library(FactoMineR)

#for protti sample correlation and volcano plots
library(dendextend)
library(pheatmap)
library(seriation)
library(UpSetR)

library(gprofiler2)

#for c-means clustering
library(e1071)


## --------------------------------------------------------------------------------------------------------------------------------------------------
alexis_theme <- function() {
  theme(
    # panel.border = element_rect(colour = "blue", fill = NA, linetype = 2),
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.x  = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black", face = "plain", family = "sans"),
    axis.title = element_text(colour = "black", family = "sans", size = 10),
    axis.ticks = element_line(colour = "black"),
    title = element_text(size = 8, hjust = 0.5),
    # legend at the bottom 6)
    legend.position = "right")   
}


## --------------------------------------------------------------------------------------------------------------------------------------------------
my_hier_cluster <- function (data,
                             sample,
                             grouping,
                             intensity_log2,
                             condition_order_for_colors = c("EGF0min", "EGF1min", "EGF3min", "EGF5min", "EGF15min"),
                             breaklist = seq(0.5, 1, by = 0.1),
                             # legend_breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1),
                             # legend_labels = c(0.5, 0.6, 0.7, 0.8, 0.9, 1),
                             
                             condition,
                             digestion = NULL, 
                             run_order = NULL, 
                             method = "spearman", 
                             fontsize_number = 14,
                             number_color = "white",
                             cell_height = 20, 
                             plot_title = "Pearson correlation of phosphosite intensities",
                             
                             interactive = FALSE) 
{
  
  protti_colours <- "placeholder"
  
  utils::data("protti_colours", envir = environment())
  
  viridis_colours <- "placeholder"
  
  utils::data("viridis_colours", envir = environment())
  
  
  
  
  correlation <- data %>%
    dplyr::distinct({{sample}}, {{grouping}}, {{intensity_log2}}) %>%
    tidyr::pivot_wider(names_from = { {sample}},
                       values_from = {{intensity_log2}}) %>%
    tibble::column_to_rownames(var = rlang::as_name(enquo(grouping))) %>% 
    stats::cor(method = {{method}}, use = "complete.obs")
  
  
  
  
  annotation <- data %>%
    dplyr::mutate(`:=`({{condition}}, as.character({{condition}}))) %>%
    dplyr::distinct({{sample}}, {{condition}}, {{digestion}}, {{run_order}}) %>%
    tibble::column_to_rownames(var = rlang::as_name(enquo(sample)))
  
  n_conditions <- 0
  n_digest <- 0
  n_run_ord <- 0
  conditions_colour <- c()
  digest_colours <- c()
  run_ord_colours <- c()
  
  
  
  
  if (!missing(condition)) {
    # conditions <- unique(dplyr::pull(annotation, {{condition}}))
    # n_conditions <- length(conditions)
    conditions <- {{condition_order_for_colors}}
    n_conditions <- length(conditions)
    conditions_colours <- protti_colours[1:n_conditions]
    names(conditions_colours) <- conditions
  }
  
  
  
  
  
  
  if (!missing(digestion)) {
    digest <- unique(dplyr::pull(annotation, {{digestion}}))
    n_digest <- length(digest)
    digest_colours <- protti_colours[(n_conditions + 1):(n_digest + 
                                                           n_conditions)]
    names(digest_colours) <- digest
  }
  
  
  
  
  
  if (!missing(run_order)) {
    colfunc <- grDevices::colorRampPalette(c("#0D0887", 
                                             "#2E0595", "#46039F", "#5C01A6", "#7201A8", "#8707A6", 
                                             "#9A169F", "#AC2694", "#BC3587", "#CA457A", "#D6556D", 
                                             "#E26561", "#EB7655", "#F48849", "#FA9B3D", "#FDAF31", 
                                             "#FDC527", "#F9DC24", "#F0F921"))
    run_ord <- unique(dplyr::pull(annotation, {
      {
        run_order
      }
    }))
    n_run_ord <- length(run_ord)
    run_ord_colours <- colfunc(n_run_ord)
    names(run_ord_colours) <- run_ord
  }
  
  
  
  
  annotation_colours <- list(conditions_colours, digest_colours, run_ord_colours)
  
  
  names(annotation_colours) <- c(if (!missing(condition)) {
    rlang::as_name(enquo(condition))
  } else {
    "condition"
  }, if (!missing(digestion)) {
    rlang::as_name(enquo(digestion))
  } else {
    "digestion"
  }, if (!missing(run_order)) {
    rlang::as_name(enquo(run_order))
  } else {
    "run_order"
  })
  
  
  
  
  
  #interactive -----------------------------------------------------------------------------------
  
  if (interactive == TRUE) {
    if (!requireNamespace("heatmaply", quietly = TRUE)) {
      message("Package \"heatmaply\" is needed for this function to work. Please install it.", 
              call. = FALSE)
      return(invisible(NULL))
    }
    heatmap_interactive <- heatmaply::heatmaply(correlation, 
                                                main = "Correlation based hirachical clustering of samples", 
                                                col_side_colors = annotation, col_side_palette = c(annotation_colours[[1]], 
                                                                                                   annotation_colours[[2]], annotation_colours[[3]]),
                                                display_numbers(round(correlation, 2)),
                                                k_col = NA, k_row = NA, plot_method = "plotly")
    return(heatmap_interactive)
  }
  
  
  
  
  
  #not interactive ----------------------------------------------------------------------------------
  
  if (interactive == FALSE) {
    dependency_test <- c(dendextend = !requireNamespace("dendextend", 
                                                        quietly = TRUE),
                         pheatmap = !requireNamespace("pheatmap", 
                                                      quietly = TRUE), seriation = !requireNamespace("seriation", quietly = TRUE))
    if (any(dependency_test)) {
      
      dependency_name <- names(dependency_test[dependency_test ==  TRUE])
      
      if (length(dependency_name) == 1) {
        message("Package \"", paste(dependency_name), "\" is needed for this function to work. Please install it.", call. = FALSE)
        return(invisible(NULL))
      }
      
      else {
        message("Packages \"", paste(dependency_name, collapse = "\" and \""), "\" are needed for this function to work. Please install them.", call. = FALSE)
        return(invisible(NULL))
      }
    }
    
    
    distance <- stats::dist(correlation)
    
    hierachical_clustering <- stats::hclust(distance)
    
    dendrogram <- stats::as.dendrogram(hierachical_clustering)
    
    dendrogram_row <- dendextend::seriate_dendrogram(dendrogram, 
                                                     distance, method = "OLO")
    
    dendrogram_column <- dendextend::rotate(dendrogram_row, 
                                            order = rev(labels(distance)[seriation::get_order(stats::as.hclust(dendrogram_row))]))
    
    #actual heatmap with cell labels, not interactive -------------------------------------------------------
    heatmap_static <- pheatmap::pheatmap(correlation,
                                         cluster_rows = stats::as.hclust(dendrogram_row), 
                                         cluster_cols = stats::as.hclust(dendrogram_column), 
                                         # display_numbers = TRUE,
                                         display_numbers = round(correlation,2),  #annotations within column, can change to asterisks, see above
                                          #can change df = correlation to a variable to allow either matrix or filtered matrix input, such as for asterixes, etc.
                                         number_color = {{number_color}}, #in-cell annotations
                                         fontsize_number = {{fontsize_number}}, #fontsize of in-cell annotations
                                         annotation = annotation,
                                         annotation_colors = annotation_colours,
                                         # legend_breaks = {{legend_breaks}},
                                         # legend_labels = {{legend_labels}},
                                         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length({{breaklist}})),
                                         
                                          breaks = {{breaklist}},
                                         main = {{plot_title}}, 
                                         # color = viridis_colours,
                                         silent = TRUE)
    return(heatmap_static)
  }
}



## --------------------------------------------------------------------------------------------------------------------------------------------------
human_fasta_2024 <- read_csv("raw_data/MainFig3to6_SuppFig9to15/sp_reviewed_UP000005640_2024_03_21_forR.csv")


#just to convert PTMSEA gene names to uniprot ID references
fasta_gene_names <- read_tsv(file = "raw_data/MainFig3to6_SuppFig9to15/uniprotkb_proteome_UP000005640_2024_09_04.tsv") %>%
  clean_names() %>% 
  #many gene names per uniprot reference, hope first gene name matches,
  #otherwise figure out if secondary matches improve join.
  separate(gene_names, into = c("gene", "gene2", "gene3", "gene4", "gene5", "gene6"), remove = FALSE)

#subcellular location analysis:
fasta_subcellular_location <- read_tsv("raw_data/MainFig3to6_SuppFig9to15/uniprotkb_proteome_UP000005640_2024_10_16_locations.txt") %>% clean_names() %>% 
  rename(reference = entry) %>% 
  filter(grepl("HUMAN", sequence) == FALSE) #remove titin which was annotated oddly


## --------------------------------------------------------------------------------------------------------------------------------------------------
EGFR1_pathway <- read_csv("raw_data/MainFig3to6_SuppFig9to15/EGFR1_phosphopathway_PTMSEA.csv") %>% clean_names() %>%
  left_join(y = fasta_gene_names, by = "gene", relationship = "many-to-many") %>%  #match on gene1 worked!
  select(-c(gene_names : gene6)) %>% 
  rename(reference = entry) %>% 
  mutate(
    ref = paste(reference, psite, sep = "_"),
    PTM_SEA = TRUE,
    mod_res = str_sub(psite, start = 1L, end = 1L))



## --------------------------------------------------------------------------------------------------------------------------------------------------
ppdia_ErbB_signaling <- read_csv(file = "raw_data/MainFig3to6_SuppFig9to15/Phosphopedia_ErbB_gene_sites.csv") %>% 
  clean_names()


## --------------------------------------------------------------------------------------------------------------------------------------------------
# PSP_all <- read_csv(file = "raw_data/MainFig3to6_SuppFig9to15/Phosphorylation_site_dataset/Phosphorylation_site_dataset.csv") %>% clean_names() %>%
#   filter(organism == "human") %>% 
#   mutate(mod_rsd = str_replace_all(mod_rsd, "-p", ""),
#          mod_res = str_sub(mod_rsd, end = 1L)) %>% 
#   rename(psite = mod_rsd) %>% 
#   mutate(gene_ref = paste(gene, psite, sep = "_"))

# write_csv(x = PSP_all, file = "raw_data/MainFig3to6_SuppFig9to15/PSP_all_psites.csv")

PSP_all <- read_csv ("raw_data/MainFig3to6_SuppFig9to15/PSP_all_psites.csv")


## --------------------------------------------------------------------------------------------------------------------------------------------------
EGFR1_pathway_PTMSigDB <- read_csv("raw_data/MainFig3to6_SuppFig9to15/EGFR1_phosphopathway_PTMSEA.csv") %>% clean_names() %>%
  filter(str_sub(psite, end = 1L) %in% c("S", "T", "Y")) %>% 
  # left_join(y = human_fasta_2024 %>% separate(gene, into = c("gene", "organism_del"), sep = "_", remove = TRUE), by = "gene") %>% 
  left_join(y = fasta_gene_names, by = "gene", relationship = "many-to-many") %>%  #match on gene1 worked!

  rename(reference = entry) %>%
  left_join(y = (human_fasta_2024 %>% mutate(HuFasta2024 = "yes")), by = "reference") %>%
  filter(HuFasta2024 == "yes") %>% 
  select(-c(gene_names : gene6)) %>%
  
  mutate(
    mod_residue = as.numeric(str_sub(psite, start = 2L)),
    ref = paste(reference, psite, sep = "_"),
    PTM_SEA = TRUE,
    mod_res = str_sub(psite, start = 1L, end = 1L)) %>% 
  
  #filter by mod_residue extracted from protein sequence equaling Y, then we know it is more likely the correct sequence.
  mutate(
    extracted_mod_res = str_sub(sequence, start = mod_residue, end = mod_residue) ) %>% 
  filter(extracted_mod_res == mod_res) %>%  #only keep entries where p-site index matches the expected psite residue (3362 to 914 rows, better)
  # filter(mod_res == "Y") %>%  #only considering pY sites from the PTM-SEA EGFR signature here.
  mutate(n_unique_pSTY_sites = n_distinct(gene.x, psite)) %>% 
  
  #still have some duplicated entries per p-site
  #now join to human_fasta_2024 to only keep references that could potentially be in my search database if duplicated.
  # rename(full_protein_sequence = sequence) %>% 
  rename(gene = gene.x) %>% 
  select(-gene.y) %>% 
  select(gene, reference, mod_res, mod_residue, PTM_SEA, sequence)

#so, these multiple references and protein sequences per gene name is present in the human fasta 2024 as well. I will then just collapse to one psite per gene after joining against my own data below.



## --------------------------------------------------------------------------------------------------------------------------------------------------
# ppdia_ErbB_signaling <- read_csv(file = "raw_data/MainFig3to6_SuppFig9to15/Phosphopedia_ErbB_gene_sites.csv") %>% 
#   clean_names() %>% 
#   rename(mod_residue = psite) %>%
#   mutate(PPDIA = TRUE,
#          mod_res = str_sub(mod_residue, end = 1L))


## --------------------------------------------------------------------------------------------------------------------------------------------------
PSP_EGFR <- read_csv(file = "raw_data/MainFig3to6_SuppFig9to15/PSP_EGF_all_20241219.csv") %>% 
  clean_names() %>% 
  filter(organism == "human") %>% #just human sites
  rename(reference = acc_id) %>% 
  rename(psite = mod_rsd) %>%
  filter(grepl("-", reference) ==FALSE) %>%  #remove isoforms, my data doesn't search agaisnt a fasta containing isoforms
  mutate(PSP = TRUE,
    mod_res = str_sub(psite, end = 1L),
         mod_residue = as.numeric(str_sub(psite, start = 2L))) %>% 
  filter(mod_type == "p") %>%  #also includes glyco, ubi
  # filter(mod_res %in% c("S", "T", "Y")) #only include S, T or Y residues
  # filter(mod_res == "Y") %>%  #only include Y residues
  left_join(y = human_fasta_2024 %>% distinct(reference, full_protein_sequence) %>% mutate(HuFasta2024 = "yes"), by = "reference") %>% 

  mutate(
    motif_in_full_protein_seq = case_when(
      str_detect(full_protein_sequence, str_to_upper(str_remove_all(site_7_aa, "_"))) == TRUE ~ "match",
      str_detect(full_protein_sequence, str_to_upper(str_remove_all(site_7_aa, "_"))) == FALSE ~ "no",
      TRUE ~ "what happened here")) %>% #all matches
  select(gene, reference, PSP, full_protein_sequence,motif_in_full_protein_seq, everything()) %>% 
  select(gene, reference, mod_res, mod_residue, PSP, full_protein_sequence)


## --------------------------------------------------------------------------------------------------------------------------------------------------
PSP_vs_PTMSigDB_EGF <- EGFR1_pathway_PTMSigDB %>%
  full_join(y = PSP_EGFR, by = c("gene", "mod_res", "mod_residue")) %>% 
  mutate(
    overlap = case_when(
      PSP == TRUE & PTM_SEA == TRUE ~ "match",
      TRUE ~ "not match"))

PSP_vs_PTMSigDB_EGF_matches <-    PSP_vs_PTMSigDB_EGF %>% 
  mutate(
    gene_reference_match = case_when(
      overlap == "match" & reference.x == reference.y ~ "full match",
      overlap == "match" & reference.x != reference.y ~ "soft match", #gene but not reference matches
      TRUE ~ "not matched"), 
    DB_overlap = case_when(
      PTM_SEA == TRUE & PSP == TRUE ~ "both\ndatabases",
      PSP == TRUE & is.na(PTM_SEA) ~ "PSP\nonly",
      is.na(PSP) & PTM_SEA == TRUE ~ "PTMSigDB\nonly")) %>%  
    # collape multiples of same site due to different genes and references in databases
  distinct(gene, mod_res, mod_residue, DB_overlap, PSP, PTM_SEA)

#COUNTS-------------------------------------
sums_PSP_vs_PTMSigDB_EGF_matches_plotting <- PSP_vs_PTMSigDB_EGF_matches %>%
  
  distinct(gene, mod_residue, mod_res, DB_overlap) %>% 
  group_by(DB_overlap) %>%
  mutate(n_sites = n()) %>%
  ungroup() %>% 
  distinct(DB_overlap, n_sites) %>% 
  mutate(n_sites_char = as.character(n_sites),
         overlap_num = case_when(
           DB_overlap == "both\ndatabases" ~ 1,
           DB_overlap == "PSP\nonly" ~ 2,
           DB_overlap == "PTMSigDB\nonly" ~ 3))

#PLOT -------------------------------------
plot_overlaps_PSP_vs_PTMSigDB_EGF_matches <- ggplot() +
  geom_bar(data = PSP_vs_PTMSigDB_EGF_matches,
           mapping = aes(x= DB_overlap, fill = mod_res), color = "black", size = 0.75, show.legend = FALSE) +
  geom_text(data = sums_PSP_vs_PTMSigDB_EGF_matches_plotting, aes(x = overlap_num, y = n_sites +100, label = n_sites_char), size = 5) +
  scale_fill_viridis_d(direction = 1) +
  alexis_theme() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, lineheight = 0.75, vjust = 1),
        axis.text.y = element_text( size = 12),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  scale_y_continuous(limits = c(0, 3000),breaks = c(seq(0, 3000, 500)), expand = c(0, 0)) +
  ylab("Annotated pSTY sites") +
  xlab("database overlap")

plot_overlaps_PSP_vs_PTMSigDB_EGF_matches

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/overlaps_PSP_vs_PTMSigDB_EGF_matches.png", plot = plot_overlaps_PSP_vs_PTMSigDB_EGF_matches, width = 6, height = 8, scale = 0.4)
ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/overlaps_PSP_vs_PTMSigDB_EGF_matches.pdf", plot = plot_overlaps_PSP_vs_PTMSigDB_EGF_matches, width = 6, height = 8, scale = 0.4)




## --------------------------------------------------------------------------------------------------------------------------------------------------
WP437_EGFR_node_depths <- read_csv(file = "raw_data/MainFig3to6_SuppFig9to15/EGFR_WP437_node_depths.csv") %>%  clean_names()


## --------------------------------------------------------------------------------------------------------------------------------------------------

comet <- read_csv("raw_data/MainFig3to6_SuppFig9to15/GlobalPhospho/comet/result.csv") %>%
  clean_names() %>% 
  separate(col = sample_name, into = c("experiment", "enrichment", "condition", "replicate", "R2P2_well", "MS_method", "inj_vol"), sep = "_", remove = FALSE) %>%
  filter(grepl("tc", condition) == FALSE) %>% 
  mutate(sample_id = paste(condition, replicate, sep = " ")) %>% 
  left_join(y = human_fasta_2024, by = "reference") %>% 
  rename(intensity = max_intensity_light_c2837) %>% 
  filter(reverse == FALSE) %>% 
  mutate(phospho = case_when(
    grepl("\\@", sequence) ~ "phospho",
    TRUE ~ "none"))


## --------------------------------------------------------------------------------------------------------------------------------------------------
ascore <- read_csv("raw_data/MainFig3to6_SuppFig9to15/GlobalPhospho/ascore/result.csv") %>%
  clean_names() %>% 
  separate(col = sample_name, into = c("experiment", "enrichment", "condition", "replicate", "R2P2_well", "MS_method", "inj_vol"), sep = "_", remove = FALSE) %>% 
  filter(grepl("tc", condition) == FALSE) %>%
  mutate(sample_id = paste(condition, replicate, sep = " ")) %>% 
  left_join(y = human_fasta_2024, by = "reference") %>% 
  rename(intensity = max_intensity_light_c2837) %>% 
  filter(reverse == FALSE) %>% 
  mutate(phospho = case_when(
    grepl("\\@", sequence) ~ "phospho",
    TRUE ~ "none"))


## --------------------------------------------------------------------------------------------------------------------------------------------------
ascore.site2 <- ascore %>% 
  distinct() %>%
  
  mutate(seq_no_mods = str_remove_all(ascore_sequence, "[:punct:]"),
         seq_no_mods = str_remove_all(seq_no_mods, "[n]")) %>% ##remove n from n## mark of nAc modification. 
  mutate_at(c("a_score_1", "a_score_2", "a_score_3", "position_1", "position_2", "position_3"), funs(as.numeric)) %>% 
  mutate_at(c("a_score_1", "a_score_2", "a_score_3", "position_1", "position_2", "position_3"), ~replace_na(., 0)) %>%

  select(sample_name, condition, replicate, sample_id,  reference, gene,
         ascore_sequence, num_sites, a_score_1, a_score_2, a_score_3, position_1, position_2, position_3,
         seq_no_mods, intensity, q_score_c2837, num_scans_light_c2837, charge, x_corr, proteins,
         missed_cleavages, num_sites, redundancy, full_protein_sequence, organism, protein_names, length) %>% 
 
##extract modified residue
  mutate(
    mod_res1 = str_sub(seq_no_mods, start = position_1+1L, end = position_1 + 1L),
    mod_res2 = str_sub(seq_no_mods, start = position_2+1L, end = position_2 + 1L),
    mod_res3 = str_sub(seq_no_mods, start = position_3+1L, end = position_3 + 1L))


## pivot wider to match each p-site to an ascore and filter for STY
ascore.site2_longer <- ascore.site2 %>% 
  unite(ascore_mod_res1, c(a_score_1, mod_res1, position_1), sep = "_") %>% 
  unite(ascore_mod_res2, c(a_score_2, mod_res2, position_2), sep = "_") %>%
  unite(ascore_mod_res3, c(a_score_3, mod_res3, position_3), sep = "_") %>%
  pivot_longer(cols = c("ascore_mod_res1", "ascore_mod_res2", "ascore_mod_res3")) %>% 
  separate(value, into = c("ascore", "mod_res", "mod_position"), sep = "_") %>% 
  mutate(ascore= as.numeric(ascore))
  



## --------------------------------------------------------------------------------------------------------------------------------------------------
ascore.site2_longer_stringent <- ascore.site2_longer %>% 
  filter(ascore >= 13) %>% 
  filter(grepl("[STY]", mod_res) == TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
ascore_fasta <- ascore.site2_longer_stringent %>% 
  filter(!is.na(ascore_sequence))


## --------------------------------------------------------------------------------------------------------------------------------------------------
ascore_fasta_protein_mod_loc <- ascore_fasta %>% 
  mutate(
    mod_position = as.numeric(mod_position), ##mod_position parsed from ascore output above
    z_peptide = str_replace_all(seq_no_mods, "[IL]", "Z"),
    z_protein = str_replace_all(full_protein_sequence, "[IL]", "Z"),
    z_pept_position_in_z_protein = str_locate(z_protein, pattern = z_peptide), ##get peptide position from z-substituted peptide and protein
    mod_position_in_protein = z_pept_position_in_z_protein[,"start"] + mod_position, ##extract mod position from true protein sequence
    test_mod_position = str_sub(full_protein_sequence, start = mod_position_in_protein, end = mod_position_in_protein)) %>% 
  select(test_mod_position, mod_res, mod_position, ascore_sequence, everything()) %>% 
  unite(col = "mod_protein_location", c(mod_res, mod_position_in_protein), sep = "", remove = FALSE) %>% 
  unite(col = "ref", c(reference, mod_protein_location), sep = "_", remove = FALSE) 
  
  
  #join to BCA quant for correlation analysis with plate protein yield
  # left_join(y = df_BCA, by = c("condition", "replicate", "density_numeric", "density")) 

#--------------------------------------------------
    ##tally up the number of rows that have NA vs. STY in test_mod_position
mod_pos_group <- ascore_fasta_protein_mod_loc %>% 
  select(test_mod_position) %>% 
  group_by(test_mod_position) %>% 
   summarize(
    num_phospeptides = n())
mod_pos_group

    ##all mods are on S, T or Y as hoped
  
#--------------------------------------------------


## --------------------------------------------------------------------------------------------------------------------------------------------------
ascore_stringent_fasta_precursor <- ascore_fasta_protein_mod_loc %>% 
  group_by(condition, replicate, ascore_sequence, charge, ascore) %>% 
  filter(intensity == max(intensity)) %>% 
  ungroup()


## --------------------------------------------------------------------------------------------------------------------------------------------------
#replicate
write_csv(x = ascore_stringent_fasta_precursor %>% distinct(condition, replicate, ascore_sequence, charge, ascore, intensity, reference, sample_name, gene, protein_names, redundancy, proteins, mod_position_in_protein, ref), file = "modified_data/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/distinct_pSTYpeptides_ascore13_replicate.csv", col_names = TRUE)


#condition
write_csv(x = ascore_stringent_fasta_precursor %>% distinct(condition,  ascore_sequence, charge, ascore, intensity, reference, gene, protein_names,redundancy, proteins, mod_position_in_protein, ref), file = "modified_data/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/distinct_pSTYpeptides_ascore13_condition.csv", col_names = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
#keep single PSM per precursor:

ascore_fasta_protein_mod_loc_sum_psite <- ascore_fasta_protein_mod_loc %>% 
  group_by(condition, replicate, ascore_sequence, charge, ascore) %>% 
  filter(intensity == max(intensity)) %>% 
  ungroup() %>% 




  ##SUM PRECURSOR INTENSITIES TO P-SITES + median normalize

 
  group_by(condition, replicate,  sample_id, reference, mod_res,  ref) %>% #ref = unique p-site
  
  #sum PSM intensities of confident p-sites to individual p-sites
  mutate(
    sum_intensity_precursor_to_psite = sum(intensity),
    log2_psite_qty = log2(sum_intensity_precursor_to_psite)) %>% 
  ungroup() %>% 
  
  #keep single summed intensity per p-site prior to median normalization
  distinct(condition, replicate, sample_id, reference, mod_res, ref, sum_intensity_precursor_to_psite, log2_psite_qty) %>% 
  
  
  #median normalize intensities summed to p-site
  mutate(
    global_median_intensity = median(log2_psite_qty)) %>% 
  group_by(condition, replicate) %>% 
  mutate(
    sample_median_intensity = median(log2_psite_qty)) %>% 
  ungroup() %>% 
  mutate(
    median_norm_intensity = log2_psite_qty - sample_median_intensity + global_median_intensity,
    raw_median_norm_intensity = 2^median_norm_intensity)



## --------------------------------------------------------------------------------------------------------------------------------------------------
#could call df 3repsonce, but omitting for easier integration into pre-existing code
ascore_fasta_protein_mod_loc_sum_psite <- ascore_fasta_protein_mod_loc_sum_psite %>% 
  group_by(condition, ref) %>% 
  mutate(
    n_obs_per_condition = n()  ) %>% 
  ungroup() %>% 
  group_by(ref) %>% 
  mutate(max_obs_in_any_condition = max(n_obs_per_condition)) %>% 
  ungroup() %>% 
  filter(max_obs_in_any_condition > 2) #requires 3 or more observations in one condition to be considered in all later analyses.


## --------------------------------------------------------------------------------------------------------------------------------------------------
key_sites_measured_gt3reps_min_once <- ascore_fasta_protein_mod_loc_sum_psite %>% 
  distinct(ref)


## --------------------------------------------------------------------------------------------------------------------------------------------------
ascore_fasta_protein_mod_loc <- ascore_fasta_protein_mod_loc %>%
  filter(ref %in% key_sites_measured_gt3reps_min_once$ref)


## --------------------------------------------------------------------------------------------------------------------------------------------------
enrichment_efficiency_count_df <- comet %>% 
  distinct( condition, replicate,sequence, phospho) %>% 
  group_by( condition, replicate, phospho) %>% 
  summarize(
    num_phospho = n()) %>% 
  pivot_wider(names_from = phospho, values_from = num_phospho) %>% 
  mutate(
    total_pept = none + phospho,
    individual_phosphopept_enrich_efficiency = phospho / total_pept) %>% 
  ungroup() %>% 
  group_by( condition) %>% 
  mutate(
    avg_phosphopeptide_enrichment_efficiency = mean(individual_phosphopept_enrich_efficiency)) %>% 
  ungroup() %>% 
  mutate(
    sample_id = paste(condition, replicate, sep = "_") )


##view df of phosphopeptide enrichment efficiency ratios
enrichment_efficiency_count_df


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_ppept_detections <- ggplot() +
  geom_bar(data = (enrichment_efficiency_count_df %>%
                     mutate(condition = fct_relevel(condition, "EGF0min", "EGF1min", "EGF3min", "EGF5min", "EGF15min")) %>% 
                     ungroup() %>% 
                     group_by(condition) %>% 
                     mutate(
                       mean_phospho = mean(phospho)) %>% 
                     distinct(condition, mean_phospho)),
           
           mapping = aes(x = condition, y = mean_phospho),
           stat = "identity", alpha = 1, fill = "gray80") + 
  geom_jitter(data = enrichment_efficiency_count_df,
              mapping = aes(x = condition, y = phospho, fill = replicate), shape = 1,
              width = 0.25, show.legend = FALSE) +
  # geom_hline(yintercept = 12000, linewidth = 0.5, alpha = 0.5, linetype = 2) +
  # annotate(geom = "text",x = 7.5, y = 10000, label = ">12,000\np-peptides", size = 4, lineheight = 0.75) +
  theme_bw(14) +
  theme(legend.position = "none") +
  ylab("unqiue p-peptides") +
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1)) +
  # scale_y_continuous(expand = c(0,0)) +
  expand_limits(y = c(0, 13000)) +
  alexis_theme()+
  # ylim(0, 15000) +
  scale_y_continuous(expand = c(0,0)) #removes whitespace below bars!
  ## ggtitle("Phosphopeptide enrichment efficiency")

plot_ppept_detections
ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/phospho_pept_counts_noAscore.png", plot = plot_ppept_detections, width = 10, height = 10, scale = 0.5)
ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/phospho_pept_counts_noAscore.pdf", plot = plot_ppept_detections, width = 10, height = 10, scale = 0.5)


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_ppept_enrich_efficiency <- ggplot() +
  geom_bar(data = (enrichment_efficiency_count_df %>%
                     mutate(condition = fct_relevel(condition, "EGF0min", "EGF1min", "EGF3min", "EGF5min", "EGF15min")) %>%
                     distinct(condition, avg_phosphopeptide_enrichment_efficiency)),
           mapping = aes(x = condition, y = avg_phosphopeptide_enrichment_efficiency),
           stat = "identity", alpha = 1, fill = "gray80") + 
  geom_jitter(data = enrichment_efficiency_count_df,
              mapping = aes(x = condition, y = individual_phosphopept_enrich_efficiency, fill = replicate), shape = 1,
              width = 0.25, show.legend = FALSE) +
  geom_hline(yintercept = 0.96, linewidth = 1, alpha = 0.5, linetype = 2) +
  annotate(geom = "text",x = 3, y = 0.88, label = ">96%\npurity", size = 4, lineheight = 0.75) +
  theme_bw(14) +
  theme(legend.position = "none") +
  ylab("p-/all pept (count)") +
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1)) +
  # scale_y_continuous(expand = c(0,0)) +
  expand_limits(y = c(0, 1)) +
  alexis_theme()+
  scale_y_continuous(expand = c(0,0))
  ## ggtitle("Phosphopeptide enrichment efficiency")

plot_ppept_enrich_efficiency
ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/phospho_pept_enrich_efficiency_counts.png", plot = plot_ppept_enrich_efficiency, width = 6, height = 8, scale = 0.5)
ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/phospho_pept_enrich_efficiency_counts.pdf", plot = plot_ppept_enrich_efficiency, width = 6, height = 9, scale = 0.5)


## --------------------------------------------------------------------------------------------------------------------------------------------------
enrichment_efficiency_intensity_df <- comet %>% 
  distinct(sequence,  condition, replicate, phospho, intensity) %>% 
  group_by(sequence,  condition, replicate, phospho) %>% 
  
##filter by max intensity for every peptide
  filter(intensity  == max(intensity)) %>% 
  ungroup() %>% 
  group_by(phospho, condition, replicate) %>% 
  
##summarize intensity of all peptides and just p-peptides
  summarise(all_peptides_intensity = sum(intensity)) %>% 
  pivot_wider(values_from = all_peptides_intensity, names_from = phospho) %>% 
  
##grouped mutate for taking mean of summed intensities for each replicate within each condition  
  group_by(condition) %>%   
  mutate(mean_phos = mean(phospho),
         mean_non = mean(none)) %>% 
  ungroup() %>% 
  mutate(mean_ratio = mean_phos / (mean_non + mean_phos),
         rep_ratio = phospho / (none + phospho)) %>% 
  ungroup() %>% 
  mutate(
    sample_id = paste(condition, replicate, sep = "_") )
  
  
  enrichment_efficiency_intensity_df


## --------------------------------------------------------------------------------------------------------------------------------------------------
##relative to intensity of all peptides (with and without phospho)

plot_ratio_phos_intensity <- ggplot(data = (enrichment_efficiency_intensity_df %>%
                                              mutate(condition = fct_relevel(condition, "EGF0min", "EGF1min", "EGF3min", "EGF5min", "EGF15min")))) + 
  geom_bar(aes(x = condition, y = mean_ratio),
           stat = "identity",
           position = "dodge",
           fill = "skyblue2") +
  geom_point(aes(x = condition,
                 y = rep_ratio,
                 fill = replicate),
             position = position_jitterdodge(dodge.width = 0.35), alpha = 1, size =2, shape = 1, show.legend = FALSE) +
  geom_hline(yintercept = 0.977, linewidth = 1, alpha = 0.5, linetype = 2) +
  annotate(geom = "text",x = 3, y = 0.87, label = "> 97%\npurity", size = 4) +
  alexis_theme()+
  scale_y_continuous(expand = c(0,0))+ 
  expand_limits(y = c(0, 1.005)) +
  # theme_bw(14) +
  # theme(axis.text.x = element_text(angle=90, vjust = 0.25, hjust = 1), legend.position = "none") +
  # theme(legend.position = "none") +
  ylab("p-/all pept (intensity)") + xlab("conditions")

plot_ratio_phos_intensity
ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/phospho_pept_intensity_ratio.png", plot = plot_ratio_phos_intensity, width = 6, height = 8, scale = 0.4)
ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/phospho_pept_intensity_ratio.pdf", plot = plot_ratio_phos_intensity, width = 6, height = 9, scale = 0.5)


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_psm_intensity_distribution_reps <- ggplot(data = comet %>%
                                                 mutate(condition = fct_relevel(condition, "EGF0min", "EGF1min", "EGF3min", "EGF5min", "EGF15min"))) +
  geom_boxplot(mapping = aes(x = condition, y = log2(intensity), fill = as.factor(replicate)),
               outlier.shape = 21,
               outlier.alpha = 0.3) +
  theme_bw(18) +
  theme(legend.position = "none") +

  scale_fill_brewer(palette = "PuBu") +
  theme(axis.text.x = element_text(angle=90, vjust = 0.25, hjust = 1), legend.position = "none") +
  ylab(expression(log[2]~(intensity))) + xlab("lysate and elution type") +
  ggtitle("Intensity of all PSMs")

plot_psm_intensity_distribution_reps


ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/all_psm_intensity_distribution_reps.png", plot = plot_psm_intensity_distribution_reps, width = 16, height = 16, scale = 0.4)
ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/all_psm_intensity_distribution_reps.pdf", plot = plot_psm_intensity_distribution_reps, width = 16, height = 16, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# plot_psm_intensity_distribution_reps <- ggplot(data = comet ) +
#   geom_boxplot(mapping = aes(x = condition, y = log2(intensity), fill = as.factor(replicate)),outliers = FALSE) +
#   theme_bw(10) +
#   theme(legend.position = "none") +
# 
#   scale_fill_brewer(palette = "Purples") +
#   theme(axis.text.x = element_text(angle=90, vjust = 0.25, hjust = 1), legend.position = "none") +
#   ylab(expression(log[2]~(intensity))) + xlab("lysate and elution type") +
#   ggtitle("Intensity of all PSMs")
# 
# plot_psm_intensity_distribution_reps
# 
# 
# ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_psm_intensity_distribution_reps.png", plot = plot_psm_intensity_distribution_reps, width = 10, height = 10, scale = 0.4)
# ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_psm_intensity_distribution_reps.pdf", plot = plot_psm_intensity_distribution_reps, width = 10, height = 10, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
distinct_p_sites_condition <- ascore_fasta_protein_mod_loc %>% 
  distinct(reference, mod_protein_location, mod_res, condition)
##this collapses all  replicates into a single set of unique p sites.


## --------------------------------------------------------------------------------------------------------------------------------------------------
write_csv(x = distinct_p_sites_condition, file = "modified_data/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pSTYsites_ascore13_condition.csv", col_names = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
counts_distinct_psites_condition <- distinct_p_sites_condition %>% 
  distinct() %>% 
  group_by(condition, mod_res) %>% 
  summarize(
    n_distinct_sites = n()) %>% 
  ungroup() %>% 
  group_by(condition) %>% 
  mutate(total_sites = sum(n_distinct_sites)) %>% 
  ungroup() %>% 
  mutate(
    ratio_pY_sites = n_distinct_sites / total_sites * 100) %>% 
  group_by(mod_res) %>% 
  mutate(
    average_pYsite_percent = mean(ratio_pY_sites) ) %>% 
  ungroup()

counts_distinct_psites_condition


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_distinct_p_sites_condition <- ggplot(data = distinct_p_sites_condition %>% 
                                            mutate(condition = fct_relevel(condition, "EGF0min", "EGF1min", "EGF3min", "EGF5min", "EGF15min"))) +
  geom_bar(mapping = aes(x = condition, fill = mod_res)) +
  theme_bw(18) +
  scale_fill_viridis_d() +
  ylab("unique phospho sites") +
  alexis_theme()+
  expand_limits(y = c(0, 12000)) +
  scale_y_continuous(expand = c(0,0))
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))

plot_distinct_p_sites_condition

ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/distinct_p_sites_condition.png", plot = plot_distinct_p_sites_condition, width = 8, height = 8, scale = 0.5)

## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_distinct_p_sites_condition <- ggplot(data = distinct_p_sites_condition %>% 
                                            mutate(condition = fct_relevel(condition, "EGF0min", "EGF1min", "EGF3min", "EGF5min", "EGF15min"))) +
  geom_bar(mapping = aes(x = condition, fill = mod_res)) +
  theme_bw(18) +
  scale_fill_viridis_d() +
  ylab("unique phospho sites") +
  alexis_theme()+
  expand_limits(y = c(0, 12000)) +
  scale_y_continuous(expand = c(0,0))
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))

plot_distinct_p_sites_condition

ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/distinct_p_sites_condition.png", plot = plot_distinct_p_sites_condition, width = 8, height = 8, scale = 0.5)


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_distinct_pY_sites_condition <- ggplot(data = distinct_p_sites_condition %>% 
                                            mutate(condition = fct_relevel(condition, "EGF0min", "EGF1min", "EGF3min", "EGF5min", "EGF15min")) %>% filter(mod_res == "Y")) +
  geom_bar(mapping = aes(x = condition, fill = mod_res), color = "black", size = 1) +
  theme_bw(18) +
  scale_fill_viridis_d(direction = -1) +
  ylab("unique phospho sites") +
  alexis_theme()+
  expand_limits(y = c(0, 250)) +
  scale_y_continuous(expand = c(0,0))
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))

plot_distinct_pY_sites_condition

ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_distinct_pY_sites_condition.png", plot = plot_distinct_pY_sites_condition, width = 8, height = 8, scale = 0.5)


## --------------------------------------------------------------------------------------------------------------------------------------------------
distinct_p_sites_reps <- ascore_fasta_protein_mod_loc %>% 
  ungroup() %>% 
  # filter(replicate %in% c("01", "02", "03", "1", "2", "3")) %>% #distinct(replicate) - successfully downsampled to 3 reps 
  distinct(condition, replicate, sample_id, reference, mod_protein_location, mod_res)
##this collapses all  replicates into a single set of unique p sites.


## --------------------------------------------------------------------------------------------------------------------------------------------------
write_csv(x = distinct_p_sites_reps, file = "modified_data/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pSTYsites_ascore13_reps.csv", col_names = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_distinct_p_sites_reps <- ggplot(data = distinct_p_sites_reps %>% 
                                       mutate(sample_id = fct_relevel(sample_id,
                                                                      "EGF0min rep1","EGF0min rep2","EGF0min rep3","EGF0min rep4","EGF0min rep5","EGF0min rep6",
                                                                      "EGF1min rep1","EGF1min rep2","EGF1min rep3","EGF1min rep4","EGF1min rep5","EGF1min rep6",
                                                                      "EGF3min rep1","EGF3min rep2","EGF3min rep3","EGF3min rep4","EGF3min rep5","EGF3min rep6",
                                                                      "EGF5min rep1","EGF5min rep2","EGF5min rep3","EGF5min rep4","EGF5min rep5","EGF5min rep6",
                                                                      "EGF15min rep1","EGF15min rep2","EGF15min rep3","EGF15min rep4","EGF15min rep5","EGF15min rep6")))+
  geom_bar(mapping = aes(x = sample_id, fill = mod_res)) +
  theme_bw(10) +
  scale_fill_viridis_d() +
  ylab("unique phospho sites")  +
  alexis_theme()+
  expand_limits(y = c(0, 10000)) +
  scale_y_continuous(expand = c(0,0))
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))

plot_distinct_p_sites_reps

ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/distinct_p_sites_reps.png", plot = plot_distinct_p_sites_reps, width = 10, height = 7, scale = 0.4)
ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/distinct_p_sites_reps.pdf", plot = plot_distinct_p_sites_reps, width = 15, height = 12, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
summarize_distinct_psites_replicates <- distinct_p_sites_reps %>%
  ungroup() %>% 
  group_by(condition, replicate, sample_id, mod_res) %>% 
  summarize(
    n_psites = n()) %>% 
  ungroup()

summarize_distinct_psites_replicates



#errorbar ------------------------------------------------------------
df_geom_errorbar <- summarize_distinct_psites_replicates %>%
  group_by(condition, mod_res) %>%
  mutate(
    min_psites = min(n_psites), 
    max_psites = max(n_psites), 
    mid_psites = mean(n_psites)) %>% 
  ungroup() %>% 
  distinct(condition, mod_res, min_psites, max_psites) #could also include mid psites later if needed

df_geom_errorbar


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_distinct_p_sites_condition_reps_wPTS <- ggplot() +
  geom_bar(data = (distinct_p_sites_condition %>%
                     filter(grepl("pool", condition) == FALSE) %>% 
                     mutate(condition = fct_relevel(condition, "EGF0min", "EGF1min", "EGF3min", "EGF5min", "EGF15min"))),
                   mapping = aes(x = condition, fill = mod_res),show.legend = FALSE) +
  
  geom_errorbar(data = (df_geom_errorbar %>%
                          filter(grepl("pool", condition) == FALSE) %>% 
                          mutate(condition = fct_relevel(condition, "EGF0min", "EGF1min", "EGF3min", "EGF5min", "EGF15min"))),
                mapping = aes(x = condition, ymin = min_psites, ymax = max_psites),  show.legend = FALSE, linewidth = 0.5, width = 0.25, color = "grey90") +
  
  
  geom_jitter(data = (distinct_p_sites_reps%>%
                        filter(grepl("pool", condition) == FALSE) %>% 
                        mutate(condition = fct_relevel(condition, "EGF0min", "EGF1min", "EGF3min", "EGF5min", "EGF15min"))),
              mapping = aes(x = condition, color = sample_id, fill = mod_res),
              stat = "count", shape = 1, show.legend = FALSE, size = 1, stroke =0.5,
              width = 0.35, height = 0, alpha = 0.5) +
   
  theme_bw(10) +
  facet_wrap(facets = vars(mod_res), ncol = 3) +
  scale_fill_viridis_d() +
  scale_color_manual(values = rep(c("grey70" ),60)) +
  ylab("unique phospho sites")   +
  alexis_theme()+
  scale_y_continuous(expand = c(0,0)) 
  # theme(axis.text.x = element_text(angle = 30, vjust = 0.99, hjust = 1, size = 10))

plot_distinct_p_sites_condition_reps_wPTS

ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/distinct_p_sites_condition_wPTS.png", plot = plot_distinct_p_sites_condition_reps_wPTS, width = 6, height = 6, scale = 0.4)
ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/distinct_p_sites_condition_wPTS.pdf", plot = plot_distinct_p_sites_condition_reps_wPTS, width = 6, height = 6, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_distinct_pY_sites_condition_reps_wPTS <- ggplot() +
  geom_bar(data = (distinct_p_sites_condition %>%
                     filter(grepl("pool", condition) == FALSE) %>% 
                     mutate(condition = fct_relevel(condition, "EGF0min", "EGF1min", "EGF3min", "EGF5min", "EGF15min")) %>% filter(mod_res == "Y")),
                   mapping = aes(x = condition, fill = mod_res),show.legend = FALSE, color = "black", size = 0.5, width = 0.75) +
  
  
  
  
  geom_jitter(data = (distinct_p_sites_reps%>%
                        filter(grepl("pool", condition) == FALSE) %>% 
                        mutate(condition = fct_relevel(condition, "EGF0min", "EGF1min", "EGF3min", "EGF5min", "EGF15min")) %>% filter(mod_res == "Y")),
              mapping = aes(x = condition, color = sample_id, fill = mod_res),
              stat = "count", shape = 1, show.legend = FALSE, size = 0.75, stroke =0.5,
              width = 0.15, height = 0, alpha = 0.5) +
  
  geom_errorbar(data = (df_geom_errorbar %>%
                          filter(grepl("pool", condition) == FALSE) %>% 
                          mutate(condition = fct_relevel(condition, "EGF0min", "EGF1min", "EGF3min", "EGF5min", "EGF15min")) %>% filter(mod_res == "Y")),
                mapping = aes(x = condition, ymin = min_psites, ymax = max_psites),  show.legend = FALSE, linewidth = 0.25, width = 0.25, color = "grey30") +
   
  geom_text(data = tibble(x = c(0.98, 2, 3, 4, 5), y = c(183,184,215,234,238), label = c("183","184","215","234","238")), inherit.aes = FALSE, mapping = aes(x = x, y = y+10, label = label), size = 2.6) +
  
  
  # theme_bw(10) +
  # facet_wrap(facets = vars(mod_res), ncol = 3) +
  scale_fill_viridis_d(direction = -1) +
  scale_color_manual(values = rep(c("grey30" ),60)) +
  ylab("unique pY sites")   +
  alexis_theme()+
  scale_y_continuous(expand = c(0,0), limits = c(0, 260), breaks = c(seq(0, 250, 50)))+
  theme(axis.text.x = element_text(angle = -90, family = "sans", hjust = 0, vjust = 0.5))
  # theme(axis.text.x = element_text(angle = 30, vjust = 0.99, hjust = 1, size = 10))

plot_distinct_pY_sites_condition_reps_wPTS

ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_distinct_pY_sites_condition_reps_wPTS.png", plot = plot_distinct_pY_sites_condition_reps_wPTS, width = 4, height = 6, scale = 0.4)
ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_distinct_pY_sites_condition_reps_wPTS.pdf", plot = plot_distinct_pY_sites_condition_reps_wPTS, width = 5, height = 8, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
n_pSTYsites_per_uniprot <- distinct_p_sites_condition %>%
  # filter(mod_res == "Y") %>% 
  group_by(reference) %>% 
  summarize(
    n_pSTY = n_distinct(mod_protein_location)) %>% 
  ungroup()

write_csv(x = n_pSTYsites_per_uniprot, file = "modified_data/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/n_pSTYsites_per_uniprot.csv")


## --------------------------------------------------------------------------------------------------------------------------------------------------
psite_location <- ascore_fasta_protein_mod_loc %>%
  mutate(
    protein_length = str_length(full_protein_sequence)) %>% 
   
  mutate(
    # mod_position_in_protein = as.numeric(str_sub(mod_loc, start = 2L)),
    position_ratio = mod_position_in_protein / protein_length) %>% 
  # filter(!is.na(position_ratio)) %>% 
  distinct(condition, ref, position_ratio)


#histogram------------------------------------
    #define range
breaks <- seq(min(psite_location$position_ratio), max(psite_location$position_ratio), length.out = 100)

plot_psite_position_in_protein_histogram <- ggplot(data = psite_location)  +
  geom_histogram(mapping = aes(x = position_ratio, fill = condition),binwidth =  NULL, breaks = breaks) +
  theme_bw(18) + 
  facet_grid(rows = vars(condition))

plot_psite_position_in_protein_histogram

ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_psite_position_in_protein_histogram.png", plot = plot_psite_position_in_protein_histogram, width = 20, height = 16, scale = 0.4)
ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_psite_position_in_protein_histogram.pdf", plot = plot_psite_position_in_protein_histogram, width = 20, height = 16, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
Nterm_mod_ascore_listofpsites <- ascore_fasta_protein_mod_loc %>% 
  distinct( ref, seq_no_mods, full_protein_sequence, mod_position_in_protein) %>% 
  mutate(peptide_position = str_locate(full_protein_sequence, seq_no_mods),
         
         #turn matrix into normal dataframe columns by binding columns explicitly
         start = peptide_position[,1],
         end = peptide_position[,2]) %>% 
  select(-peptide_position) %>% 
  filter(start <= 2) #get p-sites found on peptides with known N-terminal processing info measured. 

                    #we do see a number of pY sites found at the N-termini, but there is a tryptic cut site that removed the actually N-termini from what we detected in the pY pulldown. and total proteome data has poor overlap with pY pulldown.


Nterm_mod_ascore <- ascore_fasta_protein_mod_loc %>% 
  filter(ref %in% Nterm_mod_ascore_listofpsites$ref) %>% 
  mutate(
    Nterm_mod_type = case_when(
      str_sub(ascore_sequence, start = 1L, end = 1L) == "M" ~ "not clipped\nnot acetylated",
      str_sub(ascore_sequence, start = 1L, end = 3L) == "n#M" ~ "not clipped\nyes acetylated",
      str_sub(ascore_sequence, start = 1L, end = 3L) != "n#M" & str_sub(ascore_sequence, start = 1L, end = 2L) == "n#" ~ "yes clipped\nyes acetylated",
      str_sub(ascore_sequence, start = 1L, end = 3L) != "n#M" & str_sub(ascore_sequence, start = 1L, end = 1L) != "M" ~ "yes clipped\nnot acetylated",
      TRUE ~ "not classified" )) %>% # all are classified!
  distinct( ref, Nterm_mod_type) %>% 
  group_by(Nterm_mod_type) %>% 
  mutate(
    n_unique_psites = n()) %>% 
  ungroup()

x_axis_labels = c("not clipped\nnot acetylated", "yes clipped\nnot acetylated", "not clipped\nyes acetylated", "yes clipped\nyes acetylated")

plot_nterm_modtypes_condition <-  ggplot(data = Nterm_mod_ascore) +
  geom_bar(mapping = aes(x = Nterm_mod_type), width = 0.5) + 
  geom_text(data = Nterm_mod_ascore %>% distinct(Nterm_mod_type, n_unique_psites), mapping = aes (x = Nterm_mod_type, y = n_unique_psites + 10, label = paste0(n_unique_psites)), size = 6) +
  alexis_theme() +
  theme(axis.text.x = element_text(angle = -60, hjust = 0, vjust = 0.9, size = 14),
        axis.title.y = element_text(size = 16)) +
  # scale_y_continuous(expand = c(0,0), limits = c(0, 200), breaks = c(seq(0, 200, 50))) +
  scale_x_discrete(limits = x_axis_labels,breaks = x_axis_labels, labels = x_axis_labels) +
  xlab("N-terminal processing type") +
  ylab("n unique p-sites")

plot_nterm_modtypes_condition

ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_nterm_modtypes_condition.png", plot = plot_nterm_modtypes_condition, width = 12, height = 12, scale = 0.4)
ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_nterm_modtypes_condition.pdf", plot = plot_nterm_modtypes_condition, width = 12, height = 12, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
df_overlaps_reps <- ascore_fasta_protein_mod_loc %>% 
  ungroup() %>% 
  distinct(condition, replicate, sample_id, ref) %>% 
  mutate(psite_present = 1) %>% 
  pivot_wider(names_from = sample_id, values_fill = 0, values_from = psite_present, id_cols = ref)

df_overlaps_reps_less <- df_overlaps_reps %>% 
  select(-ref) %>% 
  mutate_all(.funs = as.numeric) %>% 
  as.matrix() %>% 
  as.data.frame()

upset_plot_reps <- upset(data = df_overlaps_reps_less, nsets = 30, order.by = "freq", text.scale = 3, nintersects = 30) 
upset_plot_reps


png("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/upset_plot_reps.png", width = 1000, height = 600)
print(upset_plot_reps)
dev.off()


pdf("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/upset_plot_reps.pdf", width = 12, height = 8)
print(upset_plot_reps)
dev.off()


## --------------------------------------------------------------------------------------------------------------------------------------------------
df_overlaps_condition <- ascore_fasta_protein_mod_loc %>% 
  ungroup() %>% 
  distinct(condition,  ref) %>% 
  mutate(psite_present = 1) %>% 
  pivot_wider(names_from = condition, values_fill = 0, values_from = psite_present, id_cols = ref)

df_overlaps_condition_less <- df_overlaps_condition %>% 
  select(-ref) %>% 
  mutate_all(.funs = as.numeric) %>% 
  as.matrix() %>% 
  as.data.frame()

upset_plot_condition <- upset(data = df_overlaps_condition_less, nsets = 7, order.by = "freq", text.scale = 3, nintersects = 30) 
upset_plot_condition


png("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/upset_plot_condition.png", width = 1000, height = 600)
print(upset_plot_condition)
dev.off()


pdf("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/upset_plot_condition.pdf", width = 12, height = 8)
print(upset_plot_condition)
dev.off()


## --------------------------------------------------------------------------------------------------------------------------------------------------
df_overlaps_condition_gt3reps <- ascore_fasta_protein_mod_loc_sum_psite %>% 
  ungroup() %>% 
  group_by(condition, ref) %>%
  filter(n_distinct(replicate) >= 3) %>% 
  ungroup() %>% 
  distinct(condition, ref) %>% 
  mutate(psite_present = 1) %>% 
  pivot_wider(names_from = condition, values_fill = 0, values_from = psite_present, id_cols = ref)

df_overlaps_condition_gt3reps_less <- df_overlaps_condition_gt3reps %>% 
  select(-ref) %>%
  mutate_all(.funs = as.numeric) %>% 
  as.matrix() %>% #trick for formatting fxns within protti.
  as.data.frame()

# df_overlaps_more <- cbind(df_overlaps_less, psite_ids) %>%  as.matrix()

upset_plot_psites_gt3reps <- upset(data = df_overlaps_condition_gt3reps_less, nsets = 8, order.by = "freq", text.scale = 3, nintersects = 10)
upset_plot_psites_gt3reps


png("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/upset_plot_psites_gt3reps.png", width = 2000, height = 600)
print(upset_plot_psites_gt3reps)
dev.off()


pdf("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/upset_plot_psites_gt3reps.pdf", width = 12, height = 8)
print(upset_plot_psites_gt3reps)
dev.off()



## --------------------------------------------------------------------------------------------------------------------------------------------------
df_unique_sites_genes_condition <- df_overlaps_condition_gt3reps %>% 
  separate(col = ref, into = c("reference", "psite"), sep = "_", remove = FALSE) %>% 
  left_join(y = fasta_gene_names %>% rename(reference = entry), by = "reference") %>% 
  mutate(
    condition_specific = case_when(
      EGF0min == 1 & EGF1min == 0 & EGF3min == 0 & EGF5min == 0 & EGF15min == 0 ~ "control only",
      EGF0min == 0 & EGF1min == 1 & EGF3min == 0 & EGF5min == 0 & EGF15min == 0 ~ "EGF1min only",
      EGF0min == 0 & EGF1min == 0 & EGF3min == 1 & EGF5min == 0 & EGF15min == 0 ~ "EGF3min only",
      EGF0min == 0 & EGF1min == 0 & EGF3min == 0 & EGF5min == 1 & EGF15min == 0 ~ "EGF5min only",
      EGF0min == 0 & EGF1min == 0 & EGF3min == 0 & EGF5min == 0 & EGF15min == 1 ~ "EGF15min only",
      
      EGF0min == 1 & EGF1min == 1 & EGF3min == 1 & EGF5min == 1 & EGF15min == 1 ~ "all conditions",
      TRUE ~ "combo")) %>% 
   
  full_join(y = ppdia_ErbB_signaling %>% mutate(pathway = "ErbB2_Phosphopedia"), by = c("gene"), suffix = c(".mydata", ".ppdia"), relationship = "many-to-many") %>% 
  select(gene, psite.mydata, psite.ppdia, condition_specific, EGF0min, EGF1min, EGF3min, EGF5min, EGF15min, pathway,everything()) 




## --------------------------------------------------------------------------------------------------------------------------------------------------
pheatmap_df <- df_unique_sites_genes_condition %>%
  
  #remove duplicate observed psites if multiple reference psites present
  group_by(gene, psite.mydata) %>% 
  filter(n_distinct(psite.ppdia) == 1) %>% 
  ungroup() %>% 
 
  filter(!is.na(psite.ppdia)) %>%
  mutate(gene_psite = paste(gene, "OBS:", psite.mydata, "REF:", psite.ppdia, sep = " ")) %>%
  select(gene_psite, EGF0min:EGF15min) %>% 
  mutate_at(vars(EGF0min:EGF15min), as.numeric) %>% 
  filter(!is.na(EGF0min)) 

# assign df colums to matrix values then row names:
pheatmap_matrix <- as.matrix(pheatmap_df[,2:6])

rownames(pheatmap_matrix) <- as.matrix(pheatmap_df[,1])



## --------------------------------------------------------------------------------------------------------------------------------------------------



## --------------------------------------------------------------------------------------------------------------------------------------------------
ErbB2_pprotein_heatmap <- pheatmap(pheatmap_matrix, cluster_cols = FALSE, cluster_rows = TRUE,
                                   legend = TRUE, legend_breaks = c(0, 1),
                                   cutree_rows = 6, cutree_cols = 2,
                                   fontsize_row = 9, fontsize_col = 10, angle_col = 270, border_color = TRUE)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/ErbB2_pprotein_heatmap.png", plot = ErbB2_pprotein_heatmap, height = 40, width = 15, scale = 0.3)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# longer_ErbB_signaling <- df_unique_sites_genes_condition %>% 
#   pivot_longer(cols = c(EGF0min:EGF15min), names_to = c("condition")) %>% 
#   select(gene, psite.mydata, psite.ppdia, condition, value) %>% 
#   mutate(gene_psite = paste(gene, "OBS:", psite.mydata, "REF:", psite.ppdia, sep = " ")) %>% 
#   select(gene_psite, condition, value) %>% 
#   rename(psite_detected = value)
# 
# ErbB2_pprotein_gt3reps_condition <- ggplot(longer_ErbB_signaling, aes(x = condition, y = gene_psite, fill = psite_detected)) +
#   geom_tile(color = "black", show.legend = FALSE) +
#   scale_x_discrete(position = "top", expand = c(0,0)) +
#   scale_y_discrete(expand = c(0,0)) +
#   labs(x = "", y = "", fill = "") +
#   alexis_theme()
#   # facet_grid(~group) +
#   # coord_fixed(ratio = 0.5) +
#   # theme(
#   #         legend.position = "none",
#   #         panel.spacing = unit(0, "lines"), 
#   #         strip.background = element_blank(),
#   #         strip.placement = "outside")
# 
# ErbB2_pprotein_gt3reps_condition


## --------------------------------------------------------------------------------------------------------------------------------------------------
psite_overlaps_condition <- df_overlaps_condition %>% 
  ungroup() %>% 
  pivot_longer(cols = c( EGF0min, EGF1min, EGF3min, EGF5min, EGF15min), names_to = c("condition")) %>% 
  group_by(ref) %>% 
  mutate(
    n_detections = sum(value)) %>% 
  ungroup()  
  # filter(replicate %in% c(1, 2, 3)) %>% 
  
ascore_psite_overlaps <- ascore_fasta_protein_mod_loc_sum_psite %>% 
  inner_join(y = psite_overlaps_condition %>%  distinct(ref, n_detections), by = c("ref")) %>% 
  mutate(n_detections = as.factor(n_detections))




## --------------------------------------------------------------------------------------------------------------------------------------------------
#plot
  ## intensity by condition overlaps---------------------------------------------
intensity_distribution_psite_overlaps <- ggplot(data = ascore_psite_overlaps) +
  geom_violin(mapping = aes(x = n_detections, y = median_norm_intensity, fill = condition), draw_quantiles = c(0.25, 0.5, 0.75), show.legend = FALSE) +
  theme_bw(18) +
  theme(strip.text = element_text(size = 8)) +
  facet_wrap(facets = vars(condition), nrow = 3)
intensity_distribution_psite_overlaps
ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/intensity_distribution_psite_overlaps.png", plot = intensity_distribution_psite_overlaps, width = 16, height = 14, scale = 0.4)




## --------------------------------------------------------------------------------------------------------------------------------------------------
all_psite_CVs <- qc_cvs(data = ascore_fasta_protein_mod_loc_sum_psite,
                            grouping = ref,
                            intensity = sum_intensity_precursor_to_psite,
                            condition = condition,
                            plot = FALSE)

all_psite_CVs_plot <- qc_cvs(data = ascore_fasta_protein_mod_loc_sum_psite,
                            grouping = ref,
                            intensity = sum_intensity_precursor_to_psite,
                            condition = condition,
                            plot_style = "violin",
                            plot = TRUE)

#return
all_psite_CVs
all_psite_CVs_plot

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/all_psite_CVs_plotNotNorm.png", plot = all_psite_CVs_plot, width = 18, height = 14, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
all_psite_CVs_NORM <- qc_cvs(data = ascore_fasta_protein_mod_loc_sum_psite,
                            grouping = ref,
                            intensity = raw_median_norm_intensity,
                            condition = condition,
                            plot = FALSE)

all_psite_CVs_plot_NORM <- qc_cvs(data = ascore_fasta_protein_mod_loc_sum_psite,
                            grouping = ref,
                            intensity = raw_median_norm_intensity,
                            condition = condition,
                            plot_style = "violin",
                            plot = TRUE)

#return
all_psite_CVs_NORM
all_psite_CVs_plot_NORM

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/all_psite_CVs_plot_Norm.png", plot = all_psite_CVs_plot_NORM, width = 18, height = 14, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
#check median normalization by plotting.

##NOT normalized intensities --------------------------------------------------------
check_all_psite_intensity_before_normalization <- qc_intensity_distribution(
  data = ascore_fasta_protein_mod_loc_sum_psite,
  sample = sample_id, 
  grouping = ref,
  intensity_log2 = log2_psite_qty,
  plot_style = "boxplot")

check_all_psite_intensity_before_normalization

ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/check_all_psite_intensity_before_normalization.pdf", plot = check_all_psite_intensity_before_normalization, width = 15, height = 15, scale = 0.5)

ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/check_all_psite_intensity_before_normalization.png", plot = check_all_psite_intensity_before_normalization, width = 15, height = 15, scale = 0.5)

##YES normalized intensities --------------------------------------------------------
check_all_psite_intensity_after_normalization <- qc_intensity_distribution(
  data = ascore_fasta_protein_mod_loc_sum_psite,
  sample = sample_id, 
  grouping = ref,
  intensity_log2 = median_norm_intensity,
  plot_style = "boxplot")

check_all_psite_intensity_after_normalization

ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/check_all_psite_intensity_after_normalization.pdf", plot = check_all_psite_intensity_after_normalization, width = 15, height = 15, scale = 0.5)

ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/check_all_psite_intensity_after_normalization.png", plot = check_all_psite_intensity_after_normalization, width = 15, height = 15, scale = 0.5)



## --------------------------------------------------------------------------------------------------------------------------------------------------
all_psite_correlation <- qc_sample_correlation(
  data = ascore_fasta_protein_mod_loc_sum_psite,
  sample = sample_id, 
  grouping = ref,
  intensity_log2 = median_norm_intensity,
  condition = condition,
  interactive = FALSE)


all_psite_correlation

ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/all_psite_correlation.png", plot = all_psite_correlation, width = 20, height = 20, scale = 0.4)
ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/all_psite_correlation.pdf", plot = all_psite_correlation, width = 16, height = 12, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
all_psite_my_correlation_pearson <- my_hier_cluster(
  data = ascore_fasta_protein_mod_loc_sum_psite %>% filter(!is.infinite(median_norm_intensity)),
  # condition_order_for_colors = c("EGF0min", "EGF1min", "EGF3min", "EGF5min", "EGF15min"),
  breaklist = seq(0.5, 1, by = 0.1),
  sample = sample_id, 
  grouping = ref,
  fontsize_number = 4.5,
  number_color = "black",
  method = "pearson",
  intensity_log2 = median_norm_intensity,
  plot_title = "Pearson correlation (R) of all pSTY site intensities\nfrom global phosphopeptide enrichment",
  condition = condition,
  interactive = FALSE)


all_psite_my_correlation_pearson

ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/all_psite_my_correlation_pearson.png", plot = all_psite_my_correlation_pearson, width = 20, height = 20, scale = 0.4)
ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/all_psite_my_correlation_pearson.pdf", plot = all_psite_my_correlation_pearson, width = 16, height = 12, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
global_pYsite_my_correlation_pearson <- my_hier_cluster(
  data = (ascore_fasta_protein_mod_loc_sum_psite %>% filter(!is.infinite(median_norm_intensity)) %>% filter(mod_res == "Y")),
  # condition_order_for_colors = c("EGF0min", "EGF1min", "EGF3min", "EGF5min", "EGF15min"),
  breaklist = seq(0.5, 1, by = 0.1),
  sample = sample_id, 
  grouping = ref,
  fontsize_number = 4.5,
  number_color = "black",
  method = "pearson",
  intensity_log2 = median_norm_intensity,
  condition = condition,
  plot_title = "Pearson correlation (R) of pY site intensities\nfrom global phosphopeptide enrichment",
  interactive = FALSE)


global_pYsite_my_correlation_pearson

ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/global_pYsite_my_correlation_pearson.png", plot = global_pYsite_my_correlation_pearson, width = 20, height = 20, scale = 0.4)
ggsave("output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/global_pYsite_my_correlation_pearson.pdf", plot = global_pYsite_my_correlation_pearson, width = 16, height = 12, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
ascore_med_norm_missing <- assign_missingness(
    data = ascore_fasta_protein_mod_loc_sum_psite,
    sample = sample_id, 
    grouping = ref, 
    intensity = median_norm_intensity,
    condition = condition, 
    ref_condition = "all",
    completeness_MAR = 0.7, #require 2 of 5 reps to have to be considered condition specific signature
    completeness_MNAR = 0.2 ) #any time obs is in 1 or less samples per condition, the psite is considered biologically missing.


## --------------------------------------------------------------------------------------------------------------------------------------------------
psite_diff_abundance <- ascore_med_norm_missing %>% 
  calculate_diff_abundance(
    sample = sample_id,
    condition = condition, 
    grouping = ref, 
    intensity_log2 = median_norm_intensity,
    missingness = missingness,
    comparison = comparison,
    filter_NA_missingness = TRUE,
    method = "moderated_t-test" )


## --------------------------------------------------------------------------------------------------------------------------------------------------
psite_volcano <- volcano_plot(
  data = psite_diff_abundance,
  grouping = ref,
  log2FC = diff,
  significance = pval,
  significance_cutoff = c(0.05, "adj_pval"),
  method = "significant",
  facet_by = comparison)

psite_volcano
ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/psite_volcano.png", plot = psite_volcano, width = 30, height = 18, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
percent_psites_diff_regulated <- psite_diff_abundance %>%
  distinct(comparison, ref, diff, adj_pval) %>%
  group_by(comparison) %>%
  mutate(total_psites = n_distinct(ref)) %>%
  ungroup() %>%
  mutate(
    regulated_psite = case_when(
      adj_pval <= 0.05 & abs(diff) >= 1 ~ "regulated",
      TRUE ~ "not regulated")) %>%
  group_by(comparison, regulated_psite) %>%
  mutate(n_regulated_psites = n()) %>%
  ungroup() %>%
  distinct(comparison, total_psites, regulated_psite, n_regulated_psites) %>%
  pivot_wider(id_cols = c(comparison, total_psites), names_from = regulated_psite, values_from = n_regulated_psites) %>%
  group_by(comparison) %>%
  mutate(percent_regulated = regulated / total_psites *100)
  

percent_psites_diff_regulated


## --------------------------------------------------------------------------------------------------------------------------------------------------
percent_psite_regulated <- psite_diff_abundance %>% 
  distinct(comparison, diff, adj_pval, ref) %>% 
  group_by(comparison) %>% 
  mutate(total_psites = n_distinct(ref)) %>%
  ungroup() %>%
  mutate(
    regulated_psite = case_when(
      adj_pval <= 0.05 & abs(diff) >= 1 ~ "regulated",
      TRUE ~ "not regulated"),
    regulated_psite_dir = case_when(
      regulated_psite == "regulated" & diff >= 1 ~ "up",
      regulated_psite == "regulated" & diff <= -1 ~ "down",
      regulated_psite == "not regulated" ~ "not regulated")) %>%
  group_by(comparison, regulated_psite_dir) %>% 
  mutate(n_regulated_psites = n()) %>% 
  ungroup() %>% 
  distinct(comparison, total_psites, regulated_psite, n_regulated_psites, regulated_psite_dir) %>% 
  pivot_wider(id_cols= c(comparison, total_psites), names_from = regulated_psite_dir, values_from = n_regulated_psites) %>% 
  group_by(comparison) %>% 
  mutate(percent_upregulated = up / total_psites*100,
         percent_downregulated = down / total_psites*100)
  

percent_psite_regulated


## --------------------------------------------------------------------------------------------------------------------------------------------------

# EGF15min vs. EGF0min (aka control/ untreated) ------------------------------------------------------

my_psite_volcano_plot_EGF0min_vs_15min <- ggplot() +
  geom_vline(xintercept = c(-1, 1), color = "darkgrey", alpha = 0.5, size =1, linetype = 2)+
  geom_hline(yintercept = -log10(0.05), color = "darkgrey", alpha = 0.5, size = 1, linetype = 2) +
  geom_point(data = psite_diff_abundance %>%
               filter(comparison == "EGF0min_vs_EGF15min") %>% 
               # filter(abs(diff) < 1) %>% 
               filter(adj_pval >0.05),
             mapping = aes(x = diff, y = -log10(adj_pval)), color = "grey", alpha = 0.1, shape =19) +
  geom_point(data = psite_diff_abundance %>%
               filter(comparison == "EGF0min_vs_EGF15min") %>% 
               filter(abs(diff) < 1) %>%
               filter(adj_pval <0.05),
             mapping = aes(x = diff, y = -log10(adj_pval)), color = "grey", alpha = 0.1, shape =19) +
  geom_point(data = psite_diff_abundance %>%
               filter(comparison == "EGF0min_vs_EGF15min") %>%
               filter(diff < -1) %>% 
               filter(adj_pval <= 0.05),
             mapping = aes(x = diff, y = -log10(adj_pval)), color = "blue", alpha = 0.1, shape =19) +
  geom_point(data = psite_diff_abundance %>%
               filter(comparison == "EGF0min_vs_EGF15min") %>% 
               filter(diff > 1) %>% 
               filter(adj_pval <= 0.05),
             mapping = aes(x = diff, y = -log10(adj_pval)), color = "red", alpha = 0.1, shape =19) +
  theme_bw(18) +
  facet_wrap(facets = vars(comparison)) +
  scale_x_continuous(name = "Log2(fold change)", breaks = seq(-8, 8, 2), limits = c(-8, 8)) +
  scale_y_continuous(name = "-Log10(adj. p-value)", breaks = seq(0,20,5), limits = c(0, 20))

my_psite_volcano_plot_EGF0min_vs_15min

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/my_psite_volcano_plot_EGF0min_vs_15min.png", plot = my_psite_volcano_plot_EGF0min_vs_15min, width = 10, height = 8, scale = 0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/my_psite_volcano_plot_EGF0min_vs_15min.pdf", plot = my_psite_volcano_plot_EGF0min_vs_15min, width = 10, height = 8, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
psite_diff_abundance_gene <- psite_diff_abundance %>% 
  separate(col = ref, into = c("reference", "modres_position"), sep = "_", remove = FALSE) %>% 
  mutate(mod_res = str_sub(modres_position, start  = 1L, end = 1L),
         position = as.numeric(str_sub(modres_position, start = 2L))) %>%
  left_join(y = (fasta_gene_names %>% select(entry, gene, protein_names, keyword_id) %>% rename(reference = entry)), by = "reference") %>% 
    #flip the 0 min and 15 min comparison
  mutate(diff = case_when(
    comparison == "EGF0min_vs_EGF15min" ~ -diff,
    TRUE ~ diff)) %>% 
  mutate(comparison = case_when(
    comparison == "EGF0min_vs_EGF15min" ~ "EGF15min_vs_EGF0min", 
    TRUE ~ comparison)) 


## --------------------------------------------------------------------------------------------------------------------------------------------------
# EGF15min vs. EGF0min (aka control/ untreated) ------------------------------------------------------

my_psite_volcano_plot_EGF15min_vs_0min_STYcolor <- ggplot() +
  geom_vline(xintercept = c(-1, 1), color = "darkgrey", alpha = 0.5, size =1, linetype = 2)+
  geom_hline(yintercept = -log10(0.05), color = "darkgrey", alpha = 0.5, size = 1, linetype = 2) +
  geom_point(data = psite_diff_abundance_gene %>%
               filter(comparison == "EGF15min_vs_EGF0min") %>% 
               # filter(abs(diff) < 1) %>% 
               filter(adj_pval >0.05),
             mapping = aes(x = diff, y = -log10(adj_pval)), color = "grey", alpha = 0.3, shape =19) +
  geom_point(data = psite_diff_abundance_gene %>%
               filter(comparison == "EGF15min_vs_EGF0min") %>% 
               filter(abs(diff) < 1) %>%
               filter(adj_pval <0.05),
             mapping = aes(x = diff, y = -log10(adj_pval)), color = "grey", alpha = 0.3, shape =19) +
  geom_point(data = psite_diff_abundance_gene %>%
               filter(comparison == "EGF15min_vs_EGF0min") %>%
               filter(diff < -1) %>% 
               filter(adj_pval <= 0.05),
             mapping = aes(x = diff, y = -log10(adj_pval), color = mod_res), alpha = 0.3, shape =19) +
  geom_point(data = psite_diff_abundance_gene %>%
               filter(comparison == "EGF15min_vs_EGF0min") %>% 
               filter(diff > 1) %>% 
               filter(adj_pval <= 0.01),
             mapping = aes(x = diff, y = -log10(adj_pval), color = mod_res), alpha = 0.3, shape =19) +
  facet_wrap(facets = vars(comparison, mod_res)) +
  # geom_text_repel(data = psite_diff_abundance_gene %>%
  #              filter(comparison == "EGF15min_vs_EGF0min") %>% 
  #              filter(diff >= 1) %>% 
  #              filter(adj_pval <= 0.05),
  #              mapping = aes(x = diff, y = -log10(adj_pval), label = gene),
  #              max.overlaps = 200,
  #              segment.color = NA,
  #              # position = position_stack( hjust = 6),
  #              nudge_x = 8,
  #              force_pull   = 0, # do not pull toward data points
  #              # xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),
  #              xlim = c(0, 10), ylim = c(-0.5, 21),
  #              direction = "y",
  #              # vjust = 1,
  #              hjust = "right",
  #              size = 2) +
  # geom_text_repel(data = psite_diff_abundance_gene %>%
  #              filter(comparison == "EGF15min_vs_EGF0min") %>% 
  #              filter(diff <= -1.5) %>% 
  #              filter(adj_pval <= 0.01),
  #              mapping = aes(x = diff, y = -log10(adj_pval), label = gene),
  #              max.overlaps = 200,
  #              segment.color = NA,
  #              # position = position_stack( hjust = 6),
  #              nudge_x = -8,
  #              direction = "y",
  #              force_pull   = 0, # do not pull toward data points
  #              # xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),
  #              xlim = c(-10, 0), ylim = c(-0.5, 21),
  #              # vjust = 1,
  #              hjust = "left",
  #              size = 2) +
  alexis_theme() +
  coord_cartesian(clip = 'off') +
  
  scale_x_continuous(name = "Log2(fold change)", breaks = seq(-8, 8, 2), limits = c(-8, 8)) +
  scale_y_continuous(name = "-Log10(adj. p-value)", breaks = seq(0,20,5), limits = c(0, 20))

my_psite_volcano_plot_EGF15min_vs_0min_STYcolor

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/my_psite_volcano_plot_EGF15min_vs_0min_STYcolor.png", plot = my_psite_volcano_plot_EGF15min_vs_0min_STYcolor, width = 15, height = 10, scale = 0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/my_psite_volcano_plot_EGF15min_vs_0min_STYcolor.pdf", plot = my_psite_volcano_plot_EGF15min_vs_0min_STYcolor, width = 15, height = 10, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# EGF5min vs. EGF0min (aka control/ untreated) ------------------------------------------------------

my_psite_volcano_plot_EGF5min_vs_0min_STYcolor <- ggplot() +
  geom_vline(xintercept = c(-1, 1), color = "darkgrey", alpha = 0.5, size =1, linetype = 2)+
  geom_hline(yintercept = -log10(0.05), color = "darkgrey", alpha = 0.5, size = 1, linetype = 2) +
  geom_point(data = psite_diff_abundance_gene %>%
               filter(comparison == "EGF5min_vs_EGF0min") %>% 
               # filter(abs(diff) < 1) %>% 
               filter(adj_pval >0.05),
             mapping = aes(x = diff, y = -log10(adj_pval)), color = "grey", alpha = 0.3, shape =19) +
  geom_point(data = psite_diff_abundance_gene %>%
               filter(comparison == "EGF5min_vs_EGF0min") %>% 
               filter(abs(diff) < 1) %>%
               filter(adj_pval <0.05),
             mapping = aes(x = diff, y = -log10(adj_pval)), color = "grey", alpha = 0.3, shape =19) +
  geom_point(data = psite_diff_abundance_gene %>%
               filter(comparison == "EGF5min_vs_EGF0min") %>%
               filter(diff < -1) %>% 
               filter(adj_pval <= 0.05),
             mapping = aes(x = diff, y = -log10(adj_pval), color = mod_res), alpha = 0.3, shape =19) +
  geom_point(data = psite_diff_abundance_gene %>%
               filter(comparison == "EGF5min_vs_EGF0min") %>% 
               filter(diff > 1) %>% 
               filter(adj_pval <= 0.01),
             mapping = aes(x = diff, y = -log10(adj_pval), color = mod_res), alpha = 0.3, shape =19) +
  facet_wrap(facets = vars(comparison, mod_res)) +
  # geom_text_repel(data = psite_diff_abundance_gene %>%
  #              filter(comparison == "EGF5min_vs_EGF0min") %>%
  #              filter(diff >= 1) %>%
  #              filter(adj_pval <= 0.05),
  #              mapping = aes(x = diff, y = -log10(adj_pval), label = gene),
  #              max.overlaps = 200,
  #              segment.color = NA,
  #              # position = position_stack( hjust = 6),
  #              nudge_x = 8,
  #              force_pull   = 0, # do not pull toward data points
  #              # xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),
  #              xlim = c(0, 10), ylim = c(-0.5, 21),
  #              direction = "y",
  #              # vjust = 1,
  #              hjust = "right",
  #              size = 2) +
  # geom_text_repel(data = psite_diff_abundance_gene %>%
  #              filter(comparison == "EGF5min_vs_EGF0min") %>%
  #              filter(diff <= -1.5) %>%
  #              filter(adj_pval <= 0.01),
  #              mapping = aes(x = diff, y = -log10(adj_pval), label = gene),
  #              max.overlaps = 200,
  #              segment.color = NA,
  #              # position = position_stack( hjust = 6),
  #              nudge_x = -8,
  #              direction = "y",
  #              force_pull   = 0, # do not pull toward data points
  #              # xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),
  #              xlim = c(-10, 0), ylim = c(-0.5, 21),
  #              # vjust = 1,
  #              hjust = "left",
  #              size = 2) +
  alexis_theme() +
  coord_cartesian(clip = 'off') +

  scale_x_continuous(name = "Log2(fold change)", breaks = seq(-8, 8, 2), limits = c(-8, 8)) +
  scale_y_continuous(name = "-Log10(adj. p-value)", breaks = seq(0,20,5), limits = c(0, 20))

my_psite_volcano_plot_EGF5min_vs_0min_STYcolor

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/my_psite_volcano_plot_EGF5min_vs_0min_STYcolor.png", plot = my_psite_volcano_plot_EGF5min_vs_0min_STYcolor, width = 15, height = 10, scale = 0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/my_psite_volcano_plot_EGF5min_vs_0min_STYcolor.pdf", plot = my_psite_volcano_plot_EGF5min_vs_0min_STYcolor, width = 15, height = 10, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# EGF3min vs. EGF0min (aka control/ untreated) ------------------------------------------------------

my_psite_volcano_plot_EGF3min_vs_0min_STYcolor <- ggplot() +
  geom_vline(xintercept = c(-1, 1), color = "darkgrey", alpha = 0.5, size =1, linetype = 2)+
  geom_hline(yintercept = -log10(0.05), color = "darkgrey", alpha = 0.5, size = 1, linetype = 2) +
  geom_point(data = psite_diff_abundance_gene %>%
               filter(comparison == "EGF3min_vs_EGF0min") %>% 
               # filter(abs(diff) < 1) %>% 
               filter(adj_pval >0.05),
             mapping = aes(x = diff, y = -log10(adj_pval)), color = "grey", alpha = 0.3, shape =19) +
  geom_point(data = psite_diff_abundance_gene %>%
               filter(comparison == "EGF3min_vs_EGF0min") %>% 
               filter(abs(diff) < 1) %>%
               filter(adj_pval <0.05),
             mapping = aes(x = diff, y = -log10(adj_pval)), color = "grey", alpha = 0.3, shape =19) +
  geom_point(data = psite_diff_abundance_gene %>%
               filter(comparison == "EGF3min_vs_EGF0min") %>%
               filter(diff < -1) %>% 
               filter(adj_pval <= 0.05),
             mapping = aes(x = diff, y = -log10(adj_pval), color = mod_res), alpha = 0.3, shape =19) +
  geom_point(data = psite_diff_abundance_gene %>%
               filter(comparison == "EGF3min_vs_EGF0min") %>% 
               filter(diff > 1) %>% 
               filter(adj_pval <= 0.01),
             mapping = aes(x = diff, y = -log10(adj_pval), color = mod_res), alpha = 0.3, shape =19) +
  facet_wrap(facets = vars(comparison, mod_res)) +
  # geom_text_repel(data = psite_diff_abundance_gene %>%
  #              filter(comparison == "EGF3min_vs_EGF0min") %>%
  #              filter(diff >= 1) %>%
  #              filter(adj_pval <= 0.05),
  #              mapping = aes(x = diff, y = -log10(adj_pval), label = gene),
  #              max.overlaps = 200,
  #              segment.color = NA,
  #              # position = position_stack( hjust = 6),
  #              nudge_x = 8,
  #              force_pull   = 0, # do not pull toward data points
  #              # xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),
  #              xlim = c(0, 10), ylim = c(-0.5, 21),
  #              direction = "y",
  #              # vjust = 1,
  #              hjust = "right",
  #              size = 2) +
  # geom_text_repel(data = psite_diff_abundance_gene %>%
  #              filter(comparison == "EGF3min_vs_EGF0min") %>%
  #              filter(diff <= -1.5) %>%
  #              filter(adj_pval <= 0.01),
  #              mapping = aes(x = diff, y = -log10(adj_pval), label = gene),
  #              max.overlaps = 200,
  #              segment.color = NA,
  #              # position = position_stack( hjust = 6),
  #              nudge_x = -8,
  #              direction = "y",
  #              force_pull   = 0, # do not pull toward data points
  #              # xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),
  #              xlim = c(-10, 0), ylim = c(-0.5, 21),
  #              # vjust = 1,
  #              hjust = "left",
  #              size = 2) +
  alexis_theme() +
  coord_cartesian(clip = 'off') +

  scale_x_continuous(name = "Log2(fold change)", breaks = seq(-8, 8, 2), limits = c(-8, 8)) +
  scale_y_continuous(name = "-Log10(adj. p-value)", breaks = seq(0,20,5), limits = c(0, 20))

my_psite_volcano_plot_EGF3min_vs_0min_STYcolor

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/my_psite_volcano_plot_EGF3min_vs_0min_STYcolor.png", plot = my_psite_volcano_plot_EGF3min_vs_0min_STYcolor, width = 15, height = 10, scale = 0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/my_psite_volcano_plot_EGF3min_vs_0min_STYcolor.pdf", plot = my_psite_volcano_plot_EGF3min_vs_0min_STYcolor, width = 15, height = 10, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# EGF1min vs. EGF0min (aka control/ untreated) ------------------------------------------------------

my_psite_volcano_plot_EGF1min_vs_0min_STYcolor <- ggplot() +
  geom_vline(xintercept = c(-1, 1), color = "darkgrey", alpha = 0.5, size =1, linetype = 2)+
  geom_hline(yintercept = -log10(0.05), color = "darkgrey", alpha = 0.5, size = 1, linetype = 2) +
  geom_point(data = psite_diff_abundance_gene %>%
               filter(comparison == "EGF1min_vs_EGF0min") %>% 
               # filter(abs(diff) < 1) %>% 
               filter(adj_pval >0.05),
             mapping = aes(x = diff, y = -log10(adj_pval)), color = "grey", alpha = 0.3, shape =19) +
  geom_point(data = psite_diff_abundance_gene %>%
               filter(comparison == "EGF1min_vs_EGF0min") %>% 
               filter(abs(diff) < 1) %>%
               filter(adj_pval <0.05),
             mapping = aes(x = diff, y = -log10(adj_pval)), color = "grey", alpha = 0.3, shape =19) +
  geom_point(data = psite_diff_abundance_gene %>%
               filter(comparison == "EGF1min_vs_EGF0min") %>%
               filter(diff < -1) %>% 
               filter(adj_pval <= 0.05),
             mapping = aes(x = diff, y = -log10(adj_pval), color = mod_res), alpha = 0.3, shape =19) +
  geom_point(data = psite_diff_abundance_gene %>%
               filter(comparison == "EGF1min_vs_EGF0min") %>% 
               filter(diff > 1) %>% 
               filter(adj_pval <= 0.01),
             mapping = aes(x = diff, y = -log10(adj_pval), color = mod_res), alpha = 0.3, shape =19) +
  facet_wrap(facets = vars(comparison, mod_res)) +
  # geom_text_repel(data = psite_diff_abundance_gene %>%
  #              filter(comparison == "EGF1min_vs_EGF0min") %>%
  #              filter(diff >= 1) %>%
  #              filter(adj_pval <= 0.05),
  #              mapping = aes(x = diff, y = -log10(adj_pval), label = gene),
  #              max.overlaps = 200,
  #              segment.color = NA,
  #              # position = position_stack( hjust = 6),
  #              nudge_x = 8,
  #              force_pull   = 0, # do not pull toward data points
  #              # xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),
  #              xlim = c(0, 10), ylim = c(-0.5, 21),
  #              direction = "y",
  #              # vjust = 1,
  #              hjust = "right",
  #              size = 2) +
  # geom_text_repel(data = psite_diff_abundance_gene %>%
  #              filter(comparison == "EGF1min_vs_EGF0min") %>%
  #              filter(diff <= -1.5) %>%
  #              filter(adj_pval <= 0.01),
  #              mapping = aes(x = diff, y = -log10(adj_pval), label = gene),
  #              max.overlaps = 200,
  #              segment.color = NA,
  #              # position = position_stack( hjust = 6),
  #              nudge_x = -8,
  #              direction = "y",
  #              force_pull   = 0, # do not pull toward data points
  #              # xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),
  #              xlim = c(-10, 0), ylim = c(-0.5, 21),
  #              # vjust = 1,
  #              hjust = "left",
  #              size = 2) +
  alexis_theme() +
  coord_cartesian(clip = 'off') +

  scale_x_continuous(name = "Log2(fold change)", breaks = seq(-8, 8, 2), limits = c(-8, 8)) +
  scale_y_continuous(name = "-Log10(adj. p-value)", breaks = seq(0,20,5), limits = c(0, 20))

my_psite_volcano_plot_EGF1min_vs_0min_STYcolor

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/my_psite_volcano_plot_EGF1min_vs_0min_STYcolor.png", plot = my_psite_volcano_plot_EGF1min_vs_0min_STYcolor, width = 15, height = 10, scale = 0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/my_psite_volcano_plot_EGF1min_vs_0min_STYcolor.pdf", plot = my_psite_volcano_plot_EGF1min_vs_0min_STYcolor, width = 15, height = 10, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
pheatmap_diff_pval_df <- psite_diff_abundance_gene %>%
  filter(comparison %in% c("EGF15min_vs_EGF0min", "EGF5min_vs_EGF0min", "EGF3min_vs_EGF0min", "EGF1min_vs_EGF0min")) %>%
  
  # #flip the 0 min and 15 min comparison
  # mutate(diff = case_when(
  #   comparison == "EGF0min_vs_EGF15min" ~ -diff,
  #   TRUE ~ diff)) %>% 
  # mutate(comparison = case_when(
  #   comparison == "EGF0min_vs_EGF15min" ~ "EGF15min_vs_EGF0min", 
  #   TRUE ~ comparison)) %>% 
  
  
  mutate(gene_psite = paste(gene, modres_position, sep = " ")) %>% 
  filter(abs(diff) > 1) %>% 
  filter(adj_pval <= 0.05) %>% 
  pivot_wider(id_cols = c(gene_psite, protein_names, mod_res), names_from = comparison, values_from = c(diff, adj_pval), values_fill = 0) %>% 
  select(gene_psite, diff_EGF1min_vs_EGF0min, diff_EGF3min_vs_EGF0min, diff_EGF5min_vs_EGF0min, diff_EGF15min_vs_EGF0min, everything()) #%>%
  # filter(mod_res == "Y")


## --------------------------------------------------------------------------------------------------------------------------------------------------
# #add in control for range of log2FC to center colors
# control_rows_df <- data.frame(row.names = c("CTRL S0", "CTRL S100", "CTRL T0", "CTRL T100","CTRL Y0", "CTRL Y100"),
#                               diff_EGF1min_vs_EGF0min = c(-6, 6, -6, 6, -6, 6),
#                               diff_EGF3min_vs_EGF0min = c(-6, 6, -6, 6, -6, 6),
#                               diff_EGF5min_vs_EGF0min = c(-6, 6, -6, 6, -6, 6),
#                               diff_EGF15min_vs_EGF0min = c(-6, 6, -6, 6, -6, 6))


## --------------------------------------------------------------------------------------------------------------------------------------------------
#add in control for range of log2FC to center colors
control_rows_df <- data.frame(row.names = c("CTRL S0", "CTRL S100"),
                              diff_EGF1min_vs_EGF0min = c(-6, 6),
                              diff_EGF3min_vs_EGF0min = c(-6, 6),
                              diff_EGF5min_vs_EGF0min = c(-6, 6),
                              diff_EGF15min_vs_EGF0min = c(-6, 6))




# assign df columns to matrix values then row names:
# pheatmap_diff_pval_matrix <- as.matrix(pheatmap_diff_pval_matrix_wCTRL[,4:13]) #column range for all conditions

#may need to uncomment:
pheatmap_diff_pval_matrix_wCTRL <- as.matrix(pheatmap_diff_pval_df[,2:5])

rownames(pheatmap_diff_pval_matrix_wCTRL) <- as.matrix(pheatmap_diff_pval_df[,1])


pheatmap_diff_pval_matrix_wCTRL <- as.matrix(rbind(control_rows_df, pheatmap_diff_pval_matrix_wCTRL))



## --------------------------------------------------------------------------------------------------------------------------------------------------
# annotate_row_pSTY <- pheatmap_diff_pval_df %>% select(gene_psite, mod_res) %>% 
#   mutate(mod_res = case_when(
#     mod_res == "S" ~ "pSer",
#     mod_res == "T" ~ "pThr",
#     mod_res == "Y" ~ "pTyr" ))


## --------------------------------------------------------------------------------------------------------------------------------------------------
  my_colour = list(
    # sample = c(normal = "#5977ff", tumour = "#f74747"),
    mod_res = c(pSer = "purple3", pThr = "darkcyan", pTyr = "darkgoldenrod1" ))
    # random = c(random1 = "#82ed82", random2 = "#9e82ed"),
    # cluster = c(cluster1 = "#e89829", cluster2 = "#cc4ee0")


## --------------------------------------------------------------------------------------------------------------------------------------------------

signif_change_pprotein_heatmap <- pheatmap(pheatmap_diff_pval_matrix_wCTRL,
                                           # annotation_row = annotate_row_pSTY,
                                           cluster_cols = FALSE,
                                           cluster_rows = TRUE,
                                           colorRampPalette(c("purple", "white", "orange"))(48),
                                           annotation_colors = my_colour,
                                           legend_breaks = c(-6,0, 6),
                                           display_numbers = TRUE,
                                           fontsize_number = 6,
                                           cutree_rows = 6, cutree_cols = 2,
                                           fontsize_row = 8, fontsize_col = 4, angle_col = 270,
                                           border_color = "black" )
# 
ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/signif_change_pprotein_heatmap_pSTY.png", plot = signif_change_pprotein_heatmap, height = 8, width = 8, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
pheatmap_diff_pval_df_all <- psite_diff_abundance_gene %>%
  filter(comparison %in% c("EGF15min_vs_EGF0min", "EGF5min_vs_EGF0min", "EGF3min_vs_EGF0min", "EGF1min_vs_EGF0min")) %>%
  # 
  # #flip the 0 min and 15 min comparison
  # mutate(diff = case_when(
  #   comparison == "EGF0min_vs_EGF15min" ~ -diff,
  #   TRUE ~ diff)) %>% 
  # mutate(comparison = case_when(
  #   comparison == "EGF0min_vs_EGF15min" ~ "EGF15min_vs_EGF0min", 
  #   TRUE ~ comparison)) %>% 
  
  
  mutate(gene_psite = paste(gene, modres_position, sep = " ")) %>% 
  # filter(abs(diff) > 1) %>% 
  # filter(adj_pval <= 0.05) %>% 
  pivot_wider(id_cols = c(gene_psite, protein_names, mod_res), names_from = comparison, values_from = c(diff, adj_pval), values_fill = 0) %>% 
  #change NA to 0:
  mutate(
    diff_EGF1min_vs_EGF0min = ifelse(is.na(diff_EGF1min_vs_EGF0min), 0, diff_EGF1min_vs_EGF0min),
    diff_EGF3min_vs_EGF0min = ifelse(is.na(diff_EGF3min_vs_EGF0min), 0, diff_EGF3min_vs_EGF0min),
    diff_EGF5min_vs_EGF0min = ifelse(is.na(diff_EGF5min_vs_EGF0min), 0, diff_EGF5min_vs_EGF0min),
    diff_EGF15min_vs_EGF0min = ifelse(is.na(diff_EGF15min_vs_EGF0min), 0, diff_EGF15min_vs_EGF0min)) %>%  
  select(gene_psite, diff_EGF1min_vs_EGF0min, diff_EGF3min_vs_EGF0min, diff_EGF5min_vs_EGF0min, diff_EGF15min_vs_EGF0min, everything()) #%>% 
  # filter(mod_res == "Y")


# assign df colums to matrix values then row names:
# pheatmap_diff_pval_matrix <- as.matrix(pheatmap_diff_pval_df[,4:13]) #column range for all conditions
pheatmap_diff_pval_matrix_all <- as.matrix(pheatmap_diff_pval_df_all[,2:5])

rownames(pheatmap_diff_pval_matrix_all) <- as.matrix(pheatmap_diff_pval_df_all[,1])

pheatmap_diff_pval_matrix_all_wCTRL <- as.matrix(rbind(control_rows_df, pheatmap_diff_pval_matrix_all))



## --------------------------------------------------------------------------------------------------------------------------------------------------



## --------------------------------------------------------------------------------------------------------------------------------------------------
signif_change_pprotein_heatmap_all <- pheatmap(pheatmap_diff_pval_matrix_all_wCTRL,
                                           cluster_cols = FALSE,
                                           cluster_rows = TRUE,
                                           colorRampPalette(c("purple", "white", "orange"))(48),
                                   legend = TRUE, legend_breaks = c(-6,0, 6),
                                   cutree_rows = 6, cutree_cols = 4,
                                   display_numbers = TRUE,
                                           fontsize_number = 6,
                                   fontsize_row = 6, fontsize_col = 4, angle_col = 270,
                                   border_color = "black" )

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/signif_change_pprotein_heatmap_all_STY.png", plot = signif_change_pprotein_heatmap_all, height = 16, width = 7, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
ascore_final_gene <- ascore_fasta_protein_mod_loc_sum_psite %>% 
  left_join(y = (fasta_gene_names %>% select(entry, gene, protein_names, keyword_id) %>% rename(reference = entry)), by = "reference") %>% 
  separate(ref, into = c("reference2", "psite"), sep = "_", remove = FALSE) %>% 
  select(-reference2) %>% 
  mutate(
    gene_ref = paste(gene, psite, sep = " "),
    condition = as.factor(condition),
    condition = fct_relevel(condition, "EGF0min", "EGF1min",  "EGF3min", "EGF5min", "EGF15min")) %>% 
  group_by(condition, ref, gene_ref) %>% 
  mutate(
    n_obs = n() ) %>% 
  ungroup()


## --------------------------------------------------------------------------------------------------------------------------------------------------
most_common_protein_redundancy <- ascore_fasta_protein_mod_loc %>% 
  mutate(redundancy = as.numeric(redundancy)) %>% 
  # distinct(condition, replicate, ref, mod_protein_location, proteins, redundancy) %>%
  group_by(ref, mod_protein_location, proteins, redundancy) %>% 
  mutate(n_obs_each_redundancy = n()) %>% 
  ungroup() %>% 
  group_by( ref, mod_protein_location) %>% 
  filter(n_obs_each_redundancy == max(n_obs_each_redundancy)) %>% 
  ungroup() %>% 
  distinct(sample_id, condition, replicate, ref, mod_protein_location, proteins, redundancy,n_obs_each_redundancy) #%>% 
  # group_by(condition, replicate, ref) %>% 
  # mutate(distinct_redundancy = n_distinct(redundancy)) %>% 
  # ungroup()

#now with a distinct and most observed redundancy assignment per psite across all conditions & replicates, let's join to the quantitative dataframe with this protein and gene redundancy info.

global_phospho_shiny <- ascore_final_gene %>% 
  mutate(gene_ref = str_replace_all(gene_ref, " ", "_")) %>% 
  left_join(most_common_protein_redundancy) %>% 
  left_join(fasta_gene_names %>% distinct(entry, gene_names) %>% rename(reference = entry)) %>% 
  rename(gene_redundancy = gene_names) %>% 
  mutate(redundancy = as.numeric(redundancy),
         ambiguous = case_when(redundancy == 0 ~ "single protein assignment",
                               redundancy >=1 ~ paste("ambiguous protein assignment,", redundancy +1, "proteins", sep = " ")))
  




## --------------------------------------------------------------------------------------------------------------------------------------------------
#global phospho shiny defined above with most common redundancy info per site and gene names and protein assignment tag 

write_csv(x = global_phospho_shiny, file = "modified_data/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/global_pSTY_SITE_input.csv")


## --------------------------------------------------------------------------------------------------------------------------------------------------
pY_boxplot_points <- ggplot() + 
  # pY_boxplot_points <- ggplot(data = ascore_final_gene %>% filter(gene %in% c("EGFR", "MAPK1", "STAT3"))) +
  geom_boxplot(data = ascore_final_gene %>% filter(gene %in% c("EGFR")),
               mapping = aes(x = condition, y = median_norm_intensity), show.legend = FALSE, alpha = 0.3, outlier.color = NA) +
  geom_point(data = ascore_final_gene %>% filter(gene %in% c("EGFR")),
             mapping = aes(x = condition, y = median_norm_intensity, fill = sample_id), shape = 21, alpha = 0.6, show.legend = FALSE) +
  geom_text(data = ascore_final_gene %>% filter(gene %in% c("EGFR")) %>% distinct(condition, gene_ref, n_obs),
             aes(x = condition, y = 25, label = n_obs),size = 2.5, show.legend = FALSE, fontface = "bold") +
  alexis_theme() +
  facet_wrap(facets = vars(gene_ref), nrow = 3)

pY_boxplot_points


ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pY_boxplot_points_EGFR.png", plot = pY_boxplot_points, height = 12, width = 20, scale = 0.4)
ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pY_boxplot_points_EGFR.pdf", plot = pY_boxplot_points, height = 12, width = 20, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
pY_boxplot_points_MAPK1 <- ggplot() + 
  # pY_boxplot_points <- ggplot(data = ascore_final_gene %>% filter(gene %in% c("EGFR", "MAPK1", "STAT3"))) +
  geom_boxplot(data = ascore_final_gene %>% filter(gene %in% c("MAPK1")), 
               mapping = aes(x = condition, y = median_norm_intensity,), show.legend = FALSE, alpha = 0.3, outlier.color = NA) +
  geom_point(data = ascore_final_gene %>% filter(gene %in% c("MAPK1")),
             mapping = aes(x = condition, y = median_norm_intensity, fill = sample_id), shape = 21, alpha = 0.6, show.legend = FALSE) +
  
  alexis_theme() +
  facet_wrap(facets = vars(gene_ref), nrow = 3) +
  geom_text(data = ascore_final_gene %>% filter(gene %in% c("MAPK1")) %>% distinct(condition, gene_ref, n_obs),
             aes(x = condition, y = 30, label = n_obs),size = 2.5, show.legend = FALSE, fontface = "bold") 

pY_boxplot_points_MAPK1


ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pY_boxplot_points_MAPK1.png", plot = pY_boxplot_points_MAPK1, height = 8, width = 6, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# boxplot_gene_psite_function <- function(vector_gene_refs, dataframe, output_dir = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/", condition = condition, sample_id = sample_id){
# 
#   #create output directory if it doesn't exist already
#   if (!dir.exists(output_dir)){
#     dir.create(output_dir)}
#   
#   for (gene_var in vector_gene_refs) {
#     if(gene_var %in% dataframe$gene_ref) {
#       
#       #create boxplot for any sites in each gene
#       p <- ggplot() + 
#         geom_boxplot(data = dataframe %>% filter(gene_ref == gene_var),
#                mapping = aes(x = condition, y = median_norm_intensity), show.legend = FALSE, alpha = 0.3, outlier.color = NA, width = 0.5) +
#         geom_point(data = dataframe %>% filter(gene_ref == gene_var),
#                    mapping = aes(x = condition, y = median_norm_intensity, fill = sample_id), position = position_jitter( width = 0.3, height = 0), shape = 21, alpha = 0.4, size = 2, show.legend = FALSE) +
#         geom_text(data = dataframe %>%
#         filter(gene_ref == gene_var) %>%
#         distinct(condition, gene_ref, n_obs),
#         aes(x = condition, y = 34, label = n_obs),size = 4, show.legend = FALSE, fontface = "bold") +
#         alexis_theme() +
#         # theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
#         facet_wrap(facets = vars(gene_ref), nrow = 4) +
#         scale_y_continuous( expand = c(0,0))  +
#         expand_limits(y = c(10, 35)) 
#       
#       # save the plot as PNG & PDF
#       ggsave(filename = file.path(output_dir, paste0(gene_var, "_boxplot_points.png")), plot = p, width = 10, height = 12, scale = 0.4)
#       ggsave(filename = file.path(output_dir, paste0(gene_var, "_boxplot_points.pdf")), plot = p, width = 6, height = 8, scale = 0.4)
#       
#     } else { 
#       message (paste("Gene", gene_var, "not found in the dataframe."))
#       
#       }
#     }
#   }


## --------------------------------------------------------------------------------------------------------------------------------------------------
# boxplot_gene_psite_function(vector_gene_refs = (ascore_final_gene %>%distinct(gene_ref))$gene_ref,
#                             dataframe = ascore_final_gene %>%
#                               distinct(gene, gene_ref, condition, replicate,
#                                        sample_id, median_norm_intensity, n_obs) %>%
#                               mutate(
#                                 condition = fct_relevel(condition,
#                                                         "EGF0min", "EGF1min", "EGF3min", "EGF5min", "EGF15min")),
#                             output_dir = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/boxplot_points/")


## --------------------------------------------------------------------------------------------------------------------------------------------------
# boxplot_gene_psite_function <- function(vector_gene_refs, dataframe, output_dir = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/boxplot_points2/", condition = condition, sample_id = sample_id, strip_color = "gray80", strip_text_color = "black"){
# 
#   #create output directory if it doesn't exist already
#   if (!dir.exists(output_dir)){
#     dir.create(output_dir)}
#   
#   for (gene_var in vector_gene_refs) {
#     if(gene_var %in% dataframe$gene_ref) {
#       
#       #create boxplot for any sites in each gene
#       p <- ggplot() + 
#         geom_boxplot(data = dataframe %>% filter(gene_ref == gene_var),
#                mapping = aes(x = condition, y = median_norm_intensity), show.legend = FALSE, alpha = 0.3, outlier.color = NA, width = 0.5) +
#         geom_point(data = dataframe %>% filter(gene_ref == gene_var),
#                    mapping = aes(x = condition, y = median_norm_intensity),
#                    fill = "grey",
#                    mod_residue = position_jitter( width = 0.15, height = 0), shape = 21, alpha = 0.5, size = 1.5, show.legend = FALSE) +
#         geom_text(data = dataframe %>%
#         filter(gene_ref == gene_var) %>%
#           mutate(max_intensity_condition = max(median_norm_intensity)) %>% 
#         distinct(condition, gene_ref, n_obs,max_intensity_condition),
#         aes(x = condition, y = max_intensity_condition + 1, label = n_obs),size = 4, show.legend = FALSE, fontface = "bold") +
#         alexis_theme() +
#         
#         facet_wrap(facets = vars(gene_ref), nrow = 4) +
#         scale_y_continuous(name = expression("Log" [2]* "(intensity)"), expand = c(0.05,0.2))  +
#         scale_x_discrete(
#           "EGF treatment duration\n(minutes)",
#           drop = FALSE,
#           # limits = c("untreated", "1 min", "3 min", "5 min", "15 min"),
#           labels = c("0", "1", "3", "5", "15")) +
#         # theme(axis.text.x = element_text(angle = -90, hjust = 0.5, vjust = 0.5)) +
#         theme(
#           
#           strip.background = element_rect(fill = adjustcolor( strip_color, alpha.f = 0.3)),
#           strip.text = element_text(size = 16, face = "bold", family = "sans", color = strip_text_color)) #+
#         # expand_limits(y = c(15, 35)) 
#       
#       # save the plot as PNG & PDF
#       ggsave(filename = file.path(output_dir, paste0(gene_var, "_boxplot_points1.png")), plot = p, width = 8, height = 8, scale = 0.4)
#       ggsave(filename = file.path(output_dir, paste0(gene_var, "_boxplot_points1.pdf")), plot = p, width = 8, height = 8, scale = 0.4)
#       
#     } else { 
#       message (paste("Gene", gene_var, "not found in the dataframe."))
#       
#       }
#     }
# }




## --------------------------------------------------------------------------------------------------------------------------------------------------

# boxplot_gene_psite_function( vector_gene_refs = (ascore_gene_clusters_Imputation_noSD %>% filter(cluster ==1) %>% distinct(gene_ref))$gene_ref,strip_text_color = "black",
#                              dataframe = ascore_gene_clusters_Imputation_noSD %>% distinct(gene, gene_ref, condition, replicate, sample_id, median_norm_intensity, cluster, n_obs) %>% mutate(condition = fct_relevel(condition, "EGF0min", "EGF1min", "EGF3min", "EGF5min", "EGF15min")), output_dir = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/cluster1/", strip_color = "#081D58")




## --------------------------------------------------------------------------------------------------------------------------------------------------
# boxplot_gene_psite_function(vector_gene_refs = (ascore_final_gene %>%distinct(gene_ref))$gene_ref,
#                             dataframe = ascore_final_gene %>%
#                               distinct(gene, gene_ref, condition, replicate,
#                                        sample_id, median_norm_intensity, n_obs) %>%
#                               mutate(
#                                 condition = fct_relevel(condition,
#                                                         "EGF0min", "EGF1min", "EGF3min", "EGF5min", "EGF15min")),
#                             output_dir = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/boxplot_points2/")


## --------------------------------------------------------------------------------------------------------------------------------------------------
# boxplot_gene_psite_function1 <- function(vector_gene_refs, dataframe, output_dir = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/", condition = condition, sample_id = sample_id, strip_color = "gray80", strip_text_color = "black"){
# 
#   #create output directory if it doesn't exist already
#   if (!dir.exists(output_dir)){
#     dir.create(output_dir)}
#   
#   for (gene_var in vector_gene_refs) {
#     if(gene_var %in% dataframe$gene_ref) {
#       
#       #create boxplot for any sites in each gene
#       p <- ggplot() + 
#         geom_boxplot(data = dataframe %>% filter(gene_ref == gene_var),
#                mapping = aes(x = condition, y = median_norm_intensity), show.legend = FALSE, alpha = 0.3, outlier.color = NA, width = 0.5) +
#         geom_point(data = dataframe %>% filter(gene_ref == gene_var),
#                    mapping = aes(x = condition, y = median_norm_intensity),
#                    fill = "grey",
#                    mod_residue = position_jitter( width = 0.15, height = 0), shape = 21, alpha = 0.5, size = 1.5, show.legend = FALSE) +
#         geom_text(data = dataframe %>%
#         filter(gene_ref == gene_var) %>%
#           mutate(max_intensity_condition = max(median_norm_intensity)) %>% 
#         distinct(condition, gene_ref, n_obs,max_intensity_condition),
#         aes(x = condition, y = max_intensity_condition + 1, label = n_obs),size = 4, show.legend = FALSE, fontface = "bold") +
#         alexis_theme() +
#         
#         facet_wrap(facets = vars(gene_ref), nrow = 4) +
#         scale_y_continuous(name = expression("Log" [2]* "(intensity)"), expand = c(0.1,0.2))  +
#         scale_x_discrete(
#           "EGF treatment duration\n(minutes)",
#           drop = FALSE,
#           # limits = c("untreated", "1 min", "3 min", "5 min", "15 min"),
#           labels = c("0", "1", "3", "5", "15")) +
#         # theme(axis.text.x = element_text(angle = -90, hjust = 0.5, vjust = 0.5)) +
#         theme(
#           
#           strip.background = element_rect(fill = adjustcolor( strip_color, alpha.f = 0.3)),
#           strip.text = element_text(size = 16, face = "bold", family = "sans", color = strip_text_color)) #+
#         # expand_limits(y = c(15, 35)) 
#       
#       # save the plot as PNG & PDF
#       ggsave(filename = file.path(output_dir, paste0(gene_var, "_boxplot_points.png")), plot = p, width = 8, height = 8, scale = 0.4)
#       ggsave(filename = file.path(output_dir, paste0(gene_var, "_boxplot_points.pdf")), plot = p, width = 8, height = 8, scale = 0.4)
#       
#     } else { 
#       message (paste("Gene", gene_var, "not found in the dataframe."))
#       
#       }
#     }
# }




## --------------------------------------------------------------------------------------------------------------------------------------------------

CST_EGFR_proteins <- c("P00533","P42336","Q9Y243","P42345","P23443","P15498","P63000","Q13153","Q13233","P45983","P05412","P01100","Q05397","P56945","P16333","P60953","P15941","P35222","P56539","P49023","P12830","P42226","P46108","P00519","Q9UPR0","P0DP24","Q16566","Q9UQM7","P22681","P62993","Q96J02","O14964","Q07889","O14807","P04049","Q02750","P01112","P28482","P05771","P12931","P01116","P52333","P07947","P43405","Q9UQC2","Q06124","Q99704","Q7L591","P42224","Q18PE1","P31751","P31749","P53779","P45984","Q6FG41","Q14289","P51636","Q03135","P52630","Q14765","P0DP25","Q13557","Q96NX5","P0DP23","Q13555","Q13554","P05129","P23458","Q6PKX4","Q9P104","Q07890","O60496","Q8TEW6","P36507","O60674","P17252","Q9NQ66","P19174","P51178"," Q4KWH8","P16885","O75038","Q00722","Q9UJM3","P52735","Q9UKW4","P15498","P16220","Q02930","O43889","P05771","P17252","P05129", "P19419", "P41970", "P28324", "P42336", "O00443", "O00750", "P23443")


## --------------------------------------------------------------------------------------------------------------------------------------------------

# boxplot_gene_psite_function1( vector_gene_refs = (ascore_final_gene %>%
#                                                    filter(reference %in% c("P19419", "P41970", "P28324")) %>% 
#                                                     # filter(grepl("STAT", gene) == TRUE) %>% 
#                                                    # filter(reference %in% CST_EGFR_proteins) %>%
#                                                    distinct(gene_ref))$gene_ref,strip_text_color = "black",
#                              dataframe = ascore_final_gene %>% distinct(gene, gene_ref, condition, replicate, sample_id, median_norm_intensity, n_obs) %>% mutate(condition = fct_relevel(condition, "EGF0min", "EGF1min", "EGF3min", "EGF5min", "EGF15min")), output_dir = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/EGF_pathway/", strip_color = "lightblue")




## --------------------------------------------------------------------------------------------------------------------------------------------------
PCA_df <- ascore_fasta_protein_mod_loc_sum_psite %>%
  # filter(!is.na(missingness)) %>% 
  filter(median_norm_intensity != -Inf)


all_psite_scree  <- qc_pca(
  data = PCA_df,
  sample = sample_id,
  grouping = ref,
  intensity = median_norm_intensity,
  condition = condition,
  digestion = NULL,
  plot_style = "scree")

all_psite_scree


## --------------------------------------------------------------------------------------------------------------------------------------------------
all_psite_PC_1_2  <- qc_pca(
  data = PCA_df,
  sample = sample_id,
  grouping = ref,
  intensity = median_norm_intensity,
  condition = condition,
  digestion = NULL,
  plot_style = "pca",
  components = c("PC1", "PC2"))

all_psite_PC_1_2

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/all_psite_PC_1_2.png", plot = all_psite_PC_1_2, width = 24, height = 12, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
all_psite_PC_1_3  <- qc_pca(
  data = PCA_df,
  sample = sample_id,
  grouping = ref,
  intensity = median_norm_intensity,
  condition = condition,
  digestion = NULL,
  plot_style = "pca",
  components = c("PC1", "PC3"))

all_psite_PC_1_3

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/all_psite_PC_1_3.png", plot = all_psite_PC_1_3, width = 24, height = 12, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
ascore_EGFR1_PTMSEA <- ascore_fasta_protein_mod_loc_sum_psite %>%
  separate(ref, into = c("reference2", "psite"), sep = "_", remove = FALSE) %>% 
  select(-reference2) %>% 
  left_join(y = fasta_gene_names %>% rename(reference = entry) %>% select(reference, gene, protein_names, keyword_id), by = "reference") %>% 
  full_join(y = EGFR1_pathway , by = c("gene", "reference", "mod_res", "psite"), relationship = "many-to-many") %>% 
mutate(PTM_SEA = case_when(
  is.na(PTM_SEA) ~ FALSE, 
  TRUE ~ TRUE)) %>%   
mutate(
    overlap = case_when(
    !is.na(condition) & PTM_SEA == TRUE ~ "both",
    PTM_SEA == TRUE ~ "literature",
    PTM_SEA == FALSE  ~ "my data")) %>%
  distinct(overlap, gene, psite, mod_res)
  


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_PTMSEA_EGFR_pathway_modres <- ggplot(data = EGFR1_pathway %>% distinct(gene, psite, mod_res)) +
  geom_bar(mapping = aes(x = mod_res, fill = mod_res), color = "black", size = 0.5, show.legend = FALSE) +
  alexis_theme() +
  theme(axis.text.x = element_text(angle = 0, size = 12, hjust = 0.5)) +
  scale_fill_viridis_d() +
  ggtitle(label = "PTM-SEA: EGFR pathway \np-site distribution") +
  ylab("unique gene + p-site")
  

plot_PTMSEA_EGFR_pathway_modres

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/literature_PTMSEA_EGFR_pathway_modres_distribution.png", plot = plot_PTMSEA_EGFR_pathway_modres, width = 6, height = 6, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
EGFR1_overlaps <- ascore_EGFR1_PTMSEA %>% 
  distinct(overlap, gene, psite, mod_res)

#plot overlaps count of psites in my data vs. PTMSEA EGFR1 pathway
EGFR1_psite_overlaps <- ggplot(data = EGFR1_overlaps ) +
  geom_bar(mapping = aes(x = overlap, fill = mod_res), color = "black", size = 0.5) +
  scale_fill_viridis_d(direction = 1) +
  alexis_theme() +
  scale_y_continuous(expand = c(0,0)) +
  # expand_limits(y = c(0, 500)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  ylab("unqiue p-sites")


EGFR1_psite_overlaps


ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/EGFR1_psite_overlaps_pSTY.png", plot = EGFR1_psite_overlaps,
       scale = 0.4, width = 8, height = 6)



## --------------------------------------------------------------------------------------------------------------------------------------------------

ascore_EGFR1_PTMSEA_protein <- ascore_fasta_protein_mod_loc_sum_psite %>% 
  left_join(y = fasta_gene_names %>% rename(reference = entry) %>% select(reference, gene, protein_names, keyword_id), by = "reference") %>%
  full_join(y = EGFR1_pathway %>% select(gene, PTM_SEA) , by = c("gene"), relationship = "many-to-many") %>%
  
  mutate(PTM_SEA = case_when(
  is.na(PTM_SEA) ~ FALSE, 
  TRUE ~ TRUE)) %>% 
  
  mutate(overlap = case_when(
    !is.na(condition) & PTM_SEA == TRUE ~ "both",
    is.na(condition) & PTM_SEA == TRUE ~ "literature",
    !is.na(condition) & PTM_SEA == FALSE  ~ "my data"))

EGFR1_overlaps_protein <- ascore_EGFR1_PTMSEA_protein %>% 
  distinct(overlap, gene)

#plot overlaps count of psites in my data vs. PTMSEA EGFR1 pathway
EGFR1_overlaps_protein <- ggplot(data = EGFR1_overlaps_protein) +
  geom_bar(mapping = aes(x = overlap)) +
  alexis_theme() +
  scale_y_continuous(expand = c(0,0)) +
  expand_limits(y = c(0, 4000)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  ylab("unqiue p-proteins")


EGFR1_overlaps_protein


ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/EGFR1_overlaps_protein.png", plot = EGFR1_overlaps_protein,
       scale = 0.4, width = 6, height = 8)

#plot
# plot_intensity_psites_EGFR1 <- ggplot(data = EGFR1_overlaps) +
#   geom_boxplot(mapping = aes())


## --------------------------------------------------------------------------------------------------------------------------------------------------
EGFR1_overlaps <- ascore_EGFR1_PTMSEA %>% 
  distinct(overlap, gene, psite, mod_res)

#plot overlaps count of psites in my data vs. PTMSEA EGFR1 pathway
EGFR1_psite_overlaps <- ggplot(data = EGFR1_overlaps) +
  geom_bar(mapping = aes(x = overlap, fill = mod_res)) +
  scale_fill_viridis_d() +
  alexis_theme() +
  scale_y_continuous(expand = c(0,0)) +
  # expand_limits(y = c(0, 12000)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  ylab("unqiue p-sites")
EGFR1_psite_overlaps
ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/EGFR1_psite_overlaps.png", plot = EGFR1_psite_overlaps,
       scale = 0.4, width = 6, height = 8)

#plot
plot_intensity_psites_EGFR1 <- ggplot(data = EGFR1_overlaps) +
  geom_boxplot(mapping = aes())


## --------------------------------------------------------------------------------------------------------------------------------------------------
psite_diff_abundance_gene %>%
  filter(comparison %in% c("EGF0min_vs_EGF15min", "EGF5min_vs_EGF0min", "EGF3min_vs_EGF0min", "EGF1min_vs_EGF0min")) %>%
  
  #flip the 0 min and 15 min comparison
  mutate(diff = case_when(
    comparison == "EGF0min_vs_EGF15min" ~ -diff,
    TRUE ~ diff)) %>% 
  mutate(comparison = case_when(
    comparison == "EGF0min_vs_EGF15min" ~ "EGF15min_vs_EGF0min", 
    TRUE ~ comparison)) 


## --------------------------------------------------------------------------------------------------------------------------------------------------
PTMNav_psite_input <- psite_diff_abundance_gene %>%
  
  #select only comparisons relative to 0 min:
  filter(comparison %in% c("EGF0min_vs_EGF15min",
                           "EGF5min_vs_EGF0min",
                           "EGF3min_vs_EGF0min",
                           "EGF1min_vs_EGF0min")) %>%
  
  #flip direction of fold change and label for 15 min vs. 0 min:
  mutate(diff = case_when(
    comparison == "EGF0min_vs_EGF15min" ~ -diff,
    TRUE ~ diff)) %>% 
  mutate(comparison = case_when(
    comparison == "EGF0min_vs_EGF15min" ~ "EGF15min_vs_EGF0min", 
    TRUE ~ comparison)) %>% 
  
  distinct(reference, modres_position, comparison, diff, adj_pval) %>%
  filter(!is.na(reference)) %>% 
  mutate(
    regulation = case_when(
      adj_pval <= 0.05 & diff > 0 ~ "up",
      adj_pval <= 0.05 & diff < 0 ~ "down",
      adj_pval > 0.05 ~ "") )


colnames(PTMNav_psite_input) <- c("Protein IDs", "p-site", "Experiment", "Fold change", "adjusted pvalue", "Regulation")

write_csv(x = PTMNav_psite_input, file = "modified_data/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/PTMNav_psite_input.csv")


## --------------------------------------------------------------------------------------------------------------------------------------------------
locations_FASTA_df <- fasta_subcellular_location %>% 
   
  mutate(
    extracellular_TM = case_when(
      grepl("extracellular|transmembrane", gene_ontology_cellular_component) == TRUE ~ TRUE,
      TRUE ~ FALSE),
    cell_membrane = case_when(
      grepl("cell membrane|extracellular|transmembrane|plasma", gene_ontology_cellular_component) == TRUE ~ TRUE,
      TRUE ~ FALSE),
    nucleus = case_when(
      grepl("nucleus|nucle", gene_ontology_cellular_component) == TRUE ~ TRUE,
      TRUE ~ FALSE),
    mitochondria = case_when(
      grepl("mitochondr", gene_ontology_cellular_component) == TRUE ~ TRUE,
      TRUE ~ FALSE),
    cytoplasm = case_when(
      grepl("cytoplasm|cytosol", gene_ontology_cellular_component) == TRUE ~ TRUE,
      TRUE ~ FALSE),
    ER_Golgi = case_when(
      grepl("ndoplasmic|olgi", gene_ontology_cellular_component) == TRUE ~ TRUE,
      TRUE ~ FALSE)) %>% 
  mutate(comparison = "FASTA")


## --------------------------------------------------------------------------------------------------------------------------------------------------
locations_FASTA_df_ratios <- fasta_subcellular_location %>% 
   
  mutate(
    extracellular_TM = case_when(
      grepl("extracellular|transmembrane", gene_ontology_cellular_component) == TRUE ~ TRUE,
      TRUE ~ FALSE),
    cell_membrane = case_when(
      grepl("cell membrane|extracellular|transmembrane|plasma", gene_ontology_cellular_component) == TRUE ~ TRUE,
      TRUE ~ FALSE),
    nucleus = case_when(
      grepl("nucleus|nucle", gene_ontology_cellular_component) == TRUE ~ TRUE,
      TRUE ~ FALSE),
    mitochondria = case_when(
      grepl("mitochondr", gene_ontology_cellular_component) == TRUE ~ TRUE,
      TRUE ~ FALSE),
    cytoplasm = case_when(
      grepl("cytoplasm|cytosol", gene_ontology_cellular_component) == TRUE ~ TRUE,
      TRUE ~ FALSE),
    ER_Golgi = case_when(
      grepl("ndoplasmic|olgi", gene_ontology_cellular_component) == TRUE ~ TRUE,
      TRUE ~ FALSE)) %>% 
  mutate(
    total_proteins = n_distinct(reference)) %>% 
  
  group_by(cell_membrane) %>% 
  mutate(
    ratio_membrane = n()/total_proteins) %>% 
  ungroup() %>% 
  mutate(comparison = "FASTA")


## --------------------------------------------------------------------------------------------------------------------------------------------------
psites_extracellular_TM_fasta <- ggplot(data = locations_FASTA_df) +
  geom_bar(mapping = aes(x = comparison, fill = extracellular_TM ), color = "black") +
  scale_fill_brewer(palette = "RdGy", direction = -1) +
  # scale_y_continuous(limits = c(0, 200), expand = c(0,0)) +
  alexis_theme()

psites_extracellular_TM_fasta

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/psites_extracellular_TM_fasta.png", plot = psites_extracellular_TM_fasta, scale= 0.4, width = 6, height = 6)

psites_extracellular_TM_fasta_fill <- ggplot(data = locations_FASTA_df) +
  geom_bar(mapping = aes(x = comparison, fill = extracellular_TM ), position = "fill", color = "black") +
  scale_fill_brewer(palette = "RdGy", direction = -1) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  alexis_theme()

psites_extracellular_TM_fasta_fill

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/psites_extracellular_TM_fasta_fill.png", plot = psites_extracellular_TM_fasta_fill, scale= 0.4, width = 6, height = 6)


## --------------------------------------------------------------------------------------------------------------------------------------------------

psites_cell_membrane_fasta <- ggplot(data = locations_FASTA_df) +
  geom_bar(mapping = aes(x = comparison, fill = cell_membrane ), color = "black") +
  scale_fill_brewer(palette = "RdGy", direction = -1) +
  # scale_y_continuous(limits = c(0, 200), expand = c(0,0)) +
  alexis_theme()

psites_cell_membrane_fasta

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/psites_cell_membrane_fasta.png", plot = psites_cell_membrane_fasta, scale= 0.4, width = 6, height = 6)

psites_cell_membrane_fasta_fill <- ggplot(data = locations_FASTA_df) +
  geom_bar(mapping = aes(x = comparison, fill = cell_membrane ), position = "fill", color = "black") +
  scale_fill_brewer(palette = "RdGy", direction = -1) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  alexis_theme()

psites_cell_membrane_fasta_fill

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/psites_cell_membrane_fasta_fill.png", plot = psites_cell_membrane_fasta_fill, scale= 0.4, width = 6, height = 6)


## --------------------------------------------------------------------------------------------------------------------------------------------------

psites_cytoplasm_fasta <- ggplot(data = locations_FASTA_df) +
  geom_bar(mapping = aes(x = comparison, fill = cytoplasm ), color = "black") +
  scale_fill_brewer(palette = "RdGy", direction = -1) +
  # scale_y_continuous(limits = c(0, 200), expand = c(0,0)) +
  alexis_theme()

psites_cytoplasm_fasta

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/psites_cytoplasm_fasta.png", plot = psites_cytoplasm_fasta, scale= 0.4, width = 6, height = 6)

psites_cytoplasm_fasta_fill <- ggplot(data = locations_FASTA_df) +
  geom_bar(mapping = aes(x = comparison, fill = cytoplasm ), position = "fill", color = "black") +
  scale_fill_brewer(palette = "RdGy", direction = -1) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  alexis_theme()

psites_cytoplasm_fasta_fill

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/psites_cytoplasm_fasta_fill.png", plot = psites_cytoplasm_fasta_fill, scale= 0.4, width = 6, height = 6)


## --------------------------------------------------------------------------------------------------------------------------------------------------

psites_nucleus_fasta <- ggplot(data = locations_FASTA_df) +
  geom_bar(mapping = aes(x = comparison, fill = nucleus ), color = "black") +
  scale_fill_brewer(palette = "RdGy", direction = -1) +
  # scale_y_continuous(limits = c(0, 200), expand = c(0,0)) +
  alexis_theme()

psites_nucleus_fasta

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/psites_nucleus_fasta.png", plot = psites_nucleus_fasta, scale= 0.4, width = 6, height = 6)

psites_nucleus_fasta_fill <- ggplot(data = locations_FASTA_df) +
  geom_bar(mapping = aes(x = comparison, fill = nucleus ), position = "fill", color = "black") +
  scale_fill_brewer(palette = "RdGy", direction = -1) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  alexis_theme()

psites_nucleus_fasta_fill

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/psites_nucleus_fasta_fill.png", plot = psites_nucleus_fasta_fill, scale= 0.4, width = 6, height = 6)


## --------------------------------------------------------------------------------------------------------------------------------------------------

psites_mitochondria_fasta <- ggplot(data = locations_FASTA_df) +
  geom_bar(mapping = aes(x = comparison, fill = mitochondria ), color = "black") +
  scale_fill_brewer(palette = "RdGy", direction = -1) +
  # scale_y_continuous(limits = c(0, 200), expand = c(0,0)) +
  alexis_theme()

psites_mitochondria_fasta

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/psites_mitochondria_fasta.png", plot = psites_mitochondria_fasta, scale= 0.4, width = 6, height = 6)

psites_mitochondria_fasta_fill <- ggplot(data = locations_FASTA_df) +
  geom_bar(mapping = aes(x = comparison, fill = mitochondria ), position = "fill", color = "black") +
  scale_fill_brewer(palette = "RdGy", direction = -1) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  alexis_theme()

psites_mitochondria_fasta_fill

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/psites_mitochondria_fasta_fill.png", plot = psites_mitochondria_fasta_fill, scale= 0.4, width = 6, height = 6)


## --------------------------------------------------------------------------------------------------------------------------------------------------
annotated_df_sum_to_psite_membrane <- ascore_fasta_protein_mod_loc_sum_psite %>% 
  left_join(locations_FASTA_df_ratios, by = "reference") %>% 
  mutate(
    membrane = case_when(
      grepl("Membrane", keywords) == TRUE ~ "membrane",
      grepl("membrane", gene_ontology_cellular_component) == TRUE ~ "membrane",
      TRUE ~ "named_protein"),
    transmembrane = case_when(
      grepl("Transmembrane", keywords) == TRUE ~ "transmembrane",
      grepl("transmembrane", gene_ontology_cellular_component) == TRUE ~ "transmembrane",
      TRUE ~ "named_protein"),
    extracellular = case_when(
      grepl("extracellular", keywords) == TRUE ~ "extracellular",
      grepl("extracellular", gene_ontology_cellular_component) == TRUE ~ "extracellular",
      TRUE ~ "named_protein"))
  

#keyword: membrane -------------------------------------
annotated_psite_summary_membrane <- annotated_df_sum_to_psite_membrane %>%
  distinct( condition, replicate, sample_id, ref, reference,  protein_names, membrane) %>% 
  group_by(condition, replicate, sample_id, membrane) %>% 
  summarize(
    n_annotated_psites = n()  ) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(condition, replicate, sample_id), names_from = membrane, values_from = n_annotated_psites) %>% 
  mutate(
    sum_psite_count = membrane + named_protein,
    ratio_membrane = membrane / sum_psite_count )

#experimental protein summary -------------------------------------
annotated_protein_summary_membrane<- annotated_df_sum_to_psite_membrane %>%
  distinct( condition, replicate, sample_id,  reference,  protein_names, membrane) %>% 
  group_by(condition, replicate, sample_id, membrane) %>% 
  summarize(
    n_annotated_proteins = n()  ) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(condition, replicate, sample_id), names_from = membrane, values_from = n_annotated_proteins) %>% 
  mutate(
    total_proteins = membrane + named_protein,
    ratio_membrane_proteins = membrane / total_proteins )


#reference Sros_GOs --------------------------------------
reference_membrane_summary <- fasta_subcellular_location %>% 
  mutate(
    membrane = case_when(
       grepl("Membrane", keywords) == TRUE ~ "membrane",
       grepl("membrane", gene_ontology_cellular_component) == TRUE ~ "membrane",
      TRUE ~ "named_protein")) %>%
  group_by(membrane) %>% 
  summarize(
    n_annotated_proteins = n()) %>% 
  ungroup() %>%
  mutate(condition = "reference") %>% 
  pivot_wider(id_cols = condition, names_from = membrane, values_from = n_annotated_proteins) %>% 
  mutate(
    total_proteins = membrane + named_protein,
    ratio_membrane_proteins = membrane / total_proteins)

plot_membrane_psite_ratio <- ggplot(data = annotated_protein_summary_membrane) +
  geom_col(mapping = aes(x = sample_id, y = ratio_membrane_proteins, fill = condition), color = "black", linewidth =1, alpha = 0.5) +
  geom_hline(yintercept = round((reference_membrane_summary$ratio_membrane_proteins[1]), 2), linetype = 2, linewidth = 1, alpha = 0.5, color = "darkred") +
  alexis_theme() +
  annotate(geom = "text", x = 34, y = (reference_membrane_summary$ratio_membrane_proteins[1]) + 0.02, label = "FASTA ratio", size = 3, color = "darkred") +
  alexis_theme() #+
  # scale_y_continuous(limits = c(0, 0.25), breaks = c(seq(0, 0.25, 0.05)), expand = c(0,0), name = "ratio membrane")

plot_membrane_psite_ratio

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_membrane_psite_ratio.png", plot = plot_membrane_psite_ratio, width = 16, height = 12, scale = 0.4)


#boxplot ---------------------------------------------------------
plot_membrane_phosphoprotein_ratio_box <- ggplot(data = annotated_protein_summary_membrane %>% mutate(component = "membrane")) +
  geom_boxplot(mapping = aes(x = component, y = ratio_membrane_proteins), color = "black", linewidth =0.5, alpha = 0.5) +
  geom_jitter(mapping = aes(x = component, y = ratio_membrane_proteins), color = "black", fill = "white", size = 0.5, alpha = 0.5, shape = 21, width = 0.25) +
  geom_hline(yintercept = round((reference_membrane_summary$ratio_membrane_proteins[1]), 2), linetype = 2, linewidth = 1, alpha = 0.5, color = "darkred") +
  annotate(geom = "text", x = 1.4, y = (reference_membrane_summary$ratio_membrane_proteins[1])  , label = "FASTA", size = 2, color = "darkred") +
  theme(
    # panel.border = element_rect(colour = "blue", fill = NA, linetype = 2),
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.x  = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = -60, vjust = 0.9, hjust = 0, size = 6, family = "sans"),
    # axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5),
    # axis.text.x = element_text(angle = -60, vjust = 0.75, hjust = 0, size = 12, family = "sans"),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black", face = "plain", family = "sans", size = 14),
    axis.title = element_text(colour = "black", family = "sans", size = 14),
    axis.ticks = element_line(colour = "black"),
    title = element_text(size = 8, hjust = 0.5),
    strip.text = element_text(size= 12, family = "sans"),
    # legend at the bottom 6)
    legend.position = "right") #+
  # scale_y_continuous(limits = c(0, 0.25), breaks = c(seq(0, 0.25, 0.05)), expand = c(0,0), name = "ratio membrane")

plot_membrane_phosphoprotein_ratio_box

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_membrane_phosphoprotein_ratio_box.png", plot = plot_membrane_phosphoprotein_ratio_box, width = 4, height = 6, scale = 0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_membrane_phosphoprotein_ratio_box.pdf", plot = plot_membrane_phosphoprotein_ratio_box, width = 4, height = 6, scale = 0.4)



## --------------------------------------------------------------------------------------------------------------------------------------------------
#keyword: membrane -------------------------------------
annotated_psite_summary_transmembrane <- annotated_df_sum_to_psite_membrane %>%
  distinct( condition, replicate, sample_id, ref, reference,  protein_names, transmembrane) %>% 
  group_by(condition, replicate, sample_id, transmembrane) %>% 
  summarize(
    n_annotated_psites = n()  ) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(condition, replicate, sample_id), names_from = transmembrane, values_from = n_annotated_psites) %>% 
  mutate(
    sum_psite_count = transmembrane + named_protein,
    ratio_transmembrane = transmembrane / sum_psite_count )

#experimental protein summary -------------------------------------
annotated_protein_summary_transmembrane<- annotated_df_sum_to_psite_membrane %>%
  distinct( condition, replicate, sample_id,  reference,  protein_names, transmembrane) %>% 
  group_by(condition, replicate, sample_id, transmembrane) %>% 
  summarize(
    n_annotated_proteins = n()  ) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(condition, replicate, sample_id), names_from = transmembrane, values_from = n_annotated_proteins) %>% 
  mutate(
    total_proteins = transmembrane + named_protein,
    ratio_transmembrane_proteins = transmembrane / total_proteins )


#reference Sros_GOs --------------------------------------
reference_transmembrane_summary <- fasta_subcellular_location %>% 
  mutate(
    transmembrane = case_when(
       grepl("Transmembrane", keywords) == TRUE ~ "transmembrane",
       grepl("transmembrane", gene_ontology_cellular_component) == TRUE ~ "transmembrane",
      TRUE ~ "named_protein")) %>%
  group_by(transmembrane) %>% 
  summarize(
    n_annotated_proteins = n()) %>% 
  ungroup() %>%
  mutate(condition = "reference") %>% 
  pivot_wider(id_cols = condition, names_from = transmembrane, values_from = n_annotated_proteins) %>% 
  mutate(
    total_proteins = transmembrane + named_protein,
    ratio_transmembrane_proteins = transmembrane / total_proteins)

plot_transmembrane_psite_ratio <- ggplot(data = annotated_protein_summary_transmembrane) +
  geom_col(mapping = aes(x = sample_id, y = ratio_transmembrane_proteins, fill = condition), color = "black", linewidth =1, alpha = 0.5) +
  geom_hline(yintercept = round((reference_transmembrane_summary$ratio_transmembrane_proteins[1]), 2), linetype = 2, linewidth = 1, alpha = 0.5, color = "darkred") +
  annotate(geom = "text", x = 20, y = (reference_transmembrane_summary$ratio_transmembrane_proteins[1]) + 0.02, label = "FASTA ratio", size = 3, color = "darkred") +
  alexis_theme() +
  scale_y_continuous(limits = c(0, 0.25), breaks = c(seq(0, 0.25, 0.05)), expand = c(0,0), name = "ratio transmembrane")

plot_transmembrane_psite_ratio

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_transmembrane_psite_ratio.png", plot = plot_transmembrane_psite_ratio, width = 16, height = 12, scale = 0.4)


#boxplot ---------------------------------------------------------
plot_transmembrane_phosphoprotein_ratio_box <- ggplot(data = annotated_protein_summary_transmembrane %>% mutate(component = "transmembrane")) +
  geom_boxplot(mapping = aes(x = component, y = ratio_transmembrane_proteins), color = "black", linewidth =0.5, alpha = 0.5) +
  geom_jitter(mapping = aes(x = component, y = ratio_transmembrane_proteins), color = "black", fill = "white", size = 0.5, alpha = 0.5, shape = 21, width = 0.25) +
  geom_hline(yintercept = round((reference_transmembrane_summary$ratio_transmembrane_proteins[1]), 2), linetype = 2, linewidth = 1, alpha = 0.5, color = "darkred") +
  annotate(geom = "text", x = 1.1, y = (reference_transmembrane_summary$ratio_transmembrane_proteins[1])  , label = "FASTA", size = 2, color = "darkred") +
  theme(
    # panel.border = element_rect(colour = "blue", fill = NA, linetype = 2),
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.x  = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = -60, vjust = 0.9, hjust = 0, size = 6, family = "sans"),
    # axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5),
    # axis.text.x = element_text(angle = -60, vjust = 0.75, hjust = 0, size = 12, family = "sans"),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black", face = "plain", family = "sans", size = 14),
    axis.title = element_text(colour = "black", family = "sans", size = 14),
    axis.ticks = element_line(colour = "black"),
    title = element_text(size = 8, hjust = 0.5),
    strip.text = element_text(size= 12, family = "sans"),
    # legend at the bottom 6)
    legend.position = "right") #+
  # scale_y_continuous(limits = c(0, 0.25), breaks = c(seq(0, 0.25, 0.05)), expand = c(0,0), name = "ratio transmembrane")

plot_transmembrane_phosphoprotein_ratio_box

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_transmembrane_phosphoprotein_ratio_box.png", plot = plot_transmembrane_phosphoprotein_ratio_box, width = 4, height = 6, scale = 0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_transmembrane_phosphoprotein_ratio_box.pdf", plot = plot_transmembrane_phosphoprotein_ratio_box, width = 4, height = 6, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
#keyword: extracellular -------------------------------------
annotated_psite_summary_extracellular <- annotated_df_sum_to_psite_membrane %>%
  distinct( condition, replicate, sample_id, ref, reference,  protein_names, extracellular) %>% 
  group_by(condition, replicate, sample_id, extracellular) %>% 
  summarize(
    n_annotated_psites = n()  ) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(condition, replicate, sample_id), names_from = extracellular, values_from = n_annotated_psites) %>% 
  mutate(
    sum_psite_count = extracellular + named_protein,
    ratio_extracellular = extracellular / sum_psite_count )

#experimental protein summary -------------------------------------
annotated_protein_summary_extracellular<- annotated_df_sum_to_psite_membrane %>%
  distinct( condition, replicate, sample_id,  reference,  protein_names, extracellular) %>% 
  group_by(condition, replicate, sample_id, extracellular) %>% 
  summarize(
    n_annotated_proteins = n()  ) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(condition, replicate, sample_id), names_from = extracellular, values_from = n_annotated_proteins) %>% 
  mutate(
    total_proteins = extracellular + named_protein,
    ratio_extracellular_proteins = extracellular / total_proteins )


#reference Sros_GOs --------------------------------------
reference_extracellular_summary <- fasta_subcellular_location %>% 
  mutate(
    extracellular = case_when(
       grepl("extracellular", keywords) == TRUE ~ "extracellular",
       grepl("extracellular", gene_ontology_cellular_component) == TRUE ~ "extracellular",
      TRUE ~ "named_protein")) %>%
  group_by(extracellular) %>% 
  summarize(
    n_annotated_proteins = n()) %>% 
  ungroup() %>%
  mutate(condition = "reference") %>% 
  pivot_wider(id_cols = condition, names_from = extracellular, values_from = n_annotated_proteins) %>% 
  mutate(
    total_proteins = extracellular + named_protein,
    ratio_extracellular_proteins = extracellular / total_proteins)

plot_extracellular_psite_ratio <- ggplot(data = annotated_protein_summary_extracellular) +
  geom_col(mapping = aes(x = sample_id, y = ratio_extracellular_proteins, fill = condition), color = "black", linewidth =1, alpha = 0.5) +
  geom_hline(yintercept = round((reference_extracellular_summary$ratio_extracellular_proteins[1]), 2), linetype = 2, linewidth = 1, alpha = 0.5, color = "darkred") +
  alexis_theme() +
  annotate(geom = "text", x = 34, y = (reference_extracellular_summary$ratio_extracellular_proteins[1]) + 0.02, label = "FASTA ratio", size = 3, color = "darkred") +
  alexis_theme()+
  scale_y_continuous(limits = c(0, 0.2), breaks = c(seq(0, 0.2, 0.05)), expand = c(0,0), name = "ratio extracellular")

plot_extracellular_psite_ratio

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_extracellular_psite_ratio.png", plot = plot_extracellular_psite_ratio, width = 16, height = 12, scale = 0.4)


#boxplot ---------------------------------------------------------
plot_extracellular_phosphoprotein_ratio_box <- ggplot(data = annotated_protein_summary_extracellular %>% mutate(component = "extracellular")) +
  geom_boxplot(mapping = aes(x = component, y = ratio_extracellular_proteins), color = "black", linewidth =0.5, alpha = 0.5) +
  geom_jitter(mapping = aes(x = component, y = ratio_extracellular_proteins), color = "black", fill = "white", size = 0.5, alpha = 0.5, shape = 21, width = 0.25) +
  geom_hline(yintercept = round((reference_extracellular_summary$ratio_extracellular_proteins[1]), 2), linetype = 2, linewidth = 1, alpha = 0.5, color = "darkred") +
  annotate(geom = "text", x = 1.4, y = (reference_extracellular_summary$ratio_extracellular_proteins[1])  , label = "FASTA", size = 2, color = "darkred") +
  theme(
    # panel.border = element_rect(colour = "blue", fill = NA, linetype = 2),
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.x  = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = -60, vjust = 0.9, hjust = 0, size = 6, family = "sans"),
    # axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5),
    # axis.text.x = element_text(angle = -60, vjust = 0.75, hjust = 0, size = 12, family = "sans"),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black", face = "plain", family = "sans", size = 14),
    axis.title = element_text(colour = "black", family = "sans", size = 14),
    axis.ticks = element_line(colour = "black"),
    title = element_text(size = 8, hjust = 0.5),
    strip.text = element_text(size= 12, family = "sans"),
    # legend at the bottom 6)
    legend.position = "right") +
  scale_y_continuous(limits = c(0, 0.2), breaks = c(seq(0, 0.2, 0.05)), expand = c(0,0), name = "ratio extracellular")

plot_extracellular_phosphoprotein_ratio_box

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_extracellular_phosphoprotein_ratio_box.png", plot = plot_extracellular_phosphoprotein_ratio_box, width = 4, height = 6, scale = 0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_extracellular_phosphoprotein_ratio_box.pdf", plot = plot_extracellular_phosphoprotein_ratio_box, width = 4, height = 6, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
annotated_df_sum_to_psite_cytoplasm <- ascore_fasta_protein_mod_loc_sum_psite %>% 
  left_join(locations_FASTA_df_ratios, by = "reference") %>% 
  mutate(
    cytoplasm = case_when(
      grepl("Cytoplasm", keywords) == TRUE ~ "cytoplasm",
      grepl("cytoplasm", gene_ontology_cellular_component) == TRUE ~ "cytoplasm",
      TRUE ~ "named_protein"))
  

#keyword: cytoplasm -------------------------------------
annotated_psite_summary_cytoplasm <- annotated_df_sum_to_psite_cytoplasm %>%
  distinct(condition, replicate, sample_id,reference, ref, protein_names, cytoplasm) %>% 
  group_by(condition, replicate, sample_id, cytoplasm) %>% 
  summarize(
    n_annotated_psites = n()  ) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(condition, replicate, sample_id), names_from = cytoplasm, values_from = n_annotated_psites) %>% 
  mutate(
    sum_psite_count = cytoplasm + named_protein,
    ratio_cytoplasm = cytoplasm / sum_psite_count )

#experimental protein summary -------------------------------------
annotated_protein_summary_cytoplasm<- annotated_df_sum_to_psite_cytoplasm %>%
  distinct(condition, replicate, sample_id,reference,  protein_names, cytoplasm) %>% 
  group_by(condition, replicate, sample_id, cytoplasm) %>% 
  summarize(
    n_annotated_proteins = n()  ) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(condition, replicate, sample_id), names_from = cytoplasm, values_from = n_annotated_proteins) %>% 
  mutate(
    total_proteins = cytoplasm + named_protein,
    ratio_cytoplasm_proteins = cytoplasm / total_proteins )


#reference Sros_GOs --------------------------------------
reference_cytoplasm_summary <- fasta_subcellular_location %>% 
  mutate(
    cytoplasm = case_when(
       grepl("Cytoplasm", keywords) == TRUE ~ "cytoplasm",
       grepl("cytoplasm", gene_ontology_cellular_component) == TRUE ~ "cytoplasm",
      TRUE ~ "named_protein")) %>%
  group_by(cytoplasm) %>% 
  summarize(
    n_annotated_proteins = n()) %>% 
  ungroup() %>%
  mutate(condition = "reference") %>% 
  pivot_wider(id_cols = condition, names_from = cytoplasm, values_from = n_annotated_proteins) %>% 
  mutate(
    total_proteins = cytoplasm + named_protein,
    ratio_cytoplasm_proteins = cytoplasm / total_proteins)

plot_cytoplasm_psite_ratio <- ggplot(data = annotated_protein_summary_cytoplasm) +
  geom_col(mapping = aes(x = sample_id, y = ratio_cytoplasm_proteins, fill = condition), color = "black", linewidth =1, alpha = 0.5) +
  geom_hline(yintercept = round((reference_cytoplasm_summary$ratio_cytoplasm_proteins[1]), 2), linetype = 2, linewidth = 1, alpha = 0.5, color = "darkred") +
  alexis_theme() +
  annotate(geom = "text", x = 34, y = (reference_cytoplasm_summary$ratio_cytoplasm_proteins[1]) + 0.02, label = "FASTA ratio", size = 3, color = "darkred") +
  alexis_theme() #+
  # scale_y_continuous(limits = c(0, 0.10), breaks = c(seq(0, 0.10, 0.05)), expand = c(0,0), name = "ratio cytoplasm")

plot_cytoplasm_psite_ratio

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_cytoplasm_psite_ratio.png", plot = plot_cytoplasm_psite_ratio, width = 16, height = 12, scale = 0.4)


#boxplot
plot_cytoplasm_phosphoprotein_ratio_box <- ggplot(data = annotated_protein_summary_cytoplasm %>% mutate(component = "cytoplasm")) +
  geom_boxplot(mapping = aes(x = component, y = ratio_cytoplasm_proteins), color = "black", linewidth =0.5, alpha = 0.5) +
  geom_jitter(mapping = aes(x = component, y = ratio_cytoplasm_proteins), color = "black", fill = "white", size = 0.5, alpha = 0.5, shape = 21, width = 0.25) +
  geom_hline(yintercept = round((reference_cytoplasm_summary$ratio_cytoplasm_proteins[1]), 2), linetype = 2, linewidth = 1, alpha = 0.5, color = "darkred") +
  annotate(geom = "text", x = 1.4, y = (reference_cytoplasm_summary$ratio_cytoplasm_proteins[1])  , label = "FASTA", size = 2, color = "darkred") +
  theme(
    # panel.border = element_rect(colour = "blue", fill = NA, linetype = 2),
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.x  = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = -60, vjust = 0.9, hjust = 0, size = 6, family = "sans"),
    # axis.text.x = element_text(angle = 0, vjust = 0.0, hjust = 0.5),
    # axis.text.x = element_text(angle = -60, vjust = 0.75, hjust = 0, size = 12, family = "sans"),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black", face = "plain", family = "sans", size = 14),
    axis.title = element_text(colour = "black", family = "sans", size = 14),
    axis.ticks = element_line(colour = "black"),
    title = element_text(size = 8, hjust = 0.5),
    strip.text = element_text(size= 12, family = "sans"),
    # legend at the bottom 6)
    legend.position = "right") +
  scale_y_continuous(limits = c(0, 0.65), breaks = c(seq(0, 0.60, 0.1)), expand = c(0,0), name = "ratio cytoplasm")

plot_cytoplasm_phosphoprotein_ratio_box

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_cytoplasm_phosphoprotein_ratio_box.png", plot = plot_cytoplasm_phosphoprotein_ratio_box, width = 4, height = 6, scale = 0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_cytoplasm_phosphoprotein_ratio_box.pdf", plot = plot_cytoplasm_phosphoprotein_ratio_box, width = 4, height = 6, scale = 0.4)



## --------------------------------------------------------------------------------------------------------------------------------------------------
annotated_df_locations <- ascore_fasta_protein_mod_loc_sum_psite %>% 
  full_join(locations_FASTA_df_ratios, by = "reference") %>% 
  mutate(
    membrane = case_when(
      grepl("Membrane", keywords) == TRUE ~ "membrane",
      grepl("membrane", gene_ontology_cellular_component) == TRUE ~ "membrane",
      TRUE ~ "named_protein"),
    transmembrane = case_when(
      grepl("Transmembrane", keywords) == TRUE ~ "transmembrane",
      grepl("transmembrane", gene_ontology_cellular_component) == TRUE ~ "transmembrane",
      TRUE ~ "named_protein"),
    extracellular = case_when(
      grepl("extracellular", keywords) == TRUE ~ "extracellular",
      grepl("extracellular", gene_ontology_cellular_component) == TRUE ~ "extracellular",
      TRUE ~ "named_protein"),
    cytoplasm = case_when(
      grepl("Cytoplasm", keywords) == TRUE ~ "cytoplasm",
      grepl("cytoplasm", gene_ontology_cellular_component) == TRUE ~ "cytoplasm",
      TRUE ~ "named_protein")) %>% 
  mutate(detected = case_when(
    is.na(condition) ~ "not_detected",
    !is.na(condition) ~ "detected"))


cytoplasm_detected <- annotated_df_locations %>% 
  filter(cytoplasm == "cytoplasm") %>% #keep only cytoplasmic proteins
  distinct(reference, detected) %>% 
  group_by(detected) %>% 
  summarize(
    n_proteins = n()  ) %>% 
  ungroup() %>% 
  mutate(location = "cytoplasm") %>% 
  pivot_wider(id_cols = location, names_from = detected, values_from = n_proteins)


membrane_detected <- annotated_df_locations %>% 
  filter(membrane == "membrane") %>% #keep only cytoplasmic proteins
  distinct(reference, detected) %>% 
  group_by(detected) %>% 
  summarize(
    n_proteins = n()  ) %>% 
  ungroup() %>% 
  mutate(location = "membrane") %>% 
  pivot_wider(id_cols = location, names_from = detected, values_from = n_proteins)

transmembrane_detected <- annotated_df_locations %>% 
  filter(transmembrane == "transmembrane") %>% #keep only cytoplasmic proteins
  distinct(reference, detected) %>% 
  group_by(detected) %>% 
  summarize(
    n_proteins = n()  ) %>% 
  ungroup() %>% 
  mutate(location = "transmembrane") %>% 
  pivot_wider(id_cols = location, names_from = detected, values_from = n_proteins)

ratios_location <- rbind(cytoplasm_detected, membrane_detected, transmembrane_detected) %>% 
  mutate(ratio_location = detected/ (detected+ not_detected))

###### plot observed out of all possible proteins detected per FASTA ###########
plot_location_ratios <- ggplot(data = ratios_location %>% filter(location != "membrane")) +
  geom_col(mapping = aes(x = location, y = ratio_location, fill = location), color = "black", linewidth = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#B8DE29FF", "#39568CFF","#287D8CFF" )) +
  alexis_theme() +
  theme(axis.text.x = element_text(angle = -60, vjust = 0.75, hjust =0)) +
  scale_y_continuous(limits = c(0, 0.3), breaks = c(seq(0, 0.3, 0.05)), name = "ratio detected", expand = c(0,0)) +
  xlab("")
  

plot_location_ratios

ggsave (filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/plot_location_ratios_forRachael.png", plot = plot_location_ratios, scale = 0.4, width = 5, height = 8)


## --------------------------------------------------------------------------------------------------------------------------------------------------
ascore_gene_noImputation <- ascore_final_gene

#likely need to modify mean and SD for imputation of global phospho relative to Y_phospho


## --------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(15)

#get SDs for condition and global to use in imputation
sd_ascore_final_gene <- ascore_gene_noImputation %>%
  filter(!is.nan((median_norm_intensity))) %>% 
  filter(!is.na(median_norm_intensity)) %>% 
  filter(!is.infinite(median_norm_intensity)) %>%
  mutate(all_sd = sd(median_norm_intensity, na.rm = TRUE)) %>% 
  # group_by(condition) %>% 
  # mutate(condition_sd = sd(median_norm_intensity, na.rm = TRUE)) %>% 
  ungroup() %>%
  group_by(condition, gene_ref) %>% 
  mutate(
    condition_gene_ref_sd = sd(median_norm_intensity, na.rm = TRUE),
    condition_gene_ref_mean_intensity = mean(median_norm_intensity, na.rm = TRUE)) %>%
  ungroup() %>%
  
  # group_by(condition) %>% #going with global average standard devation of 0.63, and not condition specific 0.61 - 0.67 for imputation.
  
  #get average within condition + psite SD
  mutate(avg_condition_gene_ref_sd = mean(condition_gene_ref_sd, na.rm = TRUE)) %>% 
  
  # ungroup() %>% #used for assessment of condition specific average standard deviations.
  distinct(condition, gene_ref, all_sd, condition_gene_ref_sd, avg_condition_gene_ref_sd, condition_gene_ref_mean_intensity)


## --------------------------------------------------------------------------------------------------------------------------------------------------
#visualize distributions of observations per p-site globally
  ##---df
observation_filter <- ascore_gene_noImputation %>%
  mutate(gene_ref = str_replace_all(gene_ref, " ", "_")) %>%
  filter(!is.na(median_norm_intensity)) %>% 
  filter(!is.infinite(median_norm_intensity)) %>% 
  group_by(gene_ref) %>% 
  mutate(
    n_observations = n())

  ##---histogram of observations per p-site globally
plot_observations <- ggplot(data = observation_filter) +
  geom_bar(mapping = aes(x = n_observations)) +
  geom_vline(xintercept = 3, linetype = 2, color = "darkblue") +
  alexis_theme()
plot_observations


## --------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(15)
#imputation using normal distributions happens here:
matrix_ascore_gene_imputed_wSD <- ascore_gene_noImputation %>%
  pivot_wider(id_cols = c("gene_ref", "condition"),
              names_from = replicate,
              values_from =  median_norm_intensity,
              values_fill = NA) %>%
  select(gene_ref, condition, rep1, rep2, rep3, rep4, rep5, rep6) %>% 
  pivot_longer(cols = c(rep1:rep6), names_to = "replicate") %>% 
  rename(median_norm_intensity = value) %>% 
  mutate(ref_condition = paste(gene_ref, condition, sep = " ")) %>% 
  filter(ref_condition != "H1 S2 EGF3min") %>% #remove the row with two intensities per psite, no protein reference. why??
  mutate(median_norm_intensity = as.character(median_norm_intensity),
         median_norm_intensity = str_replace_all(median_norm_intensity, "NULL", "0"),
         median_norm_intensity = str_replace_all(median_norm_intensity, "-Inf", "0"),
         median_norm_intensity = str_replace_all(median_norm_intensity, "0", "NA"),
         median_norm_intensity = as.numeric(median_norm_intensity)) %>% 

  left_join(y = sd_ascore_final_gene, by = c("condition", "gene_ref" )) %>%
  group_by(condition, gene_ref) %>% 
  
  mutate(
    num_missing_obs = sum(is.na(median_norm_intensity))) %>% 
  
  ungroup() %>%
  group_by(condition, gene_ref, replicate) %>% 
  mutate(
    imputed_intensity = case_when(
      !is.na(median_norm_intensity) ~ median_norm_intensity,
      
      #imputation of partially missing replicates, sd = condition SD, mean = mean of measured reps
      is.na(median_norm_intensity) & num_missing_obs < 6 ~ sample(rnorm(n = 6,
                                                                        mean = condition_gene_ref_mean_intensity,
                                                                        sd = avg_condition_gene_ref_sd), size = 1, prob = NULL),
      #imputation of psite missing in all replicates SD1= 1, mean = 10
      is.na(median_norm_intensity) & num_missing_obs  == 6 ~ sample(x = rnorm(n = 6, mean = 10, sd = 1),
                                                                    size =1, replace = FALSE, prob = NULL))) %>%
  ungroup() %>% 
  mutate(
      
      # 6 here is equivalent to adding in 'num_missing_obs' to choose vector size from normal distribution
      
      
      #add in column to backcalculate imputed measurements
      imputed_replicate = case_when(
        num_missing_obs > 0 ~ "imputed",
        num_missing_obs == 0 ~ "observed" ),
      imputed_gt_3_reps = case_when(
        num_missing_obs > 3 ~ "imputed",
        num_missing_obs <= 3 ~ "observed at least 3 reps" ))
      #can use above dataframe to look at distributions of psite missingness by plotting histogram of num_missing observations.
  


##VERY MESSY
  #separate dataframe for pivoting. seems to be used in fold change analysis below.
  matrix_ascore_gene_imputed_wSD_wider <- matrix_ascore_gene_imputed_wSD %>%
    ungroup() %>%
    distinct(condition, gene_ref, replicate,  imputed_replicate, imputed_gt_3_reps, imputed_intensity) %>%
    mutate(gene_ref = str_replace_all(gene_ref, pattern = " ", replacement = "_"),
           sample_id = paste(condition, replicate, sep = "_")) %>%

  #now identify psites missing completely at one time point but measured in others.
  pivot_wider(id_cols = c("gene_ref"),
              names_from = sample_id,
              values_from = imputed_intensity,
              values_fill = NA) %>%
    pivot_longer(cols = contains("EGF"), names_to = "sample_id") %>%
    rename(imputed_intensity = value) %>%
    group_by(gene_ref, sample_id) %>%  #to allow for different value sampling during imputation in next step
    mutate(imputed_intensity = case_when(
      !is.na(imputed_intensity) ~ imputed_intensity,
      is.na(imputed_intensity) ~ sample(x = rnorm(n = 24, mean = 10, sd = 1), size =1, replace = FALSE, prob = NULL))) %>%
    ungroup() %>%
    #chose n = 24 to account for some cases where all 4 conditions X 6 reps were not measured. sampling to size = 1 means to keep one value.
    pivot_wider(id_cols = c("gene_ref"),
              names_from = sample_id,
              values_from = imputed_intensity,
              values_fill = NA)





  # values_fill = sample(x = rnorm(n = 6, mean = 6, sd = 2.9), size =1, replace = FALSE, prob = NULL))
              # values_fill = replace_na(sample(x = rnorm(n = 6, mean = 9, sd = 3), size =1, replace = FALSE, prob = NULL)))
  
  #likely do not need the values fill argument here because I handled sampling along a normal distribution above.


## --------------------------------------------------------------------------------------------------------------------------------------------------
#to join with dataframe below for completely missing datapoints
data_completeness_wider_to_rbind <- matrix_ascore_gene_imputed_wSD %>% 
  distinct(gene_ref, condition, replicate, median_norm_intensity, imputed_intensity, imputed_gt_3_reps, imputed_replicate) %>% 
  mutate(gene_ref = str_replace_all(gene_ref, " ", "_"),
         sample_id = paste(condition, replicate, sep = "_")) %>%
  pivot_wider(id_cols = c("gene_ref"),
              names_from = sample_id,
              values_from = imputed_intensity,
              values_fill = NA) %>% 
    pivot_longer(cols = contains("EGF"), names_to = "sample_id") %>% 
    rename(imputed_intensity = value) %>% 
  mutate(imputed_replicate = case_when(
    is.na(imputed_intensity) ~ "imputed",
    !is.na(imputed_intensity) ~ "maybe observed to discard")) %>% 
  filter(imputed_replicate == "imputed") %>% 
  
  separate(sample_id, into = c("condition", "replicate"), sep = "_") %>% 
  select(condition, replicate, gene_ref, imputed_replicate)


## --------------------------------------------------------------------------------------------------------------------------------------------------
##CATEGORIZE IMPUTATION PER PSITE AND CONDITION
data_completeness <- matrix_ascore_gene_imputed_wSD %>% 
  distinct(gene_ref, condition, replicate, median_norm_intensity, imputed_intensity, imputed_gt_3_reps, imputed_replicate)


data_completeness_summary <- matrix_ascore_gene_imputed_wSD %>% 
  distinct(gene_ref, condition, replicate, median_norm_intensity, imputed_intensity, imputed_gt_3_reps, imputed_replicate) %>% 
  mutate(gene_ref = str_replace_all(gene_ref, " ", "_"),
         sample_id = paste(condition, replicate, sep = "_")) %>%
  distinct(gene_ref, condition, replicate, imputed_replicate) %>% 
  rbind(data_completeness_wider_to_rbind) %>% 
  group_by(gene_ref, condition, imputed_replicate) %>%
  
  summarize(n_observations = n()) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(gene_ref, condition), names_from = imputed_replicate, values_from = n_observations, values_fill = 0) %>% 
  mutate(percent_observed = observed / 6 * 100)


## --------------------------------------------------------------------------------------------------------------------------------------------------
#plots -----------------------------------------
# imputed_data_analysis_intensity <- ggplot(data = data_completeness) +
#   geom_boxplot(mapping = aes(x = imputed_gt_3_reps, y = imputed_intensity)) +
#   alexis_theme()
# 
# imputed_data_analysis_intensity
# 
# 
# 
# imputed_data_analysis_count <- ggplot(data = data_completeness %>% distinct(gene_ref, condition, imputed_gt_3_reps)) +
#   geom_bar(mapping = aes(x = imputed_gt_3_reps)) +
#   alexis_theme()
# 
# imputed_data_analysis_count


## --------------------------------------------------------------------------------------------------------------------------------------------------
longer_imputed_intensities_wSD <- matrix_ascore_gene_imputed_wSD_wider %>% 
  pivot_longer(cols = c(EGF3min_rep4:EGF15min_rep3)) %>% #updated column ordering, should reorder then select to hard code better
  rename(imputed_intensities = value) %>% 
  separate(col = name, into = c("condition", "replicate"), sep = "_")



longer_imputed_intensities_wSD_for_diff <- longer_imputed_intensities_wSD %>%
  mutate(sample_id = paste(condition, replicate, sep = "_")) %>% 
  assign_missingness(
    sample = sample_id,
    condition = condition, 
    grouping = gene_ref,
    intensity = imputed_intensities, 
    ref_condition = "all" )


foldchange_imputed_wSD <- longer_imputed_intensities_wSD_for_diff %>%
  calculate_diff_abundance(
    sample = sample_id,
    condition = condition,
    grouping = gene_ref,
    intensity_log2 = imputed_intensities,
    missingness = missingness,
    comparison = comparison) %>% 
  mutate(
    diff = case_when(
      comparison == "EGF0min_vs_EGF1min" ~ -diff,
      comparison == "EGF0min_vs_EGF3min" ~ -diff,
      comparison == "EGF0min_vs_EGF5min" ~ -diff,
      comparison == "EGF0min_vs_EGF15min" ~ -diff,
      TRUE ~ diff),
    comparison = case_when(
      comparison == "EGF0min_vs_EGF1min" ~ "EGF1min_vs_EGF0min",
      comparison == "EGF0min_vs_EGF3min" ~ "EGF3min_vs_EGF0min",
      comparison == "EGF0min_vs_EGF5min" ~ "EGF5min_vs_EGF0min",
      comparison == "EGF0min_vs_EGF15min" ~ "EGF15min_vs_EGF0min",
      TRUE ~ comparison))
  


##can add clusters here too!



## --------------------------------------------------------------------------------------------------------------------------------------------------
net_regulated_psites_wSD <- foldchange_imputed_wSD %>%
                filter (comparison %in% c("EGF1min_vs_EGF0min", "EGF3min_vs_EGF0min", "EGF5min_vs_EGF0min", "EGF15min_vs_EGF0min")) %>% 
                distinct(comparison, diff, adj_pval, gene_ref) %>%
                group_by(comparison) %>%
                mutate(total_psites = n_distinct(gene_ref)) %>%
                ungroup() %>%
                mutate(
                  regulated_psite = case_when(
                    adj_pval <= 0.05 & abs(diff) >= 1 ~ "regulated",
                    TRUE ~ "not regulated")) %>% 
  filter(regulated_psite == "regulated") %>% 
  distinct(gene_ref) %>% 
  summarize(
    n_regulated_psites_between_min_1_duration = n() )


## --------------------------------------------------------------------------------------------------------------------------------------------------
percent_psite_regulated_imputed <- foldchange_imputed_wSD %>%
  filter (comparison %in% c("EGF1min_vs_EGF0min", "EGF3min_vs_EGF0min", "EGF5min_vs_EGF0min", "EGF15min_vs_EGF0min")) %>% 
  distinct(comparison, diff, adj_pval, gene_ref) %>% 
  group_by(comparison) %>% 
  mutate(total_psites = n_distinct(gene_ref)) %>%
  ungroup() %>%
  mutate(
    regulated_psite = case_when(
      adj_pval <= 0.05 & abs(diff) >= 1 ~ "regulated",
      TRUE ~ "not_regulated"),
    regulated_psite_dir = case_when(
      regulated_psite == "regulated" & diff >= 1 ~ "up",
      regulated_psite == "regulated" & diff <= -1 ~ "down",
      regulated_psite == "not_regulated" ~ "not_regulated")) %>%
  group_by(comparison, regulated_psite_dir) %>% 
  mutate(n_regulated_psites = n()) %>% 
  ungroup() %>% 
  distinct(comparison, total_psites, regulated_psite, n_regulated_psites, regulated_psite_dir) %>% 
  pivot_wider(id_cols= c(comparison, total_psites), names_from = regulated_psite_dir, values_from = n_regulated_psites) %>% 
  group_by(comparison) %>% 
  mutate(percent_upregulated = up / total_psites*100,
         percent_downregulated = down / total_psites*100) %>% 
  ungroup() %>% 
  mutate(
    total_regulated_psites = up + down,
    pct_up_vs_total_regulated = up / total_regulated_psites *100,
    pct_down_vs_total_regulated = down / total_regulated_psites*100,
    percent_regulated = total_regulated_psites / total_psites*100,
    percent_not_regulated = not_regulated / total_psites*100) %>% 
  
  #annotation positions for volcano plot
  mutate(
    x_downreg_position = -12,
    x_upreg_position = 12,
    y_position = 6)
  

percent_psite_regulated_imputed


## --------------------------------------------------------------------------------------------------------------------------------------------------
psite_diff_abundance_gene_2xDB_imputed_forvolcano <- foldchange_imputed_wSD %>%
  # filter(gene_ref %in% gt3_obs_for_FC) %>% 
  separate(gene_ref, into= c("gene", "mod_residue"), sep = "_", remove = FALSE) %>% 
  rename(psite = mod_residue) %>% 
  mutate(mod_res = str_sub(psite, end = 1L),
         mod_residue = as.numeric(str_sub(psite, start = 2L))) %>% 
  filter(abs(diff) > 1) %>%
  filter (adj_pval <= 0.05) %>%
  
  mutate(my_data = TRUE) %>% 
  full_join(y = PSP_vs_PTMSigDB_EGF_matches %>% mutate(DB2x = TRUE), by = c("gene", "mod_residue", "mod_res") ) %>% 
  mutate(
    overlap_my_data = case_when(
      my_data == TRUE & DB2x ==TRUE ~ "both",
      my_data == TRUE & is.na(DB2x) ~ "this\nstudy",
      is.na(my_data) & DB2x ==TRUE ~ "DBs\nonly"),
    overlap_my_data = fct_relevel(overlap_my_data, "both", "this\nstudy", "DBs\nonly"))

#counts for annotating plot
sums_psites_overlap_DBs_my_data_imputed_forvolcano <- psite_diff_abundance_gene_2xDB_imputed_forvolcano %>%
  # filter(mod_res == "Y") %>% 
  group_by(overlap_my_data) %>%
  mutate(n_sites = n()) %>%
  ungroup() %>% 
  distinct(overlap_my_data, n_sites) %>% 
  mutate(n_sites_char = as.character(n_sites),
         overlap_num = case_when(
           overlap_my_data == "both" ~ 1,
           overlap_my_data == "this\nstudy" ~ 2,
           overlap_my_data == "DBs\nonly" ~ 3))



## --------------------------------------------------------------------------------------------------------------------------------------------------
regulated_pSTYsites_wimputation_wSD <- foldchange_imputed_wSD %>%
  filter (comparison %in% c("EGF1min_vs_EGF0min", "EGF3min_vs_EGF0min", "EGF5min_vs_EGF0min", "EGF15min_vs_EGF0min")) %>% 
  distinct(comparison, diff, adj_pval, gene_ref) %>%
  filter(adj_pval < 0.05) %>% 
  filter(diff < -1 | diff > 1)


## --------------------------------------------------------------------------------------------------------------------------------------------------
#combine all three databases of genes annotated to be EGF responsive ------------
WP437 <- c('NOS3','COX2','SHC1','STAT5A','RALA','NEDD4','ELK4','MAPK1','MAPK8','ARHGEF1','STMN1','STAT1','EPS15L1','SH3GL3','EGF','INPP5D','PTPN5','ATF1','GRB10','ERRFI1','GRIM19','MAP3K1','PTK6','RPS6KA2','CSK','PCNA','KRAS','LIMK2','RASA1','ARF6','IQGAP1','EGFR','SP1','BRAF','DOK2','PTPRR','CRK','SH2D2A','MAPK9','RPS6KA3','RIN1','CAV2','E2F1','VAV2','SOS2','PTPN12','MTOR','NCK2','MAP3K3','RALGDS','VAV3','CAMK2A','STAT3','MAPK3','SOS1','GAB2','RAP1A','JAK1','CDC42','NCOA3','ATXN2','RAB5A','RPS6KB1','MAPK14','ASAP1','SYNJ1','MAP2K1','RAF1','EPN1','NCK1','STAMBP','CAV1','STAM2','HGS','MAP2K2','CREB1','PLCG1','ERBB2','SPRY2','INPPL1','EIF4EBP1','MAP4K1','RICTOR','USP6NL','RPS6KA5','PTK2B','PIK3R2','PLD2','EPS8','RAC1','ABI1','PLSCR1','PIK3R1','MEF2D','GAB1','PLCE1','PLD1','DNM1','JUN','PIAS3','SH3GL2','RPS6KA1','ROCK1','MAP3K4','PRKCD','AP2B','AKT1','FOXO1','VAV1','PRKCZ','RALB','FOS','PRKCI','PRKCA','SRC','ITCH','ABL1','RALBP1','MAPK7','REPS2','JUND','GRB2','EPS15','USP8','AP2A1','HRAS','PDPK1','IQSEC1','PEBP1','SH3KBP1','TNK2','PTPN11','AP2S1','PRKCB','FOXO4','FOSB','PTEN','NEDD8','STAM','GJA1','CRKL','JAK2','MAP2K5','BCAR1','MAP3K2','ELK1','STAT5B','PIK3C2B','STXBP1','CBLB','CBL','CBLC','PAK1','CFL1','PTK2','AP2M1','MEF2A','PXN','MEF2C')

#turn wiki pathways into a dataframe
WP437_gene_df <- data.frame("gene" = WP437, "WP437" = TRUE)

#combine genes from all four database sources into a single dataframe
gene_3xDB_overlap <-  PSP_vs_PTMSigDB_EGF_matches %>%
              # filter(mod_res == "Y") %>% # filter for only proteins with pY observed when possible. WP437 lacks site info
              # mutate(DB1 = "PSP|PTMSigDB|PpDIA") %>% 
              full_join(y = WP437_gene_df, by = "gene") %>% 
              distinct(gene, PSP, PTM_SEA,  WP437) %>% 
  pivot_longer(cols = c(PSP, PTM_SEA,  WP437), names_to = "DB") %>%
  rename(present_in_DB = value) %>%
  filter(present_in_DB == TRUE) %>% 
  distinct(gene, present_in_DB)


#join with my data ----------------------------------------------------
gene_diff_abundance_3xDB <- foldchange_imputed_wSD %>%
  separate(gene_ref, into = c("gene", "mod_residue"), sep = "_", remove = FALSE) %>% 
  mutate(mod_res = str_sub(mod_residue, end = 1L)) %>% 
  # rename(mod_residue = modres_position) %>%
  # filter(mod_res == "Y") %>% 
  filter(abs(diff) > 1) %>%
  filter (adj_pval <= 0.05) %>%
   
  mutate(my_data = TRUE) %>%
  full_join(y = gene_3xDB_overlap %>% mutate(DB3x = TRUE), by = "gene") %>% 
  # full_join(y = PSP_vs_PTMSigDB_EGF_matches %>%
  #             filter(mod_res == "Y") %>% # filter for only proteins with pY observed when possible. WP437 lacks site info
  #             full_join(y = WP437_gene_df, by = "gene") %>% 
  #             distinct(gene) %>%
  #             mutate(DB4x = TRUE), by = c("gene") ) %>% 
  mutate(
    overlap_my_data = case_when(
      my_data == TRUE & DB3x ==TRUE ~ "both",
      my_data == TRUE & is.na(DB3x) ~ "this\nstudy",
      is.na(my_data) & DB3x ==TRUE ~ "DBs\nonly"),
    overlap_my_data = fct_relevel(overlap_my_data, "both", "this\nstudy", "DBs\nonly"))




## --------------------------------------------------------------------------------------------------------------------------------------------------
# # set.seed(15)
# # #get SDs for condition and global to use in imputation
mean_intensities_noSD <- ascore_gene_noImputation %>%
  filter(!is.nan((median_norm_intensity))) %>%
  filter(!is.na(median_norm_intensity)) %>%
  filter(!is.infinite(median_norm_intensity)) %>%
  # mutate(all_sd = sd(median_norm_intensity, na.rm = TRUE)) %>%
  # group_by(condition) %>%
  # mutate(condition_sd = sd(median_norm_intensity, na.rm = TRUE)) %>%
  # ungroup() %>%
  group_by(condition, gene_ref) %>%
  mutate(condition_gene_ref_mean_intensity = mean(median_norm_intensity, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(condition, gene_ref, condition_gene_ref_mean_intensity)


## --------------------------------------------------------------------------------------------------------------------------------------------------
#imputation using normal distributions happens here:
matrix_ascore_gene_imputed_t_2_again_noSD <- ascore_gene_noImputation %>%
  pivot_wider(id_cols = c("gene_ref", "condition"),
              names_from = replicate,
              values_from =  median_norm_intensity,
              values_fill = NA) %>%
  select(gene_ref, condition, rep1, rep2, rep3, rep4, rep5, rep6) %>% 
  pivot_longer(cols = c(rep1:rep6), names_to = "replicate") %>% 
  rename(median_norm_intensity = value) %>% 
  mutate(ref_condition = paste(gene_ref, condition, sep = " ")) %>% 
  filter(ref_condition != "H1 S2 EGF3min") %>% #remove the row with two intensities per psite, no protein reference. why??
  
  mutate(median_norm_intensity = as.character(median_norm_intensity),
         median_norm_intensity = str_replace_all(median_norm_intensity, "NULL", "0"),
         median_norm_intensity = str_replace_all(median_norm_intensity, "-Inf", "0"),
         median_norm_intensity = str_replace_all(median_norm_intensity, "0", "NA"),
         median_norm_intensity = as.numeric(median_norm_intensity)) %>%
  
  
  # select(gene_ref, condition, rep1, rep2, rep3, rep4, rep5, rep6) %>% 
  # pivot_longer(cols = c(rep1:rep6), names_to = "replicate") %>% 
  # rename(median_norm_intensity = value) %>% 
  # mutate(ref_condition = paste(gene_ref, condition, sep = " "),
  #        median_norm_intensity = case_when(is.infinite(median_norm_intensity) ~ NA,
  #                                          is.na(median_norm_intensity) ~ NA,
  #                                          TRUE ~ median_norm_intensity)) %>%
  
  
  left_join(y = mean_intensities_noSD, by = c("condition", "gene_ref" )) %>%
  group_by(condition, gene_ref) %>% 
  
  mutate(
    num_missing_obs = sum(is.na(median_norm_intensity))) %>% 
  
  ungroup() %>%
  group_by(condition, gene_ref, replicate) %>% 
  mutate(imputed_intensity = case_when(
    !is.na(median_norm_intensity) ~median_norm_intensity,
    is.na(median_norm_intensity) & num_missing_obs < 6 ~  condition_gene_ref_mean_intensity,
    is.na(median_norm_intensity) & num_missing_obs  == 6 ~ 10),
    imputed_replicate = case_when(
      is.na(median_norm_intensity) ~ "imputed",
      !is.na(median_norm_intensity) ~ "observed")) %>%  #going with bottom 5% of signal still seems too high; instead use min signal 
  ungroup() %>% 
  mutate(
      
      # 6 here is equivalent to adding in 'num_missing_obs' to choose vector size from normal distribution
      
      
      #add in column to backcalculate imputed measurements
      imputed_gt_3_reps = case_when(
        num_missing_obs > 3 ~ "imputed",
        num_missing_obs <= 3 ~ "observed at least 3 reps" ))
      #can use above dataframe to look at distributions of psite missingness by plotting histogram of num_missing observations.
  



##PIVOT WIDER TO COMPARE TRENDS ACROSS ALL TIME POINTS AND REPS
##FEEDS INTO ZSCORE AND CLUSTERING
  #separate dataframe for pivoting. 
  matrix_ascore_gene_imputed_t_3_noSD <- matrix_ascore_gene_imputed_t_2_again_noSD %>%
    ungroup() %>% 
    distinct(condition, gene_ref, replicate,  imputed_replicate, imputed_gt_3_reps, imputed_intensity) %>% 
    mutate(gene_ref = str_replace_all(gene_ref, pattern = " ", replacement = "_"),
           sample_id = paste(condition, replicate, sep = "_")) %>% 
    
  #now identify psites missing completely at one time point but measured in others.
  pivot_wider(id_cols = c("gene_ref"),
              names_from = sample_id,
              values_from = imputed_intensity,
              values_fill = NA) %>% 
    pivot_longer(cols = contains("EGF"), names_to = "sample_id") %>% 
    rename(imputed_intensity = value) %>% 
    group_by(gene_ref, sample_id) %>%  #to allow for different value sampling during imputation in next step
    mutate(imputed_intensity = case_when(
      !is.na(imputed_intensity) ~ imputed_intensity,
      is.na(imputed_intensity) ~ 10)) %>% #impute to minimum signal for conditions lacking observations at all
    ungroup() %>% 
    #chose n = 24 to account for some cases where all 4 conditions X 6 reps were not measured. sampling to size = 1 means to keep one value.
    pivot_wider(id_cols = c("gene_ref"),
              names_from = sample_id,
              values_from = imputed_intensity,
              values_fill = NA) %>% 
    #KEEP ONLY SIGNIFICANTLY CHANGING pSTY SITES PER DIFF ABUNDANCE WITH SD
    filter(gene_ref %in% ((regulated_pSTYsites_wimputation_wSD %>% distinct(gene_ref))$gene_ref)) #from section 16
    
              
  
  
  
  
  
  # values_fill = sample(x = rnorm(n = 6, mean = 6, sd = 2.9), size =1, replace = FALSE, prob = NULL))
              # values_fill = replace_na(sample(x = rnorm(n = 6, mean = 9, sd = 3), size =1, replace = FALSE, prob = NULL)))
  
  #likely do not need the values fill argument here because I handled sampling along a normal distribution above.


## --------------------------------------------------------------------------------------------------------------------------------------------------


#to join with dataframe below for completely missing datapoints
data_completeness_wider_to_rbind_noSD <- matrix_ascore_gene_imputed_t_2_again_noSD %>% 
  distinct(gene_ref, condition, replicate, median_norm_intensity, imputed_intensity, imputed_gt_3_reps, imputed_replicate) %>% 
  mutate(gene_ref = str_replace_all(gene_ref, " ", "_"),
         sample_id = paste(condition, replicate, sep = "_")) %>%
  pivot_wider(id_cols = c("gene_ref"),
              names_from = sample_id,
              values_from = imputed_intensity,
              values_fill = NA) %>% 
    pivot_longer(cols = contains("EGF"), names_to = "sample_id") %>% 
    rename(imputed_intensity = value) %>% 
  mutate(imputed_replicate = case_when(
    is.na(imputed_intensity) ~ "imputed",
    !is.na(imputed_intensity) ~ "maybe observed to discard")) %>% 
  filter(imputed_replicate == "imputed") %>% 
  
  separate(sample_id, into = c("condition", "replicate"), sep = "_") %>% 
  select(condition, replicate, gene_ref, imputed_replicate)


##CATEGORIZE IMPUTATION PER PSITE AND CONDITION
data_completeness_noSD <- matrix_ascore_gene_imputed_t_2_again_noSD %>% 
  distinct(gene_ref, condition, replicate, median_norm_intensity, imputed_intensity, imputed_gt_3_reps, imputed_replicate)


data_completeness_summary_noSD <- matrix_ascore_gene_imputed_t_2_again_noSD %>% 
  distinct(gene_ref, condition, replicate, median_norm_intensity, imputed_intensity, imputed_gt_3_reps, imputed_replicate) %>% 
  mutate(gene_ref = str_replace_all(gene_ref, " ", "_"),
         sample_id = paste(condition, replicate, sep = "_")) %>%
  distinct(gene_ref, condition, replicate, imputed_replicate) %>% 
  rbind(data_completeness_wider_to_rbind) %>% 
  group_by(gene_ref, condition, imputed_replicate) %>%
  
  summarize(n_observations = n()) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(gene_ref, condition), names_from = imputed_replicate, values_from = n_observations, values_fill = 0) %>% 
  mutate(percent_observed = observed / 6 * 100) 






#plots -----------------------------------------
imputed_data_analysis_intensity_noSD <- ggplot(data = data_completeness_noSD) +
  geom_boxplot(mapping = aes(x = imputed_gt_3_reps, y = imputed_intensity)) +
  alexis_theme()

imputed_data_analysis_intensity_noSD



imputed_data_analysis_count_noSD <- ggplot(data = data_completeness_noSD %>% distinct(gene_ref, condition, imputed_gt_3_reps)) +
  geom_bar(mapping = aes(x = imputed_gt_3_reps)) +
  alexis_theme()

imputed_data_analysis_count_noSD


## --------------------------------------------------------------------------------------------------------------------------------------------------
#9 is added if no p-sites were detected in all 6 replicates for a given time point for a specific psite.

colnames_for_zscore_matrix_noSD <- colnames(matrix_ascore_gene_imputed_t_3_noSD %>% select(-gene_ref))


cmeans_clustering_input_matrix_reps_noSD <- as.matrix(matrix_ascore_gene_imputed_t_3_noSD[,2:31])

rownames(cmeans_clustering_input_matrix_reps_noSD) <- as.matrix(matrix_ascore_gene_imputed_t_3_noSD[,1])

zscore_cmeans_clustering_matrix_reps_noSD <- apply(cmeans_clustering_input_matrix_reps_noSD, 1, scale) %>% t()
colnames(zscore_cmeans_clustering_matrix_reps_noSD) <- colnames_for_zscore_matrix_noSD

zscore_cmeans_clustering_matrix_reps_noSD[is.nan(zscore_cmeans_clustering_matrix_reps_noSD)] <- 0


## --------------------------------------------------------------------------------------------------------------------------------------------------
# #3 clusters ------------------------------
kmm3 <- kmeans(zscore_cmeans_clustering_matrix_reps_noSD, 3, nstart = 50, iter.max = 15)
# kmm3
# # Within cluster sum of squares by cluster (pY):
# # [1] 14178.881  7229.269  7366.664
# #  (between_SS / total_SS =  39.9 %)

# Within cluster sum of squares by cluster (pSTY):
# [1] 55865.62 54269.62 46332.60
#  (between_SS / total_SS =  35.2 %)
# 
# 
# #4 clusters ------------------------------
kmm4 <- kmeans(zscore_cmeans_clustering_matrix_reps_noSD, 4, nstart = 50, iter.max = 15)
# kmm4
# # Within cluster sum of squares by cluster: (pY)
# # [1]  3700.014 12782.701  6553.047  2956.340
# # (between_SS / total_SS =  45.7 %)

# Within cluster sum of squares by cluster: (pSTY)
# [1] 37713.80 24402.42 44951.31 25177.06
#  (between_SS / total_SS =  45.2 %)
# 
# #4 clusters is decent. Even 5 could work but 5 or 6 clusters begins showing slowing of explaning data variability.
# 
# #5 clusters ------------------------------
kmm5 <- kmeans(zscore_cmeans_clustering_matrix_reps_noSD, 5, nstart = 50, iter.max = 15)
# kmm5
# # Within cluster sum of squares by cluster(pY):
# # [1] 3258.073 4052.241 6157.704 8575.211 2217.156
# #  (between_SS / total_SS =  49.3 %)

# Within cluster sum of squares by cluster(pSTY):
# [1] 23158.46 20826.68 22589.17 20854.32 23747.82
#  (between_SS / total_SS =  54.0 %)
# 
# #6 clusters ------------------------------
kmm6 <- kmeans(zscore_cmeans_clustering_matrix_reps_noSD, 6, nstart = 50, iter.max = 15)
# kmm6
# # Within cluster sum of squares by cluster (pY):
# # [1] 3047.057 6414.381 6253.864 2119.976 2044.199 3026.539
# #  (between_SS / total_SS =  52.1 %)

# Within cluster sum of squares by cluster (pSTY):
# [1] 20788.498  4673.733  9967.534 23230.378 21256.672 19738.241
#  (between_SS / total_SS =  58.7 %)


## --------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(15)
fcm_result_reps <- cmeans(zscore_cmeans_clustering_matrix_reps_noSD, centers = 4, m = 1.2)

# fcm_result <- cmeans(cmeans_clustering_input_matrix, centers = 6, m = 2)

# print(fcm_result_reps)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Extract the cluster membership matrix
membership_values <- fcm_result_reps$membership

#show the membership values with simple print
# print(membership_values)

#optionally, convert the matrix to a dataframe for easier viewing
membership_df <- as.data.frame(membership_values)

#optionally, I can add clusters to the dataframe
membership_df$cluster <- apply(membership_values, 1, which.max) #get the cluster with highest membership

membership_df$cluster <- as.factor(membership_df$cluster) #convert to double or factor if desired

#update row and column names for dataframe
membership_df <- membership_df %>% 
  rownames_to_column(var = "ref")
colnames(membership_df) <- c("ref", "membership_cluster_1","membership_cluster_2","membership_cluster_3","membership_cluster_4", "cluster")

#pivot dataframe to have one correlation per psite
membership_df_long <- membership_df %>%
  pivot_longer(cols= c("membership_cluster_1","membership_cluster_2","membership_cluster_3","membership_cluster_4"), names_to = "potential_clusters", values_to = "membership_correlation")

#filter to keep only most correlated cluster per p-site
membership_df_long_less <- membership_df_long %>%
  mutate(potential_clusters = as.factor(str_sub(potential_clusters, start = -1L))) %>% 
  group_by(ref) %>% 
  filter(cluster == potential_clusters) %>% ungroup() %>% 
  select(-potential_clusters)

#now I can join the cluster membership correlation (membership_df_long_less) to the graph below for better coloration!


## --------------------------------------------------------------------------------------------------------------------------------------------------
data_clustered_reps <- data.frame(zscore_cmeans_clustering_matrix_reps_noSD)
# data_clustered <- data.frame(cmeans_clustering_input_matrix)

#assign to most correlated cluster
data_clustered_reps$cluster <- as.factor(apply(fcm_result_reps$membership, 1, which.max))

data_clustered_reps$gene_ref <- rownames(data_clustered_reps)

data_clustered_reps_df <- as.data.frame(data_clustered_reps) %>%
  pivot_longer(cols = c(EGF3min_rep1:EGF15min_rep6), names_to = c("timepoint")) %>%
  separate(timepoint, into =c("condition", "replicate"), sep = "_") %>%

  mutate(minute = case_when(
    condition == "EGF15min" ~ 15,
    condition == "EGF5min" ~ 5,
    condition == "EGF3min" ~ 3,
    condition == "EGF1min" ~ 1,
    condition == "EGF0min" ~ 0)) %>%
  select(-condition) %>%
  rename(mean_condition_intensity = value) %>%
  # group_by(gene_ref) %>%
  # mutate(
  #   ref_time_0 = mean_condition_intensity[minute == 0],
  #   fold_change_time0 = mean_condition_intensity - ref_time_0) %>%
  # ungroup() %>%
  # group_by(cluster) %>%
  separate(gene_ref, into = c("gene", "psite" ), sep = "_", remove = FALSE) %>%
  rename(cluster_unordered = cluster) %>% #rename to allow reordering safely
  mutate(
    cluster = case_when(
      cluster_unordered == 1 ~ 3,
      cluster_unordered == 2 ~ 2,
      cluster_unordered == 3 ~ 1,
      cluster_unordered == 4 ~ 4))


 data_clustered_reps_df <- data_clustered_reps_df %>% 
  left_join(y = membership_df_long_less %>% rename(gene_ref = ref) %>% rename(cluster_unordered = cluster), by = c("cluster_unordered", "gene_ref")) %>% 
  mutate(membership_correlation = as.factor(membership_correlation))


## --------------------------------------------------------------------------------------------------------------------------------------------------

cluster_membership_plot1_reps <- ggplot((data_clustered_reps_df ), aes(x = minute, y = mean_condition_intensity, color = membership_correlation )) + 
  scale_color_viridis_d() +
  geom_line(size = 0.5, show.legend = FALSE) +
 alexis_theme() +
  facet_wrap(facets = vars(cluster))

cluster_membership_plot1_reps

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/zscore_scaled_cluster_membership_plot_1_reps.png", plot = cluster_membership_plot1_reps, width = 16, height = 10, scale = 0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/zscore_scaled_cluster_membership_plot_1_reps.pdf", plot = cluster_membership_plot1_reps, width = 16, height = 10, scale = 0.4)

#-------------------------------------------------------------------------------------------
#make just a legend to save for image.
legend <- ggplot(data = data.frame(cluster_membership = c(0.3,0.27, 0.98,1), mean_condition_intensity = c(0, 0.25, 0.5,1), minute = c("a", "b", "c", "d") )) +
  geom_line(mapping = aes(x = minute, y = mean_condition_intensity , color = cluster_membership)) +
  scale_color_viridis_c() + 
  theme(legend.ticks = element_line(linewidth = 0.25, color = "white"),
        legend.ticks.length = unit(1.5, "mm"))
        # legend.frame = element_rect(color = "black"))

legend
ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/zscore_scaled_cluster_membership_plot_1_reps_legend.png", plot = legend, width = 3, height = 6, scale = 0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/zscore_scaled_cluster_membership_plot_1_reps_legend.pdf", plot = legend, width = 3, height = 6, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------

#mean intensity of clusters -------------------------------------
cluster_membership_plot2 <- ggplot((data_clustered_reps_df %>%
                                     group_by(cluster, minute) %>%
                                     summarize(mean_intensity = mean(mean_condition_intensity) + 1, count_obs = n()) %>%
                                      ungroup()), aes(x = minute, y = mean_intensity, color = count_obs )) +
  scale_color_viridis_c() +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  geom_point(size = 2, show.legend = TRUE) +
  geom_line(size = 2, show.legend = TRUE, alpha = 0.4) +
  geom_point(size = 6, alpha = 0.4, show.legend = TRUE) +
  alexis_theme() +
  scale_y_continuous(limits = c(-1, 2.5), expand = c(0,0), breaks = c(seq(0, 2, 1))) +
  facet_wrap(facets = vars(cluster))

cluster_membership_plot2

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/zscore_scaled_cluster_membership_plot_2b.png", plot = cluster_membership_plot2, width = 16, height = 10, scale = 0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/zscore_scaled_cluster_membership_plot_2b.pdf", plot = cluster_membership_plot2, width = 16, height = 10, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------

#mean intensity of clusters -------------------------------------
cluster_membership_plot2c <- ggplot((data_clustered_reps_df %>%
                                     group_by(cluster, minute) %>%
                                     summarize(mean_intensity = mean(mean_condition_intensity) + 1, count_obs = n()) %>%
                                      ungroup()), aes(x = minute, y = mean_intensity, color = as.factor(cluster))) +
  scale_color_manual(values = c("#081D58","#225EA8","#41B6C4", "#C7E9B4")) +
  # scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  # geom_point(size = 2, show.legend = FALSE) +
  geom_line(size = 2, show.legend = FALSE, alpha = 0.4) +
  geom_point(size = 6, alpha = 0.4, show.legend = FALSE) +
  alexis_theme() +
  scale_y_continuous(limits = c(-1, 2.5), expand = c(0,0), breaks = c(seq(0, 2, 1))) +
  facet_wrap(facets = vars(cluster)) +
  theme(
    # strip.background = element_rect(fill = cluster, alpha = 0.7),
        strip.text = element_text(face = "bold", family = "sans", size = 14))

cluster_membership_plot2c

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/zscore_scaled_cluster_membership_plot_2c.png", plot = cluster_membership_plot2c, width = 16, height = 10, scale = 0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/zscore_scaled_cluster_membership_plot_2c.pdf", plot = cluster_membership_plot2c, width = 16, height = 10, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
#custom labeller function
# my_labeller <-  function (variable, value) {
#   if (variable == "cluster") {
#     paste("Cluster", value, sep = " ") } else {value} }

my_labeller <- function(labels) {
  # labels is a data frame with columns = facet variables
  if ("cluster" %in% names(labels)) {
    labels$cluster <- paste("Cluster", labels$cluster)
  }
  labels
}


## --------------------------------------------------------------------------------------------------------------------------------------------------

#mean intensity of clusters -------------------------------------
cluster_membership_plot2c_cluster1 <- ggplot((data_clustered_reps_df %>%
                                                group_by(cluster, minute) %>%
                                                summarize(mean_intensity = mean(mean_condition_intensity) + 1, count_obs = n()) %>%
                                                ungroup() %>%
                                                mutate(text_position_x = 12, text_position_y = 2) %>%
                                                filter(cluster == 1)),
                                             aes(x = minute, y = mean_intensity)) +
  # scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  # geom_point(size = 2, show.legend = FALSE) +
  geom_line(size = 2, show.legend = FALSE, alpha = 0.8, color = "#081D58") +
  geom_point(size = 4, alpha = 1, show.legend = FALSE, color = "#081D58") +
  alexis_theme() +
  scale_y_continuous(limits = c(-1, 2.5), expand = c(0,0), breaks = c(seq(0, 2, 1))) +
  geom_text(mapping = aes(x = text_position_x, y = text_position_y, label = paste("n = ", count_obs, "\nobservations" ,sep = "")),
            size = 3, color= "black", lineheight = 0.75) +
  # geom_text(tibble(x = 12, y = -0.15, label = "L.O.D", cluster = 1), mapping = aes(x = x, y = y, label = label), 
  #           color = "gray30", inherit.aes = FALSE, size= 3) +
  facet_wrap(facets = vars(cluster), labeller = my_labeller) +
  theme(
    strip.background = element_rect(fill =adjustcolor( "#081D58", alpha.f = 0.3)),
        strip.text = element_text(face = "bold", family = "sans", size = 14)) +
  ylab("avg. scaled intensity") +
  xlab("EGF (minute)")

cluster_membership_plot2c_cluster1

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/zscore_scaled_cluster_membership_plot_2c_cluster1.png", plot = cluster_membership_plot2c_cluster1, width = 6, height = 6, scale = 0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/zscore_scaled_cluster_membership_plot_2c_cluster1.pdf", plot = cluster_membership_plot2c_cluster1, width = 6, height = 6, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------


#mean intensity of clusters -------------------------------------
cluster_membership_plot2c_cluster2 <- ggplot((data_clustered_reps_df %>%
                                                group_by(cluster, minute) %>%
                                                summarize(mean_intensity = mean(mean_condition_intensity) + 1, count_obs = n()) %>%
                                                ungroup() %>%
                                                mutate(text_position_x = 12, text_position_y = 2) %>%
                                                filter(cluster == 2)),
                                             aes(x = minute, y = mean_intensity)) +
  # scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  # geom_point(size = 2, show.legend = FALSE) +
  geom_line(size = 2, show.legend = FALSE, alpha = 0.8, color = "#225EA8") +
  geom_point(size = 4, alpha = 1, show.legend = FALSE, color = "#225EA8") +
  alexis_theme() +
  scale_y_continuous(limits = c(-1, 2.5), expand = c(0,0), breaks = c(seq(0, 2, 1))) +
  geom_text(mapping = aes(x = text_position_x, y = text_position_y, label = paste("n = ", count_obs, "\nobservations" ,sep = "")),
            size = 3, color= "black", lineheight = 0.75) +
  # geom_text(tibble(x = 12, y = -0.15, label = "L.O.D", cluster = 2), mapping = aes(x = x, y = y, label = label), 
  #           color = "gray30", inherit.aes = FALSE, size= 3) +
  facet_wrap(facets = vars(cluster), labeller = my_labeller) +
  theme(
    strip.background = element_rect(fill =adjustcolor( "#225EA8", alpha.f = 0.3)),
        strip.text = element_text(face = "bold", family = "sans", size = 14)) +
  ylab("avg. scaled intensity") +
  xlab("EGF (minute)")

cluster_membership_plot2c_cluster2

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/zscore_scaled_cluster_membership_plot_2c_cluster2.png", plot = cluster_membership_plot2c_cluster2, width = 6, height = 6, scale = 0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/zscore_scaled_cluster_membership_plot_2c_cluster2.pdf", plot = cluster_membership_plot2c_cluster2, width = 6, height = 6, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------

#mean intensity of clusters -------------------------------------
cluster_membership_plot2c_cluster3 <- ggplot((data_clustered_reps_df %>%
                                                group_by(cluster, minute) %>%
                                                summarize(mean_intensity = mean(mean_condition_intensity) + 1, count_obs = n()) %>%
                                                ungroup() %>%
                                                mutate(text_position_x = 12, text_position_y = 2) %>%
                                                filter(cluster == 3)),
                                             aes(x = minute, y = mean_intensity)) +
  # scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  # geom_point(size = 2, show.legend = FALSE) +
  geom_line(size = 2, show.legend = FALSE, alpha = 0.8, color = "#41B6C4") +
  geom_point(size = 4, alpha = 1, show.legend = FALSE, color = "#41B6C4") +
  alexis_theme() +
  scale_y_continuous(limits = c(-1, 2.5), expand = c(0,0), breaks = c(seq(0, 2, 1))) +
  geom_text(mapping = aes(x = text_position_x, y = text_position_y - 1, label = paste("n = ", count_obs, "\nobservations" ,sep = "")),
            size = 3, color= "black", lineheight = 0.75) +
  # geom_text(tibble(x = 12, y = -0.15, label = "L.O.D", cluster = 3), mapping = aes(x = x, y = y, label = label), 
  #           color = "gray30", inherit.aes = FALSE, size= 3) +
  facet_wrap(facets = vars(cluster), labeller = my_labeller) +
  theme(
    strip.background = element_rect(fill =adjustcolor( "#41B6C4", alpha.f = 0.3)),
        strip.text = element_text(face = "bold", family = "sans", size = 14)) +
  ylab("avg. scaled intensity") +
  xlab("EGF (minute)")

cluster_membership_plot2c_cluster3

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/zscore_scaled_cluster_membership_plot_2c_cluster3.png", plot = cluster_membership_plot2c_cluster3, width = 6, height = 6, scale = 0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/zscore_scaled_cluster_membership_plot_2c_cluster3.pdf", plot = cluster_membership_plot2c_cluster3, width = 6, height = 6, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------

#mean intensity of clusters -------------------------------------
cluster_membership_plot2c_cluster4 <- ggplot((data_clustered_reps_df %>%
                                                group_by(cluster, minute) %>%
                                                summarize(mean_intensity = mean(mean_condition_intensity) + 1, count_obs = n()) %>%
                                                ungroup() %>%
                                                mutate(text_position_x = 12, text_position_y = 2) %>%
                                                filter(cluster == 4)),
                                             aes(x = minute, y = mean_intensity)) +
  # scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  # geom_point(size = 2, show.legend = FALSE) +
  geom_line(size = 2, show.legend = FALSE, alpha = 0.8, color = "#C7E9B4") +
  geom_point(size = 4, alpha = 1, show.legend = FALSE, color = "#C7E9B4") +
  alexis_theme() +
  scale_y_continuous(limits = c(-1, 2.5), expand = c(0,0), breaks = c(seq(0, 2, 1))) +
  geom_text(mapping = aes(x = text_position_x, y = text_position_y, label = paste("n = ", count_obs, "\nobservations" ,sep = "")),
            size = 3, color= "black", lineheight = 0.75) +
  # geom_text(tibble(x = 12, y = -0.15, label = "L.O.D", cluster = 4), mapping = aes(x = x, y = y, label = label), 
  #           color = "gray30", inherit.aes = FALSE, size= 3) +
  facet_wrap(facets = vars(cluster), labeller = my_labeller) +
  theme(
    strip.background = element_rect(fill =adjustcolor( "#C7E9B4", alpha.f = 0.5)),
        strip.text = element_text(face = "bold", family = "sans", size = 14)) +
  ylab("avg. scaled intensity") +
  xlab("EGF (minute)")

cluster_membership_plot2c_cluster4

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/zscore_scaled_cluster_membership_plot_2c_cluster4.png", plot = cluster_membership_plot2c_cluster4, width = 6, height = 6, scale = 0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/zscore_scaled_cluster_membership_plot_2c_cluster4.pdf", plot = cluster_membership_plot2c_cluster4, width = 6, height = 6, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
ascore_gene_clusters_Imputation_noSD <- ascore_gene_noImputation %>% 
  left_join(y = data_clustered_reps_df %>%
              mutate(condition = case_when(
                minute == 0 ~ "EGF0min",
                minute == 1 ~ "EGF1min",
                minute == 3 ~ "EGF3min",
                minute == 5 ~ "EGF5min",
                minute == 15 ~ "EGF15min")) %>%
              select(-c( minute,  gene_ref)),
            by = c("psite", "gene",  "condition", "replicate")) #%>% 
  # mutate(median_norm_intensity = case_when(
  #                 is.infinite(median_norm_intensity) ~ 8,
  #                 TRUE ~ median_norm_intensity))


## --------------------------------------------------------------------------------------------------------------------------------------------------
conversion <- ascore_gene_clusters_Imputation_noSD %>% 
  distinct(gene_ref, ref, cluster, mod_res) %>% 
  mutate(gene_ref = str_replace_all(gene_ref, pattern = " ", replacement = "_"))


## --------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(15)


#annotation matrix
zscore_cmeans_clustering_matrix_reps_ordered_df <- zscore_cmeans_clustering_matrix_reps_noSD %>% 
  as.data.frame() %>% 
  select(contains("0min"), contains("1min"), contains("3min"), contains("EGF5min"), contains("EGF15min")) %>% 
  rownames_to_column(var = "gene_ref") %>% 
  left_join(y = data_clustered_reps_df, by = "gene_ref") %>%
  
  #add in EGF node depths here via left join column = depth
  left_join(y = WP437_EGFR_node_depths %>%
              mutate(depth = as.character(depth)),
            by = c("gene")) %>% 
  
  #add in previous database annotations
  left_join(y = gene_diff_abundance_3xDB %>%
              filter(overlap_my_data %in% c("this\nstudy", "both")) %>%
              mutate(overlap_my_data = case_when(
                overlap_my_data == "this\nstudy" ~ "this study",
                overlap_my_data == "both" ~ "both",
                TRUE ~ overlap_my_data)) %>%
              distinct(gene, overlap_my_data),
            by = "gene") %>%
  
  # #add in previous database annotations
  # left_join(y = (psite_diff_abundance_gene_2xDB_imputed_forvolcano %>%
  #             filter(overlap_my_data %in% c("this\nstudy", "both")) %>%
  #             #formatting for legend
  #             mutate(overlap_my_data = case_when(
  #               overlap_my_data == "this\nstudy" ~ "this study",
  #               overlap_my_data == "both" ~ "both", 
  #               TRUE ~ overlap_my_data)) %>% 
  #             distinct(gene,  overlap_my_data)),
  #           by = c("gene")) %>% 
  # select(gene_ref, gene, contains("min"), everything()) %>% 
  arrange(-desc(cluster), overlap_my_data) %>%
  
  mutate(
    mod_res = str_sub(psite, end = 1L)) %>% 
  distinct(gene, gene_ref, cluster, overlap_my_data) #all pY and depth is hard to interpret in figure. Therefore, just overlap_my_data w/ 2 colors.
  # distinct(gene_ref, cluster, mod_res, depth, overlap_my_data)

#pull out annotations from columns I want
row_annotations <- zscore_cmeans_clustering_matrix_reps_ordered_df %>%
  # left_join(y = data_completeness_summary %>%
  #             pivot_wider(id_cols = gene_ref, names_from = condition, values_from = observed) %>%
  #             distinct(gene_ref, EGF0min, EGF1min, EGF3min, EGF5min, EGF15min),
  #           by = "gene_ref") %>% 
  distinct( gene_ref, cluster, overlap_my_data) %>%
  # distinct(gene_ref, cluster, mod_res, depth, overlap_my_data) %>%
  select( cluster, overlap_my_data) %>%    #get rid of depth and mod_res annotations
  arrange(-desc(cluster), overlap_my_data) %>%
  # select( cluster, mod_res, depth) %>%
  as.matrix()

#turn to matrix to embed row names that match to pheatmap matrix, but return to dataframe to please function.
rownames(row_annotations) <- zscore_cmeans_clustering_matrix_reps_ordered_df$gene_ref
row_annotations <- as.data.frame(row_annotations)

#row annotation colors


my_colour = list(
    mod_res = c(S = "#440154FF", T = "#22A884FF", Y = "#FDE725FF"),
    cluster = c('1' = "#081D58", '2' = "#225EA8", '3'= "#41B6C4", '4' = "#C7E9B4"),
    # EGF15min = c("6" = "#003C30", "5" = "#01665E", "4"= "#35978F", "3" = "#80CDC1", "2" = "#D9D9D9", "1" = "#F5F5F5", "0" = "#F7F7F7"),
    # EGF5min = c("6" = "#003C30", "5" = "#01665E", "4"= "#35978F", "3" = "#80CDC1", "2" = "#D9D9D9", "1" = "#F5F5F5", "0" = "#F7F7F7"),
    # EGF3min = c("6" = "#003C30", "5" = "#01665E", "4"= "#35978F", "3" = "#80CDC1", "2" = "#D9D9D9", "1" = "#F5F5F5", "0" = "#F7F7F7"),
    # EGF1min = c("6" = "#003C30", "5" = "#01665E", "4"= "#35978F", "3" = "#80CDC1", "2" = "#D9D9D9", "1" = "#F5F5F5", "0" = "#F7F7F7"),
    # EGF0min = c("6" = "#003C30", "5" = "#01665E", "4"= "#35978F", "3" = "#80CDC1", "2" = "#D9D9D9", "1" = "#F5F5F5", "0" = "#F7F7F7"),
    depth = c("0" = "green", "1" = "green4","2" =  "darkgreen","3" =  "yellow","4" =  "yellow4","5" =  "orange", "6" = "darkorange2", "7" = "darkgoldenrod" ,"8" = "darkorange4" ,"9" = "red4", "10" = "darkmagenta" ,"11" = "#891BF2","12" =  "#7918F3","13" =  "#6A15F5","14" =  "#5B12F6" ,"15" = "#4C0FF7", "16" = "#3C0CF9", "17" =  "#1E06FC"),
    overlap_my_data = c("both" = "#22A884FF", "this study" = "#D95F02"))

# depth_colors <- colorRampPalette(c("grey", "purple", "blue"))
# depth_colors(22)

# reorder columns
zscore_cmeans_clustering_matrix_reps_ordered_df1 <- zscore_cmeans_clustering_matrix_reps_noSD %>% 
  as.data.frame() %>% 
  select(contains("0min"), contains("1min"), contains("3min"), contains("EGF5min"), contains("EGF15min")) %>%
  rownames_to_column(var = "gene_ref") %>% 
  left_join((conversion %>% distinct(gene_ref, cluster)), by = "gene_ref") %>%
  
  #pull out gene from gene_ref for joining
  separate(gene_ref, into = c("gene", "psite"), sep = "_", remove = FALSE) %>%
  
  #add in previous database annotations joined by gene level
  # left_join(y = gene_diff_abundance_gene_2xDB_imputed %>%
  left_join(y = gene_diff_abundance_3xDB %>%
              filter(overlap_my_data %in% c("this\nstudy", "both")) %>%
              mutate(overlap_my_data = case_when(
                overlap_my_data == "this\nstudy" ~ "this study",
                overlap_my_data == "both" ~ "both",
                TRUE ~ overlap_my_data)) %>%
              distinct(gene, overlap_my_data),
            by = "gene") %>%

  
  # #add in previous database annotations
  # left_join(y = psite_diff_abundance_gene_2xDB_imputed_forvolcano %>%
  #             filter(overlap_my_data %in% c("this\nstudy", "both")) %>%
  #             mutate(overlap_my_data = case_when(
  #               overlap_my_data == "this\nstudy" ~ "this study",
  #               overlap_my_data == "both" ~ "both",
  #               TRUE ~ overlap_my_data)) %>%
  #             distinct(gene,  overlap_my_data),
  #           by = c("gene")) %>%
  arrange(-desc(cluster), overlap_my_data) %>% 
  select(gene_ref, contains("min"), everything())
  

#for input into heatmap with fuzzy cmeans clustering orders (soft cluster)
zscore_cmeans_clustering_matrix_reps_ordered_matrix <- as.matrix(zscore_cmeans_clustering_matrix_reps_ordered_df1[,2:31])

rownames(zscore_cmeans_clustering_matrix_reps_ordered_matrix) <- zscore_cmeans_clustering_matrix_reps_ordered_df1[,1]



zscore_pheatmap <- pheatmap(zscore_cmeans_clustering_matrix_reps_ordered_matrix , cluster_cols = F, cluster_rows = F,
         cutree_rows = 4, cutree_cols = 5, colorRampPalette(c("#4575B4","#91BFDB",  "ghostwhite",  "#FC8D59",  "#D73027" ))(80),
         # clustering_distance_rows = "euclidean",
         
         fontsize_row = 2, annotation_row = row_annotations, annotation_colors = my_colour,
         treeheight_row = 0) 

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pheatmap_imputed_zscore_geneDBcomparison.pdf", plot = zscore_pheatmap , width = 20, height = 60, scale =0.8)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pheatmap_imputed_zscore_geneDBcomparison.png", plot = zscore_pheatmap , width = 8, height = 20, scale = 0.8, limitsize = FALSE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(15)


#annotation matrix
zscore_cmeans_clustering_matrix_reps_ordered_df <- zscore_cmeans_clustering_matrix_reps_noSD %>% 
  as.data.frame() %>% 
  select(contains("0min"), contains("1min"), contains("3min"), contains("EGF5min"), contains("EGF15min")) %>% 
  rownames_to_column(var = "gene_ref") %>% 
  left_join(y = data_clustered_reps_df, by = "gene_ref") %>%
  
  #add in EGF node depths here via left join column = depth
  left_join(y = WP437_EGFR_node_depths %>%
              mutate(depth = as.character(depth)),
            by = c("gene")) %>% 
  
  # #add in previous database annotations
  # left_join(y = gene_diff_abundance_3xDB %>%
  #             filter(overlap_my_data %in% c("this\nstudy", "both")) %>%
  #             mutate(overlap_my_data = case_when(
  #               overlap_my_data == "this\nstudy" ~ "this study",
  #               overlap_my_data == "both" ~ "both", 
  #               TRUE ~ overlap_my_data)) %>% 
  #             distinct(gene, overlap_my_data),
  #           by = "gene") %>% 
  
  #add in previous database annotations
  left_join(y = (psite_diff_abundance_gene_2xDB_imputed_forvolcano %>%
              filter(overlap_my_data %in% c("this\nstudy", "both")) %>% #this will keep only psites that are significantly regulated according to fold change tests using standard deviation during imputation.
              #formatting for legend
              mutate(overlap_my_data = case_when(
                overlap_my_data == "this\nstudy" ~ "this study",
                overlap_my_data == "both" ~ "both", 
                TRUE ~ overlap_my_data)) %>% 
              distinct(gene, gene_ref,mod_residue, overlap_my_data)),
            by = c("gene", "gene_ref")) %>% 
  select(gene_ref, contains("min"), everything()) %>% 
  arrange(-desc(cluster), overlap_my_data) %>% 
  
  mutate(
    mod_res = str_sub(psite, end = 1L)) %>% 
  distinct(gene_ref, cluster, overlap_my_data) #all pY and depth is hard to interpret in figure. Therefore, just overlap_my_data w/ 2 colors.
  # distinct(gene_ref, cluster, mod_res, depth, overlap_my_data)

#pull out annotations from columns I want
row_annotations <- zscore_cmeans_clustering_matrix_reps_ordered_df %>%
  # left_join(y = data_completeness_summary %>%
  #             pivot_wider(id_cols = gene_ref, names_from = condition, values_from = observed) %>%
  #             distinct(gene_ref, EGF0min, EGF1min, EGF3min, EGF5min, EGF15min),
  #           by = "gene_ref") %>% 
  distinct(gene_ref, cluster, overlap_my_data) %>%
  # distinct(gene_ref, cluster, mod_res, depth, overlap_my_data) %>%
  select( cluster, overlap_my_data) %>%    #get rid of depth and mod_res annotations
  # select( cluster, mod_res, depth) %>%
  as.matrix()

#turn to matrix to embed row names that match to pheatmap matrix, but return to dataframe to please function.
rownames(row_annotations) <- zscore_cmeans_clustering_matrix_reps_ordered_df$gene_ref
row_annotations <- as.data.frame(row_annotations)

#row annotation colors


my_colour = list(
    mod_res = c(S = "#440154FF", T = "#22A884FF", Y = "#FDE725FF"),
    cluster = c('1' = "#081D58", '2' = "#225EA8", '3'= "#41B6C4", '4' = "#C7E9B4"),
    # EGF15min = c("6" = "#003C30", "5" = "#01665E", "4"= "#35978F", "3" = "#80CDC1", "2" = "#D9D9D9", "1" = "#F5F5F5", "0" = "#F7F7F7"),
    # EGF5min = c("6" = "#003C30", "5" = "#01665E", "4"= "#35978F", "3" = "#80CDC1", "2" = "#D9D9D9", "1" = "#F5F5F5", "0" = "#F7F7F7"),
    # EGF3min = c("6" = "#003C30", "5" = "#01665E", "4"= "#35978F", "3" = "#80CDC1", "2" = "#D9D9D9", "1" = "#F5F5F5", "0" = "#F7F7F7"),
    # EGF1min = c("6" = "#003C30", "5" = "#01665E", "4"= "#35978F", "3" = "#80CDC1", "2" = "#D9D9D9", "1" = "#F5F5F5", "0" = "#F7F7F7"),
    # EGF0min = c("6" = "#003C30", "5" = "#01665E", "4"= "#35978F", "3" = "#80CDC1", "2" = "#D9D9D9", "1" = "#F5F5F5", "0" = "#F7F7F7"),
    depth = c("0" = "green", "1" = "green4","2" =  "darkgreen","3" =  "yellow","4" =  "yellow4","5" =  "orange", "6" = "darkorange2", "7" = "darkgoldenrod" ,"8" = "darkorange4" ,"9" = "red4", "10" = "darkmagenta" ,"11" = "#891BF2","12" =  "#7918F3","13" =  "#6A15F5","14" =  "#5B12F6" ,"15" = "#4C0FF7", "16" = "#3C0CF9", "17" =  "#1E06FC"),
    overlap_my_data = c("both" = "#22A884FF", "this study" = "#D95F02"))

# depth_colors <- colorRampPalette(c("grey", "purple", "blue"))
# depth_colors(22)

# reorder columns
zscore_cmeans_clustering_matrix_reps_ordered_df1 <- zscore_cmeans_clustering_matrix_reps_noSD %>% 
  as.data.frame() %>% 
  select(contains("0min"), contains("1min"), contains("3min"), contains("EGF5min"), contains("EGF15min")) %>%
  rownames_to_column(var = "gene_ref") %>% 
  left_join((conversion %>% distinct(gene_ref, cluster)), by = "gene_ref") %>%
  
  #pull out gene from gene_ref for joining
  separate(gene_ref, into = c("gene", "psite"), sep = "_", remove = FALSE) %>%
  
  # #add in previous database annotations joined by gene level
  # left_join(y = gene_diff_abundance_gene_3xDB_imputed %>%
  # # left_join(y = gene_diff_abundance_3xDB %>%
  #             filter(overlap_my_data %in% c("this\nstudy", "both")) %>%
  #             mutate(overlap_my_data = case_when(
  #               overlap_my_data == "this\nstudy" ~ "this study",
  #               overlap_my_data == "both" ~ "both", 
  #               TRUE ~ overlap_my_data)) %>% 
  #             distinct(gene, overlap_my_data),
  #           by = "gene") %>% 

  
  #add in previous database annotations
  left_join(y = psite_diff_abundance_gene_2xDB_imputed_forvolcano %>%
              filter(overlap_my_data %in% c("this\nstudy", "both")) %>% #also keeps only sig chg sites w/ Imputation using SD
              mutate(overlap_my_data = case_when(
                overlap_my_data == "this\nstudy" ~ "this study",
                overlap_my_data == "both" ~ "both", 
                TRUE ~ overlap_my_data)) %>% 
              distinct(gene, gene_ref,mod_residue, overlap_my_data),
            by = c("gene", "gene_ref")) %>% 
  select(gene_ref, contains("min"), everything()) %>% 
  arrange(-desc(cluster), overlap_my_data) 

#for input into heatmap with fuzzy cmeans clustering orders (soft cluster)
zscore_cmeans_clustering_matrix_reps_ordered_matrix <- as.matrix(zscore_cmeans_clustering_matrix_reps_ordered_df1[,2:31])

rownames(zscore_cmeans_clustering_matrix_reps_ordered_matrix) <- zscore_cmeans_clustering_matrix_reps_ordered_df1[,1]



zscore_pheatmap <- pheatmap(zscore_cmeans_clustering_matrix_reps_ordered_matrix , cluster_cols = F, cluster_rows = F,
         cutree_rows = 4, cutree_cols = 5, colorRampPalette(c("#4575B4","#91BFDB",  "ghostwhite",  "#FC8D59",  "#D73027" ))(80),
         # clustering_distance_rows = "euclidean",
         
         fontsize_row = 2, annotation_row = row_annotations, annotation_colors = my_colour,
         treeheight_row = 0) 

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pheatmap_imputed_zscore_psite.pdf", plot = zscore_pheatmap , width = 20, height = 60, scale =0.8)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pheatmap_imputed_zscore_psite.png", plot = zscore_pheatmap , width = 8, height = 20, scale = 0.8, limitsize = FALSE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(15)
# 
# #vector of psites only observed 3 or more times across all 30 measurments
# well_observed_psites <- (psites_obs_gt_3 %>%
#   distinct(gene_ref, well_observed) %>%
#   filter(well_observed == TRUE))$gene_ref

#annotation matrix
zscore_cmeans_clustering_matrix_reps_ordered_df2 <- zscore_cmeans_clustering_matrix_reps_noSD %>% 
  as.data.frame() %>% 
  select(contains("0min"), contains("1min"), contains("3min"), contains("EGF5min"), contains("EGF15min")) %>% 
  rownames_to_column(var = "gene_ref") %>% 
  left_join(y = data_clustered_reps_df, by = "gene_ref") %>%
  mutate(membership_correlation = as.factor(str_sub(membership_correlation, end = 4L)))%>% #reduce decimals
  
  #add in EGF node depths here via left join column = depth
  left_join(y = WP437_EGFR_node_depths %>% mutate(depth = as.character(depth)), by = c("gene")) %>% 
  
  # #add in previous database annotations
  # left_join(y = gene_diff_abundance_4xDB %>%
  #             filter(overlap_my_data %in% c("this\nstudy", "both")) %>%
  #             mutate(overlap_my_data = case_when(
  #               overlap_my_data == "this\nstudy" ~ "this study",
  #               overlap_my_data == "both" ~ "both", 
  #               is.na(overlap_my_data) ~ "n.s.",
  #               TRUE ~ overlap_my_data)) %>% 
  #             distinct(gene, overlap_my_data),
  #           by = "gene") %>% 
#add in previous database annotations
  left_join(y = gene_diff_abundance_3xDB %>%
              filter(overlap_my_data %in% c("this\nstudy", "both")) %>%
              mutate(overlap_my_data = case_when(
                overlap_my_data == "this\nstudy" ~ "this study",
                overlap_my_data == "both" ~ "both",
                TRUE ~ overlap_my_data)) %>%
              distinct(gene, overlap_my_data),
            by = "gene") %>%
  
  arrange(-desc(cluster), overlap_my_data) %>%
  
  mutate(
    mod_res = str_sub(psite, end = 1L)) %>% 
  distinct(gene_ref, cluster, membership_correlation,
           # mod_res, depth, #remove these, all pY and depth is too difficult to parse visually.
           overlap_my_data)

#pull out annotations from columns I want
row_annotations2 <- zscore_cmeans_clustering_matrix_reps_ordered_df2 %>%
  left_join(y = data_completeness_summary_noSD %>%
              pivot_wider(id_cols = gene_ref, names_from = condition, values_from = observed) %>%
              distinct(gene_ref, EGF0min, EGF1min, EGF3min, EGF5min, EGF15min) %>% 
              ungroup() %>% 
              group_by(gene_ref) %>% 
              mutate(sum_obs = as.character(sum(EGF0min, EGF1min, EGF3min, EGF5min, EGF15min))) %>% 
              ungroup(),
            by = "gene_ref") %>% 
  distinct(gene_ref, cluster, overlap_my_data, membership_correlation,
           # mod_res, depth,
           EGF0min, EGF1min, EGF3min, EGF5min, EGF15min, sum_obs) %>% 
  # filter(gene_ref %in% well_observed_psites) %>% 
  select( cluster, overlap_my_data, membership_correlation,
          # mod_res, depth,
          EGF15min, EGF5min, EGF3min, EGF1min, EGF0min, sum_obs) %>%
  as.matrix()

#turn to matrix to embed row names that match to pheatmap matrix, but return to dataframe to please function.
rownames(row_annotations2) <- zscore_cmeans_clustering_matrix_reps_ordered_df2$gene_ref
row_annotations2 <- as.data.frame(row_annotations2) %>% 
  mutate(membership_correlation = as.numeric(membership_correlation))

#row annotation colors
my_colour2 = list(
    mod_res = c("S" = "#440154FF", "T" = "#22A884FF", "Y" = "#FDE725FF"),
    cluster = c('1' = "#081D58", '2' = "#225EA8", '3'= "#41B6C4", '4' = "#C7E9B4"),
    EGF15min = c("6" = "#003C30", "5" = "#01665E", "4"= "#35978F", "3" = "#80CDC1", "2" = "#D9D9D9", "1" = "#F5F5F5", "0" = "#F7F7F7"),
    EGF5min = c("6" = "#003C30", "5" = "#01665E", "4"= "#35978F", "3" = "#80CDC1", "2" = "#D9D9D9", "1" = "#F5F5F5", "0" = "#F7F7F7"),
    EGF3min = c("6" = "#003C30", "5" = "#01665E", "4"= "#35978F", "3" = "#80CDC1", "2" = "#D9D9D9", "1" = "#F5F5F5", "0" = "#F7F7F7"),
    EGF1min = c("6" = "#003C30", "5" = "#01665E", "4"= "#35978F", "3" = "#80CDC1", "2" = "#D9D9D9", "1" = "#F5F5F5", "0" = "#F7F7F7"),
    EGF0min = c("6" = "#003C30", "5" = "#01665E", "4"= "#35978F", "3" = "#80CDC1", "2" = "#D9D9D9", "1" = "#F5F5F5", "0" = "#F7F7F7"),
    depth = c("0" = "green", "1" = "green4","2" =  "darkgreen","3" =  "yellow","4" =  "yellow4","5" =  "orange", "6" = "darkorange2", "7" = "darkgoldenrod" ,"8" = "darkorange4" ,"9" = "red4", "10" = "darkmagenta" ,"11" = "#891BF2","12" =  "#7918F3","13" =  "#6A15F5","14" =  "#5B12F6" ,"15" = "#4C0FF7", "16" = "#3C0CF9", "17" =  "#1E06FC"),
    overlap_my_data = c("both" = "#22A884FF", "this study" = "#D95F02", "n.s." = "white"),  #based on imputed data if in study.
    sum_obs = c("0" = "grey100", "1" = "grey98", "2" = "grey95", "3" = "grey90", "4" = "grey85","5" = "grey80","6" = "grey75","7" = "grey70","8" = "grey65","9" = "grey64","10"= "grey63","11" = "grey62","12" = "grey61","13" = "grey60","14" = "grey59","15" = "grey58","16" = "grey57","17" = "grey56","18" = "grey55","19" = "grey50","20" = "grey45","21" = "grey40","22" = "grey35","23" = "grey30","24" = "grey25","25" = "grey20","26" = "grey15","27" = "grey10","28" = "grey9","29" = "grey8","30" = "grey7"),
    
    membership_correlation = c("0.35 "= "#FF0000","0.36 "= "#FF0707","0.37 "= "#FF0F0F","0.38 "= "#FF1717","0.39 "= "#FF1F1F","0.4 "= "#FF2727","0.41 "= "#FF2F2F","0.42 "= "#FF3737","0.43 "= "#FF3F3F","0.44 "= "#FF4747","0.45 "= "#FF4F4F","0.46 "= "#FF5757","0.47 "= "#FF5F5F","0.48 "= "#FF6767","0.49 "= "#FF6F6F","0.5 "= "#FF7777","0.51 "= "#FF7F7F","0.52 "= "#FF8787","0.53 "= "#FF8F8F","0.54 "= "#FF9797","0.55 "= "#FF9F9F","0.56 "= "#FFA7A7","0.57 "= "#FFAFAF","0.58 "= "#FFB7B7","0.59 "= "#FFBFBF","0.6 "= "#FFC7C7","0.61 "= "#FFCFCF","0.62 "= "#FFD7D7","0.63 "= "#FFDFDF","0.64 "= "#FFE7E7","0.65 "= "#FFEFEF","0.66 "= "#FFF7F7","0.67 "= "#FFFFFF","0.68 "= "#F7FAF7","0.69 "= "#EFF5EF","0.7 "= "#E7F0E7","0.71 "= "#DFEBDF","0.72 "= "#D7E6D7","0.73 "= "#CFE1CF","0.74 "= "#C7DDC7","0.75 "= "#BFD8BF","0.76 "= "#B7D3B7","0.77 "= "#AFCEAF","0.78 "= "#A7C9A7","0.79 "= "#9FC49F","0.8 "= "#97C097","0.81 "= "#8FBB8F","0.82 "= "#87B687","0.83 "= "#7FB17F","0.84 "= "#77AC77","0.85 "= "#6FA76F","0.86 "= "#67A267","0.87 "= "#5F9E5F","0.88 "= "#579957","0.89 "= "#4F944F","0.9 "= "#478F47","0.91 "= "#3F8A3F","0.92 "= "#378537","0.93 "= "#2F812F","0.94 "= "#277C27","0.95 "= "#1F771F","0.96 "= "#177217","0.97 "= "#0F6D0F","0.98 "= "#076807","0.99 "= "#006400"))
    
    
    # membership_correlation = c("0.35 "= "#FFFFFF","0.36 "= "#FBFCFB","0.37 "= "#F7FAF7","0.38 "= "#F3F7F3","0.39 "= "#EFF5EF","0.4 "= "#EBF2EB","0.41 "= "#E7F0E7","0.42 "= "#E3EEE3","0.43 "= "#DFEBDF","0.44 "= "#DBE9DB","0.45 "= "#D7E6D7","0.46 "= "#D3E4D3","0.47 "= "#CFE1CF","0.48 "= "#CBDFCB","0.49 "= "#C7DDC7","0.5 "= "#C3DAC3","0.51 "= "#BFD8BF","0.52 "= "#BBD5BB","0.53 "= "#B7D3B7","0.54 "= "#B3D0B3","0.55 "= "#AFCEAF","0.56 "= "#ABCCAB","0.57 "= "#A7C9A7","0.58 "= "#A3C7A3","0.59 "= "#9FC49F","0.6 "= "#9BC29B","0.61 "= "#97C097","0.62 "= "#93BD93","0.63 "= "#8FBB8F","0.64 "= "#8BB88B","0.65 "= "#87B687","0.66 "= "#83B383","0.67 "= "#7FB17F","0.68 "= "#7BAF7B","0.69 "= "#77AC77","0.7 "= "#73AA73","0.71 "= "#6FA76F","0.72 "= "#6BA56B","0.73 "= "#67A267","0.74 "= "#63A063","0.75 "= "#5F9E5F","0.76 "= "#5B9B5B","0.77 "= "#579957","0.78 "= "#539653","0.79 "= "#4F944F","0.8 "= "#4B924B","0.81 "= "#478F47","0.82 "= "#438D43","0.83 "= "#3F8A3F","0.84 "= "#3B883B","0.85 "= "#378537","0.86 "= "#338333","0.87 "= "#2F812F","0.88 "= "#2B7E2B","0.89 "= "#277C27","0.9 "= "#237923","0.91 "= "#1F771F","0.92 "= "#1B741B","0.93 "= "#177217","0.94 "= "#137013","0.95 "= "#0F6D0F","0.96 "= "#0B6B0B","0.97 "= "#076807","0.98 "= "#036603","0.99 "= "#006400"))
    
    # membership_correlation = c("0.35 "= "#BEBEBE","0.36 "= "#BDBFBD","0.37 "= "#BDC1BD","0.38 "= "#BDC2BD","0.39 "= "#BCC4BC","0.4 "= "#BCC5BC","0.41 "= "#BCC7BC","0.42 "= "#BBC8BB","0.43 "= "#BBCABB","0.44 "= "#BBCBBB","0.45 "= "#BACDBA","0.46 "= "#BACEBA","0.47 "= "#BAD0BA","0.48 "= "#B9D1B9","0.49 "= "#B9D3B9","0.5 "= "#B9D4B9","0.51 "= "#B9D6B9","0.52 "= "#B8D7B8","0.53 "= "#B8D9B8","0.54 "= "#B8DAB8","0.55 "= "#B7DCB7","0.56 "= "#B7DDB7","0.57 "= "#B7DFB7","0.58 "= "#B6E0B6","0.59 "= "#B6E2B6","0.6 "= "#B6E3B6","0.61 "= "#B5E5B5","0.62 "= "#B5E6B5","0.63 "= "#B5E8B5","0.64 "= "#B4E9B4","0.65 "= "#B4EBB4","0.66 "= "#B4ECB4","0.67 "= "#B4EEB4","0.68 "= "#B6EEB3","0.69 "= "#B8EFB2","0.7 "= "#BBEFB1","0.71 "= "#BDF0B0","0.72 "= "#BFF0AF","0.73 "= "#C2F1AE","0.74 "= "#C4F1AE","0.75 "= "#C6F2AD","0.76 "= "#C9F2AC","0.77 "= "#CBF3AB","0.78 "= "#CDF3AA","0.79 "= "#D0F4A9","0.8 "= "#D2F4A9","0.81 "= "#D4F5A8","0.82 "= "#D7F5A7","0.83 "= "#D9F6A6","0.84 "= "#DBF7A5","0.85 "= "#DEF7A4","0.86 "= "#E0F8A3","0.87 "= "#E2F8A3","0.88 "= "#E5F9A2","0.89 "= "#E7F9A1","0.9 "= "#E9FAA0","0.91 "= "#ECFA9F","0.92 "= "#EEFB9E","0.93 "= "#F0FB9E","0.94 "= "#F3FC9D","0.95 "= "#F5FC9C","0.96 "= "#F7FD9B","0.97 "= "#FAFD9A","0.98 "= "#FCFE99","0.99 "= "#FFFF99"))
    
    
#     membership_correlation =  c("0.35 "= "#BEBEBE",
# "0.36 "= "#BDB9BF","0.37 "= "#BCB4C1","0.38 "= "#BBAFC2","0.39 "= "#BAAAC4","0.40 "= "#B9A5C5","0.41 "= "#B8A0C7","0.42 "= "#B79BC8","0.43 "= "#B696CA","0.44 "= "#B591CC","0.45 "= "#B48CCD","0.46 "= "#B387CF","0.47 "= "#B282D0","0.48 "= "#B17DD2","0.49 "= "#B078D3","0.50 "="#AF73D5","0.51 "= "#AF6FD7","0.52 "= "#AE6AD8","0.53 "= "#AD65DA","0.54 "= "#AC60DB","0.55 "= "#AB5BDD","0.56 "= "#AA56DE","0.57 "="#A951E0","0.58 "= "#A84CE1","0.59 "= "#A747E3","0.60 "= "#A642E5","0.61 "= "#A53DE6","0.62 "= "#A438E8","0.63 "= "#A333E9","0.64 "="#A22EEB","0.65 "= "#A129EC","0.66 "= "#A024EE","0.67 "= "#A020F0","0.68 "= "#9B1FF0","0.69 "= "#961EF0","0.70 "= "#911DF1","0.71 "="#8C1CF1","0.72 "= "#871BF2","0.73 "= "#821AF2","0.74 "= "#7D19F3","0.75 "= "#7818F3","0.76 "= "#7317F4","0.77 "= "#6E16F4","0.78 "="#6915F5","0.79 "= "#6414F5","0.80 "= "#5F13F6","0.81 "= "#5912F6","0.82 "= "#5511F7","0.83 "= "#5010F7","0.84 "= "#4B0FF7","0.85 "="#460EF8","0.86 "= "#410DF8","0.87 "= "#3C0CF9","0.88 "= "#370BF9","0.89 "= "#310AFA","0.90 "= "#2C09FA","0.91 "= "#2808FB","0.92 "="#2307FB","0.93 "= "#1E06FC","0.94 "= "#1805FC","0.95 "= "#1304FD","0.96 "= "#0F03FD","0.97 "= "#0902FE","0.98 "= "#0501FE","0.99 "="#0000FF"))


# reorder columns
zscore_cmeans_clustering_matrix_reps_ordered_df12 <- zscore_cmeans_clustering_matrix_reps_noSD %>% 
  as.data.frame() %>% 
  select(contains("0min"), contains("1min"), contains("3min"), contains("EGF5min"), contains("EGF15min")) %>%
  rownames_to_column(var = "gene_ref") %>% 
  left_join((conversion %>% distinct(gene_ref, cluster, mod_res)), by = "gene_ref") %>%

  left_join(y = (data_completeness_summary_noSD %>%
              pivot_wider(id_cols = gene_ref, names_from = condition, values_from = observed) %>%
              distinct(gene_ref, EGF0min, EGF1min, EGF3min, EGF5min, EGF15min) %>% 
              ungroup() %>% 
              group_by(gene_ref) %>% 
              mutate(sum_obs = sum(EGF0min, EGF1min, EGF3min, EGF5min, EGF15min)) %>% 
              ungroup() %>% 
              distinct(gene_ref, sum_obs)),
            by = "gene_ref") %>%
  #also organize by overlap_my_data. To do so, need to add in the overlap annotations. must pull out gene from gene_ref first
  separate(gene_ref, into = c("gene", "psite"), sep = "_", remove = FALSE) %>% 
  
  # left_join(y = psite_diff_abundance_gene_2xDB_imputed_forvolcano %>%
  #             filter(overlap_my_data %in% c("this\nstudy", "both")) %>%
  #             #formatting how the legend looks
  #             mutate(overlap_my_data = case_when(
  #               overlap_my_data == "this\nstudy" ~ "this study",
  #               overlap_my_data == "both" ~ "both", 
  #               overlap_my_data == "n.s." ~ "n.s.", 
  #               TRUE ~ overlap_my_data)) %>% 
  #             distinct(gene, gene_ref,mod_residue, overlap_my_data),
  #           by = c("gene", "gene_ref")) %>% 
  
  
  #add in previous database annotations joined by gene level
  # left_join(y = gene_diff_abundance_gene_2xDB_imputed %>%
  left_join(y = gene_diff_abundance_3xDB %>%
              filter(overlap_my_data %in% c("this\nstudy", "both")) %>%
              mutate(overlap_my_data = case_when(
                overlap_my_data == "this\nstudy" ~ "this study",
                overlap_my_data == "both" ~ "both",
                TRUE ~ overlap_my_data)) %>%
              distinct(gene, overlap_my_data),
            by = "gene") %>%
  left_join(y = data_clustered_reps_df %>% select(gene_ref, membership_correlation), by = "gene_ref") %>%
  mutate(membership_correlation = as.character(str_sub(membership_correlation, end = 4L)))%>% #reduce decimals
  
  
  select(gene_ref, contains("min"), everything()) %>% 
  arrange(-desc(cluster), overlap_my_data, desc(sum_obs)) %>% 
  mutate(sum_obs = as.character(sum_obs))

#for input into heatmap with fuzzy cmeans clustering orders (soft cluster)
zscore_cmeans_clustering_matrix_reps_ordered_matrix2 <- as.matrix(zscore_cmeans_clustering_matrix_reps_ordered_df12[,2:31])

rownames(zscore_cmeans_clustering_matrix_reps_ordered_matrix2) <- zscore_cmeans_clustering_matrix_reps_ordered_df12[,1]



zscore_pheatmap2 <- pheatmap(zscore_cmeans_clustering_matrix_reps_ordered_matrix2 ,legend_labels = "z-score", cluster_cols = F, cluster_rows = F,
         cutree_rows = 4, cutree_cols = 5, colorRampPalette(c("#4575B4","#91BFDB",  "ghostwhite",  "#FC8D59",  "#D73027" ))(80),
         # clustering_distance_rows = "euclidean",
         
         fontsize_row = 2, annotation_row = row_annotations2, annotation_colors = my_colour2,
         treeheight_row = 0)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pheatmap_imputed_zscore_withcompleteness.pdf", plot = zscore_pheatmap2 , width = 20, height = 60, scale =0.8)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pheatmap_imputed_zscore_withcompletenesslegend.pdf", plot = zscore_pheatmap2 , width = 8, height = 20, scale =0.8)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pheatmap_imputed_zscore_withcompleteness.png", plot = zscore_pheatmap2 , width = 10, height = 20, scale = 0.8, limitsize = FALSE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
write_csv(x = zscore_cmeans_clustering_matrix_reps_ordered_df2, file = "modified_data/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/zscore_cmeans_clustering_matrix_reps_ordered_df2.csv")


## --------------------------------------------------------------------------------------------------------------------------------------------------
# brewer.pal(12, "YlGnBu")
# # "#FFFFCC" "#A1DAB4" "#41B6C4" "#225EA8"
# "#FFFFD9" "#EDF8B1" "#C7E9B4" "#7FCDBB" "#41B6C4" "#1D91C0" "#225EA8" "#253494" "#081D58"
# 
# "#C7E9B4""#41B6C4""#225EA8""#081D58"


## --------------------------------------------------------------------------------------------------------------------------------------------------

# # membership_colors <- colorRampPalette(c("white","lightyellow", "cornsilk2" ,"thistle2",  "plum", "orchid", "darkorchid", "purple", "slateblue", "royalblue", "mediumblue",  "darkblue"))
# 
# # membership_colors <- colorRampPalette(c("darkblue", "royalblue", "slateblue", "purple", "orchid", "plum", "thistle",
# #                                         "darkseagreen", "darkseagreen2",  "yellowgreen","greenyellow", "#CCFF66", "#FFFF99", "yellow"))
# 
# membership_colors <- colorRampPalette(c("red", "white", "darkgreen"))
# membership_color_vector<- membership_colors(65)
# membership_color_vector
# #___________________________
# 
# membership_correlation_colors <- row_annotations2 %>%
#   distinct(membership_correlation) %>% 
#   mutate(membership_correlation_numeric = as.numeric(membership_correlation)) %>% 
#   arrange(membership_correlation_numeric) %>% 
#   mutate(membership_correlation = paste0("\"",membership_correlation, " \"= " )) %>% 
#   select(-membership_correlation_numeric)
# membership_correlation_colors$color <- paste0("\"", membership_color_vector, "\"")
# 
# vector_membership_annotation_colors <- membership_correlation_colors %>% 
#   mutate(color_to_corr = str_c(membership_correlation, color, sep = "")) %>% 
#   select(color_to_corr)
# rownames(vector_membership_annotation_colors) <- NULL
# 
# vector_membership_annotation_colors <- vector_membership_annotation_colors %>%
#   mutate(color_vec =  str_flatten(color_to_corr, collapse = ","))




## --------------------------------------------------------------------------------------------------------------------------------------------------
##PIVOT WIDER TO COMPARE TRENDS ACROSS ALL TIME POINTS, just condition using mean intensity
##FEEDS INTO ZSCORE AND CLUSTERING
  #separate dataframe for pivoting. 
  matrix_ascore_gene_imputed_df3_meancondition_noSD <- matrix_ascore_gene_imputed_t_2_again_noSD %>%
    ungroup() %>% 
    distinct(condition, gene_ref, imputed_gt_3_reps, condition_gene_ref_mean_intensity) %>% 
    mutate(gene_ref = str_replace_all(gene_ref, pattern = " ", replacement = "_")) %>% 
    
  #now identify psites missing completely at one time point but measured in others.
  pivot_wider(id_cols = c("gene_ref"),
              names_from = condition,
              values_from = condition_gene_ref_mean_intensity,
              values_fill = NA) %>% 
    pivot_longer(cols = contains("EGF"), names_to = "condition") %>% 
    rename(condition_gene_ref_mean_intensity = value) %>% 
    group_by(gene_ref, condition) %>%  #to allow for different value sampling during imputation in next step
    mutate(condition_gene_ref_mean_intensity = case_when(
      !is.na(condition_gene_ref_mean_intensity) ~ condition_gene_ref_mean_intensity,
      is.na(condition_gene_ref_mean_intensity) ~ 10)) %>% #impute to minimum signal for conditions lacking observations at all
    ungroup() %>% 
    #chose n = 24 to account for some cases where all 4 conditions X 6 reps were not measured. sampling to size = 1 means to keep one value.
    pivot_wider(id_cols = c("gene_ref"),
              names_from = condition,
              values_from = condition_gene_ref_mean_intensity,
              values_fill = NA) 


    ### intensity of 10 was added if no p-sites were detected in all 6 replicates for a given time point for a specific psite.


#z-score (scale) ----------------------------------------------------------------------------------
colnames_for_zscore_matrix_condition_noSD <- colnames(matrix_ascore_gene_imputed_df3_meancondition_noSD %>% select(-gene_ref))


cmeans_clustering_input_matrix_condition_noSD <- as.matrix(matrix_ascore_gene_imputed_df3_meancondition_noSD[,2:6])

rownames(cmeans_clustering_input_matrix_condition_noSD) <- as.matrix(matrix_ascore_gene_imputed_df3_meancondition_noSD[,1])

zscore_cmeans_clustering_matrix_condition_noSD <- apply(cmeans_clustering_input_matrix_condition_noSD, 1, scale) %>% t()
colnames(zscore_cmeans_clustering_matrix_condition_noSD) <- colnames_for_zscore_matrix_condition_noSD

zscore_cmeans_clustering_matrix_condition_noSD[is.nan(zscore_cmeans_clustering_matrix_condition_noSD)] <- 0


## --------------------------------------------------------------------------------------------------------------------------------------------------
#keep only observations with 3 or more measurements per psite
zscore_cmeans_clustering_condition_df_noSD <- zscore_cmeans_clustering_matrix_condition_noSD %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene_ref") #%>% 
  # left_join((psites_obs_gt_3 %>% distinct(gene_ref, well_observed)), by = "gene_ref") %>% 
  # filter(well_observed == TRUE)

#return to matrix after keeping only >=3 observations per psite
zscore_cmeans_clustering_matrix_condition_noSD <- as.matrix(zscore_cmeans_clustering_condition_df_noSD[,2:6])
rownames(zscore_cmeans_clustering_matrix_condition_noSD) <- as.matrix(zscore_cmeans_clustering_condition_df_noSD[,1])



## --------------------------------------------------------------------------------------------------------------------------------------------------
PSP_EGFR_pY_ref_df <- PSP_EGFR %>%
  # rename(psite = mod_residue) %>%
  mutate(psite = paste0(mod_res, mod_residue)) %>% 
  mutate(gene_ref = paste(gene, psite, sep = "_")) 
  

# reorder columns
zscore_EGFR_df1_condition_noSD <- zscore_cmeans_clustering_matrix_condition_noSD %>% 
  as.data.frame() %>% 
  select(contains("0min"), contains("1min"), contains("3min"), contains("EGF5min"), contains("EGF15min")) %>%
  rownames_to_column(var = "gene_ref") %>% 
  separate(gene_ref, into = c("gene", "psite"), sep = "_", remove = FALSE) %>%
  # filter(gene_ref %in% well_observed_psites) %>%
  # filter(grepl("_Y", gene_ref) == TRUE) %>% 
  # filter(gene %in% ((PSP_EGFR )$gene))
  filter(gene %in% PSP_EGFR_pY_ref_df$gene) %>% 
  arrange(gene)
  # filter(gene %in% ((PSP_EGFR %>% filter(mod_res == "Y") %>% rename(psite = mod_residue))$gene))
  # left_join((conversion %>% distinct(gene_ref, cluster)), by = "gene_ref") %>% 
  # arrange(-desc(cluster)) 

#for input into heatmap with fuzzy cmeans clustering orders (soft cluster)
zscore_EGFR_matrix_condition_noSD <- as.matrix(zscore_EGFR_df1_condition_noSD[,4:8])

rownames(zscore_EGFR_matrix_condition_noSD) <- zscore_EGFR_df1_condition_noSD[,1]


#plot heatmap
zscore_EGFR_pheatmap_condition_noSD <- pheatmap(zscore_EGFR_matrix_condition_noSD , cluster_cols = F, cluster_rows = T,
         cutree_rows = 3, cutree_cols = 5, colorRampPalette(c("#4575B4","#91BFDB",  "ghostwhite",  "#FC8D59",  "#D73027" ))(80),
         clustering_distance_rows = "minkowski",
         
         fontsize_row = 7.5,
         # annotation_row = row_annotations, annotation_colors = my_colour,
         treeheight_row = 1)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pheatmap_imputed_zscore_EGFR_condition_PSP_gene_noSD.pdf", plot = zscore_EGFR_pheatmap_condition_noSD , width = 12, height = 60, scale =0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pheatmap_imputed_zscore_EGFR_condition_PSP_gene_noSD.png", plot = zscore_EGFR_pheatmap_condition_noSD , width = 5, height = 40, scale = 0.5, limitsize = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
PSP_EGFR_pY_ref_df <- PSP_EGFR %>%
  mutate(psite = paste0(mod_res, mod_residue)) %>%
  mutate(gene_ref = paste(gene, psite, sep = "_")) #%>% 
  # filter(mod_res == "Y")

#somehow this code relies on gene names from ascore_gene_fasta, and doesn't play well with human_fasta_2024 genenames. Can I make a gene conversion that satisfies all analyses?
  

# reorder columns
zscore_EGFR_df1_condition_noSD <- zscore_cmeans_clustering_matrix_condition_noSD %>% 
  as.data.frame() %>% 
  select(contains("0min"), contains("1min"), contains("3min"), contains("EGF5min"), contains("EGF15min")) %>%
  rownames_to_column(var = "gene_ref") %>% 
  separate(gene_ref, into = c("gene", "psite"), sep = "_", remove = FALSE) %>%
  # filter(gene_ref %in% well_observed_psites) %>%
  # filter(grepl("_Y", gene_ref) == TRUE) %>% 
  # filter(gene %in% ((PSP_EGFR )$gene))
  filter(gene_ref %in% PSP_EGFR_pY_ref_df$gene_ref) %>% #filtering for psite level matching here
  arrange(gene, gene_ref)
  # filter(gene %in% ((PSP_EGFR %>% filter(mod_res == "Y") %>% rename(psite = mod_residue))$gene))
  # left_join((conversion %>% distinct(gene_ref, cluster)), by = "gene_ref") %>% 
  # arrange(-desc(cluster)) 

#for input into heatmap with fuzzy cmeans clustering orders (soft cluster)
zscore_EGFR_matrix_condition_noSD <- as.matrix(zscore_EGFR_df1_condition_noSD[,4:8])

rownames(zscore_EGFR_matrix_condition_noSD) <- zscore_EGFR_df1_condition_noSD[,1]


#plot heatmap
zscore_EGFR_pheatmap_condition_noSD <- pheatmap(zscore_EGFR_matrix_condition_noSD , cluster_cols = F, cluster_rows = F,
         cutree_rows = 3, cutree_cols = 5, colorRampPalette(c("#4575B4","#91BFDB",  "ghostwhite",  "#FC8D59",  "#D73027" ))(80),
         clustering_distance_rows = "minkowski",
         fontsize_col = 18,
         fontsize_row = 10,
         # annotation_row = row_annotations, annotation_colors = my_colour,
         treeheight_row = 1)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pheatmap_imputed_zscore_EGFR_condition_PSP_generef_noSD.pdf", plot = zscore_EGFR_pheatmap_condition_noSD , width = 8, height = 40, scale =0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pheatmap_imputed_zscore_EGFR_condition_PSP_generef_noSD.png", plot = zscore_EGFR_pheatmap_condition_noSD , width = 8, height = 40, scale = 0.4, limitsize = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
PTMSigDB_EGFR_pY_ref_df <- EGFR1_pathway_PTMSigDB %>%
  mutate(psite = paste0(mod_res, mod_residue)) %>%
  mutate(gene_ref = paste(gene, psite, sep = "_")) 
  

# reorder columns
zscore_EGFR_df1_condition_PTMSigDB_noSD <- zscore_cmeans_clustering_matrix_condition_noSD %>% 
  as.data.frame() %>% 
  select(contains("0min"), contains("1min"), contains("3min"), contains("EGF5min"), contains("EGF15min")) %>%
  rownames_to_column(var = "gene_ref") %>% 
  separate(gene_ref, into = c("gene", "psite"), sep = "_", remove = FALSE) %>%
  # filter(grepl("_Y", gene_ref) == TRUE) %>% 
  filter(gene %in% PTMSigDB_EGFR_pY_ref_df$gene) %>% 
  arrange(gene)
  # filter(gene %in% ((PTMSigDB_EGFR %>% filter(mod_res == "Y") %>% rename(psite = mod_residue))$gene))
  # left_join((conversion %>% distinct(gene_ref, cluster)), by = "gene_ref") %>% 
  # arrange(-desc(cluster)) 

#for input into heatmap with fuzzy cmeans clustering orders (soft cluster)
zscore_EGFR_matrix_condition_PTMSigDB_noSD <- as.matrix(zscore_EGFR_df1_condition_PTMSigDB_noSD[,4:8])

rownames(zscore_EGFR_matrix_condition_PTMSigDB_noSD) <- zscore_EGFR_df1_condition_PTMSigDB_noSD[,1]


#plot heatmap
zscore_EGFR_pheatmap_condition_PTMSigDB_noSD <- pheatmap(zscore_EGFR_matrix_condition_PTMSigDB_noSD , cluster_cols = F, cluster_rows = T,
         cutree_rows = 3, cutree_cols = 5, colorRampPalette(c("#4575B4","#91BFDB",  "ghostwhite",  "#FC8D59",  "#D73027" ))(80),
         clustering_distance_rows = "minkowski",
         
         fontsize_row = 7.5,
         # annotation_row = row_annotations, annotation_colors = my_colour,
         treeheight_row = 1)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pheatmap_imputed_zscore_EGFR_condition_PTMSigDB_gene_noSD.pdf", plot = zscore_EGFR_pheatmap_condition_PTMSigDB_noSD , width = 12, height = 60, scale =0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pheatmap_imputed_zscore_EGFR_condition_PTMSigDB_gene_noSD.png", plot = zscore_EGFR_pheatmap_condition_PTMSigDB_noSD , width = 5, height = 40, scale = 0.5, limitsize = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
PTMSigDB_EGFR_pY_ref_df <- EGFR1_pathway_PTMSigDB %>%
  mutate(psite = paste0(mod_res, mod_residue)) %>%
  mutate(gene_ref = paste(gene, psite, sep = "_"))
  

# reorder columns
zscore_EGFR_df1_condition_PTMSigDB_gene_ref_noSD <- zscore_cmeans_clustering_matrix_condition_noSD %>% 
  as.data.frame() %>% 
  select(contains("0min"), contains("1min"), contains("3min"), contains("EGF5min"), contains("EGF15min")) %>%
  rownames_to_column(var = "gene_ref") %>% 
  separate(gene_ref, into = c("gene", "psite"), sep = "_", remove = FALSE) %>%
  # filter(gene_ref %in% well_observed_psites) %>%
  # filter(grepl("_Y", gene_ref) == TRUE) %>% 
  # filter(gene %in% ((PTMSigDB_EGFR )$gene))
  filter(gene_ref %in% PTMSigDB_EGFR_pY_ref_df$gene_ref) %>% 
  arrange(gene, gene_ref)
  # filter(gene %in% ((PTMSigDB_EGFR %>% filter(mod_res == "Y") %>% rename(psite = mod_residue))$gene))
  # left_join((conversion %>% distinct(gene_ref, cluster)), by = "gene_ref") %>% 
  # arrange(-desc(cluster)) 

#for input into heatmap with fuzzy cmeans clustering orders (soft cluster)
zscore_EGFR_matrix_condition_PTMSigDB_gene_ref_noSD <- as.matrix(zscore_EGFR_df1_condition_PTMSigDB_gene_ref_noSD[,4:8])

rownames(zscore_EGFR_matrix_condition_PTMSigDB_gene_ref_noSD) <- zscore_EGFR_df1_condition_PTMSigDB_gene_ref_noSD[,1]


#plot heatmap
zscore_EGFR_pheatmap_condition_PTMSigDB_gene_ref_noSD <- pheatmap(zscore_EGFR_matrix_condition_PTMSigDB_gene_ref_noSD , cluster_cols = F, cluster_rows = F,
         cutree_rows = 3, cutree_cols = 5, colorRampPalette(c("#4575B4","#91BFDB",  "ghostwhite",  "#FC8D59",  "#D73027" ))(80),
         clustering_distance_rows = "minkowski",
         fontsize_col = 18,
         fontsize_row = 10,
         # annotation_row = row_annotations, annotation_colors = my_colour,
         treeheight_row = 1)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pheatmap_imputed_zscore_EGFR_condition_PTMSigDB_generef_noSD.pdf", plot = zscore_EGFR_pheatmap_condition_PTMSigDB_gene_ref_noSD , width = 8, height = 40, scale =0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pheatmap_imputed_zscore_EGFR_condition_PTMSigDB_generef_noSD.png", plot = zscore_EGFR_pheatmap_condition_PTMSigDB_gene_ref_noSD , width = 8, height = 60, scale = 0.4, limitsize = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------


CST_EGFR_proteins <- c("P00533","P42336","Q9Y243","P42345","P23443","P15498","P63000","Q13153","Q13233","P45983","P05412","P01100","Q05397","P56945","P16333","P60953","P15941","P35222","P56539","P49023","P12830","P42226","P46108","P00519","Q9UPR0","P0DP24","Q16566","Q9UQM7","P22681","P62993","Q96J02","O14964","Q07889","O14807","P04049","Q02750","P01112","P28482","P05771","P12931","P01116","P52333","P07947","P43405","Q9UQC2","Q06124","Q99704","Q7L591","P42224","Q18PE1","P31751","P31749","P53779","P45984","Q6FG41","Q14289","P51636","Q03135","P52630","Q14765","P0DP25","Q13557","Q96NX5","P0DP23","Q13555","Q13554","P05129","P23458","Q6PKX4","Q9P104","Q07890","O60496","Q8TEW6","P36507","O60674","P17252","Q9NQ66","P19174","P51178"," Q4KWH8","P16885","O75038","Q00722","Q9UJM3","P52735","Q9UKW4","P15498","P16220","Q02930","O43889","P05771","P17252","P05129", "P19419", "P41970", "P28324", "P42336", "O00443", "O00750", "P23443", "P46734")

CST_gene_reference_conversion <- ascore_gene_clusters_Imputation_noSD %>% 
  filter(reference %in% CST_EGFR_proteins) %>% 
  distinct(reference, gene)
  

# reorder columns
zscore_EGFR_df1_condition_noSD_CST <- zscore_cmeans_clustering_matrix_condition_noSD %>% 
  as.data.frame() %>% 
  select(contains("0min"), contains("1min"), contains("3min"), contains("EGF5min"), contains("EGF15min")) %>%
  rownames_to_column(var = "gene_ref") %>% 
  separate(gene_ref, into = c("gene", "psite"), sep = "_", remove = FALSE) %>%
  # filter(gene_ref %in% well_observed_psites) %>%
  # filter(grepl("_Y", gene_ref) == TRUE) %>% 
  # filter(gene %in% ((PSP_EGFR )$gene))
  filter(gene %in% CST_gene_reference_conversion$gene) %>% 
  arrange(gene, gene_ref)
  # filter(gene %in% ((PSP_EGFR %>% filter(mod_res == "Y") %>% rename(psite = mod_residue))$gene))
  # left_join((conversion %>% distinct(gene_ref, cluster)), by = "gene_ref") %>% 
  # arrange(-desc(cluster)) 

#for input into heatmap with fuzzy cmeans clustering orders (soft cluster)
zscore_EGFR_matrix_condition_noSD_CST <- as.matrix(zscore_EGFR_df1_condition_noSD_CST[,4:8])

rownames(zscore_EGFR_matrix_condition_noSD_CST) <- zscore_EGFR_df1_condition_noSD_CST[,1]


#plot heatmap
zscore_EGFR_pheatmap_condition_noSD_CST <- pheatmap(zscore_EGFR_matrix_condition_noSD_CST , cluster_cols = F, cluster_rows = F,
         cutree_rows = 3, cutree_cols = 5, colorRampPalette(c("#4575B4","#91BFDB",  "ghostwhite",  "#FC8D59",  "#D73027" ))(80),
         clustering_distance_rows = "minkowski",
         
         fontsize_row = 7.5,
         # annotation_row = row_annotations, annotation_colors = my_colour,
         treeheight_row = 1)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/EGF_pathway/pheatmap_imputed_zscore_EGFR_condition_CST_gene_noSD.pdf", plot = zscore_EGFR_pheatmap_condition_noSD_CST , width = 12, height = 60, scale =0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/EGF_pathway/pheatmap_imputed_zscore_EGFR_condition_CST_gene_noSD.png", plot = zscore_EGFR_pheatmap_condition_noSD_CST , width = 5, height = 40, scale = 0.5, limitsize = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
WP437 <- c('NOS3','COX2','SHC1','STAT5A','RALA','NEDD4','ELK4','MAPK1','MAPK8','ARHGEF1','STMN1','STAT1','EPS15L1','SH3GL3','EGF','INPP5D','PTPN5','ATF1','GRB10','ERRFI1','GRIM19','MAP3K1','PTK6','RPS6KA2','CSK','PCNA','KRAS','LIMK2','RASA1','ARF6','IQGAP1','EGFR','SP1','BRAF','DOK2','PTPRR','CRK','SH2D2A','MAPK9','RPS6KA3','RIN1','CAV2','E2F1','VAV2','SOS2','PTPN12','MTOR','NCK2','MAP3K3','RALGDS','VAV3','CAMK2A','STAT3','MAPK3','SOS1','GAB2','RAP1A','JAK1','CDC42','NCOA3','ATXN2','RAB5A','RPS6KB1','MAPK14','ASAP1','SYNJ1','MAP2K1','RAF1','EPN1','NCK1','STAMBP','CAV1','STAM2','HGS','MAP2K2','CREB1','PLCG1','ERBB2','SPRY2','INPPL1','EIF4EBP1','MAP4K1','RICTOR','USP6NL','RPS6KA5','PTK2B','PIK3R2','PLD2','EPS8','RAC1','ABI1','PLSCR1','PIK3R1','MEF2D','GAB1','PLCE1','PLD1','DNM1','JUN','PIAS3','SH3GL2','RPS6KA1','ROCK1','MAP3K4','PRKCD','AP2B','AKT1','FOXO1','VAV1','PRKCZ','RALB','FOS','PRKCI','PRKCA','SRC','ITCH','ABL1','RALBP1','MAPK7','REPS2','JUND','GRB2','EPS15','USP8','AP2A1','HRAS','PDPK1','IQSEC1','PEBP1','SH3KBP1','TNK2','PTPN11','AP2S1','PRKCB','FOXO4','FOSB','PTEN','NEDD8','STAM','GJA1','CRKL','JAK2','MAP2K5','BCAR1','MAP3K2','ELK1','STAT5B','PIK3C2B','STXBP1','CBLB','CBL','CBLC','PAK1','CFL1','PTK2','AP2M1','MEF2A','PXN','MEF2C')
# reorder columns
zscore_EGFR_WP437_df1_condition_noSD <- zscore_cmeans_clustering_matrix_condition_noSD %>% 
  as.data.frame() %>% 
  select(contains("0min"), contains("1min"), contains("3min"), contains("EGF5min"), contains("EGF15min")) %>%
  rownames_to_column(var = "gene_ref") %>% 
  separate(gene_ref, into = c("gene", "psite"), sep = "_", remove = FALSE) %>%
  mutate(psite_num = as.numeric(str_sub(psite, start = 2L))) %>%
  # select(everything(), psite_num) %>% 
  # filter(gene_ref %in% well_observed_psites) %>% 
  # filter(grepl("_Y", gene_ref) == TRUE) %>%
  # filter(gene %in% ((PSP_EGFR_WP437 )$gene))
  filter(gene %in% WP437) %>% 
  left_join(data_clustered_reps_df %>%distinct(gene_ref, cluster), by = "gene_ref") %>%
  filter(!is.na(cluster)) %>% 
  left_join(y = WP437_EGFR_node_depths %>% distinct(), by = "gene") %>%
  distinct() %>%
  arrange(depth, gene, psite_num)
  # arrange(gene)
  # left_join((conversion %>% distinct(gene_ref, cluster)), by = "gene_ref") %>% 
  # arrange(-desc(cluster)) 

#for input into heatmap with fuzzy cmeans clustering orders (soft cluster)
zscore_EGFR_WP437_matrix_condition_noSD <- as.matrix(zscore_EGFR_WP437_df1_condition_noSD[,4:8])

rownames(zscore_EGFR_WP437_matrix_condition_noSD) <- zscore_EGFR_WP437_df1_condition_noSD[,1]


#annotations -----------------------------------------

#annotation matrix
annotation_WP437_depths_df_condition_noSD <- zscore_EGFR_WP437_df1_condition_noSD %>%
  # left_join(data_clustered_reps_df %>% select(gene_ref, cluster), by = "gene_ref") %>% 
  filter(!is.na(cluster)) %>% 
  mutate(depth = as.character(depth)) %>% 
  # mutate(
  #     mod_res = str_sub(psite, end = 1L)) %>% 
    distinct(gene_ref,
             cluster,
             # mod_res,
             depth)

#get the annotion contents into a matrix without the gene_ref ids, which will be added later as rownames
  row_annotations_WP437_matrix_df_condition_noSD <- annotation_WP437_depths_df_condition_noSD%>% 
    filter(!is.na(cluster)) %>% 
    # distinct(gene_ref, cluster, mod_res, depth) %>%
    select(cluster,
           # mod_res,
           depth) %>% 
    as.matrix()

  
      #assign row names to be gene_ref
rownames(row_annotations_WP437_matrix_df_condition_noSD) <- annotation_WP437_depths_df_condition_noSD$gene_ref
    
#return matrix with row names into a dataframe
  #INPUT THIS FOR ROW ANNOTATIONS IN PHEATMAP
row_annotations_WP437_matrix_df_condition_noSD <- as.data.frame(row_annotations_WP437_matrix_df_condition_noSD)
  
    #colors for heatmap
#row annotation colors
my_colour_WP437_condition = list(
    # mod_res = c(S = "#440154FF", T = "#22A884FF", Y = "#FDE725FF"),
    cluster = c('1' = "#081D58", '2' = "#225EA8", '3'= "#41B6C4", '4' = "#C7E9B4"),
    # EGF15min = c("6" = "#003C30", "5" = "#01665E", "4"= "#35978F", "3" = "#80CDC1", "2" = "#D9D9D9", "1" = "#F5F5F5", "0" = "#F7F7F7"),
    # EGF5min = c("6" = "#003C30", "5" = "#01665E", "4"= "#35978F", "3" = "#80CDC1", "2" = "#D9D9D9", "1" = "#F5F5F5", "0" = "#F7F7F7"),
    # EGF3min = c("6" = "#003C30", "5" = "#01665E", "4"= "#35978F", "3" = "#80CDC1", "2" = "#D9D9D9", "1" = "#F5F5F5", "0" = "#F7F7F7"),
    # EGF1min = c("6" = "#003C30", "5" = "#01665E", "4"= "#35978F", "3" = "#80CDC1", "2" = "#D9D9D9", "1" = "#F5F5F5", "0" = "#F7F7F7"),
    # EGF0min = c("6" = "#003C30", "5" = "#01665E", "4"= "#35978F", "3" = "#80CDC1", "2" = "#D9D9D9", "1" = "#F5F5F5", "0" = "#F7F7F7"),
    depth = c("0" = "chocolate4", "1" = "#9E0142","2" =  "#D53E4F","3" =  "#F46D43","4" =  "#FDAE61","5" =  "#FEE08B", "6" = "#FFFFBF", "7" = "#E6F598" ,"8" = "#ABDDA4" ,"9" = "#66C2A5", "10" = "#3288BD" ,"11" = "lightskyblue","12" =  "deepskyblue","13" =  "blue","14" =  "#5B12F6" ,"15" = "lightslateblue", "16" = "blueviolet", "17" =  "violet"))
# 

# "#C7E9B4""#41B6C4""#225EA8""#081D58"




#plot heatmap -------------------------------------------
zscore_EGFR_WP437_pheatmap_condition <- pheatmap(zscore_EGFR_WP437_matrix_condition_noSD ,
                   cluster_cols = F,
                   cluster_rows = F,
                 # cutree_rows = 6,
                 # cutree_cols = 5,
                 colorRampPalette(c("#4575B4","#91BFDB",  "ghostwhite",  "#FC8D59",  "#D73027" ))(80),
                 # clustering_distance_rows = "minkowski",
                 
                 fontsize_row = 10,
                 annotation_row = row_annotations_WP437_matrix_df_condition_noSD,
                 annotation_colors = my_colour_WP437_condition,
                 treeheight_row = 1)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pheatmap_imputed_zscore_EGFR_WP437_condition.pdf", plot = zscore_EGFR_WP437_pheatmap_condition , width = 10, height = 60, scale =0.4)

ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pheatmap_imputed_zscore_EGFR_WP437_condition.png", plot = zscore_EGFR_WP437_pheatmap_condition , width = 7, height = 40, scale = 0.5, limitsize = FALSE)

# #taller, less wide
# ggsave(filename = "output/MainFig3to6_SuppFig9to15/GlobalPhosphoSiteLevel/pheatmap_imputed_zscore_EGFR_WP437_condition_taller.pdf", plot = zscore_EGFR_WP437_pheatmap_condition , width = 10, height = 60, scale =0.4)

