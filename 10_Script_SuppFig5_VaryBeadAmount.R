## --------------------------------------------------------------------------------------------------------------------------------------------------
alexis_theme <- function() {
  theme(
    # panel.border = element_rect(colour = "blue", fill = NA, linetype = 2),
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.x  = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0,  hjust = 0.5, size = 16),
    # axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black", face = "plain", family = "sans"),
    axis.title = element_text(colour = "black", family = "sans", size = 16),
    axis.ticks = element_line(colour = "black"),
    axis.text.y = element_text(size = 16, hjust = 1),
    plot.title = element_text(size=10),
    # legend at the bottom 6)
    legend.position = "right")   
}


## --------------------------------------------------------------------------------------------------------------------------------------------------
my_qc_cvs <- function (data,
                       grouping,
                       condition,
                       intensity,
                       plot = TRUE,
                       plot_style = "density",
                       max_cv = 200,
                       xlab = "condition",
                       showlegend = FALSE,
                       yintercept = 25){
  
  
  protti_colours <- "placeholder"
  utils::data("protti_colours", envir = environment())
  
  
#-----------------------------------------------------------------------------------  
  if (plot == FALSE) {
    
    
    if (max(dplyr::pull(data, {{intensity}}), na.rm = TRUE) < 1000) {
      stop(strwrap("Please backtransform your data or use raw values.\nThe function does not handle log2 transformed data.", 
                   prefix = "\n", initial = ""))  }
    
    
    result <- data %>% dplyr::distinct({{grouping}}, {{condition}}, {{intensity}}) %>%
      tidyr::drop_na({{intensity}}) %>%
      dplyr::group_by({{grouping}}) %>%
      dplyr::mutate(cv_combined = (stats::sd({{intensity}})/mean({{intensity}})) * 100) %>%
      dplyr::group_by({{condition}}, {{grouping}}) %>%
      dplyr::mutate(cv = (stats::sd({{intensity}})/mean({{intensity}})) * 100) %>%
      dplyr::distinct({{condition}}, {{grouping}}, .data$cv_combined, .data$cv) %>%
      tidyr::drop_na() %>%
      dplyr::group_by({{condition}}) %>%
      dplyr::mutate(median_cv = stats::median(.data$cv)) %>% 
      dplyr::ungroup() %>%
      dplyr::mutate(median_cv_combined = stats::median(.data$cv_combined)) %>% 
      dplyr::select(-{{grouping}}, -c("cv_combined", "cv")) %>%
      dplyr::distinct() 
    return(result)  }  
  
  ##-----------------------------------------------------------------------------------   
  if (plot == TRUE) {
    if (max(dplyr::pull(data, {{intensity}}), na.rm = TRUE) < 1000) {
      stop(strwrap("Please backtransform your data or use raw values.\nThe function does not handle log2 transformed data.", 
                   prefix = "\n", initial = ""))}
    
    result <- data %>%
      dplyr::distinct({{grouping}}, {{condition}}, {{intensity}}) %>%
      tidyr::drop_na({{intensity}}) %>%
      
      dplyr::group_by({{grouping}}) %>%
      dplyr::mutate(
        cv_combined = (stats::sd({{intensity}})/mean({{intensity}})) * 100) %>%
      
      
      dplyr::group_by({{condition}}, {{grouping}}) %>%
      dplyr::mutate(cv = (stats::sd({{intensity}})/mean({{intensity}})) * 100) %>%
      dplyr::ungroup() %>%
      
      dplyr::distinct({{condition}}, {{grouping}}, .data$cv_combined, .data$cv) %>%
      tidyr::drop_na() %>%
      tidyr::pivot_longer(cols = starts_with("cv"), names_to = "type", values_to = "values") %>%
      dplyr::mutate(type = ifelse(.data$type =="cv", {{condition}}, "all")) %>%
      dplyr::mutate(type = forcats::fct_relevel(as.factor(.data$type),"all")) %>%
      dplyr::select(-{{condition}}) %>%
      dplyr::group_by(.data$type) %>%
      dplyr::mutate(median = stats::median(.data$values)) %>% 
      dplyr::distinct()
    
    if (max(result$values) > max_cv) {
      cv_too_high <- result %>%
        dplyr::filter(.data$values >max_cv) %>%
        nrow()
      
      warning(paste(cv_too_high), " values were exluded from the plot (CV > ",max_cv, " %)")  }
    
    
    
    ##-----------------------------------------------------------------------------------      
    if (plot_style == "boxplot") {
      plot <- ggplot2::ggplot(result) +
        ggplot2::geom_boxplot(aes(x = .data$type, y = .data$values, fill = .data$type), na.rm = TRUE, size = 1, show.legend = showlegend,
                              alpha   = 0.5, outlier.alpha = 0.5, outlier.shape = 21) + 
        ggplot2::labs(
          # title = "Coefficients of variation",
          y = "Coefficient of variation [%]", fill = "Condition") + 
        ggplot2::geom_hline(yintercept = {{yintercept}}, linetype = 1, size = 0.75, alpha = 0.3) +
        ggplot2::scale_y_continuous(limits = c(0, max_cv)) + 
        # scale_fill_brewer(palette = "Dark2") +
        scale_color_manual(values = c("black","black","black","black")) +
  scale_fill_manual(values = c(rep("grey100", 1), rep("grey90", 1), rep("grey70", 1),rep("grey50",1),rep("grey30", 1))) +
        # ggplot2::scale_fill_manual(values = c("grey",protti_colours)) +
        alexis_theme() +
        xlab({{xlab}})
        # ggplot2::theme(plot.title = ggplot2::element_text(size = 20),
        #                axis.title.x = ggplot2::element_text(size = 15), 
        #                axis.text.y = ggplot2::element_text(size = 15),
        #                axis.text.x = ggplot2::element_text(size = 12,angle = 75, hjust = 1),
        #                axis.title.y = ggplot2::element_text(size = 15), 
        #                legend.title = ggplot2::element_text(size = 15),
        #                legend.text = ggplot2::element_text(size = 15))
      return(plot)
    }
    
    ##-----------------------------------------------------------------------------------    
    if (plot_style == "density") {
      plot <- ggplot2::ggplot(result) +
        ggplot2::geom_density(ggplot2::aes(x = .data$values, col = .data$type), size = 1, na.rm = TRUE, show.legend = showlegend) + 
        ggplot2::labs(
          # title = "Coefficients of variation",
          x = "Coefficient of variation [%]", y = "Density", color = "Condition") +
        ggplot2::scale_x_continuous(limits = c(0,max_cv)) +
        geom_vline(data = dplyr::distinct(result,  .data$median, .data$type),
                   ggplot2::aes(xintercept = median, col = .data$type),
                   size = 1,
                   linetype = "dashed", 
                   show.legend = FALSE) +
        scale_fill_brewer(palette = "Dark2") +
        # ggplot2::scale_color_manual(values = c("grey",protti_colours)) +
        alexis_theme() +
        xlab(xlab)
        # ggplot2::theme(plot.title = ggplot2::element_text(size = 20), 
        #                axis.title.x = ggplot2::element_text(size = 15),
        #                axis.text.y = ggplot2::element_text(size = 15),
        #                axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust = 1),
        #                axis.title.y = ggplot2::element_text(size = 15),
        #                legend.title = ggplot2::element_text(size = 15),
        #                legend.text = ggplot2::element_text(size = 15))
      
      return(plot)
      
    }
    
    ##-----------------------------------------------------------------------------------
    if (plot_style == "violin") {
      
      plot <- ggplot2::ggplot(result, aes(x = .data$type, 
                                          y = .data$values, color = .data$type)) +
        ggplot2::geom_violin(na.rm = TRUE, size = 1) + 
        ggplot2::geom_boxplot(width = 0.15, fill = "white", na.rm = TRUE, alpha = 0.6, size = 0.75, outlier.color = NA) +
        ggplot2::labs(
          # title = "Coefficients of variation",
                      x = "", y = "Coefficient of variation [%]",
                      fill = "Condition") +
        ggplot2::geom_hline(yintercept = {{yintercept}}, linetype = 1, size = 0.75, alpha = 0.3) +
        ggplot2::scale_y_continuous(limits = c(0, max_cv)) +
        # ggplot2::scale_fill_manual(values = c("grey",  protti_colours)) +
        # scale_fill_brewer(palette = "Dark2") +
        scale_color_manual(values = c("black","black","black","black")) +
  scale_fill_manual(values = c(rep("grey90", 1), rep("grey70", 1),rep("grey50",1),rep("grey30", 1))) +
        alexis_theme() +
        xlab(xlab)
        # ggplot2::theme_bw() +
        # ggplot2::theme(plot.title = ggplot2::element_text(size = 20),
        #                axis.title.x = ggplot2::element_text(size = 15),
        #                axis.text.y = ggplot2::element_text(size = 15), 
        #                axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust = 1), axis.title.y = ggplot2::element_text(size = 15),
        #                legend.title = ggplot2::element_text(size = 15), 
        #                legend.text = ggplot2::element_text(size = 15))
      
      return(plot)
      
    }
  }
}



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


## --------------------------------------------------------------------------------------------------------------------------------------------------
human_fasta_2024 <- read_csv("raw_data/SuppFig5/sp_reviewed_UP000005640_2024_03_21_forR.csv") 



## --------------------------------------------------------------------------------------------------------------------------------------------------

comet <- read_csv("raw_data/SuppFig5/comet/result.csv") %>%
  clean_names() %>%
  separate(sample_name, into = c("exp", "lysate","frxn", "condition", "well_replicate", "inj_vol"), sep = "_", remove = FALSE) %>% 
  mutate(
    condition = str_sub(condition, end = -7L), #just keep bead slurry volume, not "Halo"
    well = str_sub(well_replicate, start = 1L, end = 3L),
    replicate = as.numeric(str_sub(well_replicate, start = -1L)), 
    replicate = case_when(
      is.na(replicate) ~ 0,
      TRUE ~ as.numeric(replicate))) %>% 
  left_join(y = human_fasta_2024, by = c("reference")) %>% 
  
  ##rename column to 'intensity'
  rename(intensity = max_intensity_light_c2837) %>% 
  
  ##keep only forward hits
  filter(reverse == FALSE) %>% 
  
  ##classify phospho vs. no none
  mutate(phospho = case_when(
    grepl("\\@", sequence) ~ "phospho",
    TRUE~"none"))
  
  
  


## --------------------------------------------------------------------------------------------------------------------------------------------------

  ascore <- read_csv("raw_data/SuppFig5/ascore/result.csv") %>%
  clean_names() %>%
  separate(sample_name, into = c("exp", "lysate","frxn", "condition", "well_replicate", "inj_vol"), sep = "_", remove = FALSE) %>% 
  mutate(
    condition = str_sub(condition, end = -7L), #just keep bead slurry volume, not "Halo"
    well = str_sub(well_replicate, start = 1L, end = 3L),
    replicate = as.numeric(str_sub(well_replicate, start = -1L)), 
    replicate = case_when(
      is.na(replicate) ~ 0,
      TRUE ~ as.numeric(replicate))) %>% 
  
  
  ##join to FASTA
  left_join(y = human_fasta_2024, by = c("reference")) %>% 
  
  ##rename column to 'intensity'
  rename(intensity = max_intensity_light_c2837) %>% 
  
  ##keep only forward hits
  filter(reverse == FALSE) %>% 
  
  ##classify phospho vs. no none
  mutate(phospho = case_when(
    grepl("\\@", sequence) ~ "phospho",
    TRUE~"none"))
  
sample_key <- ascore %>% 
  distinct(raw_file_name, condition, replicate)
sample_key


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

set.seed(11)
plot_ppept_enrich_efficiency <- ggplot() +
  geom_bar(data = (enrichment_efficiency_count_df %>%
                     distinct(condition, avg_phosphopeptide_enrichment_efficiency)),
           mapping = aes(x = condition, y = avg_phosphopeptide_enrichment_efficiency),
           stat = "identity", alpha = 0.7, show.legend = FALSE) + 
  geom_jitter(data = enrichment_efficiency_count_df,
              mapping = aes(x = condition, y = individual_phosphopept_enrich_efficiency, fill = replicate), shape = 1, show.legend = FALSE, size = 1, width = 0.2) +
  geom_hline(yintercept = 0.96, linetype = 5, linewidth = 0.75, alpha = 0.3) +
  # annotate(geom = "text", x = 2.5, y = 0.94, label = "Avg 96%", size = 4) +
  # theme_bw(18) +
  # theme(legend.position = "none") +
  ylab("phospho / all peptides (count)") +
  xlab(expression(paste("Bead slurry volume (", mu, "L", ")"))) +
  expand_limits(y = c(0, 1)) +
  alexis_theme()+
  scale_y_continuous(expand = c(0,0),
                     breaks = c(0.00, 0.25, 0.50, 0.75, 0.96, 1.00),
                     labels = c("0.00", "0.25", "0.50", "0.75", "0.96", "1.00")) #removes whitespace below bars! custom ticks
  ## ggtitle("Phosphopeptide enrichment efficiency")

plot_ppept_enrich_efficiency
ggsave("output/SuppFig5/phospho_pept_enrich_efficiency_counts.png", plot = plot_ppept_enrich_efficiency, width = 5, height = 8, scale = 0.4)
ggsave("output/SuppFig5/phospho_pept_enrich_efficiency_counts.pdf", plot = plot_ppept_enrich_efficiency, width = 5, height = 8, scale = 0.4)


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

plot_ratio_phos_intensity <- ggplot(data = (enrichment_efficiency_intensity_df)) + 
  geom_bar(aes(x = condition, y = mean_ratio),
           stat = "identity",
           position = "dodge",
           fill = "skyblue2", show.legend = FALSE) +
  geom_point(aes(x = condition,
                 y = rep_ratio,
                 fill = replicate),
             position = position_jitterdodge(dodge.width = 0.2), alpha = 1, size =1, shape = 1, show.legend = FALSE)  +
geom_hline(yintercept = 0.995, linetype = 5, linewidth = 0.75, alpha = 0.3) +
  # annotate(geom = "text", x = 2.5, y = 0.94, label = "Avg 96%", size = 4) +
  # theme_bw(18) +
  # theme(legend.position = "none") +
  ylab("phospho / all peptides (count)") +
  expand_limits(y = c(0, 1)) +
  alexis_theme()+
  scale_y_continuous(expand = c(0,0),
                     breaks = c(0.00, 0.25, 0.50, 0.75, 0.995, 1.00),
                     labels = c("0.00", "0.25", "0.50", "0.75", "0.995", "1.00")) +#removes whitespace below bars! custom ticks
  ylab("phospho / all peptides (intensity)") + xlab(expression(paste("Bead slurry volume (", mu, "L", ")")))

plot_ratio_phos_intensity
ggsave("output/SuppFig5/phospho_pept_intensity_ratio.png", plot = plot_ratio_phos_intensity, width = 5, height = 8, scale = 0.4)
ggsave("output/SuppFig5/phospho_pept_intensity_ratio.pdf", plot = plot_ratio_phos_intensity, width = 5, height = 8, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_psm_intensity_distribution_reps <- ggplot(data = comet) +
  geom_hline(yintercept = 19.25, linewidth = 0.5, alpha = 0.3, linetype = 1) +
  geom_boxplot(mapping = aes(x = condition, y = log2(intensity), fill = as.factor(replicate)),
               outlier.shape = 21,
               outlier.alpha = 0.3,
               show.legend = FALSE,
               alpha = 0.5) +
  scale_fill_brewer(palette= "Dark2") +
  alexis_theme() +
  # theme_bw(18) +
  theme(axis.text.x = element_text(angle=0, hjust = 0.5), legend.position = "none") +
  scale_y_continuous(expand = c(0,0),
                     breaks = c(8, 16, 19.25, 24, 32),
                     labels = c("8", "16", "19.25", "24", "32")) +#removes whitespace below bars! custom ticks
  expand_limits(y = c(8, 34.5)) +
  ylab(expression(log[2]~(intensity))) + xlab(expression(paste("Bead slurry volume (", mu, "L", ")"))) +
  ggtitle("Intensity of all PSMs")

plot_psm_intensity_distribution_reps


ggsave("output/SuppFig5/all_psm_intensity_distribution_reps.png", plot = plot_psm_intensity_distribution_reps, width = 8, height = 8, scale = 0.4)

ggsave("output/SuppFig5/all_psm_intensity_distribution_reps.pdf", plot = plot_psm_intensity_distribution_reps, width = 8, height = 8, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
ascore_peptide_index <- ascore %>% 
  distinct() %>%
  
  mutate(seq_no_mods = str_remove_all(ascore_sequence, "[:punct:]"),
         seq_no_mods = str_remove_all(seq_no_mods, "[n]")) %>% ##remove n from n## mark of nAc modification. 
  mutate_at(c("a_score_1", "a_score_2", "a_score_3", "position_1", "position_2", "position_3"), funs(as.numeric)) %>% 
  mutate_at(c("a_score_1", "a_score_2", "a_score_3", "position_1", "position_2", "position_3"), ~replace_na(., 0)) %>%

  select(sample_name,  condition, replicate,  reference,
         ascore_sequence, num_sites, a_score_1, a_score_2, a_score_3, position_1, position_2, position_3,
         seq_no_mods, intensity, q_score_c2837, num_scans_light_c2837, charge, x_corr, 
         missed_cleavages, num_sites, redundancy) %>% 
 
##extract modified residue
  mutate(
    mod_res1 = str_sub(seq_no_mods, start = position_1+1L, end = position_1 + 1L),
    mod_res2 = str_sub(seq_no_mods, start = position_2+1L, end = position_2 + 1L),
    mod_res3 = str_sub(seq_no_mods, start = position_3+1L, end = position_3 + 1L)) %>% 


## pivot wider to match each p-site to an ascore and filter for STY

  unite(ascore_mod_res1, c(a_score_1, mod_res1, position_1), sep = "_") %>% 
  unite(ascore_mod_res2, c(a_score_2, mod_res2, position_2), sep = "_") %>%
  unite(ascore_mod_res3, c(a_score_3, mod_res3, position_3), sep = "_") %>%
  pivot_longer(cols = c("ascore_mod_res1", "ascore_mod_res2", "ascore_mod_res3")) %>% 
  separate(value, into = c("ascore", "mod_res", "mod_position"), sep = "_") %>% 
  mutate(ascore= as.numeric(ascore))




## --------------------------------------------------------------------------------------------------------------------------------------------------
ascore_peptide_index_stringent <- ascore_peptide_index %>% 
  filter(ascore >= 13) %>% 
  filter(grepl("[STY]", mod_res) == TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
  #join to FASTA
ascore_stringent_fasta <- ascore_peptide_index_stringent %>% 
  left_join(y = human_fasta_2024, by = c("reference")) %>% 
  
  #keep detected proteins only
  filter(!is.na(ascore_sequence)) %>% 
  
  mutate(
    mod_position = as.numeric(mod_position), ##mod_position parsed from ascore output, not Comet_Extended
    z_peptide = str_replace_all(seq_no_mods, "[IL]", "Z"),
    z_protein = str_replace_all(full_protein_sequence, "[IL]", "Z"),
    z_pept_position_in_z_protein = str_locate(z_protein, pattern = z_peptide), ##get peptide position from z-substituted peptide and protein
    mod_position_in_protein = z_pept_position_in_z_protein[,"start"] + mod_position, ##extract mod position from true protein sequence
    test_mod_position = str_sub(full_protein_sequence, start = mod_position_in_protein, end = mod_position_in_protein)) %>% 
  select(test_mod_position, mod_res, mod_position, ascore_sequence, everything()) %>% 
  unite(col = "mod_protein_location", c(mod_res, mod_position_in_protein), sep = "", remove = FALSE) %>% 
  unite(col = "ref", c(reference, mod_protein_location), sep = "_", remove = FALSE)
  





#index to protein (account for equal I and L during search)




## --------------------------------------------------------------------------------------------------------------------------------------------------
    ##tally up the number of rows that have NA vs. STY in test_mod_position
mod_pos_group <- ascore_stringent_fasta %>% 
  select(test_mod_position) %>% 
  group_by(test_mod_position) %>% 
   summarize(
    num_phospeptides = n())
mod_pos_group


## --------------------------------------------------------------------------------------------------------------------------------------------------
ascore_stringent_fasta_precursor <- ascore_stringent_fasta %>% 
  group_by(condition, replicate, ascore_sequence, charge, ascore) %>% 
  filter(intensity == max(intensity)) %>% 
  ungroup()


## --------------------------------------------------------------------------------------------------------------------------------------------------
#replicate
write_csv(x = ascore_stringent_fasta_precursor %>% distinct(condition, replicate, ascore_sequence, charge, ascore, intensity, reference, sample_name, gene, protein_names, mod_position_in_protein, ref), file = "modified_data/SuppFig5/distinct_pSTYpeptides_ascore13_replicate.csv", col_names = TRUE)


#condition
write_csv(x = ascore_stringent_fasta_precursor %>% distinct(condition,  ascore_sequence, charge, ascore, intensity, reference, gene, protein_names, mod_position_in_protein, ref), file = "modified_data/SuppFig5/distinct_pSTYpeptides_ascore13_condition.csv", col_names = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
##SUM PRECURSOR INTENSITIES TO P-SITES + median normalize

 ascore_stringent_fasta_psite <- ascore_stringent_fasta_precursor %>% 
  group_by(condition, replicate,  reference, mod_res,  ref) %>% #ref = unique p-site
  
  #sum PSM intensities of confident p-sites to individual p-sites
  mutate(
    sum_intensity_precursor_to_psite = sum(intensity),
    log2_psite_qty = log2(sum_intensity_precursor_to_psite)) %>% 
  ungroup() %>% 
  
  #keep single summed intensity per p-site prior to median normalization
  distinct(condition, replicate, reference, mod_res, ref, sum_intensity_precursor_to_psite, log2_psite_qty) 
  
  


## --------------------------------------------------------------------------------------------------------------------------------------------------
ascore_stringent_fasta_psite <- ascore_stringent_fasta_psite %>% 
  #median normalize intensities summed to p-site
  mutate(
    global_median_intensity = median(log2_psite_qty)) %>% 
  group_by(condition, replicate) %>% 
  mutate(
    sample_median_intensity = median(log2_psite_qty)) %>% 
  ungroup() %>% 
  mutate(
    median_norm_intensity = log2_psite_qty - sample_median_intensity + global_median_intensity,
    raw_median_norm_intensity = 2^median_norm_intensity) %>% 
  mutate(sample_id = paste(condition, replicate, sep = "_"))


## --------------------------------------------------------------------------------------------------------------------------------------------------
distinct_p_sites_condition <- ascore_stringent_fasta %>% 
  distinct(reference, mod_protein_location, mod_res, condition)
##this collapses all  replicates into a single set of unique p sites.



## --------------------------------------------------------------------------------------------------------------------------------------------------
write_csv(x = distinct_p_sites_condition, file = "modified_data/SuppFig5/distinct_pSTYsites_ascore13_condition.csv", col_names = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_distinct_p_sites_condition <- ggplot(data = distinct_p_sites_condition) +
  geom_bar(mapping = aes(x = condition, fill = mod_res)) +
  geom_text(data = tibble(x = c(1, 2, 3, 4),y = c(4094, 4323, 4505, 4477)), mapping = aes( x = x, y = y -150, label = paste0(y, " pY")),
            color = "black",
            size = 4,
            fontface = "bold",
            inherit.aes = FALSE) +
  alexis_theme()+
  scale_y_continuous(expand = c(0,0)) + #removes whitespace below bars!
  scale_fill_viridis_d() +
  ylab("unique phospho sites") +
  xlab(expression(paste("Bead slurry volume (", mu, "L", ")"))) 

plot_distinct_p_sites_condition

ggsave("output/SuppFig5/distinct_p_sites_condition.png", plot = plot_distinct_p_sites_condition, width = 12, height = 8, scale = 0.4)
ggsave("output/SuppFig5/distinct_p_sites_condition.pdf", plot = plot_distinct_p_sites_condition, width = 12, height = 8, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
summarize_distinct_psites_condition <- distinct_p_sites_condition %>%
  ungroup() %>% 
  group_by(condition, mod_res) %>% 
  summarize(
    n_psites = n())

summarize_distinct_psites_condition


## --------------------------------------------------------------------------------------------------------------------------------------------------
distinct_p_sites_reps <- ascore_stringent_fasta %>% 
  ungroup() %>% 
  distinct(reference, mod_protein_location, mod_res,  condition, replicate) %>% 
  mutate(
    sample_id = paste0(condition, "\n", replicate))
##this collapses all  replicates into a single set of unique p sites.


## --------------------------------------------------------------------------------------------------------------------------------------------------
write_csv(x = distinct_p_sites_reps, file = "modified_data/SuppFig5/pSTYsites_ascore13_reps.csv", col_names = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_distinct_p_sites_reps <- ggplot(data = distinct_p_sites_reps) +
  geom_bar(mapping = aes(x = sample_id, fill = mod_res)) +
  scale_fill_viridis_d() +
  ylab("unique phospho sites") +
  geom_text(data = tibble(x = c(1, 2, 3, 4, 5, 6 ,7, 8, 9, 10, 11, 12, 13, 14, 15, 16 ), y = c(2944, 3054, 3504, 3347, 3591, 3561, 3519, 3397, 3662, 3756, 3573, 3675, 3430, 3743, 3687, 3527)), mapping = aes( x = x, y = y -600, label = y),
            angle = -90,
            color = "black",
            size = 4.5,
            fontface = "bold",
            inherit.aes = FALSE) +
  alexis_theme()+
  scale_y_continuous(expand = c(0,0)) + #removes whitespace below bars!
  scale_fill_viridis_d() +
  ylab("unique phospho sites") +
  xlab(expression(paste("Bead slurry volume (", mu, "L", ")")))

plot_distinct_p_sites_reps

ggsave("output/SuppFig5/distinct_p_sites_reps.png", plot = plot_distinct_p_sites_reps, width = 16, height = 8, scale = 0.4)
ggsave("output/SuppFig5/distinct_p_sites_reps.pdf", plot = plot_distinct_p_sites_reps, width = 16, height = 8, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
summarize_distinct_psites_replicates <- distinct_p_sites_reps %>%
  ungroup() %>% 
  group_by(condition, replicate, sample_id, mod_res) %>% 
  summarize(
    n_psites = n())

summarize_distinct_psites_replicates


## --------------------------------------------------------------------------------------------------------------------------------------------------
# tibble(x = c(1, 2, 3, 4, 5, 6 ,7, 8, 9, 10, 11, 12 ), y = c(2944, 3054, 3504, 3347, 3591, 3561, 3519, 3397, 3662, 3756, 3573, 3675, 3430, 3743, 3687, 3527))


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
                     mutate(condition = fct_relevel(condition, "20", "40", "60", "80")) %>%
                     filter(mod_res == "Y")),
                   mapping = aes(x = condition, fill = mod_res),show.legend = FALSE, color = "black", size = 0.5) +
  
  
  
  
  geom_jitter(data = (distinct_p_sites_reps%>%
                        filter(grepl("pool", condition) == FALSE) %>% 
                        mutate(condition = fct_relevel(condition, "20", "40", "60", "80")) %>%
                        filter(mod_res == "Y")),
              mapping = aes(x = condition, color = sample_id, fill = mod_res),
              stat = "count", shape = 1, show.legend = FALSE, size = 1, stroke =0.5,
              width = 0.02, height = 0, alpha = 0.5) +
  
  geom_errorbar(data = (df_geom_errorbar %>%
                          filter(grepl("pool", condition) == FALSE) %>% 
                          mutate(condition = fct_relevel(condition, "20", "40", "60", "80")) %>%
                          filter(mod_res == "Y")),
                mapping = aes(x = condition, ymin = min_psites, ymax = max_psites),  show.legend = FALSE, linewidth = 0.25, width = 0.25, color = "grey30") +
   
  # theme_bw(10) +
  # facet_wrap(facets = vars(mod_res), ncol = 3) +
  scale_fill_viridis_d(direction = -1) +
  scale_color_manual(values = rep(c("grey30" ),60)) +
  ylab("unique phospho sites")   +
  alexis_theme()+
  scale_y_continuous(expand = c(0,0)) +
  # expand_limits(y = c(0, 1200)) +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0, size = 12))

plot_distinct_p_sites_condition_reps_wPTS

ggsave("output/SuppFig5/distinct_pY_sites_condition_wPTS.png", plot = plot_distinct_p_sites_condition_reps_wPTS, width = 6, height = 8, scale = 0.4)
ggsave("output/SuppFig5/distinct_pY_sites_condition_wPTS.pdf", plot = plot_distinct_p_sites_condition_reps_wPTS, width = 6, height = 8, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
# need each row to be a protein id and each column to be a cell type or sample.
df_psite_overlaps_reps <- distinct_p_sites_reps %>% 
  mutate(
    ref = paste(reference, mod_protein_location, sep = "_")) %>% 
  distinct(sample_id, ref) %>% 
  mutate(
    psite_present = 1) %>% 
  pivot_wider(
    names_from = sample_id, values_fill = 0, values_from = psite_present, id_cols = ref)

#unselect some columns
df_psite_overlaps_reps_less <- df_psite_overlaps_reps %>% 
  select(-ref) %>% 
  mutate_all(.funs = as.numeric) %>% 
  as.matrix() %>% #trick for formatting fxns within protti.
  as.data.frame()

#upset plot
upset_plot_psites_reps <- upset(data = df_psite_overlaps_reps_less, nsets = 16, order.by = "freq", text.scale = 1.5)
upset_plot_psites_reps

#save plot

png("output/SuppFig5/upset_plot_psites_reps.png", width = 800, height = 600)
print(upset_plot_psites_reps)
dev.off()


pdf("output/SuppFig5/upset_plot_psites_reps.pdf", width = 12, height = 8)
print(upset_plot_psites_reps)
dev.off()


## --------------------------------------------------------------------------------------------------------------------------------------------------
# need each row to be a protein id and each column to be a cell type or sample.
df_psite_overlaps_condition <- distinct_p_sites_condition %>% 
  mutate(
    ref = paste(reference, mod_protein_location, sep = "_")) %>% 
  distinct(condition, ref) %>% 
  mutate(
    psite_present = 1) %>% 
  pivot_wider(
    names_from = condition, values_fill = 0, values_from = psite_present, id_cols = ref)

#unselect some columns
df_psite_overlaps_condition_less <- df_psite_overlaps_condition %>% 
  select(-ref) %>% 
  mutate_all(.funs = as.numeric) %>% 
  as.matrix() %>% #trick for formatting fxns within protti.
  as.data.frame()

#upset plot
upset_plot_psites_condition <- upset(data = df_psite_overlaps_condition_less, nsets = 16, order.by = "freq", text.scale = 1.5)
upset_plot_psites_condition

#save plot

png("output/SuppFig5/upset_plot_psites_condition.png", width = 800, height = 600)
print(upset_plot_psites_condition)
dev.off()


pdf("output/SuppFig5/upset_plot_psites_condition.pdf", width = 12, height = 8)
print(upset_plot_psites_condition)
dev.off()


## --------------------------------------------------------------------------------------------------------------------------------------------------
# need each row to be a protein id and each column to be a cell type or sample.
df_psite_overlaps_reps_pY <- distinct_p_sites_reps %>%
  ungroup() %>% 
  filter(mod_res == "Y") %>% 
  mutate(
    ref = paste(reference, mod_protein_location, sep = "_")) %>% 
  distinct(sample_id, ref) %>% 
  mutate(
    psite_present = 1) %>% 
  pivot_wider(
    names_from = sample_id, values_fill = 0, values_from = psite_present, id_cols = ref)

#unselect some columns
df_psite_overlaps_reps_less_pY <- df_psite_overlaps_reps_pY %>% 
  select(-ref) %>% 
  mutate_all(.funs = as.numeric) %>% 
  as.matrix() %>% #trick for formatting fxns within protti.
  as.data.frame()

#upset plot
upset_plot_psites_reps_pY <- upset(data = df_psite_overlaps_reps_less_pY, nsets = 16, order.by = "freq", text.scale = 1.5)
upset_plot_psites_reps_pY

#save plot

png("output/SuppFig5/upset_plot_psites_reps_pY.png", width = 800, height = 600)
print(upset_plot_psites_reps_pY)
dev.off()


pdf("output/SuppFig5/upset_plot_psites_reps_pY.pdf", width = 12, height = 8)
print(upset_plot_psites_reps_pY)
dev.off()


## --------------------------------------------------------------------------------------------------------------------------------------------------
# need each row to be a protein id and each column to be a cell type or sample.
df_psite_overlaps_condition_pY <- distinct_p_sites_condition %>%
  ungroup() %>% 
  filter(mod_res == "Y") %>% 
  mutate(
    ref = paste(reference, mod_protein_location, sep = "_")) %>% 
  distinct(condition, ref) %>% 
  mutate(
    psite_present = 1) %>% 
  pivot_wider(
    names_from = condition, values_fill = 0, values_from = psite_present, id_cols = ref)

#unselect some columns
df_psite_overlaps_condition_less_pY <- df_psite_overlaps_condition_pY %>% 
  select(-ref) %>% 
  mutate_all(.funs = as.numeric) %>% 
  as.matrix() %>% #trick for formatting fxns within protti.
  as.data.frame()

#upset plot
upset_plot_psites_condition_pY <- upset(data = df_psite_overlaps_condition_less_pY, nsets = 16, order.by = "freq", text.scale = 1.5)
upset_plot_psites_condition_pY

#save plot

png("output/SuppFig5/upset_plot_psites_condition_pY.png", width = 800, height = 600)
print(upset_plot_psites_condition_pY)
dev.off()


pdf("output/SuppFig5/upset_plot_psites_condition_pY.pdf", width = 12, height = 8)
print(upset_plot_psites_condition_pY)
dev.off()


## --------------------------------------------------------------------------------------------------------------------------------------------------
psite_reps_CVs_NOTnorm <- my_qc_cvs(
  data = ascore_stringent_fasta_psite, 
  grouping = ref, 
  condition = condition,
  intensity = sum_intensity_precursor_to_psite,
  plot = TRUE, 
  plot_style = "boxplot",
  xlab = expression(paste("Bead slurry volume (", mu, "L", ")")),
  max_cv = 100)

psite_reps_CVs_NOTnorm

ggsave(filename = "output/SuppFig5/CVs_not_normalized.png", plot = psite_reps_CVs_NOTnorm, width = 8, height = 10, scale = 0.4)
ggsave(filename = "output/SuppFig5/CVs_not_normalized.pdf", plot = psite_reps_CVs_NOTnorm, width = 8, height = 10, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
psite_reps_CVs_YESnorm <- my_qc_cvs(
  data = ascore_stringent_fasta_psite, 
  grouping = ref, 
  condition = condition,
  intensity = raw_median_norm_intensity,
  plot = TRUE, 
  xlab = expression(paste("Bead slurry volume (", mu, "L", ")")),
  plot_style = "boxplot",
  max_cv = 100)

psite_reps_CVs_YESnorm

ggsave(filename = "output/SuppFig5/CVs_YES_normalized.png", plot = psite_reps_CVs_YESnorm, width = 8, height = 10, scale = 0.4)
ggsave(filename = "output/SuppFig5/CVs_YES_normalized.pdf", plot = psite_reps_CVs_YESnorm, width = 8, height = 10, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
psite_reps_CVs_NOTnorm_pY <- my_qc_cvs(
  data = ascore_stringent_fasta_psite %>% ungroup() %>% filter(mod_res =="Y"), 
  grouping = ref, 
  condition = condition,
  intensity = sum_intensity_precursor_to_psite,
  plot = TRUE, 
  plot_style = "boxplot",
  xlab = expression(paste("Bead slurry volume (", mu, "L", ")")),
  max_cv = 100)

psite_reps_CVs_NOTnorm_pY

ggsave(filename = "output/SuppFig5/CVs_not_normalized_pY.png", plot = psite_reps_CVs_NOTnorm_pY, width = 8, height = 10, scale = 0.4)
ggsave(filename = "output/SuppFig5/CVs_not_normalized_pY.pdf", plot = psite_reps_CVs_NOTnorm_pY, width = 8, height = 10, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_pY_intensity_distribution_reps <- ggplot(data = ascore_stringent_fasta_psite %>% ungroup() %>% filter(mod_res =="Y")) +
  # geom_hline(yintercept = 19.25, linewidth = 0.5, alpha = 0.3, linetype = 1) +
  geom_boxplot(mapping = aes(x = condition, y = log2_psite_qty, color = as.factor(replicate), fill = condition),
               outliers = FALSE,
               outlier.shape = 21,
               outlier.alpha = 0.3,
               show.legend = FALSE,
               alpha = 0.5) +
  # scale_fill_brewer(palette= "Dark2") +
  scale_color_manual(values = c("black","black","black","black")) +
  scale_fill_manual(values = c(rep("grey90", 1), rep("grey70", 1),rep("grey50",1),rep("grey30", 1))) +
  alexis_theme() +
  # theme_bw(18) +
  theme(axis.text.x = element_text(angle=0, hjust = 0.5), legend.position = "none") +
  scale_y_continuous(expand = c(0,0),
                     limits = c(12, 28),
                     breaks = c(seq(12, 28, 4)))+
                     # labels = c("8", "16", "19.25", "24", "32")) +#removes whitespace below bars! custom ticks
  # expand_limits(y = c(8, 34.5)) +
  ylab(expression(log[2]~(intensity))) + xlab(expression(paste("Bead slurry volume (", mu, "L", ")"))) #+
  # ggtitle("Intensity of pY sites")
 
plot_pY_intensity_distribution_reps


ggsave("output/SuppFig5/plot_pY_intensity_distribution_reps.png", plot = plot_pY_intensity_distribution_reps, width = 8, height = 8, scale = 0.4)

ggsave("output/SuppFig5/plot_pY_intensity_distribution_reps.pdf", plot = plot_pY_intensity_distribution_reps, width = 8, height = 8, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
psite_reps_CVs_YESnorm_pY <- my_qc_cvs(
  data = ascore_stringent_fasta_psite %>% ungroup() %>% filter(mod_res =="Y"), 
  grouping = ref, 
  condition = condition,
  intensity = raw_median_norm_intensity,
  plot = TRUE, 
  xlab = expression(paste("Bead slurry volume (", mu, "L", ")")),
  plot_style = "boxplot",
  max_cv = 100)

psite_reps_CVs_YESnorm_pY

ggsave(filename = "output/SuppFig5/CVs_YES_normalized_pY.png", plot = psite_reps_CVs_YESnorm_pY, width = 8, height = 10, scale = 0.4)
ggsave(filename = "output/SuppFig5/CVs_YES_normalized_pY.pdf", plot = psite_reps_CVs_YESnorm_pY, width = 8, height = 10, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
psite_correlation <- qc_sample_correlation(
  data = ascore_stringent_fasta_psite,
  sample = sample_id, 
  grouping = ref, 
  intensity_log2 = median_norm_intensity,
  condition = condition,
  # method = "pearson",
  interactive =  FALSE)

psite_correlation


ggsave(filename = "output/SuppFig5/psite_correlation.png", plot = psite_correlation, width = 14, height = 12, scale = 0.4)
ggsave(filename = "output/SuppFig5/psite_correlation.pdf", plot = psite_correlation, width = 14, height = 12, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
psite_correlation_pY <- qc_sample_correlation(
  data = ascore_stringent_fasta_psite %>% ungroup() %>% filter(mod_res =="Y"),
  sample = sample_id, 
  grouping = ref, 
  intensity_log2 = median_norm_intensity,
  condition = condition,
  # method = "pearson",
  interactive =  FALSE)

psite_correlation_pY


ggsave(filename = "output/SuppFig5/psite_correlation_pY.png", plot = psite_correlation_pY, width = 14, height = 12, scale = 0.4)
ggsave(filename = "output/SuppFig5/psite_correlation_pY.pdf", plot = psite_correlation_pY, width = 14, height = 12, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------



## --------------------------------------------------------------------------------------------------------------------------------------------------



## --------------------------------------------------------------------------------------------------------------------------------------------------



## --------------------------------------------------------------------------------------------------------------------------------------------------



## --------------------------------------------------------------------------------------------------------------------------------------------------



## --------------------------------------------------------------------------------------------------------------------------------------------------



## --------------------------------------------------------------------------------------------------------------------------------------------------



## --------------------------------------------------------------------------------------------------------------------------------------------------



## --------------------------------------------------------------------------------------------------------------------------------------------------



## --------------------------------------------------------------------------------------------------------------------------------------------------



## --------------------------------------------------------------------------------------------------------------------------------------------------



## --------------------------------------------------------------------------------------------------------------------------------------------------



## --------------------------------------------------------------------------------------------------------------------------------------------------



## --------------------------------------------------------------------------------------------------------------------------------------------------



## --------------------------------------------------------------------------------------------------------------------------------------------------



## --------------------------------------------------------------------------------------------------------------------------------------------------


