## --------------------------------------------------------------------------------------------------------------------------------------------------
alexis_theme <- function() {
  theme(
    # panel.border = element_rect(colour = "blue", fill = NA, linetype = 2),
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.x  = element_blank(),
    panel.grid.minor = element_blank(),
    # axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
    axis.text.x = element_text(angle = -90, vjust = 1, hjust = 0),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black", face = "plain", family = "sans"),
    axis.title = element_text(colour = "black", family = "sans"),
    axis.ticks = element_line(colour = "black"),
    plot.title = element_text(size=10),
    # legend at the bottom 6)
    legend.position = "right")   
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
Yeast_spurious <- read_csv(file = "raw_data/SuppFig7/Alex_yeast-and-spurious_20210219_forR.csv") %>%
  select(reference, full_protein_sequence) %>%
  mutate(organism = "Yeast") 



## --------------------------------------------------------------------------------------------------------------------------------------------------
##read in your data here, without ascore which filters for only phosphopeptides
comet <- read_csv("raw_data/SuppFig7/comet/result.csv") %>% 
  clean_names() %>%
  filter(grepl("gp", sample_name) == FALSE) %>% 
  filter(!raw_file_name %in% c("e07294", "e07281", "x09331", "x09271" ) ) %>% 
  mutate(condition = paste(raw_file_name, sample_name, sep = " "),
         replicate = 1) %>% 

  left_join(y = Yeast_spurious, by = c("reference")) %>% 
  
  ##rename column to 'intensity'
  rename(intensity = max_intensity_light_c2837) %>% 
  
  ##keep only forward hits
  filter(reverse == FALSE) %>% 
  
  ##classify phospho vs. no none
  mutate(phospho = case_when(
    grepl("\\@", sequence) ~ "phospho",
    TRUE~"none")) %>% 
  filter(replicate != 4) #dropped due to precipitate in old, and balanced replicates

  
  
  

## --------------------------------------------------------------------------------------------------------------------------------------------------
# sample_key <- comet %>% 
#   distinct(condition, replicate, raw_file_name)

sample_key <- comet %>% 
  distinct(raw_file_name, sample_name)


## --------------------------------------------------------------------------------------------------------------------------------------------------
ascore <- read_csv("raw_data/SuppFig7/ascore/result.csv") %>% 
  clean_names() %>%
  filter(grepl("gp", sample_name) == FALSE) %>%
  filter(!raw_file_name %in% c("e07294", "e07281", "x09331", "x09271" ) ) %>%
  mutate(condition =  paste(raw_file_name, sample_name, sep = " "),
         replicate = 1) %>%

  
  mutate(
    reference = case_when(
      reference == "P63185-V394A" ~ "P63185",
      reference == "P04049-D468N" ~ "P04049",
      TRUE ~ as.character(reference))) %>% 
  
  
  left_join(y = Yeast_spurious, by = c("reference")) %>% 
  
  ##rename column to 'intensity'
  rename(intensity = max_intensity_light_c2837) %>% 
  
  ##keep only forward hits
  filter(reverse == FALSE) %>% 
  
  ##classify phospho vs. no none
  mutate(phospho = case_when(
    grepl("\\@", sequence) ~ "phospho",
    TRUE~"none"))
  

# sample_key <- ascore %>% distinct(raw_file_name, condition, replicate)
# 
# sample_key


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
##relative to intensity of all peptides (with and without phospho)
set.seed(690)
plot_ratio_phos_count <- ggplot() + 
  geom_bar(data = (enrichment_efficiency_count_df %>%
                     distinct(condition, avg_phosphopeptide_enrichment_efficiency)),
           mapping = aes(x = condition, y = avg_phosphopeptide_enrichment_efficiency),
           stat = "identity",
           position = "dodge",
           fill = "grey30",
           alpha = 0.3, width = 0.8) +
  geom_point(data = enrichment_efficiency_count_df, mapping = aes(x = condition,
                 y = individual_phosphopept_enrich_efficiency,
                 fill = replicate),
             position = position_jitterdodge(dodge.width = 0.25), alpha = 1, size =1, shape = 1, show.legend =  FALSE)  +
  ylab("p- /all peptide (intensity)") +
  xlab ("bead age") +

  
  alexis_theme()+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.02), breaks = c(seq(0, 1, 0.25))) + #removes whitespace below bars!
  theme(axis.text.x = element_text(angle = -60, vjust = 0.1, hjust = 0.15, size = 10),
        axis.title = element_text(hjust = 0.5, size = 10))
  

plot_ratio_phos_count
ggsave("output/SuppFig7/plot_ratio_phos_count_ratio_style2.png", plot = plot_ratio_phos_count, width = 3, height = 6, scale = 0.4)
ggsave("output/SuppFig7/plot_ratio_phos_count_ratio_style2.pdf", plot = plot_ratio_phos_count, width = 4, height = 6, scale = 0.4)


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
set.seed(690)
plot_ratio_phos_intensity <- ggplot(data = (enrichment_efficiency_intensity_df)) + 
  geom_bar(aes(x = condition, y = mean_ratio),
           stat = "identity",
           position = "dodge",
           fill = "skyblue2",
           alpha = 0.3, width = 0.8) +
  geom_point(aes(x = condition,
                 y = rep_ratio,
                 fill = replicate),
             position = position_jitterdodge(dodge.width = 0.25), alpha = 1, size =1, shape = 1, show.legend =  FALSE)  +
  # theme(axis.text.x = element_text(angle=90, vjust = 0.25, hjust = 1), legend.position = "none") +
  # theme(legend.position = "none") +
  ylab("p- /all peptide (intensity)") +
  xlab ("bead age") +
  alexis_theme()+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.02), breaks = c(seq(0, 1, 0.25))) + #removes whitespace below bars!
  theme(axis.text.x = element_text(angle = -60, vjust = 0.1, hjust = 0.15, size = 10),
        axis.title = element_text(hjust = 0.5, size = 10))
  

plot_ratio_phos_intensity
ggsave("output/SuppFig7/phospho_pept_intensity_ratio_style2.png", plot = plot_ratio_phos_intensity, width = 3, height = 6, scale = 0.4)
ggsave("output/SuppFig7/phospho_pept_intensity_ratio_style2.pdf", plot = plot_ratio_phos_intensity, width = 4, height = 6, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
##relative to intensity of all peptides (with and without phospho)
set.seed(690)
plot_ratio_phos_intensity_oneperbeadlot <- ggplot(data = (enrichment_efficiency_intensity_df %>% 
                                                            mutate(raw_file_name = str_sub(condition, start = 1L, end = 6L)) %>% 
                                          filter(raw_file_name %in% c( "e09546", "e09595", "e13550", "e13551", "x06094", "x09014", "x09013", "x09116",  "x10147")) %>% 
                                            filter(raw_file_name %in% c( "e09546", "e09595", "e13550", "e13551", "x06094", "x09014", "x09013", "x09116",  "x10147")) %>% 
                                            mutate(raw_file_name = fct_relevel(raw_file_name, "e09546", "x06094","x10147", "e09595",  "x09013", "x09014", "x09116")) %>% 
                                            filter(!raw_file_name %in% c("e13550", "e13551")))) +
  geom_bar(aes(x = raw_file_name, y = mean_ratio),
           stat = "identity",
           position = "dodge",
           fill = "skyblue2",
           alpha = 0.3, width = 0.8) +
  geom_point(aes(x = raw_file_name,
                 y = rep_ratio,
                 fill = replicate),
             position = position_jitterdodge(dodge.width = 0.25), alpha = 1, size =1, shape = 1, show.legend =  FALSE)  +
  geom_hline(mapping = aes(yintercept = 0.98), linetype = 2, linewidth = 0.5, color = "grey60", alpha = 0.8) +
  annotate(geom = "text",x = 5, y = 0.93, label = ">98% purity", color = "grey60") +
  # theme(axis.text.x = element_text(angle=90, vjust = 0.25, hjust = 1), legend.position = "none") +
  # theme(legend.position = "none") +
  ylab("p- / all peptide (intensity)") +
  xlab ("bead lot") +
  alexis_theme()+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.02), breaks = c(seq(0, 1, 0.25))) + #removes whitespace below bars!
  theme(axis.text.x = element_text(angle = -80, vjust = 0.3, hjust = 0, size = 10),
        axis.title = element_text(hjust = 0.5, size = 10))
  

plot_ratio_phos_intensity_oneperbeadlot
ggsave("output/SuppFig7/phospho_pept_intensity_ratio_style2_oneperbeadlot.png", plot = plot_ratio_phos_intensity_oneperbeadlot, width = 8, height = 6, scale = 0.4)
ggsave("output/SuppFig7/phospho_pept_intensity_ratio_style2_oneperbeadlot.pdf", plot = plot_ratio_phos_intensity_oneperbeadlot, width = 6, height = 6, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
##relative to intensity of all peptides (with and without phospho)
set.seed(690)
plot_ratio_phos_intensity_oneperbeadlot4 <- ggplot(data = (enrichment_efficiency_intensity_df %>% 
                                                            mutate(raw_file_name = str_sub(condition, start = 1L, end = 6L)) %>% 
                                          filter(raw_file_name %in% c(  "e13550", "e13551",  "x09013", "x09014", "x09116")))) +
  geom_bar(aes(x = raw_file_name, y = mean_ratio),
           stat = "identity",
           position = "dodge",
           fill = "skyblue2",
           alpha = 0.3, width = 0.8) +
  geom_point(aes(x = raw_file_name,
                 y = rep_ratio,
                 fill = replicate),
             position = position_jitterdodge(dodge.width = 0.25), alpha = 1, size =1, shape = 1, show.legend =  FALSE)  +
  geom_hline(mapping = aes(yintercept = 0.98), linetype = 2, linewidth = 0.5, color = "grey60", alpha = 0.8) +
  annotate(geom = "text",x = 4, y = 0.93, label = ">98% purity", color = "grey60") +
  # theme(axis.text.x = element_text(angle=90, vjust = 0.25, hjust = 1), legend.position = "none") +
  # theme(legend.position = "none") +
  ylab("p- / all peptide (intensity)") +
  xlab ("bead lot") +
  alexis_theme()+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.02), breaks = c(seq(0, 1, 0.25))) + #removes whitespace below bars!
  theme(axis.text.x = element_text(angle = -80, vjust = 0.3, hjust = 0, size = 10),
        axis.title = element_text(hjust = 0.5, size = 10))
  

plot_ratio_phos_intensity_oneperbeadlot4
ggsave("output/SuppFig7/phospho_pept_intensity_ratio_style2_oneperbeadlot4.png", plot = plot_ratio_phos_intensity_oneperbeadlot4, width = 8, height = 6, scale = 0.4)
ggsave("output/SuppFig7/phospho_pept_intensity_ratio_style2_oneperbeadlot4.pdf", plot = plot_ratio_phos_intensity_oneperbeadlot4, width = 4, height = 6, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_psm_intensity_distribution_reps <- ggplot(data = comet) +
  geom_boxplot(mapping = aes(x = condition, y = log2(intensity), fill = as.factor(replicate)),
               outlier.shape = 21,
               outlier.alpha = 0.3,
               outlier.size = 0.4) +
  alexis_theme() +
  # theme(legend.position = "none") +
  scale_fill_brewer(palette = "Greys", name = "replicate") +
  ylab(expression(log[2]~(intensity))) + xlab("bead age") +
  theme(axis.title = element_text(hjust = 0.5, vjust = 1,  family = "sans", size = 14),
    axis.text.x = element_text(angle=0, vjust = 1, hjust = 0.5, family = "sans", size = 14),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, family = "sans", size = 12),
        legend.title = element_text(hjust = 0.5, face = "bold", family = "sans")) +
  theme(axis.text.x = element_text(angle = -60, vjust = 0.1, hjust = 0.15, size = 12),
        axis.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("Intensity of all PSMs")

plot_psm_intensity_distribution_reps


ggsave("output/SuppFig7/all_psm_intensity_distribution_reps.png", plot = plot_psm_intensity_distribution_reps, width = 6, height = 7, scale = 0.4)
ggsave("output/SuppFig7/all_psm_intensity_distribution_reps.pdf", plot = plot_psm_intensity_distribution_reps, width = 8, height = 6, scale = 0.4)


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
  left_join(y = Yeast_spurious, by = c("reference")) %>% 
  
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
write_csv(x = ascore_stringent_fasta_precursor %>% distinct(condition, replicate, ascore_sequence, charge, ascore, intensity, reference, sample_name,   mod_position_in_protein, ref), file = "modified_data/SuppFig7/distinct_pSTYpeptides_ascore13_replicate.csv", col_names = TRUE)


#condition
write_csv(x = ascore_stringent_fasta_precursor %>% distinct(condition,  ascore_sequence, charge, ascore, intensity, reference,   mod_position_in_protein, ref), file = "modified_data/SuppFig7/distinct_pSTYpeptides_ascore13_condition.csv", col_names = TRUE)


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
write_csv(x = distinct_p_sites_condition, file = "modified_data/SuppFig7/distinct_pSTYsites_ascore13_condition.csv", col_names = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_distinct_p_sites_condition <- ggplot(data = distinct_p_sites_condition) +
  geom_bar(mapping = aes(x = condition, fill = mod_res)) +
  theme_bw(12) +
  scale_fill_viridis_d() +
  ylab("unique phospho sites") +
  theme(axis.text.x = element_text(angle = -60, vjust = 0.9, hjust = 0))

plot_distinct_p_sites_condition

ggsave("output/SuppFig7/distinct_p_sites_condition.png", plot = plot_distinct_p_sites_condition, width = 16, height = 30, scale = 0.5)
ggsave("output/SuppFig7/distinct_p_sites_condition.pdf", plot = plot_distinct_p_sites_condition, width = 16, height = 30, scale = 0.5)


## --------------------------------------------------------------------------------------------------------------------------------------------------
summarize_distinct_psites_condition <- distinct_p_sites_condition %>%
  ungroup() %>% 
  group_by(condition, mod_res) %>% 
  summarize(
    n_psites = n())

summarize_distinct_psites_condition


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_distinct_p_sites_condition_one_per_bead_lot <- ggplot(data = distinct_p_sites_condition %>% 
                                            mutate(raw_file_name = str_sub(condition, start = 1L, end = 6L)) %>% 
                                          filter(raw_file_name %in% c( "e09546", "e09595", "e13550", "e13551", "x06094", "x09014", "x09013", "x09116",  "x10147"))) +
  geom_bar(mapping = aes(x = raw_file_name, fill = mod_res), show.legend = FALSE) +
  theme_bw(14) +
  scale_fill_viridis_d() +
  ylab("unique phospho sites") +
  theme(axis.text.x = element_text(angle = -80, vjust = 0.9, hjust = 0),
        axis.text.y = element_text(size = 18)) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(0, 8000, 1000)), limits = c(0, 8000))

plot_distinct_p_sites_condition_one_per_bead_lot

ggsave("output/SuppFig7/distinct_p_sites_condition_one_per_bead_lot.png", plot = plot_distinct_p_sites_condition_one_per_bead_lot, width = 8, height = 6, scale = 0.5)
ggsave("output/SuppFig7/distinct_p_sites_condition_one_per_bead_lot.pdf", plot = plot_distinct_p_sites_condition_one_per_bead_lot, width = 12, height = 20, scale = 0.5)


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_distinct_p_sites_condition_one_per_bead_lot_ordered <- ggplot(data = distinct_p_sites_condition %>% 
                                            mutate(raw_file_name = str_sub(condition, start = 1L, end = 6L)) %>% 
                                          filter(raw_file_name %in% c( "e09546", "e09595", "e13550", "e13551", "x06094", "x09014", "x09013", "x09116",  "x10147")) %>% 
                                            mutate(raw_file_name = fct_relevel(raw_file_name, "e09546", "x06094","x10147", "e09595",  "x09013", "x09014", "x09116")) %>% 
                                            filter(!raw_file_name %in% c("e13550", "e13551"))) +
  geom_bar(mapping = aes(x = raw_file_name, fill = mod_res), show.legend = FALSE) +
  alexis_theme() +
  scale_fill_viridis_d() +
  ylab("unique phospho sites") +
  theme(axis.text.x = element_text(angle = -80, vjust = 0.9, hjust = 0),
        axis.text.y = element_text(size = 18)) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(0, 8000, 1000)), limits = c(0, 8000))

plot_distinct_p_sites_condition_one_per_bead_lot_ordered

ggsave("output/SuppFig7/distinct_p_sites_condition_one_per_bead_lot_ordered.png", plot = plot_distinct_p_sites_condition_one_per_bead_lot_ordered, width = 8, height = 6, scale = 0.5)
ggsave("output/SuppFig7/distinct_p_sites_condition_one_per_bead_lot_ordered.pdf", plot = plot_distinct_p_sites_condition_one_per_bead_lot_ordered, width = 6, height = 6, scale = 0.5)


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_distinct_p_sites_condition_one_per_bead_lot_4 <- ggplot(data = distinct_p_sites_condition %>% 
                                            mutate(raw_file_name = str_sub(condition, start = 1L, end = 6L)) %>% 
                                          filter(raw_file_name %in% c(  "e13550", "e13551",  "x09013", "x09014", "x09116"))) +
  geom_bar(mapping = aes(x = raw_file_name, fill = mod_res), show.legend = FALSE) +
  alexis_theme() +
  scale_fill_viridis_d() +
  ylab("unique phospho sites") +
  theme(axis.text.x = element_text(angle = -80, vjust = 0.9, hjust = 0),
        axis.text.y = element_text(size = 18)) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(0, 4000, 1000)), limits = c(0, 4000))

plot_distinct_p_sites_condition_one_per_bead_lot_4

ggsave("output/SuppFig7/distinct_p_sites_condition_one_per_bead_lot4.png", plot = plot_distinct_p_sites_condition_one_per_bead_lot_4, width = 6, height = 6, scale = 0.4)
ggsave("output/SuppFig7/distinct_p_sites_condition_one_per_bead_lot4.pdf", plot = plot_distinct_p_sites_condition_one_per_bead_lot_4, width = 6, height = 6, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
distinct_p_sites_reps <- ascore_stringent_fasta %>% 
  ungroup() %>% 
  distinct(reference, mod_protein_location, mod_res,  condition, replicate) %>% 
  mutate(
    sample_id = paste(condition, replicate, sep = "_"))
##this collapses all  replicates into a single set of unique p sites.


## --------------------------------------------------------------------------------------------------------------------------------------------------
write_csv(x = distinct_p_sites_reps, file = "modified_data/SuppFig7/pSTYsites_ascore13_reps.csv", col_names = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_distinct_p_sites_reps <- ggplot(data = distinct_p_sites_reps) +
  geom_bar(mapping = aes(x = sample_id, fill = mod_res)) +
  theme_bw(18) +
  scale_fill_viridis_d() +
  ylab("unique phospho sites") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1, size = 6))

plot_distinct_p_sites_reps

ggsave("output/SuppFig7/distinct_p_sites_reps.png", plot = plot_distinct_p_sites_reps, width = 15, height = 15, scale = 0.5)
ggsave("output/SuppFig7/distinct_p_sites_reps.pdf", plot = plot_distinct_p_sites_reps, width = 15, height = 15, scale = 0.5)


## --------------------------------------------------------------------------------------------------------------------------------------------------
summarize_distinct_psites_replicates <- distinct_p_sites_reps %>%
  ungroup() %>% 
  group_by(condition, replicate, sample_id, mod_res) %>% 
  summarize(
    n_psites = n())

summarize_distinct_psites_replicates


## --------------------------------------------------------------------------------------------------------------------------------------------------
distinct_p_sites_total <- ascore_stringent_fasta %>% 
  distinct(reference, mod_position_in_protein, mod_res) %>%
  mutate(filter_level = "Ascore >= 13")


## --------------------------------------------------------------------------------------------------------------------------------------------------
distinct_pY_sites_total_summary <- distinct_p_sites_total %>% 
  group_by(mod_res) %>% 
  summarize(
    n_distinct_sites = n()) %>% 
  ungroup()

distinct_pY_sites_total_summary


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_distinct_p_sites_total <- ggplot(data = distinct_p_sites_total) +
  geom_bar(mapping = aes(x = filter_level, fill = mod_res)) +
  theme_bw(18) +
  scale_fill_viridis_d() +
  ylab("unique phospho sites") +
  alexis_theme()+
  scale_y_continuous(expand = c(0,0), limits = c(0, 5000), breaks = c(seq(0, 5000, 1000))) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5))

plot_distinct_p_sites_total

ggsave("output/SuppFig7/distinct_p_sites_total.png", plot = plot_distinct_p_sites_total, width = 6, height = 8, scale = 0.4)
ggsave("output/SuppFig7/distinct_p_sites_total.pdf", plot = plot_distinct_p_sites_total, width = 8, height = 8, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
distinct_p_sites_condition <- ascore_stringent_fasta %>% 
  distinct(reference, mod_protein_location, mod_res, condition)
##this collapses all  replicates into a single set of unique p sites.


## --------------------------------------------------------------------------------------------------------------------------------------------------
write_csv(x = distinct_p_sites_condition, file = "modified_data/SuppFig7/pSTYsites_ascore13_condition.csv", col_names = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_distinct_p_sites_condition <- ggplot(data = distinct_p_sites_condition ) +
  geom_bar(mapping = aes(x = condition, fill = mod_res)) +
  theme_bw(18) +
  scale_fill_viridis_d(name = "mod\nsite") +
  ylab("unique phospho sites") +
  alexis_theme()+
  scale_y_continuous(expand = c(0,0), limits = c(0, 5000), breaks = c(seq(0, 5000, 1000))) +
  xlab("peptide concentration \n (0.5 mg input)")+
  theme(axis.text.x = element_text(angle = -60, vjust = 0.15, hjust = 0.1, ),
        axis.title.x = element_text(size = 15))

plot_distinct_p_sites_condition

ggsave("output/SuppFig7/distinct_p_sites_condition.png", plot = plot_distinct_p_sites_condition, width = 10, height = 12, scale = 0.4)
ggsave("output/SuppFig7/distinct_p_sites_condition.pdf", plot = plot_distinct_p_sites_condition, width = 8, height = 8, scale = 0.4)


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
distinct_p_sites_reps <- ascore_stringent_fasta %>% 
  ungroup() %>% 
  mutate(sample_id = paste(condition, replicate, sep = "\n")) %>% 
  distinct(condition, replicate, sample_id, reference, mod_protein_location, mod_res)
##this collapses all  replicates into a single set of unique p sites.


## --------------------------------------------------------------------------------------------------------------------------------------------------
write_csv(x = distinct_p_sites_reps, file = "modified_data/SuppFig7/pSTYsites_ascore13_reps.csv", col_names = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------
plot_distinct_p_sites_reps <- ggplot(data = distinct_p_sites_reps) + 
  geom_bar(mapping = aes(x = sample_id, fill = mod_res)) +
  theme_bw(10) +
  scale_fill_viridis_d() +
  ylab("unique phospho sites")  +
  xlab("1 mg peptide input") +
  alexis_theme()+
 
  theme(axis.text.x = element_text(angle = -60, vjust = 0.75, hjust = 0))

plot_distinct_p_sites_reps

ggsave("output/SuppFig7/distinct_p_sites_reps.png", plot = plot_distinct_p_sites_reps, width = 8, height = 10, scale = 0.4)
ggsave("output/SuppFig7/distinct_p_sites_reps.pdf", plot = plot_distinct_p_sites_reps, width = 8, height = 10, scale = 0.4)


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
  geom_bar(data = distinct_p_sites_condition,
             
                   mapping = aes(x = condition, fill = mod_res),show.legend = FALSE) +
  
  geom_errorbar(data = (df_geom_errorbar), 
                mapping = aes(x = condition, ymin = min_psites, ymax = max_psites),  show.legend = FALSE, linewidth = 0.5, width = 0.25, color = "grey90") +
  
  
  geom_jitter(data = (distinct_p_sites_reps), 
              mapping = aes(x = condition, color = sample_id, fill = mod_res),
              stat = "count", shape = 1, show.legend = FALSE, size = 1, stroke =0.5,
              width = 0.35, height = 0, alpha = 0.5) +
   
  theme_bw(10) +
  facet_wrap(facets = vars(mod_res), ncol = 3) +
  scale_fill_viridis_d() +
  scale_color_manual(values = rep(c("grey70" ),60)) +
  ylab("unique phospho sites")   +
  alexis_theme()+
  scale_y_continuous(expand = c(0,0), limits = c(0, 8000), breaks = c(seq(0, 8000, 1000)))

plot_distinct_p_sites_condition_reps_wPTS

ggsave("output/SuppFig7/distinct_p_sites_condition_wPTS.png", plot = plot_distinct_p_sites_condition_reps_wPTS, width = 6, height = 12, scale = 0.4)
ggsave("output/SuppFig7/distinct_p_sites_condition_wPTS.pdf", plot = plot_distinct_p_sites_condition_reps_wPTS, width = 6, height = 12, scale = 0.4)


## --------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(1333)
plot_distinct_p_sites_condition_reps_wPTS <- ggplot() +
  geom_bar(data = (distinct_p_sites_condition %>%
                     filter(mod_res == "Y")),
                   mapping = aes(x = condition, fill = mod_res),
           show.legend = FALSE, color = "black", linewidth = 0.5, width = 0.75) +
  
   geom_errorbar(data = (df_geom_errorbar %>%
                           filter(mod_res == "Y")),
                mapping = aes(x = condition, ymin = min_psites, ymax = max_psites),  show.legend = FALSE, linewidth = 0.5, width = 0.2, color = "grey30") +
  
  
  geom_jitter(data = (distinct_p_sites_reps%>%
                        filter(mod_res == "Y")),
              mapping = aes(x = condition, color = sample_id, fill = mod_res),
              stat = "count", shape = 21, show.legend = FALSE, size = 1.2, stroke =0.5, fill = "white",
              width = 0.1, height = 0, alpha = 0.6) +
  geom_text(data = tibble(x = c(1, 2),y = c(4140, 4192)), mapping = aes( x = x, y = y +150, label = y),
            color = "black",
            size = 3,
            fontface = "bold",
            inherit.aes = FALSE) +
 
  scale_fill_viridis_d(direction = -1) +
  scale_color_manual(values = rep(c("grey30" ),60)) +
  ylab("unique pY sites")   +
  xlab("condition") +
  alexis_theme()+
  scale_y_continuous(expand = c(0,0), limits = c(0, 5000), breaks = c(seq(0, 5000, 1000))) +
  theme(axis.text.x = element_text(angle = -60, vjust = 0.75, hjust = 0, size = 10),
        axis.text.y = element_text( size = 10))

plot_distinct_p_sites_condition_reps_wPTS

ggsave("output/SuppFig7/distinct_pY_sites_condition_wPTS.png", plot = plot_distinct_p_sites_condition_reps_wPTS, width = 4, height = 6, scale = 0.4)
ggsave("output/SuppFig7/distinct_pY_sites_condition_wPTS.pdf", plot = plot_distinct_p_sites_condition_reps_wPTS, width = 6, height = 6, scale = 0.4)

