# merge all data for all studies

library(tidyverse)
library(DECIPHER)
library(pheatmap)
library(phyloseq)
library(microbiome)
library(MatrixGenerics)
library("optparse")
library(cowplot)
library(vegan)
library(purrr)
library(future)
library(gridExtra)
library(viridis)
library(ggpubr)
library(dplyr)

source("./dada_2_functions.R")
source("./setlists.R")

add_asv_alpha_div <- function(ps_object) {
  ps_object <- subset_taxa(ps_object, genus != "Mitochondria")
  if(min(sample_sums(ps_object)) < 2000) {
    ps_low_counts <- ps_object %>% prune_samples(samples = (sample_sums(.) < 2000))
    ps_raref <- ps_object %>% rarefy_even_depth(., sample.size = 2000) %>%
      merge_phyloseq(.,ps_low_counts)
  } else {
    ps_raref <- ps_object %>% rarefy_even_depth(., sample.size = 2000)
  }
  matching_taxa <- ps_raref@tax_table %>% data.frame() %>% filter(grepl("unclassified", genus)) %>% rownames()
  raref_otu_table <- ps_raref@otu_table %>% data.frame()
  unassigned_taxa_present_per_sample_raref <- apply(raref_otu_table[,matching_taxa], 1, function(x) sum(x > 0))
  ps_raref@sam_data$n_unassigned_asvs_sample_raref <- unassigned_taxa_present_per_sample_raref
  alpha.div_gen <- microbiome::alpha(ps_raref, index = c("Observed", "Shannon"))
  ps_raref@sam_data$shannon_asv <- alpha.div_gen$diversity_shannon
  ps_raref@sam_data$richness_asv <- alpha.div_gen$observed
  mdat_raref <- ps_raref@sam_data %>% data.frame() %>% select(shannon_asv, richness_asv, 
                                                              n_unassigned_asvs_sample_raref,
                                                              run_accession)
  ps_object@sam_data <- ps_object@sam_data %>% data.frame() %>%
    mutate(shannon_asv = NULL,
           richness_asv = NULL, 
           n_unassigned_asvs_sample_raref = NULL) %>%
    left_join(., mdat_raref, by = "run_accession") %>%
    mutate(rownames_1 = run_accession) %>%
    column_to_rownames("rownames_1") %>%
    sample_data()
  return(ps_object)
}


merged_set <- "all"

merging_step <- F
rarecurve_step <- F
diversity_step <- F
ordination_step <- F
permanova_step <- F


n_cores <- 2

option_list = list(
  make_option(c("-t", "--threads"), type="numeric", default=NULL))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$threads)){
  n_cores <- T
} else {
  n_cores <- opt$threads
}


# read phyloseq objects
# study_names <- c("kristensen_2020", "roswall_2021", "hill_2017", "wampach_2018")#, "reyman_2019", "hesla_2014", "iszatt_2019", "raman_2019")
study_names <- set_list[[merged_set]]
names(study_names) <- study_names
# load phyloseq objects
ps_list <- lapply(study_names, function(x){readRDS(list.files(paste("/fast/AG_Forslund/rob/mm_index/study_data", x, sep = "/"),
                                                   pattern =  "phyloseq.*\\d\\.rds",
                                                   full.names = T))})
# lapply(study_names, function(x){list.files(paste("../study_data", x, sep = "/"),
#                                                    pattern =  "phyloseq.*\\.rds",
#                                                    full.names = T)})
ps_list <- lapply(study_names, function(x) {ps_list[[x]]@phy_tree <- NULL; return(ps_list[[x]])})
lapply(ps_list, function(x) {x %>% tax_table() %>% nrow()}) %>% unlist %>% sum() # total taxa number
lapply(ps_list, function(x) {x %>% sample_data() %>% nrow()}) %>% unlist %>% sum() # total samples number

# get ASV level richness
ps_list_asv_div <- lapply(ps_list, add_asv_alpha_div)
# change ASV names to be compatible with BLAST output
ps_list_names <- lapply(ps_list_asv_div, function(x){taxa_names(x) <- paste0(unique(x@sam_data$study), "_", taxa_names(x)); return(x)})
rm(ps_list, ps_list_asv_div)
options(future.globals.maxSize = 25000 * 1024^2)

# data.frame(second_names = test_names, first_names = taxa_names(test)) %>% View()
start_time <- Sys.time()
print(names(ps_list_names)) %>% sort

if(merging_step){
  ps_list_genus <- lapply(ps_list_names, function(x) subset_taxa(x, genus != "Mitochondria") %>% aggregate_taxa(., level = "genus"))
  ps_list_family <- lapply(ps_list_names, function(x) subset_taxa(x, genus != "Mitochondria") %>% aggregate_taxa(., level = "family"))
  ps_list_class <- lapply(ps_list_names, function(x) subset_taxa(x, genus != "Mitochondria") %>% aggregate_taxa(., level = "class"))
  final_ps_genus <- merge_multiple_clusters(ps_list = ps_list_genus,
                                             consensus_seq = F,
                                             cores = n_cores,
                                             simple_merge = T)
  final_ps_family <- merge_multiple_clusters(ps_list = ps_list_family,
                                            consensus_seq = F,
                                            cores = n_cores,
                                            simple_merge = T)
  final_ps_class <- merge_multiple_clusters(ps_list = ps_list_class,
                                             consensus_seq = F,
                                             cores = n_cores,
                                             simple_merge = T)
  saveRDS(final_ps_genus, paste0("../data/", merged_set, "_phyloseq_raw_genus.rds"))
  saveRDS(final_ps_family, paste0("../data/", merged_set, "_phyloseq_raw_family.rds"))
  saveRDS(final_ps_class, paste0("../data/", merged_set, "_phyloseq_raw_class.rds"))
} else {
  final_ps_genus <- readRDS(paste0("../data/", merged_set, "_phyloseq_raw_genus.rds"))
  final_ps_family <- readRDS(paste0("../data/", merged_set, "_phyloseq_raw_family.rds"))
  final_ps_class <- readRDS(paste0("../data/", merged_set, "_phyloseq_raw_class.rds"))
}
# final_ps_class@sam_data <- final_ps_class@sam_data %>%
#   data.frame() %>% select(subject_ID, sample_ID, run_accession, 
#                           lifestyle, age, country, Instrument, X,
#                           geographic_location_.latitude.,
#                           geographic_location_.longitude.,
#                           n_unassigned_asvs_study, n_total_asvs_study,
#                           n_unassigned_asvs_sample, n_total_asvs_sample,
#                           shannon_asv, richness_asv, n_unassigned_asvs_sample_raref,
#                           read_count, study_accession, study, birthmode, health) %>%
#   sample_data()

print(sort(final_ps_genus@sam_data$study %>% unique))
end_time <- Sys.time()
print("Total time for merging step:")
print(end_time - start_time)
rm(ps_list_names)



# rarefaction analysis and filtering ###########################################
otu_matrix <- t(otu_table(final_ps_genus))
class(otu_matrix) <- "matrix"
# raredat <- rarecurve(otu_matrix, step = 10)

start_time <- Sys.time()
if(rarecurve_step) {
  raredat <- rarecurve(otu_matrix, step = 50, tidy = T)
  sample_data(final_ps_genus)$sample_sum <- sample_sums(final_ps_genus)
  sample_data(final_ps_genus) <- sample_data(final_ps_genus) %>%
    data.frame %>%
    mutate(age_cat = cut(age, breaks=c(-1, 7, 30, 60, 180, 365, 730, Inf),
                         labels=c("0-7", "7-30", "30-60", "60-180", "180-365", "365-730", "over 730")),
           cutoff = case_when(study == "blanton_2016" ~ 5000,
                              study == "bockulich_2016" ~ 2500,
                              study == "vatanen_2018" ~ 2500,
                              study == "roswall_2021" ~ 5000,
                              study == "hill_2017" ~ 5000,
                              study == "sprockett_2020" ~ 5000,
                              study == "pannaraj_2017" ~ 1000,
                              study == "lim_2015" ~ 2500,
                              study == "subramanian_2014" ~ 700,
                              study == "wampach_2018" ~ 2000,
                              study == "gehrig_2019" ~ 2000,
                              study == "kristensen_2020" ~ 2500,
                              study == "stokholm_2018" ~ 2000,
                              study == "raman_2019" ~ 2500,
                              study == "reyman_2019" ~ 5000,
                              study == "kamngona_2019" ~ 5000,
                              study == "kortekangas_2020" ~ 5000,
                              study == "muinck_2018" ~ 10000,
                              study == "beller_2021" ~ 2500,
                              study == "bender_2016" ~ 5000,
                              study == "davis_2017" ~ 2500,
                              study == "morandini_2023" ~ 2000)) %>%
    sample_data()
  
  raredat <- left_join(raredat , sample_data(final_ps_genus), by = c("Site" = "run_accession")) %>%
    select(Site, Sample, Species, X, subject_ID, sample_ID, age, country, study, 
           lifestyle, sample_sum, age_cat, cutoff)
  
  n_reads <- data.frame(n_reads = sample_sums(final_ps_genus), 
                        names = sample_names(final_ps_genus), 
                        study = final_ps_genus@sam_data$study,
                        age_cat = final_ps_genus@sam_data$age_cat,
                        cutoff = final_ps_genus@sam_data$cutoff)
  save(raredat, n_reads, file = "../data/rarefaction_curves.RData")
    
} else {
  load("../data/rarefaction_curves.RData")
}
rm(otu_matrix)
end_time <- Sys.time()
print("Total time for rarecurve step:")
print(end_time - start_time)


p1 <- ggplot(raredat %>% filter(study %in% c("raman_2019", "pannaraj_2017"))) +
  theme_bw() +
  # scale_fill_continuous() +
  # geom_density(data = n_reads, aes(x = n_reads)) +
  scale_color_discrete(guide = "none") +
  geom_line(aes(x = Sample, y = Species, group = Site, color = subject_ID)) +
  # geom_line(aes(x = Sample, y = dens)) +
  geom_vline(data = plyr::ddply(raredat %>% filter(study %in% c("raman_2019", "pannaraj_2017")), c("age_cat", "study"),
                                summarize, min_sample_size = min(sample_sum)),
             aes(xintercept=min_sample_size)) +
  geom_vline(aes(xintercept = cutoff), linetype = "dashed", color = "red") +
  facet_grid(age_cat ~ study) +
  xlim(0, 25000) +
  # ylim(0,500) +
  labs(title="Rarefaction curves") + xlab("Sequenced Reads") + ylab('ASVs Detected')


p2 <- ggplot(n_reads %>% filter(study %in% c("raman_2019", "pannaraj_2017")), aes(x=n_reads, color = age_cat)) +
  theme_bw() +
  geom_freqpoly(binwidth = 250) +
  xlim(0, 25000) +
  scale_y_log10()+
  geom_vline(aes(xintercept = cutoff), linetype = "dashed", color = "red") +
  facet_grid(age_cat~study) +
  ylab("N_samples")
p <- plot_grid(p1, p2, rel_heights = c(7/8, 1/8), nrow = 2, ncol = 1)
ggsave(paste0("../merged_data/all/figures/", merged_set, "_rarefaction_genus_plot.pdf"),
       device = "pdf", dpi = 320, plot = p, width = 49, height = 20)

rm(raredat, p1, p2, p)
gc()
# filter samples after rarefaction analysis
# filter criteria for each study: wampach_2018: 5.000, the rest 10.000, hill_2017: 15.00
final_ps_genus <- final_ps_genus %>% subset_samples(sample_sum >= cutoff)
final_ps_family <- final_ps_family %>% subset_samples(sample_sum >= cutoff)

# remove missing values:
final_ps_genus <- final_ps_genus %>% subset_samples(!is.na(age))
final_ps_family <- final_ps_family %>% subset_samples(!is.na(age))

saveRDS(final_ps_genus, paste0("../data/", merged_set, "_phyloseq_rf_filter_genus.rds"))
saveRDS(final_ps_family, paste0("../data/", merged_set, "_phyloseq_rf_filter_family.rds"))
