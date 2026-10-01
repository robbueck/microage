# plot PCoAs
version
setwd("/fast/AG_Forslund/rob/mm_index/R_scripts/")
library(tidyverse)
library(phyloseq)
library(microbiome)
library("optparse")
library(cowplot)
library(vegan)
library(purrr)
library(future)
library(gridExtra)
library(viridis)
library(ggpubr)
library(cowplot)
library("reshape2")
library(lme4)
library(lmtest)
library(usedist)
library(parallel)
library(rstatix)
library(ggprism)
library(gamlss)
library(gratia)
library(gamlss.ggplots)
library(ape)
library(furrr)



source("/fast/AG_Forslund/rob/mm_index/R_scripts/pcoa_functions.R")
source("/fast/AG_Forslund/rob/mm_index/R_scripts/setlists.R")
options(future.globals.maxSize = 25000 * 1024^2)

# Switches #######################
ordination_step <- T
ordination_genus_step <- T
ordination_genus_rare_step <- T
adonis_genus_step_1 <- F
adonis_genus_step_2 <- F
adonis_genus_step_3 <- F
adonis_genus_step_4 <- F
resampling_step <- F

n_cores <- 5

option_list = list(
  make_option(c("-t", "--threads"), type="numeric", default=NULL))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$threads)){
  n_cores <- T
} else {
  n_cores <- opt$threads
}
# # use rtk to extrapolate alpha diversity and counts for the low read count samples
# raw_ps <- readRDS("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_phyloseq_rf_filter_family.rds")
# raw_ps <- subset_samples(raw_ps, age <= 730) %>%
#   prune_samples(samples = (sample_sums(.) != 0)) %>%
#   subset_samples(., age > 0)
# otu_matrix <- raw_ps@otu_table@.Data
# 
# ps_low_counts <- raw_ps %>% prune_samples(samples = (sample_sums(.) < 2000))
# ps_object_family_raref <- raw_ps %>% rarefy_even_depth(., sample.size = 2000) %>%
#   merge_phyloseq(.,ps_low_counts)
# 
# # remove low abundant taxa:
# remove_taxa <- ps_object_family_raref %>%
#   prevalence(., detection = 10)
# remove_taxa <- names(remove_taxa[remove_taxa < 0.0002])
# 
# ps_object_family_raref <- ps_object_family_raref %>%
#   merge_taxa2_fixed(.,taxa = remove_taxa, name = "Other")
# 
# ps_object_family_comp <- raw_ps %>%
#   merge_taxa2_fixed(.,taxa = remove_taxa, name = "Other") %>%
#   microbiome::transform(transform = "compositional")
# 
# 
# ps_object_family_raref@sam_data %>%
#   ggplot(., aes(x = age, color = lifestyle)) +
#   geom_density() +
#   theme(
#     axis.title.x = element_text(size = 16),
#     axis.text.x = element_text(size = 14),
#     axis.title.y = element_text(size = 16),
#     axis.text.y = element_text(size = 14),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#     legend.text = element_text(size = 15),
#     legend.title = element_text(size = 15))
# ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/age_confounding.pdf", width = 11, height = 5)
# 
# 
# all_studies <- set_list$all
# studies_here <- unique(ps_object_family_raref@sam_data$study) %in% all_studies
# 
# 
# a_div_s <- diversity_plots(ps_object = ps_object_family_raref, x = "age", colour = "lifestyle",
#                 measures = c("diversity_shannon"), alpha = 0.3, size = 0.5) +
#   # scale_color_manual(values = turbo(length(all_studies))[studies_here], breaks = all_studies[studies_here] )+
#   ylab("Shannon diversity (Famliy)") +
#   xlab("age [days]") +
#   guides(color=guide_legend(override.aes = list(size = 4, alpha = 1))) +
#   geom_smooth() +
#   stat_cor(method = "spearman",
#            # aes(colour = NULL),
#            cor.coef.name = "spearman",
#            size = 6,
#            show.legend = F) +
#   theme(
#     axis.title.x = element_text(size = 16),
#     axis.text.x = element_text(size = 14),
#     axis.title.y = element_text(size = 16),
#     axis.text.y = element_text(size = 14),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#     legend.text = element_text(size = 15),
#     legend.title = element_text(size = 15))
# ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/a_div_family_present.pdf", width = 11, height = 6, plot = a_div_s)
# ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/a_div_family_present.png", width = 11, height = 6, plot = a_div_s)
# 
# 
# 
# 
# lm_full <- lm(value ~ age + lifestyle + study, data = a_div_s$data)
# anova(lm_full) %>% rstatix::eta_squared()



# alpha div genus data ##############################################################
ps_object_genus_raw <- readRDS("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_phyloseq_rf_filter_genus.rds")
ps_object_genus_raw <- subset_samples(ps_object_genus_raw, age <= 730 & age > 1) %>%
  prune_samples(samples = (sample_sums(.) != 0)) 

# filter as in RF modelling:
ps_object_comp_rf_comp <- ps_object_genus_raw %>%
  microbiome::transform(transform = "compositional") %>%
  filter_taxa(function(x) mean(x[x > 0]) > 5e-5, TRUE) %>%  # mean abundance cutoff, only counting abundances > 0
  prune_samples(samples = (sample_sums(.) != 0))
ps_object_rf <- ps_object_genus_raw %>% prune_taxa(taxa_names(ps_object_comp_rf_comp),.)
diversities <- estimate_richness(ps_object_rf, measures = c("Observed", "Shannon"))
ps_object_rf@sam_data$diversity_shannon_rf <- diversities$Shannon
ps_object_rf@sam_data$diversity_observed_rf <- diversities$Observed

mdat_rf <- ps_object_rf@sam_data %>% data.frame()
# does sequencing depth affect diversity in unrarefied data?
mdat_rf %>%
  # group_by(study, lifestyle) %>%
  # summarize(sample_sum = mean(sample_sum),
  #           diversity_shannon = mean(diversity_shannon)) %>%
  ggplot(., aes(x = sample_sum, y = diversity_shannon_rf, color = lifestyle)) +
  geom_point(, size = 0.1)
mdat_rf %>%
  ggplot(., aes(x = age, y = log(sample_sum), color = lifestyle)) +
  geom_point(, size = 0.1)

mdat_rf %>% select(sample_sum, diversity_shannon_rf) %>%
  cor(., method = "spearman")
mdat_rf %>%
  filter(lifestyle == "industrialized") %>%
  select(sample_sum, diversity_shannon_rf) %>%
  cor(., method = "spearman")
mdat_rf %>%
  filter(lifestyle == "non_industrialized") %>%
  select(sample_sum, diversity_shannon_rf) %>%
  cor(., method = "spearman")
# alpha-diversity and sequencing depth do not correlate, saturation is achieved
mdat_rf %>% select(sample_sum, age) %>%
  cor(., method = "spearman")
# read depth correlates slightly negatively with age
# age correlates strongly positively with shannon and not sample_sum, so probably read depth has no important effect on alpha diversity


otu_matrix <- ps_object_genus_raw@otu_table@.Data

ps_low_counts <- ps_object_genus_raw %>% prune_samples(samples = (sample_sums(.) < 2000))
ps_object_genus_raref <- ps_object_genus_raw %>% rarefy_even_depth(., sample.size = 2000) %>%
  merge_phyloseq(.,ps_low_counts)

# remove low abundant taxa:
remove_taxa_genus <- ps_object_genus_raw %>%
  prevalence(., detection = 10)
remove_taxa_genus <- names(remove_taxa_genus[remove_taxa_genus < 5e-5])

ps_object_genus_raref <- ps_object_genus_raref %>%
  merge_taxa2_fixed(.,taxa = remove_taxa_genus, name = "Other")

ps_object_genus_raw <- ps_object_genus_raw %>%
  merge_taxa2_fixed(.,taxa = remove_taxa_genus, name = "Other") %>%
  prune_taxa(taxa_names(ps_object_genus_raref),.)

ps_object_genus_comp <- ps_object_genus_raw %>%
  microbiome::transform(transform = "compositional")

ps_object_genus_clr <- ps_object_genus_raw %>%
  microbiome::transform(transform = "clr")

# check taxa of interest
# selected_taxa <- c("Treponema", "Succinivibrio", "Prevotella",
#                    "Lactobacillus", "Bifidobacterium")
# 
# abundance_data <- ps_object_genus_comp %>% 
#   psmelt() %>%
#   select(lifestyle, Abundance, OTU, age) %>%
#   filter(OTU %in% selected_taxa) %>%
#   mutate(Lifestyle = ifelse(lifestyle == "industrialized", yes = "Industrialized", no = "Non-industrialized"),
#          month = ceiling(age / 30)) 
# abundance_data %>%
#     ggplot(., aes(x = age, y = Abundance, color = lifestyle, fill = lifestyle)) +
#     geom_point(alpha = 0.1, size = 0.3) +
#     geom_smooth(method = "loess", alpha = 0) +
#     facet_wrap(~OTU, scales = "free_y") +
#     xlab("Age [Months]") +
#     ylab("Prevalence")+
#     ggtitle("B)")
    


ps_object_genus_raref@sam_data$lifestyle = factor(ps_object_genus_raref@sam_data$lifestyle, levels = c("non_industrialized", "industrialized"))
alpha.div_gen <- microbiome::alpha(ps_object_genus_raref, index = c("Observed", "Shannon")) %>%
  rownames_to_column()
write.csv(alpha.div_gen, "/fast/AG_Forslund/rob/mm_index/merged_data/all/a_divs.csv")
meta_df <- sample_data(ps_object_genus_raref) %>% 
  data.frame()

meta_df <- left_join(meta_df, alpha.div_gen, by = c("run_accession" = "rowname"))
a_div_s_genus <- ggplot(meta_df,# %>% filter(study %in% c("bender_2016")),
                        aes(y = diversity_shannon, x = age, color = lifestyle)) +
  geom_point(alpha = 0.4, size = 0.3) +
  ylab("Shannon diversity (Genus)") +
  xlab("age [days]") +
  geom_smooth(method = "loess") +
  # stat_cor(# aes(colour = NULL),
  #          method = "spearman",
  #          cor.coef.name = "spearman",
  #          size = 6,
  #          show.legend = F) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)) 
ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/a_div_genus_present.pdf",
       width = 11, height = 6, plot = a_div_s_genus)
ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/a_div_genus_present.png", 
       width = 11, height = 6, plot = a_div_s_genus)

# does sequencing depth influence diversity?
meta_df %>%
  # group_by(study, lifestyle) %>%
  # summarize(sample_sum = mean(sample_sum),
  #           diversity_shannon = mean(diversity_shannon)) %>%
  ggplot(., aes(x = age, y = sample_sum, color = lifestyle)) +
  geom_point(, size = 0.1) 
  # xlim(0, 100000) +
  geom_smooth(method = "lm")
meta_df %>% select(sample_sum, diversity_shannon) %>%
  cor(., method = "spearman")
meta_df %>%
  filter(lifestyle == "industrialized") %>%
  select(sample_sum, diversity_shannon) %>%
  cor(., method = "spearman")
meta_df %>%
  filter(lifestyle == "non_industrialized") %>%
  select(sample_sum, diversity_shannon) %>%
  cor(., method = "spearman")
# alpha-diversity and sequencing depth do not correlate, saturation is achieved
meta_df %>% select(sample_sum, age) %>%
  cor(., method = "spearman")
# read depth correlates slightly negatively with age
# age correlates strongly positively with shannon, so probably read depth has no important effect on alpha diversity
left_join(meta_df, mdat_rf %>% select(run_accession, diversity_shannon_rf), by = "run_accession") %>%
  ggplot(.,aes(x = diversity_shannon, y = diversity_shannon_rf)) +
  geom_point() +
  facet_wrap(~lifestyle)


## study, individual, sample rarefaction ##############
# studies
print("study rarefaction")
subsample_fn <- function(n_cat, df, cat = "study") {
  stds <- sample(unique(df[[cat]]), n_cat)
  df %>% filter(.data[[cat]] %in% stds) %>%
    select(-cat) %>% 
    colSums(. > 10) %>% 
    `>`(1) %>%
    sum()
}

mdat_all <- data.frame(ps_object_genus_raref@sam_data) %>%
  select(run_accession, lifestyle, study, subject_ID, sample_ID, age, sample_sum) %>%
  mutate(Lifestyle = case_when(lifestyle == "industrialized" ~ "Industrialized",
                               lifestyle == "non_industrialized" ~ "Non-industrialized"),
         .keep = "unused")
         
otu_table_lfst <- data.frame(ps_object_genus_raref@otu_table %>% t()) %>%
  rownames_to_column("run_accession") %>%
  left_join(mdat_all) %>%
  column_to_rownames("run_accession")

n_taxa_industrialized_studies_rarefy <- map_dfr(rep(1:13, each = 10), function(x) {
  otu_table_lfst %>%
    filter(Lifestyle == "Industrialized") %>%
    select(-Lifestyle, -subject_ID, -sample_ID, -age) %>%
    subsample_fn(x, ., cat = "study") %>%
    data.frame(n_taxa = ., n_studies = x, Lifestyle = "Industrialized")
})

n_taxa_non_industrialized_studies_rarefy <- lapply(rep(1:8, each = 10), function(x) {
  otu_table_lfst %>%
    filter(Lifestyle == "Non-industrialized") %>%
    select(-Lifestyle, -subject_ID, -sample_ID, -age) %>%
    subsample_fn(x, ., cat = "study") %>%
    data.frame(n_taxa = ., n_studies = x, Lifestyle = "Non-industrialized")
})

bind_rows(n_taxa_industrialized_studies_rarefy, n_taxa_non_industrialized_studies_rarefy) %>% 
  ggplot(.,aes(x = n_studies, y = n_taxa, color = Lifestyle)) +
  geom_point() +
  geom_smooth()
# gamlss-test:
n_taxa_studies_full <- gamlss(n_taxa ~ pb(n_studies) + Lifestyle, 
       data = bind_rows(n_taxa_industrialized_studies_rarefy,
                        n_taxa_non_industrialized_studies_rarefy) %>% filter(n_studies <= 8),
       trace = F)
n_taxa_studies_reduced <- gamlss(n_taxa ~ pb(n_studies), 
                              data = bind_rows(n_taxa_industrialized_studies_rarefy,
                                               n_taxa_non_industrialized_studies_rarefy) %>% filter(n_studies <= 8),
                              trace = F)
LR.test(n_taxa_studies_reduced, n_taxa_studies_full)
1 - deviance(n_taxa_studies_full) / deviance(n_taxa_studies_reduced)

# for individuals
print("individual rarefaction")
otu_table_lfst %>% select(Lifestyle, subject_ID) %>% table %>% `>`(0) %>% rowSums()
n_taxa_industrialized_subjects_rarefy <- map_dfr(rep(1:18*100, each = 10), function(x) {
  otu_table_lfst %>%
    filter(Lifestyle == "Industrialized") %>%
    select(-Lifestyle, -study, -sample_ID, -age) %>%
    subsample_fn(x, ., cat = "subject_ID") %>%
    data.frame(n_taxa = ., n_subjects = x, Lifestyle = "Industrialized")
})

n_taxa_non_industrialized_subjects_rarefy <- lapply(rep(1:9*100, each = 10), function(x) {
  otu_table_lfst %>%
    filter(Lifestyle == "Non-industrialized") %>%
    select(-Lifestyle, -study, -sample_ID, -age) %>%
    subsample_fn(x, ., cat = "subject_ID") %>%
    data.frame(n_taxa = ., n_subjects = x, Lifestyle = "Non-industrialized")
})

bind_rows(n_taxa_industrialized_subjects_rarefy, n_taxa_non_industrialized_subjects_rarefy) %>% 
  ggplot(.,aes(x = n_subjects, y = n_taxa, color = Lifestyle)) +
  geom_point() +
  geom_smooth()
# gamlss-test:
n_taxa_subjects_full <- gamlss(n_taxa ~ pb(n_subjects) + Lifestyle, 
                              data = bind_rows(n_taxa_industrialized_subjects_rarefy,
                                               n_taxa_non_industrialized_subjects_rarefy) %>% filter(n_subjects <= 900),
                              trace = F)
n_taxa_subjects_reduced <- gamlss(n_taxa ~ pb(n_subjects), 
                                 data = bind_rows(n_taxa_industrialized_subjects_rarefy,
                                                  n_taxa_non_industrialized_subjects_rarefy) %>% filter(n_subjects <= 900),
                                 trace = F)
LR.test(n_taxa_subjects_reduced, n_taxa_subjects_full)
1 - deviance(n_taxa_subjects_full) / deviance(n_taxa_subjects_reduced)


# for samples
print("sample rarefaction")
otu_table_lfst$Lifestyle %>% table
n_taxa_industrialized_samples_rarefy <- map_dfr(rep(1:109*100, each = 10), function(x) {
  otu_table_lfst %>%
    filter(Lifestyle == "Industrialized") %>%
    select(-Lifestyle, -study, -subject_ID, -age) %>%
    subsample_fn(x, ., cat = "sample_ID") %>%
    data.frame(n_taxa = ., n_samples = x, Lifestyle = "Industrialized")
})

n_taxa_non_industrialized_samples_rarefy <- lapply(rep(1:38*100, each = 10), function(x) {
  otu_table_lfst %>%
    filter(Lifestyle == "Non-industrialized") %>%
    select(-Lifestyle, -study, -subject_ID, -age) %>%
    subsample_fn(x, ., cat = "sample_ID") %>%
    data.frame(n_taxa = ., n_samples = x, Lifestyle = "Non-industrialized")
})

bind_rows(n_taxa_industrialized_samples_rarefy, n_taxa_non_industrialized_samples_rarefy) %>% 
  ggplot(.,aes(x = n_samples, y = n_taxa, color = Lifestyle)) +
  geom_point() +
  geom_smooth()
# gamlss-test:
n_taxa_samples_full <- gamlss(n_taxa ~ pb(n_samples) + Lifestyle, 
                               data = bind_rows(n_taxa_industrialized_samples_rarefy,
                                                n_taxa_non_industrialized_samples_rarefy) %>% filter(n_samples <= 4100),
                               trace = F)
n_taxa_samples_reduced <- gamlss(n_taxa ~ pb(n_samples), 
                                  data = bind_rows(n_taxa_industrialized_samples_rarefy,
                                                   n_taxa_non_industrialized_samples_rarefy) %>% filter(n_samples <= 4100),
                                  trace = F)
LR.test(n_taxa_samples_reduced, n_taxa_samples_full)
1 - deviance(n_taxa_samples_full) / deviance(n_taxa_samples_reduced)

## downsampling to equal age distributions #####################################
print("downsampling to equal age distributions")
# dfrm <- otu_tbl_df
# vl <- "study"
get_list <- function(vl, dfrm, study_prob = F) {
  lst <- lapply(c(industrialized = "Industrialized", non_industrialized = "Non-industrialized"),
                function(x) {dfrm %>%
                    filter(Lifestyle == x) %>%
                    pull(vl) %>% unique}) # get all ids for that lifestyle
  len <- lapply(lst, length) %>% unlist %>% min(c(.), na.rm = T)
  if(study_prob) {  #downsample the category to more equal study sizes
    probs <- dfrm %>% group_by(study) %>%
      mutate(prob = 1/length(unique(!!sym(vl)))) %>%
      ungroup() %>%
      select(prob, !!sym(vl)) %>%
      distinct() %>%
      column_to_rownames(vl)
    # print(lst)
    # lapply(lst, function(x) print(probs[x,]))
    final_list <- lapply(lst, function(x) sample(x, len, prob = probs[x,]))
  } else {
    final_list <- lapply(lst, function(x) sample(x, len))
  }
  return(final_list)
}
downsampling_age_dist <- function(otu_tbl_df) {
  stds <- get_list("study", otu_tbl_df)
  df_red <- otu_tbl_df %>%
    filter(study %in% stds$non_industrialized | study %in% stds$industrialized) # subsample studies
  sbjcts <- get_list("subject_ID", otu_tbl_df, study_prob = T)
  df_red_sb <- df_red %>%
    filter(subject_ID %in% sbjcts$non_industrialized | subject_ID %in% sbjcts$industrialized) # subsample individuals
  df_red_y <- filter(df_red_sb, age <= 365) # same distribuition for above and below one year respectively
  df_red_o <- filter(df_red_sb, age > 365)
  smpls_y <- get_list("sample_ID", df_red_y, study_prob = T)
  df_red_y_smp <- df_red_y %>%
    filter(sample_ID %in% smpls_y$non_industrialized | sample_ID %in% smpls_y$industrialized)
  if(length(table(df_red_o$Lifestyle)) > 1) {
    smpls_o <- get_list("sample_ID", df_red_o, study_prob = T)
    df_red_o_smp <- df_red_o %>%
      filter(Lifestyle == "non_industrialized" | sample_ID %in% smpls_o$industrialized) %>%
      filter(Lifestyle == "industrialized" | sample_ID %in% smpls_o$non_industrialized) # subsample samples
    df_red_smp <- rbind(df_red_y_smp, df_red_o_smp)
  } else {
    print("No old samples found for one Lifestyle")
    df_red_smp <- df_red_y_smp
  }
  df_red_smp <- df_red_smp %>%
    select(where(~ !is.numeric(.x) || sum(.x, na.rm = TRUE) > 0))
  return(df_red_smp)
}
unique_taxa_per_ls <- function(df) {
  df %>% select(-c("study", "subject_ID", "sample_ID", "age", "sample_sum")) %>%
    group_by(Lifestyle) %>%
    summarise(across(everything(), ~ sum(.))) %>%
    column_to_rownames("Lifestyle") %>%
    `==` (0) %>%
    rowSums()
}
plan(multisession, workers = 8)
n_ls_specific_taxa <- future_map_dfr(1:1000, function(x) {
  otu_table_lfst %>%
    downsampling_age_dist() %>%
    unique_taxa_per_ls()
}, .progress = T)
n_ls_specific_taxa <- n_ls_specific_taxa %>%
  mutate(`non-Industrialized` = Industrialized,
         Industrialized = `Non-industrialized`) %>%
  select(-`Non-industrialized`) %>%
  pivot_longer(c("Industrialized", "non-Industrialized"), names_to = "Lifestyle",
               values_to = "count")
n_ls_specific_taxa %>% 
  ggplot(., aes(x = Lifestyle, y = count, color = Lifestyle)) +
  geom_violin() +
  ylim(0,NA)

# non-rarefied data:
# otu_table_lfst_raw <- data.frame(ps_object_genus_raw@otu_table %>% t()) %>%
#   rownames_to_column("run_accession") %>%
#   left_join(mdat_all) %>%
#   column_to_rownames("run_accession")
# n_ls_specific_taxa_raw <- future_map_dfr(1:1000, function(x) {
#   otu_table_lfst_raw %>%
#     downsampling_age_dist() %>%
#     unique_taxa_per_ls()
# }, .progress = T)
# n_ls_specific_taxa_raw <- n_ls_specific_taxa_raw %>%
#   mutate(abc = Industrialized,
#          Industrialized = `Non-industrialized`,
#          `Non-industrialized` = abc) %>%
#   select(-abc)
# n_ls_specific_taxa_raw %>% pivot_longer(c("Industrialized", "Non-industrialized"), names_to = "Lifestyle",
#                                     values_to = "taxa") %>%
#   ggplot(., aes(x = Lifestyle, y = taxa, color = Lifestyle)) +
#   geom_violin() +
#   ylim(0,NA)
# 
# bind_rows(n_ls_specific_taxa %>% mutate(run = "rarefied"),
#           n_ls_specific_taxa_raw %>% mutate(run = "raw")) %>%
#   pivot_longer(c("Industrialized", "Non-industrialized"), names_to = "Lifestyle",
#                values_to = "taxa") %>%
#   ggplot(., aes(x = Lifestyle, y = taxa, color = Lifestyle)) +
#   geom_violin() +
#   ylim(0,NA) +
#   facet_wrap(~run)


## check unclassified counts: #################
mdat_unassigned <- data.frame(ps_object_genus_raref@sam_data)
ggplot(mdat_unassigned, aes(y = n_unassigned_asvs_sample/n_total_asvs_sample, x = lifestyle)) +
  geom_violin()

ggplot(mdat_unassigned, aes(y = n_unassigned_asvs_sample/n_total_asvs_sample, x = age, color = lifestyle)) +
  geom_point(size = 0.2, alpha = 0.5) +
  geom_smooth(method = "loess")


ggplot(mdat_unassigned, aes(y = shannon_asv, x = lifestyle)) +
  geom_boxplot()
mdat_unassigned %>% group_by(lifestyle) %>% summarize(n_taxa = mean(n_total_asvs_sample))

ggplot(mdat_unassigned, aes(y = shannon_asv, x = read_count, color = lifestyle)) +
  geom_point() +
  geom_smooth()

ggplot(mdat_unassigned, aes(y = richness_asv, x = age, color = lifestyle)) +
  geom_point(size = 0.2, alpha = 0.5) +
  labs(y = "Observed ASVs") +
  geom_smooth(method = "loess")


ggplot(mdat_unassigned, aes(y = n_unassigned_asvs_sample_raref/richness_asv, x = age, color = lifestyle)) +
  geom_point(size = 0.2, alpha = 0.5) +
  labs(y = "Observed ASVs") +
  geom_smooth(method = "loess")


mdat_unassigned %>% select(study, lifestyle, n_total_asvs_study, n_unassigned_asvs_study) %>%
  distinct() %>%
  ggplot(., aes(y = n_unassigned_asvs_study, x = lifestyle)) +
  geom_jitter()
## model unannotated asvs ######################
mdat_unassigned <- mdat_unassigned %>%
  mutate(lifestyle_industrialized = ifelse(lifestyle == "industrialized", yes = 1, no = 0),
         study = as.factor(study),
         percentage_unassigned = n_unassigned_asvs_sample_raref/richness_asv)
full_model_unas <- gamlss(percentage_unassigned ~ pb(age) + lifestyle_industrialized + random(study),
                     data = mdat_unassigned %>% select(percentage_unassigned, age, lifestyle_industrialized, study),
                     trace = F)
intercept_model_unas <- gamlss(percentage_unassigned ~ pb(age) + random(study),
                          data = mdat_unassigned %>% select(percentage_unassigned, age, lifestyle_industrialized, study),
                          trace = F)
lr_test_res_unas <- LR.test(intercept_model_unas, full_model_unas, print = F)
# explained variance: 
1 - deviance(full_model_unas) / deviance(intercept_model_unas)


## model alpha div ########################################
# non linear age term
meta_df <- meta_df %>%
  mutate(lifestyle_industrialized = ifelse(lifestyle == "industrialized", yes = 1, no = 0),
         study = as.factor(study))
full_model <- gamlss(diversity_shannon ~ pb(age) + lifestyle_industrialized + random(study),
                     data = meta_df %>% select(diversity_shannon, age, lifestyle_industrialized, study),
                     trace = F)
# is the path different between lifestyles different
# traj_model <- gamlss(diversity_shannon ~ pb(age) + lifestyle_industrialized + random(study),
#                      data = meta_df %>% select(diversity_shannon, age, lifestyle_industrialized, study),
#                      trace = F)
# LR.test(traj_model, full_model)
intercept_model <- gamlss(diversity_shannon ~ pb(age) + random(study),
                     data = meta_df %>% select(diversity_shannon, age, lifestyle_industrialized, study),
                     trace = F)
lr_test_res <- LR.test(intercept_model, full_model, print = F)
1 - deviance(full_model) / deviance(intercept_model)

# linear age term
linear_full_model <- lmer(diversity_shannon ~ age * lifestyle_industrialized + (1|study),
                          data = meta_df %>% select(diversity_shannon, age, lifestyle_industrialized, study))
# is the path different between lifestyles different
linear_traj_model <- lmer(diversity_shannon ~ age + lifestyle_industrialized + (1|study),
                          data = meta_df %>% select(diversity_shannon, age, lifestyle_industrialized, study))
lrtest(linear_traj_model, linear_full_model)
linear_intercept_model <- lmer(diversity_shannon ~ age + age:lifestyle_industrialized + (1|study),
                               data = meta_df %>% select(diversity_shannon, age, lifestyle_industrialized, study))
test <- lrtest(linear_intercept_model, linear_full_model)

lm_full <- lm(diversity_shannon ~ age + lifestyle + study, data = meta_df )
te <- anova(lm_full) %>% rstatix::eta_squared()
lm_red <- lm(diversity_shannon ~ age + study, data = meta_df )
lrtest(lm_red, lm_full)

# remove each study once:
studies <- meta_df$study %>% levels
names(studies) <- studies
get_gamlss_res <- function(dt, st) {
  print(st)
  dt = dt %>% filter(!study %in% st) %>% select(diversity_shannon, age, lifestyle_industrialized, study)

  f_model <- gamlss(diversity_shannon ~ pb(age) + lifestyle_industrialized + random(study),
                       data = dt,
                       trace = F)
  # non_lin_age_model <- gamlss(diversity_shannon ~ age + lifestyle_industrialized + random(study),
  #                   data = dt,
  #                   trace = F)
  ls_model <- gamlss(diversity_shannon ~ pb(age)  + random(study),
                            data = dt,
                            trace = F)
  # age_model <- gamlss(diversity_shannon ~ lifestyle_industrialized  + random(study),
  #                       data = dt,
  #                       trace = F)
  # lr_test_res_age <- LR.test(age_model, f_model, print = F)
  # lr_test_res_age_non_lin <- LR.test(non_lin_age_model, f_model, print = F)
  lr_test_res_ls <- tryCatch(
    LR.test(ls_model, f_model, print = F),
    error = function(e){
      print(e)
      list(p.val = NA, df = NA)
    })
  linear_fm <- lm(diversity_shannon ~ age + lifestyle_industrialized + study,
                            data = dt)
  linear_intm <- lm(diversity_shannon ~ age + study,
                                 data = dt)
  linear_res <- lrtest(linear_intm, linear_fm)
  e2 <- anova(linear_fm) %>% rstatix::eta_squared()

  return(list(intercetp_gam = f_model$mu.coefficients %>% t() %>% data.frame() %>% pull(lifestyle_industrialized),
              exp_var_ls = 1 - deviance(f_model) / deviance(ls_model),
              # exp_var_st = 1 - deviance(f_model) / deviance(study_model),
              # exp_var_age = 1 - deviance(f_model) / deviance(age_model),
              # exp_var_age_non_lin = 1 - deviance(f_model) / deviance(non_lin_age_model),
              p_val_gam_ls = lr_test_res_ls$p.val,
              # p_val_gam_age = lr_test_res_age$p.val,
              # p_val_gam_age_non_lin = lr_test_res_age_non_lin$p.val,
              df_gam_ls = lr_test_res_ls$df,
              # df_gam_age = lr_test_res_age$df,
              # df_gam_age_non_lin = lr_test_res_age_non_lin$df,
              p_val_lin = linear_res$`Pr(>Chisq)`[2],
              e2_lin = e2[2],
              exp_var_ls_lin = 1 - deviance(linear_fm) / deviance(linear_intm),
              intercept_lin = summary(linear_fm)$coefficients[3,1]))
}
rest_test <- mclapply(c(complete = "complete", studies),
                      function(x) get_gamlss_res(dt = meta_df, st = x),
                      mc.cores = 7)
rest_test_df <- do.call(rbind, lapply(rest_test, as.data.frame)) %>%
  as.data.frame()

only_raman  <- get_gamlss_res(dt = meta_df, st = "raman_2019")

### for revision ###############
# within raman
dt_ram = meta_df %>% filter(study == "raman_2019", lifestyle == "non_industrialized") %>% 
  mutate(country_b_sa = country %in% c("BANGLADESH", "SOUTH_AFRICA")) %>%
  select(diversity_shannon, age, lifestyle_industrialized, country_b_sa)

f_model_ram <- gamlss(diversity_shannon ~ pb(age) + country_b_sa,
                  data = dt_ram,
                  trace = F)
ls_model_ram <- gamlss(diversity_shannon ~ pb(age),
                   data = dt_ram,
                   trace = F)
LR.test(ls_model_ram, f_model_ram, print = F)
1 - deviance(f_model_ram) / deviance(ls_model_ram)

# within Bangladesh
dt_bangl = meta_df %>% filter(country == "BANGLADESH") %>% 
  mutate(raman = study == "raman_2019") %>%
  select(diversity_shannon, age, raman, study)

#gehrig
f_model_bang_g <- gamlss(diversity_shannon ~ pb(age) + raman,
                      data = dt_bangl %>% filter(study != "subramanian_2014"),
                      trace = F)
ls_model_bang_g <- gamlss(diversity_shannon ~ pb(age),
                       data = dt_bangl %>% filter(study != "subramanian_2014"),
                       trace = F)
LR.test(ls_model_bang_g, f_model_bang_g, print = F)
1 - deviance(f_model_bang_g) / deviance(ls_model_bang_g)

#subramanian
f_model_bang_s <- gamlss(diversity_shannon ~ pb(age) + raman,
                         data = dt_bangl %>% filter(study != "gehrig_2019"),
                         trace = F)
ls_model_bang_s <- gamlss(diversity_shannon ~ pb(age),
                          data = dt_bangl %>% filter(study != "gehrig_2019"),
                          trace = F)
LR.test(ls_model_bang_s, f_model_bang_s, print = F)
1 - deviance(f_model_bang_s) / deviance(ls_model_bang_s)


# check with rank transformation
# rest_test_rank <- mclapply(c(complete = "complete", studies),
#                       function(x) get_gamlss_res(dt = meta_df %>%
#                                                    mutate(age = rank(age),
#                                                           diversity_shannon = rank(diversity_shannon)),
#                                                  st = x),
#                       mc.cores = 7)
# rest_test_df_rank <- do.call(rbind, lapply(rest_test_rank, as.data.frame)) %>%
#   as.data.frame()



# ggplot(rest_test_df, aes(x = e2_lin, y = exp_var_ls_lin)) +
#   geom_point()

# check sliding window
# when is alpha diversity significantly different betweeen lifestyles?
sliding_window_size <- 60
alpha_div_sliding_window <- map_dfr(seq(0, max(meta_df$age) - sliding_window_size, by = 7), function(x) {
  meta_df_filt <- meta_df %>%
    filter(age < x + sliding_window_size, age >= x)# sliding window
  if (min(table(meta_df_filt$lifestyle)) < 20) {return(tibble())}
  res <- meta_df_filt %>%
    rstatix::wilcox_test(diversity_shannon ~ lifestyle,
                         detailed = T,
                         alternative = "greater",
                         ref.group = "industrialized") %>%
    mutate(sliding_window = x + sliding_window_size/2,
           n_total = n1 + n2)
  mean_shannon <- meta_df_filt %>% group_by(lifestyle) %>% 
    summarize(mean_shannon = mean(diversity_shannon)) 
  res$delta_shannon <- mean_shannon$mean_shannon[mean_shannon$lifestyle == "industrialized"] - mean_shannon$mean_shannon[mean_shannon$lifestyle != "industrialized"]
  return(res)
}) %>%   rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  mutate(set = "Complete dataset")
alpha_div_sliding_window %>% filter(sliding_window <= 200) %>% pull(delta_shannon) %>% mean()
alpha_div_sliding_window %>% filter(sliding_window > 200) %>% pull(delta_shannon) %>% mean()

alpha_div_sliding_window_no_raman <- map_dfr(seq(0, max(meta_df$age) - sliding_window_size, by = 7), function(x) {
  meta_df_filt <- meta_df %>%
    filter(age < x + sliding_window_size, age >= x,
           study != "raman_2019")# sliding window
  if (min(table(meta_df_filt$lifestyle)) < 20) {return(tibble())}
  res <- meta_df_filt %>%
    rstatix::wilcox_test(diversity_shannon ~ lifestyle,
                         detailed = T,
                         alternative = "greater",
                         ref.group = "industrialized") %>%
    mutate(sliding_window = x + sliding_window_size/2,
           n_total = n1 + n2)
  mean_shannon <- meta_df_filt %>% group_by(lifestyle) %>% 
    summarize(mean_shannon = mean(diversity_shannon)) 
  res$delta_shannon <- mean_shannon$mean_shannon[mean_shannon$lifestyle == "industrialized"] - mean_shannon$mean_shannon[mean_shannon$lifestyle != "industrialized"]
  return(res)
}) %>% rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  mutate(set = "No Raman et al. 2019")
alpha_div_sliding_window_no_raman %>% filter(sliding_window <= 200) %>% pull(delta_shannon) %>% mean()
alpha_div_sliding_window_no_raman %>% filter(sliding_window > 200) %>% pull(delta_shannon) %>% mean()

# bind_rows(alpha_div_sliding_window, alpha_div_sliding_window_no_raman) %>%
#   ggplot(aes(x = sliding_window, y = estimate, color = p.adj < 0.05,
#              size = n_total)) +
#   geom_point() +
#   geom_line(color = "grey", size = 0.4, alpha = 0.5) +
#   facet_wrap(~set)


save(a_div_s_genus, meta_df, alpha_div_sliding_window, alpha_div_sliding_window_no_raman,
     n_taxa_industrialized_studies_rarefy, n_taxa_non_industrialized_studies_rarefy,
     n_taxa_industrialized_subjects_rarefy, n_taxa_non_industrialized_subjects_rarefy,
     n_taxa_industrialized_samples_rarefy, n_taxa_non_industrialized_samples_rarefy,
     n_ls_specific_taxa,
     file = "/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/a_div_genus_present.RData")


# beta div genus #####################################################################
start_time <- Sys.time()
load("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_pcoa_bray_genus.RData")

print("Ordination step genus")
if(ordination_genus_step){
  # run for genus data only:
  load("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_pcoa_bray_genus.RData")
  dist.bray_genus <- phyloseq::distance(ps_object_genus_raref, method = "bray")
  pcoa_bray_genus <- pcoa(dist.bray_genus)
  
  # aitchison dists:
  print("get aitch dists")
  raw_count_matrix <- ps_object_genus_raw %>% otu_table %>% as.data.frame %>% as.matrix %>% t()
  count_matrix_clr <- ps_object_genus_clr %>% otu_table %>% as.data.frame %>% as.matrix %>% t() # aitchison from phyloseq

  dist.aitch.genus <- vegdist(count_matrix_clr, method = "euclidean")

  print("calc pcoas")
  pcoa_aitch_genus_1 <- prcomp(count_matrix_clr, center = T, scale. = F)
  print("calc pcoa 2")
  pcoa_aitch_genus <- pcoa(dist.aitch.genus)
  print("calc pcoa 3")
  # pcoa_aitch_genus_robust <- pcoa(dist.aitch.genus_robust)
  print("calc pcoa 4")
  # pcoa_aitch_genus_adapt <- pcoa(dist.aitch.genus_adapt)
  print("mantel tests")
  # mant_bray_aitch <- mantel(dist.bray_genus, dist.aitch.genus, parallel = 18, method = "spearman", permutations=999)
  # mant_rob_def <- mantel(dist.aitch.genus_robust, dist.aitch.genus, parallel = 18, method = "spearman", permutations=999)
  # mant_def_adapt <- mantel(dist.aitch.genus, dist.aitch.genus_adapt, parallel = 18, method = "spearman", permutations=999)
  # mant_adapt_rob <- mantel(dist.aitch.genus_adapt, dist.aitch.genus_robust, parallel = 18, method = "spearman", permutations=999)

  # pcoas
  print("adonis2")
  # fit_adonis_genus <- adonis2(dist.bray_genus ~ age + study + lifestyle, # + subject_ID,
  #                             data = meta_df, by="margin", na.action = na.omit, parallel = n_cores)
  # fit_adonis_genus_aitch <- adonis2(dist.aitch.genus ~ age + study + lifestyle, # + subject_ID,
  #                                   data = meta_df, by="margin", na.action = na.omit, parallel = n_cores)
  print("permanova 3")
  # fit_adonis_genus_aitch_null <- adonis2(dist.aitch.genus ~ age + lifestyle + study, # + subject_ID,
  #                                       data = meta_df, by=NULL, na.action = na.omit, parallel = n_cores)
  save(pcoa_bray_genus, dist.bray_genus, dist.aitch.genus, pcoa_aitch_genus,
       dist.aitch.genus_robust, dist.aitch.genus_adapt,
       pcoa_aitch_genus_1, pcoa_aitch_genus_robust,
       pcoa_aitch_genus_adapt, mant_bray_aitch,
       mant_rob_def, mant_adapt_rob, 
       fit_adonis_genus, fit_adonis_genus_aitch, fit_adonis_genus_aitch_seq,
       fit_adonis_genus_aitch_null, 
       file = "/fast/AG_Forslund/rob/mm_index/merged_data/all/all_pcoa_bray_genus.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_pcoa_bray_genus.RData")
}

print("Ordination step genus rarefied")
if(ordination_genus_rare_step){
  # aitchison dists:
  print("get aitch dists")
  count_matrix_clr_rare_h <- ps_object_genus_raw %>%
    rarefy_even_depth(., sample.size = 2000) %>%
    microbiome::transform(transform = "clr") %>%
    otu_table %>% as.data.frame %>% as.matrix %>% t() # aitchison from phyloseq
  print("calc pcoas")
  pcoa_aitch_genus_rare_h <- prcomp(count_matrix_clr_rare_h, center = T, scale. = F)
  
  count_matrix_clr_rare_l <- ps_object_genus_raw %>%
    rarefy_even_depth(., sample.size = 700) %>%
    microbiome::transform(transform = "clr") %>%
    otu_table %>% as.data.frame %>% as.matrix %>% t() # aitchison from phyloseq
  print("calc pcoas")
  pcoa_aitch_genus_rare_l <- prcomp(count_matrix_clr_rare_l, center = T, scale. = F)
  
  count_matrix_clr_rare_vh <- ps_object_genus_raw %>%
    rarefy_even_depth(., sample.size = 10000) %>%
    microbiome::transform(transform = "clr") %>%
    otu_table %>% as.data.frame %>% as.matrix %>% t() # aitchison from phyloseq
  print("calc pcoas")
  pcoa_aitch_genus_rare_vh <- prcomp(count_matrix_clr_rare_vh, center = T, scale. = F)
  

  save(count_matrix_clr_rare_h, pcoa_aitch_genus_rare_h, count_matrix_clr_rare_l, pcoa_aitch_genus_rare_l,
       count_matrix_clr_rare_vh, pcoa_aitch_genus_rare_vh,
       file = "/fast/AG_Forslund/rob/mm_index/merged_data/all/all_pcoa_bray_genus_rare.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_pcoa_bray_genus_rare.RData")
}

# set up parallel environment for adonis
print("finished pcoa start adonis runs on cores:")
print(n_cores)
library(doParallel)
print("B")
# cl <- makeCluster(n_cores, rscript_args = c("--no-init-file", "--no-site-file","--no-environ"))
#for linux   .... cl <- makePSOCKcluster(2)
cl <- makeForkCluster(n_cores)
# cl <- makePSOCKcluster(rep("localhost", n_cores))
print("C")
registerDoParallel(cl)

if(adonis_genus_step_1){
  print("permanova 1")
  fit_adonis_genus_aitch_seq <- adonis2(dist.aitch.genus ~ age + lifestyle + study, # + subject_ID,
                                        data = meta_df, by="terms", na.action = na.omit,
                                        parallel = n_cores)
  save(fit_adonis_genus_aitch_seq,
       file = "/fast/AG_Forslund/rob/mm_index/merged_data/all/all_permanova_genus_1.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_permanova_genus_1.RData")
}

if(adonis_genus_step_2){
  print("permanova 2 extended")
  fit_adonis_genus_aitch_seq_2 <- adonis2(dist.aitch.genus ~ age + lifestyle + study + sample_sum + diversity_shannon, # + subject_ID,
                                        data = meta_df, by="margin", na.action = na.omit,
                                        parallel = n_cores,
                                        permutations = 999)
  save(fit_adonis_genus_aitch_seq_2,
       file = "/fast/AG_Forslund/rob/mm_index/merged_data/all/all_permanova_genus_2.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_permanova_genus_2.RData")
}

if(adonis_genus_step_3){
  print("permanova sex")
  filter_samples_sex <- ps_object_genus_clr@sam_data %>% data.frame() %>% filter(!is.na(sex)) %>% rownames()
  count_matrix_clr_sex <- ps_object_genus_clr %>%
    prune_samples(filter_samples_sex,.) %>%
    otu_table %>% as.data.frame %>% as.matrix %>% t() # aitchison from phyloseq
  dist.aitch.genus_sex <- vegdist(count_matrix_clr_sex, method = "euclidean")
  meta_df_sex <- ps_object_genus_clr@sam_data %>% data.frame() %>% filter(!is.na(sex)) %>%
    left_join(meta_df %>% select(diversity_shannon, run_accession), by = "run_accession")

  fit_adonis_genus_aitch_seq_sex <- adonis2(dist.aitch.genus_sex ~ age + lifestyle + study + sample_sum + diversity_shannon + sex,
                                          data = meta_df_sex, by="margin", na.action = na.omit,
                                          parallel = n_cores,
                                          permutations = 999)
  save(fit_adonis_genus_aitch_seq_sex,
       file = "/fast/AG_Forslund/rob/mm_index/merged_data/all/all_permanova_genus_sex.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_permanova_genus_sex.RData")
}

if(adonis_genus_step_4){
  print("permanova 4 rarefied")
  dist.aitch.genus_rare_h <- vegdist(count_matrix_clr_rare_h, method = "euclidean")
  meta_df_rare <- meta_df %>% filter(run_accession %in% labels(dist.aitch.genus_rare_h))
  
  fit_adonis_genus_aitch_seq_rare <- adonis2(dist.aitch.genus_rare_h ~ age + lifestyle + study, # + subject_ID,
                                        data = meta_df_rare, by="terms", na.action = na.omit,
                                        parallel = n_cores)
  save(fit_adonis_genus_aitch_seq_rare,
       file = "/fast/AG_Forslund/rob/mm_index/merged_data/all/all_permanova_genus_rare_4.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_permanova_genus_rare_4.RData")
}
stopCluster(cl)


ps_object_genus_clr %>% otu_table %>% as.data.frame %>% as.matrix %>% t() %>%
  write.csv("counts_for_maya.csv")
end_time <- Sys.time()
print("Total time for beta diversity analyses:")
print(end_time - start_time)

## subsampling genus ##################
# create subsamples of 7 studies on both sides, then take the same amount of samples
# from both lifestyles (min(industrialized, non_industrialized)) and then the same amount of samples above and below 1 year of age
# to reduce sample imbalance

# subsample studies
# check between lifestyle comparisons over time, if lifestyle microbiomes become more distant over time
set.seed(130)

# print("resampling step")
# if(resampling_step) {
#   pdf("/fast/AG_Forslund/rob/mm_index/merged_data/all/subsample_pcoa_genus.pdf")
#   lapply(1:100, function(x) repeated_subsampling(df = meta_df,
#                                                 dist_obj = dist.bray_genus,
#                                                 ps_obj = ps_object_genus_raref,
#                                                 max_studies = 5,
#                                                 max_subj = 300,
#                                                 max_samples = 1000,
#                                                 stats_output = F,
#                                                 plots = T))
#   dev.off()
#   pdf("/fast/AG_Forslund/rob/mm_index/merged_data/all/subsample_pcoa_genus_aitch.pdf")
#   lapply(1:100, function(x) repeated_subsampling(df = meta_df,
#                                                  dist_obj = dist.aitch.genus,
#                                                  ps_obj = ps_object_genus_clr,
#                                                  max_studies = 5,
#                                                  max_subj = 300,
#                                                  max_samples = 1000,
#                                                  stats_output = F,
#                                                  plots = T))
#   dev.off()
#   pdf("/fast/AG_Forslund/rob/mm_index/merged_data/all/subsample_pcoa_genus_aitch_robust.pdf")
#   lapply(1:100, function(x) repeated_subsampling(df = meta_df,
#                                                  dist_obj = dist.aitch.genus_robust,
#                                                  ps_obj = ps_object_genus_clr,
#                                                  max_studies = 5,
#                                                  max_subj = 300,
#                                                  max_samples = 1000,
#                                                  stats_output = F,
#                                                  plots = T))
#   dev.off()
#   pdf("/fast/AG_Forslund/rob/mm_index/merged_data/all/subsample_pcoa_genus_aitch_adapt.pdf")
#   lapply(1:100, function(x) repeated_subsampling(df = meta_df,
#                                                  dist_obj = dist.aitch.genus_adapt,
#                                                  ps_obj = ps_object_genus_clr,
#                                                  max_studies = 5,
#                                                  max_subj = 300,
#                                                  max_samples = 1000,
#                                                  stats_output = F,
#                                                  plots = T))
#   dev.off()
#   
#   
#   resampling_results <- mclapply(1:100, function(x) repeated_subsampling(df = meta_df,
#                                                                       dist_obj = dist.bray_genus,
#                                                                       ps_obj = ps_object_genus_raref,
#                                                                       max_studies = 5,
#                                                                       max_subj = 300,
#                                                                       max_samples = 1000,
#                                                                       stats_output = T,
#                                                                       plots = F),
#                                  mc.cores = 8)
#   resampling_results_aitch <- mclapply(1:100, function(x) repeated_subsampling(df = meta_df,
#                                                                          dist_obj = dist.aitch.genus,
#                                                                          ps_obj = ps_object_genus_clr,
#                                                                          max_studies = 5,
#                                                                          max_subj = 300,
#                                                                          max_samples = 1000,
#                                                                          stats_output = T,
#                                                                          plots = F),
#                                  mc.cores = 8)
#   resampling_results_aitch_robust <- mclapply(1:100, function(x) repeated_subsampling(df = meta_df,
#                                                                                dist_obj = dist.aitch.genus_robust,
#                                                                                ps_obj = ps_object_genus_raref,
#                                                                                max_studies = 5,
#                                                                                max_subj = 300,
#                                                                                max_samples = 1000,
#                                                                                stats_output = T,
#                                                                                plots = F),
#                                        mc.cores = 8)
#   resampling_results_aitch_adapt <- mclapply(1:100, function(x) repeated_subsampling(df = meta_df,
#                                                                                       dist_obj = dist.aitch.genus_adapt,
#                                                                                       ps_obj = ps_object_genus_raref,
#                                                                                       max_studies = 5,
#                                                                                       max_subj = 300,
#                                                                                       max_samples = 1000,
#                                                                                       stats_output = T,
#                                                                                       plots = F),
#                                               mc.cores = 8)
#   
#   saveRDS(resampling_results, "/fast/AG_Forslund/rob/mm_index/merged_data/all/resampling_genus_beta_div.rds")
#   saveRDS(resampling_results_aitch, "/fast/AG_Forslund/rob/mm_index/merged_data/all/resampling_genus_beta_div_aitch.rds")
#   saveRDS(resampling_results_aitch_robust, "/fast/AG_Forslund/rob/mm_index/merged_data/all/resampling_genus_beta_div_aitch_robust.rds")
#   saveRDS(resampling_results_aitch_adapt, "/fast/AG_Forslund/rob/mm_index/merged_data/all/resampling_genus_beta_div_aitch_adapt.rds")
# } else {
#   resampling_results <- readRDS("/fast/AG_Forslund/rob/mm_index/merged_data/all/resampling_genus_beta_div.rds")
#   resampling_results_aitch <- readRDS("/fast/AG_Forslund/rob/mm_index/merged_data/all/resampling_genus_beta_div_aitch.rds")
#   resampling_results_aitch_robust <- readRDS("/fast/AG_Forslund/rob/mm_index/merged_data/all/resampling_genus_beta_div_aitch_robust.rds")
#   resampling_results_aitch_adapt <- readRDS("/fast/AG_Forslund/rob/mm_index/merged_data/all/resampling_genus_beta_div_aitch_adapt.rds")
# }
# print("Finished resampling step")
# p_vals <- map_dfr(resampling_results, "within_ls") %>% mutate(dist = "bray")
# p_vals_aitch <- map_dfr(resampling_results_aitch, "within_ls") %>% mutate(dist = "aitch")
# p_vals_aitch_robust <- map_dfr(resampling_results_aitch_robust, "within_ls") %>% mutate(dist = "aitch_robust")
# p_vals_aitch_adapt <- map_dfr(resampling_results_aitch_adapt, "within_ls") %>% mutate(dist = "aitch_adapt")
# p_vals_combined <- bind_rows(p_vals, p_vals_aitch, p_vals_aitch_adapt, p_vals_aitch_robust)
# 
# ggplot(p_vals_combined, aes(x = estimate, y = -log(p+0.000000001), color = dist)) +
#   geom_point() +
#   # geom_density() +
#   geom_hline(yintercept = -log(0.05)) +
#   geom_vline(xintercept = 0) +
#   facet_wrap(~dist, scales = "free_x")
# ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/dists_resampling.pdf")
# 
# p_vals_combined %>%
#   pivot_longer(c("non_industrialized_mean", "industrialized_mean"),
#                names_to = "lifestyle", values_to = "mean") %>%
#   pivot_longer(c("non_industrialized_median", "industrialized_median"),
#                names_to = "lifestyle2", values_to = "median") %>%
#   ggplot(., aes(y = mean, x = lifestyle)) +
#   geom_boxplot() +
#   facet_wrap(~dist, scales = "free_y")
# p_vals_combined %>%
#   ggplot(., aes(x = non_industrialized_mean, y = industrialized_mean)) +
#   geom_point() +
#   geom_abline(slope = 1) +
#   theme_classic() +
#   facet_wrap(~dist, scales = "free")
# 
# 
# 
# # between lifestyle distances from resampling
# inter_ls_dists <- map_dfr(resampling_results, "between_ls_res")
# inter_ls_dists %>%
#   ggplot(aes(x = age_cat_1, y = mean_inter_ls)) +
#   geom_boxplot()
# 
# 
# # compare distances between/within the two groups #############
# dist_df_g <- get_dist_df(dist_obj = dist.bray_genus,
#                          mdata = meta_df[,c("run_accession", "subject_ID", "age", "age_cat", "study", "lifestyle")])
# 
# # between lifestyle distances:
# dist_df_g %>% filter(!same_ls) %>% summarise(mean_inter_ls = mean(distance), median_inter_ls = median(distance))
# 
# dist_df_g_red <- dist_df_g[c(T, rep(F, 400)),]
# inter_ls_dist <- dist_df_g_red %>% filter(!same_ls) %>% summarise(mean_inter_ls = mean(distance), median_inter_ls = median(distance))
# 
# df_p_val_dist_g <- dist_df_g_red %>%
#   filter(same_ls) %>%
#   rstatix::wilcox_test(distance ~ lifestyle_1, detailed = T) %>%
#   rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
#   rstatix::add_significance(p.col = "p.adj") %>%
#   rstatix::add_xy_position(x = "lifestyle_1", dodge = 0.8)   # important for positioning!
#     # mutate(lifestyle_1 = T,
#     # y.position = 1.1)
#   
# dist_df_g_red %>% 
#   filter(same_ls) %>%
#   group_by(lifestyle_1) %>% 
#   summarise(mean = mean(distance), median = median(distance)) %>% 
#   pivot_longer(cols = c(mean, median), names_to = "stat", values_to = "value") %>%
#   unite("lifestyle_stat", lifestyle_1, stat, sep = "_")%>%
#   pivot_wider(names_from = lifestyle_stat, values_from = value) %>%
#   cbind(df_p_val_dist_g,., inter_ls_dist)
#   
# # pivot_wider(names_from = lifestyle_1)
# dist_df_g_red %>%
#   # ggplot(aes(x = distance, color = lifestyle_1)) +
#   ggplot(aes(y = distance, x = lifestyle_1, color = lifestyle_1)) +
#   # geom_boxplot(aes(x = distance, y = lifestyle_1, color = lifestyle_1), alpha = 0.3) +
#   geom_violin() +
#   geom_boxplot(alpha = 0.3) +
#   # add_pvalue(df_p_val_dist_g,
#   #            label = "{p.adj.signif}",
#   #            # step.group.by = "variation",
#   #            step.increase = 0.05,
#   #            tip.length = 0.01,
#   #            # bracket.nudge.y = 0.02,
#   #            xmin = "xmin",
#   #            xmax = "xmax",
#   #            show.legend = FALSE) +
#   theme_classic() +
#   theme(text = element_text(size = 18),
#         legend.position= c(0.2, 0.8),
#         legend.title = element_text(size=10),
#         legend.text = element_text(size = 8))
# 
# dist_df_g %>%
#   ggplot(aes(x = distance, fill = lifestyle_1)) +
#   geom_density(alpha = 0.5) +
#   # facet_wrap(~lifestyle_1) +
#   theme_classic()

# plot PCoAs #######################

p2_g_aitch_1 <- pcoa_aitch_genus_1$x %>% as.data.frame() %>%
  select(PC1, PC2) %>%
  rownames_to_column("run_accession") %>%
  left_join(meta_df, by = "run_accession") %>%
  ggplot(.,aes(x = PC1, y = PC2, color = lifestyle)) +
  geom_point(alpha = 0.8, size = 0.3) +
  ggtitle("CLR-transformed raw counts")

# with rarefied counts:
p2_g_aitch_rare_h <- pcoa_aitch_genus_rare_h$x %>% as.data.frame() %>%
  select(PC1, PC2) %>%
  rownames_to_column("run_accession") %>%
  left_join(meta_df, by = "run_accession") %>%
  ggplot(.,aes(x = PC1, y = PC2, color = lifestyle)) +
  geom_point(alpha = 0.8, size = 0.3) +
  ggtitle("Rarefied to 2000, CLR-transformed")
p2_g_aitch_rare_l <- pcoa_aitch_genus_rare_l$x %>% as.data.frame() %>%
  select(PC1, PC2) %>%
  rownames_to_column("run_accession") %>%
  left_join(meta_df, by = "run_accession") %>%
  ggplot(.,aes(x = PC1, y = PC2, color = lifestyle)) +
  geom_point(alpha = 0.8, size = 0.3) +
  ggtitle("Rarefied to 700, CLR-transformed")
grid.arrange(p2_g_aitch_1, p2_g_aitch_rare_h, p2_g_aitch_rare_l)
p2_g_aitch_rare_vh <- pcoa_aitch_genus_rare_vh$x %>% as.data.frame() %>%
  select(PC1, PC2) %>%
  rownames_to_column("run_accession") %>%
  left_join(meta_df, by = "run_accession") %>%
  ggplot(.,aes(x = PC1, y = PC2, color = lifestyle)) +
  geom_point(alpha = 0.8, size = 0.3) +
  ggtitle("Rarefied to 10000, CLR-transformed")
grid.arrange(p2_g_aitch_1, p2_g_aitch_rare_h, p2_g_aitch_rare_l, p2_g_aitch_rare_vh)


pcoa_aitch_genus_rare_l$x[,1:2] %>% data.frame() %>% rownames_to_column() %>%
  left_join(., pcoa_aitch_genus_rare_h$x[,1:2] %>% data.frame() %>% rownames_to_column(), by = "rowname", suffix = c("_rare_l", "rare_h")) %>%
  left_join(., pcoa_aitch_genus_1$x[,1:2] %>% data.frame() %>% rownames_to_column(), by = "rowname") %>%
  left_join(., pcoa_aitch_genus_rare_vh$x[,1:2] %>% data.frame() %>% rownames_to_column(), by = "rowname", suffix = c("", "rare_vh")) %>%
  select(-rowname) %>%
  cor(, method = "spearman", use = "pairwise.complete.obs") %>%
  abs() %>%
  pheatmap::pheatmap(display_numbers = T)

# correlation with age and lifestyle:
pcoa_aitch_genus_rare_l$x[,1:2] %>% data.frame() %>% rownames_to_column() %>%
  left_join(., pcoa_aitch_genus_rare_h$x[,1:2] %>% data.frame() %>% rownames_to_column(), by = "rowname", suffix = c("_rare_l", "rare_h")) %>%
  left_join(., pcoa_aitch_genus_1$x[,1:2] %>% data.frame() %>% rownames_to_column(), by = "rowname") %>%
  left_join(., pcoa_aitch_genus_rare_vh$x[,1:2] %>% data.frame() %>% rownames_to_column(), by = "rowname", suffix = c("", "rare_vh")) %>%
  left_join(., meta_df %>% select(run_accession, lifestyle, age, sample_sum, diversity_shannon),
            by = c("rowname" = "run_accession")) %>%
  mutate(lifestyle_industrialized = ifelse(lifestyle == "industrialized", yes = 1, no = 0), .keep = "unused") %>%
  select(-rowname) %>%
  cor(method = "spearman", use = "pairwise.complete.obs") %>%
  abs() %>%
  pheatmap::pheatmap(display_numbers = T)

# drivers:
drivers_df_1 <- pcoa_aitch_genus_1$rotation %>% as.data.frame() %>%
  rownames_to_column("Taxon") %>%
  select(PC1, PC2, Taxon) %>%
  mutate(PC1 = PC1,
         PC2 = PC2)
top_drivers_1 <- bind_rows(drivers_df_1 %>% slice_max(abs(PC1), n = 10) %>% mutate(top_10 = "PC1"),
                            drivers_df_1 %>% slice_max(abs(PC2), n = 10) %>% mutate(top_10 = "PC2")) %>%
  mutate(Taxon = factor(Taxon),
         Taxon = reorder(Taxon, PC1, FUN = max, na.rm = T))
  
top_drivers_1 %>% 
  pivot_longer(cols = -c("Taxon", "top_10"), values_to = "loadings", names_to = "PC") %>%
  ggplot(., aes(y = Taxon, x = loadings, fill = PC)) +
  geom_col(position = "dodge") +
  facet_wrap(~PC) +
  theme_classic()

# for rarefied versions:
drivers_df_rare_h <- pcoa_aitch_genus_rare_h$rotation %>% as.data.frame() %>%
  rownames_to_column("Taxon") %>%
  select(PC1, PC2, Taxon) %>%
  mutate(PC1 = PC1,
         PC2 = PC2)

drivers_df_rare_l <- pcoa_aitch_genus_rare_l$rotation %>% as.data.frame() %>%
  rownames_to_column("Taxon") %>%
  select(PC1, PC2, Taxon) %>%
  mutate(PC1 = PC1,
         PC2 = PC2)

drivers_df_rare_vh <- pcoa_aitch_genus_rare_vh$rotation %>% as.data.frame() %>%
  rownames_to_column("Taxon") %>%
  select(PC1, PC2, Taxon) %>%
  mutate(PC1 = PC1,
         PC2 = PC2)


drivers_all <- full_join(drivers_df_1, drivers_df_rare_h, by = "Taxon", suffix = c("", "_rare_h")) %>%
  full_join(., drivers_df_rare_l, by = "Taxon", suffix = c("", "_rare_l")) %>%
  full_join(., drivers_df_rare_vh, by = "Taxon", suffix = c("", "_rare_vh"))
top_drivers_all <- bind_rows(drivers_all %>% slice_max(abs(PC1), n = 10),
                             drivers_all %>% slice_max(abs(PC2), n = 10),
                             drivers_all %>% slice_max(abs(PC1_rare_h), n = 10),
                             drivers_all %>% slice_max(abs(PC2_rare_h), n = 10),
                             drivers_all %>% slice_max(abs(PC1_rare_l), n = 10),
                             drivers_all %>% slice_max(abs(PC2_rare_l), n = 10),
                             drivers_all %>% slice_max(abs(PC1_rare_vh), n = 10),
                             drivers_all %>% slice_max(abs(PC2_rare_vh), n = 10)) %>%
  distinct() %>%
  mutate(Taxon = factor(Taxon),
         Taxon = reorder(Taxon, PC1, FUN = max, na.rm = T))
top_drivers_all %>% 
  pivot_longer(cols = -Taxon, values_to = "loadings", names_to = "PC") %>%
  ggplot(., aes(y = Taxon, x = loadings, fill = PC)) +
  geom_col(position = "dodge") +
  facet_wrap(~PC) +
  theme_classic()

drivers_all %>% select(-Taxon) %>% cor(, method = "spearman", use = "pairwise.complete.obs") %>% abs() %>% 
  pheatmap::pheatmap(display_numbers = T)


pcoa_metadata <- ps_object_genus_clr@sam_data %>% data.frame %>%
  left_join(., alpha.div_gen, by = c("run_accession" = "rowname"))

save(pcoa_metadata, pcoa_aitch_genus_1, top_drivers_1, pcoa_aitch_genus_rare_h, pcoa_metadata,
     file = "/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/all_pcoa_genus_combined.RData")

# tempted ######################################################################
# library(tempted)
# # only run on samples with more than 2 time points
# multiple_samples <- table(ps_object_genus_clr@sam_data$subject_ID) %>%
#   data.frame() %>%
#   filter(Freq > 2) %>%
#   left_join(., data.frame(ps_object_genus_clr@sam_data), by = c("Var1" = "subject_ID")) %>%
#   pull(run_accession)
# 
# ps_object_genus_clr_no_single <- ps_object_genus_clr %>%
#   prune_samples(samples = multiple_samples)
#   
# res_tempted <- tempted_all(t(ps_object_genus_clr_no_single@otu_table@.Data),
#                            ps_object_genus_clr_no_single@sam_data$age,
#                            ps_object_genus_clr_no_single@sam_data$subject_ID,
#                            threshold=0.95,
#                            transform="none",
#                            r=4,
#                            smooth=1e-5,
#                            pct_ratio=0.1,
#                            pct_aggregate=1)
# per_subject_mdat <- ps_object_genus_clr_no_single %>% 
#   sample_data() %>%
#   data.frame() %>%
#   select(study, subject_ID, lifestyle, birthmode) %>%
#   distinct()
# 
# # sample plot
# sample_plot_data <- res_tempted$A_hat %>%
#   data.frame %>%
#   rownames_to_column() %>%
#   left_join(.,per_subject_mdat, by = c("rowname" = "subject_ID"))
# ggplot(sample_plot_data, aes(x=PC2, y=PC3, color=lifestyle)) + 
#   geom_point() +
#   labs(title='subject loading')
# 
# sample_plot_data %>%
#   pivot_longer(., cols = starts_with("PC"), names_to = "component", values_to = "values") %>%
#   ggplot(., aes(x=component, y=values, color=lifestyle)) + 
#   geom_boxplot() +
#   labs(title='subject loading')
# 
# # temporal loadings
# plot_time_loading(res_tempted, r=4) + 
#   geom_line(size=1.5) + 
#   labs(title='temporal loadings', x='days')
# 
# # feature loadings
# ggplot(as.data.frame(res_tempted$B_hat), aes(x=PC2, y=PC3)) + 
#   geom_point() + 
#   labs(x='Component 1', y='Component 2', title='feature loading')
# 
# # log ratio of top features
# plot_metafeature(res_tempted$metafeature_ratio, per_subject_mdat[,c("subject_ID", "lifestyle")]) +
#   xlab("Days of Life")
# # interestingly: divergence in PC1-3 mainly later in life
# 
# # subject trajectories
# plot_metafeature(res_tempted$metafeature_aggregate, per_subject_mdat[,c("subject_ID", "lifestyle")]) +
#   xlab("Days of Life")
# 
# # plot samples:
# tab_feat_obs <- res_tempted$metafeature_aggregate
# colnames(tab_feat_obs)[2] <- 'subject_ID'
# colnames(tab_feat_obs)[3] <- 'age'
# tab_feat_obs <- merge(tab_feat_obs, data.frame(ps_object_genus_clr_no_single@sam_data)[,c("subject_ID" ,"age", "lifestyle")])
# reshape_feat_obs <- reshape(tab_feat_obs, 
#                             idvar=c("subject_ID","age") , 
#                             v.names=c("value"), timevar="PC",
#                             direction="wide")
# colnames(reshape_feat_obs) <-   sub(".*value[.]", "",  colnames(reshape_feat_obs))
# ggplot(data=reshape_feat_obs, aes(x=PC1, y=PC2, 
#                                                        color=lifestyle)) +
#   geom_point() +
#   # scale_color_gradient(low = "#2b83ba", high = "#d7191c") + 
#   labs(x='Component 1', y='Component 2', color='Day')

# prevalence ###################################################################
# prevalence of taxa across individuals, studies and lifestyles:
# prevalence_subject <- ps_object_family_comp%>%
#   otu_table() %>%
#   t() %>%
#   data.frame() %>%
#   cbind(meta_df[,c("subject_ID", "lifestyle", "study")]) %>%
#   select(-study, -lifestyle, -Unknown, -Mitochondria) %>%
#   group_by(subject_ID) %>%
#   summarize(across(everything(), ~ mean(.>0)))
#
# # how much detected in less or equal than x% of individuals
# prevalence_subject %>% select(-subject_ID,) %>%  `>`(0) %>% colMeans() %>% sort %>% `<=` (0.01) %>% sum#densityplot()
# # how much detected in exact n individuals
# prevalence_subject %>% select(-subject_ID) %>%  `>`(0) %>% colSums() %>% sort(., decreasing = T) %>% `==` (1) %>% sum()
#
#
# prevalence_study <- ps_object_family_comp%>%
#   otu_table() %>%
#   t() %>%
#   data.frame() %>%
#   cbind(meta_df[,c("subject_ID", "lifestyle", "study")]) %>%
#   select(-subject_ID, -lifestyle, -Unknown, -Mitochondria) %>%
#   group_by(study) %>%
#   summarize(across(everything(), ~ mean(.>0)))
#
# prevalence_study %>% select(-study,) %>%  `>`(0) %>% colMeans() %>% sort #densityplot()
# prevalence_study %>% select(-study) %>%  `>`(0) %>% colSums() %>% sort %>% `==` (19) %>% sum()
#
#
#
# prevalence_lifestyle <- ps_object_family_comp%>%
#   otu_table() %>%
#   t() %>%
#   data.frame() %>%
#   cbind(meta_df[,c("subject_ID", "lifestyle", "study")]) %>%
#   select(-subject_ID, -study, -Unknown, -Mitochondria) %>%
#   group_by(lifestyle) %>%
#   summarize(across(everything(), ~ mean(.>0))) %>%  # get prevalence
#   t() %>%
#   data.frame %>%
#   `colnames<-`(.[1, ]) %>%
#   .[-1, ] %>%
#   rownames_to_column(.,var = "Family") %>%
#   pivot_longer(., c("non_industrialized", "industrialized"), names_to = "lifestyle", values_to = "prevalence") %>%
#   mutate(Family = gsub("X.Eubacterium..coprostanoligenes.group", "E.coprostanoligenes_group", Family),
#          prevalence = as.numeric(prevalence))
#
# abundance_lifestyle <- ps_object_family_comp%>%
#   otu_table() %>%
#   t() %>%
#   data.frame() %>%
#   cbind(meta_df[,c("subject_ID", "lifestyle", "study")]) %>%
#   select(-subject_ID, -study, -Unknown, -Mitochondria) %>%
#   group_by(lifestyle) %>%
#   summarize(across(everything(), ~ mean(.[.>0]))) %>% # get mean abundance
#   t() %>%
#   data.frame %>%
#   `colnames<-`(.[1, ]) %>%
#   .[-1, ] %>%
#   rownames_to_column(.,var = "Family") %>%
#   pivot_longer(., c("non_industrialized", "industrialized"), names_to = "lifestyle", values_to = "abundance") %>%
#   mutate(Family = gsub("X.Eubacterium..coprostanoligenes.group", "E.coprostanoligenes_group", Family),
#          abundance = as.numeric(abundance))
#

# prevalence_lifestyle %>% select(-lifestyle,) %>%  `>`(0) %>% colMeans() %>% sort# %>% densityplot()
# prevalence_lifestyle %>% select(-lifestyle) %>%  `>`(0) %>% colSums() %>% sort(., decreasing = T) %>% `==` (2) %>% sum()
# prevalence_lifestyle_fam <- prevalence_lifestyle %>% select(-lifestyle) %>%  `>`(0) %>% colSums()
# prevalence_lifestyle[,c(T, prevalence_lifestyle_fam == 1)] %>% column_to_rownames(var = "lifestyle") %>% `>` (0) %>% rowSums()
# # only 2 non-industrialized specific families, but lifestyle specific families have low prevalence in general
# prevalence_lifestyle %>% column_to_rownames(var = "lifestyle") %>% t %>% data.frame %>% densityplot(data = .,  ~ industrialized)
# prevalence_lifestyle[,c(T, prevalence_lifestyle_fam == 1)] %>% column_to_rownames(var = "lifestyle") %>% t %>% data.frame %>% densityplot(data = .,  ~ industrialized)

# important families ###############################################################
# prevalence of important taxa in subject/studies/lifestyles
# important_features <- readRDS("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/family_nested_cv_important_taxa_lifestyle.rds") %>%
#   filter(measure_type == "taxon") %>%
#   group_by(taxon, lifestyle, in_lit) %>%
#   summarise(importance = mean(importance)) %>%
#   ungroup %>%
#   mutate(lifestyle_importance = lifestyle,
#          .keep = "unused")
# 
# importance_prevalence_lifestyle <- left_join(abundance_lifestyle, prevalence_lifestyle, by = c("Family", "lifestyle")) %>%
#   full_join(important_features, ., by = c("taxon" = "Family"))
# 
# 
# importance_prevalence_lifestyle %>%
#   filter(!is.na(importance)) %>%
#   ggplot(., aes(x = importance, y = prevalence, color = lifestyle)) +
#   geom_point() +
#   facet_grid(~lifestyle_importance) +
#   ylim(0,1) +
#   theme_minimal()
# 
# importance_prevalence_lifestyle %>%
#   group_by(taxon, lifestyle, abundance, prevalence) %>%
#   summarise(lifestyle_importance = case_when(length(lifestyle_importance) > 1 ~ "both",
#                                              "industrialized" %in% lifestyle_importance ~ "industrialized",
#                                              "non_industrialized" %in% lifestyle_importance ~ "non_industrialized")) %>% 
#   ggplot(., aes(x = prevalence, y = abundance, color = lifestyle_importance)) +
#   geom_point() +
#   facet_grid(~lifestyle)
# ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/figures/abundance_prevalence_important_features.pdf")
# 
# # which taxa differ between the lifestyles?
# diff_imp_taxa <- importance_prevalence_lifestyle %>%
#   group_by(taxon, lifestyle, abundance, prevalence) %>%
#   summarise(lifestyle_importance = case_when(length(lifestyle_importance) > 1 ~ "both",
#                                              "industrialized" %in% lifestyle_importance ~ "industrialized",
#                                              "non_industrialized" %in% lifestyle_importance ~ "non_industrialized")) %>% 
#   filter(lifestyle_importance != "both")
# 
# ps_object_family_comp%>%
#   otu_table() %>%
#   t() %>%
#   data.frame() %>%
#   select(unique(diff_imp_taxa$taxon)) %>%
#   cbind(meta_df[,c("subject_ID", "lifestyle", "study", "age")]) %>%
#   pivot_longer(.,cols = unique(diff_imp_taxa$taxon), names_to = "taxon", values_to = "abundance") %>%
#   ggplot(., aes(x = age, y = (abundance), fill = lifestyle, color = lifestyle, group = lifestyle)) +
#   geom_point(size = 0.3, alpha = 0.2) +
#   geom_smooth() +
#   geom_smooth(color = "black", size = 0.2) +
#   facet_wrap(~taxon, scales = "free_y") +
#   theme_minimal()
# # facet_grid(~lifestyle)
# ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/figures/lifestyle_diff_important_features.pdf")
# 
# # and which are the same?
# sim_imp_taxa <- importance_prevalence_lifestyle %>%
#   group_by(taxon, lifestyle, abundance, prevalence) %>%
#   summarise(lifestyle_importance = case_when(length(lifestyle_importance) > 1 ~ "both",
#                                              "industrialized" %in% lifestyle_importance ~ "industrialized",
#                                              "non_industrialized" %in% lifestyle_importance ~ "non_industrialized")) %>% 
#   filter(lifestyle_importance == "both")
# ps_object_family_comp%>%
#   otu_table() %>%
#   t() %>%
#   data.frame() %>%
#   mutate(E.coprostanoligenes_group = `X.Eubacterium..coprostanoligenes.group`, .keep = "unused") %>%
#   select(unique(sim_imp_taxa$taxon)) %>%
#   cbind(meta_df[,c("subject_ID", "lifestyle", "study", "age")]) %>%
#   pivot_longer(.,cols = unique(sim_imp_taxa$taxon), names_to = "taxon", values_to = "abundance") %>%
#   ggplot(., aes(x = age, y = (abundance), color = lifestyle, group = lifestyle)) +
#   geom_point(size = 0.3, alpha = 0.2) +
#   geom_smooth() +
#   geom_smooth(color = "black", size = 0.2) +
#   facet_wrap(~taxon, scales = "free_y") +
#   theme_minimal()
# # facet_grid(~lifestyle)
# ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/figures/lifestyle_common_important_features.pdf")
# 
# # features important in common model:
# important_features_all <- readRDS("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/family_nested_cv_important_taxa_no_ls.rds") %>%
#   filter(!taxon %in% c("Shannon", "Observed")) %>%
#   group_by(taxon) %>%
#   summarise(importance = mean(importance)) %>%
#   ungroup
# 
# ps_object_family_comp%>%
#   otu_table() %>%
#   t() %>%
#   data.frame() %>%
#   mutate(E.coprostanoligenes_group = `X.Eubacterium..coprostanoligenes.group`, .keep = "unused") %>%
#   select(unique(important_features_all$taxon)) %>%
#   cbind(meta_df[,c("subject_ID", "lifestyle", "study", "age")]) %>%
#   pivot_longer(.,cols = unique(important_features_all$taxon), names_to = "taxon", values_to = "abundance") %>%
#   ggplot(., aes(x = age, y = (abundance), color = lifestyle, group = lifestyle)) +
#   geom_point(size = 0.3, alpha = 0.2) +
#   geom_smooth(color = "black", size = 0.2) +
#   geom_smooth() +
#   facet_wrap(~taxon, scales = "free_y") +
#   theme_minimal()
# ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/figures/lifestyle_important_features_no_ls.pdf")



# prevalence genus #############################################################
# prevalence of taxa across individuals, studies and lifestyles:
prevalence_overall_g <- ps_object_genus_comp%>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  sapply(.,function(x) sum(x > 0.0))

prevalence_overall_g_raref <- ps_object_genus_raref %>%
  microbiome::transform(transform = "compositional") %>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  sapply(.,function(x) sum(x > 0.0))


## per subject #####################
prevalence_subject_g <- ps_object_genus_comp%>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  cbind(meta_df[,c("subject_ID", "lifestyle", "study")]) %>%
  select(-study, -lifestyle) %>%
  group_by(subject_ID) %>% 
  summarize(across(everything(), ~ mean(.>0.01)))

prevalence_subject_g_raref <- ps_object_genus_raref%>%
  microbiome::transform(transform = "compositional") %>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  cbind(meta_df[,c("subject_ID", "lifestyle", "study")]) %>%
  select(-study, -lifestyle) %>%
  group_by(subject_ID) %>% 
  summarize(across(everything(), ~ mean(.>0.01)))


# how much detected in less or equal than x% of individuals
prevalence_subject_g %>% select(-subject_ID,) %>%  `>`(0) %>% colMeans() %>% sort %>% `<=` (0.01) %>% sum#densityplot()
prevalence_subject_g_raref %>% select(-subject_ID,) %>%  `>`(0) %>% colMeans() %>% sort %>% `<=` (0.01) %>% sum#densityplot()
# how much is detected in at least x% of individuals
prevalence_subject_g %>% select(-subject_ID,) %>%  `>`(0) %>% colMeans() %>% sort %>% `>=` (0.9) %>% sum#densityplot()
prevalence_subject_g_raref %>% select(-subject_ID,) %>%  `>`(0) %>% colMeans() %>% sort %>% `>=` (0.9) %>% sum#densityplot()
# how much detected in exact n individuals
prevalence_subject_g %>% select(-subject_ID) %>%  `>`(0) %>% colSums() %>% sort(., decreasing = T) %>% `==` (1) %>% sum()
prevalence_subject_g_raref %>% select(-subject_ID) %>%  `>`(0) %>% colSums() %>% sort(., decreasing = T) %>% `==` (1) %>% sum()


prevalence_study_g <- ps_object_genus_comp%>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  cbind(meta_df[,c("subject_ID", "lifestyle", "study")]) %>%
  select(-subject_ID, -lifestyle) %>%
  group_by(study) %>% 
  summarize(across(everything(), ~ mean(.>0)))

prevalence_study_g_raref <- ps_object_genus_raref %>%
  microbiome::transform(transform = "compositional") %>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  cbind(meta_df[,c("subject_ID", "lifestyle", "study")]) %>%
  select(-subject_ID, -lifestyle) %>%
  group_by(study) %>% 
  summarize(across(everything(), ~ mean(.>0)))


# genera detected once in each study:
prevalence_study_g %>% select(-study,) %>%  `>`(0) %>% colMeans() %>% sort %>% `==` (1) %>% sum#densityplot()
prevalence_study_g_raref %>% select(-study,) %>%  `>`(0) %>% colMeans() %>% sort %>% `==` (1) %>% sum#densityplot()

prevalence_study_g %>% select(-study,) %>%  `>`(0) %>% colMeans() %>% sort #%>% densityplot()
prevalence_study_g_raref %>% select(-study,) %>%  `>`(0) %>% colMeans() %>% sort #%>% densityplot()
prevalence_study_g %>% select(-study,) %>%  `>`(0) %>% colMeans() %>% sort %>% `>=` (1) %>% sum#densityplot()
prevalence_study_g_raref %>% select(-study,) %>%  `>`(0) %>% colMeans() %>% sort %>% `>=` (1) %>% sum#densityplot()

prevalence_lifestyle_g <- ps_object_genus_comp%>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  cbind(meta_df[,c("subject_ID", "lifestyle", "study")]) %>%
  select(-subject_ID, -study) %>%
  group_by(lifestyle) %>% 
  summarize(across(everything(), ~ mean(.>0.0)))  # get prevalence

prevalence_lifestyle_g_raref <- ps_object_genus_raref%>%
  microbiome::transform(transform = "compositional") %>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  cbind(meta_df[,c("subject_ID", "lifestyle", "study")]) %>%
  select(-subject_ID, -study) %>%
  group_by(lifestyle) %>% 
  summarize(across(everything(), ~ mean(.>0.0)))  # get prevalence



save(prevalence_lifestyle_g, prevalence_lifestyle_g_raref,
     file = "/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/prevalence_lifestyle.RData")


# in both lifestyles:
prevalence_lifestyle_g %>% select(-lifestyle,) %>%  `>`(0) %>% colSums() %>% sort %>% `==` (2) %>% sum#densityplot()
prevalence_lifestyle_g %>% select(-lifestyle,) %>%  `>`(0) %>% rowSums()

# only industrialized:
prevalence_lifestyle_g %>% column_to_rownames("lifestyle") %>% t() %>% 
  data.frame() %>%
  filter(industrialized > 0, non_industrialized == 0) %>% dim()
prevalence_lifestyle_g %>% column_to_rownames("lifestyle") %>% t() %>% 
  data.frame() %>%
  filter(industrialized > 0, non_industrialized == 0) %>% arrange(industrialized) 
prevalence_lifestyle_g %>% column_to_rownames("lifestyle") %>% t() %>% 
  data.frame() %>%
  filter(industrialized > 0, non_industrialized == 0) %>% rownames() %>% sort()

# only non-industrialized:
prevalence_lifestyle_g %>% column_to_rownames("lifestyle") %>% t() %>% 
  data.frame() %>%
  filter(industrialized == 0, non_industrialized > 0) %>% dim()
prevalence_lifestyle_g %>% column_to_rownames("lifestyle") %>% t() %>% 
  data.frame() %>%
  filter(industrialized == 0, non_industrialized > 0) %>% arrange(non_industrialized)
prevalence_lifestyle_g %>% column_to_rownames("lifestyle") %>% t() %>% 
  data.frame() %>%
  filter(industrialized == 0, non_industrialized > 0) %>% rownames() %>% sort()



prevalence_lifestyle_g_long <- ps_object_genus_comp %>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  cbind(meta_df[,c("subject_ID", "lifestyle", "study")]) %>%
  select(-subject_ID, -study) %>%
  group_by(lifestyle) %>% 
  summarize(across(everything(), ~ mean(.>0))) %>%  # get prevalence
  t() %>%
  data.frame %>%
  `colnames<-`(.[1, ]) %>%
  .[-1, ] %>% 
  rownames_to_column(.,var = "Genus") %>% 
  pivot_longer(., c("non_industrialized", "industrialized"), names_to = "lifestyle", values_to = "prevalence") %>%
  mutate(Genus = gsub("X.Eubacterium..coprostanoligenes.group", "E.coprostanoligenes_group", Genus),
         prevalence = as.numeric(prevalence))

ps_object_genus_comp %>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  cbind(meta_df[,c("subject_ID", "lifestyle", "study")]) %>%
  select(lifestyle, unclassified_Bacteria) %>%
  group_by(lifestyle) %>% 
  summarize(across(everything(), ~ mean(.)))  # get prevalence

# important genera ###############################################################
important_features_g <- readRDS("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/genus_nested_cv_important_taxa_lifestyle.rds") %>%
  filter(measure_type == "taxon") %>%
  group_by(taxon, lifestyle) %>%
  summarise(importance = mean(importance)) %>%
  ungroup %>%
  mutate(lifestyle_importance = lifestyle,
         .keep = "unused")

# importance_prevalence_lifestyle_g <- left_join(abundance_lifestyle_g, prevalence_lifestyle_g, by = c("Genus", "lifestyle")) %>%
#   full_join(important_features_g, ., by = c("taxon" = "Genus"))
# 
# 
# importance_prevalence_lifestyle_g %>%
#   filter(!is.na(importance)) %>%
#   # pivot_longer(., cols = c("abundance", "prevalence"), names_to = "name", values_to = "value") %>% 
#   ggplot(., aes(x = importance, y = prevalence, color = lifestyle)) +
#   geom_point(alpha = 0.8) +
#   facet_wrap(~lifestyle_importance) +
#   ylim(0,1) +
#   theme_minimal()
# ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/figures/abundance_prevalence_importantce_genus.pdf")
# ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/figures/abundance_prevalence_importantce_genus.png")
# 
# importance_prevalence_lifestyle_g %>%
#   group_by(taxon, lifestyle, abundance, prevalence) %>%
#   summarise(lifestyle_importance = case_when(length(lifestyle_importance) > 1 ~ "both",
#                                              "industrialized" %in% lifestyle_importance ~ "industrialized",
#                                              "non_industrialized" %in% lifestyle_importance ~ "non_industrialized")) %>% 
#   ggplot(., aes(x = prevalence, y = log(abundance), color = lifestyle_importance)) +
#   geom_point() +
#   facet_grid(~lifestyle) +
#   theme_classic()
# ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/figures/abundance_prevalence_important_features_genus.pdf")
# ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/figures/abundance_prevalence_important_features_genus.png")

# which taxa differ between the lifestyles?
# diff_imp_taxa_g <- importance_prevalence_lifestyle_g %>%
#   group_by(taxon, lifestyle, abundance, prevalence) %>%
#   summarise(lifestyle_importance = case_when(length(lifestyle_importance) > 1 ~ "both",
#                                              "industrialized" %in% lifestyle_importance ~ "industrialized",
#                                              "non_industrialized" %in% lifestyle_importance ~ "non_industrialized")) %>% 
#   filter(lifestyle_importance != "both")
# 
# ps_object_genus_comp%>%
#   otu_table() %>%
#   t() %>%
#   data.frame() %>%
#   select(unique(diff_imp_taxa_g$taxon)) %>%
#   cbind(meta_df[,c("subject_ID", "lifestyle", "study", "age")]) %>%
#   pivot_longer(.,cols = unique(diff_imp_taxa_g$taxon), names_to = "taxon", values_to = "abundance") %>%
#   filter(taxon %in% c("Lactobacillus", "Monoglobus", "Prevotella", "Catenibacterium", "Rothia", "Dialister")) %>%
#   ggplot(., aes(x = age, y = (abundance), fill = lifestyle, color = lifestyle, group = lifestyle)) +
#   geom_point(size = 0.3, alpha = 0.2) +
#   geom_smooth() +
#   geom_smooth(color = "black", size = 0.2) +
#   facet_wrap(~taxon, scales = "free_y") +
#   theme_minimal()
# # facet_grid(~lifestyle)
# ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/figures/lifestyle_diff_important_features_genus.pdf")
# ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/figures/lifestyle_diff_important_features_genus.png")

# and which are the same?
# sim_imp_taxa_g <- importance_prevalence_lifestyle_g %>%
#   group_by(taxon, lifestyle, abundance, prevalence) %>%
#   summarise(lifestyle_importance = case_when(length(lifestyle_importance) > 1 ~ "both",
#                                              "industrialized" %in% lifestyle_importance ~ "industrialized",
#                                              "non_industrialized" %in% lifestyle_importance ~ "non_industrialized")) %>% 
#   filter(lifestyle_importance == "both")
# ps_object_genus_comp%>%
#   otu_table() %>%
#   t() %>%
#   data.frame() %>%
#   select(unique(sim_imp_taxa_g$taxon)) %>%
#   cbind(meta_df[,c("subject_ID", "lifestyle", "study", "age")]) %>%
#   pivot_longer(.,cols = unique(sim_imp_taxa_g$taxon), names_to = "taxon", values_to = "abundance") %>%
#   filter(taxon %in% c("Bifidobacterium", "Faecalibacterium", "Ruminococcus", "Lachnospira", "Staphylococcus")) %>%
#   ggplot(., aes(x = age, y = (abundance), color = lifestyle, group = lifestyle)) +
#   geom_point(size = 0.3, alpha = 0.2) +
#   geom_smooth() +
#   geom_smooth(color = "black", size = 0.2) +
#   facet_wrap(~taxon, scales = "free_y") +
#   theme_minimal()
# # facet_grid(~lifestyle)
# ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/figures/lifestyle_common_important_features_genus.pdf")
# ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/figures/lifestyle_common_important_features_genus.png")
# 
# # for presentation
# prensent_imp_taxa_g <- importance_prevalence_lifestyle_g %>%
#   group_by(taxon, lifestyle, abundance, prevalence) %>%
#   summarise(lifestyle_importance = case_when(length(lifestyle_importance) > 1 ~ "both",
#                                              "industrialized" %in% lifestyle_importance ~ "industrialized",
#                                              "non_industrialized" %in% lifestyle_importance ~ "non_industrialized")) %>% 
#   filter(lifestyle_importance == "both")
# ps_object_genus_comp%>%
#   otu_table() %>%
#   t() %>%
#   data.frame() %>%
#   select(unique(c("Bifidobacterium", "Prevotella", "Ruminococcus", "Lactobacillus", "Staphylococcus"))) %>%
#   cbind(meta_df[,c("subject_ID", "lifestyle", "study", "age")]) %>%
#   pivot_longer(.,cols = unique(c("Bifidobacterium", "Prevotella", "Ruminococcus", "Lactobacillus")), names_to = "taxon", values_to = "abundance") %>%
#   # filter(taxon %in% c("Bifidobacterium", "Prevotella", "Ruminococcus", "Lactobacillus", "Staphylococcus")) %>%
#   ggplot(., aes(x = age, y = abundance, color = lifestyle, group = lifestyle)) +
#   geom_point(size = 0.3, alpha = 0.2) +
#   geom_smooth() +
#   geom_smooth(color = "black", size = 0.2) +
#   facet_wrap(~taxon, scales = "free_y") +
#   theme_minimal()
#   # facet_grid(~lifestyle)
# ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/figures/lifestyle_important_features_genus_present.pdf",
#        height = 5, width = 6)




# features important in common model:
important_features_all_g <- readRDS("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/genus_nested_cv_important_taxa_no_ls.rds") %>%
  filter(!taxon %in% c("Shannon", "Observed")) %>%
  group_by(taxon) %>%
  summarise(importance = mean(importance)) %>%
  ungroup

ps_object_genus_comp%>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  select(unique(important_features_all_g$taxon)) %>%
  cbind(meta_df[,c("subject_ID", "lifestyle", "study", "age")]) %>%
  pivot_longer(.,cols = unique(important_features_all_g$taxon), names_to = "taxon", values_to = "abundance") %>%
  ggplot(., aes(x = age, y = (abundance), color = lifestyle, group = lifestyle)) +
  geom_point(size = 0.3, alpha = 0.2) +
  geom_smooth(color = "black", size = 0.2) +
  geom_smooth() +
  facet_wrap(~taxon, scales = "free_y") +
  theme_minimal()
ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/figures/lifestyle_important_features_genus_no_ls.pdf")






# unknown counts ###############################################################
# how about the Unknown counts in industrialized vs non-industrialized:
ps_object_family_comp%>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  cbind(meta_df[,c("subject_ID", "lifestyle", "study")]) %>%
  select(lifestyle, Unknown) %>%
  filter(Unknown > 0) %>%
  group_by(lifestyle) %>% 
  summarize(across(everything(), ~ mean(.)))
