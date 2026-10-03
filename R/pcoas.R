# plot PCoAs
version
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



source("./pcoa_functions.R")
options(future.globals.maxSize = 25000 * 1024^2)

# Switches #######################
ordination_step <- F
ordination_genus_step <- F
ordination_genus_rare_step <- F
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


# alpha div genus data ##############################################################
ps_object_genus_raw <- readRDS("../data/all_phyloseq_rf_filter_genus.rds")
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


ps_object_genus_raref@sam_data$lifestyle = factor(ps_object_genus_raref@sam_data$lifestyle, levels = c("non_industrialized", "industrialized"))
alpha.div_gen <- microbiome::alpha(ps_object_genus_raref, index = c("Observed", "Shannon")) %>%
  rownames_to_column()
meta_df <- sample_data(ps_object_genus_raref) %>% 
  data.frame()

meta_df <- left_join(meta_df, alpha.div_gen, by = c("run_accession" = "rowname"))



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


# check sliding window
# when is alpha diversity significantly different between lifestyles?
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

# meta_df <- meta_df %>% mutate(birthmode = ifelse(study == "roswall_2021", yes = NA, no = birthmode)) 
# this information is not publicly available

save(meta_df, alpha_div_sliding_window, alpha_div_sliding_window_no_raman,
     n_taxa_industrialized_studies_rarefy, n_taxa_non_industrialized_studies_rarefy,
     n_taxa_industrialized_subjects_rarefy, n_taxa_non_industrialized_subjects_rarefy,
     n_taxa_industrialized_samples_rarefy, n_taxa_non_industrialized_samples_rarefy,
     n_ls_specific_taxa,
     file = "../data/a_div_genus_present.RData")


# beta div genus #####################################################################
start_time <- Sys.time()
print("Ordination step genus")
# aitchison dists:

count_matrix_clr <- ps_object_genus_clr %>% otu_table %>% as.data.frame %>% as.matrix %>% t() # aitchison from phyloseq
dist.aitch.genus <- vegdist(count_matrix_clr, method = "euclidean")

if(ordination_genus_step){
  # run for genus data only:
  print("calc pcoas")
  pcoa_aitch_genus_1 <- prcomp(count_matrix_clr, center = T, scale. = F)
  save(pcoa_aitch_genus_1,
       file = "../data/all_pcoa_bray_genus.RData")
} else {
  load("../data/all_pcoa_bray_genus.RData")
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
  

  save(pcoa_aitch_genus_rare_h,
       file = "../data/all_pcoa_bray_genus_rare.RData")
} else {
  load("../data/all_pcoa_bray_genus_rare.RData")
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
       file = "../data/all_permanova_genus_1.RData")
} else {
  load("../data/all_permanova_genus_1.RData")
}


if(adonis_genus_step_2){
  print("permanova 4 rarefied")
  dist.aitch.genus_rare_h <- vegdist(count_matrix_clr_rare_h, method = "euclidean")
  meta_df_rare <- meta_df %>% filter(run_accession %in% labels(dist.aitch.genus_rare_h))
  
  fit_adonis_genus_aitch_seq_rare <- adonis2(dist.aitch.genus_rare_h ~ age + lifestyle + study, # + subject_ID,
                                        data = meta_df_rare, by="terms", na.action = na.omit,
                                        parallel = n_cores)
  save(fit_adonis_genus_aitch_seq_rare,
       file = "../data/all_permanova_genus_rare_4.RData")
} else {
  load("../data/all_permanova_genus_rare_4.RData")
}
stopCluster(cl)

end_time <- Sys.time()
print("Total time for beta diversity analyses:")
print(end_time - start_time)

## correlation of PCAs with age and lifestyle: ########################
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

## drivers: ######################
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

save(pcoa_aitch_genus_1, top_drivers_1, pcoa_aitch_genus_rare_h,
     file = "../large_files/all_pcoa_genus_combined.RData")



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
     file = "../data/prevalence_lifestyle.RData")


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