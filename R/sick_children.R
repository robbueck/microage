rm(list = ls())
library(tidyverse)
library(ggplot2)
library(caret)
library(phyloseq)
library(microbiome)
library(gamlss)
library(gamlss.ggplots)
library(future)
library(furrr)
library(mgcv)
library(ggprism)
library(fastshap)
library(ranger)
library(vegan)
library(ape)
library(stats)

source("./regression_functions.R")

# Switches #####################
shap_malnurished_step <- F
pcoa_mal_step <- F
shap_preterms <- F
pcoa_preterms_step <- F
set.seed(123)
# load models ###################
load("../data/final_no_ls_genus.RData")
load("../data/all_nested_cv_genus_no_ls.RData")
load("../data/final_industrialized_genus.RData")
load("../data/all_nested_cv_genus_industrialized.RData")
load("../data/final_nonindustrialized_genus.RData")
load("../data/all_nested_cv_genus_nonindustrialized.RData")
load("../data/inter_lifestyle_pred_genus.RData")

combined_gehrig <- readRDS("../data/genus_data_no_ls_gehrig_2019.rds")
combined_sub <- readRDS("../data/genus_data_no_ls_subramanian_2014.rds")
non_industrialized_gehrig <- readRDS("../data/genus_data_nonindustrialized_gehrig_2019.rds")
non_industrialized_sub <- readRDS("../data/genus_data_nonindustrialized_subramanian_2014.rds")
# shap-values
load("../large_files/imp_shap_prev.RData")
shap_vals_combined <- bind_rows(shap_model_noindustrialized %>% mutate(model = "non-Industrialized"), 
                                shap_model_industrialized %>% mutate(model = "Industrialized"))

# Malnurishment ###########################


# load malnurished studies

ps_sub <- readRDS("../data/phyloseq_subramanian_2014_malnurished.rds")
ps_sub@sam_data$sample_sum <- sample_sums(ps_sub)
ps_gehr <- readRDS("../data/phyloseq_gehrig_2019_malnurished.rds")
ps_gehr@sam_data$sample_sum <- sample_sums(ps_gehr)
merged_ps_genus <- merge_phyloseq(ps_sub %>% aggregate_taxa(level = "genus"),
                            ps_gehr %>% aggregate_taxa(level = "genus")) %>%
  subset_samples(., health == "SAM")

  
mdat_malnurished_all <- bind_rows(data.frame(merged_ps_genus@sam_data), lifestyle_predictions_genus %>% 
                                    mutate(health = "healthy"))
# remove duplicated samples from gehrig for now
# mdat_malnurished_all <- mdat_malnurished_all[-duplicated(mdat_malnurished_all$run_accession),] 
## create data for model ###################
filtered_cdf_combined_sub <- create_caret_df(ps_sub %>% aggregate_taxa(level = "genus"),
                                    transformation = "compositional", only_multiple_samples = F,
                                    additional_cols = c("Observed", "Shannon"),
                                    filter_features = names(combined_sub$rf1$finalModel$variable.importance))
filtered_cdf_combined_gehr <- create_caret_df(ps_gehr %>% aggregate_taxa(level = "genus"),
                                             transformation = "compositional", only_multiple_samples = F,
                                             additional_cols = c("Observed", "Shannon"),
                                             filter_features = names(combined_gehrig$rf1$finalModel$variable.importance))


filtered_cdf_industrialized <- create_caret_df(merged_ps_genus,
                                     transformation = "compositional", only_multiple_samples = F,
                                     additional_cols = c("Observed", "Shannon"),
                                     filter_features = names(final_model_genus_industrialized$rf1$finalModel$variable.importance))

filtered_cdf_noindustrialized_sub <- create_caret_df(ps_sub %>% aggregate_taxa(level = "genus"),
                                    transformation = "compositional", only_multiple_samples = F,
                                    additional_cols = c("Observed", "Shannon"),
                                    filter_features = names(non_industrialized_sub$rf1$finalModel$variable.importance))

filtered_cdf_noindustrialized_gehr <- create_caret_df(ps_gehr %>% aggregate_taxa(level = "genus"),
                                           transformation = "compositional", only_multiple_samples = F,
                                           additional_cols = c("Observed", "Shannon"),
                                           filter_features = names(non_industrialized_gehrig$rf1$finalModel$variable.importance))

## apply models ################
preds_all_sub <- filtered_cdf_combined_sub$metadata %>%
  mutate(predicted_age = predict(combined_sub$rf1, newdata = filtered_cdf_combined_sub$features))
preds_all_gehr <- filtered_cdf_combined_gehr$metadata %>%
  mutate(predicted_age = predict(combined_gehrig$rf1, newdata = filtered_cdf_combined_gehr$features)) 

preds_all <- bind_rows(preds_all_sub %>% mutate(weight = as.numeric(weight)), preds_all_gehr) %>%
  select(c("age", "study", "subject_ID", "sample_ID", "health", "predicted_age", "WHZ", "WAZ", "HAZ")) %>%
  mutate(lifestyle = "non_industrialized") %>%
  bind_rows(.,genus_no_ls_nested_cv_preds %>% select(any_of(c("age", "study", "rf1", "subject_ID", "sample_ID", "lifestyle",
                                                       "WHZ", "WAZ", "HAZ"))) %>%
              mutate(predicted_age = rf1,
                     health = "healthy")) %>%
  mutate(training_set = "combined")

preds_industrialized <- filtered_cdf_industrialized$metadata %>%
  mutate(predicted_age = predict(final_model_genus_industrialized$rf1, newdata = filtered_cdf_industrialized$features),
         lifestyle = "non_industrialized") %>%
  select(c("age", "study", "subject_ID", "sample_ID", "health", "predicted_age", "lifestyle", "WHZ", "WAZ", "HAZ")) %>%
  bind_rows(.,genus_industrialized_nested_cv_preds %>% select(any_of(c("age", "study", "rf1", "subject_ID", "sample_ID", "lifestyle",
                                                         "WHZ", "WAZ", "HAZ"))) %>%
              mutate(predicted_age = rf1,
                     health = "healthy")) %>%
  mutate(training_set = "industrialized")


preds_noindustrialized_sub <- filtered_cdf_noindustrialized_sub$metadata %>%
  mutate(predicted_age = predict(non_industrialized_sub$rf1, newdata = filtered_cdf_noindustrialized_sub$features))
preds_noindustrialized_gehr <- filtered_cdf_noindustrialized_gehr$metadata %>%
  mutate(predicted_age = predict(non_industrialized_gehrig$rf1, newdata = filtered_cdf_noindustrialized_gehr$features)) 

preds_noindustrialized <- bind_rows(preds_noindustrialized_sub %>% mutate(weight = as.numeric(weight)), preds_noindustrialized_gehr) %>%
  select(c("age", "study", "subject_ID", "sample_ID", "health", "predicted_age", "WHZ", "WAZ", "HAZ")) %>%
  mutate(lifestyle = "non_industrialized") %>%
  bind_rows(.,genus_nonindustrialized_nested_cv_preds %>% select(any_of(c("age", "study", "rf1", "subject_ID", "sample_ID", "lifestyle",
                                                            "WHZ", "WAZ", "HAZ"))) %>%
              mutate(predicted_age = rf1,
                     health = "healthy")) %>%
  mutate(training_set = "non_industrialized")


all_preds_combined <- bind_rows(preds_all, preds_industrialized, preds_noindustrialized) %>%
  bind_rows(lifestyle_predictions_genus %>%   # cross-lifestyle predictions
              mutate(predicted_age = pred,
                     health = "healthy",
                     training_set = ifelse(lifestyle == "industrialized", 
                                           yes = "non_industrialized",
                                           no = "industrialized"))) %>%
  mutate(lifestyle_health = paste(lifestyle, health, sep = "_"),
         pred = predicted_age,
         healthy = ifelse(health == "healthy", yes = 1, no = 0),
         lifestyle_industrialized = ifelse(lifestyle == "industrialized", yes = 1, no = 0),
         study = factor(study)) %>%
  mutate(age_cat = cut(age, breaks = 0:112 * 7) %>% as.character()) %>%
  dplyr::group_by(age_cat, training_set) %>%
  dplyr::mutate(group_median = mean(predicted_age),
         group_sd = sd(predicted_age),
         n_healthy = sum(healthy),
         n_sick = sum(!healthy),
         ratio_h_s = log(n_healthy / n_sick, base = 2)) %>%
  ungroup() %>%
  mutate(MAZ = (predicted_age - group_median) / group_sd) %>%
  filter(age <= 600, age > 190, health != "MAM") # malnurished children are older


all_preds_combined <- all_preds_combined %>%
  filter(health == "healthy") %>%
  group_by(age_cat, training_set) %>% 
  slice_sample(n = 50) %>%
  summarize(group_median_downsampled = mean(predicted_age),
            group_sd_downsampled = sd(predicted_age)) %>%
  left_join(all_preds_combined, ., by = c("age_cat", "training_set")) %>%
  mutate(MAZ_downsampled = (predicted_age - group_median_downsampled) / group_sd_downsampled)


## plots #################
# because gam comparisons are difficult, just create the plots for GAMs and then make box-plots from MAZ
# to show the reduced age
all_preds_combined %>% 
  # filter(study %in% c("subramanian_2014", "gehrig_2019")) %>%
  # filter((!grepl("s4|s5|s6|s7|s8|s9|s10|s11|s12|s13|s14|s15|s16|s17|s18", sample_ID) & !grepl("S02|S03|S04|02G|03G|04G", sample_ID)) | health == "healthy") %>% 
  # filter(grepl("S01|S02|01G|02G", sample_ID) | health == "healthy") %>% 
  ggplot(., aes(x = age, y = predicted_age, color = lifestyle_health)) +
  geom_point(size = 0.3) +
  geom_smooth(method = "loess") + 
  xlim(0, NA) +
  facet_wrap(~training_set) +
  theme_classic()

all_preds_combined %>% 
  filter(study %in% c("subramanian_2014", "gehrig_2019")) %>%
  ggplot(., aes(x = age, y = predicted_age, color = lifestyle_health)) +
  geom_point(size = 0.3) +
  geom_smooth(method = "loess") + 
  xlim(0, NA) +
  facet_wrap(~training_set) +
  theme_classic()


all_preds_combined %>% group_by(study, health) %>%
  filter(training_set == "combined") %>%
  # filter( | health == "healthy") %>%
  # filter((!grepl("s4|s5|s6|s7|s8|s9|s10|s11|s12|s13|s14|s15|s16|s17|s18", sample_ID) & !grepl("S03|S04|03G|04G", sample_ID)) | health == "healthy") %>% 
  filter((grepl("s1|s2", sample_ID) | grepl("S01|01G", sample_ID)) | health == "healthy") %>%
  dplyr::summarize(n_samples = n(),
            n_indiv = length(unique(subject_ID)))


means_maz_mal <- all_preds_combined %>% 
  group_by(training_set, lifestyle_health) %>%
  dplyr::summarize(mean = mean(MAZ))


# only for subramanian and gehrig:
pvals_maz_mal_subset <- all_preds_combined %>%
  filter(study %in% c("subramanian_2014", "gehrig_2019")) %>%
  group_by(training_set) %>%
  rstatix::wilcox_test(MAZ ~ lifestyle_health) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>%
  rstatix::add_xy_position(x = "training_set") %>%
  left_join(.,means_maz_mal, by = c("training_set" = "training_set",
                                    "group1" = "lifestyle_health")) %>%
  left_join(.,means_maz_mal, by = c("training_set" = "training_set",
                                    "group2" = "lifestyle_health")) %>%
  mutate(MAZ_diff = abs(mean.x - mean.y)) %>%
  select(MAZ_diff, everything())

all_preds_combined %>%
  filter(study %in% c("gehrig_2019", "subramanian_2014")) %>%
  ggplot(., aes(x = training_set, y = MAZ)) +
  geom_boxplot(aes(color = lifestyle_health)) +
  add_pvalue(pvals_maz_mal_subset,
             label = "{p.adj.signif}",
             label.size = 4.5, 
             tip.length = 0.01,
             xmin = "xmin",
             xmax = "xmax",
             show.legend = FALSE) +
  ggtitle("Delayed maturation in different models") +
  theme_classic()


all_preds_combined %>% 
  ggplot(., aes(x = age, y = MAZ, color = lifestyle_health)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm") +
  xlim(0, NA) +
  ggtitle("Delayed maturation in different models") +
  facet_wrap(~training_set) +
  theme_classic()


## PCoA ########################
ps_object_genus_raw <- readRDS("../data/all_phyloseq_rf_filter_genus.rds")
ps_health_mal_genus <-  merge_phyloseq(merged_ps_genus %>% subset_samples(., health != "MAM"),
                                       ps_object_genus_raw %>% subset_samples(., age <= 600 & age > 190)) %>%
  prune_samples(samples = (sample_sums(.) != 0)) %>%
  subset_samples(., age > 0)

remove_taxa_genus <- ps_health_mal_genus %>%
  prevalence(., detection = 10)
remove_taxa_genus <- names(remove_taxa_genus[remove_taxa_genus < 5e-5])

ps_health_mal_genus <- ps_health_mal_genus %>%
  merge_taxa2_fixed(.,taxa = remove_taxa_genus, name = "Other")

ps_health_mal_genus_clr <- ps_health_mal_genus %>%
  microbiome::transform(transform = "clr")

health_mal_count_matrix_clr <- ps_health_mal_genus_clr %>% otu_table %>% as.data.frame %>% as.matrix %>% t()
ps_low_counts_health_mal <- ps_health_mal_genus %>% prune_samples(samples = (sample_sums(.) < 2000))
ps_object_health_mal_genus_raref <- ps_health_mal_genus %>% rarefy_even_depth(., sample.size = 2000) %>%
  merge_phyloseq(.,ps_low_counts_health_mal)


if(pcoa_mal_step) {
  alpha.div_gen_health_mal <- microbiome::alpha(ps_object_health_mal_genus_raref, index = c("Observed", "Shannon")) %>%
    rownames_to_column() %>%
    left_join(., mdat_malnurished_all, by = c("rowname" = "run_accession")) %>%
    mutate(lifestyle_health = paste(lifestyle, health, sep = "_"))

  pcoa_aitch_genus_mal <- prcomp(health_mal_count_matrix_clr)
  dist.aitch.genus_mal <- vegdist(health_mal_count_matrix_clr, method = "euclidean")
  fit_adonis_genus_aitch_mal <- adonis2(dist.aitch.genus_mal ~ age + lifestyle + health + study, # + subject_ID,
                                        data = mdat_malnurished_all[names(dist.aitch.genus_mal),],
                                        by="terms", na.action = na.omit,
                                            parallel = 10)
  
  save(pcoa_aitch_genus_mal, mdat_malnurished_all, fit_adonis_genus_aitch_mal, alpha.div_gen_health_mal,
       file = "../data/malnurished_pcoa.RData")
} else {
  load("../data/malnurished_pcoa.RData")
}
explained_var_mal <- ((pcoa_aitch_genus_mal$sdev^2/sum(pcoa_aitch_genus_mal$sdev^2)) * 100) %>%
  round(.,digits = 1)

mdat_malnurished_pca <- mdat_malnurished_all[rownames(pcoa_aitch_genus_mal$x),] %>%
  bind_cols(pcoa_aitch_genus_mal$x[,c(1, 2)])
mdat_malnurished_pca %>%
  mutate(health_lifestyle = paste(lifestyle, health)) %>%
  ggplot(., aes(x = PC1, y = PC2, color = health_lifestyle)) +
  geom_point(alpha = 0.8, size = 0.3) +
  scale_size_manual(values = c(0.3, 5)) +
  xlab(paste0("PC 1 (", explained_var_mal[1], "%)")) +
  ylab(paste0("PC 2 (", explained_var_mal[2], "%)")) +
  # scale_color_manual(values = fixed_colors) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
## GAMs ###########################
plan(multisession, workers = 5)  # You can change 'workers' to any number of cores you want to use

run_gamlss_tests <- function(data, var_name, filter_criteria, trace = F, outcome = "predicted_age") {
  cat("Fitting models for conditions: ", var_name, " and ", filter_criteria, "\n")
  formula_full <- as.formula(paste(outcome, "~ pb(age) +", var_name, "+ random(study)"))
  formula_intercept <- as.formula(paste(outcome, "~ pb(age) + random(study)"))
  filtered_data <- data %>%
    filter(!!rlang::parse_expr(filter_criteria)) %>%
    select(outcome, age, !!sym(var_name), study)
  # print("model formulars are:")
  # print(formula_full)
  # print(formula_traj)
  # print(formula_intercept)
  # Fit the models
  full_model <- gamlss(formula_full,
                       data = filtered_data)
  intercept_model <- gamlss(formula_intercept,
                            data = filtered_data)

  intercept_LR_test <- tryCatch(
    LR.test(intercept_model, full_model, print = F),
    error = function(e){
      print(e)
      list(p.val = NA, df = NA)
    })
  
  # list(full_m = full_model, traj_m = traj_model, int_m = intercept_model)
  data.frame(test = var_name,
             subset = filter_criteria,
             p_intercept = intercept_LR_test$p.val,
             df_gam_ls = intercept_LR_test$df,
             exp_var_intercept = 1 - deviance(full_model) / deviance(intercept_model),
             # p_traj = traj_LR_test$p.val,
             intercept = full_model$mu.coefficients[[3]],
             # shape_param = full_model$mu.coefficients[[5]],
             aic_full = full_model$aic,
             aic_intercept = intercept_model$aic)
}

# test <- run_gamlss_tests(all_preds_combined, var_name = "healthy", filter_criteria = "training_set == 'combined'")
testing_conditions <- data.frame(filter_conditions = c("training_set == 'combined' & health == 'healthy'", "training_set == 'combined'", "training_set == 'combined' & lifestyle == 'non_industrialized'",
                                                       "training_set == 'non_industrialized' & health == 'healthy'", "training_set == 'non_industrialized'", "training_set == 'non_industrialized' & lifestyle == 'non_industrialized'",
                                                       "training_set == 'industrialized' & health == 'healthy'", "training_set == 'industrialized'", "training_set == 'industrialized' & lifestyle == 'non_industrialized'"),
                       test_conditions = c("lifestyle_industrialized", "healthy", "healthy",
                                           "lifestyle_industrialized", "healthy", "healthy",
                                           "lifestyle_industrialized", "healthy", "healthy"))

results_df <- pmap_dfr(
  list(testing_conditions$filter_conditions, testing_conditions$test_conditions),
  ~ run_gamlss_tests(all_preds_combined, ..2, ..1, outcome = "predicted_age"))

testing_conditions_a_div <- data.frame(filter_conditions = c("health == 'healthy'", "lifestyle == 'non_industrialized'",
                                                             "lifestyle == 'industrialized'"),
                                 test_conditions = c("lifestyle_industrialized", "healthy", "healthy"))

alpha.div_gen_health_mal$lifestyle_industrialized <- alpha.div_gen_health_mal$lifestyle == "industrialized"
alpha.div_gen_health_mal$healthy <- alpha.div_gen_health_mal$health == "healthy"
alpha.div_gen_health_mal$study <- factor(alpha.div_gen_health_mal$study)
results_a_div <- pmap_dfr(
  list(testing_conditions_a_div$filter_conditions, testing_conditions_a_div$test_conditions),
  ~ run_gamlss_tests(alpha.div_gen_health_mal, ..2, ..1, outcome = "Shannon"))




## SHAP ###########################
# run shap analysis on mlanurished children for all models
# combined
# need to include the healthy ones also in the calculation, because they influence 
# the calculation of the shap-values, because this is based on permutation of existing abundances.

if(shap_malnurished_step) {
  shap_res_combined <- get_shap_long(final_model_genus_no_ls$rf1,
                                         test_data = bind_rows(filtered_cdf_combined_sub$features,
                                                               filtered_cdf_combined_gehr$features,
                                                               genus_train_data$features),
                                         features = intersect(colnames(filtered_cdf_combined_sub$features), 
                                                              colnames(filtered_cdf_combined_gehr$features)) %>%
                                           intersect(., colnames(genus_train_data$features)) %>%
                                           intersect(., unique(shap_vals_combined$taxon)),
                                     normalize_ab_values = F) %>%
  mutate(model = "combined")

  # industrialized
  shap_res_industrialized <- get_shap_long(final_model_genus_industrialized$rf1,
                                 test_data = bind_rows(filtered_cdf_industrialized$features,
                                                       genus_train_data$features),
                                 features = intersect(colnames(filtered_cdf_industrialized$features),
                                                      colnames(genus_train_data_industrialized$features)) %>%
                                   intersect(.,unique(shap_vals_combined$taxon)),
                                 normalize_ab_values = F) %>%
    mutate(model = "Industrialized")
  
  # non-industrialized
  shap_res_nonindustrialized <- get_shap_long(final_model_genus_nonindustrialized$rf1,
                                    test_data = bind_rows(filtered_cdf_noindustrialized_sub$features, 
                                                          filtered_cdf_noindustrialized_gehr$features,
                                                          genus_train_data$features),
                                    features = intersect(colnames(filtered_cdf_noindustrialized_sub$features), 
                                                         colnames(filtered_cdf_noindustrialized_gehr$features)) %>%
                                      intersect(., colnames(genus_train_data$features)) %>%
                                      intersect(.,unique(shap_vals_combined$taxon)),
                                    normalize_ab_values = F) %>%
    mutate(model = "non-Industrialized")
  
  save(shap_res_industrialized, shap_res_nonindustrialized,
       file = "../large_files//malnurished_shap_res.RData")
} else {
  load("../large_files/malnurished_shap_res.RData")
}
print("Finished malnurished SHAP step")
# combine shap values:
shap_mdat_mal <- bind_rows(shap_res_industrialized, shap_res_nonindustrialized) %>%
  left_join(.,mdat_malnurished_all %>% select(run_accession, age, lifestyle, health), by = "run_accession") %>%
  filter(age <= 600, age > 190, health != "MAM") %>% # malnurished children are older
  mutate(age_cat = cut(age, breaks = 0:112 * 7) %>% as.character(),
         lifestyle_health = paste(lifestyle, health, sep = "_"),
         display_taxon = gsub(".*unclassified_", "Uncl_ ",taxon)) %>%
  group_by(age_cat, model, display_taxon) %>%
  dplyr::mutate(group_median = mean(shap_value),
                group_sd = sd(shap_value)) %>%
  ungroup() %>%
  mutate(SVZ = (shap_value - group_median) / group_sd) %>%
  filter(model == "non-Industrialized", lifestyle_health != "industrialized_healthy")

shap_mdat_mal %>%
  filter(display_taxon %in% c("Bifidobacterium", "Faecalibacterium", "Lactobacillus")) %>%
  ggplot(., aes(x = age, y = shap_value, color = lifestyle_health)) +
  geom_point(size = 0.2, alpha = 0.3) +
  geom_smooth() +
  facet_wrap(~taxon) +
  theme_classic()

shap_mdat_mal %>%
  filter(taxon %in% c("Bifidobacterium", "Faecalibacterium", "Lactobacillus")) %>%
  ggplot(., aes(x = age, y = ab_value, color = lifestyle_health)) +
  geom_point(size = 0.2, alpha = 0.3) +
  geom_smooth() +
  facet_wrap(~display_taxon) +
  theme_classic()



df_counts <- shap_mdat_mal %>%
  # filter(display_taxon %in% c("Bifidobacterium", "Faecalibacterium", "Lactobacillus")) %>%
  count(display_taxon, lifestyle_health) %>%
  tidyr::complete(display_taxon, lifestyle_health, fill = list(n = 0)) %>%
  filter(n == 0) %>% pull(display_taxon) %>%
  unique()
means_svz_mal <- shap_mdat_mal %>% 
  filter(!display_taxon %in% df_counts) %>%
  group_by(display_taxon, health) %>%
  dplyr::summarize(mean = mean(shap_value))

##  heatmaps  #################
pvals_svz_mal_age <- shap_mdat_mal %>%
  mutate(age_bin = floor(ceiling(age) / 30)) %>%
  group_by(display_taxon, age_bin) %>%
  rstatix::wilcox_test(shap_value ~ health) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1))

interesting_taxa_mal <- pvals_svz_mal_age %>% group_by(display_taxon) %>%
  summarize(n_signif = mean(p.adj.signif != "ns")) %>% #pull(n_signif) %>%table
  filter(n_signif >= 0.8) %>% pull(display_taxon)


pvals_svz_mal_age %>% group_by(display_taxon) %>%
  summarize(n_sign = mean(p.adj.signif != "ns")) %>% pull(n_sign) %>% table

shap_mdat_mal %>%
  mutate(age_bin = floor(ceiling(age) / 30)) %>%
  left_join(.,pvals_svz_mal_age,
            by = c("display_taxon" = "display_taxon", "age_bin" = "age_bin")) %>%
  mutate(p.adj.signif = ifelse(p.adj.signif == "ns", yes = NA, no = p.adj.signif)) %>%
  group_by(health, display_taxon, age_bin) %>%
  summarize(mean_shap = mean(shap_value),
            significance = unique(p.adj.signif)) %>%
  ggplot(., aes(x = age_bin, y = display_taxon, fill = mean_shap)) +
  geom_tile(color = "white") +                   # White borders between tiles
  geom_text(aes(label = significance), color = "black", size = 2) +  # significance stars
  scale_fill_viridis_c(
    option = "H", 
    name = "mean SVZ") + 
  facet_wrap(~health) +
  theme_classic()

# heatmap mean shap diff:
mean_shap_diff_mal <- shap_mdat_mal %>%
  mutate(age_bin = floor(ceiling(age) / 30)) %>%
  left_join(.,pvals_svz_mal_age,
            by = c("display_taxon" = "display_taxon", "age_bin" = "age_bin")) %>%
  mutate(p.adj.signif = ifelse(p.adj.signif == "ns", yes = NA, no = p.adj.signif)) %>%
  group_by(display_taxon, age_bin) %>%
  summarize(mean_shap_healthy = mean(shap_value[health == "healthy"], na.rm = T),
            mean_shap_sam = mean(shap_value[health == "SAM"], na.rm = T),
            diff_shap = mean_shap_sam - mean_shap_healthy,
            significance = unique(p.adj.signif))# %>% pull(diff_SVZ) %>% sum()

high_effect_taxa_mal <- mean_shap_diff_mal %>% filter(abs(diff_shap) > 10) %>% pull(display_taxon) %>% unique

mean_shap_diff_mal %>%
  filter(display_taxon %in% intersect(interesting_taxa_mal, high_effect_taxa_mal)) %>%
  ggplot(., aes(x = age_bin, y = display_taxon, fill = diff_shap)) +
  geom_tile(color = "white") +                   # White borders between tiles
  geom_text(aes(label = significance), color = "black", size = 2) +  # significance stars
  scale_fill_gradient2(name = "mean SVZ diff",
                      mid = "white",
                      low = c("#30123BFF", "#4662D7FF", "#36AAF9FF", "#1AE4B6FF"),
                      high = c("#FABA39FF", "#F66B19FF", "#CB2A04FF", "#7A0403FF")) +
  xlab("Chronological age [Months]") +
  theme_classic()


# preterm children #######################################
## load data #############
ps_gibson <- readRDS("../data/phyloseq_gibson_2016.rds")
ps_ryan <- readRDS("../data/phyloseq_ryan_2019.rds")
ps_kamdar <- readRDS("../data/phyloseq_kamdar_2020.rds")

merged_ps_preterm <- merge_phyloseq(ps_gibson %>% aggregate_taxa(level = "genus"),
                                    # ps_rao %>% aggregate_taxa(level = "genus"),
                                    ps_ryan %>% aggregate_taxa(level = "genus"),
                                    # ps_gregory %>% aggregate_taxa(level = "genus")
                                    ps_kamdar %>% aggregate_taxa(level = "genus"))
mdat_preterm_all <- bind_rows(data.frame(merged_ps_preterm@sam_data), lifestyle_predictions_genus %>% 
                                    mutate(health = "healthy"))

## create data for model ###################
filtered_cdf_combined_preterm <- create_caret_df(merged_ps_preterm,
                                             transformation = "compositional", only_multiple_samples = F,
                                             additional_cols = c("Observed", "Shannon"),
                                             filter_features = names(final_model_genus_no_ls$rf1$finalModel$variable.importance))
filtered_cdf_industrialized_pretern <- create_caret_df(merged_ps_preterm,
                                             transformation = "compositional", only_multiple_samples = F,
                                             additional_cols = c("Observed", "Shannon"),
                                             filter_features = names(final_model_genus_industrialized$rf1$finalModel$variable.importance))

filtered_cdf_noindustrialized_preterm <- create_caret_df(merged_ps_preterm,
                                           transformation = "compositional", only_multiple_samples = F,
                                           additional_cols = c("Observed", "Shannon"),
                                           filter_features = names(final_model_genus_nonindustrialized$rf1$finalModel$variable.importance))

## apply models ################
preds_comb_preterm <- filtered_cdf_combined_preterm$metadata %>%
  mutate(predicted_age = predict(final_model_genus_no_ls$rf1, newdata = filtered_cdf_combined_preterm$features)) %>%
  select(c("age", "study", "subject_ID", "sample_ID", "health", "predicted_age")) %>%
  mutate(lifestyle = "industrialized") %>%
  bind_rows(.,genus_no_ls_nested_cv_preds %>% select(c("age", "study", "rf1", "subject_ID", "sample_ID", "lifestyle")) %>%
              mutate(predicted_age = rf1,
                     health = "healthy")) %>%
  mutate(training_set = "combined")

preds_industrialized_preterm <- filtered_cdf_industrialized_pretern$metadata %>%
  mutate(predicted_age = predict(final_model_genus_industrialized$rf1, newdata = filtered_cdf_industrialized_pretern$features),
         lifestyle = "industrialized") %>%
  select(c("age", "study", "subject_ID", "sample_ID", "health", "predicted_age", "lifestyle")) %>%
  bind_rows(.,genus_industrialized_nested_cv_preds %>% select(c("age", "study", "rf1", "subject_ID", "sample_ID", "lifestyle")) %>%
              mutate(predicted_age = rf1,
                     health = "healthy")) %>%
  mutate(training_set = "industrialized")

preds_noindustrialized_preterm <- filtered_cdf_noindustrialized_preterm$metadata %>%
  mutate(predicted_age = predict(final_model_genus_nonindustrialized$rf1, newdata = filtered_cdf_noindustrialized_preterm$features)) %>%
  select(c("age", "study", "subject_ID", "sample_ID", "health", "predicted_age")) %>%
  mutate(lifestyle = "industrialized") %>%
  bind_rows(.,genus_no_ls_nested_cv_preds %>% select(c("age", "study", "rf1", "subject_ID", "sample_ID", "lifestyle")) %>%
              mutate(predicted_age = rf1,
                     health = "healthy")) %>%
  mutate(training_set = "non_industrialized")


all_preds_combined_preterm <- bind_rows(preds_comb_preterm, preds_industrialized_preterm, preds_noindustrialized_preterm) %>%
  bind_rows(lifestyle_predictions_genus %>% # add cross-lifestyle predictions
              mutate(predicted_age = pred,
                     health = "healthy",
                     training_set = ifelse(lifestyle == "industrialized", 
                                           yes = "non_industrialized",
                                           no = "industrialized"))) %>%
  filter(age < 88) %>%
  mutate(lifestyle_health = paste(lifestyle, health, sep = "_"),
         pred = predicted_age,
         healthy = ifelse(health == "healthy", yes = 1, no = 0),
         lifestyle_industrialized = ifelse(lifestyle == "industrialized", yes = 1, no = 0),
         study = factor(study)) %>%
  mutate(age_cat = cut(age, breaks = 0:15 * 7)) %>%
  group_by(age_cat, training_set) %>%
  dplyr::mutate(group_median = mean(predicted_age),
         group_sd = sd(predicted_age),
         n_healthy = sum(healthy),
         n_sick = sum(!healthy),
         ratio_h_s = log(n_healthy / n_sick, base = 2)) %>%
  ungroup() %>%
  mutate(MAZ = (predicted_age - group_median) / group_sd)
         
all_preds_combined_preterm <- all_preds_combined_preterm %>%
  filter(health == "healthy") %>%
  group_by(age_cat, training_set) %>% 
  slice_sample(n = 50) %>%
  summarize(group_median_downsampled = mean(predicted_age),
            group_sd_downsampled = sd(predicted_age)) %>%
  left_join(all_preds_combined_preterm, ., by = c("age_cat", "training_set")) %>%
  mutate(MAZ_downsampled = (predicted_age - group_median_downsampled) / group_sd_downsampled)
  

## plots ############

ggplot(all_preds_combined_preterm, aes(x = age, y = predicted_age, color = lifestyle_health)) +
  geom_point(size = 0.3) +
  geom_smooth() + 
  facet_wrap(~training_set)

means_maz_preterm <- all_preds_combined_preterm %>% 
  group_by(training_set, lifestyle_health) %>%
  summarize(mean = mean(MAZ))

pvals_preterm_maz <- all_preds_combined_preterm %>%
  group_by(training_set) %>%
  rstatix::wilcox_test(MAZ ~ lifestyle_health) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>%
  rstatix::add_xy_position(x = "training_set") %>%
  left_join(.,means_maz_preterm, by = c("training_set" = "training_set",
                                    "group1" = "lifestyle_health")) %>%
  left_join(.,means_maz_preterm, by = c("training_set" = "training_set",
                                    "group2" = "lifestyle_health")) %>%
  mutate(MAZ_diff = abs(mean.x - mean.y)) %>%
  select(MAZ_diff, everything())

ggplot(all_preds_combined_preterm,
       aes(x = training_set, y = MAZ)) +
  # geom_smooth() +
  # geom_point(size = 0.3) +
  geom_boxplot(aes(color = lifestyle_health)) +
  add_pvalue(pvals_preterm_maz,
             label = "{p.adj.signif}",
             label.size = 4.5, 
             tip.length = 0.01,
             xmin = "xmin",
             xmax = "xmax",
             show.legend = FALSE) +
  theme_classic()



## PCoA ########################
ps_health_preterm_genus <-  merge_phyloseq(merged_ps_preterm,
                                       ps_object_genus_raw %>% subset_samples(., age < 88)) %>%
  prune_samples(samples = (sample_sums(.) != 0)) %>%
  subset_samples(., age > 0)

remove_taxa_genus_preterm <- ps_health_preterm_genus %>%
  prevalence(., detection = 10)
remove_taxa_genus_preterm <- names(remove_taxa_genus_preterm[remove_taxa_genus_preterm < 5e-5])

ps_health_preterm_genus <- ps_health_preterm_genus %>%
  merge_taxa2_fixed(.,taxa = remove_taxa_genus_preterm, name = "Other")

ps_health_preterm_genus_clr <- ps_health_preterm_genus %>%
  microbiome::transform(transform = "clr")

health_preterm_count_matrix_clr <- ps_health_preterm_genus_clr %>% otu_table %>% as.data.frame %>% as.matrix %>% t()

ps_low_counts_health_preterm <- ps_health_preterm_genus %>% prune_samples(samples = (sample_sums(.) < 2000))
ps_object_health_preterm_genus_raref <- ps_health_preterm_genus %>% rarefy_even_depth(., sample.size = 2000) %>%
  merge_phyloseq(.,ps_low_counts_health_preterm)

if(pcoa_preterms_step) {
  alpha.div_gen_health_preterm <- microbiome::alpha(ps_object_health_preterm_genus_raref, index = c("Observed", "Shannon")) %>%
    rownames_to_column() %>%
    left_join(., mdat_preterm_all, by = c("rowname" = "run_accession")) %>%
    mutate(lifestyle_health = paste(lifestyle, health, sep = "_"))
  
  pcoa_aitch_genus_preterm <- prcomp(health_preterm_count_matrix_clr)
  dist.aitch.genus_preterm <- vegdist(health_preterm_count_matrix_clr, method = "euclidean")
  fit_adonis_genus_aitch_preterm <- adonis2(dist.aitch.genus_preterm ~ age + lifestyle + health + study, # + subject_ID,
                                        data = mdat_preterm_all[names(dist.aitch.genus_preterm),], by="terms", na.action = na.omit,
                                        parallel = 10)
  
  save(pcoa_aitch_genus_preterm, mdat_preterm_all, fit_adonis_genus_aitch_preterm, alpha.div_gen_health_preterm,
       file = "../data/preterm_pcoa.RData")
} else {
  load("../data/preterm_pcoa.RData")
}

explained_var_pre <- ((pcoa_aitch_genus_preterm$sdev^2/sum(pcoa_aitch_genus_preterm$sdev^2)) * 100) %>%
  round(.,digits = 1)

mdat_preterm_pca <- mdat_preterm_all[rownames(pcoa_aitch_genus_preterm$x),] %>%
  bind_cols(pcoa_aitch_genus_preterm$x[,c(1, 2)]) %>%
  filter(!is.na(health))
mdat_preterm_pca %>%
  mutate(health_lifestyle = paste(lifestyle, health)) %>%
  ggplot(., aes(x = PC1, y = PC2, color = health_lifestyle)) +
  geom_point(alpha = 0.8, size = 0.3) +
  scale_size_manual(values = c(0.3, 5)) +
  xlab(paste0("PC 1 (", explained_var_pre[1], "%)")) +
  ylab(paste0("PC 2 (", explained_var_pre[2], "%)")) +
  # scale_color_manual(values = fixed_colors) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

## GAMs ##########################################
results_df_preterm <- pmap_dfr(
  list(testing_conditions$filter_conditions, testing_conditions$test_conditions),
  ~ run_gamlss_tests(all_preds_combined_preterm, ..2, ..1, outcome = "predicted_age"))

alpha.div_gen_health_preterm$lifestyle_industrialized <- alpha.div_gen_health_preterm$lifestyle == "industrialized"
alpha.div_gen_health_preterm$healthy <- alpha.div_gen_health_preterm$health == "healthy"
alpha.div_gen_health_preterm$study <- factor(alpha.div_gen_health_preterm$study)
results_a_div_preterm <- pmap_dfr(
  list(testing_conditions_a_div$filter_conditions, testing_conditions_a_div$test_conditions),
  ~ run_gamlss_tests(alpha.div_gen_health_preterm, ..2, ..1, outcome = "Shannon"))

## SHAP-analysis ####################
if(shap_preterms) {
  shap_res_combined_preterm <- get_shap_long(final_model_genus_no_ls$rf1,
                                         test_data = bind_rows(genus_train_data$features, filtered_cdf_combined_preterm$features),
                                         features = intersect(colnames(filtered_cdf_combined_preterm$features), 
                                                              colnames(genus_train_data$features)) %>%
                                           intersect(.,unique(shap_vals_combined$taxon)),
                                         normalize_ab_values = F) %>%
    mutate(model = "combined")
  # industrialized
  # need to include the healthy ones also in the calculation, because they influence 
  # the calculation of the shap-values, because this is based on permutation of existing abundances.
  shap_res_industrialized_preterm <- get_shap_long(final_model_genus_industrialized$rf1,
                                 test_data = bind_rows(genus_train_data$features, filtered_cdf_industrialized_pretern$features),
                                 # features = test_taxa) %>%
                                 features = intersect(colnames(filtered_cdf_industrialized_pretern$features),
                                                      colnames(genus_train_data$features)) %>%
                                   intersect(.,unique(shap_vals_combined$taxon)),
                                 normalize_ab_values = F) %>%
    mutate(model = "Industrialized")
  # non-industrialized
  shap_res_nonindustrialized_preterm <- get_shap_long(final_model_genus_nonindustrialized$rf1,
                                        test_data = bind_rows(genus_train_data$features, filtered_cdf_noindustrialized_preterm$features),
                                        # features = test_taxa) %>%
                                        features = intersect(colnames(filtered_cdf_noindustrialized_preterm$features),
                                                             colnames(genus_train_data$features)) %>%
                                          intersect(., unique(shap_vals_combined$taxon)),
                                        normalize_ab_values = F) %>%
    mutate(model = "non-Industrialized")
  save(shap_res_combined_preterm, shap_res_industrialized_preterm, shap_res_nonindustrialized_preterm,
       file = "../large_files/preterm_shap_res.RData")
} else {
  load("../large_files/preterm_shap_res.RData")
}
print("Finished preterm SHAP step")

# combine shap values:
shap_mdat_preterm <- shap_res_industrialized_preterm %>% #bind_rows(shap_res_industrialized_preterm, shap_res_nonindustrialized_preterm) %>%
  left_join(.,mdat_preterm_all %>% select(run_accession, age, lifestyle, health, study), by = "run_accession") %>%
  filter(age <= 87, ) %>% # preterms are younger
  mutate(age_bin = (floor(age / 14) * 14) %>% as.character(),
         lifestyle_health = paste(lifestyle, health, sep = "_"),
         display_taxon = gsub(".*unclassified_", "Uncl_ ",taxon)) %>%
  filter(model == "Industrialized", lifestyle_health != "non_industrialized_healthy") %>%
  group_by(age_bin, model, display_taxon) %>%
  dplyr::mutate(group_median = mean(shap_value),
                group_sd = sd(shap_value)) %>%
  ungroup() %>%
  mutate(SVZ = (shap_value - group_median) / group_sd)


shap_mdat_preterm %>%
  filter(taxon %in% c("Bifidobacterium", "Faecalibacterium", "Lactobacillus", "Staphylococcus")) %>%
  ggplot(., aes(x = age, y = shap_value, color = lifestyle_health)) +
  geom_point(size = 0.2, alpha = 0.3) +
  geom_smooth() +
  facet_wrap(~display_taxon, scale = "free_y") +
  theme_classic()


means_svz_preterm <- shap_mdat_preterm %>% 
  group_by(display_taxon, health) %>%
  dplyr::summarize(mean = mean(shap_value))

pvals_svz_preterm <- shap_mdat_preterm %>%
  # filter(display_taxon %in% c("Bifidobacterium", "Faecalibacterium", "Lactobacillus")) %>%
  group_by(display_taxon) %>%
  rstatix::wilcox_test(shap_value ~ health) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>%
  rstatix::add_xy_position(x = "health") %>%
  left_join(.,means_svz_preterm, by = c(#"model" = "model",
    "display_taxon" = "display_taxon",
    "group1" = "health")) %>%
  left_join(.,means_svz_preterm, by = c(#"model" = "model",
    "display_taxon" = "display_taxon",
    "group2" = "health")) %>%
  mutate(shap_diff = abs(mean.x - mean.y)) %>%
  select(shap_diff, everything())

pvals_svz_preterm_filter <- pvals_svz_preterm %>% 
  filter(#group1 != "industrialized_healthy",
    #group2 != "industrialized_healthy",
    shap_diff > 0.2,
    mean.x > mean.y)

shap_mdat_preterm %>%
  filter(display_taxon %in% pvals_svz_preterm_filter$display_taxon) %>%
  # filter(display_taxon %in% c("Observed", "Streptococcus", "Gemella")) %>%
  ggplot(., aes(x = health, y = shap_value)) +
  geom_boxplot(aes(color = health)) +
  add_pvalue(pvals_svz_preterm_filter,
             label = "{p.adj.signif}",
             label.size = 4.5, 
             tip.length = 0.01,
             xmin = "xmin",
             xmax = "xmax",
             show.legend = FALSE) +
  facet_wrap(~display_taxon) +
  theme_classic()

##  heatmaps  #################
pvals_svz_preterm_age <- shap_mdat_preterm %>%
  mutate(age_bin = floor(age / 7),
         age_bin = ifelse(age_bin == 12, yes = 11, no = age_bin)) %>% # last bin contains too few values for tests
  group_by(display_taxon, age_bin) %>%
  rstatix::wilcox_test(shap_value ~ health) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1))

interesting_taxa_pret <- pvals_svz_preterm_age %>% group_by(display_taxon) %>%
  summarize(n_signif = mean(p.adj.signif != "ns")) %>% #pull(n_signif) %>%table
  filter(n_signif >= 0.1) %>% pull(display_taxon)


shap_mdat_preterm %>%
  mutate(age_bin = floor(age / 7)) %>%
  left_join(.,pvals_svz_preterm_age,
            by = c("display_taxon" = "display_taxon", "age_bin" = "age_bin")) %>%
  filter(display_taxon %in% interesting_taxa_pret) %>%
  mutate(p.adj.signif = ifelse(p.adj.signif == "ns", yes = NA, no = p.adj.signif)) %>%
         # health = ifelse(health == "preterm", yes = paste(health, study, sep = "_"), no = health)) %>%
  group_by(health, display_taxon, age_bin) %>%
  summarize(mean_shap = mean(shap_value),
            significance = unique(p.adj.signif)) %>%
  ggplot(., aes(x = age_bin, y = display_taxon, fill = mean_shap)) +
  geom_tile(color = "white") +                   # White borders between tiles
  geom_text(aes(label = significance), color = "black", size = 2) +  # significance stars
  scale_fill_viridis_c(
    option = "H", 
    name = "mean SHAP") + 
  facet_wrap(~health) +
  theme_classic()


# check non industrialized model:
bind_rows(shap_res_industrialized_preterm, shap_res_nonindustrialized_preterm) %>%
  left_join(.,mdat_preterm_all %>% select(run_accession, age, lifestyle, health), by = "run_accession") %>%
  filter(model == "non-Industrialized") %>%
  filter(age <= 87, ) %>% # preterms are younger
  mutate(age_cat = cut(age, breaks = 0:112 * 7) %>% as.character(),
         lifestyle_health = paste(lifestyle, health, sep = "_"),
         display_taxon = gsub(".*unclassified_", "Uncl_ ",taxon)) %>%
  group_by(age_cat, model, display_taxon) %>%
  dplyr::mutate(group_median = mean(shap_value),
                group_sd = sd(shap_value)) %>%
  ungroup() %>%
  mutate(SVZ = (shap_value - group_median) / group_sd) %>%
  mutate(age_bin = floor(age / 7)) %>%
  group_by(health, display_taxon, age_bin, lifestyle, model) %>%
  summarize(mean_shap = mean(shap_value)) %>%
  ggplot(., aes(x = age_bin, y = display_taxon, fill = mean_shap)) +
  geom_tile(color = "white") +                   # White borders between tiles
  scale_fill_viridis_c(
    option = "H", 
    name = "mean SHAP") + 
  facet_grid(~health + lifestyle + model) +
  theme_classic()


# heatmap mean shap diff:
mean_shap_diff_preterm <- shap_mdat_preterm %>%
  mutate(age_bin = floor(age / 7),
         age_bin = ifelse(age_bin == 12, yes = 11, no = age_bin)) %>% # last bin contains too few values for tests
  left_join(.,pvals_svz_preterm_age,
            by = c("display_taxon" = "display_taxon", "age_bin" = "age_bin")) %>%
  mutate(p.adj.signif = ifelse(p.adj.signif == "ns", yes = NA, no = p.adj.signif),
         health_study = ifelse(health == "preterm", yes = paste(health, study, sep = "_"), no = health)) %>%
  group_by(display_taxon, age_bin) %>%
  summarize(mean_shap_healthy = mean(shap_value[health == "healthy"], na.rm = T),
            mean_shap_preterm = mean(shap_value[health == "preterm"], na.rm = T),
            diff_shap = mean_shap_preterm - mean_shap_healthy,
            significance = unique(p.adj.signif)) # %>%  # pull(diff_SVZ) %>% mean() # need to check this again, should be positive
high_effect_taxa <- mean_shap_diff_preterm %>% filter(abs(diff_shap) > 2) %>% pull(display_taxon) %>% unique
mean_shap_diff_preterm %>%
  filter(display_taxon %in% intersect(interesting_taxa_pret, high_effect_taxa)) %>%
  ggplot(., aes(x = age_bin, y = display_taxon, fill = diff_shap)) +
  geom_tile(color = "white") +                   # White borders between tiles
  geom_text(aes(label = significance), color = "black", size = 2) +  # significance stars
  xlab("Chronological age [Weeks]") +
  scale_fill_gradient2(name = "mean SVZ diff",
                       mid = "white",
                       low = c("#30123BFF", "#4662D7FF", "#36AAF9FF", "#1AE4B6FF"),
                       high = c("#FABA39FF", "#F66B19FF", "#CB2A04FF", "#7A0403FF")) +
  # facet_wrap(~study) +
  theme_classic()


shap_mdat_preterm %>%
  filter(display_taxon %in% interesting_taxa_pret) %>%
  ggplot(., aes(x = age, y = SVZ, color = lifestyle_health)) +
  geom_point(size = 0.2, alpha = 0.3) +
  geom_smooth() +
  facet_wrap(~display_taxon) +
  theme_classic()

save(shap_mdat_preterm, shap_mdat_mal, file = "../large_files/sick_shap_data.RData")