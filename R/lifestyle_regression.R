# create model on industrialized or non-industrialized studies and predict age in the respective other set of studies
setwd("/fast/AG_Forslund/rob/mm_index/publication_R_scripts/R/")
library(caret)
library(caretEnsemble)
library(tidyverse)
library(phyloseq)
library(microbiome)
library("optparse")
library(doParallel)
library(rsample)
library(ggpmisc)
library(viridis)
library(ggplot2)
library(gbm)
library(kknn)
library(gridExtra)
library(tictoc)
library(gridExtra)
library(furrr)
library(rstatix)
library(ggprism)
library(Boruta)
library(lme4)
library(lmtest)
# library(MatchIt)
source("./alt_models.R")
source("./regression_functions.R")
source("./lifestyle_regression_functions.R")

# Switches ###############
dataset_nested_cv_genus_step <- F
dataset_nested_cv_family_step <- F
dataset_nested_cv_genus_2months_bins_step <- F
shap_step <- F
ks_step <- F


n_cores <- 5

option_list = list(
  make_option(c("-t", "--threads"), type="numeric", default=NULL))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$threads)){
  n_cores <- T
  n_cores <- 1
} else {
  n_cores <- opt$threads
}

merged_set <- "all"

cl <- makePSOCKcluster(ceiling(n_cores/2))
registerDoParallel(cl)
getDoParWorkers()
set.seed(825)
future::plan(multisession, workers = ceiling(n_cores/2))


# genus data ###################################################################

ps_object_genus_raw <- readRDS("../data/all_phyloseq_rf_filter_genus.rds") %>%
  subset_samples(age <= 730 & age > 1)
metadata <- ps_object_genus_raw %>%
  sample_data() %>%
  data.frame

print("genus data industrialized datasets")
tic()
if(dataset_nested_cv_genus_step){
  lifestyle_predictions_genus <- c("industrialized", "non_industrialized") %>% 
    future_map_dfr(~ get_predictions_lifestyle(. ,ps = ps_object_genus_raw))
  save(lifestyle_predictions_genus, file = "../data/inter_lifestyle_pred_genus.RData")
} else {
  load("../data/inter_lifestyle_pred_genus.RData")
}
toc()

# compare with prediction from same lifestyle:
load("../data/all_nested_cv_genus_industrialized.RData")
genus_industrialized_nested_cv_preds_long <- genus_industrialized_nested_cv_preds %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "industrialized")

industrialized_lm_list_genus <- get_lm_list(pred_df = genus_industrialized_nested_cv_preds_long,
                                     grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "industrialized")

load("../data/all_nested_cv_genus_nonindustrialized.RData")
genus_nonindustrialized_nested_cv_preds_long <- genus_nonindustrialized_nested_cv_preds %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "non_industrialized")
nonindustrialized_lm_list_genus <- get_lm_list(pred_df = genus_nonindustrialized_nested_cv_preds_long, grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "non_industrialized")


# default model:
load("../data/all_nested_cv_genus_no_ls.RData")
genus_nested_cv_preds_long <- genus_no_ls_nested_cv_preds %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "combined")
all_lm_list_genus <- get_lm_list(pred_df = genus_nested_cv_preds_long, grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "combined")


lifestyle_predictions_genus <- lifestyle_predictions_genus %>%
  mutate(model_name = "rf1") %>%
  mutate(training_set = ifelse(lifestyle == "industrialized", yes = "non_industrialized", no = "industrialized"))
lifestyle_lm_list_genus <- get_lm_list(pred_df = lifestyle_predictions_genus, grouping = c("study", "lifestyle", "model_name")) %>%
  mutate(training_set = ifelse(lifestyle == "industrialized", yes = "non_industrialized", no = "industrialized")) %>%
  rbind(., industrialized_lm_list_genus, nonindustrialized_lm_list_genus, all_lm_list_genus)

lifestyle_lm_list_genus <- lifestyle_lm_list_genus %>%
  mutate(training_set = factor(training_set, levels = c("non_industrialized", "combined", "industrialized")),
         lifestyle = factor(lifestyle, levels = c("industrialized", "non_industrialized")))


# stat test (stay with R2, RSME performs worse)
df_p_val_lifestyle_genus <- lifestyle_lm_list_genus %>%
  # lifestyle_lm_list_genus %>%
  filter(study != "bender_2016") %>%
  arrange(study) %>%
  mutate(training_set = factor(training_set)) %>%
  rstatix::group_by(lifestyle) %>%
  rstatix::wilcox_test(R2 ~ training_set, paired = T) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>% 
  rstatix::add_xy_position(x = "lifestyle", dodge = 0.8) 

ggplot(lifestyle_lm_list_genus, aes(x=lifestyle, y = R2)) +
  geom_boxplot(aes(fill = training_set)) +
  # ylim(0, 1) +
  xlab("Test set") +
  add_pvalue(df_p_val_lifestyle_genus,
             label = "{p.adj.signif}",
             # step.group.by = "variation",
             step.increase = 0.05,
             tip.length = 0.01,
             # bracket.nudge.y = 0.02,
             xmin = "xmin", 
             xmax = "xmax",
             show.legend = FALSE) +
  ylim(0,NA)

# family model #############################################
ps_object_family_raw <- readRDS("../data/all_phyloseq_rf_filter_family.rds") %>%
  subset_samples(., age <= 730 & age > 1) 

print("family data industrialized datasets")
tic()
if(dataset_nested_cv_family_step){
  lifestyle_predictions_family <- c("industrialized", "non_industrialized") %>% 
    future_map_dfr(~ get_predictions_lifestyle(. ,ps = ps_object_family_raw))
  save(lifestyle_predictions_family,
       file = "../data/inter_lifestyle_pred_family.RData")
} else {
  load("../data/inter_lifestyle_pred_family.RData")
}
toc()
# compare with prediction from same lifestyle:
load("../data/all_nested_cv_family_industrialized.RData")
family_industrialized_nested_cv_preds_long <- family_industrialized_nested_cv_preds %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "industrialized")

industrialized_lm_list_family <- get_lm_list(pred_df = family_industrialized_nested_cv_preds_long,
                                            grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "industrialized")

load("../data/all_nested_cv_family_nonindustrialized.RData")
family_nonindustrialized_nested_cv_preds_long <- family_nonindustrialized_nested_cv_preds %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "non_industrialized")
nonindustrialized_lm_list_family <- get_lm_list(pred_df = family_nonindustrialized_nested_cv_preds_long, 
                                               grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "non_industrialized")

# load combined model:
load("../data/all_nested_cv_family_no_ls.RData")
family_nested_cv_preds_long <- family_no_ls_nested_cv_preds %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "combined")
all_lm_list_family <- get_lm_list(pred_df = family_nested_cv_preds_long, 
                                 grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "combined")

# merge model outputs
lifestyle_predictions_family <- lifestyle_predictions_family %>%
  mutate(model_name = "rf1") %>%
  mutate(training_set = ifelse(lifestyle == "industrialized", yes = "non_industrialized", no = "industrialized"))
lifestyle_lm_list_family <- get_lm_list(pred_df = lifestyle_predictions_family, grouping = c("study", "lifestyle", "model_name")) %>%
  mutate(training_set = ifelse(lifestyle == "industrialized", yes = "non_industrialized", no = "industrialized")) %>%
  rbind(., industrialized_lm_list_family, nonindustrialized_lm_list_family, all_lm_list_family)

lifestyle_lm_list_family <- lifestyle_lm_list_family %>%
  mutate(training_set = factor(training_set, levels = c("non_industrialized", "combined", "industrialized")),
         lifestyle = factor(lifestyle, levels = c("industrialized", "non_industrialized")))

# stat test (stay with R2, RSME performs worse)
df_p_val_lifestyle_family <- lifestyle_lm_list_family %>%
  # lifestyle_lm_list_genus %>%
  filter(study != "bender_2016") %>%
  arrange(study) %>%
  mutate(training_set = factor(training_set)) %>%
  rstatix::group_by(lifestyle) %>%
  rstatix::wilcox_test(R2 ~ training_set, paired = T) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>% 
  rstatix::add_xy_position(x = "lifestyle", dodge = 0.8) 

ggplot(lifestyle_lm_list_family, aes(x=lifestyle, y = R2)) +
  geom_boxplot(aes(fill = training_set)) +
  # ylim(0, 1) +
  xlab("Test set") +
  add_pvalue(df_p_val_lifestyle_family,
             label = "{p.adj.signif}",
             # step.group.by = "variation",
             step.increase = 0.05,
             tip.length = 0.01,
             # bracket.nudge.y = 0.02,
             xmin = "xmin", 
             xmax = "xmax",
             show.legend = FALSE) +
  ylim(0,NA)

save(lifestyle_lm_list_genus, lifestyle_lm_list_family,
     file = "../data/lifestyle_perfomance_genus.RData")


## comparison predictions genus vs family level #############################
combined_family <- bind_rows(family_nested_cv_preds_long, lifestyle_predictions_family,
                             family_nonindustrialized_nested_cv_preds_long, 
                             family_industrialized_nested_cv_preds_long) %>% 
  select(lifestyle, study, training_set, run_accession, age, pred)
combined_genus <- bind_rows(genus_nested_cv_preds_long, lifestyle_predictions_genus,
                           genus_nonindustrialized_nested_cv_preds_long, 
                           genus_industrialized_nested_cv_preds_long) %>%
  select(lifestyle, study, training_set, run_accession, age, pred)
combined_genus_family <- left_join(combined_family, combined_genus,
          by = c("lifestyle", "study", "training_set", "run_accession", "age"),
          suffix = c("_family", "_genus"))
ggplot(combined_genus_family, aes(x = pred_genus, y = pred_family, color = lifestyle)) +
  geom_point(size = 0.3, alpha = 0.2) +
  geom_smooth(method = "lm") +
  facet_wrap(~training_set) +
  coord_equal()
combined_genus_family %>%
  group_by(lifestyle, training_set) %>%
  rstatix::cor_test(pred_family, pred_genus, method = "spearman")






## different feature importance #################################################
industrialized_studies <- lifestyle_lm_list_genus %>% 
  filter(lifestyle == "industrialized") %>%
  pull(study) %>% unique
non_industrialized_studies <- lifestyle_lm_list_genus %>% 
  filter(lifestyle == "non_industrialized") %>%
  pull(study) %>% unique

combined_importance_nonindustrialized_genus <- map(non_industrialized_studies, 
                                            function(x) get_feature_importance(x, lst = "nonindustrialized")) %>%
  purrr::reduce(rbind) %>%
  mutate(lifestyle = "non_industrialized") %>%
  group_by(taxon) %>%
  filter(n() > 1) %>% # only take those that occur in more than 2 models
  ungroup()
combined_importance_industrialized_genus <- map(industrialized_studies, function(x) get_feature_importance(x, lst = "industrialized")) %>%
  purrr::reduce(rbind) %>%
  mutate(lifestyle = "industrialized") %>%
  group_by(taxon) %>%
  filter(n() > 1) %>% # only take those that occur in more than 2 models
  ungroup()

bind_rows(combined_importance_nonindustrialized_genus, combined_importance_industrialized_genus) %>%
  group_by(ID, lifestyle) %>%
  dplyr::summarize(count = n()) %>%
  ggplot(., aes( x = lifestyle, y = count)) +
  geom_boxplot() +
  geom_jitter() +
  ylim(0, 100)

## important taxa stats #################
# n important samples 
important_union <- unique(c(combined_importance_nonindustrialized_genus$taxon, combined_importance_industrialized_genus$taxon))
important_industrialized <- unique(combined_importance_industrialized_genus$taxon) 
length(important_industrialized)
important_non_industrialized <- unique(combined_importance_nonindustrialized_genus$taxon) 
length(important_non_industrialized)
# shared:
length(intersect(important_industrialized, important_non_industrialized))
# specific to one:
important_industrialized[!important_industrialized %in% important_non_industrialized] %>% length
important_non_industrialized[!important_non_industrialized %in% important_industrialized] %>% length

# relative abundance of important taxa in data:
for_caret_list <- create_caret_df(ps_object = ps_object_genus_raw, transformation = "compositional",
                                  mean_ab_cutoff = 5e-5, study_prevalence_cutoff = 1,
                                  additional_cols = c("lifestyle", "study"),
                                  only_multiple_samples = F)

# combined importances
important_taxa_long <- for_caret_list$features %>% 
  select(important_union[!important_union %in% c("Observed", "Shannon")], lifestyle) %>%
  rownames_to_column %>%
  pivot_longer(., cols = -c(rowname, lifestyle),
               names_to = "taxon",
               values_to = "value")
important_taxa_long %>% group_by(rowname) %>%
  summarize(rel_abund = sum(value)) %>%
  pull(rel_abund) %>% mean
important_taxa_long %>% group_by(rowname, lifestyle) %>%
  summarize(rel_abund = sum(value)) %>%
  group_by(lifestyle) %>%
  summarize(mean_total = mean(rel_abund))
# only industial important
important_taxa_long %>% 
  filter(taxon %in% important_industrialized[!important_industrialized %in% c("Observed", "Shannon")]) %>%
  group_by(rowname, lifestyle) %>%
  summarize(rel_abund = sum(value)) %>%
  group_by(lifestyle) %>%
  summarize(mean_total = mean(rel_abund))
# only non-industrial important
important_taxa_long %>% 
  filter(taxon %in% important_non_industrialized[!important_non_industrialized %in% c("Observed", "Shannon")]) %>%
  group_by(rowname, lifestyle) %>%
  summarize(rel_abund = sum(value)) %>%
  group_by(lifestyle) %>%
  summarize(mean_total = mean(rel_abund))

# are taxa only important in one lifestyle also present in another lifestyle?
important_taxa_long %>% 
  # filter(lifestyle == "non_industrialized") %>% 
  filter(taxon %in% important_industrialized[!important_industrialized %in% important_non_industrialized]) %>%
  mutate(value = value > 5e-5) %>%
  group_by(taxon, lifestyle) %>%
  summarize(prevalence = mean(value)) %>%
  ggplot(.,aes(x = lifestyle, y = prevalence)) +
  geom_boxplot() +
  geom_jitter()
important_taxa_long %>% 
  filter(lifestyle == "industrialized") %>% 
  filter(taxon %in% important_non_industrialized[!important_non_industrialized %in% important_industrialized]) %>%
  mutate(value = value > 5e-5) %>%
  group_by(taxon) %>%
  summarize(prevalence = sum(value)) %>% pull(prevalence) %>% table


# add prevalence:
prevalence_lifestyle <- for_caret_list$features %>%
  # mutate(study_lifestyle = paste(lifestyle, study.1, sep = "xxx")) %>%
  group_by(lifestyle, study.1) %>% 
  summarize(across(everything(), ~ mean(.>0)), .groups = "drop") %>%  # get prevalence, mean over studies
  select(-study.1) %>%
  group_by(lifestyle) %>%
  summarize(across(everything(), ~ mean(.)), .groups = "drop") %>%
  t() %>%
  data.frame %>%
  `colnames<-`(.[1, ]) %>%
  .[-1, ] %>% 
  rownames_to_column(.,var = "taxon") %>% 
  pivot_longer(., -taxon, names_to = "lifestyle", values_to = "prevalence") %>%
  mutate(prevalence = as.numeric(prevalence) * 100)# %>%

combined_importance_prev_lifestyles_genus <- rbind(combined_importance_industrialized_genus, combined_importance_nonindustrialized_genus) %>%
  complete(taxon, lifestyle) %>%
  mutate(lifestyle = factor(lifestyle, levels = c("industrialized", "non_industrialized")),
         model = lifestyle,
         measure_type = ifelse(taxon %in% c("Shannon", "Observed"), yes = "a_div", no = "taxon")) %>%
  left_join(., prevalence_lifestyle, by = c("taxon", "lifestyle"))

# shorten taxon names
short_taxa_names <- combined_importance_prev_lifestyles_genus %>% select(taxon) %>%
  distinct() %>%
  mutate(display_taxon = taxon %>% gsub(".*unclassified_", "Uncl_ ",.) %>%
           make.unique() %>%
           gsub("X[.]", "", .) %>%
           gsub("[.]", " ", .) %>%
           gsub("  ", " ", .) %>%
           gsub("_", ".", .) %>%
           gsub("Family XIII AD3011 group", "Anaerovoracaceae Fam. XIII AD3011", .) %>%
           gsub("UCG 002", "Oscillospiraceae UCG 002", .))

combined_importance_prev_lifestyles_genus <- combined_importance_prev_lifestyles_genus %>% 
  left_join(.,short_taxa_names, by = "taxon") %>%
  # reorder according to importance
  mutate(display_taxon = reorder(display_taxon, importance, FUN = max),
         display_taxon = factor(display_taxon, levels = unique(c("Shannon", "Observed", levels(display_taxon)))),
         lifestyle = factor(lifestyle, levels = c("non_industrialized", "industrialized"))
  )

combined_importance_prev_lifestyles_genus %>%
  group_by(taxon, model) %>%
  summarize(prevalence = mean(prevalence, na.rm = T),
            importance = mean(importance, na.rm = T)) %>%
  ggplot(., aes(x = prevalence, y = importance, color = model)) +
  geom_point() +
  geom_smooth() +
  xlim(0,100)+
  ggpubr::stat_cor(# aes(colour = NULL),
    method = "spearman",
    cor.coef.name = "spearman",
    size = 6,
    show.legend = F)

save(combined_importance_prev_lifestyles_genus, 
     file = "../data/imp_prevs.RData")
  

## shap analysis genus ################################################################
# combined model
load("../data/final_no_ls_genus.RData")
# final_model_genus_no_ls
# shap_out_genus_no_ls
# industrialized model
load("../data/final_industrialized_genus.RData")
# final_model_genus_industrialized
# shap_out_genus_industrialized
# non industrialized model
load("../data/final_nonindustrialized_genus.RData")
# final_model_genus_nonindustrialized
# shap_out_genus_nonindustrialized
# genus_train_data_non_industrialized

# important taxa:
print("run shap analysis")
if(shap_step){
  shap_model_all <- get_shap_long(final_model_genus_no_ls$rf1,
                                  test_data = genus_train_data$features,
                                  features = unique(combined_importance_prev_lifestyles_genus$taxon))
  print("done shap all")
  shap_model_industrialized <- get_shap_long(final_model_genus_industrialized$rf1, 
                                   test_data = genus_train_data$features,
                                   features = unique(combined_importance_prev_lifestyles_genus$taxon))
  print("done shap industrialized")
  shap_model_noindustrialized <- get_shap_long(final_model_genus_nonindustrialized$rf1,
                                     test_data = genus_train_data$features,
                                     features = unique(combined_importance_prev_lifestyles_genus$taxon))
  print("done shap non industrialized")
  save(shap_model_industrialized, shap_model_noindustrialized,
       file = "../data/lifestyle_shap_data.RData")
} else {
  load("../data/lifestyle_shap_data.RData")
}

## shap vs abundances #################

merged_shap_models <- rbind(shap_model_industrialized %>% mutate(model = "industrialized"), 
                           shap_model_noindustrialized %>% mutate(model = "non_industrialized")) %>%
  left_join(., short_taxa_names, by = "taxon") %>%
  left_join(., metadata, by = "run_accession")

shap_cors <- merged_shap_models %>%
  group_by(display_taxon, model) %>%  
  # lifestyle  maybe only check for all data combined, because it's more a property of the model than of the data?
  summarise(cor_test = list(cor.test(ab_value, shap_value, method = "spearman", use = "pairwise.complete.obs")),
            display_taxon = unique(display_taxon),
            .groups = "drop") %>%
  mutate(R2 = map_dbl(cor_test, ~ .x$estimate),
         p_value = map_dbl(cor_test, ~ .x$p.value),
         p_adj = p.adjust(p_value, method = "bonferroni"),
         p_adj = ifelse(p_adj == 0,
                         yes = min(p_adj[p_adj!=0], na.rm = T), # remove 0 p-values
                         no = p_adj),
         # keep the order of the first figure
         display_taxon = factor(display_taxon, 
                                levels = levels(combined_importance_prev_lifestyles_genus$display_taxon)))

# correlation between shap and abund for each taxon in different models based on different data
shap_cors %>%
  mutate(R2 = case_when(is.na(R2) ~ 0,
                        p_adj > 0.05 ~ 0,
                        .default = R2)) %>%
  select(-cor_test, -p_value, -p_adj) %>%
  pivot_wider(names_from = "model", values_from = R2) %>%
  filter(industrialized != 0 | non_industrialized != 0) %>%
  # filter(!taxon %in% c("X.Eubacterium..coprostanoligenes.group", "UCG.002", "Observed",
  #                      "Shannon", "Lachnospiraceae.ND3007.group", "Intestinibacter",
  #                      "Erysipelatoclostridium", "Family.XIII.AD3011.group",
  #                      "Anaerostipes", "Oscillibacter", "Catenibacterium", "uncl_Firmicutes",
  #                      "X.Eubacterium..eligens.group", "Erysipelotrichaceae.UCG.003")) %>%
  ggplot(., aes(x = industrialized, y = non_industrialized)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black") +
  # Set axis limits and expand to ensure 0 is included
  # scale_x_continuous(expand = c(0, 0)) +
  # scale_y_continuous(expand = c(0, 0)) +
  theme_classic()


### mixed models #################
# linear model to test for lifestyle effect in modeling shap and abundance:
# run ks test between taxa within one model

if(ks_step) {
  model_res_shap_ab <- lapply(unique(merged_shap_models$taxon) %>% sort, 
                              function(x) ks_test_lifestyle(tx = x, 
                                                            df = merged_shap_models))
  taxa_combs <- combn(unique(merged_shap_models$taxon), m = 2, simplify = F)
  options(future.globals.maxSize = 8000 * 1024^2)
  plan(multisession, workers = 8)
  test_ks_industrialized <- map_dfr(taxa_combs, 
                                 function(x) ks_test_per_model(merged_shap_models %>% 
                                                                 filter(model == "industrialized"), x))
  test_ks_noindustrialized <- map_dfr(taxa_combs, 
                                   function(x) ks_test_per_model(merged_shap_models %>% 
                                                                   filter(model == "non_industrialized"), x))
  save(model_res_shap_ab, test_ks_industrialized, test_ks_noindustrialized, 
       file = "../data/ks_test_shap_res.RData")
} else {load("../data/ks_test_shap_res.RData")}


between_taxa_comp <- bind_rows(test_ks_industrialized %>% mutate(model = "industrialized",
                                                       comparison = "between_taxa"),
                               test_ks_noindustrialized %>%
                                 mutate(model = "non_industrialized",
                                        comparison = "between_taxa"))
percentile <- quantile(between_taxa_comp$D_ks, 0.05)


model_res_shap_ab_df <- bind_rows(model_res_shap_ab) %>%
  rstatix::adjust_pvalue(p.col = "p_ks", output.col = "q_ks", method = "fdr") %>%
  mutate(comparison = "between_models",
         empirical_pval = mean(abs(between_taxa_comp$D_ks) >= abs(D_ks)))

ks_res_all <- bind_rows(model_res_shap_ab_df,between_taxa_comp)

ggplot(ks_res_all, aes(x = D_ks, color = comparison)) +
  geom_density() +
  theme_classic()

# combine shap and feature importance:
shap_and_imps <- bind_rows(combined_importance_prev_lifestyles_genus %>% mutate(type = "importance"),
                           shap_cors %>% mutate(type = "shap") %>% filter(!is.na(R2))) %>%
  mutate(model = paste0(model, "_model"))
facet_labels <- shap_and_imps %>%
  distinct(model) %>%
  pull(model) %>%
  setNames(., .)  # Keep only the `model` as the facet label

save(combined_importance_prev_lifestyles_genus, shap_cors, short_taxa_names, 
     shap_model_industrialized, shap_model_noindustrialized,
     file = "../data/imp_shap_prev.RData")



# plots for indiviidual taxa
rbind(shap_model_industrialized, shap_model_noindustrialized) %>%
  filter(taxon %in% c("Faecalibacterium", "Bifidobacterium")) %>%
  left_join(., metadata, by = "run_accession") %>%
  group_by(taxon) %>%
  summarise(mean = mean(ab_value),
            std = sd(ab_value),
            abs_max = max(abs(ab_value))) %>%
  left_join(left_join(rbind(shap_model_industrialized, shap_model_noindustrialized) %>%
                        filter(taxon %in% c("Faecalibacterium", "Bifidobacterium")), metadata, by = "run_accession"), ., by = c("taxon")) %>%
  mutate(feature_value_z = (ab_value - mean) / std) %>%
  # filter(ab_value > 5e-5) %>%
  # filter(!taxon %in% c("")) %>%
  ggplot(., aes(y = feature_value_z, x = shap_value)) +
  geom_point(alpha = 0.5, size = .4) +# geom_bin2d(aes(fill = ..ndensity..), bins = 100)+
  # scale_fill_continuous(type = "viridis")+
  # geom_violin() +
  theme_classic() +
  geom_smooth() +
  facet_wrap(~taxon, scales = "free") +
  theme(panel.grid.major.y = element_line(color = "grey",
                                          size = 0.3,
                                          linetype = 2),
        strip.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16))


rbind(shap_model_industrialized, shap_model_noindustrialized) %>% 
  filter(taxon == c("unclassified_Enterobacteriaceae")) %>%
  left_join(., metadata, by = "run_accession") %>%
  group_by(taxon) %>%
  summarise(mean = mean(ab_value),
            std = sd(ab_value),
            abs_max = max(abs(ab_value))) %>%
  left_join(left_join(rbind(shap_model_industrialized, shap_model_noindustrialized) %>%
                        filter(taxon == c("unclassified_Enterobacteriaceae")), metadata, by = "run_accession"), ., by = c("taxon")) %>%
  mutate(feature_value_z = (ab_value - mean) / std) %>%
  # filter(ab_value > 5e-5) %>%
  # filter(!taxon %in% c("")) %>%
  ggplot(., aes(y = feature_value_z, x = shap_value)) +
  geom_point(alpha = 0.5, size = .4) +# geom_bin2d(aes(fill = ..ndensity..), bins = 100)+
  # scale_fill_continuous(type = "viridis")+
  # geom_violin() +
  theme_classic() +
  geom_smooth(method = "lm") +
  facet_wrap(~lifestyle , scales = "free") +
  theme(panel.grid.major.y = element_line(color = "grey",
                                          size = 0.3,
                                          linetype = 2),
        strip.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16))



# rarefied data ###############################
if(dataset_nested_cv_genus_step){
  lifestyle_predictions_genus_rare_h <- c("industrialized", "non_industrialized") %>% 
    future_map_dfr(~ get_predictions_lifestyle(. ,ps = ps_object_genus_raw %>% rarefy_even_depth(., sample.size = 2000)))
  save(lifestyle_predictions_genus_rare_h, file = "../data/inter_lifestyle_pred_genus_rare_h.RData")
} else {
  load("../data/inter_lifestyle_pred_genus_rare_h.RData")
}



# compare with prediction from same lifestyle:
load("../data/all_nested_cv_genus_industrialized_rare_h.RData")
genus_industrialized_nested_cv_rare_h_preds_long <- genus_industrialized_nested_cv_rare_h_preds %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "industrialized")

industrialized_lm_list_genus_rare_h <- get_lm_list(pred_df = genus_industrialized_nested_cv_rare_h_preds_long,
                                                   grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "industrialized")

load("../data/all_nested_cv_genus_nonindustrialized_rare_h.RData")
genus_nonindustrialized_nested_cv_rare_h_preds_long <- genus_nonindustrialized_nested_cv_rare_h_preds %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "non_industrialized")
nonindustrialized_lm_list_genus_rare_h <- get_lm_list(pred_df = genus_nonindustrialized_nested_cv_rare_h_preds_long, grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "non_industrialized")


# default model:
load("../data/all_nested_cv_genus_no_ls_rare_h.RData")
genus_nested_cv_rare_h_preds_long <- genus_no_ls_nested_rare_h_cv_preds %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "combined")
all_lm_list_genus_rare_h <- get_lm_list(pred_df = genus_nested_cv_rare_h_preds_long, grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "combined")


lifestyle_predictions_genus_rare_h <- lifestyle_predictions_genus_rare_h %>%
  mutate(model_name = "rf1") %>%
  mutate(training_set = ifelse(lifestyle == "industrialized", yes = "non_industrialized", no = "industrialized"))
lifestyle_lm_list_genus_rare_h <- get_lm_list(pred_df = lifestyle_predictions_genus_rare_h, grouping = c("study", "lifestyle", "model_name")) %>%
  mutate(training_set = ifelse(lifestyle == "industrialized", yes = "non_industrialized", no = "industrialized")) %>%
  rbind(., industrialized_lm_list_genus_rare_h, nonindustrialized_lm_list_genus_rare_h, all_lm_list_genus_rare_h)

lifestyle_lm_list_genus_rare_h <- lifestyle_lm_list_genus_rare_h %>%
  mutate(training_set = factor(training_set, levels = c("non_industrialized", "combined", "industrialized")),
         lifestyle = factor(lifestyle, levels = c("industrialized", "non_industrialized")))

save(lifestyle_lm_list_genus_rare_h,
     file = "../data/lifestyle_perfomance_genus_raref_h.RData")