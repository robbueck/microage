# create model on western or non-western studies and predict age in the respective other set of studies
setwd("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/")
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
library(vroom)
source("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/alt_models.R")
source("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/regression_functions.R")
source("/fast/AG_Forslund/rob/mm_index/R_scripts/setlists.R")

dataset_nested_cv_picrust_step <- T
shap_step <- T

# shap analysis:
pfun <- function(object, newdata) {
  require(ranger)
  predict(object, data = newdata)$predictions
}

get_shap_long <- function(model, test_data = NULL, features = NULL) {
  if(!is.null(features)) {
    if(all(features == "important")) {
      print("runnin shap analysis on the 15 most important features")
      features <- model$finalModel$variable.importance %>% sort %>% tail(., n = 15) %>% names
    }
  }
  if(is.null(test_data)) {
    test_data <- model$trainingData %>%
      select(-.outcome)
  }
  shap <- fastshap::explain(
    model$finalModel,
    X = test_data,
    pred_wrapper = pfun,
    nsim = 10,
    shap_only = F,
    feature_names = features,
    parallel = T)
  shap_long <- shap$shapley_values %>%
    data.frame() %>%
    mutate(sample_ID = 1:nrow(.)) %>%
    pivot_longer(cols = -sample_ID,
                 names_to = "taxon",
                 values_to = "shap_value")
  ab_long <- shap$feature_values %>%
    data.frame() %>%
    mutate(across(everything(), ~ . / max(.)))%>%
    mutate(sample_ID = 1:nrow(.)) %>%
    rownames_to_column(var = "run_accession") %>%
    pivot_longer(cols = -c(sample_ID, run_accession),
                 names_to = "taxon",
                 values_to = "ab_value")
  # ab-value: feature values divided by the maximum value of that feature => max = 1
  ab_shap_long <- left_join(shap_long, ab_long) 
  return(ab_shap_long)
}


get_predictions_lifestyle <- function(test_set, ps){
  print(test_set)
  oldDF <- as(sample_data(ps), "data.frame") 
  trainDF <- subset(oldDF, lifestyle != test_set)
  # %>%
  #   subset(., is.na(antibiotics_before) | !antibiotics_before) %>%
  #   subset(., is.na(antibiotics_one_week_before) | !antibiotics_one_week_before) %>%
  #   subset(., is.na(antibiotics_any) | !antibiotics_any)
  
  testDF <- subset(oldDF, lifestyle == test_set)
  ps_train <- ps
  ps_test <- ps
  rm(ps)
  sample_data(ps_train) <- sample_data(trainDF)
  sample_data(ps_test) <- sample_data(testDF)
  for_caret_list_train <- create_caret_df(ps_object = ps_train, transformation = "compositional",
                                          mean_ab_cutoff = 5e-5, study_prevalence_cutoff = 1,
                                          additional_cols = NULL,
                                          only_multiple_samples = F)
  for_caret_list_test <- create_caret_df(ps_object = ps_test, transformation = "compositional",
                                         additional_cols = NULL,
                                         only_multiple_samples = F,
                                         filter_features = colnames(for_caret_list_train$features))
  
  # subset features in test set:
  common_features <- colnames(for_caret_list_test$features) %in% colnames(for_caret_list_train$features)
  for_caret_list_test$features <- for_caret_list_test$features[,common_features]
  
  splits <- group_vfold_cv(for_caret_list_train$metadata, group = "study")
  caret_split <- rsample2caret(splits)
  tc_grouped <- trainControl(method = "cv",
                             savePredictions = "final",
                             selectionFunction = "oneSE", # function to select the best model
                             index = caret_split$index,
                             indexOut = caret_split$indexOut,
                             preProcOptions = list(uniqueCut = 5),
                             allowParallel = T)
  rfGrid <- expand.grid(mtry = c(3, 10, 20),
                        num.trees = c(200),
                        min.node.size = c(10, 20),
                        splitrule = c("variance"),
                        replace = F)
  rf_model <- train(x = for_caret_list_train$features, y = for_caret_list_train$metadata$age,
                    method = adapt_ranger, 
                    metric="Rsquared",
                    importance = "permutation",
                    trControl = tc_grouped,
                    preProcess = c("nzv"),
                    tuneGrid = rfGrid)
  predictions <- for_caret_list_test$metadata %>%
    mutate(pred = predict(rf_model, newdata = for_caret_list_test$features))
  return(predictions)
}

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

# picrust data ###################################################################
ps_object_family_raw <- readRDS("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_phyloseq_rf_filter_family.rds") %>%
  subset_samples(., age <= 730 & age > 1)

meta_df_picrust <- sample_data(ps_object_family_raw) %>% data.frame
# change rownames for bockulich
rownames(meta_df_picrust)[meta_df_picrust$study == "bockulich_2016"] <- 
  paste0("X", meta_df_picrust$sample_ID[meta_df_picrust$study == "bockulich_2016"])

picrust_ko_table <- vroom("/fast/AG_Forslund/rob/mm_index/merged_data/all/merged_picrust_ec.tsv", delim = "\t", col_names = T) %>%
  column_to_rownames(.,var = "function") %>%
  t() %>%
  as.data.frame()
picrust_ko_table <- picrust_ko_table[rownames(picrust_ko_table) %in% rownames(meta_df_picrust),] %>%
  round()

# create phyloseq object
picrust_ps <- phyloseq(otu_table(t(picrust_ko_table), taxa_are_rows = T), sample_data(meta_df_picrust)) %>%
  subset_samples(., age <= 730 & age > 1)

all_studies_picrust <- unique(picrust_ps@sam_data$study)

print("picrust data western datasets")
tic()
if(dataset_nested_cv_picrust_step){
  lifestyle_predictions_picrust <- c("westernized", "non_westernized") %>% 
    future_map_dfr(~ get_predictions_lifestyle(. ,ps = ps_object_picrust_raw))
  save(lifestyle_predictions_picrust, 
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/inter_lifestyle_pred_picrust.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/inter_lifestyle_pred_picrust.RData")
}
toc()
# compare with prediction from same lifestyle:
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_picrust_western.RData")
picrust_western_nested_cv_preds_long <- picrust_western_nested_cv_preds %>%
  pivot_longer(cols = c("rf1", "lasso", "glmnet"), names_to = "model_name", values_to = "pred" )
western_lm_list_picrust <- get_lm_list(pred_df = picrust_western_nested_cv_preds_long, grouping = c("study", "lifestyle")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "westernized")

load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_picrust_non_western.RData")
picrust_non_western_nested_cv_preds_long <- picrust_non_western_nested_cv_preds %>%
  pivot_longer(cols = c("rf1", "lasso", "glmnet"), names_to = "model_name", values_to = "pred" )
nonwestern_lm_list_picrust <- get_lm_list(pred_df = picrust_non_western_nested_cv_preds_long, grouping = c("study", "lifestyle")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "non_westernized")

# default model:
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_picrust.RData")
picrust_nested_cv_preds_long <- picrust_nested_cv_preds %>%
  pivot_longer(cols = c("rf1", "lasso", "glmnet"), names_to = "model_name", values_to = "pred" )
all_lm_list_picrust <- get_lm_list(pred_df = picrust_nested_cv_preds_long, grouping = c("study", "lifestyle")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "combined")


lifestyle_predictions_picrust <- lifestyle_predictions_picrust %>%
  mutate(model_name = "rf1")
lifestyle_lm_list_picrust <- get_lm_list(pred_df = lifestyle_predictions_picrust, grouping = c("study", "lifestyle")) %>%
  mutate(training_set = ifelse(lifestyle == "westernized", yes = "non_westernized", no = "westernized")) %>%
  rbind(., western_lm_list_picust, nonwestern_lm_list_picrust, all_lm_list_picrust)

lifestyle_lm_list_picrust <- lifestyle_lm_list_picrust %>%
  mutate(training_set = factor(training_set, levels = c("westernized", "combined", "non_westernized")),
         lifestyle = factor(lifestyle, levels = c("westernized", "non_westernized")))
# stat test
df_p_val_lifestyle_picrust <- lifestyle_lm_list_picrust %>%
  arrange(study) %>%
  mutate(training_set = factor(training_set)) %>%
  rstatix::group_by(lifestyle) %>%
  rstatix::wilcox_test(R2 ~ training_set, paired = T) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>% 
  rstatix::add_xy_position(x = "lifestyle", dodge = 0.8) 

ggplot(lifestyle_lm_list_picrust, aes(x=lifestyle, y = R2)) +
  geom_boxplot(aes(fill = training_set)) +
  ylim(0, 1) +
  xlab("Test set") +
  add_pvalue(df_p_val_lifestyle_picrust,
             label = "{p.adj.signif}",
             # step.group.by = "variation",
             step.increase = 0.05,
             tip.length = 0.01,
             # bracket.nudge.y = 0.02,
             xmin = "xmin", 
             xmax = "xmax",
             show.legend = FALSE) +
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
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/lifestyle_perfomance_picrust.pdf", width = 10)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/lifestyle_perfomance_picrust.png", width = 10)

# different feature importance #################################################
get_feature_importance <- function(study, lst = "nonwestern") {
  train_object <- readRDS(paste0("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/nested_cv_dataset_models/",
                                 "picrust_data_", lst, "_", study, ".rds"))
  var_imp <- train_object$rf1$finalModel$variable.importance %>%
    sort(., decreasing = T) %>%
    head(15)
  var_imp <- (var_imp / max(var_imp)) *100
  # Extract variable importance data
  data.frame(importance = var_imp,
             taxon = names(var_imp),
             ID = study)
}
western_studies <- lifestyle_lm_list_picrust %>% 
  filter(lifestyle == "westernized") %>%
  pull(study) %>% unique
non_western_studies <- lifestyle_lm_list_picrust %>% 
  filter(lifestyle == "non_westernized") %>%
  pull(study) %>% unique

combined_importance_nonwestern_picrust <- map(non_western_studies, function(x) get_feature_importance(x, lst = "nonwestern")) %>%
  purrr::reduce(rbind) %>%
  mutate(lifestyle = "non_westernized")
combined_importance_western_picrust <- map(western_studies, function(x) get_feature_importance(x, lst = "western")) %>%
  purrr::reduce(rbind) %>%
  mutate(lifestyle = "westernized")


combined_importance_lifestyles_picrust <- rbind(combined_importance_western_picrust, combined_importance_nonwestern_picrust) %>%
  mutate(lifestyle = factor(lifestyle, levels = c("westernized", "non_westernized")),
         # in_lit = ifelse(taxon %in% in_literature, yes = "red", no = "black"),
         taxon = reorder(taxon, importance))

saveRDS(combined_importance_lifestyles_picrust,
        "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/picrust_nested_cv_important_taxa_lifestyle.rds")

ggplot(combined_importance_lifestyles_picrust, aes(x = taxon, y = importance)) +
  # geom_line() + 
  geom_boxplot() +
  ylim(0,100) +
  ylab("Relative importance") +
  facet_wrap(~ lifestyle, nrow = 2, drop = T) +
  # geom_text(category, y = 0, size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "grey",
                                          size = 0.3,
                                          linetype = 2),
        strip.text = element_text(size = 14),
        plot.margin = margin(0.2,0.2,0.2,1, "cm"))
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/picrust_nested_cv_importance_lifestyle.pdf")
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/picrust_nested_cv_importance_lifestyle.png")



# shap analysis picrust ################################################################
# combined model
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_picrust.RData")
# final_model_pc, shap_out_pc
# western model
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_picrust_western.RData")
# final_model_pc_western, shap_out_pc_western
# non western model
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_picrust_non_western.RData")
# final_model_pc_non_western, shap_out_pc_non_western

# important features:
if(shap_step){
  shap_model_all <- get_shap_long(final_model_pc$rf1,
                                  features = levels(unique(combined_importance_lifestyles_picrust$taxon)))
  shap_model_west <- get_shap_long(final_model_pc_western$rf1, 
                                   test_data = final_model_pc$rf1$trainingData %>%
                                     select(-.outcome),
                                   features = levels(unique(combined_importance_lifestyles_picrust$taxon)))
  shap_model_nowest <- get_shap_long(final_model_pc_non_western$rf1,
                                     test_data = final_model_pc$rf1$trainingData %>%
                                       select(-.outcome),
                                     features = levels(unique(combined_importance_lifestyles_picrust$taxon)))
  save(shap_model_all, shap_model_west, shap_model_nowest,
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/lifestyle_shap_data_picrust.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/lifestyle_shap_data_picrust.RData")
}


shap_model_all %>%
  # filter(ab_value > 5e-5) %>%
  # filter(!taxon %in% c("")) %>%
  filter(!taxon %in% c("Observed", "Shannon")) %>%
  left_join(., metadata, by = "run_accession") %>%
  ggplot(., aes(y = taxon, x = shap_value)) +
  # geom_jitter(size = 0.1,alpha = 0.5) +
  geom_bin2d(aes(fill = ..ndensity..), bins = 100)+
  scale_fill_continuous(type = "viridis")+
  # geom_violin() +
  geom_vline(xintercept = 0) +
  theme_classic() +
  facet_wrap(~lifestyle) +
  theme(panel.grid.major.y = element_line(color = "grey",
                                          size = 0.3,
                                          linetype = 2))+
  ggtitle("Combined_model")
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_combined_model_picrust.pdf")
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_combined_model_picrust.png")


shap_model_all %>%
  left_join(., metadata, by = "run_accession") %>%
  group_by(taxon, lifestyle) %>% 
  summarise(mean = mean(ab_value),
            std = sd(ab_value),
            abs_max = max(abs(ab_value))) %>%
  left_join(left_join(shap_model_all, metadata, by = "run_accession"), ., by = c("taxon", "lifestyle")) %>%
  mutate(feature_value_z = ab_value / abs_max) %>%
  # filter(ab_value > 5e-5) %>%
  # filter(!taxon %in% c("")) %>%
  ggplot(., aes(y = taxon, x = shap_value, color = feature_value_z)) +
  geom_jitter(size = 0.0002,alpha = 0.3) +
  # geom_bin2d(aes(fill = ..ndensity..), bins = 100)+
  # scale_fill_continuous(type = "viridis")+
  # geom_violin() +
  geom_vline(xintercept = 0) +
  theme_classic() +
  facet_wrap(~lifestyle) +
  scale_colour_gradient(low = "blue", high = "red", na.value = NA) +
  theme(panel.grid.major.y = element_line(color = "grey",
                                          size = 0.3,
                                          linetype = 2))+
  ggtitle("Combined_model_picrust")



shap_model_all %>%
  left_join(., metadata, by = "run_accession") %>%
  group_by(taxon, lifestyle) %>% 
  summarise(R2 = cor(ab_value, shap_value)) %>%
  filter(!taxon %in% c("Observed", "Shannon")) %>%
  ggplot(., aes(y = taxon, x = R2, fill = lifestyle, group = lifestyle)) +
  geom_bar(position="dodge", stat="identity") +
  theme_classic() +
  theme(panel.grid.major.y = element_line(color = "grey",
                                          size = 0.3,
                                          linetype = 2))+
  ggtitle("Combined_model_picrust")





shap_model_west %>%
  # filter(abs(shap_value) > 0.01,
  #        ab_value > 5e-5) %>%
  filter(!taxon %in% c("Observed", "Shannon")) %>%
  left_join(., metadata, by = "run_accession") %>%
  ggplot(., aes(y = taxon, x = shap_value)) +
  # geom_jitter(size = 0.1, alpha = 0.5) +
  geom_bin2d(aes(fill = ..ndensity..), bins = 100)+
  scale_fill_continuous(type = "viridis")+
  # geom_violin() +
  geom_vline(xintercept = 0) +
  theme_classic() +
  facet_wrap(~lifestyle) +
  theme(panel.grid.major.y = element_line(color = "grey",
                                          size = 0.3,
                                          linetype = 2))+
  ggtitle("Westernized model picrust")
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_western_model_picrust.pdf")
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_western_model_picrust.png")



shap_model_nowest %>%
  # filter(abs(shap_value) > 0.1,
  #        ab_value > 5e-3) %>%
  filter(!taxon %in% c("Observed", "Shannon")) %>%
  left_join(., metadata, by = "run_accession") %>%
  ggplot(., aes(y = taxon, x = shap_value)) +
  # geom_violin() +
  geom_bin2d(aes(fill = ..ndensity..), bins = 100)+
  scale_fill_continuous(type = "viridis")+
  # geom_jitter(size = 0.4, alpha = 0.1) +
  # geom_boxplot(alpha = 0.3) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  facet_wrap(~lifestyle) +
  theme(panel.grid.major.y = element_line(color = "grey",
                                                  size = 0.3,
                                                  linetype = 2))+
  ggtitle(" Non-westernized Model picrust")
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_non_western_model_picrust.pdf")
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_non_western_model_picrust.png")

shap_model_west$model <- "westernized_model"
shap_model_nowest$model <- "non_westernized_model"
rbind(shap_model_west, shap_model_nowest) %>%
  filter(!taxon %in% c("Observed", "Shannon")) %>%
  left_join(., metadata, by = "run_accession") %>%
  mutate(lifestyle = ifelse(lifestyle == "westernized", yes = "westernized_samples", no = "non_westernized_samples")) %>%
  ggplot(., aes(y = taxon, x = shap_value)) +
  geom_bin2d(aes(fill = ..ndensity..), bins = 100)+
  scale_fill_continuous(type = "viridis")+
  # geom_jitter(size = 0.4, alpha = 0.1) +
  # geom_boxplot(alpha = 0.3) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  facet_wrap(~lifestyle + model) +
  theme(panel.grid.major.y = element_line(color = "grey",
                                          size = 0.3,
                                          linetype = 2))
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_lifestyle_models_picrust.pdf")
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_lifestyle_models_picrust.png")



sha_value_cor <- rbind(shap_model_west, shap_model_nowest) %>%
  left_join(., metadata, by = "run_accession") %>%
  group_by(taxon, lifestyle, model) %>% 
  summarise(R2 = cor(ab_value, shap_value, method = "spearman")) %>%
  filter(!taxon %in% c("X.Eubacterium..eligens.group", "UCG.002", "Observed",
                       "Shannon", "Lachnospiraceae.ND3007.group", "Intestinibacter",
                       "Erysipelatoclostridium", "Family.XIII.AD3011.group",
                       "Lachnospiraceae.ND3007.group","Anaerostipes", "Oscillibacter")) %>%
  ggplot(., aes(y = taxon, x = R2, fill = lifestyle, group = lifestyle)) +
  geom_bar(position="dodge", stat="identity") +
  theme_classic() +
  facet_grid(~model) +
  xlab("Spearman relative feature abundances vs. SHAP-values") +
  theme(panel.grid.major.y = element_line(color = "grey",
                                          size = 0.3,
                                          linetype = 2),
        strip.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_blank())
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_lifestyle_models_cor_picrust.pdf",
       width = 9, height = 10, plot = sha_value_cor)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_lifestyle_models_cor_picrust.png")




shap_vs_value <- rbind(shap_model_west, shap_model_nowest) %>%
  left_join(., metadata, by = "run_accession") %>%
  group_by(taxon, lifestyle) %>% 
  summarise(mean = mean(ab_value),
            std = sd(ab_value),
            abs_max = max(abs(ab_value))) %>%
  left_join(left_join(rbind(shap_model_west, shap_model_nowest), metadata, by = "run_accession"), ., by = c("taxon", "lifestyle")) %>%
  mutate(feature_value_z = ab_value / abs_max) %>%
  filter(!taxon %in% c("X.Eubacterium..eligens.group", "UCG.002", "Observed",
                       "Shannon", "Lachnospiraceae.ND3007.group", "Intestinibacter",
                       "Erysipelatoclostridium", "Family.XIII.AD3011.group",
                       "Lachnospiraceae.ND3007.group", "Anaerostipes", "Oscillibacter")) %>%
  # filter(ab_value > 5e-5) %>%
  # filter(!taxon %in% c("")) %>%
  ggplot(., aes(y = taxon, x = shap_value, color = feature_value_z)) +
  geom_jitter(size = 0.0002,alpha = 0.3) +
  # geom_bin2d(aes(fill = ..ndensity..), bins = 100)+
  # scale_fill_continuous(type = "viridis")+
  # geom_violin() +
  geom_vline(xintercept = 0) +
  theme_classic() +
  facet_wrap(~model) +
  scale_colour_gradient(low = "blue", high = "red", na.value = NA) +
  theme(panel.grid.major.y = element_line(color = "grey",
                                          size = 0.3,
                                          linetype = 2),
        strip.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_blank())
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_lifestyle_models_ab_value_picrust.pdf")
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_lifestyle_models_ab_value_picrusts.png",
       width = 9, height = 10,
       plot = shap_vs_value)



rbind(shap_model_west, shap_model_nowest) %>%
  filter(taxon %in% c("Faecalibacterium", "Bifidobacterium")) %>%
  left_join(., metadata, by = "run_accession") %>%
  group_by(taxon) %>%
  summarise(mean = mean(ab_value),
            std = sd(ab_value),
            abs_max = max(abs(ab_value))) %>%
  left_join(left_join(rbind(shap_model_west, shap_model_nowest) %>%
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
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_correlation_picrust.pdf", width = 16, height = 8)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_correlation_picrust.png",
       width = 9, height = 10)


