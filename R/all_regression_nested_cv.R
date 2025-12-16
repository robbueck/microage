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
library(furrr)
library(vroom)
library(ranger)
library(purrr)
library(fastshap)
library(foreach)
library(ggpubr)
library(Boruta)
library(gamlss)
library(lmtest)

source("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/alt_models.R")
source("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/regression_functions.R")
source("/fast/AG_Forslund/rob/mm_index/R_scripts/setlists.R")


# Switches ################################
nested_cv_family_step <- F
final_model_family_step <- F
nested_cv_extra_metadata_step <- F
final_model_family_extra_metadata_step <- F
nested_cv_family_no_ls_step <- F
final_model_family_no_ls_step <- F
nested_cv_family_industrialized_step <- F
final_model_family_industrialized_step <- F
nested_cv_family_nonindustrialized_step <- F
final_model_family_nonindustrialized_step <- F
# genus
nested_cv_genus_no_ls_step <- F
final_model_genus_no_ls_step <- F
nested_cv_genus_industrialized_step <- F
final_model_genus_industrialized_step <- F
nested_cv_genus_nonindustrialized_step <- F
final_model_genus_nonindustrialized_step <- F
# rarefaction
nested_cv_genus_no_ls_rarefied_step <- F
final_model_genus_no_ls_rarefied_step <- F
nested_cv_genus_industrialized_rarefied_step <- F
final_model_genus_industrialized_rarefied_step <- F
nested_cv_genus_nonindustrialized_rarefied_step <- F
final_model_genus_nonindustrialized_rarefied_step <- F

nested_cv_genus_no_ls_rarefied_l_step <- T
final_model_genus_no_ls_rarefied_l_step <- T
nested_cv_genus_industrialized_rarefied_l_step <- T
final_model_genus_industrialized_rarefied_l_step <- T
nested_cv_genus_nonindustrialized_rarefied_l_step <- T
final_model_genus_nonindustrialized_rarefied_l_step <- T



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
  important_features <-
    model$finalModel$variable.importance %>% sort %>% tail(., n = 15)
  return(ab_shap_long)
}



get_final_model <- function(ps, transform = "compositional",
                            extra_cols = c("Observed", "Shannon", "lifestyle_industrialized"),
                            prefilter = F,
                            rf_depth = 2000,
                            pre_processing = "nzv",
                            pre_processing_option = list(uniqueCut = 5)) {
  for_caret_list_train <- create_caret_df(ps_object = ps, transformation = transform,
                                          mean_ab_cutoff = 5e-5, study_prevalence_cutoff = 2,
                                          prevalence_in_study_cutoff = 5,
                                          additional_cols = extra_cols,
                                          rf_depth = rf_depth,
                                          only_multiple_samples = F)
  if(prefilter) {
    family_cleaned_cor <- cor(for_caret_list_train$features)
    highlyCorrelated_family_cleaned <- findCorrelation(family_cleaned_cor, cutoff=0.8)
    for_caret_list_train$features <- for_caret_list_train$features[,-highlyCorrelated_family_cleaned]
  }
  splits <- group_vfold_cv(for_caret_list_train$metadata, group = "study")
  caret_split <- rsample2caret(splits)
  tc_grouped <- trainControl(method = "cv",
                             savePredictions = "final",
                             selectionFunction = "oneSE", # function to select the best model
                             index = caret_split$index,
                             indexOut = caret_split$indexOut,
                             preProcOptions = pre_processing_option,
                             returnData = T,
                             allowParallel = T)
  rf_model_list <- caretList(y=for_caret_list_train$metadata$age, x=for_caret_list_train$features,
                             metric="Rsquared",
                             trControl=tc_grouped,
                             preProcess = pre_processing,
                             methodList=c("lasso"),
                             tuneList = list(
                               rf1=caretModelSpec(method=adapt_ranger, tuneGrid=expand.grid(mtry = c(3, 10),
                                                                                            num.trees = c(200),
                                                                                            min.node.size = c(5, 10),
                                                                                            splitrule = c("variance"),
                                                                                            replace = F),
                                                  importance = "permutation")
                             )
  )
  return(rf_model_list)
}

get_predictions_nested_cv <- function(test_set, ps, prefix = "default",
                                      transform = "compositional",
                                      extra_cols = c("Observed", "Shannon", "lifestyle_industrialized"),
                                      prefilter = F,
                                      pre_processing = "nzv",
                                      rf_depth = 2000,
                                      pre_processing_option = list(uniqueCut = 5)){
  print(test_set)
  oldDF <- as(sample_data(ps), "data.frame") 
  trainDF <- subset(oldDF, study != test_set)
  testDF <- subset(oldDF, study == test_set)
  ps_train <- ps
  ps_test <- ps
  rm(ps)
  sample_data(ps_train) <- sample_data(trainDF)
  sample_data(ps_test) <- sample_data(testDF)
  for_caret_list_train <- create_caret_df(ps_object = ps_train, transformation = transform,
                                          mean_ab_cutoff = 5e-5, study_prevalence_cutoff = 2,
                                          prevalence_in_study_cutoff = 5,
                                          rf_depth = rf_depth,
                                          additional_cols = extra_cols,
                                          only_multiple_samples = F)
  if(prefilter) {
    family_cleaned_cor <- cor(for_caret_list_train$features)
    highlyCorrelated_family_cleaned <- findCorrelation(family_cleaned_cor, cutoff=0.8)
    for_caret_list_train$features <- for_caret_list_train$features[,-highlyCorrelated_family_cleaned]
  }
  for_caret_list_test <- create_caret_df(ps_object = ps_test, transformation = transform,
                                         additional_cols = extra_cols,
                                         only_multiple_samples = F,
                                         rf_depth = rf_depth,
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
                             preProcOptions = pre_processing_option,
                             allowParallel = T)
  possibleError <- try(
    expr = {
      rf_model_list <- caretList(y=for_caret_list_train$metadata$age, x=for_caret_list_train$features,
                                 metric="Rsquared",
                                 trControl=tc_grouped,
                                 preProcess = pre_processing,
                                 methodList=c("lasso"),
                                 tuneList = list(
                                   rf1=caretModelSpec(method=adapt_ranger, tuneGrid=expand.grid(mtry = c(3, 10),
                                                                                                num.trees = c(50, 200),
                                                                                                min.node.size = c(5, 10),
                                                                                                splitrule = c("variance"),
                                                                                                replace = F),
                                                      num.threads = 10,
                                                      importance = "permutation")
                                 )
      )
      write_rds(rf_model_list,
                file = paste0("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/nested_cv_dataset_models/",
                              prefix, "_", test_set, ".rds"))
      # Boruta feature importance
      print("Run Boruta")
      sel_feat <- Boruta(x=for_caret_list_train$features, y = for_caret_list_train$metadata$age,
                         num.trees = rf_model_list$rf1$bestTune$num.trees, 
                         mtry = rf_model_list$rf1$bestTune$mtry, num.threads = ceiling(n_cores/5))
      write_rds(sel_feat,
                file = paste0("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/nested_cv_dataset_models/",
                              prefix, "_", test_set, "_boruta_res.rds"))
      
      predictions <- cbind(for_caret_list_test$metadata,
                           predict(rf_model_list, newdata = for_caret_list_test$features))
    }
  )
  if(inherits(possibleError, "try-error")){
    message("No model fit possible for: ", test_set)
    predictions <- mutate(for_caret_list_test$metadata,
                          lasso = NA,
                          rf1 = NA)
  }
  
  return(predictions)
}

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
set.seed(825)
future::plan(multisession, workers = ceiling(n_cores/2))



# family data with extra metadata ##############################################
ps_object_family_raw <- readRDS("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_phyloseq_rf_filter_family.rds") %>%
  subset_samples(., age <= 730 & age > 1) 

all_studies <- unique(ps_object_family_raw@sam_data$study)


print("family data with extra metadata")

sample_data(ps_object_family_raw) <- sample_data(ps_object_family_raw) %>%
  data.frame() %>%
  mutate(birthmode_vd = case_when(birthmode == "v" ~ 1,
                                  birthmode == "c" ~ 0,
                                  is.na(birthmode) ~ 3),
         sex_f = case_when(sex == "f" ~ 1,
                           sex == "female" ~ 1,
                           sex == "m" ~ 0,
                           is.na(sex) ~ 3)) %>%
  sample_data()

tic()
if(nested_cv_extra_metadata_step){
  family_nested_extra_d_cv_preds <- all_studies %>% 
    future_map_dfr(~ get_predictions_nested_cv(. ,ps = ps_object_family_raw,
                                               prefix = "family_data_extra_d",
                                               extra_cols = c("Observed", "Shannon",
                                                              "lifestyle_industrialized",
                                                              "sex_f", "birthmode_vd"),
                                               transform = "compositional"), progress = F)
  family_nested_extra_d_cv_preds_long <- family_nested_extra_d_cv_preds %>%
    pivot_longer(cols = c("rf1", "lasso"), names_to = "model_name", values_to = "pred")
  ps_nested_cv_extra_d <- modelplots(all_results = family_nested_extra_d_cv_preds_long)
  save(family_nested_extra_d_cv_preds, ps_nested_cv_extra_d,
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_family_extra_d.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_family_extra_d.RData")
}


ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/family_extra_d_nested_cv_point_binary.pdf",
       plot = ps_nested_cv_extra_d$correlation)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/family_extra_d_nested_cv_heatmap_binary.pdf",
       plot = ps_nested_cv_extra_d$model_heatmap)



if(final_model_family_extra_metadata_step){
  final_model_family_extra_d <- get_final_model(ps = ps_object_family_raw,
                                        extra_cols = c("Observed", "Shannon",
                                                       "lifestyle_industrialized",
                                                       "sex_f", "birthmode_vd"))
  shap_out_family_extra_d <- get_shap_long(final_model_family_extra_d$rf1, features = "important")
  shap_out_family_extra_d %>%
    filter(abs(shap_value) > 0.01,
           ab_value > 5e-5) %>%
    ggplot(., aes(y = taxon, x = shap_value, color = ab_value)) +
    geom_jitter(size = 0.5) +
    geom_violin() +
    geom_vline(xintercept = 0) +
    theme_classic()
  ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/family_extra_d_all_shap.pdf")
  save(final_model_family_extra_d, shap_out_family_extra_d,
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_all_family_extra_d.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_all_family_extra_d.RData")
}

# shap_out_family_extra_d %>%
#   left_join(.,data.frame(ps_object_family_raw@sam_data), by = "run_accession") %>%
#   filter(abs(shap_value) > 0.01,
#          ab_value > 5e-5) %>%
#   ggplot(., aes(y = taxon, x = shap_value, color = lifestyle)) +
#   geom_violin() +
#   geom_vline(xintercept = 0) +
#   theme_classic()
# 
# 
# ps_object_family_raw %>% sample_data() %>%
#   data.frame() %>%
#   ggplot(., aes(x = age, color = as.factor(birthmode_vd))) +
#   geom_histogram() +
#   facet_grid(~lifestyle, scales = "free_x")

# family data without lifestyle ################################################

studies_all <- unique(ps_object_family_raw@sam_data$study)
print("family data nolifestyle samples")
tic()
if(nested_cv_family_no_ls_step){
  family_no_ls_nested_cv_preds <- studies_all %>% 
    future_map_dfr(~ get_predictions_nested_cv(. ,ps = ps_object_family_raw,
                                               prefix = "family_data_no_ls",
                                               transform = "compositional",
                                               extra_cols = c("Observed", "Shannon")), progress = F)
  family_no_ls_nested_cv_preds_long <- family_no_ls_nested_cv_preds %>%
    pivot_longer(cols = c("rf1", "lasso"), names_to = "model_name", values_to = "pred" )
  ps_nested_cv_no_ls <- modelplots(all_results = family_no_ls_nested_cv_preds_long)
  save(family_no_ls_nested_cv_preds, ps_nested_cv_no_ls,
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_family_no_ls.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_family_no_ls.RData")
}
toc()


if(final_model_family_no_ls_step){
  final_model_family_no_ls <- get_final_model(ps = ps_object_family_raw, extra_cols = c("Observed", "Shannon"))
  shap_out_family_no_ls <- get_shap_long(final_model_family_no_ls$rf1, features = "important")
  shap_out_family_no_ls %>%
    filter(abs(shap_value) > 0.01,
           ab_value > 5e-5) %>%
    ggplot(., aes(y = taxon, x = shap_value, color = ab_value)) +
    geom_jitter(size = 0.5) +
    geom_violin() +
    geom_vline(xintercept = 0) +
    theme_classic()
  ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/family_no_ls_shap.pdf")
  save(final_model_family_no_ls, shap_out_family_no_ls, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_no_ls_family.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_no_ls_family.RData")
}

shap_out_family_no_ls %>%
  filter(abs(shap_value) > 0.01,
         ab_value > 5e-5) %>% 
  ggplot(., aes(y = taxon, x = shap_value, color = ab_value)) +
  geom_jitter(size = 0.5) +
  geom_violin() +
  geom_vline(xintercept = 0) +
  theme_classic()



ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/family_nested_cv_point_no_ls.pdf",
       plot = ps_nested_cv_no_ls$correlation)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/family_nested_cv_heatmap_no_ls.pdf",
       plot = ps_nested_cv_no_ls$model_heatmap)

# check feature importances:
create_feature_importance_df <- function(study) {
  train_object <- readRDS(paste0("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/nested_cv_dataset_models/",
                                 "family_data_no_ls_", study, ".rds"))
  var_imp <- train_object$rf1$finalModel$variable.importance %>%
    sort(., decreasing = T) %>%
    head(15)
  var_imp <- (var_imp / max(var_imp)) *100
  # Extract variable importance data
  data.frame(importance = var_imp,
             taxon = names(var_imp),
             ID = study)
}
# combined_importance_no_ls <- map(studies_all, create_feature_importance_df) %>%
#   purrr::reduce(rbind)  # Combine the individual plots
# 
# combined_importance_no_ls <- combined_importance_no_ls %>%
#   mutate(taxon = gsub("X.Eubacterium..coprostanoligenes.group", "E.coprostanoligenes_group", taxon))
# saveRDS(combined_importance_no_ls,
#         "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/family_nested_cv_important_taxa_no_ls.rds")



# family data for industrialized populations only specific model ##################

ps_object_family_raw_industrialized <- ps_object_family_raw %>%
  subset_samples(lifestyle == "industrialized")
industrialized_studies <- unique(ps_object_family_raw_industrialized@sam_data$study)

print("family data only for industrialized samples")
tic()
if(nested_cv_family_industrialized_step){
  family_industrialized_nested_cv_preds <- industrialized_studies %>% 
    future_map_dfr(~ get_predictions_nested_cv(. ,ps = ps_object_family_raw_industrialized,
                                               prefix = "family_data_industrialized",
                                               transform = "compositional",
                                               extra_cols = c("Observed", "Shannon")), progress = F)
  family_industrialized_nested_cv_preds_long <- family_industrialized_nested_cv_preds %>%
    pivot_longer(cols = c("rf1", "lasso"), names_to = "model_name", values_to = "pred" )
  ps_nested_cv_industrialized <- modelplots(all_results = family_industrialized_nested_cv_preds_long)
  save(family_industrialized_nested_cv_preds, ps_nested_cv_industrialized,
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_family_industrialized.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_family_industrialized.RData")
}
toc()


if(final_model_family_industrialized_step){
  final_model_family_industrialized <- get_final_model(ps = ps_object_family_raw_industrialized, extra_cols = c("Observed", "Shannon"))
  shap_out_family_industrialized <- get_shap_long(final_model_family_industrialized$rf1, features = "important")
  shap_out_family_industrialized %>%
    filter(abs(shap_value) > 0.01,
           ab_value > 5e-5) %>%
    ggplot(., aes(y = taxon, x = shap_value, color = ab_value)) +
    geom_jitter(size = 0.5) +
    geom_violin() +
    geom_vline(xintercept = 0) +
    theme_classic()
  ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/family_industrialized_shap.pdf")
  save(final_model_family_industrialized, shap_out_family_industrialized, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_industrialized_family.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_industrialized_family.RData")
}


ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/family_nested_cv_point_industrialized.pdf",
       plot = ps_nested_cv_industrialized$correlation)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/family_nested_cv_heatmap_industrialized.pdf",
       plot = ps_nested_cv_industrialized$model_heatmap)

# check feature importances:
create_feature_importance_df <- function(study) {
  train_object <- readRDS(paste0("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/nested_cv_dataset_models/",
                                 "family_data_industrialized_", study, ".rds"))
  var_imp <- train_object$rf1$finalModel$variable.importance %>%
    sort(., decreasing = T) %>%
    head(15)
  var_imp <- (var_imp / max(var_imp)) *100
  # Extract variable importance data
  data.frame(importance = var_imp,
             taxon = names(var_imp),
             ID = study)
}
# combined_importance_industrialized <- map(industrialized_studies, create_feature_importance_df) %>%
#   purrr::reduce(rbind)  # Combine the individual plots
#   
# combined_importance_industrialized <- combined_importance_industrialized %>%
#   mutate(taxon = gsub("X.Eubacterium..coprostanoligenes.group", "E.coprostanoligenes_group", taxon),
#          taxon = gsub("Bacteria_Firmicutes_Clostridia_unclassified_Clostridia_unclassified_Clostridia", 
#                       "unclassified_Clostridia", taxon),
#          lifestyle = "industrialized")
# 
# in_literature <- c("Ruminococcaceae", "Lachnospiraceae", "Oscillospiraceae", "Bifidobacteriaceae",
#                    "Peptostreptococcaceae", "Erysipelotrichaceae", "Enterobacteriaceae", "Bacteroidaceae")
# in_lit_col <- c("black", "red", "black", "black", "black", "black", "red", "black", 
#                 "black", "black", "red", "red", "red", "black", "black", 
#                 "red", "black", "red", "red", "red")
# ggplot(combined_importance_industrialized, aes(x = reorder(taxon, importance), y = importance)) +
#   # geom_line() + 
#   geom_boxplot() +
#   ylim(0,100) +
#   ylab("Relative importance") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14, colour = in_lit_col),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(size = 16),
#         axis.text.y = element_text(size = 14),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank())
# 
# # important industrialized features:
# c("Alistipes", "Anaerostipes", "Bacteroides", "Bifidobacteriales",
#   "Bifidobacterium", "Blautia", "Butyricicoccus", "Clostridia",
#   "Clostridioides_difficile", "Clostridiales", "Clostridiales_bacterium_VE202-13",
#   "Clostridium", "Eubacterium", "Enterobacter", "Escherichia",
#   "Erysipelatoclostridium", "Erysipelotrichaceae", "Faecalibacterium",
#   "Faecalibacterium_prausnitzii", "Firmicutes", "Flavonifractor",
#   "Gordonibacter", "Hungatella", "Intestinibacter", "Lachnospiraceae",
#   "Lachnospiraceae_bacterium_3-1", "Oscillibacter", "Oscillospiraceae",
#   "Propionibacterium", "Roseburia_inulinivorans_DSM_16841",
#   "Ruminococcaceae", "Ruminococcus", "Sellimonas", "Streptococcus",
#   "Subdoligranulum", "Terrabacteria_group", "Turicibacter", "Veilonella")
# # as families:
# c("Bacteroidaceae", "Bifidobacteriaceae", "Clostridiaceae", "Eggerthellaceae",
#   "Enterobacteriaceae", "Erysipelotrichaceae", "Lachnospiraceae", "Oscillospiraceae",
#   "Peptostreptococcaceae", "Propionibacteriaceae", "Ruminococcaceae", "Streptococcaceae",
#   "Veillonellaceae")

# family data for non industrialized populaitons only #############################

ps_object_family_raw_nonindustrialized <- ps_object_family_raw %>%
  subset_samples(lifestyle != "industrialized")
nonindustrialized_studies <- unique(ps_object_family_raw_nonindustrialized@sam_data$study)


print("family data only for non-industrialized samples")
tic()
if(nested_cv_family_nonindustrialized_step){
  family_nonindustrialized_nested_cv_preds <- nonindustrialized_studies %>% 
    future_map_dfr(~ get_predictions_nested_cv(. ,ps = ps_object_family_raw_nonindustrialized,
                                               prefix = "family_data_nonindustrialized",
                                               transform = "compositional",
                                               extra_cols = c("Observed", "Shannon")), progress = F)
  family_nonindustrialized_nested_cv_preds_long <- family_nonindustrialized_nested_cv_preds %>%
    pivot_longer(cols = c("rf1", "lasso"), names_to = "model_name", values_to = "pred" )
  ps_nested_cv_nonindustrialized <- modelplots(all_results = family_nonindustrialized_nested_cv_preds_long)
  save(family_nonindustrialized_nested_cv_preds, ps_nested_cv_nonindustrialized,
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_family_nonindustrialized.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_family_nonindustrialized.RData")
}
toc()

if(final_model_family_nonindustrialized_step){
  final_model_family_nonindustrialized <- get_final_model(ps = ps_object_family_raw_nonindustrialized, extra_cols = c("Observed", "Shannon"))
  shap_out_family_nonindustrialized <- get_shap_long(final_model_family_nonindustrialized$rf1, features = "important")
  shap_out_family_nonindustrialized %>%
    filter(abs(shap_value) > 0.01,
           ab_value > 5e-5) %>%
    ggplot(., aes(y = taxon, x = shap_value, color = ab_value)) +
    geom_jitter(size = 0.5) +
    geom_violin() +
    geom_vline(xintercept = 0) +
    theme_classic()
  ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/family_nonindustrialized_shap.pdf")
  save(final_model_family_nonindustrialized, shap_out_family_nonindustrialized, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_nonindustrialized_family.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_nonindustrialized_family.RData")
}

ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/family_nested_cv_point_nowest.pdf",
       plot = ps_nested_cv_nonindustrialized$correlation)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/family_nested_cv_heatmap_nowest.pdf",
       plot = ps_nested_cv_nonindustrialized$model_heatmap)

# check feature importances:
create_feature_importance_df_nw <- function(study) {
  train_object <- readRDS(paste0("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/nested_cv_dataset_models/",
                                 "family_data_nonindustrialized_", study, ".rds"))
  var_imp <- train_object$rf1$finalModel$variable.importance %>%
    sort(., decreasing = T) %>%
    head(15)
  var_imp <- (var_imp / max(var_imp)) *100
  # Extract variable importance data
  data.frame(importance = var_imp,
             taxon = names(var_imp),
             ID = study)
}
# combined_importance_nonindustrialized <- map(nonindustrialized_studies, create_feature_importance_df_nw) %>%
#   purrr::reduce(rbind)  # Combine the individual plots
# 
# combined_importance_nonindustrialized <- combined_importance_nonindustrialized %>%
#   mutate(taxon = gsub("X.Eubacterium..coprostanoligenes.group", "E.coprostanoligenes_group", taxon),
#          taxon = gsub("Bacteria_Firmicutes_Clostridia_unclassified_Clostridia_unclassified_Clostridia", 
#                       "unclassified_Clostridia", taxon),
#          lifestyle = "non_industrialized")
# 
# in_literature <- c("Bacteroidaceae", "Bifidobacteriaceae", "Eggerthellaceae", "Enterobacteriaceae",
#                    "Erysipelotrichaceae", "Lachnospiraceae", "Lactobacillaceae", "Oscillospiraceae",
#                    "Pasteurellaceae", "Peptostreptococcaceae", "Prevotellaceae", "Ruminococcaceae")
# in_lit_col <- c("black", "black", "black", "red", "black", "red", "black", "red", "black", "black", "black", "black",
#                 "red", "black", "red", "red", "red", "red", "red", "red",
#                 "black", "black", "black", "red", "black", "red")
# 
# combined_importance_lifestyles <- rbind(combined_importance_industrialized, combined_importance_nonindustrialized) %>%
#   mutate(lifestyle = factor(lifestyle, levels = c("industrialized", "non_industrialized")),
#          measure_type = ifelse(taxon %in% c("Shannon", "Observed"), yes = "a_div", no = "taxon"),
#          in_lit = ifelse(taxon %in% in_literature, yes = "red", no = "black"),
#          taxon = reorder(taxon, importance),
#          taxon = factor(taxon, levels = unique(c("Shannon", "Observed", levels(taxon)))))
# 
# saveRDS(combined_importance_lifestyles, "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/family_nested_cv_important_taxa_lifestyle.rds")
# 
# ggplot(combined_importance_lifestyles, aes(x = taxon, y = importance)) +
#   # geom_line() + 
#   geom_boxplot() +
#   ylim(0,100) +
#   ylab("Relative importance") +
#   facet_wrap(~ lifestyle, nrow = 2, drop = T) +
#   # geom_text(category, y = 0, size = 3) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14, colour = in_lit_col),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(size = 16),
#         axis.text.y = element_text(size = 14),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major.x = element_line(color = "grey",
#                                           size = 0.3,
#                                           linetype = 2),
#         strip.text = element_text(size = 14),
#         plot.margin = margin(0.2,0.2,0.2,1, "cm"))
# ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/family_nested_cv_importance_lifestyle.pdf")
# 
# # non industrialized:
# c("Bacteroides", "Bifidobacterium", "Blautia", "Catenibacterium",
#   "Clostridium", "Dialister", "Dorea", "Eubacterium", "E. coli",
#   "Faecalibacterium", "Haemophilus", "Intestinibacter", "Lactobacillus",
#   "Pasteurellaceae", "Prevotella", "Ruminococcus", "Ruminococcaceae",
#   "Streptococcus", "Subdoligranulum", "Veilonella", "Weissella")
# 
# # combined:
# c("Bacteroidaceae", "Bifidobacteriaceae", "Clostridiaceae", "Eggerthellaceae",
#   "Enterobacteriaceae", "Erysipelotrichaceae", "Lachnospiraceae", "Oscillospiraceae",
#   "Peptostreptococcaceae", "Propionibacteriaceae", "Ruminococcaceae", "Streptococcaceae",
#   "Veillonellaceae", "Lactobacillaceae", "Pasteurellaceae", "Prevotellaceae", 
#   "Leuconostocaceae")


# genus data without lifestyle #################################################

print("genus data")
ps_object_genus_raw <- readRDS("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_phyloseq_rf_filter_genus.rds") %>%
  subset_samples(., age <= 730 & age > 1)
all_studies <- unique(ps_object_genus_raw@sam_data$study)


studies_all <- unique(ps_object_genus_raw@sam_data$study)
print("genus data nolifestyle samples")
tic()
if(nested_cv_genus_no_ls_step){
  genus_no_ls_nested_cv_preds <- studies_all %>% 
    future_map_dfr(~ get_predictions_nested_cv(. ,ps = ps_object_genus_raw,
                                               prefix = "genus_data_no_ls",
                                               transform = "compositional",
                                               extra_cols = c("Observed", "Shannon")))
  genus_no_ls_nested_cv_preds_long <- genus_no_ls_nested_cv_preds %>%
    pivot_longer(cols = c("rf1", "lasso"), names_to = "model_name", values_to = "pred" )
  ps_nested_cv_genus_no_ls <- modelplots(all_results = genus_no_ls_nested_cv_preds_long)
  save(genus_no_ls_nested_cv_preds, ps_nested_cv_genus_no_ls,
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_no_ls.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_no_ls.RData")
}
toc()


if(final_model_genus_no_ls_step){
  final_model_genus_no_ls <- get_final_model(ps = ps_object_genus_raw, extra_cols = c("Observed", "Shannon"))
  genus_train_data <- create_caret_df(ps_object = ps_object_genus_raw, transformation = "compositional",
                                      mean_ab_cutoff = 5e-5, study_prevalence_cutoff = 2,
                                      prevalence_in_study_cutoff = 5,
                                      additional_cols = c("Observed", "Shannon"),
                                      only_multiple_samples = F)
  shap_out_genus_no_ls <- get_shap_long(final_model_genus_no_ls$rf1, 
                                  features = "important",
                                  test_data = genus_train_data$features)
  shap_out_genus_no_ls %>%
    filter(abs(shap_value) > 0.01,
           ab_value > 5e-5) %>%
    ggplot(., aes(y = taxon, x = shap_value, color = ab_value)) +
    geom_jitter(size = 0.5) +
    geom_violin() +
    geom_vline(xintercept = 0) +
    theme_classic()
  ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/genus_no_ls_shap.pdf")
  save(final_model_genus_no_ls, shap_out_genus_no_ls, genus_train_data, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_no_ls_genus.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_no_ls_genus.RData")
}

ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/genus_nested_cv_point_no_ls.pdf",
       plot = ps_nested_cv_genus_no_ls$correlation)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/genus_nested_cv_heatmap_no_ls.pdf",
       plot = ps_nested_cv_genus_no_ls$model_heatmap)


(age_prediction_dotplot <- genus_no_ls_nested_cv_preds %>% 
  select(age, rf1, lifestyle, study) %>%
  # filter(study != "raman_2019") %>%
  ggplot(.,aes(x = age, y = rf1, color = lifestyle)) +
  geom_point(alpha = 0.2, size = 0.7) +
  theme_classic() +
  stat_cor(# aes(colour = NULL),
    method = "spearman",
    cor.coef.name = "spearman",
    size = 6,
    show.legend = F) +
  xlab("True age [days]") +
  ylab("Predicted age [days]") +
  # geom_smooth(size = 2) +
  geom_smooth(method = "loess") +
  theme(text = element_text(size = 20),
        legend.position=c(0.8, 0.15),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 13))   )
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/genus_nested_cv_no_ls_point_present.pdf", 
       width = 10, plot = age_prediction_dotplot)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/genus_nested_cv_no_ls_point_present.png",
       width = 10, plot = age_prediction_dotplot)
save(genus_no_ls_nested_cv_preds, age_prediction_dotplot,
     file = "/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/genus_nested_cv_no_ls.RData")

cor(genus_no_ls_nested_cv_preds$age, genus_no_ls_nested_cv_preds$rf1, method = "spearman")
genus_no_ls_nested_cv_preds <- genus_no_ls_nested_cv_preds %>%
  mutate(lifestyle_industrialized = ifelse(lifestyle == "industrialized", yes = 1, no = 0),
         study = as.factor(study))
full_model <- gamlss(rf1 ~ pb(age) * lifestyle_industrialized + random(study),
                     data = genus_no_ls_nested_cv_preds %>% select(age, lifestyle_industrialized, study, rf1),
                     trace = F)
# is the path different between lifestyles different
traj_model <- gamlss(rf1 ~ pb(age) + lifestyle_industrialized + random(study),
                     data = genus_no_ls_nested_cv_preds %>% select(age, lifestyle_industrialized, study, rf1),
                     trace = F)
LR.test(traj_model, full_model)
# is the intercept different?
intercept_model <- gamlss(rf1 ~ pb(age) + pb(age):lifestyle_industrialized + random(study),
                          data = genus_no_ls_nested_cv_preds %>% select(age, lifestyle_industrialized, study, rf1),
                          trace = F)
LR.test(intercept_model, full_model)

# do the same with linear models and rank transformed predictions
full_model_lin <- lmer(rank(rf1) ~ age + lifestyle_industrialized + (1|study),
                     data = genus_no_ls_nested_cv_preds %>% select(age, lifestyle_industrialized, study, rf1))
# is the intercept different?
intercept_model_lin <- lmer(rank(rf1) ~ age + (1|study),
                          data = genus_no_ls_nested_cv_preds %>% select(age, lifestyle_industrialized, study, rf1))
test <- lrtest(intercept_model_lin, full_model_lin)


# check effect of studies:
studies <- genus_no_ls_nested_cv_preds$study %>% levels
names(studies) <- studies
get_gamlss_res <- function(dt, st) {
  print(st)
  dt = dt %>% filter(study != st) %>% select(rf1, age, lifestyle_industrialized, study)
  # f_model <- gamlss(rf1 ~ pb(age) + lifestyle_industrialized + random(study),
  #                   data = dt,
  #                   trace = F)
  # non_lin_age_model <- gamlss(rf1 ~ age + lifestyle_industrialized + random(study),
  #                             data = dt,
  #                             trace = F)
  # ls_model <- gamlss(rf1 ~ pb(age)  + random(study),
  #                    data = dt,
  #                    trace = F)
  # age_model <- gamlss(rf1 ~ lifestyle_industrialized  + random(study),
  #                     data = dt,
  #                     trace = F)
  # lr_test_res_age <- LR.test(age_model, f_model, print = F)
  # lr_test_res_age_non_lin <- LR.test(non_lin_age_model, f_model, print = F)
  f_model <- lmer(rank(rf1) ~ age + lifestyle_industrialized + (1|study),
                    data = dt)
  ls_model <- lmer(rank(rf1) ~ age  + (1|study),
                     data = dt)
  age_model <- lmer(rank(rf1) ~ lifestyle_industrialized  + (1|study),
                      data = dt)
  lr_test_res_age <- lrtest(age_model, f_model)
  lr_test_res_ls <- tryCatch(
    lrtest(ls_model, f_model),
    error = function(e){
      print(e)
      list(p.val = NA, df = NA)
    })
  linear_fm <- lm(rank(rf1) ~ age + lifestyle_industrialized + study,
                  data = dt)
  linear_intm <- lm(rank(rf1) ~ age + study,
                    data = dt)
  linear_res <- lrtest(linear_intm, linear_fm)
  e2 <- anova(linear_fm) %>% rstatix::eta_squared()
  return(list(intercetp_lmer = data.frame(summary(full_model_lin)$coefficients)$Estimate[1],
              exp_var_ls = 1 - deviance(f_model) / deviance(ls_model),
              exp_var_age = 1 - deviance(f_model) / deviance(age_model),
              coeff_ls = summary(f_model)$coefficients[3, 1],
              coeff_age = summary(f_model)$coefficients[2, 1],
              p_val_gam_ls = lr_test_res_ls$`Pr(>Chisq)`[2], 
              p_val_gam_age = lr_test_res_age$`Pr(>Chisq)`[2], 
              df_gam_ls = lr_test_res_ls$Df[2], 
              df_gam_age = lr_test_res_age$Df[2], 
              p_val_lin = linear_res$`Pr(>Chisq)`[2], 
              e2_lin = e2[2],
              exp_var_ls_lin = 1 - deviance(linear_fm) / deviance(linear_intm),
              intercept_lin = summary(linear_fm)$coefficients[3,1]))
}
# rest_test <- mclapply(c(complete = "complete", studies), 
#                       function(x) get_gamlss_res(dt = genus_no_ls_nested_cv_preds, st = x))
# rest_test_df <- do.call(rbind, lapply(rest_test, as.data.frame)) %>%
#   as.data.frame()





# SHAP analysis:
genus_no_ls_nested_cv_preds %>%
  select(run_accession, lifestyle, study, rf1, age) %>%
  left_join(shap_out_genus_no_ls,., by = "run_accession") %>%
  mutate(lifestyle = gsub("transitional", "non_industrialized", lifestyle)) %>%
  filter(!is.na(lifestyle)) %>%
  filter(abs(shap_value) > 0.01,
         ab_value > 5e-5) %>% 
  ggplot(., aes(y = taxon, x = shap_value, color = lifestyle)) +
  # geom_jitter(size = 0.5) +
  geom_violin() +
  geom_vline(xintercept = 0) +
  theme_classic()


library(rstatix)
library(ggprism)
df_p_val_shap <- genus_no_ls_nested_cv_preds %>%
  select(run_accession, lifestyle, study, rf1, age) %>%
  left_join(shap_out_genus_no_ls,., by = "run_accession") %>%
  mutate(lifestyle = gsub("transitional", "non_industrialized", lifestyle)) %>%
  filter(!is.na(lifestyle)) %>%
  filter(ab_value > 5e-5) %>% 
  rstatix::group_by(taxon) %>%
  rstatix::wilcox_test(shap_value ~ lifestyle) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj") %>%
  rstatix::add_xy_position(x = "taxon", dodge = 0.8)  # important for positioning!


genus_no_ls_nested_cv_preds %>%
  select(run_accession, lifestyle, study, rf1, age) %>%
  left_join(shap_out_genus_no_ls,., by = "run_accession") %>%
  mutate(lifestyle = gsub("transitional", "non_industrialized", lifestyle)) %>%
  filter(!is.na(lifestyle)) %>%
  filter(ab_value > 5e-5) %>% 
  ggplot(., aes(x = taxon, y = shap_value))+
  geom_violin(aes(fill = lifestyle)) +
  geom_boxplot(aes(fill = lifestyle), alpha = 0.5) +
  geom_hline(yintercept = 0)+
  add_pvalue(df_p_val_shap,
             label = "{p.adj.signif}",
             # step.increase = 0.025,
             # tip.length = 0.01,
             xmin = "xmin",
             xmax = "xmax",
             show.legend = FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_lifestyle.pdf")
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_lifestyle.png")


# MAZ and shap:
genus_no_ls_nested_cv_preds <- genus_no_ls_nested_cv_preds %>%
  mutate(age_group = group_age(age, breaks=seq(0, 735, 7)))

genus_no_ls_maz_shap <- genus_no_ls_nested_cv_preds %>%
  group_by(age_group) %>%
  summarize(group_median = median(rf1),
            group_sd = sd(rf1)) %>%
  left_join(., genus_no_ls_nested_cv_preds, by = "age_group") %>%
  mutate(MAZ = (rf1 - group_median) / group_sd) %>%
  left_join(shap_out_genus_no_ls, ., by = "run_accession")


genus_no_ls_maz_shap %>%
  filter(ab_value > 5e-5, !is.na(MAZ)) %>% 
  ggplot(., aes(y = taxon, x = shap_value, color = MAZ)) +
  # geom_jitter(size = 0.5) +
  # geom_violin() +
  geom_jitter(alpha = 0.3, size = 0.1) +
  geom_vline(xintercept = 0) +
  facet_wrap(~lifestyle)+
  scale_color_gradientn(name = "MAZ", colors = c("red", "green")) +
  theme_classic()


genus_no_ls_maz_shap %>%
  filter(ab_value > 5e-5, !is.na(MAZ)) %>% 
  filter(taxon %in% c("Prevotella",
                      "Faecalibacterium",
                      "Ruminococcus",
                      "Subdoligranulum")) %>%
  ggplot(., aes(x = shap_value, y = MAZ, color = taxon)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_smooth(method = "lm") +
  # scale_color_gradientn(name = "MAZ", colors = c("red", "green")) +
  theme_classic()
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_maz_taxon.pdf")
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_maz_taxon.png")



# check feature importances:
create_feature_importance_df <- function(study) {
  train_object <- readRDS(paste0("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/nested_cv_dataset_models/",
                                 "genus_data_no_ls_", study, ".rds"))
  var_imp <- train_object$rf1$finalModel$variable.importance %>%
    sort(., decreasing = T) %>%
    head(15)
  var_imp <- (var_imp / max(var_imp)) *100
  # Extract variable importance data
  data.frame(importance = var_imp,
             taxon = names(var_imp),
             ID = study)
}
combined_importance_genus_no_ls <- map(studies_all, create_feature_importance_df) %>%
  purrr::reduce(rbind)  # Combine the individual plots
combined_importance_genus_no_ls %>%
  ggplot(., aes(x = taxon, y = importance)) +
  geom_boxplot() +  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
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


saveRDS(combined_importance_genus_no_ls,
        "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/genus_nested_cv_important_taxa_no_ls.rds")



# genus data for industrialized populaitons only specific model ##################

ps_object_genus_raw_industrialized <- ps_object_genus_raw %>%
  subset_samples(lifestyle == "industrialized")
industrialized_studies <- unique(ps_object_genus_raw_industrialized@sam_data$study)

print("genus data only for industrialized samples")
tic()
if(nested_cv_genus_industrialized_step){
  genus_industrialized_nested_cv_preds <- industrialized_studies %>% 
    future_map_dfr(~ get_predictions_nested_cv(. ,ps = ps_object_genus_raw_industrialized,
                                               prefix = "genus_data_industrialized",
                                               transform = "compositional",
                                               extra_cols = c("Observed", "Shannon")))
  genus_industrialized_nested_cv_preds_long <- genus_industrialized_nested_cv_preds %>%
    pivot_longer(cols = c("rf1", "lasso"), names_to = "model_name", values_to = "pred" )
  ps_nested_cv_genus_industrialized <- modelplots(all_results = genus_industrialized_nested_cv_preds_long)
  save(genus_industrialized_nested_cv_preds, ps_nested_cv_genus_industrialized,
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_industrialized.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_industrialized.RData")
}
toc()


if(final_model_genus_industrialized_step){
  final_model_genus_industrialized <- get_final_model(ps = ps_object_genus_raw_industrialized, extra_cols = c("Observed", "Shannon"))
  genus_train_data_industrialized <- create_caret_df(ps_object = ps_object_genus_raw_industrialized, transformation = "compositional",
                                      mean_ab_cutoff = 5e-5, study_prevalence_cutoff = 2,
                                      prevalence_in_study_cutoff = 5,
                                      additional_cols = c("Observed", "Shannon"),
                                      only_multiple_samples = F)
  shap_out_genus_industrialized <- get_shap_long(final_model_genus_industrialized$rf1, 
                                        features = "important",
                                        test_data = genus_train_data_industrialized$features)
  shap_out_genus_industrialized %>%
    filter(abs(shap_value) > 0.01,
           ab_value > 5e-5) %>%
    ggplot(., aes(y = taxon, x = shap_value, color = ab_value)) +
    geom_jitter(size = 0.5) +
    geom_violin() +
    geom_vline(xintercept = 0) +
    theme_classic()
  ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/genus_industrialized_shap.pdf")
  
  save(final_model_genus_industrialized, shap_out_genus_industrialized, genus_train_data_industrialized, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_industrialized_genus.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_industrialized_genus.RData")
}


ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/genus_nested_cv_point_industrialized.pdf",
       plot = ps_nested_cv_genus_industrialized$correlation)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/genus_nested_cv_heatmap_industrialized.pdf",
       plot = ps_nested_cv_genus_industrialized$model_heatmap)

# check feature importances:
create_feature_importance_df_genus <- function(study) {
  train_object <- readRDS(paste0("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/nested_cv_dataset_models/",
                                 "genus_data_industrialized_", study, ".rds"))
  var_imp <- train_object$rf1$finalModel$variable.importance %>%
    sort(., decreasing = T) %>%
    head(15)
  var_imp <- (var_imp / max(var_imp)) *100
  # Extract variable importance data
  data.frame(importance = var_imp,
             taxon = names(var_imp),
             ID = study)
}
combined_importance_genus_industrialized <- map(industrialized_studies, create_feature_importance_df_genus) %>%
  purrr::reduce(rbind)  %>%
  mutate(taxon = gsub("X.Eubacterium..coprostanoligenes.group", "E.coprostanoligenes_group", taxon),
         lifestyle = "industrialized")

ggplot(combined_importance_genus_industrialized, aes(x = reorder(taxon, importance), y = importance)) +
  # geom_line() + 
  geom_boxplot() +
  ylim(0,100) +
  ylab("Relative importance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# genus data for non industrialized populations only #############################

ps_object_genus_raw_nonindustrialized <- ps_object_genus_raw %>%
  subset_samples(lifestyle != "industrialized")
nonindustrialized_studies <- unique(ps_object_genus_raw_nonindustrialized@sam_data$study)


print("genus data only for non-industrialized samples")
tic()
if(nested_cv_genus_nonindustrialized_step){
  genus_nonindustrialized_nested_cv_preds <- nonindustrialized_studies %>% 
    future_map_dfr(~ get_predictions_nested_cv(. ,ps = ps_object_genus_raw_nonindustrialized,
                                               prefix = "genus_data_nonindustrialized",
                                               transform = "compositional",
                                               extra_cols = c("Observed", "Shannon")))
  genus_nonindustrialized_nested_cv_preds_long <- genus_nonindustrialized_nested_cv_preds %>%
    pivot_longer(cols = c("rf1", "lasso"), names_to = "model_name", values_to = "pred" )
  ps_nested_cv_genus_nonindustrialized <- modelplots(all_results = genus_nonindustrialized_nested_cv_preds_long)
  save(genus_nonindustrialized_nested_cv_preds, ps_nested_cv_genus_nonindustrialized,
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_nonindustrialized.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_nonindustrialized.RData")
}
toc()

if(final_model_genus_nonindustrialized_step){
  final_model_genus_nonindustrialized <- get_final_model(ps = ps_object_genus_raw_nonindustrialized, extra_cols = c("Observed", "Shannon"))
  genus_train_data_non_industrialized <- create_caret_df(ps_object = ps_object_genus_raw_nonindustrialized, transformation = "compositional",
                                              mean_ab_cutoff = 5e-5, study_prevalence_cutoff = 2,
                                              prevalence_in_study_cutoff = 5,
                                              additional_cols = c("Observed", "Shannon"),
                                              only_multiple_samples = F)
  shap_out_genus_nonindustrialized <- get_shap_long(final_model_genus_nonindustrialized$rf1, 
                                          features = "important",
                                          test_data = genus_train_data_non_industrialized$features)
  shap_out_genus_nonindustrialized %>%
    filter(abs(shap_value) > 0.01,
           ab_value > 5e-5) %>%
    ggplot(., aes(y = taxon, x = shap_value, color = ab_value)) +
    geom_jitter(size = 0.5) +
    geom_violin() +
    geom_vline(xintercept = 0) +
    theme_classic()
  ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/genus_nonindustrialized_shap.pdf")
  
  save(final_model_genus_nonindustrialized, shap_out_genus_nonindustrialized, genus_train_data_non_industrialized,
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_nonindustrialized_genus.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_nonindustrialized_genus.RData")
}

ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/genus_nested_cv_point_nowest.pdf",
       plot = ps_nested_cv_genus_nonindustrialized$correlation)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/genus_nested_cv_heatmap_nowest.pdf",
       plot = ps_nested_cv_genus_nonindustrialized$model_heatmap)

# check feature importances:
create_feature_importance_df_nw_genus <- function(study) {
  train_object <- readRDS(paste0("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/nested_cv_dataset_models/",
                                 "genus_data_nonindustrialized_", study, ".rds"))
  var_imp <- train_object$rf1$finalModel$variable.importance %>%
    sort(., decreasing = T) %>%
    head(15)
  var_imp <- (var_imp / max(var_imp)) *100
  # Extract variable importance data
  data.frame(importance = var_imp,
             taxon = names(var_imp),
             ID = study)
}
combined_importance_nonindustrialized_genus <- map(nonindustrialized_studies, create_feature_importance_df_nw_genus) %>%
  purrr::reduce(rbind) %>%
  mutate(taxon = gsub("X.Eubacterium..coprostanoligenes.group", "E.coprostanoligenes_group", taxon),
         lifestyle = "non_industrialized")

combined_importance_lifestyles_genus <- rbind(combined_importance_genus_industrialized, combined_importance_nonindustrialized_genus) %>%
  mutate(lifestyle = factor(lifestyle, levels = c("industrialized", "non_industrialized")),
         measure_type = ifelse(taxon %in% c("Shannon", "Observed"), yes = "a_div", no = "taxon"),
         # in_lit = ifelse(taxon %in% in_literature, yes = "red", no = "black"),
         taxon = reorder(taxon, importance),
         taxon = factor(taxon, levels = unique(c("Shannon", "Observed", levels(taxon)))))

saveRDS(combined_importance_lifestyles_genus,
        "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/genus_nested_cv_important_taxa_lifestyle.rds")

combined_importance_lifestyles_genus %>% 
  mutate(taxon = gsub("Bacteria_Firmicutes_Clostridia_unclassified_Clostridia_unclassified_Clostridia_unclassified_Clostridia",
                      "unclassified_Clostridia", taxon)) %>% 
  ggplot(., aes(x = taxon, y = importance)) +
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
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/genus_nested_cv_importance_lifestyle.pdf")
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/genus_nested_cv_importance_lifestyle.png")




# genus data without lifestyle rarefied #################################################

studies_all <- unique(ps_object_genus_raw@sam_data$study)
ps_object_genus_rare_h <- ps_object_genus_raw %>% rarefy_even_depth(., sample.size = 2000)
print("genus data nolifestyle samples rare_h")
tic()
if(nested_cv_genus_no_ls_rarefied_step){
  genus_no_ls_nested_rare_h_cv_preds <- studies_all %>% 
    future_map_dfr(~ get_predictions_nested_cv(. ,ps = ps_object_genus_rare_h,
                                               prefix = "genus_data_no_ls_rare",
                                               transform = "compositional",
                                               extra_cols = c("Observed", "Shannon")))
  save(genus_no_ls_nested_rare_h_cv_preds,
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_no_ls_rare_h.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_no_ls_rare_h.RData")
}
toc()


if(final_model_genus_no_ls_rarefied_step){
  final_model_genus_no_ls_rare_h <- get_final_model(ps = ps_object_genus_rare_h, extra_cols = c("Observed", "Shannon"))
  genus_train_data_rare_h <- create_caret_df(ps_object = ps_object_genus_rare_h, transformation = "compositional",
                                      mean_ab_cutoff = 5e-5, study_prevalence_cutoff = 2,
                                      prevalence_in_study_cutoff = 5,
                                      additional_cols = c("Observed", "Shannon"),
                                      only_multiple_samples = F)
  shap_out_genus_no_ls_rare_h <- get_shap_long(final_model_genus_no_ls_rare_h$rf1, 
                                        features = "important",
                                        test_data = genus_train_data_rare_h$features)
  save(final_model_genus_no_ls_rare_h, shap_out_genus_no_ls_rare_h, genus_train_data_rare_h, 
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_no_ls_genus_rare_h.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_no_ls_genus_rare_h.RData")
}

# genus data for industrialized populations only specific model ##################

ps_object_genus_rare_h_industrialized <- ps_object_genus_rare_h %>%
  subset_samples(lifestyle == "industrialized")
industrialized_studies_rare_h <- unique(ps_object_genus_rare_h_industrialized@sam_data$study)

print("genus data only for industrialized samples")
tic()
if(nested_cv_genus_industrialized_rarefied_step){
  genus_industrialized_nested_cv_rare_h_preds <- industrialized_studies_rare_h %>% 
    future_map_dfr(~ get_predictions_nested_cv(. ,ps = ps_object_genus_rare_h_industrialized,
                                               prefix = "genus_data_industrialized_rare_h",
                                               transform = "compositional",
                                               extra_cols = c("Observed", "Shannon")))
  save(genus_industrialized_nested_cv_rare_h_preds,
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_industrialized_rare_h.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_industrialized_rare_h.RData")
}
toc()


if(final_model_genus_industrialized_rarefied_step){
  final_model_genus_industrialized_rare_h <- get_final_model(ps = ps_object_genus_rare_h_industrialized, extra_cols = c("Observed", "Shannon"))
  genus_train_data_industrialized_rare_h <- create_caret_df(ps_object = ps_object_genus_rare_h_industrialized, transformation = "compositional",
                                                     mean_ab_cutoff = 5e-5, study_prevalence_cutoff = 2,
                                                     prevalence_in_study_cutoff = 5,
                                                     additional_cols = c("Observed", "Shannon"),
                                                     only_multiple_samples = F)
  shap_out_genus_industrialized_rare_h <- get_shap_long(final_model_genus_industrialized_rare_h$rf1, 
                                                 features = "important",
                                                 test_data = genus_train_data_industrialized_rare_h$features)
  save(final_model_genus_industrialized_rare_h, shap_out_genus_industrialized_rare_h, genus_train_data_industrialized_rare_h,
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_industrialized_genus_rare_h.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_industrialized_genus_rare_h.RData")
}

# genus data for non industrialized populations only #############################

ps_object_genus_rare_h_nonindustrialized <- ps_object_genus_rare_h %>%
  subset_samples(lifestyle != "industrialized")
nonindustrialized_studies_rare_h <- unique(ps_object_genus_rare_h_nonindustrialized@sam_data$study)


print("genus data only for non-industrialized samples _rare_h")
tic()
if(nested_cv_genus_nonindustrialized_rarefied_step){
  genus_nonindustrialized_nested_cv_rare_h_preds <- nonindustrialized_studies_rare_h %>% 
    future_map_dfr(~ get_predictions_nested_cv(. ,ps = ps_object_genus_rare_h_nonindustrialized,
                                               prefix = "genus_data_nonindustrialized_rare_h",
                                               transform = "compositional",
                                               extra_cols = c("Observed", "Shannon")))
    save(genus_nonindustrialized_nested_cv_rare_h_preds,
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_nonindustrialized_rare_h.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_nonindustrialized_rare_h.RData")
}
toc()

if(final_model_genus_nonindustrialized_rarefied_step){
  final_model_genus_nonindustrialized_rare_h <- get_final_model(ps = ps_object_genus_rare_h_nonindustrialized, extra_cols = c("Observed", "Shannon"))
  genus_train_data_non_industrialized_rare_h <- create_caret_df(ps_object = ps_object_genus_rare_h_nonindustrialized, transformation = "compositional",
                                                         mean_ab_cutoff = 5e-5, study_prevalence_cutoff = 2,
                                                         prevalence_in_study_cutoff = 5,
                                                         additional_cols = c("Observed", "Shannon"),
                                                         only_multiple_samples = F)
  shap_out_genus_nonindustrialized_rare_h <- get_shap_long(final_model_genus_nonindustrialized_rare_h$rf1, 
                                                    features = "important",
                                                    test_data = genus_train_data_non_industrialized_rare_h$features)
  save(final_model_genus_nonindustrialized_rare_h, shap_out_genus_nonindustrialized_rare_h, genus_train_data_non_industrialized_rare_h,
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_nonindustrialized_genus_rare_h.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_nonindustrialized_genus_rare_h.RData")
}





# genus data without lifestyle rarefied low #################################################

studies_all <- unique(ps_object_genus_raw@sam_data$study)
ps_object_genus_rare_l <- ps_object_genus_raw %>% rarefy_even_depth(.) %>%
  filter_taxa(function(x) mean(x) > 0, TRUE)  # mean abundance cutoff, only counting abundances > 0
  
print("genus data nolifestyle samples rare_l")
tic()
if(nested_cv_genus_no_ls_rarefied_l_step){
  genus_no_ls_nested_rare_l_cv_preds <- studies_all %>% 
    map_dfr(~ get_predictions_nested_cv(. ,ps = ps_object_genus_rare_l,
                                               prefix = "genus_data_no_ls_rare_l",
                                               transform = "compositional",
                                        rf_depth = min(sample_sums(ps_object_genus_rare_l)),
                                               extra_cols = c("Observed", "Shannon")))
  save(genus_no_ls_nested_rare_l_cv_preds,
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_no_ls_rare_l.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_no_ls_rare_l.RData")
}
toc()


if(final_model_genus_no_ls_rarefied_l_step){
  final_model_genus_no_ls_rare_l <- get_final_model(ps = ps_object_genus_rare_l, extra_cols = c("Observed", "Shannon"),
                                                    rf_depth = min(sample_sums(ps_object_genus_rare_l)))
  genus_train_data_rare_l <- create_caret_df(ps_object = ps_object_genus_rare_l, transformation = "compositional",
                                             rf_depth = min(sample_sums(ps_object_genus_rare_l)),
                                             mean_ab_cutoff = 5e-5, study_prevalence_cutoff = 2,
                                             prevalence_in_study_cutoff = 5,
                                             additional_cols = c("Observed", "Shannon"),
                                             only_multiple_samples = F)
  shap_out_genus_no_ls_rare_l <- get_shap_long(final_model_genus_no_ls_rare_l$rf1, 
                                               features = "important",
                                               test_data = genus_train_data_rare_l$features)
  save(final_model_genus_no_ls_rare_l, shap_out_genus_no_ls_rare_l, genus_train_data_rare_l, 
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_no_ls_genus_rare_l.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_no_ls_genus_rare_l.RData")
}

# genus data for industrialized populations only specific model ##################

ps_object_genus_rare_l_industrialized <- ps_object_genus_rare_l %>%
  subset_samples(lifestyle == "industrialized")
industrialized_studies_rare_l <- unique(ps_object_genus_rare_l_industrialized@sam_data$study)

print("genus data only for industrialized samples")
tic()
if(nested_cv_genus_industrialized_rarefied_l_step){
  genus_industrialized_nested_cv_rare_l_preds <- industrialized_studies_rare_l %>% 
    future_map_dfr(~ get_predictions_nested_cv(. ,ps = ps_object_genus_rare_l_industrialized,
                                               prefix = "genus_data_industrialized_rare_l",
                                               transform = "compositional",
                                               rf_depth = min(sample_sums(ps_object_genus_rare_l)),
                                               extra_cols = c("Observed", "Shannon")))
  save(genus_industrialized_nested_cv_rare_l_preds,
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_industrialized_rare_l.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_industrialized_rare_l.RData")
}
toc()


if(final_model_genus_industrialized_rarefied_l_step){
  final_model_genus_industrialized_rare_l <- get_final_model(ps = ps_object_genus_rare_l_industrialized,
                                                             extra_cols = c("Observed", "Shannon"),
                                                             rf_depth = min(sample_sums(ps_object_genus_rare_l)))
  genus_train_data_industrialized_rare_l <- create_caret_df(ps_object = ps_object_genus_rare_l_industrialized, transformation = "compositional",
                                                            mean_ab_cutoff = 5e-5, study_prevalence_cutoff = 2,
                                                            prevalence_in_study_cutoff = 5,
                                                            rf_depth = min(sample_sums(ps_object_genus_rare_l)),
                                                            additional_cols = c("Observed", "Shannon"),
                                                            only_multiple_samples = F)
  shap_out_genus_industrialized_rare_l <- get_shap_long(final_model_genus_industrialized_rare_l$rf1, 
                                                        features = "important",
                                                        test_data = genus_train_data_industrialized_rare_h$features)
  save(final_model_genus_industrialized_rare_l, shap_out_genus_industrialized_rare_l, genus_train_data_industrialized_rare_l,
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_industrialized_genus_rare_l.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_industrialized_genus_rare_l.RData")
}

# genus data for non industrialized populations only #############################

ps_object_genus_rare_l_nonindustrialized <- ps_object_genus_rare_l %>%
  subset_samples(lifestyle != "industrialized")
nonindustrialized_studies_rare_l <- unique(ps_object_genus_rare_l_nonindustrialized@sam_data$study)


print("genus data only for non-industrialized samples _rare_l")
tic()
if(nested_cv_genus_nonindustrialized_rarefied_l_step){
  genus_nonindustrialized_nested_cv_rare_l_preds <- nonindustrialized_studies_rare_l %>% 
    future_map_dfr(~ get_predictions_nested_cv(. ,ps = ps_object_genus_rare_l_nonindustrialized,
                                               prefix = "genus_data_nonindustrialized_rare_l",
                                               transform = "compositional",
                                               rf_depth = min(sample_sums(ps_object_genus_rare_l)),
                                               extra_cols = c("Observed", "Shannon")))
  save(genus_nonindustrialized_nested_cv_rare_l_preds,
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_nonindustrialized_rare_l.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_nonindustrialized_rare_l.RData")
}
toc()

if(final_model_genus_nonindustrialized_rarefied_l_step){
  final_model_genus_nonindustrialized_rare_l <- get_final_model(ps = ps_object_genus_rare_l_nonindustrialized,
                                                                extra_cols = c("Observed", "Shannon"),
                                                                rf_depth = min(sample_sums(ps_object_genus_rare_l)))
  genus_train_data_non_industrialized_rare_l <- create_caret_df(ps_object = ps_object_genus_rare_l_nonindustrialized, transformation = "compositional",
                                                                mean_ab_cutoff = 5e-5, study_prevalence_cutoff = 2,
                                                                prevalence_in_study_cutoff = 5,
                                                                rf_depth = min(sample_sums(ps_object_genus_rare_l)),
                                                                additional_cols = c("Observed", "Shannon"),
                                                                only_multiple_samples = F)
  shap_out_genus_nonindustrialized_rare_l <- get_shap_long(final_model_genus_nonindustrialized_rare_l$rf1, 
                                                           features = "important",
                                                           test_data = genus_train_data_non_industrialized_rare_l$features)
  save(final_model_genus_nonindustrialized_rare_l, shap_out_genus_nonindustrialized_rare_l, 
       genus_train_data_non_industrialized_rare_l,
       file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_nonindustrialized_genus_rare_l.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_nonindustrialized_genus_rare_l.RData")
}



# compare runs #################################################################

family_nested_cv_preds$variation <- "two_years"
family_no_ls_nested_cv_preds$variation <- "family_no_ls"
family_industrialized_nested_cv_preds$variation <- "industrialized"
family_industrialized_nested_cv_preds$lifestyle_industrialized <- 1
genus_nested_cv_preds$variation <- "genus"
family_nonindustrialized_nested_cv_preds$variation <- "non_industrialized"
family_nonindustrialized_nested_cv_preds$lifestyle_industrialized <- 0
picrust_nested_cv_preds$variation <- "picrust"
picrust_industrialized_nested_cv_preds$variation <- "picrust_industrialized"
picrust_pathway_nested_cv_preds$variation <- "pathways"
picrust_gbm_nested_cv_preds$variation <- "gbms"
picrust_gmm_nested_cv_preds$variation <- "gmms"


all_models <- bind_rows(family_nested_cv_preds,
                    family_no_ls_nested_cv_preds,
                    family_industrialized_nested_cv_preds,
                    genus_nested_cv_preds,
                    family_nonindustrialized_nested_cv_preds,
                    picrust_nested_cv_preds,
                    picrust_industrialized_nested_cv_preds,
                    picrust_gbm_nested_cv_preds,
                    picrust_pathway_nested_cv_preds,
                    picrust_gmm_nested_cv_preds)

all_models_long <- all_models %>%
  pivot_longer(cols = c("rf1", "lasso"), names_to = "model_name", values_to = "pred" )


all_lm_list <- get_lm_list(pred_df = all_models_long, grouping = c("study", "variation")) %>%
  mutate(variation = factor(variation,
                            levels = c("two_years", "binary", "months", "pathways",
                                       "gbms", "genus", "no_ab", "industrialized", "one_year", "non_industrialized",
                                       "family_no_ls", "picrust", "picrust_industrialized", "gmms")),
         model_name = factor(model_name, levels = c("rf1", "lasso")))

all_studies <- set_list$all
studies_here <- unique(all_lm_list$study) %in% all_studies

ggplot(all_lm_list, aes(x = variation, y = R2))+
  geom_boxplot() +
  geom_line(aes(group=study, color = study), size=0.4, alpha=0.5) +
  geom_point(aes(group=study, color = study), size=1, alpha=0.9) +
  ylim(0, 1) +
  scale_color_manual(values = turbo(length(all_studies))[studies_here], breaks = all_studies[studies_here] )+
  facet_wrap(~model_name) + 
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1),)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/variations_x_models.pdf")


# influence of binary data on different models:
binary_lm_list <- all_lm_list %>% 
  filter(variation %in% c("binary", "two_years")) %>%
  mutate(variation = gsub("two_years", "compositional", variation))
         # model_name = factor(model_name),
         # model_name_variation = paste(model_name, variation, sep = "_"))

res_wilc <- pairwise.wilcox.test(binary_lm_list$R2, binary_lm_list$model_name_variation)

custom_colors <- c("rf1" = "#0000ff", "lasso" = "#00ff00")

library(rstatix)
library(ggprism)

df_p_val1 <- binary_lm_list %>%
  arrange(study) %>%
  rstatix::group_by(variation) %>%
  rstatix::wilcox_test(R2 ~ model_name, ref.group = "rf1", alternative = "greater", paired = T) %>%
  filter(variation == "compositional") %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj") %>% 
  rstatix::add_xy_position() %>%
  mutate(xmin = xmin + 0.2,
          xmax = xmax + 0.2)

df_p_val2 <- binary_lm_list %>%
  arrange(study) %>%
  rstatix::group_by(model_name) %>%
  rstatix::wilcox_test(R2 ~ variation, ref.group = "compositional", paired = T) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj") %>%
  rstatix::add_xy_position(x = "model_name", dodge = 0.8)  # important for positioning!

ggplot(binary_lm_list, aes(x = model_name, y = R2))+
  geom_boxplot(aes(fill = variation)) +
  ylim(0, 1.05) +
  add_pvalue(df_p_val1,
             label = "{p.adj.signif}",
             step.increase = 0.025,
             tip.length = 0.01,
             xmin = "xmin",
             xmax = "xmax",
             show.legend = FALSE) +
  add_pvalue(df_p_val2, 
             xmin = "xmin", 
             xmax = "xmax",
             bracket.nudge.y = -0.76,
             label = "{p.adj.signif}",
             tip.length = -0.01) +
  theme_minimal() +
  theme(axis.text.x = element_text(colour = custom_colors[c(1, 5, 4)]))

ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/binary_model.pdf", width = 14, height = 8)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/binary_model.png", width = 7, height = 4)



# influence of different data
variations_lm_list <- all_lm_list %>% 
  filter(variation != "binary") %>%
  filter(model_name == "rf1") %>%
  filter(!variation %in% c("months", "no_ab", "one_year", "industrialized", "non_industrialized", "picrust_industrialized", "family_no_ls")) %>%
  mutate(variation = droplevels(variation)) %>%
  mutate(variation = factor(variation, levels = c("two_years", "genus", "pathways", "picrust", "gbms", "gmms")))
# model_name = factor(model_name),
# model_name_variation = paste(model_name, variation, sep = "_"))

df_p_val_all <- variations_lm_list %>%
  arrange(study) %>%
  # filter(variation != "industrialized", variation != "non_industrialized") %>%
  mutate(variation = droplevels(variation)) %>%
  # rstatix::group_by(variation) %>%
  rstatix::wilcox_test(R2 ~ variation, ref.group = "genus") %>%
  # rstatix::wilcox_test(R2 ~ variation, ref.group = "two_years", alternative = "less", paired = T) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj") %>% 
  rstatix::add_xy_position()

df_p_val_all_2 <- variations_lm_list %>%
  filter(variation == "industrialized"| variation == "non_industrialized" | variation == "two_years") %>%
  mutate(variation = droplevels(variation)) %>% 
  rstatix::wilcox_test(R2 ~ variation, ref.group = "genus", alternative = "less") %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj") %>% 
  rstatix::add_xy_position() %>%
  mutate(xmax = xmax + 5)

ggplot(variations_lm_list, aes(x = variation, y = R2))+
  geom_boxplot(aes(fill = variation)) +
  # ylim(0, 1) +
  add_pvalue(df_p_val_all,
             label = "{p.adj.signif}",
             # step.group.by = "variation",
             step.increase = 0.045,
             tip.length = 0.01,
             # bracket.nudge.y = 0.02,
             xmin = "xmin",
             xmax = "xmax",
             show.legend = FALSE) +
  theme_classic() #+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/variations.pdf")
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/variations.png")

# correlation between predictions from different data:
taxon_functional_wide <- all_models_long %>%
  filter(model_name == "rf1",
         !variation %in% c("months", "no_ab", "family_no_ls", "one_year", "non_industrialized", "industrialized")) %>%
  pivot_wider(names_from = variation, values_from = pred, id_cols = c("sample_ID", "run_accession", "age")) 

plot(taxon_functional_wide %>% select(-sample_ID, -run_accession))
cor_plot(taxon_functional_wide %>% select(-sample_ID, -run_accession) %>% cor(.,use = "pairwise.complete.obs"), label = T)
plot(taxon_functional_wide$genus, taxon_functional_wide$picrust)
# Idea: maybe train a model on gbms and taxonomy

