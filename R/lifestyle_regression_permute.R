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
library(gridExtra)
library(tictoc)
library(gridExtra)
library(furrr)
library(rstatix)
library(ggprism)
library(permute)
source("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/alt_models.R")
source("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/regression_functions.R")
source("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/lifestyle_regression_functions.R")
source("/fast/AG_Forslund/rob/mm_index/R_scripts/setlists.R")

reshuffle_genus_step <- T


# get the list of instances downsampled to a specific size
get_list <- function(vl, dfrm, max = NA, study_prob = F) {
  lst <- lapply(c(westernized = "westernized", non_westernized = "non_westernized"),
                function(x) {dfrm %>%
                    filter(lifestyle == x) %>%
                    pull(vl) %>% unique}) # get all ids for that lifestyle
  len <- lapply(lst, length) %>% unlist %>% min(c(., max), na.rm = T)
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


# splits the dataset into test and training set
# traing set gets downsampled so that both lifestyles have the same amount
# of studies/subjects and samples
get_test_train_data <- function(test_set, ps, extra_cols = c("Observed", "Shannon"),
                                max_samples = NA, max_studies = NA, max_subj = NA) {
  # spit into test/train data
  oldDF <- as(sample_data(ps), "data.frame")
  trainDF <- subset(oldDF, study != test_set)
  testDF <- subset(oldDF, study == test_set)
  ps_train <- ps
  ps_test <- ps
  rm(ps, oldDF)
  # subsample train set to equal amount of studies/individuals/samples per lifestyle
  stds <- get_list("study", trainDF, max = max_studies)
  df_red <- trainDF %>%
    filter(lifestyle == "non_westernized" | study %in% stds$westernized) # subsample studies
  sbjcts <- get_list("subject_ID", trainDF, max = max_subj, study_prob = T)
  df_red_sb <- df_red %>%
    filter(lifestyle == "non_westernized" | subject_ID %in% sbjcts$westernized) %>%
    filter(lifestyle == "westernized" | subject_ID %in% sbjcts$non_westernized) # subsample individuals
  df_red_y <- filter(df_red_sb, age <= 365) # same distribuition for above and below one year respectively
  df_red_o <- filter(df_red_sb, age > 365)
  smpls_y <- get_list("sample_ID", df_red_y, max = max_samples / 2, study_prob = T)
  df_red_y_smp <- df_red_y %>%
    filter(lifestyle == "non_westernized" | sample_ID %in% smpls_y$westernized) %>%
    filter(lifestyle == "westernized" | sample_ID %in% smpls_y$non_westernized)
  if(length(table(df_red_o$lifestyle)) > 1) {
    smpls_o <- get_list("sample_ID", df_red_o, max = max_samples / 2, study_prob = T)
    df_red_o_smp <- df_red_o %>%
      filter(lifestyle == "non_westernized" | sample_ID %in% smpls_o$westernized) %>%
      filter(lifestyle == "westernized" | sample_ID %in% smpls_o$non_westernized) # subsample samples
    df_red_smp <- rbind(df_red_y_smp, df_red_o_smp)
  } else {
    print("No old samples found for one lifestyle")
    df_red_smp <- df_red_y_smp
  }
  sample_data(ps_train) <- sample_data(df_red_smp)
  sample_data(ps_test) <- sample_data(testDF)
  return(list(test = ps_test, train = ps_train))
}


get_predictions <- function(ps_test, ps_train, extra_cols = c("Observed", "Shannon"),
                            test_set_name = NULL, run_importances_boruta = T,
                            run_importances_shap = F){
  tic()
  print(test_set_name)
  for_caret_list_train <- create_caret_df(ps_object = ps_train, transformation = "compositional",
                                          mean_ab_cutoff = 5e-5, study_prevalence_cutoff = 1,
                                          additional_cols = extra_cols,
                                          only_multiple_samples = F)
  for_caret_list_test <- create_caret_df(ps_object = ps_test, transformation = "compositional",
                                         additional_cols = extra_cols,
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
  rfGrid <- expand.grid(mtry = c(10),
                        num.trees = c(200),
                        min.node.size = c(5),
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
  if(run_importances_boruta) {
    # get boruta importance:
    print("finished model, run Boruta")
    selected_features <- get_boruta_important_features(train_data = for_caret_list_train,
                                                       model = rf_model,
                                                       study = test_set_name)
  } else {selected_features = NULL}
  if(run_importances_shap){
    # run shap
    print("Get SHAP values")
    shap_out <- get_shap_long(rf_model,
                              test_data = for_caret_list_train$features,
                              features = unique(selected_features$importances$taxon))
  }else {shap_out = NULL}
  toc()
  return(list(preds = predictions,
              boruta = selected_features, 
              shap = shap_out, 
              test_set = test_set_name,
              training_data = for_caret_list_train$features))
}


west_non_west_preds <- function(test_set, ps, extra_cols = c("Observed", "Shannon"),
                                max_samples = NA, max_studies = NA, max_subj = NA) {
  print(test_set)
  test_train_set <- get_test_train_data(test_set, ps, extra_cols,
                                        max_samples = max_samples,
                                        max_studies = max_studies,
                                        max_subj = max_subj)
  # combined
  print("run combined model")
  combined_res <- get_predictions(ps_test = test_train_set$test, 
                                  ps_train = test_train_set$train,
                                  extra_cols = extra_cols,
                                  test_set_name = test_set,
                                  run_importances_boruta = F,
                                  run_importances_shap = F) 
  combined_preds <- combined_res$preds %>%
    select(sample_ID, subject_ID, study, lifestyle, age, pred) %>%
    mutate(pred_combined = pred, .keep = "unused")
  
  # western
  print("run western model")
  test_train_set_west <- test_train_set
  test_train_set_west$train <- test_train_set_west$train %>% subset_samples(lifestyle == "westernized")
  west_res <- get_predictions(ps_test = test_train_set_west$test, 
                              ps_train = test_train_set_west$train,
                              extra_cols = extra_cols,
                              test_set_name = test_set) 
  west_preds <- west_res$preds %>%
    select(sample_ID, subject_ID, study, lifestyle, age, pred) %>%
    mutate(pred_western = pred, .keep = "unused")
  
  # non western
  print("run non western model")
  test_train_set_non_west <- test_train_set
  test_train_set_non_west$train <- test_train_set_non_west$train %>% subset_samples(lifestyle == "non_westernized")
  non_west_res <- get_predictions(ps_test = test_train_set_non_west$test,
                                  ps_train = test_train_set_non_west$train,
                                  extra_cols = extra_cols,
                                  test_set_name = test_set) 
  non_west_preds <- non_west_res$preds %>%
    select(sample_ID, subject_ID, study, lifestyle, age, pred) %>%
    mutate(pred_non_western = pred, .keep = "unused")
  # combine
  all_preds <- full_join(west_preds, non_west_preds) %>%
    full_join(.,combined_preds)
  cat("finished with ", test_set, "\n")
  return(list(preds = all_preds, ft_imps = list(west = west_res, non_west = non_west_res)))
}


get_permute_preds <- function(ps, blocks = "subject_ID", round = NULL) {
  # create permutation
  cat("Permuting round: ", round, "\n")
  mtdt <- ps %>%
    sample_data() %>%
    data.frame 
  shuffle_data <- mtdt %>%
    select(blocks, lifestyle) %>%  # assign each individual a new lifestyle definition randomly
    distinct()
  rwnms <- rownames(mtdt)
  shuffle_data$lifestyle_shuffled <- shuffle_data$lifestyle[shuffle(nrow(shuffle_data))]
  mtdt <- left_join(mtdt, shuffle_data)
  rownames(mtdt) <- rwnms
  sample_data(ps) <- sample_data(mtdt)
  # get predictions
  print("start cv")
  preds <- unique(mtdt$study) %>%
    future_map(~ west_non_west_preds(. ,ps = ps,
                                     extra_cols = c("Observed", "Shannon")), .progress = T)
  return(preds)
  }


# n_cores <- 5

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

cl <- makePSOCKcluster(ceiling(n_cores/4))
registerDoParallel(cl)
getDoParWorkers()
set.seed(825)
future::plan(multisession, workers = ceiling(n_cores/4))

# genus data ###################################################################

ps_object_genus_raw <- readRDS("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_phyloseq_rf_filter_genus.rds") %>%
  subset_samples(age <= 730 & age > 1)
metadata <- ps_object_genus_raw %>%
  sample_data() %>%
  data.frame 

print("start reshuffling step")
tic()
if(reshuffle_genus_step){
  # permuted_preds_genus <- get_permute_preds(ps = ps_object_genus_raw, blocks = "study")
  permuted_preds_genus_list <- future_map(1:50, ~ get_permute_preds(ps = ps_object_genus_raw,
                                                                    blocks = "subject_ID",
                                                                    round = .x))
  save(permuted_preds_genus_list, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/inter_intra_lifestyle_genus_permuted.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/inter_intra_lifestyle_genus_permuted.RData")
}
toc()
# create correlation list
get_lm_list_from_list <- function(x) {
  x %>% 
    pivot_longer(cols = c("pred_western", "pred_non_western"), names_to = "training_set", values_to = "pred") %>%
    get_lm_list(pred_df = ., grouping = c("study", "lifestyle", "training_set"))
}

lifestyle_lm_perm_list <- lapply(permuted_preds_genus_list, get_lm_list_from_list)
lifestyle_lm_perm <- bind_rows(lifestyle_lm_perm_list, .id = "source")


# for shufffled data, with permutations
combined_lm_list <- lifestyle_lm_perm %>% 
  bind_rows(., lifestyle_lm_list) %>%
  filter(!training_set %in% c("non_westernized", "westernized")) %>%
  group_by(study, lifestyle, training_set) %>%
  summarize(Coefficient = mean(Coefficient),
            R2 = mean(R2),
            RSME = mean(RSME)) %>%
  ungroup

df_p_val_lifestyle_perm <- combined_lm_list %>%
  arrange(study) %>%
  mutate(training_set = factor(training_set)) %>%
  rstatix::group_by(lifestyle) %>%
  rstatix::wilcox_test(R2 ~ training_set, paired = T) %>%  
  # take care, paired is not possible between permuted and not permuted
  # as both sets contain different studies within the same lifestyle
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>% 
  # rstatix::add_significance(p.col = "p", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>% 
  rstatix::add_xy_position(x = "lifestyle", dodge = 0.8) 


combined_lm_list %>%
  ggplot(., aes(x=lifestyle, y = R2)) +
  geom_boxplot(aes(fill = training_set)) +
  xlab("Test set") +
  add_pvalue(df_p_val_lifestyle_perm,
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


combined_lm_list %>%
  group_by(lifestyle, training_set) %>%
  summarize(mean_r2 = mean(R2))


# check differences to the combined model
combined_lm_list_diffs <- lifestyle_lm_perm %>% 
  bind_rows(., lifestyle_lm_list) %>%
  filter(training_set != "combined") %>%
  left_join(.,lifestyle_lm_list %>% 
              filter(training_set == "combined") %>%
              select(-model_name, -training_set),
            by = c("study", "lifestyle"),
            suffix = c("_subset", "_combined")) %>%
  mutate(diff_to_combined = Coefficient_subset - Coefficient_combined)


df_p_val_lifestyle_perm_diffs <- combined_lm_list_diffs %>%
  group_by(study, lifestyle, training_set) %>%
  summarise(diff_to_combined = mean(diff_to_combined)) %>%
  arrange(study) %>%
  mutate(training_set = factor(training_set)) %>%
  rstatix::group_by(lifestyle) %>%
  rstatix::wilcox_test(diff_to_combined ~ training_set, paired = T) %>%  
  rstatix::add_xy_position(x = "lifestyle", dodge = 0.8) %>%
  # filter(!(group1 == "non_westernized" & group2 == "pred_western")) %>%
  # filter(!(group1 == "pred_non_western" & group2 == "westernized")) %>%
  # filter(!(group1 == "non_westernized" & group2 == "westernized")) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>% 
  rstatix::add_significance(p.col = "p", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1))
  

combined_lm_list_diffs %>%
  group_by(study, lifestyle, training_set) %>%
  summarise(diff_to_combined = mean(diff_to_combined)) %>%
  ggplot(., aes(x=lifestyle, y = diff_to_combined)) +
  # geom_violin(aes(fill = training_set)) +
  geom_boxplot(aes(fill = training_set)) +
  add_pvalue(df_p_val_lifestyle_perm_diffs,
             label = "{p.adj.signif}",
             # step.group.by = "variation",
             step.increase = 0.05,
             tip.length = 0.01,
             # bracket.nudge.y = 0.02,
             xmin = "xmin",
             xmax = "xmax",
             show.legend = FALSE) +
  theme_classic()
  

