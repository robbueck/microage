# create model on industrialized or non-industrialized studies and predict age in the respective other set of studies
# downsample each lifestyle to eaqual sizes
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
library(Boruta)
library(gamlss)
source("./alt_models.R")
source("./regression_functions.R")
source("./lifestyle_regression_functions.R")

# Switches ######################################
resample_genus_step <- F
resample_genus_sick_step <- F

# get the list of instances downsampled to a specific size
get_list <- function(vl, dfrm, max = NA, study_prob = F) {
  lst <- lapply(c(industrialized = "industrialized", non_industrialized = "non_industrialized"),
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
  trainDF <- subset(oldDF, study != test_set & health == "healthy")
  testDF <- subset(oldDF, study == test_set | health != "healthy")
  ps_train <- ps
  ps_test <- ps
  rm(ps, oldDF)
  # subsample train set to equal amount of studies/individuals/samples per lifestyle
  stds <- get_list("study", trainDF, max = max_studies)
  df_red <- trainDF %>%
    filter(lifestyle == "non_industrialized" | study %in% stds$industrialized) # subsample studies
  sbjcts <- get_list("subject_ID", trainDF, max = max_subj, study_prob = T)
  df_red_sb <- df_red %>%
    filter(lifestyle == "non_industrialized" | subject_ID %in% sbjcts$industrialized) %>%
    filter(lifestyle == "industrialized" | subject_ID %in% sbjcts$non_industrialized) # subsample individuals
  df_red_y <- filter(df_red_sb, age <= 365) # same distribuition for above and below one year respectively
  df_red_o <- filter(df_red_sb, age > 365)
  smpls_y <- get_list("sample_ID", df_red_y, max = max_samples / 2, study_prob = T)
  df_red_y_smp <- df_red_y %>%
    filter(lifestyle == "non_industrialized" | sample_ID %in% smpls_y$industrialized) %>%
    filter(lifestyle == "industrialized" | sample_ID %in% smpls_y$non_industrialized)
  if(length(table(df_red_o$lifestyle)) > 1) {
    smpls_o <- get_list("sample_ID", df_red_o, max = max_samples / 2, study_prob = T)
    df_red_o_smp <- df_red_o %>%
      filter(lifestyle == "non_industrialized" | sample_ID %in% smpls_o$industrialized) %>%
      filter(lifestyle == "industrialized" | sample_ID %in% smpls_o$non_industrialized) # subsample samples
    df_red_smp <- rbind(df_red_y_smp, df_red_o_smp)
  } else {
    print("No old samples found for one lifestyle")
    df_red_smp <- df_red_y_smp
  }
  sample_data(ps_train) <- sample_data(df_red_smp)
  sample_data(ps_test) <- sample_data(testDF)
  return(list(test = ps_test, train = ps_train))
}

indust_non_indust_preds <- function(test_set, ps, extra_cols = c("Observed", "Shannon"),
                                max_samples = NA, max_studies = NA, max_subj = NA,
                                run_importances = T) {
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
                                  run_importances = F)
  combined_preds <- combined_res$preds %>%
    select(sample_ID, subject_ID, study, lifestyle, age, pred) %>%
    mutate(pred_combined = pred, .keep = "unused")
  
  # industrialized
  print("run industrialized model")
  test_train_set_industrialized <- test_train_set
  test_train_set_industrialized$train <- test_train_set_industrialized$train %>% subset_samples(lifestyle == "industrialized")
  industrialized_res <- get_predictions(ps_test = test_train_set_industrialized$test, 
                              ps_train = test_train_set_industrialized$train,
                              extra_cols = extra_cols,
                              run_importances = run_importances,
                              test_set_name = test_set)
  industrialized_preds <- industrialized_res$preds %>%
    select(sample_ID, subject_ID, study, lifestyle, age, pred) %>%
    mutate(pred_industrialized = pred, .keep = "unused")
  
  # non industrialized
  print("run non industrialized model")
  test_train_set_non_industrialized <- test_train_set
  test_train_set_non_industrialized$train <- test_train_set_non_industrialized$train %>% subset_samples(lifestyle == "non_industrialized")
  non_industrialized_res <- get_predictions(ps_test = test_train_set_non_industrialized$test,
                                  ps_train = test_train_set_non_industrialized$train,
                                  extra_cols = extra_cols,
                                  run_importances = run_importances,
                                  test_set_name = test_set)
  non_industrialized_preds <- non_industrialized_res$preds %>%
    select(sample_ID, subject_ID, study, lifestyle, age, pred) %>%
    mutate(pred_non_industrialized = pred, .keep = "unused")
  # combine
  all_preds <- full_join(industrialized_preds, non_industrialized_preds) %>%
    full_join(.,combined_preds)
  cat("finished with ", test_set, "\n")
  return(list(preds = all_preds, ft_imps = list(industrialized = industrialized_res, non_industrialized = non_industrialized_res)))
}


get_predictions <- function(ps_test, ps_train, extra_cols = c("Observed", "Shannon"),
                            test_set_name = NULL, run_importances = T){
  tic()
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
  if(run_importances) {
    # get boruta importance:
    print("finished model, run Boruta")
    selected_features <- get_boruta_important_features(train_data = for_caret_list_train,
                                                       model = rf_model,
                                                       study = test_set_name)
    # run shap
    print("Get SHAP values")
    shap_out <- get_shap_long(rf_model,
                              test_data = for_caret_list_train$features,
                              features = unique(selected_features$importances$taxon))
  } else {
    shap_out = NULL
    selected_features = NULL
  }
  toc()
  return(list(preds = predictions,
              boruta = selected_features, 
              shap = shap_out, 
              test_set = test_set_name,
              training_data = for_caret_list_train$features))
}

run_all <- function(ps, n_run, run_importances = T){
  print(n_run)
  all_studies <- ps@sam_data$study %>% unique()
  all_studies <- all_studies[!all_studies %in% c( "gibson_2016", "ryan_2019", "kamdar_2020" )]
  names(all_studies) <- all_studies
  preds_imps <- all_studies %>% 
    future_map(~ indust_non_indust_preds(. ,ps = ps, run_importances = run_importances), progress = T)
  return(preds_imps)
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

cl <- makePSOCKcluster(ceiling(sqrt(n_cores)))
registerDoParallel(cl)
getDoParWorkers()
set.seed(825)
future::plan(multisession, workers = ceiling(sqrt(n_cores)))
options(ranger.num.threads = ceiling(sqrt(n_cores)))

# genus data ###################################################################

ps_object_genus_raw <- readRDS("../data/all_phyloseq_rf_filter_genus.rds") %>%
  subset_samples(age <= 730 & age > 1)
metadata <- ps_object_genus_raw %>%
  sample_data() %>%
  data.frame 
ps_object_genus_raw@sam_data$health <- "healthy"

# add sick children for prediction:
ps_sub <- readRDS("../data/phyloseq_subramanian_2014_malnurished.rds")
ps_gehr <- readRDS("../data/phyloseq_gehrig_2019_malnurished.rds")
ps_gehr@sam_data$sample_sum <- sample_sums(ps_gehr)
ps_sub@sam_data$sample_sum <- sample_sums(ps_sub)
merged_ps_genus <- merge_phyloseq(ps_sub %>% aggregate_taxa(level = "genus"),
                                  ps_gehr %>% aggregate_taxa(level = "genus")) %>%
  subset_samples(., health == "SAM")
ps_gibson <- readRDS("../data/phyloseq_gibson_2016.rds")
ps_ryan <- readRDS("../data/phyloseq_ryan_2019.rds")
ps_kamdar <- readRDS("../data/phyloseq_kamdar_2020.rds")

merged_ps_preterm <- merge_phyloseq(ps_gibson %>% aggregate_taxa(level = "genus"),
                                    ps_ryan %>% aggregate_taxa(level = "genus"),
                                    ps_kamdar %>% aggregate_taxa(level = "genus"))

ps_object_genus_raw_healthy_sick <- merge_phyloseq(merged_ps_genus, merged_ps_preterm, ps_object_genus_raw)

if(resample_genus_sick_step){
  # permuted_preds_genus <- get_permute_preds(ps = ps_object_genus_raw, blocks = "study")
  downsampled_genus_preds_sick <- future_map(1:50, ~ run_all(ps = ps_object_genus_raw_healthy_sick, 
                                                             run_importances = F,
                                                             n_run = .x))
  save(downsampled_genus_preds_sick, file = "../data/inter_intra_lifestyle_genus_downsampled_sick.RData")
} else {
  load("../data/inter_intra_lifestyle_genus_downsampled_sick.RData")
}

get_preds <- function(x) {
  x <- lapply(x, `[[`, "preds") %>% bind_rows()
}
downsampled_genus_df <- lapply(downsampled_genus_preds_sick, get_preds) %>%
  bind_rows(., .id = "source") %>%
  group_by(sample_ID, subject_ID, study, lifestyle, age) %>%
  summarize(pred_industrialized = mean(pred_industrialized, na.rm = T),
            pred_non_industrialized = mean(pred_non_industrialized, na.rm = T),
            pred_combined = mean(pred_combined, na.rm = T),
            .groups = "drop")
save(downsampled_genus_df, 
     file = "..//data/downsampled_genus_df.RData")


# downsampling healthy only #########################################
# requires lots of memory and time
tic()
if(resample_genus_step){
  # permuted_preds_genus <- get_permute_preds(ps = ps_object_genus_raw, blocks = "study")
  downsampled_genus_preds_imps <- future_map(1:50, ~ run_all(ps = ps_object_genus_raw, n_run = .x))
  save(downsampled_genus_preds_imps, file = "../data/inter_intra_lifestyle_genus_downsampled.RData")
} else {
  load("../data/inter_intra_lifestyle_genus_downsampled.RData")
}
toc()

# create correlation list
get_lm_list_from_list <- function(x) {
  x <- lapply(x, `[[`, "preds") %>% bind_rows()
  x %>% 
    pivot_longer(cols = c("pred_industrialized", "pred_non_industrialized"), names_to = "training_set", values_to = "pred") %>%
    get_lm_list(pred_df = ., grouping = c("study", "lifestyle", "training_set"))
}

lifestyle_lm_down <- lapply(downsampled_genus_preds_imps, get_lm_list_from_list) %>%
  bind_rows(., .id = "source") %>%
  group_by(study, lifestyle, training_set) %>%
  summarize(Coefficient = mean(Coefficient),
            R2 = mean(R2),
            RSME = mean(RSME)) %>%
  ungroup


# default model:
load("../data/all_nested_cv_genus_no_ls.RData")
genus_nested_cv_preds_long <- genus_no_ls_nested_cv_preds %>%
  pivot_longer(cols = c("rf1", "lasso"), names_to = "model_name", values_to = "pred" )
all_lm_list <- get_lm_list(pred_df = genus_nested_cv_preds_long, grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "combined")

lifestyle_lm_list <- bind_rows(lifestyle_lm_down, all_lm_list)

df_p_val_lifestyle_genus <- lifestyle_lm_list %>%
  arrange(study) %>%
  mutate(training_set = factor(training_set)) %>%
  rstatix::group_by(lifestyle) %>%
  rstatix::wilcox_test(R2 ~ training_set, paired = T) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>% 
  rstatix::add_xy_position(x = "lifestyle", dodge = 0.8) 



# number of important features ###############
# check if the number of important features in industrialized and non-industrialized models
# differ in the downsampling analysis too
important_fts_resampling <- downsampled_genus_preds_imps %>%
  map2_dfr(., seq_along(.),
           function(main_elem, resample_idx) {
             map2_dfr(main_elem, names(main_elem),
                      function(second_elem, study) {
                        map2_dfr(second_elem$ft_imps, names(second_elem$ft_imps),
                                 function(third_elem, lifestyle) {
                                   data.frame(
                                     resample_idx = resample_idx,
                                     study = study,
                                     lifestyle = lifestyle,
                                     finalDecision = third_elem$boruta$boruta_out$finalDecision,
                                     taxon_name = names(third_elem$boruta$boruta_out$finalDecision)
                                   )
                                 })
                      })
           })

important_fts_resampling_ls_filter <- important_fts_resampling %>%
  mutate(model_lifestyle = ifelse(lifestyle == "industrialized", yes = "industrialized", no = "non_industrialized"),
         .keep = "unused") %>%
  left_join(., metadata %>% select(study, lifestyle) %>% distinct, by = "study") #%>%
  # filter(model_lifestyle == lifestyle) # check only models where validation lifestyle = training lifestyle

# each study was used the same amount as a validation set
important_fts_resampling_ls_filter %>% 
  rownames_to_column() %>%
  select(resample_idx, study, lifestyle) %>% 
  distinct() %>% 
  pull(study) %>% table

n_imp_fts_ls_rs_filter <- important_fts_resampling_ls_filter %>% 
  filter(finalDecision == "Confirmed") %>%
  group_by(resample_idx, lifestyle, study) %>%
  summarize(counts = n()) %>% 
  ungroup()
n_imp_fts_ls_rs_filter %>% rstatix::wilcox_test(counts ~ lifestyle) # significant difference between lifestyles
  
# difference in number of important features
ggplot(n_imp_fts_ls_rs_filter,aes(x = lifestyle, y = counts)) +
  geom_boxplot() +
  geom_jitter()

n_imp_fts_ls_rs_filter  %>% 
  group_by(lifestyle) %>%
  summarize(mean = mean(counts))

n_imp_fts <- important_fts_resampling_ls_filter %>% 
  filter(finalDecision == "Confirmed") %>%
  group_by(resample_idx, model_lifestyle, taxon_name) %>%
  summarize(counts = n()) %>% 
  ungroup() %>%
  filter(counts > 2) %>% # from each resampling, remove all taxas that were important in only two or less models
  group_by(resample_idx, model_lifestyle) %>%
  summarize(counts = length(unique(taxon_name))) #%>%
  # mutate(counts = ifelse(lifestyle == "industrialized", yes = counts/13, no = counts/8))

n_imp_fts %>%
  ggplot(., aes(x = model_lifestyle, y = counts)) +
  geom_boxplot() +
  geom_jitter() +
  ylim(0, NA)

# because I don't use the same resampling per study
n_imp_fts_per_study <- important_fts_resampling_ls_filter %>% 
  filter(finalDecision == "Confirmed") %>%
  group_by(resample_idx, model_lifestyle, study) %>%
  summarize(counts = length(unique(taxon_name))) %>% 
  ungroup()

n_imp_fts_per_study %>%
  ggplot(., aes(x = model_lifestyle, y = counts)) +
  geom_boxplot() +
  geom_jitter() +
  ylim(0, NA)

# means in each lifestyle:
n_imp_fts_per_study %>%
  group_by(model_lifestyle) %>%
  summarize(mean_counts = mean(counts))


# importance vs prevalence in models ###############
# get prevalences from shap-output
shap_res_resampling <- downsampled_genus_preds_imps %>%
  map2_dfr(., seq_along(.),
           function(main_elem, resample_idx) {
             map2_dfr(main_elem, names(main_elem),
                      function(second_elem, study) {
                        map2_dfr(second_elem$ft_imps, names(second_elem$ft_imps),
                                 function(third_elem, lifestyle) {
                                   third_elem$shap %>%
                                     mutate(resample_idx = resample_idx)
                                 })
                      })
           }) %>%
  left_join(., metadata %>% select(run_accession, study, lifestyle))

# get prevalences in each resampling set
prevs_per_resampling <- shap_res_resampling %>%
  group_by(resample_idx, lifestyle, taxon) %>%
  summarize(prev = mean(ab_value > 0)) %>%
  ungroup() %>%
  mutate(lifestyle = ifelse(lifestyle == "non_industrialized", yes = "non_industrialized", no = "industrialized"))

# get importances
importances_res_resampling <- downsampled_genus_preds_imps %>%
  map2_dfr(., seq_along(.),
           function(main_elem, resample_idx) {
             map2_dfr(main_elem, names(main_elem),
                      function(second_elem, study) {
                        map2_dfr(second_elem$ft_imps, names(second_elem$ft_imps),
                                 function(third_elem, lifestyle) {
                                   # imp_fts <- third_elem$boruta$boruta_out$finalDecision
                                   # imp_fts <- names(imp_fts)[imp_fts == "Confirmed"]
                                   third_elem$boruta$importances %>%
                                     # filter(taxon %in% imp_fts) %>%
                                     mutate(resample_idx = resample_idx,
                                            lifestyle = lifestyle)
                                 })
                      })
           })

# compare importances with prevalences:
imp_prevs <- left_join(importances_res_resampling, prevs_per_resampling) %>%
  left_join(., important_fts_resampling %>% 
              filter(finalDecision == "Confirmed") %>%
              group_by(resample_idx, lifestyle, taxon_name) %>%
              summarize(counts = n()) %>%
              mutate(taxon = taxon_name))

imp_prevs %>%
  filter(counts > 1) %>%
  group_by(taxon, lifestyle) %>%
  summarize(prev = mean(prev),
            importance = mean(importance)) %>%
  ggplot(., aes(x = prev, y = importance, color = lifestyle)) +
  geom_point() +
  geom_smooth(method = "loess") +
  # facet_wrap(~resample_idx) +
  ggpubr::stat_cor(# aes(colour = NULL),
    method = "spearman",
    cor.coef.name = "spearman",
    size = 3,
    show.legend = F)


cors_per_ls_resample <- imp_prevs %>%
  filter(counts > 1) %>%
  group_by(taxon, lifestyle, resample_idx) %>%
  summarize(prev = mean(prev),
            importance = mean(importance)) %>%
  group_by(lifestyle, resample_idx) %>%
  rstatix::cor_test(., vars = c("prev", "importance"), method = "spearman")


# total number of taxa in each resampling set ##################################
n_fts_resampling <- downsampled_genus_preds_imps %>%
  map2_dfr(., seq_along(.),
           function(main_elem, resample_idx) {
             map2_dfr(main_elem, names(main_elem),
                      function(second_elem, study) {
                        map2_dfr(second_elem$ft_imps, names(second_elem$ft_imps),
                                 function(third_elem, lifestyle) {
                                   taxa_names <- third_elem$training_data %>%
                                     select(-Observed, -Shannon) %>%
                                     summarise(across(everything(), ~mean(.x > 0, na.rm = TRUE))) %>%
                                     keep(~ .x > 0) %>%
                                     names(.) %>%
                                     unique()
                                   data.frame(taxa_names = taxa_names, 
                                              resample_idx = resample_idx,
                                              study = study,
                                              lifestyle = lifestyle)
                                 })
                      })
           })
print("A")
# number of taxa shared/unique between lifestyles
n_taxa_per_lifestyle <- n_fts_resampling %>%
  group_by(study, resample_idx, taxa_names) %>%
  summarise(lifestyle_n = n_distinct(lifestyle), .groups = "drop",
            lifestyle_i = sum(lifestyle == "industrialized"),
            lifestyle_ni = sum(lifestyle == "non_industrialized")) %>%
  group_by(study, resample_idx) %>%
  summarize(unique_taxa = sum(lifestyle_n ==1), .groups = "drop",
            shared_taxa = sum(lifestyle_n == 2),
            industrialized_taxa = sum(lifestyle_ni == 0),
            non_industrialized_taxa = sum(lifestyle_i == 0))
print("B")
# unique taxa per lifestyle
unique_per_lifestyle <- n_taxa_per_lifestyle %>%
  pivot_longer(cols = c("industrialized_taxa", "non_industrialized_taxa"),
               names_to = "category", values_to = "count") %>%
  mutate(category = factor(category, levels = c("non_industrialized_taxa", "industrialized_taxa")))

unique_per_lifestyle %>% rstatix::wilcox_test(count ~ category, paired = T)
ggplot(unique_per_lifestyle,aes(x  = category, y = count)) +
  geom_boxplot() +
  geom_jitter() +
  ylim(0,NA)

n_taxa_per_lifestyle %>%
  pivot_longer(cols = c("industrialized_taxa", "non_industrialized_taxa"),
               names_to = "category", values_to = "count") %>%
  rstatix::wilcox_test(count ~ category, paired = T)
summary(n_taxa_per_lifestyle)

# number of taxa per lifestyle:
ft_counts <- n_fts_resampling %>%
  mutate(lifestyle = ifelse(lifestyle == "industrialized", yes = "industrialized", no = "non_industrialized"),
         .keep = "unused") %>%
  group_by(resample_idx, study, lifestyle) %>%
  summarize(n_taxa = length(unique(taxa_names))) %>%
  ungroup()

# significant difference between lifestyles
ft_counts %>% rstatix::wilcox_test(n_taxa ~ lifestyle, paired = T)
# effect size:
ft_counts %>% 
  group_by(lifestyle) %>%
  summarize(mean_taxa = mean(n_taxa))

ggplot(ft_counts,aes(x = lifestyle, y = n_taxa))+
  geom_boxplot() +
  geom_jitter() +
  ylim(0,NA)
            
n_imps_n_total <- left_join(n_imp_fts_ls_rs_filter, ft_counts, by = c("resample_idx", "lifestyle", "study"))
ggplot(n_imps_n_total, aes(x = counts, y = n_taxa, color = lifestyle)) +
  geom_point(alpha = 0.3)+
  ylim(0,NA) +
  xlim(0,NA) +
  geom_smooth()


save(lifestyle_lm_list, n_imp_fts_per_study,
     file = "../data/downsampling_data.RData")
