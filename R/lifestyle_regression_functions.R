library(Boruta)


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
  features <- features[features %in% colnames(test_data)]
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

# get boruta importances from a dataset
get_boruta_important_features <- function(train_data, model, study){
  boruta_res <- Boruta(x=train_data$features, y = train_data$metadata$age,
                       num.trees = model$bestTune$num.trees, 
                       mtry = model$bestTune$mtry)
  imp_ftrs <- boruta_res$finalDecision %>% grep("Confirmed", value = T,.) %>% names
  var_imp <- model$finalModel$variable.importance %>%
    { .[names(.) %in% imp_ftrs] } %>%
    { (. / max(.)) * 100 }
  res_df <- data.frame(importance = var_imp,
                       taxon = names(var_imp),
                       ID = study) %>%
    arrange(-importance)# %>% head(n = 20)
  return(list(boruta_out = boruta_res, importances = res_df))
}


# shap analysis can be done also on the traing data
# https://stats.stackexchange.com/questions/615290/if-feature-importance-is-only-computed-based-on-training-set-does-it-mean-one-s
# load model without respective study
# test_set <- "wampach_2018"
# prefix <- "genus_data_no_ls"
# data <- bind_cols(genus_train_data$features, genus_train_data$metadata %>% select(study))
# get_shap_per_study <- function(test_set, prefix, data){
#   # load model and boruta data
#   model <- read_rds(file = paste0("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/nested_cv_dataset_models/",
#                           prefix, "_", test_set, ".rds"))
#   model <- model$rf1
#   boruta_res <- read_rds(file = paste0("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/nested_cv_dataset_models/",
#                           prefix, "_", test_set, "_boruta_res.rds"))
#   important_features <- boruta_res$finalDecision %>% grep("Confirmed", value = T,.) %>% names
#   # subset data
#   test_data <- test_data %>%
#     filter(study != test_set)
#   important_features <- important_features[important_features %in% colnames(test_data)]
#   # run shap
# }

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
                                          additional_cols = c("Observed", "Shannon"),
                                          only_multiple_samples = F)
  for_caret_list_test <- create_caret_df(ps_object = ps_test, transformation = "compositional",
                                         additional_cols = c("Observed", "Shannon"),
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


# run model in bins
run_in_bins <- function(ps_obj, cutoffs, prefix) {
  print(cutoffs)
  oldDF <- as(sample_data(ps_obj), "data.frame") 
  bin_DF <- subset(oldDF, age >= cutoffs[[1]] & age < cutoffs[2])
  ps_obj_bin <- ps_obj
  sample_data(ps_obj_bin) <- sample_data(bin_DF)
  stds_all <- unique(ps_obj_bin@sam_data$study)
  res <- c("industrialized", "non_industrialized") %>% 
    future_map_dfr(~ get_predictions_lifestyle(. ,ps = ps_obj))
  res$interval <- cutoffs[1]
  res$n_studies <- length(stds_all)
  return(res)
}


# read feature importances
get_feature_importance <- function(study, lst = "nonindustrialized") {
  train_object <- readRDS(paste0("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/nested_cv_dataset_models/",
                                 "genus_data_", lst, "_", study, ".rds"))
  boruta_object <- readRDS(paste0("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/nested_cv_dataset_models/",
                                  "genus_data_", lst, "_", study, "_boruta_res.rds"))
  imp_ftrs <- boruta_object$finalDecision %>% grep("Confirmed", value = T,.) %>% names
  # get, filter and normalize importances
  var_imp <- train_object$rf1$finalModel$variable.importance %>%
    { .[names(.) %in% imp_ftrs] } %>%
    { (. / max(.)) * 100 }
  data.frame(importance = var_imp,
             taxon = names(var_imp),
             ID = study) %>%
    arrange(-importance)# %>% head(n = 20)
}


# linear model to test for lifestyle effect in modeling shap and abundance:
ks_test_lifestyle <- function(tx, df) {
  print(tx)
  dr_flt <- df %>% filter(taxon == tx) %>%
    # filter(model == lifestyle) %>%
    select(shap_value, ab_value, model, study, subject_ID, age) 
  # # subsample to equal age dist:
  # print("match it")
  # m_out <- matchit(model_industrialized ~ age, 
  #         data = dr_flt,
  #         method = "cem", k2k = T, m.order = "farthest")
  # plot(m_out, type = "density", interactive = F,
  #      which.xs = ~age)
  # plot(summary(m_out))
  # # plot(m_out, type = "jitter", interactive = FALSE)
  # dr_flt_match <- match.data(m_out, data = dr_flt,
  #            drop.unmatched = T)
  print("run ks test")
  ks_test_res <- ks.test(shap_value ~ model, data = dr_flt)
  # return(ks_test_res)
  return(list(taxon = tx,
              # p_intercept = intercept_test$`Pr(>Chisq)`[2],
              # p_slope = slope_test$`Pr(>Chisq)`[2],
              ls_ratio = table(dr_flt$model)[1] / table(dr_flt$model)[2],
              n_samples = nrow(dr_flt),
              p_ks = ks_test_res$p.value,
              D_ks = ks_test_res$statistic))
}

# run ks test between taxa within one model:
ks_test_per_model <- function(data, taxa) {
  print(taxa)
  data_grouped <- data %>% filter(taxon %in% taxa)
  ks_test_res <- ks.test(shap_value ~ taxon, data = data_grouped)
  # return(ks_test_res)
  return(list(taxon1 = taxa[1],
              taxon2 = taxa[2],
              n_samples = nrow(data_grouped),
              p_ks = ks_test_res$p.value,
              D_ks = ks_test_res$statistic))
  return(ks_results)
}


# check empirical p-values:
get_p_value <- function(obs_value, null_dist) {
  # Calculate the proportion of null values greater than or equal to the observed value
  p_value <- mean(abs(null_dist) >= abs(obs_value))
  return(p_value)
}
