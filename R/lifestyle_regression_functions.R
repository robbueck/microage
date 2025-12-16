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
