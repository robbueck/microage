# functions to run the regression scripts
library(caret)
library(rsample)
library(tidyverse)
library(ranger)
library(caret)
library(caretEnsemble)
library(ggplot2)
library(gbm)
library(lme4)

# fix merge_taxa2 naming
merge_taxa2_fixed <- function (x, taxa = NULL, pattern = NULL, name = "Merged")
{
  if (is.null(taxa) && is.null(pattern)) {
    return(x)
  }
  if (!is.null(pattern)) {
    if (!is.null(taxa)) {
      mytaxa <- taxa
    }
    else {
      mytaxa <- taxa(x)
    }
    if (length(grep(pattern, mytaxa)) == 0) {
      return(x)
    }
    mytaxa <- mytaxa[grep(pattern, mytaxa)]
  }
  else if (is.null(taxa)) {
    mytaxa <- taxa(x)
  }
  else {
    mytaxa <- taxa
  }
  x2 <- phyloseq::merge_taxa(x, mytaxa, 1)
  mytaxa <- gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", mytaxa))
  taxa_names(x2) <- gsub(mytaxa[[1]], name, taxa_names(x2))
  tax_table(x2)[1, ] <- rep(name, ncol(tax_table(x2)))
  x2
}



# function to run the regression methods
# k: The number of partitions of the data set. If left as NULL (the default),
# v will be set to the number of unique values in the grouping variable, creating "leave-one-group-out" splits.
# outdated, not used anymore
run_age_cv <- function(data_set, meta_df, group, k = NULL, methods = c("adapt_ranger", "enet", "lasso", "glmnet", "gbm"),
                       rfGrid = NULL, xgbGrid = NULL){
  # create subsets for cv
  set.seed(100)
  splits <- group_vfold_cv(meta_df, group = group, v = k)
  caret_split <- rsample2caret(splits)
  # create subsets for cv
  tc_grouped <- trainControl(method = "cv", number = 5, 
                             savePredictions = "final",
                             selectionFunction = "oneSE", # function to select the best model
                             index = caret_split$index,
                             indexOut = caret_split$indexOut,
                             allowParallel = T)
  res_list <- setNames(vector("list", length(methods)), methods)
  
  # create functions:
  run_ranger <- function(...) {
    print("random forest")
    if(is.null(rfGrid)) {
      rfGrid <-  expand.grid(mtry = (1:10)*3, 
                             num.trees = c(20, 50, 200, 1000), 
                             min.node.size = (1:5)*4,
                             splitrule = c("variance"))
    }
    return(train(x = data_set, y = meta_df$age,
                 method = adapt_ranger, 
                 importance = "permutation",
                 trControl = tc_grouped,
                 tuneGrid = rfGrid)
           )
  }
  run_enet <- function(){
    print("elastic net")
    return(train(x = data_set, y = meta_df$age,
                 method = "enet", 
                 trControl = tc_grouped,
                 verbose = FALSE)
    )
  }
  run_lasso <- function(){
    print("lasso regression")
    return(train(x = data_set, y = meta_df$age, method = "lasso", trControl = tc_grouped))
  }
  run_glmnet <- function(){
    print("glmnet")
    return(train(x = data_set, y = meta_df$age,method = "glmnet", trControl = tc_grouped))
  }
  run_gbm <- function(){
    print("gradient boosting")
    return(train(x = data_set, y = meta_df$age,
                 method = "gbm", 
                 trControl = tc_grouped_subjects,
                 verbose = FALSE)
    )
  }
  run_xgbtree <- function(){
    if(is.null(xgbGrid)) {
      xgbGrid <- expand.grid(
      nrounds = seq(from = 200, to = 1000, by = 50),
      eta = c(0.025, 0.05, 0.1, 0.3),
      max_depth = c(2, 3, 4, 5, 6),
      gamma = 0,
      colsample_bytree = 1,
      min_child_weight = 1,
      subsample = 1
    )}
    xgbFit1_subjects <- train(x = data_set, y = meta_df$age,
                              method = "xgbTree",
                              tuneGrid = tune_grid_xgb_subjects,
                              trControl = tc_grouped_subjects)
  }
  
  # run functions:
  res_list <- setNames(map(methods, function(task) {
    if (task %in% c("adapt_ranger", "enet", "lasso", "glmnet", "gbm", "xgbTree")) {
      switch(task,
             "adapt_ranger" = run_ranger(),
             "enet" = run_enet(),
             "lasso" = run_lasso(),
             "glmnet" = run_glmnet(),
             "gbm" = run_gbm(),
             "xgbTree" = run_xgbtree(),
             )
    } else {
      return(NULL)  # Return NULL if the task is not found in list_total
    }
  }), methods)
  return(res_list)
}

# creates one data frame with the predictions and some metadata from a list of trained models
pred_from_list <- function(model_list, metad,
                           metad_cols = c("age", "study", "lifestyle", "country", "run_accession")){
  if(class(model_list) == "caretList"){
    res_df <- map(names(model_list), ~ model_list[[.x]]$pred %>% 
                         arrange(rowIndex) %>% 
                         cbind(metad[,metad_cols], .) %>%
                         mutate(model_name = .x)) %>%
      bind_rows()
  } else if(class(model_list) == "train") {
    res_df <- model_list$pred %>%
      arrange(rowIndex) %>%
      cbind(metad[,metad_cols], .) %>%
      mutate(model_name = model_list$method)
  } else {
    stop("model_list is neither an object of lass caretList or train")
  }
  res_df$residuals <- res_df$pred - res_df$age
  if(!all(res_df$age == res_df$obs)){
    warning("observed values in model do not correspond to true age in metadata")
  }
  return(res_df)
}

# get R2 and coefficient for each study for each model from a list of models
# go with R2 on rank transformed predictions
get_lm_list <- function(pred_df, grouping){
  pred_df <- pred_df %>%
    tidyr::unite(., col = "grouping_x_model_name", grouping, sep = "xxx") %>%
    mutate(pred = rank(pred))     # rank transform the predictions
  list_models <- lme4::lmList(data = pred_df, formula = pred ~ age | grouping_x_model_name)
  # list_models_sample_sum <- lme4::lmList(data = pred_df, formula = pred ~ age * sample_sum | grouping_x_model_name)
  model_df <- data.frame(grouping_x_model_name = names(list_models),
                         Coefficient = coef(list_models)$age,
                         # R2_sample_sum = summary(list_models_sample_sum)$r.squared,
                         R2 = summary(list_models)$r.squared)
  RSME_list <- pred_df %>%
    group_by(grouping_x_model_name) %>%
    dplyr::summarize(RSME = sqrt(mean(age - pred)^2),
                     n_samples = n()) %>%
    ungroup()
  model_df <- left_join(model_df, RSME_list) %>%
    separate(., grouping_x_model_name, into = c(grouping), sep = "xxx")
  return(model_df)
}


# function to create plots:
# if a list is entered, create each plot for each element of the list
# all_studies argument is for color compatibility with other plots. provide the 
# same list as in other plots, such as each study gets the same color in different plots
# independent of the other studies present or absent
# input is either a list of caret models (model_list) or
# a data frame containing the columns age, pred, color and model_name
modelplots <- function(model_list = NULL, all_results = NULL, metadata_df = NULL, all_studies = NULL, color = "study",
                       metadata_cols = c("age", "study", "lifestyle", "country", "run_accession")){
  if(!is.null(model_list) & !is.null(all_results)){
    stop("model_list and all_results are specified, decide for one of that!")
  } else if(!is.null(model_list) & is.null(metadata_df)){
    stop("model_list is given, but metadata missing")
  } else if(!is.null(all_results) & !is.null(metadata_df)){
    warning("all_results given, metadata_df will be ignored")
  }
  res_list = list()
  if(!is.null(model_list)){
    all_results <- pred_from_list(model_list = model_list, metad = metadata_df, metad_cols = metadata_cols)
  } else if(!is.null(all_results)){
    all_results$residuals <- all_results$pred - all_results$age
  }
  coeff_df <- get_lm_list(pred_df = all_results, grouping = c(color, "model_name"))
  
    # only when color compatibility is required
  if(!is.null(all_studies)){
    studies_here <- unique(all_results[[color]]) %in% all_studies
    all_results <- all_results %>% filter(.data[[color]] %in% all_studies)
    coeff_df <- coeff_df %>% filter(.data[[color]] %in% all_studies)
  }
  # predicted vs true age
  res_list$correlation <- ggplot(all_results, aes(x = age, y = pred, color = .data[[color]])) +
    geom_point(alpha = 0.2, size = 0.7) +
    theme_classic() +
    xlab("True age [days]") +
    ylab("Predicted age [days]") +
    stat_poly_line(fm.values = T) +
    ylim(0, max(all_results$pred)) +
    facet_wrap(~model_name) +
    geom_abline(slope = 1, color = "red", linetype = 2) +
    theme(text = element_text(size = 20),
          legend.position="none", #c(0.4, 0.05),
          legend.title = element_text(size=15),
          legend.text = element_text(size = 13)) 
  # residual plot
  res_list$residuals <- ggplot(all_results, aes(x = age, y = residuals, color = study)) +
    geom_point(alpha = 0.2, size = 0.7) +
    theme_classic() +
    facet_wrap(~model_name) +
    stat_poly_line(fm.values = T) +
    geom_abline(slope = 0, color = "red", linetype = 2) +
    theme(text = element_text(size = 20),
          legend.position="none", #c(0.4, 0.05),
          legend.title = element_text(size=15),
          legend.text = element_text(size = 13)) 
  # R2 and slope for true vs predicted age
  res_list$model <- ggplot(coeff_df, aes(x = R2, y = Coefficient, color = .data[[color]]))+#, size = Coefficient)) +
    geom_point(size = 3) +
    theme_classic() +
    # expand_limits(y=0)+
    facet_wrap(~model_name) +
    # labs(size = "slope") +
    xlab("R\u00b2") +
    ylab("Slope") +
    # geom_abline(slope = 1, color = "red", linetype = 2) +
    geom_vline(aes(xintercept = 0.5), linetype = "dashed", color = "red") +
    geom_hline(aes(yintercept = 0.5), linetype = "dashed", color = "red") +
    guides(colour = "none") +
    theme(#axis.title.y=element_blank(),
          # axis.text.y=element_blank(),
          # axis.ticks.y = element_blank(),
          # axis.line.y = element_blank(),
          text = element_text(size = 20),
          legend.position=c(0.17, 0.2),
          legend.title = element_text(size=15),
          legend.text = element_text(size = 13)) 

  res_list$model_heatmap <- ggplot(coeff_df, aes(model_name, .data[[color]])) +
    geom_tile(aes(fill = R2)) +
    geom_text(aes(label = round(R2, 2))) +
    scale_fill_gradient(low = "white", high = "red", limits = c(0,1))

  if(!is.null(model_list)){
    # importance plots
    if(class(model_list) == "caretList") {
      res_list$importance <- setNames(map(names(model_list), ~ varImp(model_list[[.x]]) %>%
                                  ggplot(., top = 10) +
                                  theme_classic() +
                                  ylab("relative importance") +
                                  xlab("") +
                                  ggtitle(.x)+
                                  theme(text = element_text(size = 20),
                                        legend.position="none", #c(1.5, 0.5),
                                  ))
                            , names(model_list))
    } else if(class(model_list) == "train"){
      res_list$importance <- varImp(res_family_study$rf1) %>%
        ggplot(., top = 10) +
        theme_classic() +
        ylab("relative importance") +
        xlab("") +
        theme(text = element_text(size = 20),
              legend.position="none", #c(1.5, 0.5),
        )
    } else{stop("model_list must be either of class caretList or train")}
  }
  if(!is.null(all_studies)){
    res_list$cor_plot <- res_list$cor_plot +
      scale_color_manual(values = turbo(length(all_studies))[studies_here], breaks = all_studies[studies_here] )
    res_list$res_plot <- res_list$res_plot +
      scale_color_manual(values = turbo(length(all_studies))[studies_here], breaks = all_studies[studies_here] )
    res_list$lm_plot <- res_list$lm_plot +
      scale_color_manual(values = turbo(length(all_studies))[studies_here], breaks = all_studies[studies_here] )
  }
  # res_list <- c(res_list, correlations = cor_plot, residuals = res_plot, linear_model = lm_plot)
  return(res_list)
}

# transform the otu matrix to binary presence absence format
binary_transform_out <- function(ps){
  otu <- otu_table(ps)
  t_a_r <- taxa_are_rows(otu)
  otu <- unclass(otu)
  otu <- otu > 0
  # cast from logical to double
  storage.mode(otu) <- "double"
  otu_table(ps) <- otu_table(otu, taxa_are_rows = t_a_r)
  return(ps)
}


# funciton to preprocess a phyloseq object to be used by caret. returns a list with
# the feature-dataframe and corresponding metadata-dataframe
# stud_prevalence_cutoff: in ho many studies must a samplefeature be present
# study_prevalence_cutoff: prevalence cutoff in how many samples a taxon must be detected to be considered present in this study
# mean_ab_cutoff: mean abundance cutoff, only counting non-zero abundances
# only_multiple_samples: if only subjects with multiple samples should be contained
# filter_features: if logical, filter features according to cutoffs or return all features
# if vector: logical/numerical/character vector of features to keep
create_caret_df <- function(ps_object, transformation = "identity", mean_ab_cutoff = 5e-5, 
                            study_prevalence_cutoff = 2, # detected in how many studies
                            prevalence_in_study_cutoff = 5, # 
                            additional_cols = c("Observed", "Shannon", "lifestyle_industrialized"),
                            rf_depth = 2000,
                            only_multiple_samples = F, filter_features = T){
  cols_to_keep <- filter_features
  if(length(filter_features) > 1){
    filter_features <- F
  }
  if(transformation == "identity" & mean_ab_cutoff < 1){
    warning("No transformation, mean_ab_cutoff below 1 will have no effect")
  }
  
  if(filter_features) {
    ps_object_comp <- ps_object %>%
      microbiome::transform(transform = "compositional") %>%
      filter_taxa(function(x) mean(x[x > 0]) > mean_ab_cutoff, TRUE) %>%  # mean abundance cutoff, only counting abundances > 0
      prune_samples(samples = (sample_sums(.) != 0))
    ps_object <- ps_object %>% prune_taxa(taxa_names(ps_object_comp),.)
  }
  
  if(transformation == "binary"){
    ps_object_trans <- binary_transform_out(ps_object)
  } else {
    ps_object_trans <- ps_object %>%
      microbiome::transform(transform = transformation)
  }
  
  if(only_multiple_samples){
    # remove single sample subjects
    if(!"subject_ID" %in% colnames(ps_object_trans@sam_data)){
      stop("subject_ID not in sam_data of ps_object")
    }
    multiple_samples <- ps_object_trans@sam_data %>%
      data.frame() %>%
      filter(duplicated(subject_ID)) %>%
      select(subject_ID) %>%
      unique %>%
      deframe
    # ps_object_trans <- ps_object_trans %>% subset_samples(subject_ID %in% multiple_samples)
    oldDF <- as(sample_data(ps_object_trans), "data.frame") 
    newDF <- subset(oldDF, subject_ID %in% multiple_samples) 
    sample_data(ps_object_trans) <- sample_data(newDF) 
  }
  
  # apply filter to raw data
  ps_object <- ps_object %>% 
    prune_samples(sample_names(ps_object_trans),.)

  
  if(any(c("Observed", "Shannon") %in% additional_cols)){
    # add richness and shannon diversity to the metadata
    if(any(sample_sums(ps_object) < rf_depth)) {
      ps_low_counts <- ps_object %>% prune_samples(samples = (sample_sums(.) < rf_depth))
      ps_object_raref <- ps_object %>% rarefy_even_depth(., sample.size = rf_depth) %>%
        merge_phyloseq(.,ps_low_counts)
    } else {ps_object_raref <- ps_object %>% rarefy_even_depth(., sample.size = rf_depth)}
    diversities <- estimate_richness(ps_object_raref, measures = c("Observed", "Shannon"))
    diversities <- diversities[sample_names(ps_object_trans),]
    sample_data(ps_object_trans)$Observed <- diversities$Observed
    sample_data(ps_object_trans)$Shannon <- diversities$Shannon
  }
  rm(ps_object, ps_low_counts, ps_object_raref)
  metadf <- sample_data(ps_object_trans) %>% data.frame()
  if("lifestyle_industrialized" %in% additional_cols){
    metadf$lifestyle_industrialized <- ifelse(metadf$lifestyle == "industrialized", yes = 1, no = 0)
  }
  if(!"study" %in% colnames(ps_object_trans@sam_data)){
    stop("study not in sam_data of ps_object, but required")
  }
  data_set_caret <- ps_object_trans%>%
    otu_table() %>%
    t() %>%
    data.frame() %>%
    # select(!Unknown) %>%
    cbind(metadf[,c(additional_cols, "study", "run_accession")]) %>%
    select(-run_accession)
  
  if(filter_features){
    # prevalence filtering, only keep taxa present in multiple studies
    cols_to_keep <- data_set_caret %>% group_by(study) %>%
      summarize(across(everything(), ~ sum(. > 0))) %>% 
      select(-study) %>%
      select(where(~ sum(. >= prevalence_in_study_cutoff) >= study_prevalence_cutoff)) %>% # adjust threshold, in how many studies can a taxon must be detected
      colnames()  
    if("lifestyle_industrialized" %in% additional_cols){
      cols_to_keep <- c(cols_to_keep, "lifestyle_industrialized") %>% unique # is either 1 or 0 for all samples of one study
    }
  } else if(length(cols_to_keep) == 1) { # filter_features = F, use all features
    cols_to_keep <- colnames(data_set_caret)
  } 
  # add missing columns
  missing_cols <- cols_to_keep[!cols_to_keep %in% colnames(data_set_caret)]
  data_set_caret[,missing_cols] <- 0
  data_set_caret <- data_set_caret[,cols_to_keep]
  return(list(features = data_set_caret, metadata = metadf))
}


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



# for lodocv:
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
  if(is.null(rf_model_list$rf1$trainingData)){
    training_data <- for_caret_list_train$features
    training_data$.outcome <- for_caret_list_train$metadata$age
    rf_model_list$rf1$trainingData <- training_data
  }
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
                file = paste0("..//data/nested_cv_dataset_models/",
                              prefix, "_", test_set, ".rds"))
      # Boruta feature importance
      print("Run Boruta")
      sel_feat <- Boruta(x=for_caret_list_train$features, y = for_caret_list_train$metadata$age,
                         num.trees = rf_model_list$rf1$bestTune$num.trees, 
                         mtry = rf_model_list$rf1$bestTune$mtry, num.threads = ceiling(n_cores/5))
      write_rds(sel_feat,
                file = paste0("../data/nested_cv_dataset_models/",
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
