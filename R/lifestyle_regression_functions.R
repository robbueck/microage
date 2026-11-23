##Project: Microbiome age index
##Aim: functions for cross-lifestyle predictions
##Author: Robert Bücking

library(Boruta)

#' Get Boruta-confirmed features with their model importances
#'
#' Runs Boruta feature selection on the training data to find features that
#' are significantly relevant for predicting age, and returns the variable
#' importances of these confirmed features from a trained random forest model.
#'
#' @param train_data A list with two elements:
#'   \describe{
#'     \item{features}{A data frame or matrix of predictors (e.g. taxa), one
#'       row per sample.}
#'     \item{metadata}{A data frame with the same samples in the same order,
#'       containing the outcome column \code{age}.}
#'   }
#' @param model A random forest model trained with \code{caret::train()} using
#'   the \code{ranger} method, with variable importance enabled (e.g.
#'   \code{importance = "permutation"}). Its tuned \code{num.trees} and
#'   \code{mtry} are reused for Boruta.
#' @param study Character. An identifier (e.g. study or dataset name) added
#'   as column \code{ID} to the importance table, useful when combining
#'   results from several datasets.
#'
#' @details
#' The outcome is fixed to \code{train_data$metadata$age}. Boruta compares
#' each feature against randomly permuted "shadow" features using ranger
#' random forests. Only features with the final decision "Confirmed" are kept;
#' "Tentative" features are excluded.
#'
#' The importances in the result are not Boruta's own scores. They come from
#' \code{model$finalModel$variable.importance} and are scaled so that the most
#' important confirmed feature has the value 100.
#'
#' Boruta is random; call \code{set.seed()} before this function for
#' reproducible results.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{boruta_out}{The full \code{Boruta} result object.}
#'     \item{importances}{A data frame with the columns \code{importance}
#'       (scaled 0 to 100), \code{taxon} and \code{ID}, sorted by decreasing
#'       importance, with one row per confirmed feature.}
#'   }
#'
#' @seealso \code{\link[Boruta]{Boruta}}, \code{\link[caret]{train}}
#'
#' @importFrom Boruta Boruta
#' @importFrom dplyr arrange %>%
#'
#' @examples
#' #' set.seed(1)
#' res <- get_boruta_important_features(train_data, rf_model, study = "cohort_A")
#' head(res$importances)
#' }
#'
#' @export
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


#' Predict age for one lifestyle group from a model trained on the others
#'
#' Performs one step of a leave-one-lifestyle-out validation: trains a random
#' forest regression model for age on all samples except those of one
#' lifestyle group, then predicts age for the held-out group.
#'
#' @param test_set Character. The value of the \code{lifestyle} column whose
#'   samples are held out as test set (e.g. \code{"industrialized"}).
#' @param ps A \code{phyloseq} object. Its sample data must contain the
#'   columns \code{lifestyle}, \code{study}, \code{age}, \code{Observed} and
#'   \code{Shannon}.
#'
#' @details
#' Features are prepared with \code{create_caret_df()}: relative abundances
#' (\code{transformation = "compositional"}) plus observed richness and
#' Shannon diversity. In the training set, taxa with a mean abundance below
#' 5e-5 are removed. The test set is restricted to the features selected in
#' the training set.
#'
#' The model is tuned with grouped cross-validation, where each study forms
#' one group (\code{rsample::group_vfold_cv()}), so samples from the same
#' study are never split between training and validation folds. The tuning
#' uses the custom caret method \code{adapt_ranger} with permutation
#' importance and the grid:
#' \itemize{
#'   \item \code{mtry}: 3, 10, 20
#'   \item \code{min.node.size}: 10, 20
#'   \item \code{num.trees}: 200; \code{splitrule}: variance;
#'     \code{replace}: \code{FALSE}
#' }
#' Near-zero-variance features are removed before training. The final model
#' is selected with the one-standard-error rule on R², which picks the
#' simplest model within one SE of the best.
#'
#' Samples with a missing \code{lifestyle} are excluded from both training
#' and test set. Model training is random; call \code{set.seed()} beforehand
#' for reproducible results. \code{allowParallel = TRUE} only takes effect if
#' a parallel backend is registered (e.g. with \pkg{doParallel}).
#'
#' @return The metadata of the test samples (from \code{create_caret_df()})
#'   with an additional column \code{pred} holding the predicted age.
#'
#' @seealso \code{\link{get_boruta_important_features}},
#'   \code{\link[caret]{train}}, \code{\link[rsample]{group_vfold_cv}}
#'
#' @importFrom phyloseq sample_data sample_data<-
#' @importFrom rsample group_vfold_cv rsample2caret
#' @importFrom caret trainControl train
#' @importFrom dplyr mutate %>%
#' @importFrom methods as
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' lifestyles <- unique(sample_data(ps)$lifestyle)
#' preds <- purrr::map_dfr(lifestyles, get_predictions_lifestyle, ps = ps)
#' }
#'
#' @export
get_predictions_lifestyle <- function(test_set, ps){
  print(test_set)
  oldDF <- as(sample_data(ps), "data.frame") 
  trainDF <- subset(oldDF, lifestyle != test_set)
  
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


#' Run lifestyle predictions within an age bin
#'
#' Restricts the samples to one age range and runs
#' \code{\link{get_predictions_lifestyle}} for both lifestyle groups
#' (industrialized and non-industrialized) in parallel.
#'
#' @param ps_obj A \code{phyloseq} object with the sample data columns
#'   \code{age}, \code{study} and \code{lifestyle}.
#' @param cutoffs Numeric vector of length 2: lower (inclusive) and upper
#'   (exclusive) age limit of the bin.
#'
#' @return The combined predictions for both lifestyle groups, with the
#'   additional columns \code{interval} (lower bin limit) and \code{n_studies}
#'   (number of studies in the bin).
#'
#' @seealso \code{\link{get_predictions_lifestyle}}
#'
#' @importFrom furrr future_map_dfr
#'
#' @examples
#' \dontrun{
#' future::plan(future::multisession, workers = 2)
#' res <- run_in_bins(ps, cutoffs = c(0, 365))
#' }
#'
#' @export
run_in_bins <- function(ps_obj, cutoffs) {
  print(cutoffs)
  oldDF <- as(sample_data(ps_obj), "data.frame") 
  bin_DF <- subset(oldDF, age >= cutoffs[[1]] & age < cutoffs[2])
  ps_obj_bin <- ps_obj
  sample_data(ps_obj_bin) <- sample_data(bin_DF)
  stds_all <- unique(ps_obj_bin@sam_data$study)
  res <- c("industrialized", "non_industrialized") %>% 
    future_map_dfr(~ get_predictions_lifestyle(. ,ps = ps_obj_bin))
  res$interval <- cutoffs[1]
  res$n_studies <- length(stds_all)
  return(res)
}


#' Read Boruta-confirmed feature importances for a saved model
#'
#' Loads a saved random forest model and its Boruta result for one study,
#' and returns the importances of the Boruta-confirmed features, scaled so
#' that the most important feature has the value 100.
#'
#' @param study Character. Study identifier used in the file names.
#' @param lst Character. Lifestyle group used in the file names. Default
#'   \code{"nonindustrialized"}.
#'
#' @details
#' Reads two files from \code{../data/nested_cv_dataset_models/} (relative to
#' the working directory):
#' \code{genus_data_<lst>_<study>.rds}, which must contain the caret model as
#' element \code{rf1}, and \code{genus_data_<lst>_<study>_boruta_res.rds}.
#' Importances are taken from the model, not from Boruta.
#'
#' @return A data frame with the columns \code{importance}, \code{taxon} and
#'   \code{ID} (the study), sorted by decreasing importance.
#'
#' @seealso \code{\link{get_boruta_important_features}}
#'
#' @importFrom dplyr arrange %>%
#'
#' @examples
#' \dontrun{
#' imp <- get_feature_importance("backhed_2015", lst = "industrialized")
#' }
#'
#' @export
get_feature_importance <- function(study, lst = "non_industrialized") {
  train_object <- readRDS(paste0("../data/nested_cv_dataset_models/",
                                 "genus_data_", lst, "_", study, ".rds"))
  boruta_object <- readRDS(paste0("../data/nested_cv_dataset_models/",
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

#' Compare SHAP value distributions of one taxon between lifestyle models
#'
#' Runs a two-sample Kolmogorov-Smirnov test on the SHAP values of one taxon,
#' comparing the two levels of \code{model} (e.g. the industrialized and
#' non-industrialized models).
#'
#' @param tx Character. Name of the taxon to test.
#' @param df A data frame in long format with the columns \code{taxon},
#'   \code{shap_value}, \code{ab_value}, \code{model}, \code{study},
#'   \code{subject_ID} and \code{age}. \code{model} must have exactly two
#'   levels for the selected taxon.
#'
#' @return A list with:
#'   \describe{
#'     \item{taxon}{The tested taxon.}
#'     \item{ls_ratio}{Number of values for the first level of \code{model}
#'       divided by the number for the second (levels in alphabetical or
#'       factor order).}
#'     \item{n_samples}{Total number of values.}
#'     \item{p_ks}{P-value of the KS test.}
#'     \item{D_ks}{KS test statistic D.}
#'   }
#'
#' @importFrom dplyr filter select %>%
#' @importFrom stats ks.test
#'
#' @examples
#' \dontrun{
#' res <- purrr::map_dfr(unique(shap_df$taxon), ks_test_lifestyle, df = shap_df)
#' res$p_adj <- p.adjust(res$p_ks, method = "BH")
#' }
#'
#' @export
ks_test_lifestyle <- function(tx, df) {
  print(tx)
  dr_flt <- df %>% filter(taxon == tx) %>%
    # filter(model == lifestyle) %>%
    select(shap_value, ab_value, model, study, subject_ID, age) 
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

#' Compare SHAP value distributions between two taxa
#'
#' Runs a two-sample Kolmogorov-Smirnov test comparing the SHAP values of two
#' taxa within one model.
#'
#' @param data A data frame in long format with the columns \code{taxon} and
#'   \code{shap_value}, containing the SHAP values of a single model.
#' @param taxa Character vector of length 2 with the names of the two taxa to
#'   compare.
#'
#' @return A list with \code{taxon1}, \code{taxon2}, \code{n_samples} (total
#'   number of values), \code{p_ks} (KS p-value) and \code{D_ks} (KS test
#'   statistic D).
#'
#' @seealso \code{\link{ks_test_lifestyle}}
#'
#' @importFrom dplyr filter %>%
#' @importFrom stats ks.test
#'
#' @examples
#' \dontrun{
#' pairs <- combn(unique(shap_df$taxon), 2, simplify = FALSE)
#' res <- purrr::map_dfr(pairs, ks_test_per_model, data = shap_df)
#' }
#'
#' @export
ks_test_per_model <- function(data, taxa) {
  print(taxa)
  data_grouped <- data %>% filter(taxon %in% taxa)
  ks_test_res <- ks.test(shap_value ~ taxon, data = data_grouped)
  return(list(taxon1 = taxa[1],
              taxon2 = taxa[2],
              n_samples = nrow(data_grouped),
              p_ks = ks_test_res$p.value,
              D_ks = ks_test_res$statistic))
}


# check empirical p-values:
get_p_value <- function(obs_value, null_dist) {
  # Calculate the proportion of null values greater than or equal to the observed value
  p_value <- mean(abs(null_dist) >= abs(obs_value))
  return(p_value)
}
