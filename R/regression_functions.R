##Project: Microbiome age index
##Aim: functions to run the regression scripts
##Author: Robert Bücking

library(caret)
library(rsample)
library(tidyverse)
library(ranger)
library(caret)
library(caretEnsemble)
library(ggplot2)
library(gbm)
library(lme4)

# fix merge_taxa2 naming, obsolete with update of microbiome package
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

# shap analysis:
pfun <- function(object, newdata) {
  require(ranger)
  predict(object, data = newdata)$predictions
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


#' Calculate R² and slope of predicted versus true age per group
#'
#' Fits a linear regression of rank-transformed predicted age on true age
#' separately for each group (e.g. each combination of study and model), and
#' returns the slope, R², an RSME and the sample size per group.
#'
#' @param pred_df A data frame with one row per prediction, containing the
#'   columns \code{pred} (predicted age), \code{age} (true age) and the
#'   columns named in \code{grouping}.
#' @param grouping Character vector of column names that define the groups,
#'   e.g. \code{c("study", "model_name")}.
#'
#' @details
#' Predictions are replaced by their ranks before fitting, which makes R²
#' less sensitive to outliers and non-linear relationships. The ranks are
#' calculated across the whole data frame, not within each group. Models are
#' fitted with \code{lme4::lmList(pred ~ age | group)}.
#'
#' The grouping columns are temporarily combined into one column using
#' \code{"xxx"} as separator, so the group values must not contain this
#' string.
#'
#' @return A data frame with one row per group, containing the grouping
#'   columns, \code{Coefficient} (slope of the prediction rank on age),
#'   \code{R2}, \code{RSME} and \code{n_samples}.
#'
#' @seealso \code{\link[lme4]{lmList}}
#'
#' @importFrom tidyr unite separate
#' @importFrom dplyr mutate group_by summarize ungroup left_join n %>%
#'
#' @examples
#' \dontrun{
#' res <- get_lm_list(pred_df, grouping = c("study", "model_name"))
#' }
#'
#' @export
get_lm_list <- function(pred_df, grouping){
  pred_df <- pred_df %>%
    tidyr::unite(., col = "grouping_x_model_name", grouping, sep = "xxx") %>%
    mutate(pred_rank = rank(pred))     # rank transform the predictions
  list_models <- lme4::lmList(data = pred_df, formula = pred_rank ~ age | grouping_x_model_name)
  model_df <- data.frame(grouping_x_model_name = names(list_models),
                         Coefficient = coef(list_models)$age,
                         R2 = summary(list_models)$r.squared)
  RSME_list <- pred_df %>%
    group_by(grouping_x_model_name) %>%
    dplyr::summarize(RSME = sqrt(mean((age - pred)^2)),
                     n_samples = n()) %>%
    ungroup()
  model_df <- left_join(model_df, RSME_list) %>%
    separate(., grouping_x_model_name, into = c(grouping), sep = "xxx")
  return(model_df)
}


#' Create evaluation plots for age prediction models
#'
#' Creates a set of plots to evaluate one or more age prediction models:
#' predicted versus true age, residuals, R² and slope per group, and, if
#' models are given, variable importance. Takes either trained caret models
#' or a data frame of predictions.
#'
#' @param model_list A \code{caretList} or a single caret \code{train}
#'   object. Predictions are extracted with \code{pred_from_list()}. Requires
#'   \code{metadata_df}. Use either this or \code{all_results}.
#' @param all_results A data frame of predictions with the columns
#'   \code{age}, \code{pred}, \code{model_name}, \code{study} and the column
#'   named in \code{color}. Use either this or \code{model_list}.
#' @param metadata_df A data frame of sample metadata, required with
#'   \code{model_list}; ignored with \code{all_results}.
#' @param all_studies Character vector of all studies used across your plots.
#'   If given, each study gets the same colour in every plot, independent of
#'   which studies are present in this one. Default \code{NULL}.
#' @param color Character. Column used for colouring and grouping, e.g. for
#'   R² per group. Default \code{"study"}.
#' @param metadata_cols Character vector of metadata columns added to the
#'   predictions when using \code{model_list}.
#'
#' @details
#' R² and slope are calculated per group and model with
#' \code{\link{get_lm_list}}, i.e. on rank-transformed predictions. In the
#' R²/slope plot, dashed lines mark 0.5 on both axes.
#'
#' @return A named list of ggplot objects:
#'   \describe{
#'     \item{correlation}{Predicted versus true age per model, with a
#'       regression line per group and the identity line (red, dashed).}
#'     \item{residuals}{Residuals (predicted minus true age) versus true age.}
#'     \item{model}{R² versus slope per group, faceted by model.}
#'     \item{model_heatmap}{Heatmap of R² per group and model.}
#'     \item{importance}{Only with \code{model_list}: the ten most important
#'       features, as one plot or a list of plots per model.}
#'   }
#'
#' @seealso \code{\link{get_lm_list}}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_tile geom_text geom_abline
#'   geom_vline geom_hline facet_wrap theme_classic theme element_text xlab
#'   ylab ylim ggtitle guides scale_fill_gradient scale_color_manual
#' @importFrom ggpmisc stat_poly_line
#' @importFrom caret varImp
#' @importFrom viridisLite turbo
#' @importFrom purrr map
#' @importFrom dplyr filter %>%
#'
#' @examples
#' \dontrun{
#' plots <- modelplots(all_results = pred_df, all_studies = all_studies)
#' plots$correlation
#' }
#'
#' @export
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
    coord_cartesian(ylim = c(0, max(all_results$pred))) +
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
    xlab("R\u00b2") +
    ylab("Slope") +
    geom_vline(aes(xintercept = 0.5), linetype = "dashed", color = "red") +
    geom_hline(aes(yintercept = 0.5), linetype = "dashed", color = "red") +
    guides(colour = "none") +
    theme(text = element_text(size = 20),
          legend.position=c(0.17, 0.2),
          legend.title = element_text(size=15),
          legend.text = element_text(size = 13)) 

  res_list$model_heatmap <- ggplot(coeff_df, aes(model_name, .data[[color]])) +
    geom_tile(aes(fill = R2)) +
    geom_text(aes(label = round(R2, 2))) +
    scale_fill_gradient(low = "white", high = "red", limits = c(0,1))

  if(!is.null(model_list)){
    # importance plots
    if(inherits(model_list, "caretList")) {
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
  return(res_list)
}



#' Prepare a phyloseq object as input for caret models
#'
#' Filters and transforms the taxa of a phyloseq object, optionally adds
#' alpha diversity and lifestyle as predictors, and returns a feature table
#' and the matching metadata for use with \code{caret::train()}.
#'
#' @param ps_object A \code{phyloseq} object with taxa as rows in the OTU
#'   table. The sample data must contain \code{study} and
#'   \code{run_accession}, plus \code{subject_ID} if
#'   \code{only_multiple_samples = TRUE} and \code{lifestyle} if
#'   \code{"lifestyle_industrialized"} is in \code{additional_cols}.
#' @param transformation Character. Transformation of the abundances, passed
#'   to \code{microbiome::transform()} (e.g. \code{"identity"},
#'   \code{"compositional"}, \code{"clr"}), or \code{"binary"} for
#'   presence/absence via \code{binary_transform_out()}. Default
#'   \code{"identity"}.
#' @param mean_ab_cutoff Numeric. Taxa are kept only if their mean relative
#'   abundance, counting only samples where they are present, exceeds this
#'   value. Default \code{5e-5}.
#' @param study_prevalence_cutoff Integer. Minimum number of studies in which
#'   a taxon must be prevalent to be kept. Default 2.
#' @param prevalence_in_study_cutoff Integer. Minimum number of samples in
#'   which a taxon must be detected for it to count as prevalent in a study.
#'   Default 5.
#' @param additional_cols Character vector of extra predictors added to the
#'   features. \code{"Observed"} and \code{"Shannon"} add richness and Shannon
#'   diversity; \code{"lifestyle_industrialized"} adds 1 for industrialized
#'   and 0 for other samples. Other names are taken from the sample data.
#' @param rf_depth Integer. Rarefaction depth for calculating diversity.
#'   Default 2000.
#' @param only_multiple_samples Logical. If \code{TRUE}, only subjects with
#'   at least two samples are kept. Default \code{FALSE}.
#' @param filter_features Either a logical or a character vector:
#'   \describe{
#'     \item{\code{TRUE}}{(default) Filter taxa by mean abundance and
#'       prevalence.}
#'     \item{\code{FALSE}}{Keep all taxa.}
#'     \item{character vector}{Feature names to keep, e.g. the columns of a
#'       training set when preparing a test set. Features missing from the
#'       data are added as columns of zeros, so the result matches the
#'       training data.}
#'   }
#'
#' @details
#' The mean abundance filter is always applied to relative abundances,
#' independent of \code{transformation}, and is followed by the prevalence
#' filter on the transformed data. Both are skipped if \code{filter_features}
#' is \code{FALSE} or a vector of feature names.
#'
#' Samples with at least \code{rf_depth} reads are rarefied before calculating
#' diversity; samples with fewer reads use their raw counts. Rarefaction is
#' random, so call \code{set.seed()} beforehand for reproducible results.
#'
#' @return A list with:
#'   \describe{
#'     \item{features}{A data frame with one row per sample and one column
#'       per feature (taxa and additional predictors).}
#'     \item{metadata}{The sample data in the same order, including
#'       \code{Observed}, \code{Shannon} and \code{lifestyle_industrialized}
#'       if requested.}
#'   }
#'
#' @seealso \code{\link{get_predictions_lifestyle}}
#'
#' @importFrom phyloseq filter_taxa prune_taxa prune_samples sample_sums
#'   taxa_names sample_names sample_data otu_table rarefy_even_depth
#'   merge_phyloseq estimate_richness
#' @importFrom dplyr filter select group_by summarize across everything where
#'   %>%
#' @importFrom tibble deframe
#' @importFrom methods as
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' train <- create_caret_df(ps_train, transformation = "compositional")
#' test  <- create_caret_df(ps_test, transformation = "compositional",
#'                          filter_features = colnames(train$features))
#' }
#'
#' @export
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
  if(transformation == "compositional" & mean_ab_cutoff < 1){
    warning("transformation is not compositional, but filter will be applied to compositional values")
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
  if(taxa_are_rows(ps_object_trans)) {
    data_set_caret <- ps_object_trans%>%
      otu_table() %>%
      t() %>%
      data.frame() %>%
      cbind(metadf[,c(additional_cols, "study", "run_accession")]) %>%
      select(-run_accession)
  } else {
    data_set_caret <- ps_object_trans%>%
      otu_table() %>%
      data.frame() %>%
      cbind(metadf[,c(additional_cols, "study", "run_accession")]) %>%
      select(-run_accession)
  }
  
  if(filter_features){
    # prevalence filtering, only keep taxa present in multiple studies
    cols_to_keep <- data_set_caret %>% group_by(study) %>%
      summarize(across(everything(), ~ sum(. > 0))) %>% 
      select(-study) %>%
      select(where(~ sum(. >= prevalence_in_study_cutoff) >= study_prevalence_cutoff)) %>% # adjust threshold, in how many studies a taxon must be detected
      colnames()  
    if("lifestyle_industrialized" %in% additional_cols){
      cols_to_keep <- c(cols_to_keep, "lifestyle_industrialized") %>% unique # is either 1 or 0 for all samples of one study
    }
  } else if(length(cols_to_keep) == 1) { # filter_features = F, use all features
    cols_to_keep <- setdiff(colnames(data_set_caret), "study")
  } 
  # add missing columns
  missing_cols <- cols_to_keep[!cols_to_keep %in% colnames(data_set_caret)]
  data_set_caret[,missing_cols] <- 0
  data_set_caret <- data_set_caret[,cols_to_keep]
  return(list(features = data_set_caret, metadata = metadf))
}



#' Calculate SHAP values for a random forest model in long format
#'
#' Calculates approximate SHAP values with \pkg{fastshap} for a trained
#' caret random forest model, and returns them in long format together with
#' the scaled feature values, e.g. for SHAP summary or dependence plots.
#'
#' @param model A random forest model trained with \code{caret::train()}
#'   (ranger-based). For \code{features = "important"}, it must have been
#'   trained with variable importance enabled.
#' @param test_data A data frame of feature values for which SHAP values are
#'   calculated, with the same feature columns as the training data. If
#'   \code{NULL} (default), the model's training data is used.
#' @param nsim Integer. Passed to fastshap::explain. Default \code{10}.
#' @param features Which features to explain:
#'   \describe{
#'     \item{\code{NULL}}{(default) All features.}
#'     \item{\code{"important"}}{The 15 features with the highest variable
#'       importance in the model.}
#'     \item{character vector}{The named features; names not present in
#'       \code{test_data} are dropped.}
#'   }
#'
#' @details
#' SHAP values are approximated by Monte Carlo sampling with 10 simulations
#' (\code{nsim = 10}), so results vary between runs; call \code{set.seed()}
#' beforehand for reproducible results. Predictions are made with a function
#' \code{pfun(object, newdata)}, which must exist in the calling environment.
#' The calculation runs in parallel if a \pkg{foreach} backend is registered
#' (e.g. with \pkg{doParallel}).
#'
#' Feature values are scaled by dividing by the maximum of each feature, so
#' they range up to 1. This scaling is only meaningful for non-negative
#' features such as abundances.
#'
#' @return A long data frame with one row per sample and feature, containing:
#'   \describe{
#'     \item{sample_ID}{Row number of the sample in \code{test_data}.}
#'     \item{taxon}{Feature name.}
#'     \item{shap_value}{SHAP value of the feature for that sample.}
#'     \item{run_accession}{Row name of the sample in \code{test_data}.}
#'     \item{ab_value}{Feature value divided by the feature's maximum.}
#'   }
#'
#' @seealso \code{\link[fastshap]{explain}}, \code{\link{ks_test_lifestyle}}
#'
#' @importFrom fastshap explain
#' @importFrom dplyr select mutate across everything left_join %>%
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' shap_df <- get_shap_long(rf_model, test_data = test$features,
#'                          features = "important")
#' }
#'
#' @export
get_shap_long <- function(model, test_data = NULL, features = NULL, nsim= 10) {
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
    nsim = nsim,
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



#' Train final age prediction model
#'
#' Prepares the features of a phyloseq object with
#' \code{\link{create_caret_df}} and trains a lasso regression and a random
#' forest to predict age. Hyperparameters are tuned with grouped
#' cross-validation in which each study forms one fold.
#'
#' @param ps A \code{phyloseq} object with the training samples; see
#'   \code{\link{create_caret_df}} for the required sample data.
#' @param transform Character. Abundance transformation passed to
#'   \code{create_caret_df()}. Default \code{"compositional"}.
#' @param extra_cols Character vector of additional predictors passed to
#'   \code{create_caret_df()}. Default richness and Shannon diversity.
#' @param prefilter Logical. If \code{TRUE}, removes highly correlated
#'   features (absolute Pearson correlation > 0.8) with
#'   \code{caret::findCorrelation()} before training. Default \code{FALSE}.
#' @param rf_depth Integer. Rarefaction depth for diversity, passed to
#'   \code{create_caret_df()}. Default 2000.
#' @param pre_processing Character vector of caret pre-processing steps.
#'   Default \code{"nzv"} (remove near-zero-variance features).
#' @param pre_processing_option List of options for the pre-processing,
#'   passed to \code{preProcOptions} in \code{caret::trainControl()}.
#'   Default \code{list(uniqueCut = 5)}.
#'
#' @details
#' Features are filtered in \code{create_caret_df()} by mean abundance
#' (> 5e-5) and prevalence (detected in at least 5 samples in at least
#' 2 studies). Both models are tuned on the same leave-one-study-out splits,
#' using R² and the one-standard-error rule to select the final parameters.
#' The random forest uses the custom caret method \code{adapt_ranger} with
#' permutation importance and the grid:
#' \itemize{
#'   \item \code{mtry}: 3, 10
#'   \item \code{min.node.size}: 5, 10
#'   \item \code{num.trees}: 200; \code{splitrule}: variance;
#'     \code{replace}: \code{FALSE}
#' }
#' The training data is stored in the random forest model, so that
#' \code{\link{get_shap_long}} can use it by default.
#'
#' Training is random; call \code{set.seed()} beforehand for reproducible
#' results. \code{allowParallel = TRUE} only takes effect if a parallel
#' backend is registered.
#'
#' @return A \code{caretList} with the trained models \code{lasso} and
#'   \code{rf1}.
#'
#' @seealso \code{\link{create_caret_df}}, \code{\link{get_shap_long}},
#'   \code{\link{get_predictions_lifestyle}}
#'
#' @importFrom caret findCorrelation trainControl
#' @importFrom caretEnsemble caretList caretModelSpec
#' @importFrom rsample group_vfold_cv rsample2caret
#' @importFrom stats cor
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' models <- get_final_model(ps_train, prefilter = TRUE)
#' models$rf1$bestTune
#' }
#'
#' @export
get_final_model <- function(ps, transform = "compositional",
                            extra_cols = c("Observed", "Shannon"),
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
    family_cleaned_cor <- cor(for_caret_list_train$features, use = "pairwise.complete.obs")
    highlyCorrelated_family_cleaned <- findCorrelation(family_cleaned_cor, cutoff=0.8)
    if (length(highlyCorrelated_family_cleaned) > 0){
      for_caret_list_train$features <- for_caret_list_train$features[,-highlyCorrelated_family_cleaned]
    }
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


#' Train age prediction models for leave-one-dataset-out validation
#'
#' Prepares the features of a phyloseq object with
#' \code{\link{create_caret_df}} and trains a lasso regression and a random
#' forest to predict age. Hyperparameters are tuned with grouped
#' cross-validation in which each study forms one fold. Intended as the model
#' fitting step within leave-one-dataset-out cross-validation, where
#' \code{ps} contains all studies except the held-out one.
#'
#' @param ps A \code{phyloseq} object with the training samples; see
#'   \code{\link{create_caret_df}} for the required sample data.
#' @param transform Character. Abundance transformation passed to
#'   \code{create_caret_df()}. Default \code{"compositional"}.
#' @param extra_cols Character vector of additional predictors passed to
#'   \code{create_caret_df()}. Default richness, Shannon diversity and
#'   industrialized lifestyle.
#' @param prefilter Logical. If \code{TRUE}, removes highly correlated
#'   features (absolute Pearson correlation > 0.8) with
#'   \code{caret::findCorrelation()} before training. Default \code{FALSE}.
#' @param rf_depth Integer. Rarefaction depth for diversity, passed to
#'   \code{create_caret_df()}. Default 2000.
#' @param pre_processing Character vector of caret pre-processing steps.
#'   Default \code{"nzv"} (remove near-zero-variance features).
#' @param pre_processing_option List of options for the pre-processing,
#'   passed to \code{preProcOptions} in \code{caret::trainControl()}.
#'   Default \code{list(uniqueCut = 5)}.
#'
#' @details
#' Features are filtered in \code{create_caret_df()} by mean abundance
#' (> 5e-5) and prevalence (detected in at least 5 samples in at least
#' 2 studies). Both models are tuned on the same leave-one-study-out splits,
#' using R² and the one-standard-error rule to select the final parameters.
#' The random forest uses the custom caret method \code{adapt_ranger} with
#' permutation importance and the grid:
#' \itemize{
#'   \item \code{mtry}: 3, 10
#'   \item \code{min.node.size}: 5, 10
#'   \item \code{num.trees}: 200; \code{splitrule}: variance;
#'     \code{replace}: \code{FALSE}
#' }
#' The training data is stored in the random forest model, so that
#' \code{\link{get_shap_long}} can use it by default.
#'
#' Training is random; call \code{set.seed()} beforehand for reproducible
#' results. \code{allowParallel = TRUE} only takes effect if a parallel
#' backend is registered.
#'
#' @return A \code{caretList} with the trained models \code{lasso} and
#'   \code{rf1}.
#'
#' @seealso \code{\link{create_caret_df}}, \code{\link{get_shap_long}},
#'   \code{\link{get_predictions_lifestyle}}
#'
#' @importFrom caret findCorrelation trainControl
#' @importFrom caretEnsemble caretList caretModelSpec
#' @importFrom rsample group_vfold_cv rsample2caret
#' @importFrom stats cor
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' ps_train <- subset_samples(ps, study != "backhed_2015")
#' models <- get_final_model(ps_train, prefilter = TRUE)
#' models$rf1$bestTune
#' }
#'
#' @export
get_predictions_nested_cv <- function(test_set, ps, prefix = "default",
                                      transform = "compositional",
                                      extra_cols = c("Observed", "Shannon", "lifestyle_industrialized"),
                                      prefilter = F,
                                      pre_processing = "nzv",
                                      rf_depth = 2000,
                                      n_cores = 10,
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
    family_cleaned_cor <- cor(for_caret_list_train$features, use = "pairwise.complete.obs")
    highlyCorrelated_family_cleaned <- findCorrelation(family_cleaned_cor, cutoff=0.8)
    if (length(highlyCorrelated_family_cleaned) > 0) {
      for_caret_list_train$features <- for_caret_list_train$features[,-highlyCorrelated_family_cleaned]
    }
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
                file = paste0("../data/nested_cv_dataset_models/",
                              prefix, "_", test_set, ".rds"))
      predictions <- cbind(for_caret_list_test$metadata,
                           predict(rf_model_list, newdata = for_caret_list_test$features))
      # Boruta feature importance
      print("Run Boruta")
      sel_feat <- Boruta(x=for_caret_list_train$features, y = for_caret_list_train$metadata$age,
                         num.trees = rf_model_list$rf1$bestTune$num.trees, 
                         mtry = rf_model_list$rf1$bestTune$mtry, num.threads = ceiling(n_cores/5))
      write_rds(sel_feat,
                file = paste0("../data/nested_cv_dataset_models/",
                              prefix, "_", test_set, "_boruta_res.rds"))
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
