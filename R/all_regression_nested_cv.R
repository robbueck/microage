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
library(lme4)
library(rstatix)
library(ggprism)

source("./alt_models.R")
source("./regression_functions.R")


# Switches ################################
# family level
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

cl <- makePSOCKcluster(ceiling(n_cores/2))
registerDoParallel(cl)
set.seed(825)
future::plan(multisession, workers = ceiling(n_cores/2))


# family data without lifestyle ################################################
ps_object_family_raw <- readRDS("../data/all_phyloseq_rf_filter_family.rds") %>%
  subset_samples(., age <= 730 & age > 1) 

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
  save(family_no_ls_nested_cv_preds,
       file = "../data/all_nested_cv_family_no_ls.RData")
} else {
  load("../data/all_nested_cv_family_no_ls.RData")
}
toc()

if(final_model_family_no_ls_step){
  final_model_family_no_ls <- get_final_model(ps = ps_object_family_raw, extra_cols = c("Observed", "Shannon"))
  shap_out_family_no_ls <- get_shap_long(final_model_family_no_ls$rf1, features = "important")
  save(final_model_family_no_ls, shap_out_family_no_ls, file = "../data/final_no_ls_family.RData")
} else {
  load("../data/final_no_ls_family.RData")
}


## family data for industrialized populations only specific model ##################

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
  save(family_industrialized_nested_cv_preds,
       file = "../data/all_nested_cv_family_industrialized.RData")
} else {
  load("../data/all_nested_cv_family_industrialized.RData")
}
toc()


if(final_model_family_industrialized_step){
  final_model_family_industrialized <- get_final_model(ps = ps_object_family_raw_industrialized, extra_cols = c("Observed", "Shannon"))
  shap_out_family_industrialized <- get_shap_long(final_model_family_industrialized$rf1, features = "important")
  save(final_model_family_industrialized, shap_out_family_industrialized,
       file = "../data/final_industrialized_family.RData")
} else {
  load("../data/final_industrialized_family.RData")
}

## family data for non industrialized populations only #############################

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
  save(family_nonindustrialized_nested_cv_preds,
       file = "../data/all_nested_cv_family_nonindustrialized.RData")
} else {
  load("../data/all_nested_cv_family_nonindustrialized.RData")
}
toc()

if(final_model_family_nonindustrialized_step){
  final_model_family_nonindustrialized <- get_final_model(ps = ps_object_family_raw_nonindustrialized, extra_cols = c("Observed", "Shannon"))
  shap_out_family_nonindustrialized <- get_shap_long(final_model_family_nonindustrialized$rf1, features = "important")
  save(final_model_family_nonindustrialized, shap_out_family_nonindustrialized, 
       file = "../data/final_nonindustrialized_family.RData")
} else {
  load("../data/final_nonindustrialized_family.RData")
}

# genus data without lifestyle #################################################

print("genus data")
ps_object_genus_raw <- readRDS("../data/all_phyloseq_rf_filter_genus.rds") %>%
  subset_samples(., age <= 730 & age > 1)

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
  save(genus_no_ls_nested_cv_preds,
       file = "../data/all_nested_cv_genus_no_ls.RData")
} else {
  load("../data/all_nested_cv_genus_no_ls.RData")
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
  save(final_model_genus_no_ls,
       file = "../data/final_no_ls_genus.RData")
  save(shap_out_genus_no_ls,
       file = "../large_files/shap_no_ls_genus.RData")
} else {
  load("../data/final_no_ls_genus.RData")
  load("../large_files/shap_no_ls_genus.RData")
}



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


### effect of delivery mode:
bm_data <- genus_no_ls_nested_cv_preds %>% 
  filter(!is.na(birthmode), lifestyle == "industrialized") %>%
  mutate(c_section = birthmode == "c") %>%
  select(age, c_section, study, rf1)
ggplot(bm_data, aes(x = age, y = (rf1), color = c_section)) +
  geom_point(size = 0.2) +
  geom_smooth()

full_model_lin_bm <- lmer(rank(rf1) ~ age + c_section + (1|study),
                       data = bm_data)
# is the intercept different?
intercept_model_lin_bm <- lmer(rank(rf1) ~ age + (1|study),
                            data = bm_data)
1 - deviance(full_model_lin_bm) / deviance(intercept_model_lin_bm)
test_bm <- lrtest(intercept_model_lin_bm, full_model_lin_bm)



## genus data for industrialized populations only specific model ##################

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
  save(genus_industrialized_nested_cv_preds,
       file = "../data/all_nested_cv_genus_industrialized.RData")
} else {
  load("../data/all_nested_cv_genus_industrialized.RData")
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
  save(final_model_genus_industrialized,
       file = "../data/final_industrialized_genus.RData")
  save(shap_out_genus_industrialized,
       file = "../large_files/shap_industrialized_genus.RData")
} else {
  load("../data/final_industrialized_genus.RData")
  load("../large_files/shap_industrialized_genus.RData")
}


## genus data for non industrialized populations only #############################

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
  save(genus_nonindustrialized_nested_cv_preds,
       file = "../data/all_nested_cv_genus_nonindustrialized.RData")
} else {
  load("../data/all_nested_cv_genus_nonindustrialized.RData")
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
  save(final_model_genus_nonindustrialized, shap_out_genus_nonindustrialized,
       file = "../data/final_nonindustrialized_genus.RData")
} else {
  load("../data/final_nonindustrialized_genus.RData")
}


# genus data without lifestyle rarefied #################################################
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
       file = "../data/all_nested_cv_genus_no_ls_rare_h.RData")
} else {
  load("../data/all_nested_cv_genus_no_ls_rare_h.RData")
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
  save(final_model_genus_no_ls_rare_h, shap_out_genus_no_ls_rare_h, 
       file = "../data/final_no_ls_genus_rare_h.RData")
} else {
  load("../data/final_no_ls_genus_rare_h.RData")
}

## genus data for industrialized populations only specific model ##################

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
       file = "../data/all_nested_cv_genus_industrialized_rare_h.RData")
} else {
  load("../data/all_nested_cv_genus_industrialized_rare_h.RData")
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
  save(final_model_genus_industrialized_rare_h, shap_out_genus_industrialized_rare_h,
       file = "../data/final_industrialized_genus_rare_h.RData")
} else {
  load("../data/final_industrialized_genus_rare_h.RData")
}

## genus data for non industrialized populations only #############################

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
       file = "../data/all_nested_cv_genus_nonindustrialized_rare_h.RData")
} else {
  load("../data/all_nested_cv_genus_nonindustrialized_rare_h.RData")
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
  save(final_model_genus_nonindustrialized_rare_h, shap_out_genus_nonindustrialized_rare_h,
       file = "../data/final_nonindustrialized_genus_rare_h.RData")
} else {
  load("../data/final_nonindustrialized_genus_rare_h.RData")
}