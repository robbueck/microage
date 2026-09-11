# create model on industrialized or non-industrialized studies and predict age in the respective other set of studies
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
library(Boruta)
library(lme4)
library(lmtest)
# library(MatchIt)
source("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/alt_models.R")
source("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/regression_functions.R")
source("/fast/AG_Forslund/rob/mm_index/R_scripts/setlists.R")
source("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/lifestyle_regression_functions.R")

# Switches ###############
dataset_nested_cv_genus_step <- T
dataset_nested_cv_family_step <- F
dataset_nested_cv_genus_2months_bins_step <- F
shap_step <- F
ks_step <- F


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


# genus data ###################################################################

ps_object_genus_raw <- readRDS("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_phyloseq_rf_filter_genus.rds") %>%
  subset_samples(age <= 730 & age > 1)
metadata <- ps_object_genus_raw %>%
  sample_data() %>%
  data.frame

print("genus data industrialized datasets")
tic()
if(dataset_nested_cv_genus_step){
  lifestyle_predictions_genus <- c("industrialized", "non_industrialized") %>% 
    future_map_dfr(~ get_predictions_lifestyle(. ,ps = ps_object_genus_raw))
  save(lifestyle_predictions_genus, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/inter_lifestyle_pred_genus.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/inter_lifestyle_pred_genus.RData")
}
toc()

# compare with prediction from same lifestyle:
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_industrialized.RData")
genus_industrialized_nested_cv_preds_long <- genus_industrialized_nested_cv_preds %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "industrialized")

industrialized_lm_list_genus <- get_lm_list(pred_df = genus_industrialized_nested_cv_preds_long,
                                     grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "industrialized")
library(gamlss)
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_nonindustrialized.RData")
genus_nonindustrialized_nested_cv_preds_long <- genus_nonindustrialized_nested_cv_preds %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "non_industrialized")
nonindustrialized_lm_list_genus <- get_lm_list(pred_df = genus_nonindustrialized_nested_cv_preds_long, grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "non_industrialized")


# default model:
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_no_ls.RData")
genus_nested_cv_preds_long <- genus_no_ls_nested_cv_preds %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "combined")
all_lm_list_genus <- get_lm_list(pred_df = genus_nested_cv_preds_long, grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "combined")


lifestyle_predictions_genus <- lifestyle_predictions_genus %>%
  mutate(model_name = "rf1") %>%
  mutate(training_set = ifelse(lifestyle == "industrialized", yes = "non_industrialized", no = "industrialized"))
lifestyle_lm_list_genus <- get_lm_list(pred_df = lifestyle_predictions_genus, grouping = c("study", "lifestyle", "model_name")) %>%
  mutate(training_set = ifelse(lifestyle == "industrialized", yes = "non_industrialized", no = "industrialized")) %>%
  rbind(., industrialized_lm_list_genus, nonindustrialized_lm_list_genus, all_lm_list_genus)

lifestyle_lm_list_genus <- lifestyle_lm_list_genus %>%
  mutate(training_set = factor(training_set, levels = c("non_industrialized", "combined", "industrialized")),
         lifestyle = factor(lifestyle, levels = c("industrialized", "non_industrialized")))

# what is the effect of sample_sum on the model?
# lifestyle_lm_list_genus %>% pivot_longer(cols = c("R2", "R2_sample_sum"), names_to = "model", values_to = "R2") %>%
#   ggplot(.,aes(x = lifestyle, y = R2, color = training_set, fill = model)) +
#   geom_boxplot(alpha = 0.5) +
#   xlab("Test set") +
#   geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
#   theme_classic()
# lifestyle_lm_list_genus %>% pivot_longer(cols = c("R2", "R2_sample_sum"), names_to = "model", values_to = "R2") %>%
#   rstatix::group_by(lifestyle, training_set) %>%
#   rstatix::wilcox_test(R2 ~ model, paired = T) %>%
#   rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
#   rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1))
  

# genus_industrialized_nested_cv_preds_long
# genus_nonindustrialized_nested_cv_preds_long
# genus_nested_cv_preds_long
# lifestyle_predictions_genus
# genus_nested_cv_preds_long_combined <- bind_rows(genus_industrialized_nested_cv_preds_long, 
#                                                  genus_nonindustrialized_nested_cv_preds_long) %>%
#   bind_rows(.,genus_nested_cv_preds_long) %>%
#   bind_rows(.,lifestyle_predictions_genus)
# 
# test_combined <- genus_nested_cv_preds_long_combined %>%
#   mutate(study = factor(study),
#          subject_ID = factor(subject_ID)) %>% # if I add subject ID as random factor, all significance disappears
#   filter(!study %in% c("davis_2017", "bender_2016")) %>%
#   group_by(lifestyle, training_set, study) %>%
#   summarize(R2 = Rsq(gamlss(pred ~ pb(age) + random(subject_ID), family = NO)))
# 
# genus_nested_cv_preds_long_combined %>% group_by(study) %>% summarize(age_dist = max(age) - min(age),
#                                                                       max_age = max(age)) 
# colSums(table(genus_nested_cv_preds_long_combined$age, genus_nested_cv_preds_long_combined$study) > 0)
# test <-  genus_nonindustrialized_nested_cv_preds_long %>% filter(study == "raman_2019")
# test_model <- gamlss(pred ~ lo(~age), data = test %>% select(age, pred))
# test_non_industrialized <- gamlss(pred ~ pb(age) + random(study), data = genus_nonindustrialized_nested_cv_preds_long %>% 
#                           select(age, pred, study) %>%
#                           mutate(study = factor(study)), , control = gamlss.control(n.cyc = 100))
# test_industrialized <- gamlss(pred ~ pb(age) + random(study), data = genus_industrialized_nested_cv_preds_long %>% 
#                       select(age, pred, study) %>%
#                       mutate(study = factor(study)), , control = gamlss.control(n.cyc = 100))


# stat test (stay with R2, RSME performs worse)
df_p_val_lifestyle_genus <- lifestyle_lm_list_genus %>%
  # lifestyle_lm_list_genus %>%
  filter(study != "bender_2016") %>%
  arrange(study) %>%
  mutate(training_set = factor(training_set)) %>%
  rstatix::group_by(lifestyle) %>%
  rstatix::wilcox_test(R2 ~ training_set, paired = T) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>% 
  rstatix::add_xy_position(x = "lifestyle", dodge = 0.8) 

ls_performance <- ggplot(lifestyle_lm_list_genus, aes(x=lifestyle, y = R2)) +
  geom_boxplot(aes(fill = training_set)) +
  # ylim(0, 1) +
  xlab("Test set") +
  add_pvalue(df_p_val_lifestyle_genus,
             label = "{p.adj.signif}",
             # step.group.by = "variation",
             step.increase = 0.05,
             tip.length = 0.01,
             # bracket.nudge.y = 0.02,
             xmin = "xmin", 
             xmax = "xmax",
             show.legend = FALSE) +
  ylim(0,NA) +
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
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/lifestyle_perfomance_genus.pdf", 
       width = 10, plot = ls_performance)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/lifestyle_perfomance_genus.png", 
       width = 10, plot = ls_performance)
save(lifestyle_lm_list_genus, ls_performance,
     file = "/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/lifestyle_perfomance_genus.RData")

# family model #############################################
ps_object_family_raw <- readRDS("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_phyloseq_rf_filter_family.rds") %>%
  subset_samples(., age <= 730 & age > 1) 

print("family data industrialized datasets")
tic()
if(dataset_nested_cv_family_step){
  lifestyle_predictions_family <- c("industrialized", "non_industrialized") %>% 
    future_map_dfr(~ get_predictions_lifestyle(. ,ps = ps_object_family_raw))
  save(lifestyle_predictions_family, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/inter_lifestyle_pred_family.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/inter_lifestyle_pred_family.RData")
}
toc()
# compare with prediction from same lifestyle:
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_family_industrialized.RData")
family_industrialized_nested_cv_preds_long <- family_industrialized_nested_cv_preds %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "industrialized")

industrialized_lm_list_family <- get_lm_list(pred_df = family_industrialized_nested_cv_preds_long,
                                            grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "industrialized")

load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_family_nonindustrialized.RData")
family_nonindustrialized_nested_cv_preds_long <- family_nonindustrialized_nested_cv_preds %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "non_industrialized")
nonindustrialized_lm_list_family <- get_lm_list(pred_df = family_nonindustrialized_nested_cv_preds_long, 
                                               grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "non_industrialized")

# load combined model:
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_family_no_ls.RData")
family_nested_cv_preds_long <- family_no_ls_nested_cv_preds %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "combined")
all_lm_list_family <- get_lm_list(pred_df = family_nested_cv_preds_long, 
                                 grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "combined")

# merge model outputs
lifestyle_predictions_family <- lifestyle_predictions_family %>%
  mutate(model_name = "rf1") %>%
  mutate(training_set = ifelse(lifestyle == "industrialized", yes = "non_industrialized", no = "industrialized"))
lifestyle_lm_list_family <- get_lm_list(pred_df = lifestyle_predictions_family, grouping = c("study", "lifestyle", "model_name")) %>%
  mutate(training_set = ifelse(lifestyle == "industrialized", yes = "non_industrialized", no = "industrialized")) %>%
  rbind(., industrialized_lm_list_family, nonindustrialized_lm_list_family, all_lm_list_family)

lifestyle_lm_list_family <- lifestyle_lm_list_family %>%
  mutate(training_set = factor(training_set, levels = c("non_industrialized", "combined", "industrialized")),
         lifestyle = factor(lifestyle, levels = c("industrialized", "non_industrialized")))

# stat test (stay with R2, RSME performs worse)
df_p_val_lifestyle_family <- lifestyle_lm_list_family %>%
  # lifestyle_lm_list_genus %>%
  filter(study != "bender_2016") %>%
  arrange(study) %>%
  mutate(training_set = factor(training_set)) %>%
  rstatix::group_by(lifestyle) %>%
  rstatix::wilcox_test(R2 ~ training_set, paired = T) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>% 
  rstatix::add_xy_position(x = "lifestyle", dodge = 0.8) 

ls_performance_family <- ggplot(lifestyle_lm_list_family, aes(x=lifestyle, y = R2)) +
  geom_boxplot(aes(fill = training_set)) +
  # ylim(0, 1) +
  xlab("Test set") +
  add_pvalue(df_p_val_lifestyle_family,
             label = "{p.adj.signif}",
             # step.group.by = "variation",
             step.increase = 0.05,
             tip.length = 0.01,
             # bracket.nudge.y = 0.02,
             xmin = "xmin", 
             xmax = "xmax",
             show.legend = FALSE) +
  ylim(0,NA) +
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

save(lifestyle_lm_list_genus, ls_performance, lifestyle_lm_list_family,
     file = "/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/lifestyle_perfomance_genus.RData")


## comparison predictions genus vs family level #############################
combined_family <- bind_rows(family_nested_cv_preds_long, lifestyle_predictions_family,
                             family_nonindustrialized_nested_cv_preds_long, 
                             family_industrialized_nested_cv_preds_long) %>% 
  select(lifestyle, study, training_set, run_accession, age, pred)
combined_genus <- bind_rows(genus_nested_cv_preds_long, lifestyle_predictions_genus,
                           genus_nonindustrialized_nested_cv_preds_long, 
                           genus_industrialized_nested_cv_preds_long) %>%
  select(lifestyle, study, training_set, run_accession, age, pred)
combined_genus_family <- left_join(combined_family, combined_genus,
          by = c("lifestyle", "study", "training_set", "run_accession", "age"),
          suffix = c("_family", "_genus"))
ggplot(combined_genus_family, aes(x = pred_genus, y = pred_family, color = lifestyle)) +
  geom_point(size = 0.3, alpha = 0.2) +
  geom_smooth(method = "lm") +
  facet_wrap(~training_set) +
  coord_equal()
combined_genus_family %>%
  group_by(lifestyle, training_set) %>%
  rstatix::cor_test(pred_family, pred_genus, method = "spearman")

### per study performance #######################
left_join(lifestyle_lm_list_family, lifestyle_lm_list_genus, 
          by = c("study", "lifestyle", "model_name", "training_set"),
          suffix = c("_family", "_genus")) %>%
  ggplot(., aes(x = R2_genus, y = R2_family, color = lifestyle)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(~training_set) +
  coord_equal()

# 2 months binned model #############################
print("genus data 2 months bins")
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
intervals <- lapply(seq(0, 660, 60), \(x) c(x, min(x + 60, 720)))
tic()

if(dataset_nested_cv_genus_2months_bins_step){
  lifestyle_predictions_genus_2months_bins_list <- lapply(intervals, 
                                                   function(x) run_in_bins(ps_object_genus_raw,
                                                                           cutoffs = x))
  lifestyle_predictions_genus_2months_bins <- lifestyle_predictions_genus_2months_bins_list %>% bind_rows()
  save(lifestyle_predictions_genus_2months_bins, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/inter_lifestyle_pred_genus_2months_binned.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/inter_lifestyle_pred_genus_2months_binned.RData")
}
toc()

# compare with prediction from same lifestyle:
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_indust_no_ls_2month_bins.RData")
genus_2month_bins_industrialized_nested_cv_preds_long <- genus_indust_no_ls_nested_cv_2month_bins %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "industrialized")

industrialized_lm_list_genus_2months_bins <- get_lm_list(pred_df = genus_2month_bins_industrialized_nested_cv_preds_long,
                                             grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "industrialized")

load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_non_indust_no_ls_2month_bins.RData")
genus_non_indust_no_ls_nested_cv_2month_bins_long <- genus_non_indust_no_ls_nested_cv_2month_bins %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "non_industrialized")
nonindustrialized_lm_list_genus_2months_bins <- get_lm_list(pred_df = genus_non_indust_no_ls_nested_cv_2month_bins_long, 
                                                grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "non_industrialized")

# load combined model:
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_no_ls_2month_bins.RData")
genus_no_ls_nested_cv_2month_bins_long <- genus_no_ls_nested_cv_2month_bins %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "combined")
all_lm_list_genus_2months_bins <- get_lm_list(pred_df = genus_no_ls_nested_cv_2month_bins_long, 
                                  grouping = c("study", "lifestyle", "model_name", "interval")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "combined")
ggplot(genus_no_ls_nested_cv_2month_bins_long %>% filter(age <= 500), 
       aes(x = age, y = pred, group = interval)) +
  geom_point(size = 0.3) +
  geom_smooth(method = "lm") +
  geom_smooth(method = "lm", group = NULL) +
  coord_equal() +
  facet_wrap(~lifestyle)+
  theme_classic()


lifestyle_predictions_family <- lifestyle_predictions_family %>%
  mutate(model_name = "rf1") %>%
  mutate(training_set = ifelse(lifestyle == "industrialized", yes = "non_industrialized", no = "industrialized"))
lifestyle_lm_list_family <- get_lm_list(pred_df = lifestyle_predictions_family, grouping = c("study", "lifestyle", "model_name")) %>%
  mutate(training_set = ifelse(lifestyle == "industrialized", yes = "non_industrialized", no = "industrialized")) %>%
  rbind(., industrialized_lm_list_family, nonindustrialized_lm_list_family, all_lm_list_family)

lifestyle_lm_list_family <- lifestyle_lm_list_family %>%
  mutate(training_set = factor(training_set, levels = c("non_industrialized", "combined", "industrialized")),
         lifestyle = factor(lifestyle, levels = c("industrialized", "non_industrialized")))

# stat test (stay with R2, RSME performs worse)
df_p_val_lifestyle_family <- lifestyle_lm_list_family %>%
  # lifestyle_lm_list_genus %>%
  filter(study != "bender_2016") %>%
  arrange(study) %>%
  mutate(training_set = factor(training_set)) %>%
  rstatix::group_by(lifestyle) %>%
  rstatix::wilcox_test(R2 ~ training_set, paired = T) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>% 
  rstatix::add_xy_position(x = "lifestyle", dodge = 0.8) 

ls_performance_family <- ggplot(lifestyle_lm_list_family, aes(x=lifestyle, y = R2)) +
  geom_boxplot(aes(fill = training_set)) +
  # ylim(0, 1) +
  xlab("Test set") +
  add_pvalue(df_p_val_lifestyle_family,
             label = "{p.adj.signif}",
             # step.group.by = "variation",
             step.increase = 0.05,
             tip.length = 0.01,
             # bracket.nudge.y = 0.02,
             xmin = "xmin", 
             xmax = "xmax",
             show.legend = FALSE) +
  ylim(0,NA) +
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


# tipping point? ###############################
sliding_window_size <- 50
rf1_sliding_window <- map_dfr(seq(0, max(genus_no_ls_nested_cv_preds$age) - sliding_window_size, by = 7),
                              function(x) {
                                meta_df_filt <- genus_no_ls_nested_cv_preds %>%
                                  filter(age < x + sliding_window_size, age >= x) %>% # sliding window
                                  mutate(rf1_r = rank(rf1))
                                list_models <- lme4::lmList(data = meta_df_filt, formula = rf1_r ~ age | lifestyle)
                                model_df <- data.frame(grouping_x_model_name = names(list_models),
                                                       Coefficient = coef(list_models)$age,
                                                       R2 = summary(list_models)$r.squared,
                                                       sliding_window = x + sliding_window_size / 2) %>%
                                  mutate(sd = ifelse(grouping_x_model_name == "industrialized", 
                                                     yes = sd(meta_df_filt$rf1[meta_df_filt$lifestyle == "industrialized"]),
                                                     no = sd(meta_df_filt$rf1[meta_df_filt$lifestyle == "non_industrialized"])))
                                return(model_df)
                              })
ggplot(rf1_sliding_window, aes(x = sliding_window, y = R2, color = grouping_x_model_name)) + geom_point()
ggplot(rf1_sliding_window, aes(x = sliding_window, y = sd, color = grouping_x_model_name)) + geom_point()

rf1_expanding_window <- map_dfr(seq(30, max(genus_no_ls_nested_cv_preds$age), by = 30),
                                function(x) {
                                  meta_df_filt <- genus_no_ls_nested_cv_preds %>%
                                    filter(age > x) %>% # sliding window
                                    mutate(rf1_r = rank(rf1))
                                  list_models <- lme4::lmList(data = meta_df_filt, formula = rf1_r ~ age | lifestyle)
                                  model_df <- data.frame(grouping_x_model_name = names(list_models),
                                                         Coefficient = coef(list_models)$age,
                                                         R2 = summary(list_models)$r.squared,
                                                         expanding_window = x) %>%
                                    mutate(sd = ifelse(grouping_x_model_name == "industrialized", 
                                                       yes = sd(meta_df_filt$rf1[meta_df_filt$lifestyle == "industrialized"]),
                                                       no = sd(meta_df_filt$rf1[meta_df_filt$lifestyle == "non_industrialized"])))
                                  
                                  return(model_df)
                                })
ggplot(rf1_expanding_window, aes(x = expanding_window, y = R2, color = grouping_x_model_name)) + geom_point()
ggplot(rf1_expanding_window, aes(x = expanding_window, y = sd, color = grouping_x_model_name)) + geom_point()



## different feature importance #################################################
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
industrialized_studies <- lifestyle_lm_list_genus %>% 
  filter(lifestyle == "industrialized") %>%
  pull(study) %>% unique
non_industrialized_studies <- lifestyle_lm_list_genus %>% 
  filter(lifestyle == "non_industrialized") %>%
  pull(study) %>% unique

combined_importance_nonindustrialized_genus <- map(non_industrialized_studies, 
                                            function(x) get_feature_importance(x, lst = "nonindustrialized")) %>%
  purrr::reduce(rbind) %>%
  mutate(lifestyle = "non_industrialized") %>%
  group_by(taxon) %>%
  filter(n() > 1) %>% # only take those that occur in more than 2 models
  ungroup()
combined_importance_industrialized_genus <- map(industrialized_studies, function(x) get_feature_importance(x, lst = "industrialized")) %>%
  purrr::reduce(rbind) %>%
  mutate(lifestyle = "industrialized") %>%
  group_by(taxon) %>%
  filter(n() > 1) %>% # only take those that occur in more than 2 models
  ungroup()

bind_rows(combined_importance_nonindustrialized_genus, combined_importance_industrialized_genus) %>%
  group_by(ID, lifestyle) %>%
  dplyr::summarize(count = n()) %>%
  ggplot(., aes( x = lifestyle, y = count)) +
  geom_boxplot() +
  geom_jitter() +
  ylim(0, 100)

data.frame(counts = c(97, 91, 61,45), studies = c(13, 13, 8, 8), threshold = c(1, 3, 1, 3)) %>%
  ggplot(.,aes(x = threshold, y = counts/studies, color = threshold)) +
  geom_point()
## important taxa stats #################
# n important samples 
important_union <- unique(c(combined_importance_nonindustrialized_genus$taxon, combined_importance_industrialized_genus$taxon))
important_industrialized <- unique(combined_importance_industrialized_genus$taxon) 
length(important_industrialized)
important_non_industrialized <- unique(combined_importance_nonindustrialized_genus$taxon) 
length(important_non_industrialized)
# shared:
length(intersect(important_industrialized, important_non_industrialized))
# specific to one:
important_industrialized[!important_industrialized %in% important_non_industrialized] %>% length
important_non_industrialized[!important_non_industrialized %in% important_industrialized] %>% length

# relative abundance of important taxa in data:
for_caret_list <- create_caret_df(ps_object = ps_object_genus_raw, transformation = "compositional",
                                  mean_ab_cutoff = 5e-5, study_prevalence_cutoff = 1,
                                  additional_cols = c("lifestyle", "study"),
                                  only_multiple_samples = F)

# combined importances
important_taxa_long <- for_caret_list$features %>% 
  select(important_union[!important_union %in% c("Observed", "Shannon")], lifestyle) %>%
  rownames_to_column %>%
  pivot_longer(., cols = -c(rowname, lifestyle),
               names_to = "taxon",
               values_to = "value")
important_taxa_long %>% group_by(rowname) %>%
  summarize(rel_abund = sum(value)) %>%
  pull(rel_abund) %>% mean
important_taxa_long %>% group_by(rowname, lifestyle) %>%
  summarize(rel_abund = sum(value)) %>%
  group_by(lifestyle) %>%
  summarize(mean_total = mean(rel_abund))
# only industial important
important_taxa_long %>% 
  filter(taxon %in% important_industrialized[!important_industrialized %in% c("Observed", "Shannon")]) %>%
  group_by(rowname, lifestyle) %>%
  summarize(rel_abund = sum(value)) %>%
  group_by(lifestyle) %>%
  summarize(mean_total = mean(rel_abund))
# only non-industrial important
important_taxa_long %>% 
  filter(taxon %in% important_non_industrialized[!important_non_industrialized %in% c("Observed", "Shannon")]) %>%
  group_by(rowname, lifestyle) %>%
  summarize(rel_abund = sum(value)) %>%
  group_by(lifestyle) %>%
  summarize(mean_total = mean(rel_abund))

# are taxa only important in one lifestyle also present in another lifestyle?
important_taxa_long %>% 
  # filter(lifestyle == "non_industrialized") %>% 
  filter(taxon %in% important_industrialized[!important_industrialized %in% important_non_industrialized]) %>%
  mutate(value = value > 5e-5) %>%
  group_by(taxon, lifestyle) %>%
  summarize(prevalence = mean(value)) %>%
  ggplot(.,aes(x = lifestyle, y = prevalence)) +
  geom_boxplot() +
  geom_jitter()
important_taxa_long %>% 
  filter(lifestyle == "industrialized") %>% 
  filter(taxon %in% important_non_industrialized[!important_non_industrialized %in% important_industrialized]) %>%
  mutate(value = value > 5e-5) %>%
  group_by(taxon) %>%
  summarize(prevalence = sum(value)) %>% pull(prevalence) %>% table


vary_imp_thres <- function(thres) {
  n_noindustrialized<- map(non_industrialized_studies, function(x) get_feature_importance(x, lst = "nonindustrialized")) %>%
    purrr::reduce(rbind) %>% 
    filter(!taxon %in% c("Shannon", "Observed")) %>% 
    pull(taxon) %>% table %>% `>=` (thres) %>% sum
  n_industrialized <-   map(industrialized_studies, function(x) get_feature_importance(x, lst = "industrialized")) %>%
    purrr::reduce(rbind) %>% filter(!taxon %in% c("Shannon", "Observed")) %>% 
    pull(taxon) %>% table %>% `>=` (thres) %>% sum
  data.frame(threshold = thres, industrialized = n_industrialized, non_industrialized = n_noindustrialized)
}
vary_thresh <- lapply(c(1:14), vary_imp_thres) %>%
  bind_rows %>%
  mutate(industrialized = industrialized / 13,
         non_industrialized = non_industrialized / 8) %>%
  pivot_longer(cols = c("industrialized", "non_industrialized"), names_to = "lifestyle",
                             values_to = "n_important_features_per_study")
ggplot(vary_thresh, aes( x = threshold, y = n_important_features_per_study, color = lifestyle)) +
  geom_point()+
  ylim(0, 8) +
  geom_line()


# add prevalence:
prevalence_lifestyle <- for_caret_list$features %>%
  # mutate(study_lifestyle = paste(lifestyle, study.1, sep = "xxx")) %>%
  group_by(lifestyle, study.1) %>% 
  summarize(across(everything(), ~ mean(.>0)), .groups = "drop") %>%  # get prevalence, mean over studies
  select(-study.1) %>%
  group_by(lifestyle) %>%
  summarize(across(everything(), ~ mean(.)), .groups = "drop") %>%
  t() %>%
  data.frame %>%
  `colnames<-`(.[1, ]) %>%
  .[-1, ] %>% 
  rownames_to_column(.,var = "taxon") %>% 
  pivot_longer(., -taxon, names_to = "lifestyle", values_to = "prevalence") %>%
  mutate(prevalence = as.numeric(prevalence) * 100)# %>%

combined_importance_prev_lifestyles_genus <- rbind(combined_importance_industrialized_genus, combined_importance_nonindustrialized_genus) %>%
  complete(taxon, lifestyle) %>%
  mutate(lifestyle = factor(lifestyle, levels = c("industrialized", "non_industrialized")),
         model = lifestyle,
         measure_type = ifelse(taxon %in% c("Shannon", "Observed"), yes = "a_div", no = "taxon")) %>%
  left_join(., prevalence_lifestyle, by = c("taxon", "lifestyle"))

saveRDS(combined_importance_prev_lifestyles_genus,
        "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/genus_nested_cv_important_taxa_lifestyle.rds")

# shorten taxon names
short_taxa_names <- combined_importance_prev_lifestyles_genus %>% select(taxon) %>%
  distinct() %>%
  mutate(display_taxon = taxon %>% gsub(".*unclassified_", "Uncl_ ",.) %>%
           make.unique() %>%
           gsub("X[.]", "", .) %>%
           gsub("[.]", " ", .) %>%
           gsub("  ", " ", .) %>%
           gsub("_", ".", .) %>%
           gsub("Family XIII AD3011 group", "Anaerovoracaceae Fam. XIII AD3011", .) %>%
           gsub("UCG 002", "Oscillospiraceae UCG 002", .))

combined_importance_prev_lifestyles_genus <- combined_importance_prev_lifestyles_genus %>% 
  left_join(.,short_taxa_names, by = "taxon") %>%
  # reorder according to importance
  mutate(display_taxon = reorder(display_taxon, importance, FUN = max),
         display_taxon = factor(display_taxon, levels = unique(c("Shannon", "Observed", levels(display_taxon)))),
         lifestyle = factor(lifestyle, levels = c("non_industrialized", "industrialized"))
  )

importance_prevalence_plot <- ggplot(combined_importance_prev_lifestyles_genus,
       aes(y = display_taxon, x = importance, color = lifestyle)) +
  geom_col(data = combined_importance_prev_lifestyles_genus %>%
             distinct(display_taxon, prevalence, model),
           aes(y = display_taxon, x = prevalence, alpha = 0.5), fill = "grey", color = "white") +
  geom_boxplot() +
  xlim(0,100) +
  xlab("Relative importance / Prevalence") +
  facet_wrap(~ model, nrow = 1, drop = T) +
  theme(axis.text.x = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), #element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "grey",
                                          size = 0.3,
                                          linetype = 2),
        strip.text = element_text(size = 14),
        plot.margin = margin(0.2,0.2,0.2,1, "cm"))
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/genus_nested_cv_importance_lifestyle.pdf",
       height = 17, plot = importance_prevalence_plot)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/genus_nested_cv_importance_lifestyle.png",
       height = 17, plot = importance_prevalence_plot)

combined_importance_prev_lifestyles_genus %>%
  group_by(taxon, model) %>%
  summarize(prevalence = mean(prevalence, na.rm = T),
            importance = mean(importance, na.rm = T)) %>%
  ggplot(., aes(x = prevalence, y = importance, color = model)) +
  geom_point() +
  geom_smooth() +
  xlim(0,100)+
  ggpubr::stat_cor(# aes(colour = NULL),
    method = "spearman",
    cor.coef.name = "spearman",
    size = 6,
    show.legend = F) +
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
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/genus_nested_cv_importance_vs_prev_lifestyle.pdf")
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/genus_nested_cv_importance_vs_prev_lifestyle.png")

save(combined_importance_prev_lifestyles_genus, 
     file = "/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/imp_prevs.RData")
  

## shap analysis genus ################################################################
# combined model
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_no_ls_genus.RData")
# final_model_genus_no_ls
# shap_out_genus_no_ls
# industrialized model
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_industrialized_genus.RData")
# final_model_genus_industrialized
# shap_out_genus_industrialized
# non industrialized model
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_nonindustrialized_genus.RData")
# final_model_genus_nonindustrialized
# shap_out_genus_nonindustrialized
# genus_train_data_non_industrialized

# important taxa:
print("run shap analysis")
if(shap_step){
  shap_model_all <- get_shap_long(final_model_genus_no_ls$rf1,
                                  test_data = genus_train_data$features,
                                  features = unique(combined_importance_prev_lifestyles_genus$taxon))
  print("done shap all")
  shap_model_industrialized <- get_shap_long(final_model_genus_industrialized$rf1, 
                                   test_data = genus_train_data$features,
                                   features = unique(combined_importance_prev_lifestyles_genus$taxon))
  print("done shap industrialized")
  shap_model_noindustrialized <- get_shap_long(final_model_genus_nonindustrialized$rf1,
                                     test_data = genus_train_data$features,
                                     features = unique(combined_importance_prev_lifestyles_genus$taxon))
  print("done shap non industrialized")
  save(shap_model_all, shap_model_industrialized, shap_model_noindustrialized, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/lifestyle_shap_data.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/lifestyle_shap_data.RData")
}

## shap vs abundances #################

shap_model_all %>%
  left_join(short_taxa_names, by = "taxon") %>%
  left_join(., metadata, by = "run_accession") %>%
  group_by(display_taxon, lifestyle) %>% 
  summarise(R2 = cor(ab_value, shap_value)) %>%
  # filter(!display_taxon %in% c("Observed", "Shannon")) %>%
  ggplot(., aes(y = display_taxon, x = R2, fill = lifestyle, group = lifestyle)) +
  geom_bar(position="dodge", stat="identity") +
  theme_classic() +
  theme(panel.grid.major.y = element_line(color = "grey",
                                          size = 0.3,
                                          linetype = 2))+
  ggtitle("Combined_model")

merged_shap_models <- rbind(shap_model_industrialized %>% mutate(model = "industrialized"), 
                           shap_model_noindustrialized %>% mutate(model = "non_industrialized")) %>%
  left_join(., short_taxa_names, by = "taxon") %>%
  left_join(., metadata, by = "run_accession")

shap_cors <- merged_shap_models %>%
  group_by(display_taxon, model) %>%  
  # lifestyle  maybe only check for all data combined, because it's more a property of the model than of the data?
  summarise(cor_test = list(cor.test(ab_value, shap_value, method = "spearman", use = "pairwise.complete.obs")),
            display_taxon = unique(display_taxon),
            .groups = "drop") %>%
  mutate(R2 = map_dbl(cor_test, ~ .x$estimate),
         p_value = map_dbl(cor_test, ~ .x$p.value),
         p_adj = p.adjust(p_value, method = "bonferroni"),
         p_adj = ifelse(p_adj == 0,
                         yes = min(p_adj[p_adj!=0], na.rm = T), # remove 0 p-values
                         no = p_adj),
         # keep the order of the first figure
         display_taxon = factor(display_taxon, 
                                levels = levels(combined_importance_prev_lifestyles_genus$display_taxon)))

sha_value_cor <- shap_cors %>%
  ggplot(., aes(y = display_taxon, x = R2, fill = model)) + # fill = lifestyle, group = lifestyle
  geom_bar(position="dodge", stat="identity") +
  facet_wrap(~ model, nrow = 1, drop = T) +
  xlab("Spearman(feature abund, SHAP)") +
  ylab("") +
  theme(axis.text.x = element_text(size = 14),
        # axis.title.x = element_blank(),
        axis.title.y = element_blank(), #element_text(size = 14),
        axis.text.y = element_blank(), #element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        panel.grid.major.y = element_line(color = "grey",
                                          size = 0.3,
                                          linetype = 2),
        strip.text = element_text(size = 14),
        plot.margin = margin(0.2,0.2,0.2,1, "cm"))

ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_lifestyle_models_cor.pdf",
       plot = sha_value_cor, height = 15, width = 6)#, width = 9, height = 10)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_lifestyle_models_cor.png",
       plot = sha_value_cor, height = 15, width = 6)#, width = 9, height = 10)

# combine with importance plot:
library(ggpubr)
ggarrange(importance_prevalence_plot,sha_value_cor)


# correlation between shap and abund for each taxon in different models based on different data
sha_value_cor_point <- shap_cors %>%
  mutate(R2 = case_when(is.na(R2) ~ 0,
                        p_adj > 0.05 ~ 0,
                        .default = R2)) %>%
  select(-cor_test, -p_value, -p_adj) %>%
  pivot_wider(names_from = "model", values_from = R2) %>%
  filter(industrialized != 0 | non_industrialized != 0) %>%
  # filter(!taxon %in% c("X.Eubacterium..coprostanoligenes.group", "UCG.002", "Observed",
  #                      "Shannon", "Lachnospiraceae.ND3007.group", "Intestinibacter",
  #                      "Erysipelatoclostridium", "Family.XIII.AD3011.group",
  #                      "Anaerostipes", "Oscillibacter", "Catenibacterium", "uncl_Firmicutes",
  #                      "X.Eubacterium..eligens.group", "Erysipelotrichaceae.UCG.003")) %>%
  ggplot(., aes(x = industrialized, y = non_industrialized)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black") +
  # Set axis limits and expand to ensure 0 is included
  # scale_x_continuous(expand = c(0, 0)) +
  # scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size = 14),
        plot.margin = margin(0.2,0.2,0.2,1, "cm"))
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_lifestyle_models_cor_point.pdf",
       plot = sha_value_cor_point)#, width = 9, height = 10)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_lifestyle_models_cor_point.png",
       plot = sha_value_cor_point)#, width = 9, height = 10)

### mixed models #################
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

if(ks_step) {
  model_res_shap_ab <- lapply(unique(merged_shap_models$taxon) %>% sort, 
                              function(x) ks_test_lifestyle(tx = x, 
                                                            df = merged_shap_models))
  taxa_combs <- combn(unique(merged_shap_models$taxon), m = 2, simplify = F)
  options(future.globals.maxSize = 8000 * 1024^2)
  plan(multisession, workers = 8)
  test_ks_industrialized <- map_dfr(taxa_combs, 
                                 function(x) ks_test_per_model(merged_shap_models %>% 
                                                                 filter(model == "industrialized"), x))
  test_ks_noindustrialized <- map_dfr(taxa_combs, 
                                   function(x) ks_test_per_model(merged_shap_models %>% 
                                                                   filter(model == "non_industrialized"), x))
  save(model_res_shap_ab, test_ks_industrialized, test_ks_noindustrialized, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/ks_test_shap_res.RData")
} else {load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/ks_test_shap_res.RData")}

between_taxa_comp <- bind_rows(test_ks_industrialized %>% mutate(model = "industrialized",
                                                       comparison = "between_taxa"),
                               test_ks_noindustrialized %>%
                                 mutate(model = "non_industrialized",
                                        comparison = "between_taxa"))
percentile <- quantile(between_taxa_comp$D_ks, 0.05)

# check empirical p-values:
get_p_value <- function(obs_value, null_dist) {
  # Calculate the proportion of null values greater than or equal to the observed value
  p_value <- mean(abs(null_dist) >= abs(obs_value))
  return(p_value)
}

model_res_shap_ab_df <- bind_rows(model_res_shap_ab) %>%
  rstatix::adjust_pvalue(p.col = "p_ks", output.col = "q_ks", method = "fdr") %>%
  mutate(comparison = "between_models",
         empirical_pval = mean(abs(between_taxa_comp$D_ks) >= abs(D_ks)))

ks_res_all <- bind_rows(model_res_shap_ab_df,between_taxa_comp)

ggplot(ks_res_all, aes(x = D_ks, color = comparison)) +
  geom_density() +
  # geom_vline(xintercept = c(quantile(test_ks_industrialized$D_ks, 0.05),
  #                           mean(test_ks_industrialized$D_ks),
  #                           quantile(test_ks_industrialized$D_ks, 0.95)), color = "blue") +
  # geom_vline(xintercept = c(quantile(test_ks_noindustrialized$D_ks, 0.05),
  #                           mean(test_ks_noindustrialized$D_ks),
  #                           quantile(test_ks_noindustrialized$D_ks, 0.95)), color = "red") +
  # geom_vline(xintercept = c(quantile(model_res_shap_ab_df$D_ks, 0.05),
  #                           mean(model_res_shap_ab_df$D_ks),
  #                           quantile(model_res_shap_ab_df$D_ks, 0.95)), color = "grey") +
  theme_classic()
save(ks_res_all, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/ks_test_res.RData")

# combine shap and feature importance:
shap_and_imps <- bind_rows(combined_importance_prev_lifestyles_genus %>% mutate(type = "importance"),
                           shap_cors %>% mutate(type = "shap") %>% filter(!is.na(R2))) %>%
  mutate(model = paste0(model, "_model"))
facet_labels <- shap_and_imps %>%
  distinct(model) %>%
  pull(model) %>%
  setNames(., .)  # Keep only the `model` as the facet label

shap_imp_prev_combined_plot <- ggplot() +
  geom_col(data = shap_and_imps %>%
             filter(type == "importance") %>%
             distinct(display_taxon, prevalence, model, type),
           aes(y = display_taxon, x = prevalence), alpha = 0.5, fill = "grey", color = "white") +
  geom_boxplot(data = shap_and_imps %>%
                 filter(type =="importance"),
               aes(color = lifestyle, y = display_taxon, x = importance)) +
  geom_bar(data = shap_and_imps %>%
             filter(type == "shap"),
           aes(x = R2, y = display_taxon, fill = model),
           position="dodge", stat="identity") +
  ylab("") +
  facet_wrap(~type + model, drop = TRUE, nrow = 1, scales = "free_x") +
  # facet_grid(~model, drop model= TRUE, scales = "free_x") +
  # facet_grid(cols = vars(type, model), drop = TRUE, scales = "free_x") +
  theme(axis.text.x = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        # legend.position = "none",
        legend.position= c(0.2, 0.2),
        legend.title = element_text(size=10),
        legend.text = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "grey",
                                          size = 0.3,
                                          linetype = 2),
        strip.text = element_text(size = 14),
        strip.placement = "outside",
        plot.margin = margin(0.2,0.2,0.2,1, "cm"))

save(combined_importance_prev_lifestyles_genus, shap_cors, short_taxa_names, shap_model_industrialized, shap_model_noindustrialized,
     file = "/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/imp_shap_prev.RData")



# plots for indiviidual taxa
rbind(shap_model_industrialized, shap_model_noindustrialized) %>%
  filter(taxon %in% c("Faecalibacterium", "Bifidobacterium")) %>%
  left_join(., metadata, by = "run_accession") %>%
  group_by(taxon) %>%
  summarise(mean = mean(ab_value),
            std = sd(ab_value),
            abs_max = max(abs(ab_value))) %>%
  left_join(left_join(rbind(shap_model_industrialized, shap_model_noindustrialized) %>%
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
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_correlation.pdf", width = 16, height = 8)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_correlation.png",
       width = 9, height = 10)



rbind(shap_model_industrialized, shap_model_noindustrialized) %>% 
  filter(taxon == c("unclassified_Enterobacteriaceae")) %>%
  left_join(., metadata, by = "run_accession") %>%
  group_by(taxon) %>%
  summarise(mean = mean(ab_value),
            std = sd(ab_value),
            abs_max = max(abs(ab_value))) %>%
  left_join(left_join(rbind(shap_model_industrialized, shap_model_noindustrialized) %>%
                        filter(taxon == c("unclassified_Enterobacteriaceae")), metadata, by = "run_accession"), ., by = c("taxon")) %>%
  mutate(feature_value_z = (ab_value - mean) / std) %>%
  # filter(ab_value > 5e-5) %>%
  # filter(!taxon %in% c("")) %>%
  ggplot(., aes(y = feature_value_z, x = shap_value)) +
  geom_point(alpha = 0.5, size = .4) +# geom_bin2d(aes(fill = ..ndensity..), bins = 100)+
  # scale_fill_continuous(type = "viridis")+
  # geom_violin() +
  theme_classic() +
  geom_smooth(method = "lm") +
  facet_wrap(~lifestyle , scales = "free") +
  theme(panel.grid.major.y = element_line(color = "grey",
                                          size = 0.3,
                                          linetype = 2),
        strip.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16))
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_weird_correlation.pdf", width = 16, height = 8)




# rarefied data ###############################
if(dataset_nested_cv_genus_step){
  lifestyle_predictions_genus_rare_h <- c("industrialized", "non_industrialized") %>% 
    future_map_dfr(~ get_predictions_lifestyle(. ,ps = ps_object_genus_raw %>% rarefy_even_depth(., sample.size = 2000)))
  save(lifestyle_predictions_genus_rare_h, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/inter_lifestyle_pred_genus_rare_h.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/inter_lifestyle_pred_genus_rare_h.RData")
}
if(dataset_nested_cv_genus_step){
  lifestyle_predictions_genus_rare_l <- c("industrialized", "non_industrialized") %>% 
    future_map_dfr(~ get_predictions_lifestyle(. ,ps = ps_object_genus_raw %>% rarefy_even_depth(., sample.size = 10000)))
  save(lifestyle_predictions_genus_rare_l, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/inter_lifestyle_pred_genus_rare_l.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/inter_lifestyle_pred_genus_rare_l.RData")
}
# compare predictinos from rarefied and un-rarefied data:
lifestyle_predictions_genus %>% left_join(lifestyle_predictions_genus_rare_h %>% select(run_accession, pred),
                                          suffix = c("", "_rare_h"), by = "run_accession") %>%
  left_join(lifestyle_predictions_genus_rare_l %>% select(run_accession, pred), suffix = c("", "_rare_l"), by = "run_accession") %>%
  pivot_longer(cols = c("pred_rare_l", "pred_rare_h"), names_to = "rarefaction_type", values_to = "pred_rare") %>%
  ggplot(., aes(x = pred, y = pred_rare)) +
  geom_point() +
  facet_wrap(~lifestyle + rarefaction_type)

lifestyle_predictions_genus %>% left_join(lifestyle_predictions_genus_rare_h %>% select(run_accession, pred),
                                          suffix = c("", "_rare_h"), by = "run_accession") %>%
  left_join(lifestyle_predictions_genus_rare_l %>% select(run_accession, pred), suffix = c("", "_rare_l"), by = "run_accession") %>%
  select(pred, pred_rare_h, pred_rare_l) %>%
  cor(method = "spearman", use = "pairwise.complete.obs") %>%
  pheatmap::pheatmap(display_numbers = T)

# compare with prediction from same lifestyle:
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_industrialized_rare_h.RData")
genus_industrialized_nested_cv_rare_h_preds_long <- genus_industrialized_nested_cv_rare_h_preds %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "industrialized")

industrialized_lm_list_genus_rare_h <- get_lm_list(pred_df = genus_industrialized_nested_cv_rare_h_preds_long,
                                                   grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "industrialized")

load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_nonindustrialized_rare_h.RData")
genus_nonindustrialized_nested_cv_rare_h_preds_long <- genus_nonindustrialized_nested_cv_rare_h_preds %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "non_industrialized")
nonindustrialized_lm_list_genus_rare_h <- get_lm_list(pred_df = genus_nonindustrialized_nested_cv_rare_h_preds_long, grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "non_industrialized")


# default model:
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_no_ls_rare_h.RData")
genus_nested_cv_rare_h_preds_long <- genus_no_ls_nested_rare_h_cv_preds %>%
  pivot_longer(cols = c("rf1"), names_to = "model_name", values_to = "pred" ) %>%
  mutate(training_set = "combined")
all_lm_list_genus_rare_h <- get_lm_list(pred_df = genus_nested_cv_rare_h_preds_long, grouping = c("study", "lifestyle", "model_name")) %>%
  filter(model_name == "rf1") %>%
  mutate(training_set = "combined")


lifestyle_predictions_genus_rare_h <- lifestyle_predictions_genus_rare_h %>%
  mutate(model_name = "rf1") %>%
  mutate(training_set = ifelse(lifestyle == "industrialized", yes = "non_industrialized", no = "industrialized"))
lifestyle_lm_list_genus_rare_h <- get_lm_list(pred_df = lifestyle_predictions_genus_rare_h, grouping = c("study", "lifestyle", "model_name")) %>%
  mutate(training_set = ifelse(lifestyle == "industrialized", yes = "non_industrialized", no = "industrialized")) %>%
  rbind(., industrialized_lm_list_genus_rare_h, nonindustrialized_lm_list_genus_rare_h, all_lm_list_genus_rare_h)

lifestyle_lm_list_genus_rare_h <- lifestyle_lm_list_genus_rare_h %>%
  mutate(training_set = factor(training_set, levels = c("non_industrialized", "combined", "industrialized")),
         lifestyle = factor(lifestyle, levels = c("industrialized", "non_industrialized")))

# which studies differ between rarefied and non-rarefied?
test <- left_join(lifestyle_lm_list_genus_rare_h, lifestyle_lm_list_genus, by = c("study", "lifestyle", "model_name",
                                                                                  "training_set"), suffix = c("_rarefied", ""))
test %>% group_by(lifestyle, training_set) %>%
  summarize(R2_mean = mean(R2),
            R2_raref_mean = mean(R2_rarefied))
ggplot(test, aes(x = R2, y = R2_rarefied, color = lifestyle)) +
  geom_point() +
  geom_abline(slope = 1) +
  facet_wrap(~training_set)

# stat test (stay with R2, RSME performs worse)
df_p_val_lifestyle_genus_rare_h <- lifestyle_lm_list_genus_rare_h %>%
  # lifestyle_lm_list_genus %>%
  # filter(study != "bender_2016") %>%
  arrange(study) %>%
  mutate(training_set = factor(training_set)) %>%
  rstatix::group_by(lifestyle) %>%
  rstatix::wilcox_test(R2 ~ training_set, paired = T) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>% 
  rstatix::add_xy_position(x = "lifestyle", dodge = 0.8) 

ls_performance_rare_h <- ggplot(lifestyle_lm_list_genus_rare_h, aes(x=lifestyle, y = R2)) +
  geom_boxplot(aes(fill = training_set)) +
  # ylim(0, 1) +
  xlab("Test set") +
  add_pvalue(df_p_val_lifestyle_genus_rare_h,
             label = "{p.adj.signif}",
             # step.group.by = "variation",
             step.increase = 0.05,
             tip.length = 0.01,
             # bracket.nudge.y = 0.02,
             xmin = "xmin", 
             xmax = "xmax",
             show.legend = FALSE) +
  ylim(0,NA) +
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
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/lifestyle_perfomance_genus_raref_h.pdf", 
       width = 10, plot = ls_performance_rare_h)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/lifestyle_perfomance_genus_raref_h.png", 
       width = 10, plot = ls_performance_rare_h)
save(lifestyle_lm_list_genus_rare_h, ls_performance_rare_h,
     file = "/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/lifestyle_perfomance_genus_raref_h.RData")

bind_rows(lifestyle_lm_list_genus_rare_h %>% mutate(run = "rarefied"),
          lifestyle_lm_list_genus %>% mutate(run = "raw")) %>%
  ggplot(., aes(x=lifestyle, y = R2, color = run)) +
  geom_boxplot(aes(fill = training_set), alpha = 0.5)

bind_rows(lifestyle_lm_list_genus_rare_h %>% mutate(run = "rarefied"),
          lifestyle_lm_list_genus %>% mutate(run = "raw")) %>%
  arrange(study) %>%
  mutate(training_set = factor(training_set)) %>%
  rstatix::group_by(lifestyle, training_set) %>%
  rstatix::wilcox_test(R2 ~ run, paired = T) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1))
  
  
## different feature importance #################################################
combined_importance_nonindustrialized_genus_rare_h <- map(non_industrialized_studies, 
                                                   function(x) get_feature_importance(x, lst = "nonindustrialized_rare_h")) %>%
  purrr::reduce(rbind) %>%
  mutate(lifestyle = "non_industrialized") %>%
  group_by(taxon) %>%
  filter(n() > 1) %>% # only take those that occur in more than 2 models
  ungroup()
combined_importance_industrialized_genus_rare_h <- map(industrialized_studies, function(x) get_feature_importance(x, lst = "industrialized_rare_h")) %>%
  purrr::reduce(rbind) %>%
  mutate(lifestyle = "industrialized") %>%
  group_by(taxon) %>%
  filter(n() > 1) %>% # only take those that occur in more than 2 models
  ungroup()

bind_rows(combined_importance_nonindustrialized_genus_rare_h, combined_importance_industrialized_genus_rare_h) %>%
  group_by(ID, lifestyle) %>%
  dplyr::summarize(count = n()) %>%
  ggplot(., aes( x = lifestyle, y = count)) +
  geom_boxplot() +
  geom_jitter() +
  ylim(0, 100)

## important taxa stats #################
# n important samples 
important_union_rare_h <- unique(c(combined_importance_nonindustrialized_genus_rare_h$taxon, 
                            combined_importance_industrialized_genus_rare_h$taxon))
important_industrialized_rare_h <- unique(combined_importance_industrialized_genus_rare_h$taxon) 
length(important_industrialized_rare_h)
important_non_industrialized_rare_h <- unique(combined_importance_nonindustrialized_genus_rare_h$taxon) 
length(important_non_industrialized_rare_h)
# shared:
length(intersect(important_industrialized_rare_h, important_non_industrialized_rare_h))
# specific to one:
important_industrialized_rare_h[!important_industrialized_rare_h %in% important_non_industrialized_rare_h] %>% length
important_non_industrialized_rare_h[!important_non_industrialized_rare_h %in% important_industrialized_rare_h] %>% length





# relative abundance of important taxa in data:
for_caret_list_rare_h <- create_caret_df(ps_object = ps_object_genus_raw %>% rarefy_even_depth(., sample.size = 2000),
                                  transformation = "compositional",
                                  mean_ab_cutoff = 5e-5, study_prevalence_cutoff = 1,
                                  additional_cols = c("lifestyle", "study"),
                                  only_multiple_samples = F)



# add prevalence:
prevalence_lifestyle_rare_h <- for_caret_list_rare_h$features %>%
  # mutate(study_lifestyle = paste(lifestyle, study.1, sep = "xxx")) %>%
  group_by(lifestyle, study.1) %>% 
  summarize(across(everything(), ~ mean(.>0)), .groups = "drop") %>%  # get prevalence, mean over studies
  select(-study.1) %>%
  group_by(lifestyle) %>%
  summarize(across(everything(), ~ mean(.)), .groups = "drop") %>%
  t() %>%
  data.frame %>%
  `colnames<-`(.[1, ]) %>%
  .[-1, ] %>% 
  rownames_to_column(.,var = "taxon") %>% 
  pivot_longer(., -taxon, names_to = "lifestyle", values_to = "prevalence") %>%
  mutate(prevalence = as.numeric(prevalence) * 100)# %>%

combined_importance_prev_lifestyles_genus_rare_h <- rbind(combined_importance_industrialized_genus_rare_h, 
                                                   combined_importance_nonindustrialized_genus_rare_h) %>%
  complete(taxon, lifestyle) %>%
  mutate(lifestyle = factor(lifestyle, levels = c("industrialized", "non_industrialized")),
         model = lifestyle,
         measure_type = ifelse(taxon %in% c("Shannon", "Observed"), yes = "a_div", no = "taxon")) %>%
  left_join(., prevalence_lifestyle_rare_h, by = c("taxon", "lifestyle"))

# saveRDS(combined_importance_prev_lifestyles_genus,
#         "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/genus_nested_cv_important_taxa_lifestyle.rds")

# shorten taxon names
short_taxa_names_rare_h <- combined_importance_prev_lifestyles_genus_rare_h %>% select(taxon) %>%
  distinct() %>%
  mutate(display_taxon = taxon %>% gsub(".*unclassified_", "Uncl_ ",.) %>%
           make.unique() %>%
           gsub("X[.]", "", .) %>%
           gsub("[.]", " ", .) %>%
           gsub("  ", " ", .) %>%
           gsub("_", ".", .) %>%
           gsub("Family XIII AD3011 group", "Anaerovoracaceae Fam. XIII AD3011", .) %>%
           gsub("UCG 002", "Oscillospiraceae UCG 002", .))

combined_importance_prev_lifestyles_genus_rare_h <- combined_importance_prev_lifestyles_genus_rare_h %>% 
  left_join(.,short_taxa_names_rare_h, by = "taxon") %>%
  # reorder according to importance
  mutate(display_taxon = reorder(display_taxon, importance, FUN = max),
         display_taxon = factor(display_taxon, levels = unique(c("Shannon", "Observed", levels(display_taxon)))),
         lifestyle = factor(lifestyle, levels = c("non_industrialized", "industrialized"))
  )

(importance_prevalence_plot_rare_h <- bind_rows(combined_importance_prev_lifestyles_genus_rare_h %>% mutate(run = "rarefied"),
                                                combined_importance_prev_lifestyles_genus %>% mutate(run = "raw")) %>%
    ggplot(.,aes(y = display_taxon, x = importance, color = run)) +
    geom_col(data = bind_rows(combined_importance_prev_lifestyles_genus_rare_h %>%
                                mutate(run = "rarefied"),
                              combined_importance_prev_lifestyles_genus %>%
                                mutate(run = "raw")) %>%
               distinct(display_taxon, prevalence, model, run),
             aes(y = display_taxon, x = prevalence, fill = run), alpha = 0.5,,
             position = position_dodge(width = 0.8)) +
    geom_boxplot() +
    xlim(0,100) +
    xlab("Relative importance / Prevalence") +
    facet_wrap(~ model, nrow = 1, drop = T) +
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(), #element_text(size = 14),
          axis.text.y = element_text(size = 14),
          # legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_line(color = "grey",
                                            size = 0.3,
                                            linetype = 2),
          strip.text = element_text(size = 14),
          plot.margin = margin(0.2,0.2,0.2,1, "cm")))
# the non-rarefied version has much more important taxa
unique(c(combined_importance_prev_lifestyles_genus_rare_h$display_taxon, 
         combined_importance_prev_lifestyles_genus$display_taxon))
left_join(combined_importance_prev_lifestyles_genus_rare_h %>%
            pivot_longer(cols = c("importance", "prevalence"), names_to = "measure", values_to = "values"), 
          combined_importance_prev_lifestyles_genus %>%
            pivot_longer(cols = c("importance", "prevalence"), names_to = "measure", values_to = "values"),
          by = c("taxon", "lifestyle", "model", "display_taxon", "measure_type", "ID", "measure"),
          suffix = c("_rarefied", "_raw")) %>%
  ggplot(., aes(x = values_raw, y = values_rarefied, color = model)) +
  geom_point() +
  geom_abline(slope = 1) +
  facet_wrap(~measure, scale = "free")



# ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/genus_nested_cv_importance_lifestyle.pdf",
#        height = 17, plot = importance_prevalence_plot)
# ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/genus_nested_cv_importance_lifestyle.png",
#        height = 17, plot = importance_prevalence_plot)

# save(combined_importance_prev_lifestyles_genus,
# file = "/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/imp_prevs.RData")


## shap analysis genus ################################################################
# combined model
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_no_ls_genus.RData")
# final_model_genus_no_ls
# shap_out_genus_no_ls
# industrialized model
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_industrialized_genus.RData")
# final_model_genus_industrialized
# shap_out_genus_industrialized
# non industrialized model
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_nonindustrialized_genus.RData")
# final_model_genus_nonindustrialized
# shap_out_genus_nonindustrialized
# genus_train_data_non_industrialized

# important taxa:
print("run shap analysis")
if(shap_step){
  shap_model_all <- get_shap_long(final_model_genus_no_ls$rf1,
                                  test_data = genus_train_data$features,
                                  features = unique(combined_importance_prev_lifestyles_genus$taxon))
  print("done shap all")
  shap_model_industrialized <- get_shap_long(final_model_genus_industrialized$rf1, 
                                             test_data = genus_train_data$features,
                                             features = unique(combined_importance_prev_lifestyles_genus$taxon))
  print("done shap industrialized")
  shap_model_noindustrialized <- get_shap_long(final_model_genus_nonindustrialized$rf1,
                                               test_data = genus_train_data$features,
                                               features = unique(combined_importance_prev_lifestyles_genus$taxon))
  print("done shap non industrialized")
  save(shap_model_all, shap_model_industrialized, shap_model_noindustrialized, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/lifestyle_shap_data.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/lifestyle_shap_data.RData")
}

## shap vs abundances #################

shap_model_all %>%
  left_join(short_taxa_names, by = "taxon") %>%
  left_join(., metadata, by = "run_accession") %>%
  group_by(display_taxon, lifestyle) %>% 
  summarise(R2 = cor(ab_value, shap_value)) %>%
  # filter(!display_taxon %in% c("Observed", "Shannon")) %>%
  ggplot(., aes(y = display_taxon, x = R2, fill = lifestyle, group = lifestyle)) +
  geom_bar(position="dodge", stat="identity") +
  theme_classic() +
  theme(panel.grid.major.y = element_line(color = "grey",
                                          size = 0.3,
                                          linetype = 2))+
  ggtitle("Combined_model")

merged_shap_models <- rbind(shap_model_industrialized %>% mutate(model = "industrialized"), 
                            shap_model_noindustrialized %>% mutate(model = "non_industrialized")) %>%
  left_join(., short_taxa_names, by = "taxon") %>%
  left_join(., metadata, by = "run_accession")

shap_cors <- merged_shap_models %>%
  group_by(display_taxon, model) %>%  
  # lifestyle  maybe only check for all data combined, because it's more a property of the model than of the data?
  summarise(cor_test = list(cor.test(ab_value, shap_value, method = "spearman", use = "pairwise.complete.obs")),
            display_taxon = unique(display_taxon),
            .groups = "drop") %>%
  mutate(R2 = map_dbl(cor_test, ~ .x$estimate),
         p_value = map_dbl(cor_test, ~ .x$p.value),
         p_adj = p.adjust(p_value, method = "bonferroni"),
         p_adj = ifelse(p_adj == 0,
                        yes = min(p_adj[p_adj!=0], na.rm = T), # remove 0 p-values
                        no = p_adj),
         # keep the order of the first figure
         display_taxon = factor(display_taxon, 
                                levels = levels(combined_importance_prev_lifestyles_genus$display_taxon)))

sha_value_cor <- shap_cors %>%
  ggplot(., aes(y = display_taxon, x = R2, fill = model)) + # fill = lifestyle, group = lifestyle
  geom_bar(position="dodge", stat="identity") +
  facet_wrap(~ model, nrow = 1, drop = T) +
  xlab("Spearman(feature abund, SHAP)") +
  ylab("") +
  theme(axis.text.x = element_text(size = 14),
        # axis.title.x = element_blank(),
        axis.title.y = element_blank(), #element_text(size = 14),
        axis.text.y = element_blank(), #element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        panel.grid.major.y = element_line(color = "grey",
                                          size = 0.3,
                                          linetype = 2),
        strip.text = element_text(size = 14),
        plot.margin = margin(0.2,0.2,0.2,1, "cm"))

ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_lifestyle_models_cor.pdf",
       plot = sha_value_cor, height = 15, width = 6)#, width = 9, height = 10)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_lifestyle_models_cor.png",
       plot = sha_value_cor, height = 15, width = 6)#, width = 9, height = 10)

# combine with importance plot:
library(ggpubr)
ggarrange(importance_prevalence_plot,sha_value_cor)


# correlation between shap and abund for each taxon in different models based on different data
sha_value_cor_point <- shap_cors %>%
  mutate(R2 = case_when(is.na(R2) ~ 0,
                        p_adj > 0.05 ~ 0,
                        .default = R2)) %>%
  select(-cor_test, -p_value, -p_adj) %>%
  pivot_wider(names_from = "model", values_from = R2) %>%
  filter(industrialized != 0 | non_industrialized != 0) %>%
  # filter(!taxon %in% c("X.Eubacterium..coprostanoligenes.group", "UCG.002", "Observed",
  #                      "Shannon", "Lachnospiraceae.ND3007.group", "Intestinibacter",
  #                      "Erysipelatoclostridium", "Family.XIII.AD3011.group",
  #                      "Anaerostipes", "Oscillibacter", "Catenibacterium", "uncl_Firmicutes",
  #                      "X.Eubacterium..eligens.group", "Erysipelotrichaceae.UCG.003")) %>%
  ggplot(., aes(x = industrialized, y = non_industrialized)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black") +
  # Set axis limits and expand to ensure 0 is included
  # scale_x_continuous(expand = c(0, 0)) +
  # scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size = 14),
        plot.margin = margin(0.2,0.2,0.2,1, "cm"))
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_lifestyle_models_cor_point.pdf",
       plot = sha_value_cor_point)#, width = 9, height = 10)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/figs/shap_lifestyle_models_cor_point.png",
       plot = sha_value_cor_point)#, width = 9, height = 10)

