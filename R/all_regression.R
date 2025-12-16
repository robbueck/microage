setwd("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/")
# another idea: instead of abundances, use ratios between taxa
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
source("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/alt_models.R")
source("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/regression_functions.R")
source("/fast/AG_Forslund/rob/mm_index/R_scripts/setlists.R")

dataset_cv_family_step <- F
dataset_cv_family_one_year_step <- F
dataset_cv_family_age_month_step <- F
dataset_cv_family_bin_step <- F
dataset_cv_family_noabx_step <- F
dataset_cv_family_western_step <- F
dataset_cv_family_nonwestern_step <- F
dataset_cv_family_clr_step <- T
dataset_cv_family_ratio_step <- T
dataset_cv_genus_step <- T


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

cl <- makePSOCKcluster(n_cores)
registerDoParallel(cl)
getDoParWorkers()
set.seed(825)

#############
# family data
#############
print("family data")
ps_object_family_raw <- readRDS("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_phyloseq_rf_filter_family.rds") %>%
  subset_samples(., age <= 730 & age > 1)
  
ps_object_family_comp <- ps_object_family_raw %>%
  # filter_taxa(function(x){sum(x > 0) > 5}, TRUE) %>%  # prevalence cutoff 10 % of the smallest study
  microbiome::transform(transform = "compositional") %>%
  filter_taxa(function(x) mean(x[x > 0]) > 5e-5, TRUE) %>%  # mean abundance cutoff, only counting abundances > 0
  prune_samples(samples = (sample_sums(.) != 0))

# remove single sample subjects
multiple_samples <- ps_object_family_comp@sam_data %>%
  data.frame() %>%
  filter(duplicated(subject_ID)) %>%
  select(subject_ID) %>%
  unique %>%
  deframe
ps_object_family_comp <- ps_object_family_comp %>% subset_samples(subject_ID %in% multiple_samples)
ps_object_family_raw <- ps_object_family_raw %>% 
  prune_taxa(taxa_names(ps_object_family_comp), .) %>%
  prune_samples(sample_names(ps_object_family_comp),.) %>%
  subset_samples(subject_ID %in% multiple_samples)


# add richness and shannon diversity to the metadata
diversities <- estimate_richness(ps_object_family_raw, measures = c("Observed", "Shannon"))
sample_data(ps_object_family_comp)$observed <- diversities$Observed
sample_data(ps_object_family_comp)$shannon <- diversities$Shannon
metadata_df <- sample_data(ps_object_family_comp) %>% data.frame()
# xtabs(~ study + lifestyle, data = metadata_df)
metadata_df$lifestyle_westernized <- ifelse(metadata_df$lifestyle == "westernized", yes = 1, no = 0)

data_set_family <- ps_object_family_comp%>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  # select(!Unknown) %>%
  cbind(metadata_df[,c("observed", "shannon", "lifestyle_westernized", "study")])

# prevalence filtering, only keep taxa present in all studies
cols_to_keep <- data_set_family %>% group_by(study) %>%
  summarize(across(everything(), ~ sum(.))) %>% 
  select(-study) %>%
  select(where(~ sum(. == 0) <= 1)) %>% # adjust threshold, in how many studies can a taxon be missing
  colnames() %>%
  c(., "lifestyle_westernized") # is either 1 or 0 for all samples of one study
data_set_family <- data_set_family[,cols_to_keep]


# remove near zero variance predictors
# nzv_family=preProcess(data_set_family,method="nzv",uniqueCut = 5)
# data_set_family_cleaned <- predict(nzv_family,data_set_family) %>%
#   mutate(Eubacteriaceae = X.Eubacterium..coprostanoligenes.group, .keep = "unused")
# family_cleaned_cor <- cor(data_set_family_cleaned)
# highlyCorrelated_family_cleaned <- findCorrelation(family_cleaned_cor, cutoff=0.5)
if(any(rowSums(data_set_family) == 0)) {
  warning("One or more samples do not contain any abundance for any of the selected taxa")
}

data_set_family <- data_set_family %>%
  mutate(Eubacteriaceae = X.Eubacterium..coprostanoligenes.group, .keep = "unused")


# # create correlation plot
# cor_network <- cor(data_set_family)
# cor_network
# pheatmap::pheatmap(cor_network, treeheight_row = 0, treeheight_col = 0)
# ggsave("showcase_correlation_plot.pdf")

print("remove studies from family data")
tic()
if(dataset_cv_family_step){
  splits <- group_vfold_cv(metadata_df, group = "study", v = length(unique(metadata_df$study)))
  caret_split <- rsample2caret(splits)
  tc_grouped <- trainControl(method = "cv",
                             savePredictions = "final",
                             selectionFunction = "oneSE", # function to select the best model
                             index = caret_split$index,
                             indexOut = caret_split$indexOut,
                             preProcOptions = list(uniqueCut = 5),
                             allowParallel = T)
  res_family_study <- caretList(y=metadata_df$age, x=data_set_family,
                                metric="Rsquared",
                                trControl=tc_grouped,
                                preProcess = c("nzv"),
                                methodList=c("lasso", "glmnet"),
                                tuneList = list(
                                  enet=caretModelSpec(method="enet", verbose = F),
                                  gbm=caretModelSpec(method="gbm", verbose = F, tuneGrid=expand.grid(shrinkage=c(0.1, 0.2, 0.3),
                                                                                                     n.minobsinnode=c(5, 15),
                                                                                                     n.trees=c(100, 200),
                                                                                                     interaction.depth=c(2, 3, 4))),
                                  kknn=caretModelSpec(method = 'kknn', tuneGrid=expand.grid(kmax = c(5, 9, 13, 15, 17),
                                                                                            distance = c(1, 5, 9),
                                                                                            kernel = "optimal")),
                                  # xgb=caretModelSpec(method="xgbTree", tuneGrid=expand.grid(nrounds = seq(from = 200, to = 1000, by = 50),
                                  #                                                           eta = c(0.025, 0.05, 0.1, 0.3),
                                  #                                                           max_depth = c(2, 3, 4, 5, 6),
                                  #                                                           gamma = 0,
                                  #                                                           colsample_bytree = 1,
                                  #                                                           min_child_weight = 1,
                                  #                                                           subsample = 1)),
                                  rf1=caretModelSpec(method=adapt_ranger, tuneGrid=expand.grid(mtry = c(3, 10, 20),
                                                                                               num.trees = c(20, 50, 200, 400),
                                                                                               min.node.size = c(5, 10, 20),
                                                                                               splitrule = c("variance"),
                                                                                               replace = F),
                                                     importance = "permutation")
                                )
  )
  plots_family_study <- modelplots(res_family_study, metadata_df = metadata_df, all_studies = set_list$all, color = "study")

  save(res_family_study, plots_family_study,
       file = paste0("data/", merged_set, "_regression_models_dataset_family_cv.RData"))
} else {
  load(paste0("data/", merged_set, "_regression_models_dataset_family_cv.RData"))
}
toc()

plots_family_study$importance$rf1 +
  ggtitle("") +
  theme(
    axis.title.x = element_text(size = 24),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 24))

ggsave("figs/importance_family.pdf", width = 12, height = 6)

# track dynamics over time
importances <- varImp(res_family_study$rf1)
top_fams <- importances$importance %>% arrange(-Overall) %>% rownames(.) %>% head(10)
top_fams <- top_fams[top_fams!= "lifestyle_westernized"]
# top_fams <- c("Lachnospiraceae", "Ruminococcaceae", "Bifidobacteriaceae", "Micrococcaceae")
# top_fams <- c("Bifidobacteriaceae")
# top_fams <- "Micrococcaceae"
top_fams <- c("Lachnospiraceae", "Ruminococcaceae", "Veillonellaceae", "Micrococcaceae")

df_to_plot <- merge(metadata_df[,c("age", "study")], data_set_family[,c(top_fams, "lifestyle_westernized")], by = 0) %>%
  pivot_longer(all_of(top_fams)) %>%
  filter(!name %in% c("observed", "shannon")) 

sum(data_set_family$Lachnospiraceae != 0) / nrow(data_set_family)
sum(data_set_family$Ruminococcaceae != 0) / nrow(data_set_family)
sum(data_set_family$Veillonellaceae != 0) / nrow(data_set_family)
sum(data_set_family$Micrococcaceae != 0) / nrow(data_set_family)
df_to_plot$name <- df_to_plot$name %>%
  gsub("Lachnospiraceae", "Lachnospiraceae (0.61)", .) %>%
  gsub("Ruminococcaceae", "Ruminococcaceae (0.26)", .) %>%
  gsub("Veillonellaceae", "Veillonellaceae (0.75)", .) %>%
  gsub("Micrococcaceae", "Micrococcaceae (0.30)" ,.)
  

ggplot(df_to_plot, aes(x = age, y = value)) +
  geom_point(alpha = 0.7) + 
  facet_wrap(~name, ncol = 2) +
  ylab("relative abundance") +
  xlab("age [days]") + 
  scale_y_continuous(trans='log10') +
  theme(
    axis.title.x = element_text(size = 24),
    axis.text.x = element_text(size = 22),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 22),
    strip.text = element_text(size=22))
ggsave("figs/fams_over_time.pdf", width = 10, height = 4)

rf_plots <- modelplots(res_family_study$rf1, metadata_df = metadata_df, all_studies = set_list$all, color = "study")
rf_plots$correlations +
  theme( strip.background = element_blank(),
         strip.text.x = element_blank()
  )
ggsave("figs/age_vs_pred_all.pdf", width = 5, height = 5)
rf_plots$linear_model +
  xlim(c(0,1)) +
  ylim(c(0,1))+
  theme( strip.background = element_blank(),
         strip.text.x = element_blank(),
         legend.position=c(0.1, 0.8)
  )
ggsave("figs/model_plots_all.pdf", width = 4, height = 4)

xyplot(resamples(res_family_study[c(1,3)]))
modelCor(resamples(res_family_study))

ans = resamples(res_family_study) #resamples helps to tabularize the results
# summary(ans)
# dotplot(ans)
# stack.glm = caretStack(res_family_study[c(1, 3, 4)], method="rf", trControl=tc_grouped) #logistic
# stack.glm = caretStack(res_family_study[c(4)], method="lm", trControl=tc_grouped) #logistic
# print(stack.glm)
# ensemble_plots <- modelplots(stack.glm$ens_model, metadata_df = metadata_df, all_studies = set_list$all, color = "study")
# rf_plots <- modelplots(res_family_study$rf1, metadata_df = metadata_df, all_studies = set_list$all, color = "study")
# grid.arrange(rf_plots$linear_model, ensemble_plots$linear_model)
# 
# 
p1 <- ggplot(res_family_study$rf1$results,# %>% filter(replace),
       aes(x=num.trees,y=Rsquared,col=factor(mtry), shape = splitrule)) +
  geom_line() + geom_point() + facet_grid(~min.node.size) +
  ggtitle("replace = T")
p2 <- ggplot(res_family_study$rf1$results,# %>% filter(replace),
             aes(x=num.trees,y=Rsquared,col=factor(mtry), shape = splitrule)) +
  geom_line() + geom_point() + facet_grid(~min.node.size) +
  ggtitle("replace = T")
p2 <- ggplot(res_family_study$rf1$results %>% filter(!replace),
       aes(x=num.trees,y=Rsquared,col=factor(mtry), shape = splitrule)) +
  geom_line() + geom_point() + facet_grid(~min.node.size) +
  ggtitle("replace = F")
gridExtra::grid.arrange(p1, p2)
ggplot(res_family_study$gbm$results,
       aes(x=shrinkage,y=Rsquared,col=factor(n.trees), shape = factor(n.minobsinnode))) +
  geom_line() + geom_point() + facet_grid(~interaction.depth)
ggplot(res_family_study$enet)
ggplot(res_family_study$kknn)
ggplot(res_family_study$lasso)
ggplot(res_family_study$glmnet)



################################################################################
# family data with max one year
################################################################################
print("family data one year")
ps_object_family_one_raw <- readRDS("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_phyloseq_rf_filter_family.rds") %>%
  subset_samples(., age <= 365 & age > 1)

ps_object_family_one_comp <- ps_object_family_one_raw %>%
  # filter_taxa(function(x){sum(x > 0) > 5}, TRUE) %>%  # prevalence cutoff 10 % of the smallest study
  microbiome::transform(transform = "compositional") %>%
  filter_taxa(function(x) mean(x[x > 0]) > 5e-5, TRUE) %>%  # mean abundance cutoff, only counting abundances > 0
  prune_samples(samples = (sample_sums(.) != 0))

# remove single sample subjects
multiple_samples_one <- ps_object_family_one_comp@sam_data %>%
  data.frame() %>%
  filter(duplicated(subject_ID)) %>%
  select(subject_ID) %>%
  unique %>%
  deframe
ps_object_family_one_comp <- ps_object_family_one_comp %>% subset_samples(subject_ID %in% multiple_samples_one)
ps_object_family_one_raw <- ps_object_family_one_raw %>% 
  prune_taxa(taxa_names(ps_object_family_one_comp), .) %>%
  prune_samples(sample_names(ps_object_family_one_comp),.) %>%
  subset_samples(subject_ID %in% multiple_samples_one)


# add richness and shannon diversity to the metadata
diversities <- estimate_richness(ps_object_family_one_raw, measures = c("Observed", "Shannon"))
sample_data(ps_object_family_one_comp)$observed <- diversities$Observed
sample_data(ps_object_family_one_comp)$shannon <- diversities$Shannon
metadata_df_one <- sample_data(ps_object_family_one_comp) %>% data.frame()
# xtabs(~ study + lifestyle, data = metadata_df_one)
metadata_df_one$lifestyle_westernized <- ifelse(metadata_df_one$lifestyle == "westernized", yes = 1, no = 0)

data_set_family_one <- ps_object_family_one_comp%>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  # select(!Unknown) %>%
  cbind(metadata_df_one[,c("observed", "shannon", "lifestyle_westernized", "study")])

# prevalence filtering, only keep taxa present in all studies
cols_to_keep <- data_set_family_one %>% group_by(study) %>%
  summarize(across(everything(), ~ sum(.))) %>% 
  select(-study) %>%
  select(where(~ sum(. == 0) <= 1)) %>% # adjust threshold, in how many studies can a taxon be missing
  colnames() %>%
  c(., "lifestyle_westernized") # is either 1 or 0 for all samples of one study
data_set_family_one <- data_set_family_one[,cols_to_keep]

if(any(rowSums(data_set_family_one) == 0)) {
  warning("One or more samples do not contain any abundance for any of the selected taxa")
}


print("remove studies from family data")
tic()
if(dataset_cv_family_one_year_step){
  splits <- group_vfold_cv(metadata_df_one, group = "study", v = length(unique(metadata_df_one$study)))
  caret_split <- rsample2caret(splits)
  tc_grouped <- trainControl(method = "cv",
                             savePredictions = "final",
                             selectionFunction = "oneSE", # function to select the best model
                             index = caret_split$index,
                             indexOut = caret_split$indexOut,
                             preProcOptions = list(uniqueCut = 5),
                             allowParallel = T)
  res_family_one_study <- caretList(y=metadata_df_one$age, x=data_set_family_one,
                                metric="Rsquared",
                                trControl=tc_grouped,
                                preProcess = c("nzv"),
                                methodList=c("lasso", "glmnet"),
                                tuneList = list(
                                  enet=caretModelSpec(method="enet", verbose = F),
                                  gbm=caretModelSpec(method="gbm", verbose = F, tuneGrid=expand.grid(shrinkage=c(0.1, 0.2, 0.3),
                                                                                                     n.minobsinnode=c(5, 15),
                                                                                                     n.trees=c(100, 200),
                                                                                                     interaction.depth=c(2, 3, 4))),
                                  kknn=caretModelSpec(method = 'kknn', tuneGrid=expand.grid(kmax = c(5, 9, 13, 15, 17),
                                                                                            distance = c(1, 5, 9),
                                                                                            kernel = "optimal")),
                                  # xgb=caretModelSpec(method="xgbTree", tuneGrid=expand.grid(nrounds = seq(from = 200, to = 1000, by = 50),
                                  #                                                           eta = c(0.025, 0.05, 0.1, 0.3),
                                  #                                                           max_depth = c(2, 3, 4, 5, 6),
                                  #                                                           gamma = 0,
                                  #                                                           colsample_bytree = 1,
                                  #                                                           min_child_weight = 1,
                                  #                                                           subsample = 1)),
                                  rf1=caretModelSpec(method=adapt_ranger, tuneGrid=expand.grid(mtry = c(3, 10, 20),
                                                                                               num.trees = c(20, 50, 200, 400),
                                                                                               min.node.size = c(5, 10, 20),
                                                                                               splitrule = c("variance"),
                                                                                               replace = F),
                                                     importance = "permutation")
                                )
  )
  plots_family_one_study <- modelplots(res_family_one_study, metadata_df = metadata_df_one, all_studies = set_list$all, color = "study")
  
  save(res_family_one_study, plots_family_one_study,
       file = paste0("data/", merged_set, "_regression_models_dataset_family_one_year_cv.RData"))
} else {
  load(paste0("data/", merged_set, "_regression_models_dataset_family_one_year_cv.RData"))
}
toc()




################################################################################
# family data with age rounded to months
################################################################################
metadata_df_month <- metadata_df
metadata_df_month$age <- ceiling(metadata_df$age/30)

age_study <- xtabs(~ study + age, data = metadata_df_month) %>% as.data.frame.matrix
# View(age_study)
colSums(age_study > 0)
# xtabs(~ study + lifestyle, data = metadata_df_month)


if(any(rowSums(data_set_family) == 0)) {
  warning("One or more samples do not contain any abundance for any of the selected taxa")
}


print("remove studies from family data")
tic()
if(dataset_cv_family_age_month_step){
  splits <- group_vfold_cv(metadata_df_month, group = "study", v = length(unique(metadata_df_month$study)))
  caret_split <- rsample2caret(splits)
  tc_grouped <- trainControl(method = "cv",
                             savePredictions = "final",
                             selectionFunction = "oneSE", # function to select the best model
                             index = caret_split$index,
                             indexOut = caret_split$indexOut,
                             preProcOptions = list(uniqueCut = 5),
                             allowParallel = T)
  res_family_study_month <- caretList(y=metadata_df_month$age, x=data_set_family,
                                metric="Rsquared",
                                trControl=tc_grouped,
                                preProcess = c("nzv"),
                                methodList=c("lasso", "glmnet"),
                                tuneList = list(
                                  enet=caretModelSpec(method="enet", verbose = F),
                                  gbm=caretModelSpec(method="gbm", verbose = F, tuneGrid=expand.grid(shrinkage=c(0.1, 0.2, 0.3),
                                                                                                     n.minobsinnode=c(5, 15),
                                                                                                     n.trees=c(100, 200),
                                                                                                     interaction.depth=c(2, 3, 4))),
                                  kknn=caretModelSpec(method = 'kknn', tuneGrid=expand.grid(kmax = c(5, 9, 13, 15, 17),
                                                                                            distance = c(1, 5, 9),
                                                                                            kernel = "optimal")),
                                  # xgb=caretModelSpec(method="xgbTree", tuneGrid=expand.grid(nrounds = seq(from = 200, to = 1000, by = 50),
                                  #                                                           eta = c(0.025, 0.05, 0.1, 0.3),
                                  #                                                           max_depth = c(2, 3, 4, 5, 6),
                                  #                                                           gamma = 0,
                                  #                                                           colsample_bytree = 1,
                                  #                                                           min_child_weight = 1,
                                  #                                                           subsample = 1)),
                                  rf1=caretModelSpec(method=adapt_ranger, tuneGrid=expand.grid(mtry = c(3, 10, 20),
                                                                                               num.trees = c(20, 50, 200, 400),
                                                                                               min.node.size = c(5, 10, 20),
                                                                                               splitrule = c("variance"),
                                                                                               replace = F),
                                                     importance = "permutation")
                                )
  )
  plots_family_study_month <- modelplots(res_family_study_month, metadata_df = metadata_df_month, all_studies = set_list$all, color = "study")
  
  save(res_family_study_month, plots_family_study_month,
       file = paste0("data/", merged_set, "_regression_models_dataset_family_month_cv.RData"))
} else {
  load(paste0("data/", merged_set, "_regression_models_dataset_family_month_cv.RData"))
}
toc()



####################
# binary family data
####################
# compare to other models, if abundance gives more accuracy compared to presence-absence

data_set_family_bin <- ps_object_family_comp %>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, 0))) %>%
  select(!Unknown) %>%
  cbind(metadata_df[,c("observed", "shannon", "lifestyle_westernized", "study")])


# prevalence filtering, only keep taxa present in all studies
cols_to_keep <- data_set_family_bin %>% group_by(study) %>%
  summarize(across(everything(), ~ sum(.))) %>% 
  select(-study) %>%
  select(where(~ sum(. == 0) <= 1)) %>% # adjust threshold, in how many studies can a taxon be missing
  colnames() %>%
  c(., "lifestyle_westernized") # is either 1 or 0 for all samples of one study
data_set_family_bin <- data_set_family_bin[,cols_to_keep]

# remove near zero variance predictors
# data_set_family_bin_cleaned <- data_set_family_bin %>% select(which(colSums(.) <= 0.99 * nrow(.) | # remove those which are present almost everywhere
#                                                                       colSums(.) > 1.5 * nrow(.))) %>% # but keep the observed taxa number
#   select(which(colSums(.) >= 0.01 * nrow(.))) 


print("remove studies from binary family data")
tic()
if(dataset_cv_family_bin_step){
  splits <- group_vfold_cv(metadata_df, group = "study", v = length(unique(metadata_df$study)))
  caret_split <- rsample2caret(splits)
  tc_grouped <- trainControl(method = "cv",
                             savePredictions = "final",
                             selectionFunction = "oneSE", # function to select the best model
                             index = caret_split$index,
                             indexOut = caret_split$indexOut,
                             allowParallel = T)
  res_family_study_bin <- caretList(y=metadata_df$age, x=data_set_family_bin,
                                metric="Rsquared",
                                trControl=tc_grouped,
                                methodList=c("lasso", "glmnet"),
                                preProcess = c("nzv"),
                                tuneList = list(
                                  enet=caretModelSpec(method="enet", verbose = F),
                                  gbm=caretModelSpec(method="gbm", verbose = F, tuneGrid=expand.grid(shrinkage=c(0.1, 0.2, 0.3),
                                                                                                     n.minobsinnode=c(5, 15),
                                                                                                     n.trees=c(100, 200),
                                                                                                     interaction.depth=c(2, 3, 4))),
                                  kknn=caretModelSpec(method = 'kknn', tuneGrid=expand.grid(kmax = c(5, 9, 13, 15, 17),
                                                                                            distance = c(1, 5, 9),
                                                                                            kernel = "optimal")),
                                  # xgb=caretModelSpec(method="xgbTree", tuneGrid=expand.grid(nrounds = seq(from = 200, to = 1000, by = 50),
                                  #                                                           eta = c(0.025, 0.05, 0.1, 0.3),
                                  #                                                           max_depth = c(2, 3, 4, 5, 6),
                                  #                                                           gamma = 0,
                                  #                                                           colsample_bytree = 1,
                                  #                                                           min_child_weight = 1,
                                  #                                                           subsample = 1)),
                                  rf1=caretModelSpec(method=adapt_ranger, tuneGrid=expand.grid(mtry = c(3, 10, 20),
                                                                                               num.trees = c(5, 20, 50, 200, 400),
                                                                                               min.node.size = c(5, 10, 20),
                                                                                               splitrule = c("variance"),
                                                                                               replace = F),
                                                     importance = "permutation")
                                )
  )
  plots_family_study_bin <- modelplots(res_family_study_bin, metadata_df = metadata_df, all_studies = set_list$all, color = "study")
  
  save(res_family_study_bin, plots_family_study_bin,
       file = paste0("data/", merged_set, "_regression_models_dataset_bin_family_cv.RData"))
} else {
  load(paste0("data/", merged_set, "_regression_models_dataset_bin_family_cv.RData"))
}
toc()

################################################################################
# compositional family data without known abx samples
################################################################################

ps_object_family_comp_noabx <- ps_object_family_raw %>%
  subset_samples(., is.na(antibiotics_before) | !antibiotics_before) %>%
  subset_samples(., is.na(antibiotics_one_week_before) | !antibiotics_one_week_before) %>%
  subset_samples(., is.na(antibiotics_any) | !antibiotics_any) %>%
  microbiome::transform(transform = "compositional") %>%
  filter_taxa(function(x) mean(x[x > 0]) > 5e-5, TRUE) %>%  # mean abundance cutoff, only counting abundances > 0
  prune_samples(samples = (sample_sums(.) != 0))

# remove single sample subjects
multiple_samples <- ps_object_family_comp_noabx@sam_data %>%
  data.frame() %>%
  filter(duplicated(subject_ID)) %>%
  select(subject_ID) %>%
  unique %>%
  deframe
ps_object_family_comp_noabx <- ps_object_family_comp_noabx %>% subset_samples(subject_ID %in% multiple_samples)
ps_object_family_raw_noabx <- ps_object_family_raw %>% 
  prune_taxa(taxa_names(ps_object_family_comp_noabx), .) %>%
  prune_samples(sample_names(ps_object_family_comp_noabx),.) %>%
  subset_samples(subject_ID %in% multiple_samples)


# add richness and shannon diversity to the metadata
diversities <- estimate_richness(ps_object_family_raw_noabx, measures = c("Observed", "Shannon"))
sample_data(ps_object_family_comp_noabx)$observed <- diversities$Observed
sample_data(ps_object_family_comp_noabx)$shannon <- diversities$Shannon
metadata_df_noabx <- sample_data(ps_object_family_comp_noabx) %>% data.frame()
metadata_df_noabx$lifestyle_westernized <- ifelse(metadata_df_noabx$lifestyle == "westernized", yes = 1, no = 0)

data_set_family_noabx <- ps_object_family_comp_noabx%>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  # select(!Unknown) %>%
  cbind(metadata_df_noabx[,c("observed", "shannon", "lifestyle_westernized", "study")])

# prevalence filtering, only keep taxa present in all studies
cols_to_keep <- data_set_family_noabx %>% group_by(study) %>%
  summarize(across(everything(), ~ sum(.))) %>% 
  select(-study) %>%
  select(where(~ sum(. == 0) <= 1)) %>% # adjust threshold, in how many studies can a taxon be missing
  colnames() %>%
  c(., "lifestyle_westernized") # is either 1 or 0 for all samples of one study
data_set_family_noabx <- data_set_family_noabx[,cols_to_keep]


# remove near zero variance predictors
# nzv_family=preProcess(data_set_family,method="nzv",uniqueCut = 5)
# data_set_family_cleaned <- predict(nzv_family,data_set_family) %>%
#   mutate(Eubacteriaceae = X.Eubacterium..coprostanoligenes.group, .keep = "unused")
# family_cleaned_cor <- cor(data_set_family_cleaned)
# highlyCorrelated_family_cleaned <- findCorrelation(family_cleaned_cor, cutoff=0.5)
if(any(rowSums(data_set_family_noabx) == 0)) {
  warning("One or more samples do not contain any abundance for any of the selected taxa")
}


print("remove abx samples from family data")
tic()
if(dataset_cv_family_noabx_step){
  splits <- group_vfold_cv(metadata_df_noabx, group = "study", v = length(unique(metadata_df_noabx$study)))
  caret_split <- rsample2caret(splits)
  tc_grouped <- trainControl(method = "cv",
                             savePredictions = "final",
                             selectionFunction = "oneSE", # function to select the best model
                             index = caret_split$index,
                             indexOut = caret_split$indexOut,
                             preProcOptions = list(uniqueCut = 5),
                             allowParallel = T)
  res_family_study_noabx <- caretList(y=metadata_df_noabx$age, x=data_set_family_noabx,
                                metric="Rsquared",
                                trControl=tc_grouped,
                                preProcess = c("nzv"),
                                methodList=c("lasso", "glmnet"),
                                tuneList = list(
                                  enet=caretModelSpec(method="enet", verbose = F),
                                  gbm=caretModelSpec(method="gbm", verbose = F, tuneGrid=expand.grid(shrinkage=c(0.1, 0.2, 0.3),
                                                                                                     n.minobsinnode=c(5, 15),
                                                                                                     n.trees=c(100, 200),
                                                                                                     interaction.depth=c(2, 3, 4))),
                                  kknn=caretModelSpec(method = 'kknn', tuneGrid=expand.grid(kmax = c(5, 9, 13, 15, 17),
                                                                                            distance = c(1, 5, 9),
                                                                                            kernel = "optimal")),
                                  # xgb=caretModelSpec(method="xgbTree", tuneGrid=expand.grid(nrounds = seq(from = 200, to = 1000, by = 50),
                                  #                                                           eta = c(0.025, 0.05, 0.1, 0.3),
                                  #                                                           max_depth = c(2, 3, 4, 5, 6),
                                  #                                                           gamma = 0,
                                  #                                                           colsample_bytree = 1,
                                  #                                                           min_child_weight = 1,
                                  #                                                           subsample = 1)),
                                  rf1=caretModelSpec(method=adapt_ranger, tuneGrid=expand.grid(mtry = c(3, 10, 20),
                                                                                               num.trees = c(20, 50, 200, 400),
                                                                                               min.node.size = c(5, 10, 20),
                                                                                               splitrule = c("variance"),
                                                                                               replace = F),
                                                     importance = "permutation")
                                )
  )
  plots_family_study_noabx <- modelplots(res_family_study_noabx, metadata_df = metadata_df_noabx, all_studies = set_list$all, color = "study")
  
  save(res_family_study_noabx, plots_family_study_noabx,
       file = paste0("data/", merged_set, "_regression_models_dataset_family_cv_noabx.RData"))
} else {
  load(paste0("data/", merged_set, "_regression_models_dataset_family_cv_noabx.RData"))
}
toc()

# evaluate abx samples:
ps_object_family_comp_with_abx <- ps_object_family_raw %>%
  subset_samples(antibiotics_before | antibiotics_one_week_before | antibiotics_any)



################################################################################
# family data for westernized populaitons only specific model
################################################################################

ps_object_family_comp_western <- ps_object_family_raw %>%
  subset_samples(lifestyle == "westernized") %>%
  microbiome::transform(transform = "compositional") %>%
  filter_taxa(function(x) mean(x[x > 0]) > 5e-5, TRUE) %>%  # mean abundance cutoff, only counting abundances > 0
  prune_samples(samples = (sample_sums(.) != 0))

# remove single sample subjects
multiple_samples <- ps_object_family_comp_western@sam_data %>%
  data.frame() %>%
  filter(duplicated(subject_ID)) %>%
  select(subject_ID) %>%
  unique %>%
  deframe
ps_object_family_comp_western <- ps_object_family_comp_western %>% subset_samples(subject_ID %in% multiple_samples)
ps_object_family_raw_western <- ps_object_family_raw %>% 
  prune_taxa(taxa_names(ps_object_family_comp_western), .) %>%
  prune_samples(sample_names(ps_object_family_comp_western),.) %>%
  subset_samples(subject_ID %in% multiple_samples)


# add richness and shannon diversity to the metadata
diversities <- estimate_richness(ps_object_family_raw_western, measures = c("Observed", "Shannon"))
sample_data(ps_object_family_comp_western)$observed <- diversities$Observed
sample_data(ps_object_family_comp_western)$shannon <- diversities$Shannon
metadata_df_western <- sample_data(ps_object_family_comp_western) %>% data.frame()

data_set_family_western <- ps_object_family_comp_western%>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  # select(!Unknown) %>%
  cbind(metadata_df_western[,c("observed", "shannon", "study")])

# prevalence filtering, only keep taxa present in all studies
cols_to_keep <- data_set_family_western %>% group_by(study) %>%
  summarize(across(everything(), ~ sum(.))) %>% 
  select(-study) %>%
  select(where(~ sum(. == 0) <= 1)) %>% # adjust threshold, in how many studies can a taxon be missing
  colnames()
data_set_family_western <- data_set_family_western[,cols_to_keep]

if(any(rowSums(data_set_family_western) == 0)) {
  warning("One or more samples do not contain any abundance for any of the selected taxa")
}

print("family data only for western samples")
tic()
if(dataset_cv_family_western_step){
  splits <- group_vfold_cv(metadata_df_western, group = "study", v = length(unique(metadata_df_western$study)))
  caret_split <- rsample2caret(splits)
  tc_grouped <- trainControl(method = "cv",
                             savePredictions = "final",
                             selectionFunction = "oneSE", # function to select the best model
                             index = caret_split$index,
                             indexOut = caret_split$indexOut,
                             preProcOptions = list(uniqueCut = 5),
                             allowParallel = T)
  res_family_study_western <- caretList(y=metadata_df_western$age, x=data_set_family_western,
                                      metric="Rsquared",
                                      trControl=tc_grouped,
                                      preProcess = c("nzv"),
                                      methodList=c("lasso", "glmnet"),
                                      tuneList = list(
                                        enet=caretModelSpec(method="enet", verbose = F),
                                        gbm=caretModelSpec(method="gbm", verbose = F, tuneGrid=expand.grid(shrinkage=c(0.1, 0.2, 0.3),
                                                                                                           n.minobsinnode=c(5, 15),
                                                                                                           n.trees=c(100, 200),
                                                                                                           interaction.depth=c(2, 3, 4))),
                                        kknn=caretModelSpec(method = 'kknn', tuneGrid=expand.grid(kmax = c(5, 9, 13, 15, 17),
                                                                                                  distance = c(1, 5, 9),
                                                                                                  kernel = "optimal")),
                                        # xgb=caretModelSpec(method="xgbTree", tuneGrid=expand.grid(nrounds = seq(from = 200, to = 1000, by = 50),
                                        #                                                           eta = c(0.025, 0.05, 0.1, 0.3),
                                        #                                                           max_depth = c(2, 3, 4, 5, 6),
                                        #                                                           gamma = 0,
                                        #                                                           colsample_bytree = 1,
                                        #                                                           min_child_weight = 1,
                                        #                                                           subsample = 1)),
                                        rf1=caretModelSpec(method=adapt_ranger, tuneGrid=expand.grid(mtry = c(3, 10, 20),
                                                                                                     num.trees = c(20, 50, 200, 400),
                                                                                                     min.node.size = c(5, 10, 20),
                                                                                                     splitrule = c("variance"),
                                                                                                     replace = F),
                                                           importance = "permutation")
                                      )
  )
  plots_family_study_western <- modelplots(res_family_study_western, metadata_df = metadata_df_western, all_studies = set_list$all, color = "study")
  
  save(res_family_study_western, plots_family_study_western,
       file = paste0("data/", merged_set, "_regression_models_dataset_family_cv_western.RData"))
} else {
  load(paste0("data/", merged_set, "_regression_models_dataset_family_cv_western.RData"))
}
toc()



################################################################################
# family data for non westernized populaitons only
################################################################################

ps_object_family_comp_nonwestern <- ps_object_family_raw %>%
  subset_samples(lifestyle != "westernized") %>%
  microbiome::transform(transform = "compositional") %>%
  filter_taxa(function(x) mean(x[x > 0]) > 5e-5, TRUE) %>%  # mean abundance cutoff, only counting abundances > 0
  prune_samples(samples = (sample_sums(.) != 0))

# remove single sample subjects
multiple_samples <- ps_object_family_comp_nonwestern@sam_data %>%
  data.frame() %>%
  filter(duplicated(subject_ID)) %>%
  select(subject_ID) %>%
  unique %>%
  deframe
ps_object_family_comp_nonwestern <- ps_object_family_comp_nonwestern %>% subset_samples(subject_ID %in% multiple_samples)
ps_object_family_raw_nonwestern <- ps_object_family_raw %>% 
  prune_taxa(taxa_names(ps_object_family_comp_nonwestern), .) %>%
  prune_samples(sample_names(ps_object_family_comp_nonwestern),.) %>%
  subset_samples(subject_ID %in% multiple_samples)


# add richness and shannon diversity to the metadata
diversities <- estimate_richness(ps_object_family_raw_nonwestern, measures = c("Observed", "Shannon"))
sample_data(ps_object_family_comp_nonwestern)$observed <- diversities$Observed
sample_data(ps_object_family_comp_nonwestern)$shannon <- diversities$Shannon
metadata_df_nonwestern <- sample_data(ps_object_family_comp_nonwestern) %>% data.frame()

data_set_family_nonwestern <- ps_object_family_comp_nonwestern%>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  # select(!Unknown) %>%
  cbind(metadata_df_nonwestern[,c("observed", "shannon", "study")])

# prevalence filtering, only keep taxa present in all studies
cols_to_keep <- data_set_family_nonwestern %>% group_by(study) %>%
  summarize(across(everything(), ~ sum(.))) %>% 
  select(-study) %>%
  select(where(~ sum(. == 0) <= 1)) %>% # adjust threshold, in how many studies can a taxon be missing
  colnames()
data_set_family_nonwestern <- data_set_family_nonwestern[,cols_to_keep]

if(any(rowSums(data_set_family_nonwestern) == 0)) {
  warning("One or more samples do not contain any abundance for any of the selected taxa")
}

print("family data only for non-western samples")
tic()
if(dataset_cv_family_nonwestern_step){
  splits <- group_vfold_cv(metadata_df_nonwestern, group = "study", v = length(unique(metadata_df_nonwestern$study)))
  caret_split <- rsample2caret(splits)
  tc_grouped <- trainControl(method = "cv",
                             savePredictions = "final",
                             selectionFunction = "oneSE", # function to select the best model
                             index = caret_split$index,
                             indexOut = caret_split$indexOut,
                             preProcOptions = list(uniqueCut = 5),
                             allowParallel = T)
  res_family_study_nonwestern <- caretList(y=metadata_df_nonwestern$age, x=data_set_family_nonwestern,
                                        metric="Rsquared",
                                        trControl=tc_grouped,
                                        preProcess = c("nzv"),
                                        methodList=c("lasso", "glmnet"),
                                        tuneList = list(
                                          enet=caretModelSpec(method="enet", verbose = F),
                                          gbm=caretModelSpec(method="gbm", verbose = F, tuneGrid=expand.grid(shrinkage=c(0.1, 0.2, 0.3),
                                                                                                             n.minobsinnode=c(5, 15),
                                                                                                             n.trees=c(100, 200),
                                                                                                             interaction.depth=c(2, 3, 4))),
                                          kknn=caretModelSpec(method = 'kknn', tuneGrid=expand.grid(kmax = c(5, 9, 13, 15, 17),
                                                                                                    distance = c(1, 5, 9),
                                                                                                    kernel = "optimal")),
                                          # xgb=caretModelSpec(method="xgbTree", tuneGrid=expand.grid(nrounds = seq(from = 200, to = 1000, by = 50),
                                          #                                                           eta = c(0.025, 0.05, 0.1, 0.3),
                                          #                                                           max_depth = c(2, 3, 4, 5, 6),
                                          #                                                           gamma = 0,
                                          #                                                           colsample_bytree = 1,
                                          #                                                           min_child_weight = 1,
                                          #                                                           subsample = 1)),
                                          rf1=caretModelSpec(method=adapt_ranger, tuneGrid=expand.grid(mtry = c(3, 10, 20),
                                                                                                       num.trees = c(20, 50, 200, 400),
                                                                                                       min.node.size = c(5, 10, 20),
                                                                                                       splitrule = c("variance"),
                                                                                                       replace = F),
                                                             importance = "permutation")
                                        )
  )
  plots_family_study_nonwestern <- modelplots(res_family_study_nonwestern, metadata_df = metadata_df_nonwestern, all_studies = set_list$all, color = "study")
  
  save(res_family_study_nonwestern, plots_family_study_nonwestern,
       file = paste0("data/", merged_set, "_regression_models_dataset_family_cv_nonwestern.RData"))
} else {
  load(paste0("data/", merged_set, "_regression_models_dataset_family_cv_nonwestern.RData"))
}
toc()




##########
# family data with clr transformation
##########
# convert to total counts:
data_set_family_clr <- ps_object_family_raw %>%
  microbiome::transform(transform = "rclr") %>%
  filter_taxa(function(x){sum(x > 0) > 5}, TRUE) %>% # prevalence cutoff 10 % of the smallest study
  subset_samples(subject_ID %in% metadata_df$subject_ID) %>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  select(!Unknown) %>%
  cbind(metadata_df[,c("observed", "shannon", "lifestyle_westernized", "study")])

# prevalence filtering, only keep taxa present in all studies
cols_to_keep <- data_set_family_clr %>% group_by(study) %>%
  summarize(across(everything(), ~ sum(.))) %>% 
  select(-study) %>%
  select(where(~ sum(. == 0) <= 1)) %>% # adjust threshold, in how many studies can a taxon be missing
  colnames() %>%
  c(., "lifestyle_westernized") # is either 1 or 0 for all samples of one study
data_set_family_clr <- data_set_family_clr[,cols_to_keep]

# nzv_family_clr=preProcess(data_set_family_clr,method="nzv",uniqueCut = 5)
# data_set_family_clr_cleaned=predict(nzv_family_clr,data_set_family_clr) %>%
#   mutate(Eubacteriaceae = X.Eubacterium..coprostanoligenes.group, .keep = "unused")


data_set_family_clr <- data_set_family_clr %>%
  mutate(Eubacteriaceae = X.Eubacterium..coprostanoligenes.group, .keep = "unused")

# remove datasets for cv
print("remove datasets from family clr data")
tic()
if(dataset_cv_family_clr_step){
  # remove subjects for cv
  set.seed(100)
  splits <- group_vfold_cv(metadata_df, group = "study", v = length(unique(metadata_df$study)))
  caret_split <- rsample2caret(splits)
  tc_grouped <- trainControl(method = "cv",
                             savePredictions = "final",
                             selectionFunction = "oneSE", # function to select the best model
                             index = caret_split$index,
                             indexOut = caret_split$indexOut,
                             preProcOptions = list(uniqueCut = 5),
                             allowParallel = T)
  res_family_clr_study <- caretList(y=metadata_df$age, x=data_set_family_clr,
                                    metric="Rsquared",
                                    trControl=tc_grouped,
                                    preProcess = c("nzv"),
                                    methodList=c("lasso", "glmnet"),
                                    tuneList = list(
                                      enet=caretModelSpec(method="enet", verbose = F),
                                      gbm=caretModelSpec(method="gbm", verbose = F, tuneGrid=expand.grid(shrinkage=c(0.1, 0.2, 0.3),
                                                                                                              n.minobsinnode=c(5, 15),
                                                                                                              n.trees=c(100, 200),
                                                                                                              interaction.depth=c(2, 3, 4))),
                                      kknn=caretModelSpec(method = 'kknn', tuneGrid=expand.grid(kmax = c(5, 9, 13, 15, 17),
                                                                                                distance = c(1, 5, 9),
                                                                                                kernel = "optimal")),
                                      # xgb=caretModelSpec(method="xgbTree", tuneGrid=expand.grid(nrounds = seq(from = 200, to = 1000, by = 50),
                                      #                                                           eta = c(0.025, 0.05, 0.1, 0.3),
                                      #                                                           max_depth = c(2, 3, 4, 5, 6),
                                      #                                                           gamma = 0,
                                      #                                                           colsample_bytree = 1,
                                      #                                                           min_child_weight = 1,
                                      #                                                           subsample = 1)),
                                      rf1=caretModelSpec(method=adapt_ranger, tuneGrid=expand.grid(mtry = c(3, 10, 20),
                                                                                                        num.trees = c(20, 50, 400),
                                                                                                        min.node.size = c(5, 10),
                                                                                                        splitrule = c("variance", "extratrees"),
                                                                                                        replace = F),
                                                              importance = "permutation")
                                         )
  )
  plots_family_clr_study <- modelplots(res_family_clr_study, metadata_df = metadata_df, all_studies = set_list$all, color = "study")
  
  save(res_family_clr_study, plots_family_clr_study,
       file = paste0("data/", merged_set, "_regression_models_dataset_clr_family_cv.RData"))
} else {
  load(paste0("data/", merged_set, "_regression_models_dataset_clr_family_cv.RData"))
}
toc()

p1 <- ggplot(res_family_clr_study$rf1$results,# %>% filter(replace),
       aes(x=num.trees,y=Rsquared,col=factor(mtry), shape = splitrule)) +
  geom_line() + geom_point() + facet_grid(~min.node.size)
p2 <- ggplot(res_family_clr_study$rf1$results,# %>% filter(!replace),
       aes(x=num.trees,y=RMSE,col=factor(mtry), shape = splitrule)) +
  geom_line() + geom_point() + facet_grid(~min.node.size)
grid.arrange(p1, p2)

# track dynamics over time
importances <- varImp(res_family_clr_study$rf1)
top_fams <- importances$importance %>% arrange(-Overall) %>% rownames(.) %>% head(10)
top_fams <- top_fams[top_fams != "lifestyle_westernized"]
df_to_plot <- merge(metadata_df[,c("age", "study")], data_set_family_clr[,c(top_fams, "lifestyle_westernized")], by = 0) %>%
  pivot_longer(all_of(top_fams))
df_to_plot <- df_to_plot %>% filter(lifestyle_westernized == 0)
ggplot(df_to_plot, aes(x = age, y = value, color = study)) +
  geom_point() + 
  facet_wrap(~name)


# xyplot(resamples(res_family_clr_study[c(3,4)]))
# modelCor(resamples(res_family_clr_study))

# ans = resamples(res_family_clr_study) #resamples helps to tabularize the results
# summary(ans)
# dotplot(ans)
# stack.glm = caretStack(res_family_clr_study[c(3, 4)], method="ridge", trControl=trainControl(method = "cv", number = 10, 
#                                                                                              savePredictions = "final",
#                                                                                              selectionFunction = "oneSE", # function to select the best model
#                                                                                              allowParallel = T), metric = "Rsquared") #logistic
# stack.glm = caretStack(res_family_clr_study[c(4)], method="lm", trControl=trainControl(method = "cv", number = 10, 
#                                                                                        savePredictions = "final",
#                                                                                        selectionFunction = "oneSE", # function to select the best model
#                                                                                        allowParallel = T), metric = "Rsquared") #logistic
# ensembl <- caretEnsemble(res_family_clr_study[c(3, 4)], trControl=tc_grouped, metric = "Rsquared")
# stack_plots <- modelplots(stack.glm$ens_model, metadata_df = metadata_df, all_studies = set_list$all, color = "study")
# ensembl_plots <- modelplots(ensembl$ens_model, metadata_df = metadata_df, all_studies = set_list$all, color = "study")
# rf_plots <- modelplots(res_family_clr_study$rf1, metadata_df = metadata_df, all_studies = set_list$all, color = "study")
# grid.arrange(rf_plots$linear_model,stack_plots$linear_model, ensemble_plots$linear_model)
# cor(res_family_clr_study$rf1$pred$pred, res_family_clr_study$rf1$pred$obs, method = 'spearman')^2
# cor(stack.glm$ens_model$pred$pred, stack.glm$ens_model$pred$obs, method = 'spearman')^2
# cor(ensembl$ens_model$pred$pred, ensembl$ens_model$pred$obs, method = 'spearman')^2

# ############
# # family data with taxa ratio
# ############
print("family data with taxa ratio")

data_set_family <- ps_object_family_comp%>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  select(!Unknown) %>%
  cbind(metadata_df[,c("observed", "study")]) %>%
  select(-observed) %>% # doesn't work otherwise to add study information
  mutate(Eubacteriaceae = X.Eubacterium..coprostanoligenes.group, .keep = "unused")

# prevalence filtering, only keep taxa present in all studies
cols_to_keep <- data_set_family %>% group_by(study) %>%
  summarize(across(everything(), ~ sum(.))) %>% 
  select(-study) %>%
  select(where(~ sum(. == 0) <= 1)) %>% # adjust threshold, in how many studies can a taxon be missing
  colnames()
data_set_family <- data_set_family[,cols_to_keep]

ratio_names <- combn(names(data_set_family), m = 2, paste, collapse = "_")
data_set_family_ratios <- `colnames<-`(combn(data_set_family, m = 2, FUN = function(x) x[[1]]/x[[2]]),ratio_names)%>%
  replace(is.na(.), 0) %>%
  replace(. == Inf, 0) %>%
  cbind(metadata_df[,c("observed", "shannon", "lifestyle_westernized")]) %>%
  select_if(~ sum(.) != 0)

# remove near zero variance predictors
# nzv_family_ratios=preProcess(data_set_family_ratios,method="nzv",uniqueCut = 5)
# data_set_family_ratios_cleaned=predict(nzv_family_ratios,data_set_family_ratios)

print("remove studies from taxa ratio family data")
tic()
if(dataset_cv_family_ratio_step){
  splits <- group_vfold_cv(metadata_df, group = "study", v = length(unique(metadata_df$study)))
  caret_split <- rsample2caret(splits)
  tc_grouped <- trainControl(method = "cv",
                             savePredictions = "final",
                             selectionFunction = "oneSE", # function to select the best model
                             index = caret_split$index,
                             indexOut = caret_split$indexOut,
                             preProcOptions = list(uniqueCut = 5),
                             allowParallel = T)
  res_family_ratio_study <- caretList(y=metadata_df$age, x=data_set_family_ratios,
                                      metric="Rsquared",
                                      trControl=tc_grouped,
                                      preProcess = c("nzv"),
                                      methodList=c("lasso", "glmnet"),
                                      tuneList = list(
                                        enet=caretModelSpec(method="enet", verbose = F),
                                        gbm=caretModelSpec(method="gbm", verbose = F, tuneGrid=expand.grid(shrinkage=c(0.1, 0.2, 0.3),
                                                                                                     n.minobsinnode=c(5, 15),
                                                                                                     n.trees=c(100, 200),
                                                                                                     interaction.depth=c(2, 3, 4))),
                                        kknn=caretModelSpec(method = 'kknn', tuneGrid=expand.grid(kmax = c(5, 9, 13, 15, 17),
                                                                                                  distance = c(1, 5, 9),
                                                                                                  kernel = "optimal")),
                                        # xgb=caretModelSpec(method="xgbTree", tuneGrid=expand.grid(nrounds = seq(from = 200, to = 1000, by = 50),
                                  #                                                           eta = c(0.025, 0.05, 0.1, 0.3),
                                  #                                                           max_depth = c(2, 3, 4, 5, 6),
                                  #                                                           gamma = 0,
                                  #                                                           colsample_bytree = 1,
                                  #                                                           min_child_weight = 1,
                                  #                                                           subsample = 1)),
                                  rf1=caretModelSpec(method=adapt_ranger, tuneGrid=expand.grid(mtry = c(3, 10, 20),
                                                                                               num.trees = c(20, 50, 400),
                                                                                               min.node.size = c(5, 10),
                                                                                               splitrule = c("variance"),
                                                                                               replace = F),
                                                     importance = "permutation")
                                )
  )
  plots_family_ratio_study <- modelplots(res_family_ratio_study, metadata_df = metadata_df, all_studies = set_list$all, color = "study")

  save(res_family_ratio_study, plots_family_ratio_study,
       file = paste0("data/", merged_set, "_regression_models_dataset_ratio_family_cv.RData"))
} else {
  load(paste0("data/", merged_set, "_regression_models_dataset_ratio_family_cv.RData"))
}
toc()

# track dynamics over time
importances <- varImp(res_family_ratio_study$rf1)
top_fams <- importances$importance %>% arrange(-Overall) %>% rownames(.) %>% head(10)
df_to_plot <- merge(metadata_df[,c("age", "study")], data_set_family_ratios[,top_fams], by = 0) %>%
  pivot_longer(all_of(top_fams)) %>%
  filter(!name %in% c("observed", "shannon"))

ggplot(df_to_plot, aes(x = age, y = value, color = study)) +
  geom_point() + 
  facet_wrap(~name)



# p1 <- ggplot(res_family_ratio_study$rf1$results %>% filter(replace),
#        aes(x=num.trees,y=RMSE,col=factor(mtry), shape = splitrule)) +
#   geom_line() + geom_point() + facet_grid(~min.node.size) +
#   ggtitle("replace = T")
# p2 <- ggplot(res_family_ratio_study$rf1$results %>% filter(!replace),
#        aes(x=num.trees,y=RMSE,col=factor(mtry), shape = splitrule)) +
#   geom_line() + geom_point() + facet_grid(~min.node.size) +
#   ggtitle("replace = F")
# gridExtra::grid.arrange(p1, p2)


# xyplot(resamples(res_family_ratio_study[c(3,4)]))
# modelCor(resamples(res_family_ratio_study))
# 
# ans = resamples(res_family_ratio_study) #resamples helps to tabularize the results
# # summary(ans)
# # dotplot(ans)
# stack.glm = caretStack(res_family_ratio_study[c(3, 4)], method="ridge", trControl=trainControl(method = "cv", number = 10,
#                                                                                              savePredictions = "final",
#                                                                                              selectionFunction = "oneSE", # function to select the best model
#                                                                                              allowParallel = T),
#                        metric = "Rsquared") #logistic
# stack.glm = caretStack(res_family_ratio_study[c(4)], method="lm", trControl=trainControl(method = "cv", number = 10,
#                                                                                        savePredictions = "final",
#                                                                                        selectionFunction = "oneSE", # function to select the best model
#                                                                                        allowParallel = T),
#                        metric = "Rsquared") #logistic
# ensembl <- caretEnsemble(res_family_ratio_study[c(3, 4)], trControl= trainControl(method = "cv", number = 10,
#                                                                                   savePredictions = "final",
#                                                                                   selectionFunction = "oneSE", # function to select the best model
#                                                                                   allowParallel = T),
#                          metric = "Rsquared")
# stack_plots <- modelplots(stack.glm$ens_model, metadata_df = metadata_df, all_studies = set_list$all, color = "study")
# ensembl_plots <- modelplots(ensembl$ens_model, metadata_df = metadata_df, all_studies = set_list$all, color = "study")
# rf_plots <- modelplots(res_family_ratio_study$rf1, metadata_df = metadata_df, all_studies = set_list$all, color = "study")
# grid.arrange(rf_plots$linear_model,stack_plots$linear_model, ensemble_plots$linear_model)
# cor(res_family_ratio_study$rf1$pred$pred, res_family_ratio_study$rf1$pred$obs, method = 'spearman')^2
# cor(stack.glm$ens_model$pred$pred, stack.glm$ens_model$pred$obs, method = 'spearman')^2
# cor(ensembl$ens_model$pred$pred, ensembl$ens_model$pred$obs, method = 'spearman')^2


############
# genus data
############
print("genus data")
ps_object_genus_raw <- readRDS("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_phyloseq_rf_filter_genus.rds") %>%
  subset_samples(., age <= 730 & age > 1)


ps_object_genus_comp <- ps_object_genus_raw %>%
  microbiome::transform(transform = "compositional") %>%
  filter_taxa(function(x){sum(x > 0) > 5}, TRUE) %>% # prevalence cutoff 10 % of the smallest study
  filter_taxa(function(x) mean(x[x > 0]) > 5e-5, TRUE) %>%  # mean abundance cutoff, only counting abundances > 0
  prune_samples(samples = (sample_sums(.) != 0))

ps_object_genus_raw <- ps_object_genus_raw %>% 
  prune_taxa(taxa_names(ps_object_genus_comp), .) %>%
  prune_samples(sample_names(ps_object_genus_comp),.)


# add richness and shannon diversity to the metadata
diversities <- estimate_richness(ps_object_genus_raw, measures = c("Observed", "Shannon"))
sample_data(ps_object_genus_comp)$observed <- diversities$Observed
sample_data(ps_object_genus_comp)$shannon <- diversities$Shannon
metadata_df_genus <- sample_data(ps_object_genus_comp) %>% data.frame()
metadata_df_genus$lifestyle_westernized <- ifelse(metadata_df_genus$lifestyle == "westernized", yes = 1, no = 0)


data_set_genus <- ps_object_genus_comp%>%
  otu_table() %>%
  t() %>%
  data.frame() %>%
  select(!Unknown) %>%
  cbind(metadata_df_genus[,c("observed", "shannon", "lifestyle_westernized", "study")])

# prevalence filtering, only keep taxa present in all studies
gena_to_keep <- data_set_genus %>% group_by(study) %>%
  summarize(across(everything(), ~ sum(.))) %>% 
  select(-study) %>%
  select(where(~ sum(. == 0) <= 1)) %>% # adjust threshold, in how many studies can a taxon be missing
  colnames() %>%
  c(., "lifestyle_westernized") # is either 1 or 0 for all samples of one study
data_set_genus <- data_set_genus[,gena_to_keep]


# remove near zero variance predictors
# nzv_genus=preProcess(data_set_genus,method="nzv",uniqueCut = 5)
# data_set_genus_cleaned=predict(nzv_genus,data_set_genus)



print("remove studies from genus data")
tic()
if(dataset_cv_genus_step){
  splits <- group_vfold_cv(metadata_df_genus, group = "study", v = length(unique(metadata_df_genus$study)))
  caret_split <- rsample2caret(splits)
  tc_grouped <- trainControl(method = "cv",
                             savePredictions = "final",
                             selectionFunction = "oneSE", # function to select the best model
                             index = caret_split$index,
                             indexOut = caret_split$indexOut,
                             preProcOptions = list(uniqueCut = 5),
                             allowParallel = T)
  res_genus_study <- caretList(y=metadata_df_genus$age, x=data_set_genus,
                               metric="Rsquared",
                               trControl=tc_grouped,
                               preProcess = c("nzv"),
                               methodList=c("lasso", "glmnet"),
                               tuneList = list(
                                 enet=caretModelSpec(method="enet", verbose = F),
                                 gbm=caretModelSpec(method="gbm", verbose = F, tuneGrid=expand.grid(shrinkage=c(0.1, 0.2, 0.3),
                                                                                                     n.minobsinnode=c(5, 15),
                                                                                                     n.trees=c(100, 200),
                                                                                                     interaction.depth=c(2, 3, 4))),
                                 kknn=caretModelSpec(method = 'kknn', tuneGrid=expand.grid(kmax = c(5, 9, 13, 15, 17),
                                                                                           distance = c(1, 5, 9),
                                                                                           kernel = "optimal")),
                                 # xgb=caretModelSpec(method="xgbTree", tuneGrid=expand.grid(nrounds = seq(from = 200, to = 1000, by = 50),
                                 #                                                           eta = c(0.025, 0.05, 0.1, 0.3),
                                 #                                                           max_depth = c(2, 3, 4, 5, 6),
                                 #                                                           gamma = 0,
                                 #                                                           colsample_bytree = 1,
                                 #                                                           min_child_weight = 1,
                                 #                                                           subsample = 1)),
                                 rf1=caretModelSpec(method=adapt_ranger, tuneGrid=expand.grid(mtry = c(3, 10, 20),
                                                                                               num.trees = c(20, 50, 400),
                                                                                               min.node.size = c(5, 10),
                                                                                               splitrule = c("variance"),
                                                                                               replace = F),
                                                     importance = "permutation")
                                )
  )
  plots_genus_study <- modelplots(res_genus_study, metadata_df = metadata_df_genus, all_studies = set_list$all, color = "study")
  
  save(res_genus_study, plots_genus_study,
       file = paste0("data/", merged_set, "_regression_models_dataset_genus_cv.RData"))
} else {
  load(paste0("data/", merged_set, "_regression_models_dataset_genus_cv.RData"))
}
toc()

p1 <- ggplot(res_genus_study$rf1$results,
       aes(x=num.trees,y=RMSE,col=factor(mtry), shape = splitrule)) +
  geom_line() + geom_point() + facet_grid(~min.node.size) +
  ggtitle("replace = T")


importances <- varImp(res_genus_study$rf1)
top_fams <- importances$importance %>% arrange(-Overall) %>% rownames(.) %>% head(10)
top_fams <- top_fams[top_fams!= "lifestyle_westernized"]
df_to_plot <- merge(metadata_df[,c("age", "study")], data_set_genus[,c(top_fams, "lifestyle_westernized")], by = 0) %>%
  pivot_longer(all_of(top_fams)) %>%
  filter(!name %in% c("observed", "shannon"))

p1 <- ggplot(df_to_plot[df_to_plot$lifestyle_westernized == 1,], aes(x = age, y = value, color = study)) +
  geom_point(alpha = 0.7) + 
  facet_wrap(~name) +
  ggtitle("westernized")
p2 <- ggplot(df_to_plot[df_to_plot$lifestyle_westernized == 0,], aes(x = age, y = value, color = study)) +
  geom_point(alpha = 0.7) + 
  facet_wrap(~name) +
  ggtitle("not-westernized")
grid.arrange(p1, p2)



# other ideas: ensemble learning: combine learners. Would only make sense, if results are not correlated
