setwd("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/")
# another idea: instead of abundances, use ratios between taxa
rm(list = ls())
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
library(gridExtra)
source("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/alt_models.R")
source("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/regression_functions.R")
source("/fast/AG_Forslund/rob/mm_index/R_scripts/setlists.R")

rarefaction_all_families <-  F
rarefaction_all_genus <- F

n_cores <- 5

option_list = list(
  make_option(c("-t", "--threads"), type="numeric", default=NULL))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$threads)){
  n_cores <- T
  n_cores <- 4
} else {
  n_cores <- opt$threads
}

merged_set <- "all"

cl <- makePSOCKcluster(ceiling(n_cores/2))
registerDoParallel(cl)
getDoParWorkers()
set.seed(825)
future::plan(multisession, workers = ceiling(n_cores/2))



run_rf_rf <- function(s_list, ps){
  print(s_list)
  oldDF <- as(sample_data(ps), "data.frame") 
  newDF <- subset(oldDF, study %in% s_list) 
  sample_data(ps) <- sample_data(newDF) 
  for_caret_list <- create_caret_df(ps_object = ps, transformation = "compositional",
                                    mean_ab_cutoff = 5e-5, study_prevalence_cutoff = 2,
                                    prevalence_in_study_cutoff = 5,
                                    additional_cols = c("Observed", "Shannon", "lifestyle_industrialized"),
                                    only_multiple_samples = F)
  splits <- group_vfold_cv(for_caret_list$metadata, group = "study", v = length(s_list))
  caret_split <- rsample2caret(splits)
  tc_grouped <- trainControl(method = "cv",
                             savePredictions = "final",
                             selectionFunction = "oneSE", # function to select the best model
                             index = caret_split$index,
                             indexOut = caret_split$indexOut,
                             preProcOptions = list(uniqueCut = 5),
                             allowParallel = T)
  rfGrid <- expand.grid(mtry = c(3, 10, 20),
                        num.trees = c(50, 100),
                        min.node.size = c(10, 20),
                        splitrule = c("variance"),
                        replace = F)
  # sometimes rerunning helps:
  retry <- 1
  rf_model <- NULL
  while (retry < 4) { # try three times
    tryCatch({
      cat("building model, try number: ", retry, "\n")
      retry <- retry + 1
      rf_model <- train(x = for_caret_list$features, y = for_caret_list$metadata$age,
                        method = adapt_ranger, 
                        metric="Rsquared",
                        # importance = "permutation",
                        trControl = tc_grouped,
                        preProcess = c("nzv"),
                        tuneGrid = rfGrid)
      retry <- 4 # if successful, continue
    },
    error=function(e){
      if(retry < 4){
        message("training model for following combination was not successful, try again for:")
        print(s_list)
        print(e)
      } else{
        message("training model for following combination was not successful, skip:")
        print(s_list)
        print(e)
      }
    },
    warning = function(w){
      print(w)
    }
    )
  }
  if(is.null(rf_model)) {
    res <- data.frame(study = s_list, model_name = "failed", Coefficient = NA, R2 = NA, n_studies = i,
                      model_studies = paste(s_list, sep = "xxx", collapse = "xxx"))
  } else {
    pred_meta_df <- pred_from_list(model_list = rf_model, metad = for_caret_list$metadata)
    res <- get_lm_list(pred_meta_df, grouping = "study") %>%
      mutate(n_studies = length(s_list),
             model_studies = paste(s_list, sep = "xxx", collapse = "xxx"))
  }
  return(res)
}



# family data #############
print("family data")
ps_object_family <- readRDS("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_phyloseq_rf_filter_family.rds") %>%
  subset_samples(., age <= 730 & age > 1)
# if(any(rowSums(for_caret_list_family$features) == 0)) {
#   warning("One or more samples do not contain any abundance for any of the selected taxa")
# }
all_studies <- ps_object_family@sam_data$study %>%
  unique
n_studies <- length(all_studies)

subsample_studies <- list()
for(i in 2:n_studies){
  all_combs <- combn(all_studies, m = i, simplify = F)
  if(length(all_combs) > n_studies * 3){  # limit sample size to number of studies
    all_combs <- all_combs[sample(length(all_combs), size = n_studies * 3)]
  }
  subsample_studies <- c(subsample_studies, all_combs)
}

if(rarefaction_all_families){
  tic()
  rf_rf_df <- subsample_studies %>% 
    future_map_dfr(~ run_rf_rf(. ,ps = ps_object_family), .progress = F)
  toc()
  save(rf_rf_df, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_regression_rarefaction.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_regression_rarefaction.RData")
}

study_lifestyle <- ps_object_family %>%
  sample_data() %>%
  data.frame %>%
  select(study, lifestyle) %>%
  distinct(study, .keep_all = T)

rf_rf_df <- left_join(rf_rf_df, study_lifestyle, by = "study") %>%
  arrange(-R2) %>%
  mutate(study = factor(study, levels = unique(study)))

# c("subramanian_2014", "blanton_2016", "sprockett_2020",
#   "raman_2019", "gehrig_2019", "kamngona_2019",
#   "davis_2017", "pannaraj_2017", "vatanen_2018",
#   "wampach_2018", "hill_2017", "kristensen_2020",
#   "beller_2021", "muinck_2018", "bockulich_2016",
#   "roswall_2021", "reyman_2019", "lim_2015",        
#   "stokholm_2018")
library(ggpubr)
ggplot(rf_rf_df, aes(x = n_studies, y = R2, color = lifestyle)) +
  geom_point() +
  facet_wrap(~ study) +
  geom_hline(data = plyr::ddply(rf_rf_df, c("study"),
                                summarize, final_model = R2[which.max(n_studies)]),
             aes(yintercept=final_model, linetype = "final_model")) +
  scale_linetype_manual(name = NULL, values = "solid") +
  ylab("R2 ( pred age, true age )")+
  guides(linetype = guide_legend(order = 2),color = guide_legend(order = 1)) +
  geom_smooth(method = "lm") +
  stat_cor(method = "spearman",
           # aes(colour = NULL),
           cor.coef.name = "spearman",
           size = 3,
           color = "black",
           show.legend = F)

ggsave("figs/rarefaction_all_per_study.pdf")





ggplot(rf_rf_df, aes(x = n_studies, y = R2, color = study)) +
  geom_point() +
  facet_wrap(~lifestyle)


reml_rarefy <- lmer(R2 ~ n_studies + (1|study), data = rf_rf_df)
reml_rarefy_null <- lmer(R2 ~ (1|study), data = rf_rf_df)
res_anova <- anova(reml_rarefy, reml_rarefy_null)

summary(reml_rarefy)
plot(reml_rarefy)
qqmath(reml_rarefy)
hist(resid(reml_rarefy), probability = TRUE, breaks = 50)
curve(dnorm(x, mean = mean(resid(reml_rarefy)),
            sd = sd(resid(reml_rarefy))),
      from = -0.2, to = 0.2, add= TRUE, col = "blue")
plot(reml_rarefy, form = study ~ resid(., scale = TRUE))

# lm(R2 ~ n_studies, data = rf_rf_df)


# genus data #############
print("genus data")
ps_object_genus <- readRDS("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_phyloseq_rf_filter_genus.rds") %>%
  subset_samples(., age <= 730 & age > 1)

subsample_studies <- list()
for(i in 1:(n_studies/2)*2){
  all_combs <- combn(all_studies, m = i, simplify = F)
  if(length(all_combs) > n_studies * 3){  # limit sample size to number of studies
    all_combs <- all_combs[sample(length(all_combs), size = n_studies * 3)]
  }
  subsample_studies <- c(subsample_studies, all_combs)
}

if(rarefaction_all_genus){
  tic()
  rf_rf_df_genus <- subsample_studies %>% 
    future_map_dfr(~ run_rf_rf(. ,ps = ps_object_genus), .progress = F)
  toc()
  save(rf_rf_df_genus, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_regression_rarefaction_genus.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_regression_rarefaction_genus.RData")
}

study_lifestyle_genus <- ps_object_genus %>%
  sample_data() %>%
  data.frame %>%
  select(study, lifestyle) %>%
  distinct(study, .keep_all = T)

rf_rf_df_genus <- left_join(rf_rf_df_genus, study_lifestyle_genus, by = "study") %>%
  arrange(-R2) %>%
  mutate(study = factor(study, levels = unique(study)))

# if they are so big, how come they were not described before?

save(rf_rf_df_genus, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/rarefaction_studies.RData")
ggplot(rf_rf_df_genus,# %>%
         # group_by(study, lifestyle, n_studies) %>%
         # summarize(R2 = mean(R2)),
       aes(x = as.factor(20-n_studies), y = R2, color = lifestyle)
       ) +
  geom_point() +
  # geom_boxplot() +
  facet_wrap(~ study) +
  geom_smooth(method = "lm", aes(group = study)) +
  geom_hline(data = plyr::ddply(rf_rf_df_genus, c("study"),
                                summarize, final_model = R2[which.max(n_studies)]),
             aes(yintercept=final_model, linetype = "final_model")) +
  scale_linetype_manual(name = NULL, values = "solid") +
  ylab("R2 ( pred age, true age )")+
  guides(linetype = guide_legend(order = 2),color = guide_legend(order = 1))
ggsave("figs/rarefaction_all_per_study_genus.pdf")
ggsave("figs/rarefaction_all_per_study_genus.png")

reml_rarefy_genus <- lmer(R2 ~ n_studies + (1|study), data = rf_rf_df_genus)
reml_rarefy_genus_null <- lmer(R2 ~ (1|study), data = rf_rf_df_genus)
res_anova_genus <- anova(reml_rarefy_genus, reml_rarefy_genus_null)
