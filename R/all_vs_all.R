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
library(gridExtra)
library(tictoc)
library(furrr)
library(gridExtra)
library(igraph)
library(tidygraph)
library(ggraph)
library(rstatix)
source("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/alt_models.R")
source("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/regression_functions.R")
source("/fast/AG_Forslund/rob/mm_index/R_scripts/setlists.R")

all_vs_all_genus <-  T

n_cores <- 8

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

cl <- makePSOCKcluster(n_cores)
registerDoParallel(cl)
getDoParWorkers()
set.seed(825)
future::plan(multisession, workers = n_cores)



get_predictions_all_vs_all <- function(train_set, ps){
  print(train_set)
  oldDF <- as(sample_data(ps), "data.frame") 
  trainDF <- subset(oldDF, study_country == train_set)

  testDF <- subset(oldDF, study_country != train_set)
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
    tc_grouped <- trainControl(method = "cv",
                               number = 10,
                             savePredictions = "final",
                             selectionFunction = "oneSE", # function to select the best model
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
    mutate(pred = predict(rf_model, newdata = for_caret_list_test$features),
           train_study = train_set)
  return(predictions)
}



# genus data #############
print("genus data")
ps_object_genus <- readRDS("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_phyloseq_rf_filter_genus.rds") %>%
  subset_samples(., age <= 730 & age > 1)
ps_object_genus@sam_data$study_country <- paste(ps_object_genus@sam_data$study, ps_object_genus@sam_data$country, sep = "_")
all_studies <- ps_object_genus@sam_data$study_country %>%
  unique
n_studies <- length(all_studies)


if(all_vs_all_genus){
  tic()
  all_vs_all_df <- all_studies %>% 
    future_map_dfr(~ get_predictions_all_vs_all(. ,ps = ps_object_genus), .progress = F)
  toc()
  save(all_vs_all_df, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_vs_all_genus.RData")
} else {
  load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_vs_all_genus.RData")
}

study_ls_dict <- all_vs_all_df %>% 
  group_by(study_country) %>%
  summarise(lifestyle = unique(lifestyle),
            country = unique(country),
            n_samples = length(unique(sample_ID)),
            age_span = abs(min(age) - max(age)))

all_vs_all_R2 <- all_vs_all_df %>%
  mutate(model_name = "_",
         train_test = paste(train_study, study_country, sep = "___")) %>%
  get_lm_list(., grouping = "train_test") %>%
  separate(., train_test, into = c("train_study", "test_study"), sep = "___") %>%
  left_join(., study_ls_dict, by = c("train_study" = "study_country")) %>%
  left_join(., study_ls_dict, by = c("test_study" = "study_country"), suffix = c("_train", "_test")) %>%
  mutate(model_name = NULL)
# export to cytoscape:
write.csv(all_vs_all_R2, "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_vs_all_genus_network.csv")
write.csv(study_ls_dict, "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_vs_all_genus_network_annotation.csv")

study_graph <- graph_from_data_frame(all_vs_all_R2[,1:3], directed = T)
V(study_graph)$country <- study_ls_dict$country[match(V(study_graph)$name, study_ls_dict$study_country)]
V(study_graph)$lifestyle <- study_ls_dict$lifestyle[match(V(study_graph)$name, study_ls_dict$study_country)]
V(study_graph)$n_samples <- study_ls_dict$n_samples[match(V(study_graph)$name, study_ls_dict$study_country)]
V(study_graph)$age_span <- study_ls_dict$age_span[match(V(study_graph)$name, study_ls_dict$study_country)]
communities <- cluster_infomap(study_graph)
V(study_graph)$community <- communities$membership

V(study_graph)$color <- ifelse(V(study_graph)$lifestyle == "westernized", yes = "red", no = "blue")
E(study_graph)$width <- E(study_graph)$Coefficient * 2

write_graph(study_graph, 
            file = "/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_vs_all_genus_network.graphml", 
            format = "graphml")


plot(study_graph, edge.curved=0.5)


# Convert to tbl_graph
tbl_graph <- as_tbl_graph(study_graph)



# Plot the graph with curved edges and community layout
ggraph(tbl_graph, layout = "fr") + 
  geom_edge_arc(aes(edge_alpha = Coefficient), show.legend = TRUE) +  # Curved edges
  geom_node_point(aes(color = lifestyle), size = 3) +
  geom_node_text(aes(label = name), vjust = 1.8, size = 5) +
  theme_void() +
  theme(legend.position = "right")  # Adjust legend position as needed



# intra vs. inter lifestyle correlation
all_vs_all_R2 %>%
  mutate(same_lifestyle = lifestyle_test == lifestyle_train,
         same_country = country_train == country_test) %>%
  ggplot(., aes(x = Coefficient, color = same_lifestyle)) +
  geom_histogram()

all_vs_all_R2 %>%
  mutate(same_lifestyle = lifestyle_test == lifestyle_train,
         same_country = country_train == country_test) %>%
  rstatix::wilcox_test(Coefficient ~ same_lifestyle) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj") %>%
  rstatix::add_xy_position()




