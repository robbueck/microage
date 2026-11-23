# microage
Scripts for the publication Global analysis of infant gut microbiota revealed distinctive maturation dynamics across lifestyles (https://doi.org/10.1038/s41522-026-01130-4)

## File description:
- examples for dada2-scripts for studies consisting of one cohort: `dada_2_single_cohort.R` or multiple cohorts: `dada_2_multiple_cohorts.R`, both require `dada_2_functions.R`
- aggregate all studies and perform rarefaction filtering with `merge_all_simple.R`
- alpha and beta diversity analysis, prevalence analyses done with `pcoas.R`, requires `pcoa_functions.R`
- main models are trained with: `all_regression_nested_cv.R`
- cross-application of lifestyle specific models with: `lifestyle_regression.R`, requires `lifestyle_regression_functions.R`
- training models on downsampled lifestyle groups with matching age distributions with: `lifestyle_regression_downsample.R`         
- all analyses including SAM and preterm infants done with: `sick_children.R`
- all random forest scripts need `regression_functions.R` and `alt_models.R` (a modified version of ranger for caret to vary selected hyperparameters)
- scripts for generating all figures in the paper: `paper_figures.R` and `paper_supplement_figures.R`


## Application of age index
Generate count table following the example dada2-Rscript using the SILVA database version 138.
Using another processing pipeline like LotuS should work as well, but might introduce some bias.

### Useful files in data/:
- industrialized model: `final_industrialized_genus.RData`
- non-industrialized model: `final_nonindustrialized_genus.RData`
- training data (also works as input data example): `all_phyloseq_raw_genus.rds`

Model predictions:
- combined model predictions (predicted age in column rf1): `all_nested_cv_genus_no_ls.RData`
- predictions from industrialized model on industrialized data: `all_nested_cv_genus_industrialized.RData`
- predictions from non-industrialized model on non-industrialized data: `all_nested_cv_genus_nonindustrialized.RData`
- Cross lifestyle predictions: `inter_lifestyle_pred_genus.RData`

```
# load functions and model
source(regression_functions.R)
ps_all <- readRDS("../data/all_phyloseq_raw_genus.rds") # example taxonomic input profile
load("../data/final_industrialized_genus.RData") # change here to load different models

# Formatting input data from a phyloseq object
filtered_cdf <- create_caret_df(ps_all, transformation = "compositional", only_multiple_samples = F,
                filter_features = names(final_model_genus_industrialized$finalModel$variable.importance),
                additional_cols = c("Observed", "Shannon"))

# predict microbiome age
predicted_age = predict(final_model_genus_industrialized$rf1, newdata = filtered_cdf$features)
```
## Get SHAP values for new data
```
shap_res_combined <- get_shap_long(final_model_genus_industrialized$rf1,
                                   test_data = filtered_cdf$features
                                   features = intersect(colnames(filtered_cdf$features),
                                                        names(final_model_genus_industrialized$finalModel$variable.importance)),
                                   normalize_ab_values = F)

```

# Comments
- I didn't upload all models, further models can either be trained with `all_regression_nested_cv.R`, or contact me, then I can share them
- I will set up a zenodo repository in the next weeks to also share the combined and a preterm specific model.
- Some intermediate files are too large to be uploaded here, but can be shared on request.
- Work in progress: models for updated SILVA taxonomy v144 and shotgun data processed by metaphlan/motus.
