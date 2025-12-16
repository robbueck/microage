library(tidyverse)
library(ggplot2)
library(microbiome)
library(lme4)
library(lmtest)
library(metadeconfoundR)
library(scales)

metavars <- c("predicted_age", "predicted_age_west", "predicted_age_nowest", "age", "antibiotics_any", "antibiotics_before", "antibiotics_one_week_before", "birth_height", 
              "birth_weight", "birthmode", "country", "exclusive_breastfeeding", "family_ID", "first_solid_food",
              "food", "gestational_age", "height", "sample_ID", "sex", "study", "subject_ID", "lifestyle", "weight",
              "Shannon", "Observed", "WHZ", "WAZ", "HAZ")

load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_no_ls.RData")
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_western.RData")
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/all_nested_cv_genus_nonwestern.RData")
# get MAZ and metadata
pred_metadata_df <- genus_no_ls_nested_cv_preds %>% 
  left_join(.,genus_western_nested_cv_preds %>% select(run_accession, rf1) %>% mutate(predicted_age_west = rf1, .keep = "unused")) %>%
  left_join(.,genus_nonwestern_nested_cv_preds %>% select(run_accession, rf1) %>% mutate(predicted_age_nowest = rf1, .keep = "unused")) %>%
  mutate(predicted_age = rf1, rf1 = NULL) %>%
  select(all_of(metavars)) %>%
  mutate(age_cat_MAZ =  cut(age, breaks=c(-1, 2, 7, 20, 1:25 * 30))) %>%
  group_by(age_cat_MAZ) %>%
  # mutate(across(where(is.numeric) & !all_of(c("predicted_age", "age", "Shannon",
  #                                             "Observed", "WHZ", "WAZ", "HAZ")), ~ scale(.)[,1])) %>%
  ungroup

# check availability of metadata
available_mdata <- pred_metadata_df %>% select(study, lifestyle, sex, birthmode, antibiotics_any, 
                            antibiotics_before, antibiotics_one_week_before, food) %>%
  dplyr::mutate(sex = gsub("female", "f", sex),
                study = paste(lifestyle, study, sep = "_")) %>%
  fastDummies::dummy_cols("food", remove_first_dummy = F, remove_selected_columns = T, ignore_na = T) %>%
  fastDummies::dummy_cols(c("sex", "birthmode"), remove_first_dummy = T, 
                          remove_selected_columns = T, ignore_na = T) %>%
  select(-food_nothing, -food_unknown, -food_, food_mixed) %>%
  dplyr::mutate(food_formula = max(food_formula, food_mixed_b_f, food_mixed_b_f_s, food_mixed_f_b_s, food_mixed_f_s, na.rm = T),
         food_breast = max(food_breast, food_mixed_b_f, food_mixed_b_f_s, food_mixed_b_s, food_mixed_f_b_s, na.rm = T),
         food_solid = max(food_solid, food_mixed_b_f_s, food_mixed_b_s, food_mixed_f_b_s, food_mixed_f_s, na.rm = T),
         .keep = "unused")

mdat_summary <- bind_rows(available_mdata, available_mdata %>%
                            mutate(study = paste(lifestyle, "combined", sep = "_"))) %>% 
  select(-lifestyle) %>%
  pivot_longer(-c("study"), names_to = "variable", values_to = "value") %>%
  group_by(study, variable) %>%
  dplyr::summarise( count_1 = sum(value == 1, na.rm = TRUE),  # Count of 1s
    count_0 = sum(value == 0, na.rm = TRUE),  # Count of 0s
    non_na_count = sum(!is.na(value)) %>% log(),        # Count of non-NA values
    ratio = count_1 / (count_0 + count_1),                # Calculate ratio of 1s to 0s
    ratio = min(ratio, 1-ratio, na.rm = T),
    .groups = "drop")

mdat_summary %>% filter(non_na_count > 0, !variable %in% c("food_breast", "food_formula", "food_mixed", "food_solid")) %>%
  mutate(study = factor(study, levels = unique(c("non_westernized_combined", "westernized_combined", study)))) %>%
  ggplot(., aes(x = variable, y = study, fill = ratio, size = non_na_count)) +
  geom_point(shape = 21, color = "black", alpha = 0.5) +  # Points with black borders
  scale_fill_gradient(low = "blue", high = "red", na.value = "grey") +  # Color scale for ratio
  scale_size_continuous(range = c(3, 10)) +  # Control point size based on non-NA count
  theme_minimal() +
  labs(title = "Available_metadata",
       x = NULL,
       y = NULL,
       fill = "Ratio",
       size = "Non-NA Count (log)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  geom_text(aes(label = paste0(count_1, "\n", count_0)), size = 3, color = "black")  # Add text labels


age_mean_sd <- pred_metadata_df %>%
  group_by(age_cat_MAZ) %>%
  summarize(group_median = median(predicted_age, na.rm = T),
            group_sd = sd(predicted_age, na.rm = T),
            group_median_west = median(predicted_age_west, na.rm = T),
            group_sd_west = sd(predicted_age_west, na.rm = T),
            group_median_nowest = median(predicted_age_nowest, na.rm = T),
            group_sd_nowest = sd(predicted_age_nowest, na.rm = T))


pred_metadata_df <- pred_metadata_df %>%
  left_join(., age_mean_sd, by = "age_cat_MAZ") %>%
  mutate(MAZ = (predicted_age - group_median) / group_sd,
         MAZ_west = (predicted_age_west - group_median_west) / group_sd_west,
         MAZ_nowest = (predicted_age_nowest - group_median_nowest) / group_sd_nowest) %>%
  mutate(
    # across(where(is.numeric) & !all_of(c("MAZ", "predicted_age", "group_median", "group_sd",
    #                                           "MAZ_west", "predicted_age_west", "group_median_west", "group_sd_west",
    #                                           "MAZ_nowest", "predicted_age_nowest", "group_median_nowest", "group_sd_nowest",
    #                                           "age")), ~ rescale(.)),
         lifestyle_westernized = lifestyle == "westernized") %>%
  select(-lifestyle)
  

ggplot(pred_metadata_df, aes(x = age_cat_MAZ, y = MAZ)) +
  geom_boxplot()


# does the subject/family even matter for MAZ? yes, it does
# age_model <- lmer ( (MAZ) ~   age  + (1 | family_ID) + (1 | study),
#                      REML = F, data = pred_metadata_df)
# confint(age_model, method = "boot")


# associations with MAZ:
full_model <- lmer ( (MAZ) ~   age * antibiotics_one_week_before + (1 | family_ID) + (1 | study),
                                REML = F, data = pred_metadata_df)
no_slope_model <- lmer ( (MAZ) ~   age + antibiotics_one_week_before + (1 | family_ID) + (1 | study),
                               REML = F, data = pred_metadata_df)
no_intercept_model <- lmer ( (MAZ) ~  age + age:antibiotics_one_week_before + (1 | family_ID) + (1 | study),
                                   REML = F, data = pred_metadata_df)
lrtest (full_model, no_slope_model) 
lrtest (full_model, no_intercept_model) 
pred_metadata_df %>%
  filter(!is.na(antibiotics_before)) %>% 
         # !study %in% c("wampach_2018", "sprockett_2020", "roswall_2021"),
         # study == "bockulich_2016") %>%
  ggplot(., aes(x = age, y = MAZ, color = antibiotics_before)) +
  # geom_point(alpha = 0.3) +
  geom_smooth(method = "loess") +
  facet_wrap(~study)


# Shannon vs MAZ
full_model <- lmer ( (MAZ) ~   age * Shannon + (1 | family_ID) + (1 | study),
                     REML = F, data = pred_metadata_df)
no_slope_model <- lmer ( (MAZ) ~   age + Shannon + (1 | family_ID) + (1 | study),
                         REML = F, data = pred_metadata_df)
no_intercept_model <- lmer ( (MAZ) ~  age + age:Shannon + (1 | family_ID) + (1 | study),
                             REML = F, data = pred_metadata_df)
lrtest (full_model, no_slope_model) 
lrtest (full_model, no_intercept_model) 
pred_metadata_df %>%
  filter(!is.na(Shannon)) %>%
  ggplot(., aes(x = Shannon, y = MAZ, color = Shannon)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm") +
  facet_wrap(~study)


# check metadeconfounder
ftrs <- pred_metadata_df %>% select(MAZ, MAZ_west, MAZ_nowest, Shannon, Observed)
mvars <- pred_metadata_df %>% select(-c("MAZ", "MAZ_west", "MAZ_nowest", "Shannon", "Observed", "subject_ID",
                                        "country", "predicted_age",
                                        "group_sd", "group_median",
                                        "predicted_age_west", "group_sd_west", "group_median_west",
                                        "predicted_age_nowest", "group_sd_nowest", "group_median_nowest")) %>%
  select(sex, everything())


metad_all <- metadeconfoundR::MetaDeconfound(featureMat = ftrs %>% data.frame(),
                                             metaMat = mvars %>% data.frame(),
                                             fixedVar = c("age"),
                                             returnLong = T,
                                             # collectMods = T,
                                             DCutoff = 0.1,
                                             startStop = "naiveStop",
                                             QCutoff = 0.05,
                                             randomVar = c("study", "family_ID"),
                                             nnodes = 10)

metad_all %>%
  mutate(status = ifelse(Qs <= 0.1, yes = "OK_sd", no = "NS")) %>%
  metadeconfoundR::BuildHeatmap(., q_cutoff = 0.1, d_cutoff = 0.05) #, keepMeta = colnames(test_metad$Ps))


# inspect significant associations:
full_model_first_solid_food <- lmer ( (MAZ) ~   age * birth_weight + (1 | family_ID) + (1 | study),
                     REML = F, data = pred_metadata_df)
no_slope_model_first_solid_food <- lmer ( (MAZ) ~   age + birth_weight + (1 | family_ID) + (1 | study),
                         REML = F, data = pred_metadata_df)
no_intercept_model_first_solid_food <- lmer ( (MAZ) ~  age + age:birth_weight + (1 | family_ID) + (1 | study),
                             REML = F, data = pred_metadata_df)
lrtest (full_model_first_solid_food, no_slope_model_first_solid_food) 
lrtest (full_model_first_solid_food, no_intercept_model_first_solid_food) 
pred_metadata_df %>%
  filter(!is.na(birth_weight)) %>%
  ggplot(., aes(x = birth_weight, y = MAZ, color = birth_weight)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm") +
  facet_wrap(~study)

# Only young:
pred_metadata_df_young <- pred_metadata_df %>% filter(age < 180)
ftrs_young <- pred_metadata_df_young %>% select(MAZ, MAZ_west, MAZ_nowest, Shannon, Observed) %>% data.frame()
mvars_young <- pred_metadata_df_young %>% select(-c("MAZ", "MAZ_west", "MAZ_nowest", "Shannon", "Observed", "subject_ID",
                                                    "country", "predicted_age",
                                                    "group_sd", "group_median",
                                                    "predicted_age_west", "group_sd_west", "group_median_west",
                                                    "predicted_age_nowest", "group_sd_nowest", "group_median_nowest")) %>%
  select(sex, everything()) %>%
  data.frame()

metad_young <- metadeconfoundR::MetaDeconfound(featureMat = ftrs_young,
                                             metaMat = mvars_young,
                                             fixedVar = c("age"),
                                             returnLong = T,
                                             # collectMods = T,
                                             DCutoff = 0.1,
                                             startStop = "naiveStop",
                                             QCutoff = 0.05,
                                             randomVar = c("study", "family_ID"),
                                             nnodes = 10)

metad_young %>%
  mutate(status = ifelse(Qs <= 0.1, yes = "OK_sd", no = "NS")) %>%
  metadeconfoundR::BuildHeatmap(., q_cutoff = 0.1, d_cutoff = 0.05) #, keepMeta = colnames(test_metad$Ps))

pred_metadata_df_young %>% filter(!is.na(exclusive_breastfeeding)) %>%
  ggplot(., aes( x = exclusive_breastfeeding, y = MAZ, color = study)) +
  geom_point() +
  geom_smooth()
# Only western:
pred_metadata_df_west <- pred_metadata_df %>% filter(lifestyle_westernized == 1)
ftrs_west <- pred_metadata_df_west %>% select(MAZ, MAZ_west, MAZ_nowest, Shannon, Observed) %>% data.frame()
mvars_west <- pred_metadata_df_west %>% select(-c("MAZ", "MAZ_west", "MAZ_nowest", "Shannon", "Observed", "subject_ID",
                                                  "country", "predicted_age",
                                                  "group_sd", "group_median",
                                                  "predicted_age_west", "group_sd_west", "group_median_west",
                                                  "predicted_age_nowest", "group_sd_nowest", "group_median_nowest")) %>%
  select(sex, everything()) %>%
  data.frame()


metad_west <- metadeconfoundR::MetaDeconfound(featureMat = ftrs_west,
                                               metaMat = mvars_west,
                                               fixedVar = c("age"),
                                               returnLong = T,
                                               # collectMods = T,
                                               DCutoff = 0.1,
                                               startStop = "naiveStop",
                                               QCutoff = 0.05,
                                               randomVar = c("study", "family_ID"),
                                               nnodes = 10)

metad_west %>%
  mutate(status = ifelse(Qs <= 0.1, yes = "OK_sd", no = "NS")) %>%
  metadeconfoundR::BuildHeatmap(., q_cutoff = 0.1, d_cutoff = 0.05) #, keepMeta = colnames(test_metad$Ps))



pred_metadata_df %>%
  mutate(birthmode_v = birthmode == "v",
         sex_m = sex == "m") %>%
  fastDummies::dummy_columns(select_columns = "food", remove_selected_columns = T) %>%
  pivot_longer(cols = c("antibiotics_any", "antibiotics_before", "antibiotics_one_week_before",
               "birthmode_v", "sex_m"), names_to = "variable",
               values_to = "value") %>%
  mutate(value = as.factor(value)) %>%
  ggplot(.,aes(x = age, y = predicted_age, color = value)) +
  geom_point(size = 0.1, alpha = 0.3) +
  geom_smooth(method = "loess", se = T) +
  facet_wrap(~variable + lifestyle_westernized, scales = "free_y")


pred_metadata_df %>%
  mutate(birthmode_v = birthmode == "v",
         sex_m = sex == "m") %>%
  fastDummies::dummy_columns(select_columns = "food", remove_selected_columns = T) %>%
  pivot_longer(cols = c("antibiotics_any", "antibiotics_before", "antibiotics_one_week_before",
                        "birthmode_v", "sex_m"), names_to = "variable",
               values_to = "value") %>%
  mutate(value = as.factor(value)) %>%
  ggplot(.,aes(x = age, y = MAZ, color = value)) +
  geom_point(size = 0.1, alpha = 0.3) +
  geom_smooth(method = "loess", se = T) +
  facet_wrap(~variable + lifestyle_westernized, scales = "free_y")


pred_metadata_df_west %>% 
  filter(!is.na(antibiotics_before)) %>%
  ggplot(aes(x = age, y = MAZ, color = antibiotics_before)) +
  geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~study)

pred_metadata_df_west %>% 
  filter(!is.na(antibiotics_before), !study %in% c("roswall_2021", "wampach_2018"), study == "muinck_2018") %>% 
  ggplot(aes(x = age, y = predicted_age, color  = subject_ID %in% c("10_muinck_2018", "11_muinck_2018"), group = subject_ID)) +
  geom_point(size = 0.3) +
  geom_smooth(method = "loess") 

pred_metadata_df_west %>% 
  mutate(antibiotics_before = ifelse(subject_ID == "11_muinck_2018", yes = T, no = antibiotics_before)) %>%
  filter(!is.na(antibiotics_before), !study %in% c("roswall_2021", "wampach_2018", "beller_2021")) %>% 
  ggplot(aes(x = age, y = predicted_age, color = antibiotics_before)) +
  geom_point(size = 0.3) +
  geom_smooth(method = "loess") +
  facet_wrap(~study)



pred_metadata_df_west %>% 
  filter(!is.na(birthmode), !study %in% c("")) %>% 
  ggplot(aes(x = age, y = MAZ, color = birthmode)) +
  geom_point(size = 0.3) +
  geom_smooth(method = "loess") +
  facet_wrap(~study)

pred_metadata_df_west %>% 
  filter(!is.na(birthmode), !study %in% c("")) %>%
  ggplot(aes(x = age, y = Shannon, color = birthmode)) +
  geom_point(size = 0.3) +
  geom_smooth(method = "loess") +
  facet_wrap(~study)


# Malnutrition #######
pred_metadata_df %>% 
  pivot_longer(cols = c("WAZ", "WHZ", "HAZ"), names_to = "score_name", values_to = "zscore_value") %>%
  filter(!is.na(zscore_value), !study %in% c("")) %>%
  ggplot(aes(x = zscore_value, y = MAZ, color = study)) +
  geom_point(size = 0.3) +
  geom_smooth(method = "loess") +
  geom_vline(xintercept = -2) +
  facet_wrap(~score_name + study)



