library(phyloseq)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(microbiome)
library(broom)
library(mgcv)
library(gamm4)
library(gratia)
library(gamlss)
library(gamlss.ggplots)

setwd("/fast/AG_Forslund/rob/mm_index/R_scripts/")




ps_object_genus_raw <- readRDS("/fast/AG_Forslund/rob/mm_index/merged_data/all/all_phyloseq_rf_filter_genus.rds") %>%
  subset_samples(., age <= 730 & age > 1)
mdat <- sample_data(ps_object_genus_raw)
samples_to_keep <- mdat %>%
  group_by(subject_ID, study) %>%
  summarise(samples = list(sample_ID),
            n_samples = length(sample_ID),
            tps_min = min(age),
            tps_max = max(age)) %>%
  filter(n_samples > 2) %>%
  filter(tps_min < 30 & tps_max > 180) %>% # has samples above and below half a year
  pull(samples, name = "subject_ID")

ps_object_genus_comp <- ps_object_genus_raw %>% 
  microbiome::transform(., transform = "compositional") %>%
  aggregate_rare(., detection = 1e-4, prevalence = 0.005, level = "genus") 

ps_object_genus_raref_comp <- ps_object_genus_raw %>% 
  rarefy_even_depth(., sample.size = 2000) %>%
  microbiome::transform(., transform = "compositional") %>%
  aggregate_rare(., detection = 1e-4, prevalence = 0.005, level = "genus") 


# get abundance per individual, timepoint and taxon:
abundance_long <- ps_object_genus_comp %>%
  psmelt() %>%
  select(subject_ID, age, sample_ID, Abundance, OTU, study, lifestyle,
         age_cat)

abundance_long_raref <- ps_object_genus_raref_comp %>%
  psmelt() %>%
  select(subject_ID, age, sample_ID, Abundance, OTU, study, lifestyle,
         age_cat)

  
# prevalence filter, each taxon must be detected 10 times in 2 studies at least
prev_per_study <- abundance_long %>%
  filter(Abundance > 1e-3) %>%
  group_by(OTU, study) %>%
  summarize(prev = n())

# when each taxon was detected first and last time in a subject.
# biased, as some subjects to don't span the entire 2 years
first_last_timepoint <- abundance_long %>%
  filter(subject_ID %in% names(samples_to_keep)) %>%
  # filter(Abundance > 1e-3) %>%
  # filter(OTU %in% (prev_per_study %>% # maybe change to top 50 abundant OTUs?
  #                    group_by(OTU) %>%
  #                    filter(sum(prev > 10) >= 2) %>%
  #                    pull(OTU) %>%
  #                    unique)) %>%
  group_by(subject_ID, OTU) %>%
  summarize(first = min(age),
            last = max(age),
            lifestyle = unique(lifestyle),
            study = unique(study)) %>%
  filter(!OTU %in% c("Unknown", "Other")) %>%
  ungroup()


first_last_timepoint %>%
  # left_join(., mdat %>% select(subje, sepsis_lons_binary, antibiotics_neonatal) %>% distinct()) %>%
  group_by(OTU) %>%
  mutate(OTU = paste(OTU, length(OTU))) %>%
  # filter(length(OTU) >= 10) %>%
  ungroup %>%
  ggplot(., aes(x = first, y = last, color = lifestyle)) +
  # geom_jitter() +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = c(30,180)) +
  geom_hline(yintercept = c(30,180))+
  facet_wrap(~OTU)
ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/first_last_time_colonizers.pdf",
       width = 10, height = 5.5)

# better check overall prevalence
prev_per_month <- abundance_long %>%
  filter(OTU %in% (prev_per_study %>% # maybe change to top 50 abundant OTUs?
                     group_by(OTU) %>% 
                     filter(sum(prev > 10) >= 2) %>% 
                     pull(OTU) %>% 
                     unique)) %>%
  mutate(month = ceiling(age / 30)) %>%
  group_by(OTU, month, lifestyle) %>%
  summarize(prev = mean(Abundance > 1e-3)) %>%
  ungroup()

prev_per_month_raref <- abundance_long_raref %>%
  filter(OTU %in% (prev_per_study %>% # maybe change to top 50 abundant OTUs?
                     group_by(OTU) %>% 
                     filter(sum(prev > 10) >= 2) %>% 
                     pull(OTU) %>% 
                     unique)) %>%
  mutate(month = ceiling(age / 30)) %>%
  group_by(OTU, month, lifestyle) %>%
  summarize(prev = mean(Abundance > 1e-3)) %>%
  ungroup()


# all taxa
ggplot(prev_per_month %>% filter(OTU %in% c("Treponema", "Succinivibrio", "Prevotella")), aes(x = month, y = prev, color = lifestyle)) +
  geom_point() +
  facet_wrap(~OTU)

  # check each taxon per lifestyle, what is the direction over time
# run per lifestyle
slopes <- prev_per_month %>%
  mutate(prev = case_when(prev == 0 ~ prev + 0.0000001,
                          prev == 1 ~ prev - 0.0000001,
                          .default = prev)) %>% 
  group_by(OTU, lifestyle) %>%
  mutate(month_centered = month - month[which.max(prev)]) %>% # shift the max/min to the max of prev, would be at month 0 otherwise
  group_modify(~ {
    model <- gamlss(prev ~ month, data = ., family = BE) # + poly(month_centered, 2)
    res_df <- summary(model) %>% 
      data.frame %>%
      rownames_to_column("factor") %>%
      filter(!factor %in% c("X.Intercept.", "X.Intercept..1"))
    res_df
  }) %>%
  ungroup %>%
  mutate(p_adj = p.adjust(Pr...t..)) %>% # should I do that for all p-values or per term?
  pivot_wider(names_from = "factor", values_from = -c("OTU", "factor", "lifestyle")) %>%
  mutate(pattern = case_when((p_adj_month < 0.05) & Estimate_month < 0 ~ "Early",
                             (p_adj_month < 0.05) & Estimate_month > 0 ~ "Late",
                             (p_adj_month > 0.05) ~ "None"))

# early colonizers:
prev_per_month %>%
  left_join(., slopes) %>%
  filter(p_adj_month < 0.05,
         Estimate_month < 0) %>%
  ggplot(., aes(x = month, y = prev, color = lifestyle)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~OTU, scales = "free_y")
ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/early_colonizers.png")

# late colonizers
prev_per_month %>%
  left_join(., slopes) %>%
  filter(p_adj_month < 0.05,
         Estimate_month > 0) %>%
  ggplot(., aes(x = month, y = prev, color = lifestyle)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~OTU, scales = "free_y")
ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/late_colonizers.png")

# transient colonizers/no pattern
prev_per_month %>%
  left_join(., slopes) %>%
  filter(p_adj_month > 0.05) %>% #, p_adj_poly.month_centered..2.2 < 0.05) %>%
  ggplot(., aes(x = month, y = prev, color = lifestyle)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~OTU, scales = "free_y")
ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/transient_colonizers.png")

save(prev_per_month, prev_per_month_raref, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/prevalences.RData")
save(ps_object_genus_comp, ps_object_genus_raref_comp, file = "/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/phyloseq_supplement.RData")

# compare slopes:
library(lmtest)
full_model_prev <- lm(prev ~ month * lifestyle, data = prev_per_month %>% filter(OTU == "Prevotella"))
slope_model_prev <- lm(prev ~ month + lifestyle, data = prev_per_month %>% filter(OTU == "Prevotella"))
lrtest(slope_model_prev, full_model_prev)

full_model_lac <- lm(prev ~ month * lifestyle, data = prev_per_month %>% filter(OTU == "Lactobacillus"))
slope_model_lac <- lm(prev ~ month + lifestyle, data = prev_per_month %>% filter(OTU == "Lactobacillus"))
lrtest(slope_model_lac, full_model_lac)

# early vs late colonizers volcano:
slopes_to_plot <- slopes %>%
  mutate(pattern = case_when(p_adj_month < 0.05 & Estimate_month > 0 ~ "late",
                             p_adj_month < 0.05 & Estimate_month < 0 ~ "early",
                             .default = "NA")) 
  
ggplot(slopes_to_plot, aes(y = -log(p_adj_month), x = Estimate_month, color = lifestyle, shape = pattern, group = OTU)) +
  geom_line(color = "black", size = 0.1, data = slopes_to_plot %>% filter(p_adj_month < 0.05)) +
  geom_point() +
  geom_hline(yintercept = -log(0.05), linetype = "dashed", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.6) +
  theme_classic()
ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/early_late_linear_volcano.pdf")
ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/early_late_linear_volcano.png")



patterns_df <- slopes %>%
  group_by(OTU) %>%
  summarize(pattern = case_when(min(p_adj_month) > 0.05 ~ "none", # significant nowhere
                                # significant in only one ls
                                max(p_adj_month) > 0.05 & Estimate_month[which.min(p_adj_month)] > 0 ~ paste0("late_", lifestyle[which.min(p_adj_month)]), 
                                max(p_adj_month) > 0.05 & Estimate_month[which.min(p_adj_month)] < 0 ~ paste0("early_", lifestyle[which.min(p_adj_month)]),
                                # significant in both 
                                max(p_adj_month) <= 0.05 & prod(Estimate_month) < 0 ~ "bidirectional", # different dirs
                                max(p_adj_month) <= 0.05 & prod(Estimate_month) > 0 & sum(Estimate_month) > 0 ~ "late_both", # positive slope
                                max(p_adj_month) <= 0.05 & prod(Estimate_month) > 0 & sum(Estimate_month) < 0 ~ "early_both") # negative slope
  )
table(patterns_df$pattern)






# check lifestyle specificity #########################
# run models as below with gamlss, maybe compare results from there with results from here comparing different trajectories
difference_in_ls <- slopes %>%
  filter(`p_adj_month` <= 0.05) %>%
  group_by(OTU) %>%
  summarize(n_ls = length(unique(lifestyle)),
            ls_diff = prod(estimate_month)) %>%
  filter(n_ls == 1 | ls_diff < 0) # only one ls signif, or different directions

# lifestyle specific colonizers:
prev_per_month %>%
  filter(OTU %in% difference_in_ls$OTU) %>%
  ggplot(., aes(x = month, y = prev, color = lifestyle)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~OTU, scales = "free_y")
ggsave("/fast/AG_Forslund/rob/mm_index/merged_data/all/lifestyle_specific_colonizers.png")

# check numbers, how many are early/late/consistently
patterns_df <- slopes %>%
  group_by(OTU) %>%
  summarize(pattern = case_when(min(p_adj_month) > 0.05 ~ "none", # significant nowhere
                                # significant in only one ls
                                max(p_adj_month) > 0.05 & estimate_month[which.min(p_adj_month)] > 0 ~ paste0("late_", lifestyle[which.min(p_adj_month)]), 
                                max(p_adj_month) > 0.05 & estimate_month[which.min(p_adj_month)] < 0 ~ paste0("early_", lifestyle[which.min(p_adj_month)]),
                                max(p_adj_month) <= 0.05 & prod(estimate_month) < 0 ~ "bidirectional", # significant in both, but different dirs
                                max(p_adj_month) <= 0.05 & prod(estimate_month) > 0 & sum(estimate_month) > 0 ~ "late_both", # positive slope
                                max(p_adj_month) <= 0.05 & prod(estimate_month) > 0 & sum(estimate_month) < 0 ~ "early_both")
            )

################################################################################

# check lifestyle effect in model
slopes_lifestyle_effect <- prev_per_month %>%
  group_by(OTU) %>%
  do({
    model <- lm(rank(prev) ~ month * lifestyle, data = .)  # Fit model for each group
    tidy(model) %>%                      # Extract model summary
      filter(term != "(Intercept)") 
  }) %>%
  ungroup %>% 
  mutate(p_adj = p.adjust(p.value)) %>%
  pivot_wider(id_cols = "OTU", names_from = "term", values_from = -c("OTU", "term"))

# does month:lifestyle agree with lifestyle specific differences?
combined_analyses <- slopes %>%
  pivot_wider(id_cols = OTU, names_from = lifestyle, values_from = -c("OTU", "lifestyle")) %>%
  left_join(.,slopes_lifestyle_effect) %>%
  mutate(diff_lifestyle = estimate_month_industrialized - estimate_month_non_industrialized)
ggplot(combined_analyses, aes(x = diff_lifestyle, y = `estimate_month:lifestyleindustrialized`)) +
  geom_point()


# lifestyle effect is significant in those taxa:
prev_per_month %>%
  left_join(., slopes_lifestyle_effect) %>%
  filter(p_adj_month <= 0.05) %>%
  filter(`p_adj_month:lifestyleindustrialized` <= 0.05) %>%
  filter(`estimate_month:lifestyleindustrialized` < 0) %>%
  ggplot(., aes(x = month, y = prev, color = lifestyle)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~OTU, scales = "free_y")

quadratic_effect <- prev_per_month %>%
  group_by(OTU) %>%
  do({
    model <- lm(rank(prev) ~ I(month^2) + lifestyle, data = .)  # Fit model for each group
    tidy(model) %>%                      # Extract model summary
      filter(term != "(Intercept)") 
  }) %>%
  ungroup %>% 
  mutate(p_adj = p.adjust(p.value)) %>%
  pivot_wider(id_cols = "OTU", names_from = "term", values_from = -c("OTU", "term"))

linear_quadratic_effects <- left_join(slopes_lifestyle_effect, quadratic_effect, by = "OTU")
plot(linear_quadratic_effects$estimate_month, linear_quadratic_effects$`estimate_I(month^2)`)

# quadratic effects per lifestyle
# quadratic term is significant for taxa, where linear term is also significant, no improvement
quadratic_effect_per_ls <- prev_per_month %>%
  group_by(OTU, lifestyle) %>%
  do({
    model <- lm(rank(prev) ~ month + I(month^2), data = .)  # Fit model for each group
    tidy(model) %>%                      # Extract model summary
      filter(term != "(Intercept)") 
  }) %>%
  ungroup %>% 
  pivot_wider(id_cols = c("OTU", "lifestyle"), names_from = "term", values_from = -c("OTU", "lifestyle", "term")) %>%
  mutate(`p.adjust_I(month^2)` = p.adjust(`p.value_I(month^2)`))

# # transient colonizers
# prev_per_month %>%
#   left_join(., slopes) %>%
#   filter(`p_adj_I(month^2)` <= 0.05) %>%
#   ggplot(., aes(x = month, y = prev, color = lifestyle)) +
#   geom_point() +
#   geom_smooth() +
#   facet_wrap(~OTU, scales = "free_y")

  
# early vs late vs transient colonizers:
slopes %>%
  left_join(quadratic_effect_per_ls) %>%
  mutate(pattern = case_when(p_adj_month < 0.05 & estimate_month > 0 ~ "late",
                             # p_adj_month > 0.05 & `p.adjust_I(month^2)` < 0.05 ~ "transient",
                             p_adj_month < 0.05 & estimate_month < 0 ~ "early",
                             .default = "NA")) %>%
  ggplot(., aes(y = -log(`p_adj_month`), x = `estimate_month`, color = lifestyle, shape = pattern)) +
  geom_point() +
  # facet_wrap(~lifestyle) +
  theme_classic() +
  geom_hline(yintercept = -log(0.05))


# GAMs #########################################################################

taxa_names <- unique(prev_per_month$OTU)
names(taxa_names) <- taxa_names
test_data <- prev_per_month %>% filter(OTU == "Lactobacillus") %>%
  mutate(lifestyle_industrialized = as.numeric(lifestyle == "industrialized"),
         lifestyle = as.factor(lifestyle))
test_gam <- gam(prev ~ s(month,lifestyle_industrialized, k = 10), data = test_data,
                method="REML")
test_lmer <- lmer(prev ~ month + (1|lifestyle), data = test_data)
draw(test_gam)

coef(test_gam)
fixef(test_lmer)
variance_comp(test_gam)
summary(test_lmer)$varcor


test_data_pred <- with(test_data,
                       expand.grid(month=seq(min(month), max(month), length=100),
                                   lifestyle_industrialized=unique(lifestyle_industrialized)))
test_data_pred <- cbind(test_data_pred,
                        predict(test_gam, 
                                test_data_pred, 
                                se.fit=TRUE, 
                                type="response"))
ggplot(test_data_pred, aes(x = month, y = fit, group = lifestyle_industrialized,
                    colour = lifestyle_industrialized)) +
  geom_line() +
  geom_point(data = test_data, aes(y = prev)) +
  facet_wrap(~ lifestyle_industrialized)




prev_per_month <- prev_per_month %>%
  mutate(lifestyle_industrialized = as.numeric(lifestyle == "industrialized"))

gams_lfs <- prev_per_month %>%
  group_by(OTU) %>%
  do({
    model <- gam(prev ~ s(month, bs = "tp") + s(lifestyle_industrialized, bs = "re"), data = .,
                 method="REML", family="gaussian")  # Fit model for each group
    # tidy(model) %>%                      # Extract model summary
    #   filter(term != "(Intercept)") 
  }) %>%
  ungroup %>% 
  mutate(p_adj = p.adjust(p.value)) %>%
  pivot_wider(id_cols = "OTU", names_from = "term", values_from = -c("OTU", "term"))
draw(gams_lfs$Lactobacillus)

# with gamlss:
"Bifidobacterium"
"Lactobacillus"
test_data <- prev_per_month %>% filter(OTU == "Lactobacillus") %>%
  mutate(lifestyle_industrialized = as.numeric(lifestyle == "industrialized"),
         lifestyle = as.factor(lifestyle)) %>% select(month, prev, lifestyle_industrialized, lifestyle) %>% data.frame
test_data$prev[test_data$prev == 1] <- 0.999999999
fitsum<-gamlss::gamlss(prev ~ month + lifestyle_industrialized, data = test_data)
fitsum_be <-gamlss::gamlss(prev ~ month + lifestyle_industrialized, data = test_data, trace = FALSE, family = BE)
fitsum_pb<-gamlss::gamlss(prev ~ pb(month) + lifestyle_industrialized, data = test_data, trace = FALSE)
fitsum_pb_r<-gamlss::gamlss(prev ~ pb(month) + random(lifestyle), data = test_data, trace = FALSE)
fitsum_pb_be<-gamlss::gamlss(prev ~ pb(month) + lifestyle_industrialized, data = test_data, trace = FALSE, family = BE)
fitsum_pb_be_s<-gamlss::gamlss(prev ~ pb(month) + lifestyle_industrialized, sigma.formula = ~ lifestyle_industrialized + month, data = test_data, trace = FALSE, family = BE)

fitsum_pb_be_by<-gamlss::gamlss(prev ~ pb(month) * lifestyle_industrialized, data = test_data, family = BE, trace = FALSE)
fitsum_pb_be_shape<-gamlss::gamlss(prev ~ pb(month) + lifestyle_industrialized, data = test_data, family = BE, trace = FALSE)
fitsum_pb_be_int<-gamlss::gamlss(prev ~ pb(month) + pb(month):lifestyle_industrialized, data = test_data, family = BE, trace = FALSE)
LR.test(fitsum_pb_be_shape, fitsum_pb_be_by) # shape
LR.test(fitsum_pb_be_int, fitsum_pb_be_by) # intercept
# Extract fitted values
test_data$fitted_interaction <- fitted(fitsum_pb_be_by)
test_data$fitted_no_interaction <- fitted(fitsum_pb_be_red)

# Plot the curves for the interaction model
ggplot(test_data, aes(x = month, y = prev, color = as.factor(lifestyle_industrialized))) +
  geom_point(alpha = 0.5) +  # Raw data
  geom_line(aes(y = fitted_interaction), size = 1) +  # Fitted curves
  labs(color = "Lifestyle Industrialized") +
  ggtitle("Fitted Curves for Interaction Model")
library(gratia)
library(mgcv)
gam_fit <- gam(prev ~ s(month, by = lifestyle_industrialized, k=20, bs="tp", m=8) + lifestyle_industrialized, data = test_data, family = betar())
draw(gam_fit)
smooth_estimates(gam_fit)
derivatives <- derivatives(gam_fit, term = "s(month)", partial_match = TRUE)
fitted_values(gam_fit, data = test_data) |>
  mutate(fz = factor(lifestyle_industrialized)) |>
  ggplot(aes(x = month, y = .fitted, colour = fz, group = fz)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = fz, colour = NULL),
              alpha = 0.2
  ) +
  geom_line() +
  labs(
    title = "Fitted values from model",
    y = expression(hat(y)), colour = "z", fill = "z"
  )


plot(fitsum)
plot(fitsum_pb_be)
plot(fitsum_pb_be_by)
resid_index(fitsum)
resid_index(fitsum_pb_be)
resid_index(fitsum_pb_be_by)
resid_wp(fitsum)
resid_wp(fitsum_pb_be)
resid_wp(fitsum_pb_be_by)
moment_bucket(fitsum_pb, fitsum, fitsum_be, fitsum_pb_be, fitsum_pb_be_by)

fitted_terms(fitsum, data = test_data, rug = T, ylim = "free")
fitted_terms(fitsum_pb)
fitted_terms(fitsum_pb_be_by)
fitted_terms(fitsum_pb_be, data = test_data, rug = T, ylim = "free")
pe_param_grid(fitsum_pb_be_by, c("month", "lifestyle_industrialized"), filled = T)
pe_param(fitsum_pb_be_by, "month")
pe_param(fitsum_pb_be_by, "lifestyle_industrialized")
resp_mu(fitsum)
resp_mu(fitsum_pb)
resp_mu(fitsum_pb_be_by)
resp_mu(fitsum_pb_be)



GAIC(fitsum_pb, fitsum, fitsum_be, fitsum_pb_be, fitsum_pb_be_by,
     fitsum_pb_r, fitsum_pb_be_s)
model_GAIC_lollipop(fitsum_pb, fitsum, fitsum_be, fitsum_pb_be, fitsum_pb_be_by,
     fitsum_pb_r, fitsum_pb_be_s)


test_dists <- prev_per_month %>%
  group_by(OTU) %>%
  do({
    res <- fitDist(.$prev, type = c("real0to1"))  # test for each group
    res$fits %>% data.frame %>% rownames_to_column("rowname")
    }) %>%
  ungroup %>% 
  mutate(p_adj = p.adjust(p.value)) %>%
  pivot_wider(id_cols = "OTU", names_from = "term", values_from = -c("OTU", "term"))
test_dists %>%
  group_by(rowname) %>%
  summarize(mn = median(`.`)) %>% View()
test_dists %>%
  filter(!rowname %in% c("GP", "GU", "LOGNO", "LOGNO2", "PARETO2", "PARETO2o",
                         "SN1", "BCT", "WEI", "WEI2", "WEI3", "exGAUS",
                         "EGB2", "TF", "TF2", "SHASHo", "SHASHo2", "NO")) %>%
  ggplot(.,aes(x = `.`,))+# color = rowname)) +
  geom_density() +
  geom_vline(xintercept = c(0, -50)) +
  facet_wrap(~rowname) +
  xlim(-300, 100)

test <- fitDist(test_data$prev)



fitcoef<-as.data.frame(fitsum$coef.table[rownames(fitsum$coef.table)!="(Intercept)",])
fitcoef[,"varname"]<-rownames(fitcoef)
fitcoef[,"id"]<-"test"
fitcoefw<-reshape(fitcoef, idvar="id", timevar="varname", direction="wide")

np<-length(colnames(fitcoefw)[grep("Pr(>|t|)",colnames(fitcoefw))])
fitcoefw[,sub('.*\\.', 'pval.adjust.',colnames(fitcoefw)[grep("Pr(>|t|)",colnames(fitcoefw))])]<-
  apply(fitcoefw[,colnames(fitcoefw)[grep("Pr(>|t|)",colnames(fitcoefw))]],2,p.adjust,method = p.adjust.method)

test_data



