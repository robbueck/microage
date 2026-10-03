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



ps_object_genus_raw <- readRDS("../data/all_phyloseq_rf_filter_genus.rds") %>%
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


save(prev_per_month, prev_per_month_raref, file = "../data/prevalences.RData")
save(ps_object_genus_comp, ps_object_genus_raref_comp, file = "../data/phyloseq_supplement.RData")