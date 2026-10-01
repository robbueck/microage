library(ggplot2)
library(tidyverse)
library(viridis)
library(gridExtra)
library(grid)
library(ggpubr)
library(ggpattern)
library(ggpmisc)
library(patchwork)
library(cowplot)
library(rstatix)
library(ggprism)
library(ggh4x)
library(grid)
library(ggtext)
library(glue)
library(ggExtra)
library(wesanderson)
library(ggdist)


source("/fast/AG_Forslund/rob/mm_index/R_scripts/functions.R")
setwd("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures")
# define global colors for lifestyle:
fixed_colors <- c(Industrialized = "#737125", `non-Industrialized` = "#1A97C8", Combined = "#BC85A9")

# Figure 1: Methods ############################################################
## studies map #####################
load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/16S_study_locations_ls.RData")
hdi_data <- read.table("/fast/AG_Forslund/rob/mm_index/R_scripts/hdi_data.txt",
                       sep = "\t", col.names = c("region", "hdi")) 
world_coordinates <- left_join(world_coordinates, hdi_data, by = "region") %>%
  mutate(HDI = hdi, .keep = "unused")


(map_16S <- ggplot() +
    geom_map(
      data = world_coordinates, map = world_coordinates,
      aes(long, lat, map_id = region, alpha = HDI),
      color = "grey70", fill = "grey", size = 0.15)+
    geom_point(data = data_loc_comb %>% filter (Type == "16S") %>%
                 mutate(Lifestyle = ifelse(lifestyle == "industrialized", yes = "Industrialized", no = "non-Industrialized")),
               aes(x = longitude, y = latitude, size = sample_count,  fill = Lifestyle),
               alpha = 0.6,
               position = position_jitter(width=0.7),
               shape=21) +
    # coord_map("gall") +
    ylim(-50, 78) +
    xlim(-155, 175)+
    labs(size = "Samples") +
    scale_alpha_continuous(range = c(0, 0.7)) +
    scale_fill_manual(values = fixed_colors,
                      guide = guide_legend(override.aes = list(alpha = 1, size = 5, order = 2))) +
    scale_size(breaks = c(500, 1000, 2000), guide = guide_legend(order = 1)) +
    ggtitle("A)") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          # panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 18),
          plot.margin = margin(0, 0, 0, 0, "cm"),
          legend.position = c(0.14,0.41),
          # panel.background = element_rect(fill = "#BDE3FF"),    
          panel.background = element_blank(),    
          legend.background = element_blank(),
          legend.title = element_text(size=15),
          legend.text = element_text(size=15)))
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/16S_study_locations_ls.png",
       width = 19.2, height = 12.8, units = "cm", plot = map_16S)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/16S_study_locations_ls.svg",
       width = 19.2, height = 12.8, units = "cm", plot = map_16S)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/16S_study_locations_ls.pdf",
       width = 19.2, height = 12.8, units = "cm", plot = map_16S)

# Figure 2 #####################################################################
## A: Alpha div ##############
load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/a_div_genus_present.RData")
meta_df %>% filter(study == "raman_2019") %>% select(lifestyle, country) %>% table
meta_df_noraman <- meta_df %>%
  mutate(Lifestyle = ifelse(lifestyle == "industrialized",
                            yes = "Industrialized",
                            no = "non-Industrialized"),
         group = "No Raman") %>%
  filter(# !study %in% c("raman_2019", "gehrig_2019", "subramanian_2014"), 
         !country %in% c("BANGLADESH", "SOUTH_AFRICA"),
         Lifestyle == "non-Industrialized")
linetypes_1 <- c("No Raman" = "88")
(a_div_s_genus <- ggplot(meta_df %>% mutate(Lifestyle = ifelse(lifestyle == "industrialized",
                                                               yes = "Industrialized",
                                                               no = "non-Industrialized"),
                                            diversity_shannon = diversity_shannon),
                         aes(y = diversity_shannon, x = age, color = Lifestyle,
                             group = Lifestyle)) +
    geom_point(alpha = 0.1, size = 0.3) +
    geom_smooth(method = "loess", size = 1, se = T, linetype = "88", 
                aes(fill = Lifestyle),
                data = meta_df_noraman,
                color = "#1A97C8", show.legend = F) +
    geom_smooth(aes(fill = Lifestyle), method = "loess", size = 1, se = T, alpha = 0.3) +
    scale_color_manual(values = fixed_colors,
                       guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
    scale_fill_manual(values = fixed_colors) +
    labs(x = "Chronological age [days]",
         y = "Shannon diversity (genus)",
         title = NULL) +
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          # legend.position = c(0.8, 0.1),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)))
(a_div_s_genus_hist <- ggMarginal(a_div_s_genus + theme(legend.position = 'hidden'),
                                   type = "density",      
                                   margins = "x",             
                                   groupFill = TRUE,
                                  size = 5))
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/a_div_genus_present.png",
       width = 11, height = 6, plot = a_div_s_genus_hist)

## B: Alpha div difference #############################
(a_div_diff <- bind_rows(alpha_div_sliding_window, alpha_div_sliding_window_no_raman) %>%
   mutate(Samples = n_total,
          p_val = ifelse(p.adj < 0.05, yes = "< 0.05", no = "> 0.05")) %>%
   ggplot(aes(x = sliding_window, y = delta_shannon, color = set, alpha = p_val, 
              group = set,
              size = Samples)) +
   geom_point() +
   geom_line(size = 0.4, alpha = 0.5, show.legend = F) +
   geom_vline(xintercept = 226) +
   scale_color_manual(
     name = NULL,   # This also sets the legend title
     values = c(
       "Complete dataset" = "#419A97",
       "No Raman et al. 2019" = "#C04136"),
     guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
   scale_alpha_manual(name = "p-Value",
                      values = c("< 0.05" = 1,
                                 "> 0.05" = 0.2)) +
   guides(color = guide_legend(order = 1, override.aes = list(shape = 15, size = 5, alpha = 1)),
          alpha = guide_legend(oder = 2, override.aes = list(size = 5)),
          size = guide_legend(order = 3)) +
   labs(title = NULL,
        y = "\U0394 Shannon diversity [I - NI]",
        x = "Chronological age [days]") +
   theme(axis.title.x = element_text(size = 16),
         axis.text.x = element_text(size = 14),
         axis.title.y = element_text(size = 16),
         axis.text.y = element_text(size = 14),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         panel.grid.major.y = element_line(color = "grey",
                                           size = 0.3,
                                           linetype = 2),
         strip.text = element_text(size = 14),
         strip.background = element_blank(),
         legend.position = "top",
         # legend.position = c(0.65, 0.85),
         legend.text = element_text(size = 12),
         legend.title = element_text(size = 12),
         legend.direction = "vertical", legend.box = "horizontal"))


## C: Beta div ##############
load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/all_pcoa_genus_combined.RData")

explained_var <- ((pcoa_aitch_genus_1$sdev^2/sum(pcoa_aitch_genus_1$sdev^2)) * 100) %>%
  round(.,digits = 1)
plot_pcoa_data <- pcoa_aitch_genus_1$x %>% 
  as.data.frame() %>%
  select(PC1, PC2, PC3, PC4, PC5) %>%
  rownames_to_column("run_accession") %>%
  left_join(pcoa_metadata, by = "run_accession") %>%
  mutate(Lifestyle = ifelse(lifestyle == "industrialized", yes = "Industrialized",
                            no = "non-Industrialized"))

# associations between metadata and axes:
plot_pcoa_data %>% select(PC1, PC2, PC3, PC4, sample_sum, age, lifestyle, shannon_asv) %>%
  mutate(lifestyle_industrialized = ifelse(lifestyle == "industrialized", yes = 1, no = 0),
         .keep = "unused") %>%
  cor(., method = "spearman", use = "pairwise.complete.obs")


(pcoa_age <- ggplot(plot_pcoa_data, aes(x = PC1, y = PC2, color = age)) +
    geom_point(alpha = 0.8, size = 0.3) +
    scale_color_gradientn(guide = guide_colorbar(direction = "horizontal",
                                                 display = "gradient",
                                                 title.position = "top",
                                                 barwidth = 10),
                          colors = wes_palette("Zissou1", type = "continuous")) +
    ggnewscale::new_scale_color() +
    stat_ellipse(data = plot_pcoa_data %>% filter(age <= 200),
                 alpha = 1, 
                 aes(color = Lifestyle),
                 show.legend = F) +
    scale_color_manual(values = fixed_colors) +
    coord_equal() +
    annotate(geom = "text", x = min(plot_pcoa_data$PC1) * -0.005,  # Adjust position
             y = max(plot_pcoa_data$PC2) * 0.9,
             label = "R\u00B2 = 0.05; p < 0.001",
             size = 4, hjust = 1) +
    labs(title = NULL,
         x = paste0("PC 1 (", explained_var[1], "%)"),
         y = paste0("PC 2 (", explained_var[2], "%)"),
         color = "Age [days]") +
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)))

(pcoa_lifestyle <- ggplot(plot_pcoa_data, aes(x = PC1, y = PC2, color = Lifestyle)) +
    geom_point(alpha = 0.8, size = 0.3) +
    stat_ellipse(alpha = 1) +
    scale_size_manual(values = c(0.3, 5)) +
    scale_color_manual(values = fixed_colors) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    coord_equal() +
    annotate("text", x = min(plot_pcoa_data$PC1) * - 0.005,  # Adjust position
             y = max(plot_pcoa_data$PC2) * 0.9,
             label = "R\u00B2 = 0.04; p < 0.001", size = 4, hjust = 1) +
    labs(title = NULL,
         x = paste0("PC 1 (", explained_var[1], "%)"),
         y = paste0("PC 2 (", explained_var[2], "%)")) +
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.position = "none"))


pcoa_age_ls_combined <- grid.arrange(pcoa_lifestyle, pcoa_age,
                      layout_matrix = rbind(c(1, 1, 2, 2, 2)))

ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/all_pcoa_genus_combined.png",
       width = 14, height = 7, plot = pcoa_age_ls_combined)
## combined Figure 1:
legend_age <- get_legend(pcoa_age)
legend_ls <- get_legend(a_div_s_genus)
legend_diff <- get_legend(a_div_diff)
legends <- cowplot::plot_grid(legend_ls, legend_age, nrow = 1)

(fig_2_combined <- cowplot::plot_grid(cowplot::plot_grid(a_div_s_genus_hist,
                                                         # legends, 
                                                         a_div_diff, # + theme(legend.position = "none"),
                                                         rel_widths = c(1, 1),
                                                         nrow = 1,
                                                         align = "h",
                                                         axis = "b",
                                                         labels = c("A)", "B)"),
                                                         label_size = 18,
                                                         label_fontface = "plain"),
                                      NULL,
                                      cowplot::plot_grid(pcoa_lifestyle,
                                                         pcoa_age + theme(legend.position='hidden'),
                                                         rel_widths = c(1, 1),
                                                         nrow = 1,
                                                         labels = c("C)", "D)"),
                                                         label_size = 18,
                                                         label_x = c(0, -0.005),
                                                         label_fontface = "plain"),
                                      legends,
                                     rel_heights = c(0.9, 0.1, 1, 0.15), ncol = 1))

cowplot::save_plot("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/figure_2_combined.png", 
                   fig_2_combined,dpi = 900,
                   bg = "white", base_height = 14, base_width = 14)
cowplot::save_plot("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/figure_2_combined.pdf", 
                   fig_2_combined,dpi = 900,
                   bg = "white", base_height = 14, base_width = 14)

# Figure 3 #####################################################################
## A: Predicted vs true Age ######################
load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/genus_nested_cv_no_ls.RData")

genus_no_ls_nested_cv_preds_noraman <- genus_no_ls_nested_cv_preds %>%
  mutate(Lifestyle = ifelse(lifestyle == "industrialized",
                            yes = "Industrialized",
                            no = "non-Industrialized"),
         group = "No Raman") %>%
  filter(study != "raman_2019", Lifestyle == "non-Industrialized")
linetypes_1 <- c("No Raman" = "88")

# fit logistic growth curve:
logistic_model <- function(x, A, k, x0) {
  A / (1 + exp(-k * (x - x0)))
}
fit_logistic <- function(data) {
  tryCatch({
    nls(
      rf1 ~ logistic_model(age, A, k, x0),
      data = data,
      start = list(A = quantile(data$rf1, 0.95), k = 1, x0 = median(data$age))
    )
  }, error = function(e) e)
}
add_predictions <- function(data, model) {
    # Predict for all ages in this experiment
    pred <- predict(model, newdata = data.frame(age = data$age))
    data.frame(age = data$age, fitted = pred)
}

fit <- genus_no_ls_nested_cv_preds %>%
  group_by(lifestyle) %>%
  nest() %>%
  mutate(model = map(data, fit_logistic), 
         params = map(model, ~ if (!is.null(.x)) coef(.x) else NA))
plot_data <- fit %>%
  mutate(preds = map2(data, model, add_predictions)) %>%
  select(lifestyle, preds) %>%
  unnest(preds) %>%
  distinct() %>%
  mutate(Lifestyle = ifelse(lifestyle == "industrialized",
                            yes = "Industrialized",
                            no = "non-Industrialized"))

# 95% of carying capacity
# industrialized
plot_data %>% filter(Lifestyle == "Industrialized",
                     fitted >= summary(fit$model[[1]])$coefficients[1,1] * 0.90) %>%
  pull(age) %>% min

plot_data %>% filter(Lifestyle == "non-Industrialized",
                     fitted >= summary(fit$model[[2]])$coefficients[1,1] * 0.90) %>%
  pull(age) %>% min

plot_data_no_raman <- genus_no_ls_nested_cv_preds %>%
  filter(study != "raman_2019", lifestyle == "non_industrialized") %>%
  group_by(lifestyle) %>%
  nest() %>%
  mutate(model = map(data, fit_logistic),
         params = map(model, ~ if (!is.null(.x)) coef(.x) else NA)) %>%
  mutate(preds = map2(data, model, add_predictions)) %>%
  select(lifestyle, preds) %>%
  unnest(preds) %>%
  distinct() %>%
  mutate(Lifestyle = ifelse(lifestyle == "industrialized",
                            yes = "Industrialized",
                            no = "non-Industrialized"))


(age_prediction_dotplot <- genus_no_ls_nested_cv_preds %>% 
    left_join(plot_data, by = c("lifestyle", "age")) %>%
    mutate(Lifestyle = ifelse(lifestyle == "industrialized",
                              yes = "Industrialized",
                              no = "non-Industrialized")) %>%
    select(age, rf1, Lifestyle, study) %>%
    ggplot(.,aes(x = age, y = rf1, color = Lifestyle, fill = Lifestyle)) +
    geom_point(alpha = 0.2, size = 0.3) +
    stat_poly_eq(use_label("R2", "P"),
                 formula = as.formula("rank(y) ~ x"),
                 vstep = 0.1,
                 size = 4) +
    geom_line(data = plot_data, aes(y = fitted, x = age), size = 1, show.legend = F) +
    geom_line(data = plot_data_no_raman, aes(y = fitted, x = age), size = 1, linetype = "88", show.legend = F) +
    scale_color_manual(values = fixed_colors, 
                       guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
    scale_fill_manual(values = fixed_colors) +
    coord_equal() +
    labs(title = NULL,
         x = "Chronological age [days]",
         y = "Microbial age [days]") +
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = c(0.7, 0.15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)))
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/genus_nested_cv_no_ls_point_present.png", 
       width = 10, plot = age_prediction_dotplot)

(age_prediction_dotplot_hist <- ggMarginal(age_prediction_dotplot,
                                  type = "density",      
                                  margins = "x",             
                                  groupFill = TRUE))


## B: Lifestyle performance #######################
load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/lifestyle_perfomance_genus.RData")


lifestyle_lm_list_genus <- lifestyle_lm_list_genus %>%
  mutate(Lifestyle = case_when(lifestyle == "industrialized" ~ "Industrialized",
                               lifestyle == "non_industrialized" ~ "non-Industrialized",
                               lifestyle == "combined" ~ "Combined"),
         Lifestyle = factor(Lifestyle, levels = c("non-Industrialized", "Industrialized")),
         Training_set = case_when(training_set == "industrialized" ~ "Industrialized",
                                  training_set == "non_industrialized" ~ "non-Industrialized",
                                  training_set == "combined" ~ "Combined"),
         Training_set = factor(Training_set, levels = c("non-Industrialized", "Combined", "Industrialized")))
df_p_val_lifestyle_genus <- lifestyle_lm_list_genus %>%
  arrange(study) %>%
  mutate(Training_set = factor(Training_set)) %>%
  rstatix::group_by(Lifestyle) %>%
  rstatix::wilcox_test(R2 ~ Training_set, paired = T) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>% 
  rstatix::add_xy_position(x = "Lifestyle", dodge = 0.75,
                           step.increase = 0.2) %>%
  mutate(y.position = ifelse(Lifestyle == "Industrialized",
                             yes = y.position - 0.002,
                             no = y.position))

lifestyle_lm_list_genus %>% group_by(Lifestyle, Training_set) %>%
  summarize(mean = mean(R2),
            median = median(R2))

library(ggbeeswarm)
(ls_performance <- ggplot(lifestyle_lm_list_genus, aes(x=Lifestyle, y = R2)) +
    geom_boxplot(aes(fill = Training_set), outlier.shape = NA) +
    geom_quasirandom(aes(color = Training_set), dodge.width = 0.8,
                     shape = 21, width = 0.1, size = 1.5, show.legend = T) +
    geom_quasirandom(aes(fill = Training_set), dodge.width = 0.8,
                     shape = 21, width = 0.1, size = 1.5) +
    scale_color_manual(values = fixed_colors,
                       guide = guide_legend(#title.hjust = 0.5,
                         # title.position = "top",
                         # direction = "horizontal",
                         # position = "bottom",
                         override.aes = list(shape = 15, size = 5, alpha = 1))) +
    scale_fill_manual(values = fixed_colors,
                      guide = "none") +
    labs(title = NULL,
         x = "Validation set",
         y = bquote("Performance ["~R^2~"]"),
         color = "Training set",
         fill = "Training set") +
    add_pvalue(df_p_val_lifestyle_genus,
               label = "{p.adj.signif}",
               label.size = 4.5,
               tip.length = 0.01,
               xmin = "xmin",
               xmax = "xmax",
               show.legend = FALSE) +
    ylim(0, 1) +
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = c(0.7, 0.1),
          # legend.position="bottom",
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)))

ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/lifestyle_perfomance_genus.png", 
       width = 10, plot = ls_performance)

(fig_3_combined <- cowplot::plot_grid(age_prediction_dotplot_hist,
                                     ls_performance,
                                     nrow = 1,
                                     labels = c("A)", "B)"),
                                     label_size = 18,
                                     label_fontface = "plain"))

cowplot::save_plot("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/figure_3_combined.png", 
                   fig_3_combined, dpi = 900,
                   bg = "white", base_height = 7.6, base_width = 14)

cowplot::save_plot("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/figure_3_combined.pdf", 
                   fig_3_combined, dpi = 900,
                   bg = "white", base_height = 7.6, base_width = 14)


# Figure 4 ####################################################################

load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/imp_shap_prev.RData")
# select interesting Taxa:
interesting_taxa <- combined_importance_prev_lifestyles_genus %>%
  filter(!is.na(importance)) %>%
  group_by(display_taxon, model) %>%
  summarize(importance = mean(importance)) %>%
  group_by(model) %>% #filter(model == "industrialized") %>%View()
  top_n(19, importance) %>%
  pull(display_taxon) %>%
  as.character() %>%
  unique()
# add some more interesting taxa
interesting_taxa <- c(interesting_taxa, "Clostridium sensu stricto 1",
                      "Lactobacillus", "Collinsella",
                      "Observed", "Shannon") %>%
  unique
## A: Importance, Prevalence ####################
imp_prev_ls <- combined_importance_prev_lifestyles_genus %>%
  mutate(model = ifelse(model == "industrialized",
                        yes = "Industrialized",
                        no = "non-Industrialized"),
         Lifestyle = ifelse(lifestyle == "industrialized",
                            yes = "Industrialized",
                            no = "non-Industrialized"),
         filter_taxon = display_taxon,
         display_taxon = if_else(grepl("_| ", display_taxon) |
                                    display_taxon == "Shannon" |
                                    display_taxon == "Observed",
                                  true = glue("{display_taxon}"),
                                  false = glue("<i>{display_taxon}</i>")),
         display_taxon = factor(display_taxon),
         display_taxon = reorder(display_taxon, importance, FUN = max, na.rm = T),
         display_taxon = factor(display_taxon, levels = unique(c("Shannon", "Observed",
                                                                 levels(display_taxon)))))

imp_prev_ls_filter <- imp_prev_ls %>% filter(filter_taxon %in% interesting_taxa) %>%
  mutate(display_taxon = droplevels(display_taxon),
         filter_taxon = droplevels(filter_taxon))

# numer of important taxa:
taxa_vs_ls <- imp_prev_ls %>% select(Lifestyle, taxon) %>% distinct()
# total taxa
length(unique(taxa_vs_ls$taxon))
# number of important features in each lifestyle:
table(taxa_vs_ls$Lifestyle)
# number of lifestyle specific taxa
taxa_vs_ls$taxon %>% table %>% table
# per lifestyle (for interpretation, switch lifestyle labels):
table(taxa_vs_ls) %>% `==` (0) %>% rowSums()


(imp_prev_combined_plot <- ggplot() +
    geom_col(data = imp_prev_ls_filter %>%
               distinct(display_taxon, prevalence, model),
             aes(y = display_taxon, x = prevalence), alpha = 0.5, fill = "grey", color = "white") +
    geom_boxplot(data = imp_prev_ls_filter,
                 aes(color = Lifestyle, y = display_taxon, x = importance)) +
    scale_color_manual(values = fixed_colors,
                      guide = guide_legend(override.aes = list(
                        shape = 15,
                        size = 5,
                        alpha = 1))) +
    labs(title = NULL,
        x = "Feature importance",
        y = NULL) +
    facet_grid(~ model,
                scales = "free_x") +
    theme(axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 16),
      axis.text.y = element_markdown(size = 12),
      # legend.position = c(0.2, 0.2),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.grid.major.y = element_line(color = "grey",
        size = 0.3,
        linetype = 2),
      strip.text = element_text(size = 14),
      strip.background = element_blank()))

lifestyle_legend <- ggplot() +
  geom_bar(data = imp_prev_ls_filter,
           aes(x = importance, y = display_taxon, fill = Lifestyle),
           position = "dodge",
           stat = "identity") +
  scale_fill_manual(values = fixed_colors,
                    name = "Model",
                    guide = guide_legend(override.aes = list(shape = 15,
                                                             size = 5,
                                                             nrow = 1,
                                                             alpha = 1),
                                         title.position = "top", 
                                         direction = "horizontal")) +
  theme(legend.title = element_text(size = 15, margin = margin(r = 15)),
        legend.text = element_text(size = 15, margin = margin(r = 15)))
prevalence_legend <- ggplot(data.frame(a = 1, b = 2, c = "Prevalence"), aes(x = a, y = b, fill = c)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_manual(
    name = "",         # No title, or change as needed
    values = c("Prevalence" = "grey"),
    guide = guide_legend(override.aes = list(shape = 15,
                                             size = 5,
                                             nrow = 1,
                                             alpha = 1))) +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 15))

## B: Shap ####################

shap_cors_all <- shap_cors %>%
  mutate(model = ifelse(model == "industrialized",
                        yes = "Industrialized",
                        no = "non-Industrialized"),
         display_taxon_filter = display_taxon,
         display_taxon = if_else(grepl("_| ", display_taxon) |
                                   display_taxon == "Shannon" |
                                   display_taxon == "Observed",
                                 true = glue("{display_taxon}"),
                                 false = glue("<i>{display_taxon}</i>")),
         display_taxon = factor(display_taxon)) %>%
  left_join(.,imp_prev_ls %>%
              group_by(display_taxon, model) %>%
              summarize(importance = mean(importance),
                        prevalence = mean(prevalence),
                        .groups = "drop"),
            by = c("display_taxon", "model")) %>%
  mutate(R2 = case_when(is.na(R2) ~ 0,
                        is.na(importance) ~ 0,
                        p_adj > 0.05 ~ 0,
                        .default = R2))

shap_cors_filt <- shap_cors_all %>%
  filter(display_taxon_filter %in% interesting_taxa) %>%
  mutate(display_taxon = factor(display_taxon, levels = levels(imp_prev_ls$display_taxon)))


shap_vals_combined <- bind_rows(shap_model_noindustrialized %>% mutate(model = "non-Industrialized"), 
                                shap_model_industrialized %>% mutate(model = "Industrialized")) %>%
  left_join(short_taxa_names, by = "taxon") %>%
  mutate(display_taxon_filter = display_taxon,
         display_taxon = if_else(grepl("_| ", display_taxon) |
                                   display_taxon == "Shannon" |
                                   display_taxon == "Observed",
                                 true = glue("{display_taxon}"),
                                 false = glue("<i>{display_taxon}</i>")),
         display_taxon = factor(display_taxon)) %>%
  left_join(.,imp_prev_ls %>%
              group_by(display_taxon, model) %>%
              summarize(importance = mean(importance, na.rm = T),
                        prevalence = mean(prevalence, na.rm = T),
                        .groups = "drop"),
            by = c("display_taxon", "model")) %>%
  filter(display_taxon_filter %in% interesting_taxa) %>%
  mutate(display_taxon = factor(display_taxon, levels = levels(imp_prev_ls$display_taxon)),
         shap_value = ifelse(!is.na(importance), yes = shap_value, no = NA))

violin_plot_df <- shap_vals_combined %>%
  mutate(shap_range = cut(shap_value, breaks = 5000)) %>%
  group_by(shap_range, display_taxon, model) %>%
  summarize(shap_widths = n(),
            shap = mean(shap_value, na.rm = T),
            abundance = mean(ab_value, na.rm = T)) %>%
  group_by(display_taxon, model) %>%
  mutate(abundance = log(abundance / max(abundance) + 0.00001)) %>%
  mutate(shap_widths = shap_widths/max(shap_widths)) %>%
  ungroup() %>%
  mutate(display_taxon = factor(display_taxon, levels = levels(imp_prev_ls_filter$display_taxon)))


(shap_violin_plot <- ggplot() +
    geom_segment(data = violin_plot_df, 
                 aes(y = as.numeric(display_taxon) - shap_widths / 2.2,
                     yend = as.numeric(display_taxon) + shap_widths / 2.2,
                     x = shap,
                     xend = shap,
                     color = abundance),
                 size = 1) +
    scale_y_continuous(breaks = seq_along(levels(violin_plot_df$display_taxon)),
                       labels = levels(violin_plot_df$display_taxon),
                       expand = c(0.005, 0.005)) +
    scale_color_gradientn(name = "Feature value",
                          breaks = range(violin_plot_df$abundance, na.rm = TRUE),
                          labels = c("Low", "High"),
                          guide = guide_colorbar(direction = "horizontal",
                                                 display = "gradient",
                                                 title.position = "top",
                                                 barwidth = 10),
                          colors = wes_palette("Zissou1", type = "continuous")) +
    geom_vline(xintercept = 0) +
    xlab("SHAP") +
    facet_grid(~ model, scales = "free_x") +
    theme(axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size=18),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_line(color = "grey",
                                            size = 0.3,
                                            linetype = 2),
          strip.text = element_text(size = 14),
          strip.background = element_blank()))

cowplot::set_null_device("agg")
fig_4_combined <- cowplot::plot_grid(cowplot::plot_grid(imp_prev_combined_plot, 
                                                        shap_violin_plot + theme(legend.position = "none"),
                                                        ncol = 2,
                                                        labels = c("A)", "B)"),
                                                        label_size = 18,
                                                        label_fontface = "plain",
                                                        # align = "v",
                                                        # axis = "b",
                                                        rel_widths = c(1.6, 1),
                                                        label_x = c(0, -0.03)),
                                     cowplot::plot_grid(get_legend(lifestyle_legend), 
                                                        get_legend(prevalence_legend),
                                                        get_legend(shap_violin_plot), 
                                                        nrow = 1, rel_widths = c(1, 0.5, 1)),
                                     nrow = 2,
                                     rel_heights = c(1, 0.13))

cowplot::save_plot("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/figure_4_combined.png", 
                   fig_4_combined, dpi = 900,
                   bg = "white", base_height = 9, base_width = 14)

cowplot::save_plot("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/figure_4_combined.pdf", 
                   fig_4_combined, dpi = 900,
                   bg = "white", base_height = 9, base_width = 14)


# Fig 5: #######################################################################

## SHAP Cors Lifestyles ########################################################
selected_taxa <- c("Staphylococcus", "Faecalibacterium", "Prevotella",
                   "Lactobacillus", "Bifidobacterium")

shap_cors_filt_xy_plot <- shap_cors_all %>%
  select(display_taxon, display_taxon_filter, model, R2) %>%
  mutate(R2 = ifelse(R2 == 0, yes = NA, no = R2)) %>%
  pivot_wider(names_from = "model", values_from = R2) %>%
  filter(!(Industrialized == 0 & `non-Industrialized` == 0), 
         !display_taxon_filter %in% c("Shannon", "Observed")) %>%
  mutate(Classification = case_when(`non-Industrialized` > 0 & Industrialized > 0 ~ "Late colonizer",
                           `non-Industrialized` < 0 & Industrialized < 0 ~ "Early colonizer",
                           is.na(Industrialized) | is.na(`non-Industrialized`) ~ "Single lifestyle",
                           .default = "Mixed"),
         display_taxon_filter = droplevels(display_taxon_filter)) %>%
  mutate(annotation_taxon = case_when(display_taxon_filter %in% selected_taxa ~ display_taxon_filter, 
                                      .default = ""))
table(shap_cors_filt_xy_plot$Classification) %>% prop.table



shap_cors_filt_xy_plot %>%
  # filter(Classification != "Single lifestyle") %>%
  pull(Classification) %>% table
shap_cors_filt_xy_plot %>%
  filter(Classification == "Single lifestyle") %>%
  select(`non-Industrialized`, Industrialized) %>%
  mutate(non_industrialized_late = `non-Industrialized` > 0,
         industrialized_late = Industrialized > 0,
         non_industrialized_early = `non-Industrialized` < 0,
         industrialized_early = Industrialized < 0) %>%
  colSums(., na.rm = T)

(single_ls_industrialized <- shap_cors_filt_xy_plot %>% filter(is.na(`non-Industrialized`)) %>%
  ggplot(., aes(x = Industrialized, fill = Classification, shape = Classification)) +
    geom_segment(x = -1.01, xend = 1.01, y = 0, yend = 0,
                 alpha = 0.1,
                 arrow = arrow(length = unit(0.3, "cm"), type = "closed", ends = "both"),
                 color = "grey") +
    geom_point(size = 3, alpha = 0.9, fill = "#746CFF", shape = 21, y = 0) +
    xlim(-1.05, 1.05) +
    ylim(-0.05, 0.05) +
    theme_void() +
    theme(#plot.margin = margin(b = -10, t = -10, l = 0, r = 0),
          legend.position = "none"))

(single_ls_non_industrialized <- shap_cors_filt_xy_plot %>% filter(is.na(`Industrialized`)) %>%
    ggplot(., aes(y = `non-Industrialized`, fill = Classification, shape = Classification)) +
    geom_segment(y = -1.01, yend = 1.01, x = 0, xend = 0,
                 alpha = 0.1,
                 arrow = arrow(length = unit(0.3, "cm"), type = "closed", ends = "both"),
                 color = "grey") +
    geom_point(size = 3, alpha = 0.9, fill = "#746CFF", shape = 21, x = 0) +
    ylim(-1.05, 1.05) +
    xlim(-0.05, 0.05) +
    theme_void() +
    theme(axis.line.y = element_line(color = "black"),
          plot.margin = margin(l = -10, t = 0, r = -10, b = 0),
          legend.position = "none"))
    
(diff_shap_cors_both <- ggplot(shap_cors_filt_xy_plot, 
                          aes(x = Industrialized, y = `non-Industrialized`, shape = Classification, fill = Classification)) +
    geom_segment(aes(x = -1.01, xend = 1.01, y = 0, yend = 0),
                 arrow = arrow(length = unit(0.3, "cm"), type = "closed", ends = "both"),
                 color = "black") +
    geom_segment(aes(y = -1.01, yend = 1.01, x = 0, xend = 0),
                 arrow = arrow(length = unit(0.3, "cm"), type = "closed", ends = "both"),
                 color = "black") +  
    annotate("text", x = -1, y = 0.1, 
             label = "Early", hjust = 0, size = 5, fontface = "bold") +
    annotate("text", x = 1, y = 0.1, 
             label = "Late", hjust = 1, size = 5, fontface = "bold") +
    annotate("text", y = 1, x = - 0.05, 
             label = "Late", hjust = 1, size = 5, fontface = "bold", angle = 0) +
    annotate("text", y = -1, x = - 0.05, 
             label = "Early", hjust = 1, size = 5, fontface = "bold", angle = 0) +
    geom_point(size = 3, alpha = 0.9) +
    ggrepel::geom_text_repel(aes(label = annotation_taxon),
                             fontface = "italic",
                             point.size = 2,
                             size = 4,
                             seed = 123,
                             max.overlaps = Inf,
                             box.padding = 0.7) +
    scale_fill_manual(values = c(`Single lifestyle` = "#746CFF", `Late colonizer` = "#F08B3C", 
                                 Mixed = "#7A0403", `Early colonizer` = "#7CC235")) +
    scale_shape_manual(values = c(22,24,23,21)) +
    labs(title = NULL,
         x = bquote("Industrialized ["~R^2~"]"),
         y = bquote("non-Industrialized ["~R^2~"]")) +
    # coord_equal() +
    xlim(-1.05, 1.05) +
    ylim(-1.05, 1.05) +
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          plot.title = element_text(size=18),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          axis.line = element_blank(),
          legend.position = "inside",
          legend.position.inside=c(0.8, 0.2),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          # plot.margin = margin(t = -1, r = -1, b = 0.2, l = 1, "cm")
          strip.text = element_text(size = 14)))


hist_x <- axis_canvas(diff_shap_cors_both, axis = 'x') + single_ls_industrialized
hist_y <- axis_canvas(diff_shap_cors_both, axis = 'y', coord_flip = TRUE) + 
  single_ls_non_industrialized

test <- diff_shap_cors_both %>% insert_xaxis_grob(., hist_x, grid::unit(.2, "null"), position = "top") %>%
  insert_yaxis_grob(., hist_y, grid::unit(.2, "null"), position = "right") %>%
  ggdraw()

diff_shap_cors <- ggplot2::ggplotGrob(diff_shap_cors_both)
diff_shap_cors <- ggExtra:::addTopMargPlot(diff_shap_cors, top = single_ls_industrialized, size = 10)
diff_shap_cors <- ggExtra:::addRightMargPlot(diff_shap_cors, right = single_ls_non_industrialized, size = 10)
class(diff_shap_cors) <- c("ggExtraPlot", class(diff_shap_cors))

ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/diff_shap_cors_lifestyle.png", 
       width = 10, plot = diff_shap_cors)



## Prevalences ########################
load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/prevalences.RData")

selected_taxa <- c(selected_taxa) # Bacteroides "Atopobium", "Clostridium.sensu.stricto.1", "Collinsella"

shap_vals_combined %>%
  filter(display_taxon_filter %in% selected_taxa) %>%
  ggplot(.,aes(y = ab_value, x = shap_value, color = model)) +
  geom_point(size = 0.1, alpha = 0.1) +
  geom_smooth() +
  scale_color_manual(values = fixed_colors, 
                     guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1),
                                          direction = "vertical",
                                          position = "bottom",
                                          title.position = "left")) +
  facet_wrap(~display_taxon_filter)


(prevalence_plot <- prev_per_month %>%
    mutate(Lifestyle = ifelse(lifestyle == "industrialized", yes = "Industrialized", no = "non-Industrialized")) %>%
    left_join(., short_taxa_names, by = c("OTU" = "taxon")) %>% #pull(display_taxon) %>% unique %>% sort
    filter(display_taxon %in% selected_taxa) %>%
    ggplot(., aes(x = month, y = prev, color = Lifestyle, fill = Lifestyle)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", alpha = 0) +
    scale_color_manual(values = fixed_colors, 
                       guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1),
                                            direction = "vertical",
                                            position = "bottom",
                                            title.position = "left")) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    scale_fill_manual(values = fixed_colors) +
    facet_wrap(~display_taxon, nrow = 5, ncol = 1) +
    labs(title = NULL,
         x = "Age [months]",
         y = "Prevalence") +
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_blank(),
          strip.text = element_text(size = 14, face = "italic"),
          strip.background = element_blank()))

ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/selected_prevalences.png", 
       width = 10, plot = prevalence_plot)


(fig_5_combined <- cowplot::plot_grid(
  ggdraw(diff_shap_cors), 
  prevalence_plot,
  nrow = 1,
  labels = c("A)", "B)"),
  label_size = 18,
  label_fontface = "plain",
  rel_widths = c(2, 1)))

cowplot::save_plot("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/figure_5_combined.png", 
                   fig_5_combined, dpi = 900,
                   bg = "white", base_height = 8.15, base_width = 14)
cowplot::save_plot("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/figure_5_combined.pdf", 
                   fig_5_combined, dpi = 900,
                   bg = "white", base_height = 8.15, base_width = 14)


# Figure 6 ###########################################
load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/sick_age_predictions.RData")
fixed_colors <- c(Industrialized = "#737125", `non-Industrialized` = "#1A97C8", Combined = "#2A2359")
fixed_colors_sick <- fixed_colors
names(fixed_colors_sick) <- paste0(names(fixed_colors_sick), " healthy")
fixed_colors_sick <- c(fixed_colors_sick, `non-Industrialized SAM` = "#812B8C", `Industrialized preterm` = "#D9731A")

## A: Malnurished microbial age ########################

(mal_pred_age <- all_preds_combined %>% 
    mutate(lifestyle_health = gsub("_", " ", lifestyle_health) %>%
             gsub("industrialized", "Industrialized", .) %>%
             gsub("non ", "non-", .),
           training_set = case_when(training_set == "combined" ~ "Combined",
                                    training_set == "non_industrialized" ~ "non-Industrialized",
                                    training_set == "industrialized" ~ "Industrialized")) %>%
    ggplot(., aes(x = age, y = predicted_age, color = lifestyle_health, fill = lifestyle_health)) +
    geom_point(alpha = 0.2, size = 0.1) +
    geom_smooth(method = "loess") +
   scale_x_continuous(limits = c(NA, 630), breaks = c(200, 400, 600)) +
    scale_color_manual(values = fixed_colors_sick,
                       name = "Lifestyle and health",
                       guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
    scale_fill_manual(values = fixed_colors_sick,
                      name = "Lifestyle and health") +
   labs(title = NULL,
        x = "Chronological age [days]",
        y = "Microbial age [days]") +
    facet_wrap(~training_set) +
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 13),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          strip.text = element_text(size = 14),
          strip.background = element_blank(),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)))


## B: Malnurished MAZ ########################
means_maz_mal <- all_preds_combined %>% 
  group_by(training_set, lifestyle_health) %>%
  dplyr::summarize(mean = mean(MAZ))

pvals_maz_mal <- all_preds_combined %>%
  group_by(training_set) %>%
  rstatix::wilcox_test(MAZ ~ lifestyle_health) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>%
  rstatix::add_xy_position(x = "training_set", dodge = 1) %>%
  left_join(.,means_maz_mal, by = c("training_set" = "training_set",
                                    "group1" = "lifestyle_health")) %>%
  left_join(.,means_maz_mal, by = c("training_set" = "training_set",
                                    "group2" = "lifestyle_health")) %>%
  mutate(MAZ_diff = abs(mean.x - mean.y)) %>%
  select(MAZ_diff, everything())

(mal_maz_diff <- all_preds_combined %>%
    mutate(lifestyle_health = gsub("_", " ", lifestyle_health) %>%
             gsub("industrialized", "Industrialized", .) %>%
             gsub("non ", "non-", .),
           training_set = case_when(training_set == "combined" ~ "Combined",
                                    training_set == "non_industrialized" ~ "non-Industrialized",
                                    training_set == "industrialized" ~ "Industrialized")) %>%
    ggplot(., aes(x = training_set, y = MAZ)) +
    stat_eye(aes(fill = lifestyle_health), position = position_dodge(1), scale = 0.9,
                 .width = c(0, 0.5, 0.95), adjust = 1, shape = 23, point_size = 2,
                 side = "both") +
    scale_y_continuous(breaks = c(-5, -2.5, 0, 2.5, 5)) +
    scale_fill_manual(values = fixed_colors_sick,
                      name = "Lifestyle and health",
                      guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
    add_pvalue(pvals_maz_mal,
               step.increase = 0.05 ,
               label = "{p.adj.signif}",
               label.size = 4.5,
               tip.length = 0.01,
               xmin = "xmin",
               xmax = "xmax",
               show.legend = FALSE) +
    labs(title = NULL,
         x = "Training set") +
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 13),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 18),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)))

## C: Preterms microbial age ###########################
(preterms_pred_age <- all_preds_combined_preterm %>% 
   mutate(lifestyle_health = gsub("_", " ", lifestyle_health) %>%
            gsub("industrialized", "Industrialized", .) %>%
            gsub("non ", "non-", .),
          training_set = case_when(training_set == "combined" ~ "Combined",
                                   training_set == "non_industrialized" ~ "non-Industrialized",
                                   training_set == "industrialized" ~ "Industrialized")) %>%
   ggplot(., aes(x = age, y = predicted_age, color = lifestyle_health, fill = lifestyle_health)) +
   geom_point(alpha = 0.2, size = 0.1) +
   geom_smooth(method = "loess") + 
   scale_color_manual(values = fixed_colors_sick,
                      name = "Lifestyle and health",
                      guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
   scale_fill_manual(values = fixed_colors_sick,
                     name = "Lifestyle and health") +
   labs(title = NULL,
        x = "Chronological age [days]",
        y = "Microbial age [days]") +
   xlim(0, NA) +
   facet_wrap(~training_set) +
   theme(axis.title.x = element_text(size = 16),
         axis.text.x = element_text(size = 13),
         axis.title.y = element_text(size = 16),
         axis.text.y = element_text(size = 14),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         legend.position = "none",
         strip.text = element_text(size = 14),
         strip.background = element_blank(),
         legend.text = element_text(size = 15),
         legend.title = element_text(size = 15)))

## D: Preterms MAZ ###########################

means_maz_preterm <- all_preds_combined_preterm %>% 
  group_by(training_set, lifestyle_health) %>%
  dplyr::summarize(mean = mean(MAZ))

pvals_maz_preterm <- all_preds_combined_preterm %>%
  group_by(training_set) %>%
  rstatix::wilcox_test(MAZ ~ lifestyle_health) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>%
  rstatix::add_xy_position(x = "training_set", dodge = 1) %>%
  left_join(.,means_maz_preterm, by = c("training_set" = "training_set",
                                    "group1" = "lifestyle_health")) %>%
  left_join(.,means_maz_preterm, by = c("training_set" = "training_set",
                                    "group2" = "lifestyle_health")) %>%
  mutate(MAZ_diff = abs(mean.x - mean.y)) %>%
  select(MAZ_diff, everything())

(preterm_maz_diff <- all_preds_combined_preterm %>%
    mutate(lifestyle_health = gsub("_", " ", lifestyle_health) %>%
             gsub("industrialized", "Industrialized", .) %>%
             gsub("non ", "non-", .),
           training_set = case_when(training_set == "combined" ~ "Combined",
                                    training_set == "non_industrialized" ~ "non-Industrialized",
                                    training_set == "industrialized" ~ "Industrialized")) %>%
    ggplot(., aes(x = training_set, y = MAZ)) +
    # geom_boxplot(aes(fill = lifestyle_health), outliers = F) +
    stat_eye(aes(fill = lifestyle_health), position = position_dodge(1), scale = 0.9,
                 .width = c(0, 0.5, 0.95), adjust = 1, shape = 23, point_size = 2,
                 side = "both") +
    scale_fill_manual(values = fixed_colors_sick,
                      name = "Lifestyle and health",
                      guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
    add_pvalue(pvals_maz_mal,
               label = "{p.adj.signif}",
               step.increase = 0.05 ,
               # bracket.nudge.y = 3,
               label.size = 4.5,
               tip.length = 0.01,
               xmin = "xmin",
               xmax = "xmax",
               show.legend = FALSE) +
    labs(title = NULL,
         x = "Training set") +
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 13),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)))

## E: Malnurished SHAP ###############################
load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/sick_shap_data.RData")
pvals_svz_mal_age <- shap_mdat_mal %>%
  mutate(age_bin = floor(ceiling(age) / 30)) %>%
  group_by(display_taxon, age_bin) %>%
  rstatix::wilcox_test(shap_value ~ health) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1))
interesting_taxa_mal <- pvals_svz_mal_age %>% group_by(display_taxon) %>%
  summarize(n_signif = mean(p.adj.signif != "ns")) %>% #pull(n_signif) %>%table
  filter(n_signif >= 0.8) %>% pull(display_taxon)

mean_shap_diff_mal <- shap_mdat_mal %>%
  mutate(age_bin = floor(ceiling(age) / 30)) %>%
  left_join(.,pvals_svz_mal_age,
            by = c("display_taxon" = "display_taxon", "age_bin" = "age_bin")) %>%
  mutate(p.adj.signif = ifelse(p.adj.signif == "ns", yes = NA, no = p.adj.signif)) %>%
  group_by(display_taxon, age_bin) %>%
  summarize(mean_shap_healthy = mean(shap_value[health == "healthy"], na.rm = T),
            mean_shap_sam = mean(shap_value[health == "SAM"], na.rm = T),
            diff_shap = mean_shap_sam - mean_shap_healthy,
            significance = unique(p.adj.signif))# %>% pull(diff_SVZ) %>% sum()

high_effect_taxa_mal <- mean_shap_diff_mal %>% filter(abs(diff_shap) > 10) %>% pull(display_taxon) %>% unique

(heatmap_mal <- mean_shap_diff_mal %>%
    filter(display_taxon %in% intersect(interesting_taxa, high_effect_taxa_mal)) %>%
    mutate(display_taxon = if_else(grepl("_| ", display_taxon) |
                                     display_taxon == "Shannon" |
                                     display_taxon == "Observed",
                                   true = glue("{display_taxon}"),
                                   false = glue("<i>{display_taxon}</i>")) %>%
             gsub("_", ".",.),
           display_taxon = factor(display_taxon, levels = unique(c("Observed", "Shannon", display_taxon)))) %>%
    ggplot(., aes(x = age_bin, y = display_taxon, fill = diff_shap)) +
    geom_tile(color = "white") +                   # White borders between tiles
    geom_text(aes(label = significance), color = "black", size = 3, angle = 90) +  # significance stars
    scale_x_continuous(breaks = c(6, 12, 18)) +
    scale_fill_gradient2(name = "\U0394 mean\nnon-Industrialized\nSHAP",
                         mid = "white",
                         low = c("#4662D7FF", "#36AAF9FF", "#1AE4B6FF"),
                         high = c("#FABA39FF", "#F66B19FF", "#CB2A04FF")) +
    labs(title = "SAM",
         x = "Chronological age [months]") +
    theme(axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.text.y = element_markdown(size = 12),
          axis.title.y = element_blank(),
          title = element_text(size = 14),
          plot.title = element_text(size=18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_blank(),
          strip.text = element_text(size = 14)))

## F: Preterm SHAP ########################################
pvals_svz_preterm_age <- shap_mdat_preterm %>%
  mutate(age_bin = floor(age / 7),
         age_bin = ifelse(age_bin == 12, yes = 11, no = age_bin)) %>% # last bin contains too few values for tests
  group_by(display_taxon, age_bin) %>%
  rstatix::wilcox_test(shap_value ~ health) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1))

interesting_taxa_pret <- pvals_svz_preterm_age %>% group_by(display_taxon) %>%
  summarize(n_signif = mean(p.adj.signif != "ns")) %>% #pull(n_signif) %>%table
  filter(n_signif >= 0.1) %>% pull(display_taxon)

mean_shap_diff_preterm <- shap_mdat_preterm %>%
  mutate(age_bin = floor(age / 7),
         age_bin = ifelse(age_bin == 12, yes = 11, no = age_bin)) %>% # last bin contains too few values for tests
  left_join(.,pvals_svz_preterm_age,
            by = c("display_taxon" = "display_taxon", "age_bin" = "age_bin")) %>%
  mutate(p.adj.signif = ifelse(p.adj.signif == "ns", yes = NA, no = p.adj.signif),
         health_study = ifelse(health == "preterm", yes = paste(health, study, sep = "_"), no = health)) %>%
  group_by(display_taxon, age_bin) %>%
  summarize(mean_shap_healthy = mean(shap_value[health == "healthy"], na.rm = T),
            mean_shap_preterm = mean(shap_value[health == "preterm"], na.rm = T),
            diff_shap = mean_shap_preterm - mean_shap_healthy,
            significance = unique(p.adj.signif)) 

high_effect_taxa <- mean_shap_diff_preterm %>% filter(abs(diff_shap) > 3) %>% pull(display_taxon) %>% unique

(heatmap_preterm <- mean_shap_diff_preterm %>%
    filter(display_taxon %in% intersect(interesting_taxa_pret, high_effect_taxa)) %>%
    mutate(display_taxon = if_else(grepl("_| ", display_taxon) |
                                     display_taxon == "Shannon" |
                                     display_taxon == "Observed",
                                   true = glue("{display_taxon}"),
                                   false = glue("<i>{display_taxon}</i>")) %>%
             gsub("_", ".",.),
           display_taxon = factor(display_taxon, levels = unique(c("Observed", display_taxon)))) %>%
    # group_by(age_bin) %>% summarise(sum_shap = sum(diff_shap, na.rm = T))
    ggplot(., aes(x = age_bin, y = display_taxon, fill = diff_shap)) +
    geom_tile(color = "white") +                   # White borders between tiles
    geom_text(aes(label = significance), color = "black", size = 3, angle = 90) +  # significance stars
    scale_x_continuous(breaks = c(0, 4, 8, 11)) +
    scale_fill_gradient2(name = "\U0394 mean\nIndustrialized\nSHAP",
                         mid = "white",
                         low = c("#4662D7FF", "#36AAF9FF", "#1AE4B6FF"),
                         high = c("#FABA39FF", "#F66B19FF", "#CB2A04FF")) +
    labs(title = "Preterm",
         x = "Chronological age [weeks]") +
    theme(axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.text.y = element_markdown(size = 12),
          axis.title.y = element_blank(),
          title = element_text(size = 14),
          plot.title = element_text(size=18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_blank()))


(legend_plot  <- data.frame(lifestyle_health = names(fixed_colors_sick)) %>% 
    filter(lifestyle_health != "Combined healthy") %>%
    ggplot(., aes(x = 1, y = 2)) +
    geom_point(aes(color = lifestyle_health)) +
    scale_color_manual(values = fixed_colors_sick,
                       name = "Lifestyle and health",
                       guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1),
                                            nrow = 1)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)))

# combine and save plot
cowplot::set_null_device("agg") # pdf, png, cairo, agg (on max-cluster)
fig_6_combined <- cowplot::plot_grid(cowplot::plot_grid(mal_pred_age, mal_maz_diff, 
                                                        preterms_pred_age, preterm_maz_diff,
                                                        nrow = 2,
                                                        rel_widths = c(1.5, 1),
                                                        labels = c("A)", "B)", "C)", "D)"),
                                                        label_size = 18,
                                                        label_y = 1.02,
                                                        label_x = c(0, 0, -0.005, -0.01),
                                                        align = "vh",
                                                        axis = "tblr",
                                                        label_fontface = "plain"),
                                     # NULL,
                                     # cowplot::plot_grid(preterms_pred_age, preterm_maz_diff, 
                                     #                    nrow = 1,
                                     #                    rel_widths = c(1.5, 1),
                                     #                    labels = c("C)", "D)"),
                                     #                    label_size = 18,
                                     #                    label_y = 1.02,
                                     #                    label_fontface = "plain"),
                                     get_legend(legend_plot),
                                     NULL,
                                     cowplot::plot_grid(heatmap_mal, heatmap_preterm,
                                                        rel_widths = c(1, 1.2),
                                                        labels = c("E)", "F)"),
                                                        label_size = 18,
                                                        label_y = 1.02,
                                                        label_x = c(0.005, 0),
                                                        label_fontface = "plain"),
                                     ncol = 1, rel_heights = c(2, 0.2, 0.1, 1.2))


cowplot::save_plot("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/figure_6_combined.png", 
                   fig_6_combined, dpi = 900,
                   bg = "white", base_height = 11.5, base_width = 14)
cowplot::save_plot("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/figure_6_combined.pdf", 
                   fig_6_combined, dpi = 900,
                   bg = "white", base_height = 11.5, base_width = 14)