library(ggplot2)
library(tidyverse)
library(viridis)
library(gridExtra)
library(grid)
library(ggpubr)
library(ggpattern)
library(patchwork)
library(cowplot)
library(rstatix)
library(ggprism)
library(ggh4x)
library(grid)
library(ggtext)
library(glue)
library(phyloseq)
library(ggVennDiagram)
library(ggvenn)
library(ggforce)
library(wesanderson)
library(ggdist)
library(ggbeeswarm)

# to do:
# maybe: abundance over time
cowplot::set_null_device("agg")


source("/fast/AG_Forslund/rob/mm_index/R_scripts/functions.R")
setwd("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement")
# define global colors for lifestyle:
# fixed_colors <- setNames(viridisLite::turbo(3), c("Industrialized", "non-Industrialized", "Combined")) 
fixed_colors <- c(Industrialized = "#737125", `non-Industrialized` = "#1A97C8", Combined = "#BC85A9")

# Figure S1 #####################################################################
# rarefaction curves
load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/rarefaction_curves.RData")
raredat <- raredat %>%
  mutate(study = paste0(gsub("^(.)", "\\U\\1",
                             gsub("_", " et al. ", study),
                             perl=TRUE),
                        "; ",
                        cutoff)) %>%
  filter(Sample < 21000)


(rarefaction_curves <- ggplot(raredat) +
# scale_color_discrete(guide = "none") +
geom_line(aes(x = Sample,
              y = Species,
              color = sample_sum,
              # color = factor(as.numeric(factor(Site)) %% 7),
              group = Site),
          alpha = 0.5) +
geom_vline(aes(xintercept = cutoff),
           linetype = "dashed",
           color = "red", size = 0.8) +
scale_x_continuous(n.breaks = 3) +
scale_y_continuous(n.breaks = 4) +
scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous"),
                      trans = "log10",
                      breaks = c(1e+1, 1e+3, 1e+5),
                      name = "Total reads") +
guides(color = guide_colorbar(
  direction = "horizontal",
  title.position = "top",
  title.hjust = 1,
  barwidth = 12,
  barheight = 2)) +
labs(title = NULL,
     x = "Number of reads",
     y = "Number of genera detected") +
facet_wrap(~ study, scales = "free", ncol = 3) +
theme(axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      strip.text = element_text(size = 14),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      legend.position = c(0.85, 0.05),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()))

ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/figure_S1_combined_supplement.png",
       dpi = 900,
       width = 14, height = 12.72, plot = rarefaction_curves)



# Figure S2 #####################################################################
## A: Prevalence Venn diagram ###########################################
load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/prevalence_lifestyle.RData")
prevalence_lifestyle_g_bin <- prevalence_lifestyle_g %>%
  column_to_rownames("lifestyle") %>%
  t() %>%
  `>` (0) %>%
  data.frame()
venn_lists <- list(
  `non-Industrialized` = rownames(prevalence_lifestyle_g_bin)[prevalence_lifestyle_g_bin$non_industrialized],
  Industrialized = rownames(prevalence_lifestyle_g_bin)[prevalence_lifestyle_g_bin$industrialized]
)
union <- unique(unlist(venn_lists)) %>% length
intersection <- intersect(venn_lists$`non-Industrialized`, venn_lists$Industrialized) %>% length
only_non_industrial <- sum(!(venn_lists$`non-Industrialized` %in% venn_lists$Industrialized))
only_industrial <- sum(!(venn_lists$Industrialized %in% venn_lists$`non-Industrialized`))

circle_data <- data.frame(
  x0 = c(-0.65, 0.65),  # Adjust to increase/decrease overlap
  y0 = c(0, 0),
  r = c(1, 1),
  labels = c("non-Industrialized", "Industrialized"),
  counts = c(paste0(only_non_industrial, "\n", round(100* only_non_industrial / union, digits = 0), "%"),
             paste0(only_industrial, "\n", round(100* only_industrial / union, digits = 0), "%"))  # Replace with actual numbers
)

(prev_venn <- ggplot() +
    geom_circle(data = circle_data,
      aes(x0 = x0, y0 = y0, r = r, fill = labels),
      color = "black", alpha = 0.6, size = 0.5) +
    coord_fixed() +
    scale_fill_manual(values = unname(fixed_colors)) +
    # number/percentage
    annotate("text",
      x = circle_data$x0[1], y = 0,
      label = circle_data$counts[1], size = 5) +
    annotate("text",
      x = circle_data$x0[2], y = 0,
      label = circle_data$counts[2], size = 5) +
    annotate("text",
             x = 0, y = 0, 
             label = paste0(intersection, "\n", round(100* intersection / union, digits = 0), "%"), size = 4) +
    # lifestyle labels
    annotate("text",
             x = circle_data$x0[1] * 1.2, y = 1.2,
             label = circle_data$labels[1], size = 4.5) +
    annotate("text",
             x = circle_data$x0[2] * 1.2, y = 1.2,
             label = circle_data$labels[2], size = 4.5) +
    
    # ggtitle("A)") +
    theme_void() +
    theme(plot.title = element_text(size = 18),
          legend.position = "none"))
### on rarefied counts: ##############################
prevalence_lifestyle_g_bin_raref <- prevalence_lifestyle_g_raref %>%
  column_to_rownames("lifestyle") %>%
  t() %>%
  `>` (0) %>%
  data.frame()
venn_lists_raref <- list(
  `non-Industrialized` = rownames(prevalence_lifestyle_g_bin_raref)[prevalence_lifestyle_g_bin_raref$non_industrialized],
  Industrialized = rownames(prevalence_lifestyle_g_bin_raref)[prevalence_lifestyle_g_bin_raref$industrialized]
)

union_raref <- unique(unlist(venn_lists_raref)) %>% length
intersection_raref <- intersect(venn_lists_raref$`non-Industrialized`, venn_lists_raref$Industrialized) %>% length
only_non_industrial_raref <- sum(!(venn_lists_raref$`non-Industrialized` %in% venn_lists_raref$Industrialized))
only_industrial_raref <- sum(!(venn_lists_raref$Industrialized %in% venn_lists_raref$`non-Industrialized`))

circle_data_raref <- data.frame(
  x0 = c(-0.65, 0.65),  # Adjust to increase/decrease overlap
  y0 = c(0, 0),
  r = c(1, 1),
  labels = c("non-Industrialized", "Industrialized"),
  counts = c(paste0(only_non_industrial_raref, "\n", round(100* only_non_industrial_raref / union_raref, digits = 0), "%"),
             paste0(only_industrial_raref, "\n", round(100* only_industrial_raref / union_raref, digits = 0), "%"))  # Replace with actual numbers
)

(prev_venn_raref <- ggplot() +
    geom_circle(data = circle_data_raref,
                aes(x0 = x0, y0 = y0, r = r, fill = labels),
                color = "black", alpha = 0.6, size = 0.5) +
    coord_fixed() +
    scale_fill_manual(values = unname(fixed_colors)) +
    # number/percentage
    annotate("text",
             x = circle_data_raref$x0[1], y = 0,
             label = circle_data_raref$counts[1], size = 5) +
    annotate("text",
             x = circle_data_raref$x0[2], y = 0,
             label = circle_data_raref$counts[2], size = 4) +
    annotate("text",
             x = 0, y = 0, 
             label = paste0(intersection_raref, "\n", round(100* intersection_raref / union_raref, digits = 0), "%"), size = 4) +
    # lifestyle labels
    annotate("text",
             x = circle_data_raref$x0[1] * 1.2, y = 1.2,
             label = circle_data_raref$labels[1], size = 4) +
    annotate("text",
             x = circle_data_raref$x0[2] * 1.2, y = 1.2,
             label = circle_data_raref$labels[2], size = 4) +
    
    # ggtitle("A)") +
    theme_void() +
    theme(plot.title = element_text(size = 18),
          legend.position = "none"))
cowplot::plot_grid(prev_venn, prev_venn_raref)


## B: sequencing depth ################
load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/all_pcoa_genus_combined.RData")
plot_pcoa_data <- pcoa_aitch_genus_1$x %>% data.frame() %>% 
  rownames_to_column(., "run_accession") %>%
  left_join(., pcoa_metadata, by = "run_accession") %>%
  mutate(Lifestyle = ifelse(lifestyle == "industrialized", yes = "Industrialized",
                            no = "non-Industrialized"),
         Study = gsub("_", " et. al ", study) %>%
           str_replace(., "^\\w{1}", toupper))

(sequencing_depth <- ggplot(plot_pcoa_data, aes(x = Lifestyle, y = sample_sum, fill = Lifestyle)) +
    stat_eye(position = position_dodge(1), scale = 0.9,
             .width = c(0, 0.5, 0.95), adjust = 1, shape = 23, point_size = 3, show.legend = F) +
    scale_y_log10() +
    geom_hline(yintercept = 10000, linetype = "dashed") +
    labs(title = NULL,
         x = NULL,
         y = "Total read count") +
    scale_fill_manual(values = fixed_colors) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 13),
          axis.text.y = element_text(size = 11),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)))


## C: richness ##############
load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/a_div_genus_present.RData")
(a_div_s_genus_richness <- ggplot(meta_df %>% mutate(Lifestyle = ifelse(lifestyle == "industrialized",
                                                               yes = "Industrialized",
                                                               no = "non-Industrialized")),
                         aes(y = observed, x = age, color = Lifestyle, fill = Lifestyle)) +
    geom_point(alpha = 0.2, size = 0.3) +
    # geom_smooth(aes(fill = Lifestyle), method = "loess", size = 2, color = "black") +
    geom_smooth(method = "loess", size = 1) +
    scale_color_manual(values = fixed_colors, 
                       guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
    scale_fill_manual(values = fixed_colors) +
    labs(title = NULL,
         x = "Age [Days]",
         y = "Richness (Genus)") +
    # guides(color = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) + 
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = c(0.5, 0.8),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)))

ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/a_div_genus_richness_present.png",
       width = 11, height = 6, plot = a_div_s_genus_richness)

## D unassigned ASVs ##########################
(a_div_unassigned_asvs <- ggplot(meta_df %>% mutate(Lifestyle = ifelse(lifestyle == "industrialized",
                                                                        yes = "Industrialized",
                                                                        no = "non-Industrialized")),
                                  aes(y = n_unassigned_asvs_sample_raref/richness_asv, x = age, color = Lifestyle, fill = Lifestyle)) +
   geom_point(alpha = 0.2, size = 0.3) +
   # geom_smooth(aes(fill = Lifestyle), method = "loess", size = 2, color = "black") +
   geom_smooth(method = "loess", size = 1) +
   scale_color_manual(values = fixed_colors) +
   scale_fill_manual(values = fixed_colors) +
   labs(title = NULL,
        x = "Age [Days]",
        y = "Fraction of unassigned ASVs") +
   # guides(color = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) + 
   theme(axis.title.x = element_text(size = 16),
         axis.text.x = element_text(size = 14),
         axis.title.y = element_text(size = 16),
         axis.text.y = element_text(size = 14),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size = 15),
         legend.title = element_text(size = 15)))

ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/a_div_genus_richness_present.png",
       width = 11, height = 6, plot = a_div_s_genus_richness)


# (a_div_diff <- bind_rows(alpha_div_sliding_window, alpha_div_sliding_window_no_raman) %>%
#    mutate(Samples = n_total,
#           p_val = ifelse(p.adj < 0.05, yes = "< 0.05", no = "> 0.05")) %>%
#    ggplot(aes(x = sliding_window, y = estimate, color = set, alpha = p_val, 
#               group = set,
#               size = Samples)) +
#    geom_point() +
#    geom_line(size = 0.4, alpha = 0.5) +
#    geom_vline(xintercept = 226) +
#    scale_color_manual(
#      name = NULL,   # This also sets the legend title
#      values = c(
#        "Complete dataset" = "#4777EFFF",
#        "No Raman et al. 2019" = "#DB3A07FF"),
#      guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
#    scale_alpha_manual(name = "p-Value",
#                       values = c("< 0.05" = 1,
#                                  "> 0.05" = 0.2)) +
#      # name = "p-Value",   # This also sets the legend title
#      # values = c(
#      #   "< 0.05" = "#4777EFFF",
#      #   "> 0.05" = "#DB3A07FF"
#      # )) +
#    guides(color = guide_legend(order = 1, override.aes = list(shape = 15, size = 5, alpha = 1)),
#           alpha = guide_legend(oder = 2),
#           size = guide_legend(order = 3)) +
#    ylab("\U0394 Shannon Diversity [W - NW]") +
#    xlab("Chronological age [Days]") +
#    # facet_wrap(~set, nrow = 2) +
#    theme(axis.title.x = element_text(size = 16),
#         axis.text.x = element_text(size = 14),
#         axis.title.y = element_text(size = 16),
#         axis.text.y = element_text(size = 14),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         panel.grid.major.y = element_line(color = "grey",
#                                           size = 0.3,
#                                           linetype = 2),
#         strip.text = element_text(size = 14),
#         strip.background = element_blank(),
#         plot.title = element_text(size = 18),
#         legend.position = c(0.65, 0.85),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 12),
#         legend.direction = "vertical", legend.box = "horizontal"))

## E: drivers of separation ###############

metavar_loadings <- plot_pcoa_data %>% 
  mutate(lifestyle_industrialized = lifestyle == "industrialized") %>%
  select(age, lifestyle_industrialized, diversity_shannon, sample_sum, PC1, PC2) %>%
  cor(., method = "spearman") %>% data.frame %>% rownames_to_column(var = "display_taxon") %>% 
  select(display_taxon, PC1, PC2) %>% 
  filter(display_taxon %in% c("age", "lifestyle_industrialized", "diversity_shannon", "sample_sum")) %>%
  mutate(display_taxon = gsub("lifestyle_industrialized", "Industrialized lifestyle", display_taxon) %>%
           gsub("diversity_shannon", "Shannon diversity", .) %>%
           gsub("sample_sum", "Sequencing depth", .),
         display_taxon = reorder(display_taxon, abs(PC1), FUN = max, na.rm = T)) %>%
  pivot_longer(cols = c("PC1", "PC2"), names_to = "PC", values_to = "loadings") %>%
  mutate(sign = ifelse(loadings > 0, yes = "positive", no = "negative"),
         loadings = abs(loadings))
  
(pcoa_loadings <- ggplot(metavar_loadings, aes(y = display_taxon, x = loadings, fill = sign)) +
  geom_col(position = "dodge") +
  geom_vline(xintercept = 0) +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2)) +
  scale_fill_manual(values = c("#DBD547", "#6279AD"), name = NULL) +
  labs(title = NULL,
       x = NULL,
       y = NULL) +
    facet_wrap(~PC) +
  theme(axis.text.y = element_markdown(size = 12),#, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",          # move *all* legends below
        panel.grid.major.x = element_line(color = "grey",
                                          size = 0.3,
                                          linetype = 2),
        panel.grid.major.y = element_line(color = "grey",
                                          size = 0.3,
                                          linetype = 2),
        strip.text = element_text(size = 14),
        strip.background = element_blank()))


(pc_drivers <- top_drivers_1 %>%
    mutate(display_taxon = if_else(grepl("_| ", Taxon),
                                  true = glue("{Taxon}"),
                                  false = glue("<i>{Taxon}</i>")) %>%
            gsub("Bacteria_Firmicutes_Clostridia_unclassified_Clostridia_unclassified_Clostridia_unclassified_Clostridia",
                 "unclassified_Clostridia",.) %>%
            gsub("unclassified_", "Uncl. ",.),
           display_taxon = factor(display_taxon),
           display_taxon = reorder(display_taxon, abs(PC1), FUN = max, na.rm = T)) %>%
    # bind_rows(metavar_loadings) %>%
    pivot_longer(cols = c("PC1", "PC2"), names_to = "PC", values_to = "loadings") %>%
    mutate(sign = ifelse(loadings > 0, yes = "positive", no = "negative"),
           loadings = abs(loadings),
           top_10 = ifelse(PC == top_10, yes = 1, no = 0.3)) %>%
    ggplot(., aes(y = display_taxon, x = loadings, fill = sign, alpha = top_10)) +
    geom_col(position = "dodge") +
    geom_vline(xintercept = 0) +
    scale_x_continuous(breaks = c(-0.2, 0, 0.2)) +
    scale_alpha_identity() +
    scale_fill_manual(values = c("#DBD547", "#6279AD"), name = NULL,
                      guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
    labs(title = NULL,
         x = "Loadings",
         y = NULL) +
    facet_wrap(~PC) +
    theme(axis.text.y = element_markdown(size = 12),#, angle = 45, vjust = 1, hjust=1),
         axis.title.x = element_text(size = 16),
         axis.text.x = element_text(size = 12),
         axis.title.y = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         legend.position = "bottom",          # move *all* legends below
         legend.box = "horizontal",           # arrange them horizontally
         panel.grid.major.x = element_line(color = "grey",
                                           size = 0.3,
                                           linetype = 2),
         panel.grid.major.y = element_line(color = "grey",
                                           size = 0.3,
                                           linetype = 2),
         legend.text = element_text(size = 15),
         legend.title = element_text(size = 15),
         strip.text = element_blank(),
         strip.background = element_blank()))

## F: Beta div study ##############
explained_var <- (pcoa_aitch_genus_1$sdev^2 / sum(pcoa_aitch_genus_1$sdev^2)) %>%
  `*` (100) %>%
  round(.,digits = 1)

(pcoa_study <- ggplot(plot_pcoa_data, aes(x = PC1, y = PC2, fill = study, color = Lifestyle)) +
    geom_point(alpha = 0.8, size = 0.6, shape = 21, color = "transparent") +
    stat_ellipse(alpha = 0.5) +
    scale_fill_manual(values = viridisLite::turbo(20)) +
    scale_color_manual(values = fixed_colors) +
    guides(fill = guide_legend(override.aes = list(size = 3, alpha = 1, shape = 22), ncol = 2),
           color = "none") +
    coord_equal() +
    labs(title = NULL,
         x = paste0("PC 1 (", explained_var[1], "%)"),
         y = paste0("PC 2 (", explained_var[2], "%)"),
         fill = "Study") +
    annotate("text", x = min(plot_pcoa_data$PC1) * 0,  # Adjust position
             y = max(plot_pcoa_data$PC2) * 0.9,
             # label = expression(atop(R^2 ~ "= 0.07", "p < 0.001")), size = 4, hjust = 1) +
             label = "R^2 == 0.07*';'~p < 0.001", size = 4, hjust = 1, parse = T) +
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.text = element_text(size = 15, margin = margin(r = 13)),
          legend.title = element_text(size = 15)))
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/all_pcoa_genus_study.png",
       width = 14, height = 7, plot = pcoa_study)


## G: Beta div rarefied lifestyle ######################
explained_var_rare <- (pcoa_aitch_genus_rare_h$sdev^2 / sum(pcoa_aitch_genus_rare_h$sdev^2)) %>%
  `*` (100) %>%
  round(.,digits = 1)
plot_pcoa_data_rare <- pcoa_aitch_genus_rare_h$x %>% as.data.frame() %>%
  rownames_to_column("run_accession") %>%
  left_join(., pcoa_metadata, by = "run_accession") %>%
  mutate(Lifestyle = ifelse(lifestyle == "industrialized", yes = "Industrialized",
                            no = "non-Industrialized"))

(pcoa_lifestyle_rare <- ggplot(plot_pcoa_data_rare, aes(x = PC1, y = PC2, color = Lifestyle)) +
    geom_point(alpha = 0.8, size = 0.3) +
    scale_color_manual(values = fixed_colors) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    coord_equal() +
    labs(title = NULL,
         x = paste0("PC 1 (", explained_var_rare[1], "%)"),
         y = paste0("PC 2 (", explained_var_rare[2], "%)")) +
    annotate("text",
             x = max(plot_pcoa_data_rare$PC1) * 0.95,  # Adjust position
             y = max(plot_pcoa_data_rare$PC2) * 0.9,
             # label = expression(atop(R^2 ~ "= 0.07", "p < 0.001")), size = 4, hjust = 1) +
             label = "R^2 == 0.06*';'~p < 0.001", size = 4, hjust = 1, parse = T) +
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          legend.text = element_text(size = 13, margin = margin(r = 13)),
          legend.title = element_text(size = 13)))


## H: Beta div rarefied age ######################
(pcoa_age_rare <- ggplot(plot_pcoa_data_rare, aes(x = PC1, y = PC2, color = age)) +
    geom_point(alpha = 0.8, size = 0.3) +
    # scale_color_viridis_c(option = "C", guide = guide_colorbar(direction = "vertical",
    #                                                            display = "gradient",
    #                                                            title.position = "top",
    #                                                            barwidth = 1.5)) +
    scale_color_gradientn(guide = guide_colorbar(direction = "vertical",
                                                 display = "gradient",
                                                 title.position = "top",
                                                 barwidth = 1.5),
                          colors = wes_palette("Zissou1", type = "continuous"))+
   coord_equal() +
   labs(title = NULL,
        x = paste0("PC 1 (", explained_var_rare[1], "%)"),
        y = paste0("PC 2 (", explained_var_rare[2], "%)"),
        color = "Age [days]") +
   
    annotate(geom = "text", 
             x = max(plot_pcoa_data_rare$PC1) * 0.95,  # Adjust position
             y = max(plot_pcoa_data_rare$PC2) * 0.9,
             label = "R\u00B2 = 0.09; p < 0.001",
             size = 4, hjust = 1) +
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          # legend.position = c(0.83, 0.9),
          # legend.spacing.y = unit(100, 'pt'),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)))




legend_study <- get_legend(pcoa_study)
cowplot::set_null_device("agg")
(fig_S2_combined <- cowplot::plot_grid(cowplot::plot_grid(cowplot::plot_grid(cowplot::plot_grid(prev_venn, sequencing_depth,
                                                                                                rel_widths = c(1, 0.5),
                                                                                                labels = c("A)", "B)"),
                                                                                                label_size = 18,
                                                                                                label_x = c(0, -0.1),
                                                                                                label_fontface = "plain"),
                                                                             a_div_s_genus_richness,
                                                                             a_div_unassigned_asvs,
                                                          rel_heights = c(1, 1.8, 1.8), 
                                                          nrow = 3,
                                                          labels = c(NA, "C)", "D)"),
                                                          label_y = 1.05,
                                                          label_size = 18,
                                                          label_fontface = "plain"),
                                      cowplot::plot_grid(pcoa_loadings, pc_drivers, ncol = 1, rel_heights = c(0.26, 1),
                                                         align = "v", axis = "l", greedy = TRUE),
                                      nrow = 1,
                                      labels = c(NA, "E)"),
                                      label_size = 18,
                                      label_fontface = "plain"),
                                      cowplot::plot_grid(pcoa_study + theme(legend.position='hidden'),
                                                         legend_study,
                                                         nrow = 1,
                                                         labels = c("F)", ""),
                                                         label_size = 18,
                                                         label_y = 1.03,
                                                         label_fontface = "plain"),
                   cowplot::plot_grid(pcoa_lifestyle_rare,
                                      pcoa_age_rare, #pcoa_lifestyle_adapt,
                                      nrow = 1,
                                      rel_widths = c(0.8, 1),
                                      labels = c("G)", "H)"),
                                      label_size = 18,
                                      label_fontface = "plain"),
                   rel_heights = c(1.6, 1, 1), ncol = 1))
cowplot::save_plot("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/figure_S2_combined_supplement.png", 
                   fig_S2_combined, dpi = 900, bg = "white", base_height = 23, base_width = 14)


# Figure S3 Downsampling #######################
## A: Number of features per lifestyle ###########################
load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/downsampling_data.RData")
unique_per_lifestyle <- unique_per_lifestyle %>%
  mutate(Lifestyle = ifelse(category == "industrialized_taxa", yes = "Industrialized", no = "non-Industrialized"))

pvals_unique_per_lifestyle <- unique_per_lifestyle %>% 
  rstatix::wilcox_test(count ~ Lifestyle, paired = T) %>%
  rstatix::add_significance(p.col = "p", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>% 
  rstatix::add_xy_position(x = "Lifestyle", dodge = 0.75,
                           step.increase = 0.6)

(n_unique_resampling <- unique_per_lifestyle %>%
    ggplot(., aes(x = Lifestyle, y = count)) +
    stat_eye(aes(fill = Lifestyle, color = Lifestyle), position = position_dodge(1), scale = 0.9,
             .width = c(0, 0.5, 0.95), adjust = 1, shape = 23, point_size = 3) +
    stat_eye(aes(fill = Lifestyle), position = position_dodge(1), scale = 0.9,
             .width = c(0, 0.5, 0.95), adjust = 1, shape = 23, point_size = 3, show.legend = F) +
    scale_fill_manual(values = fixed_colors,
                      guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
    scale_color_manual(values = fixed_colors, 
                       guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
    ylim(0, NA) +
    add_pvalue(pvals_unique_per_lifestyle,
               label = "{p.signif}") +
    labs(title = NULL,
         x = NULL,
         y = "Unique Genera Count") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "inside",
          legend.position.inside = c(0.4, 0.2),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)))


## B, C, D: Rarefaction of studies, individuals and samples #####################
theme_raefy <- function() {
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size = 18),
    legend.position = "none",
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)
  )
}

(studies_rarefaction <- bind_rows(n_taxa_industrialized_studies_rarefy, n_taxa_non_industrialized_studies_rarefy) %>% 
    mutate(Lifestyle = ifelse(Lifestyle == "Non-industrialized", yes = "non-Industrialized", no = Lifestyle)) %>% #pull(Lifestyle) %>% table()
  ggplot(.,aes(x = n_studies, y = n_taxa, color = Lifestyle)) +
  # geom_point() +
    geom_smooth() +
    ylab("Taxa Count") +
    xlab("Number of Studies") +
  scale_color_manual(values = fixed_colors,
                    guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
  # stat_cor(method = "spearman",
  #          cor.coef.name = "spearman",
  #          size = 4,
  #          show.legend = F) +
  theme_raefy())

(subjects_rarefactin <- bind_rows(n_taxa_industrialized_subjects_rarefy, n_taxa_non_industrialized_subjects_rarefy) %>% 
    mutate(Lifestyle = ifelse(Lifestyle == "Non-industrialized", yes = "non-Industrialized", no = Lifestyle)) %>% #pull(Lifestyle) %>% table()
  ggplot(.,aes(x = n_subjects, y = n_taxa, color = Lifestyle)) +
  # geom_point() +
    geom_smooth() +
    ylab("Taxa Count") +
    xlab("Number of Subjects") +
  scale_color_manual(values = fixed_colors,
                     guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
  theme_raefy())

(samples_rarefaction <- bind_rows(n_taxa_industrialized_samples_rarefy, n_taxa_non_industrialized_samples_rarefy) %>% 
    mutate(Lifestyle = ifelse(Lifestyle == "Non-industrialized", yes = "non-Industrialized", no = Lifestyle)) %>% #pull(Lifestyle) %>% table()
  ggplot(.,aes(x = n_samples, y = n_taxa, color = Lifestyle)) +
  # geom_point() +
    geom_smooth() +
    ylab("Taxa Count") +
    xlab("Number of Samples") +
  scale_color_manual(values = fixed_colors,
                     guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
  theme_raefy())

rarefactions_cobined <- cowplot::plot_grid(studies_rarefaction, subjects_rarefactin, samples_rarefaction)
cowplot::save_plot("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/rarefaction_figures.png", 
                   rarefactions_cobined, bg = "white", base_height = 11, base_width = 7, dpi = 900)

## E: rarefaction performance #####################################
load(file = "/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/lifestyle_perfomance_genus_raref_h.RData")

lifestyle_lm_list_genus_rare_h <- lifestyle_lm_list_genus_rare_h %>%
  mutate(Lifestyle = case_when(lifestyle == "industrialized" ~ "Industrialized",
                               lifestyle == "non_industrialized" ~ "non-Industrialized"),
         Lifestyle = factor(Lifestyle, levels = c("non-Industrialized", "Industrialized")),
         Training_set = case_when(training_set == "industrialized" ~ "Industrialized",
                                  training_set == "non_industrialized" ~ "non-Industrialized",
                                  training_set == "combined" ~ "Combined"),
         Training_set = factor(Training_set, levels = c("non-Industrialized", "Combined", "Industrialized")))


df_p_val_lifestyle_genus_rare_h <- lifestyle_lm_list_genus_rare_h %>%
  arrange(study) %>%
  rstatix::group_by(Lifestyle) %>%
  rstatix::wilcox_test(R2 ~ Training_set, paired = T) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>% 
  rstatix::add_xy_position(x = "Lifestyle", dodge = 0.75,
                           step.increase = 0.6) %>%
  mutate(y.position = ifelse(Lifestyle == "Industrialized",
                             yes = y.position - 0.002,
                             no = y.position))

(downsampling_performance_rare_h <- ggplot(lifestyle_lm_list_genus_rare_h, aes(x=Lifestyle, y = R2)) +
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
         color = "Training Set",
         fill = "Training Set") +
    add_pvalue(df_p_val_lifestyle_genus_rare_h,
               label = "{p.adj.signif}",
               label.size = 4.5,
               tip.length = 0.01,
               xmin = "xmin",
               xmax = "xmax",
               show.legend = FALSE) +
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size=18),
          legend.position="none",
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)))

## F: number of important features #########################
n_imp_fts_per_study <- n_imp_fts_per_study %>%
  ungroup() %>%
  mutate(Lifestyle = ifelse(model_lifestyle == "industrialized", yes = "Industrialized", no = "non-Industrialized"))

df_p_val_n_fts <- n_imp_fts_per_study %>%
  arrange(resample_idx) %>%
  rstatix::wilcox_test(counts ~ Lifestyle, paired = T) %>%
  rstatix::add_significance(p.col = "p", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>% 
  rstatix::add_xy_position()

(n_imp_fts_resampling <- n_imp_fts_per_study %>%
    ggplot(., aes(x = Lifestyle, y = counts)) +
    stat_eye(aes(fill = Lifestyle, color = Lifestyle), position = position_dodge(1), scale = 0.9,
             .width = c(0, 0.5, 0.95), adjust = 1, shape = 23, point_size = 3) +
    stat_eye(aes(fill = Lifestyle), position = position_dodge(1), scale = 0.9,
             .width = c(0, 0.5, 0.95), adjust = 1, shape = 23, point_size = 3, show.legend = F) +
    scale_fill_manual(values = fixed_colors) +
    scale_color_manual(values = fixed_colors) +
    ylim(0, NA) +
    add_pvalue(df_p_val_n_fts,
               label = "{p.signif}") +
    labs(title = NULL,
         x = NULL,
         y = "Important Features Count") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)))
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/n_imp_fts_resampling.png",
       width = 14, height = 7, plot = n_imp_fts_resampling)



## G: Downsampling performance ##########################

lifestyle_lm_list <- lifestyle_lm_list %>%
  mutate(Lifestyle = case_when(lifestyle == "industrialized" ~ "Industrialized",
                               lifestyle == "non_industrialized" ~ "non-Industrialized",
                               lifestyle == "combined" ~ "Combined"),
         Lifestyle = factor(Lifestyle, levels = c("non-Industrialized", "Industrialized")),
         Training_set = case_when(training_set == "pred_industrialized" ~ "Industrialized",
                                  training_set == "pred_non_industrialized" ~ "non-Industrialized",
                                  training_set == "combined" ~ "Combined"),
         Training_set = factor(Training_set, levels = c("non-Industrialized", "Combined", "Industrialized")))


df_p_val_lifestyle_genus <- lifestyle_lm_list %>%
  arrange(study) %>%
  rstatix::group_by(Lifestyle) %>%
  rstatix::wilcox_test(R2 ~ Training_set, paired = T) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>% 
  rstatix::add_xy_position(x = "Lifestyle", dodge = 0.75,
                           step.increase = 0.6) %>%
  mutate(y.position = ifelse(Lifestyle == "Industrialized",
                             yes = y.position - 0.002,
                             no = y.position))

(downsampling_performance <- ggplot(lifestyle_lm_list, aes(x=Lifestyle, y = R2)) +
    # geom_boxplot(aes(fill = Training_set, color = Training_set), outlier.shape = NA) +
    geom_point(aes(color = Training_set), dodge.width = 0.8,
                     shape = 15, width = 0.1, size = 1.5, alpha = 0) +
    geom_boxplot(aes(fill = Training_set), outlier.shape = NA, show.legend = F) +
    geom_quasirandom(aes(fill = Training_set), dodge.width = 0.8,
                     shape = 21, width = 0.1, size = 1.5, show.legend = F) +
    scale_color_manual(values = fixed_colors,
                      guide = guide_legend(override.aes = list(size = 5, shape = 15, alpha = 1),
                                           title = "Training Set")) +
    scale_fill_manual(values = fixed_colors,
                      guide = "none") +
    labs(title = NULL,
         x = "Validation set",
         y = bquote("Performance ["~R^2~"]")) +
    add_pvalue(df_p_val_lifestyle_genus,
               label = "{p.adj.signif}",
               label.size = 4.5,
               tip.length = 0.01,
               xmin = "xmin",
               xmax = "xmax",
               show.legend = FALSE) +
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          # legend.position=c(0.7, 0.08),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)))
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/downsampling_performance.png",
       width = 14, height = 7, plot = downsampling_performance)


(fig_S3_combined <- cowplot::plot_grid(cowplot::plot_grid(n_unique_resampling, studies_rarefaction, subjects_rarefactin,
                                                          nrow = 1,
                                                          rel_widths = c(1, 1.5, 1.5),
                                                          labels = c("A)", "B)", "C)"),
                                                          label_size = 18,
                                                          label_x = -0.03,
                                                          label_fontface = "plain"),
                                       cowplot::plot_grid(samples_rarefaction,
                                                          downsampling_performance_rare_h,
                                                          nrow = 1,
                                                          labels = c("D)", "E)"),
                                                          label_size = 18,
                                                          label_fontface = "plain",
                                                          label_x = -0.03,
                                                          rel_widths = c(1.5, 2.5)),
                                       cowplot::plot_grid(n_imp_fts_resampling, downsampling_performance,
                                                          nrow = 1,
                                                          labels = c("F)", "G)"),
                                                          label_size = 18,
                                                          label_fontface = "plain",
                                                          label_x = -0.03,
                                                          rel_widths = c(1, 3)),
                                       ncol = 1,
                                       rel_heights = c(1, 1)))

cowplot::save_plot("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/figure_S3_combined_supplement.png", 
                   fig_S3_combined, dpi = 900, bg = "white", base_height = 11.85, base_width = 14)


# Figure S4 ####################################################################

load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/imp_shap_prev.RData")
## A: Importance, Prevalence ####################
imp_prev_ls <- combined_importance_prev_lifestyles_genus %>%
  mutate(model = ifelse(model == "industrialized",
                        yes = "Industrialized",
                        no = "non-Industrialized"),
         Lifestyle = ifelse(lifestyle == "industrialized",
                            yes = "Industrialized",
                            no = "non-Industrialized"),
         display_taxon = if_else(grepl("_| ", display_taxon) |
                                    display_taxon == "Shannon" |
                                    display_taxon == "Observed",
                                  true = glue("{display_taxon}"),
                                  false = glue("<i>{display_taxon}</i>")),
         display_taxon = factor(display_taxon),
         display_taxon = reorder(display_taxon, importance, FUN = max, na.rm = T),
         display_taxon = factor(display_taxon, levels = unique(c("Shannon", "Observed",
                                                                 levels(display_taxon)))))

(imp_prev_combined_plot <- ggplot() +
    geom_col(data = imp_prev_ls %>%
             distinct(display_taxon, prevalence, model),
             aes(y = display_taxon, x = prevalence), alpha = 0.5, fill = "grey", color = "white") +
    geom_boxplot(data = imp_prev_ls,
                 aes(color = Lifestyle, y = display_taxon, x = importance)) +
    scale_color_manual(values = fixed_colors,
                      guide = guide_legend(override.aes = list(
                        shape = 15,
                        size = 5,
                        alpha = 1))) +
    labs(title = NULL,
         x = "Feature Importance",
         y = NULL) +
    facet_wrap(~ model,
                scales = "free_x") +
    theme(axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 16),
      axis.text.y = element_markdown(size = 12),
      # legend.position = c(0.2, 0.2),
      legend.position = "none",
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.grid.major.y = element_line(color = "grey",
        size = 0.3,
        linetype = 2),
      strip.text = element_text(size = 14),
      strip.background = element_blank()))

lifestyle_legend <- ggplot() +
  geom_bar(data = imp_prev_ls,
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
  filter(display_taxon %in% imp_prev_ls$display_taxon) %>%
  left_join(.,imp_prev_ls %>%
              group_by(display_taxon, model) %>%
              summarize(importance = mean(importance),
                        prevalence = mean(prevalence),
                        .groups = "drop"),
            by = c("display_taxon", "model")) %>%
  mutate(display_taxon = factor(display_taxon, levels = levels(imp_prev_ls$display_taxon)),
         R2 = case_when(is.na(R2) ~ 0,
                        is.na(importance) ~ 0,
                        p_adj > 0.05 ~ 0,
                        .default = R2))


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
  filter(display_taxon %in% levels(imp_prev_ls$display_taxon)) %>%
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
  mutate(display_taxon = factor(display_taxon, levels = levels(imp_prev_ls$display_taxon)))


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
                       expand = c(0.001, 0.001)) +
    # scale_color_viridis_c(
    #   option = "H", 
    #   name = "Feature value",
    #   breaks = range(violin_plot_df$abundance, na.rm = TRUE),
    #   labels = c("Low", "High"),
    #   guide = guide_colorbar(direction = "horizontal",
    #                          display = "gradient",
    #                          title.position = "top",
    #                          barwidth = 10)) +
    scale_color_gradientn(name = "Feature value",
                          breaks = range(violin_plot_df$abundance, na.rm = TRUE),
                          labels = c("Low", "High"),
                          guide = guide_colorbar(direction = "horizontal",
                                                 display = "gradient",
                                                 title.position = "top",
                                                 barwidth = 10),
                          colors = wes_palette("Zissou1", type = "continuous"))+
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
(fig_4_combined <- cowplot::plot_grid(cowplot::plot_grid(imp_prev_combined_plot, 
                                                        shap_violin_plot + theme(legend.position = "none"),
                                                        ncol = 2,
                                                        labels = c("A)", "B)"),
                                                        label_size = 18,
                                                        label_fontface = "plain",
                                                        # align = "v",
                                                        # axis = "b",
                                                        rel_widths = c(1.5, 1),
                                                        label_x = c(0, -0.03)),
                                     cowplot::plot_grid(get_legend(lifestyle_legend), 
                                                        get_legend(prevalence_legend),
                                                        get_legend(shap_violin_plot), 
                                                        nrow = 1, rel_widths = c(1, 0.5, 1)),
                                     nrow = 2,
                                     rel_heights = c(1, 0.1)))

cowplot::save_plot("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/figure_S4_combined_supplement.png", 
                   fig_4_combined, dpi = 900,
                   bg = "white", base_height = 20, base_width = 12)
# Fig S5 ############################################
load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/ks_test_res.RData")

## A: Single lifestyle taxa prevalence ##################### 
# Difference in prevalence between lifestyles depending on their importance:
single_taxa <- imp_prev_ls %>%
  group_by(lifestyle, taxon) %>%
  summarise(imp = mean(importance, na.rm = T)) %>%
  filter(is.na(imp)) %>% pull(taxon) %>% unique()

singe_ls_prev_p_val <- imp_prev_ls %>%
  filter(measure_type == "taxon", taxon %in% single_taxa) %>% 
  group_by(taxon, lifestyle) %>%
  summarize(important = any(!is.na(importance)),
            prevalence = unique(prevalence),
            .groups = "drop") %>%
  rstatix::wilcox_test(prevalence ~ important) %>%
  rstatix::add_significance(p.col = "p", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>% 
  rstatix::add_xy_position()

ls_importance <- imp_prev_ls %>%
  filter(measure_type == "taxon", taxon %in% single_taxa) %>% 
  group_by(taxon, lifestyle) %>%
  summarize(important = any(!is.na(importance)),
            prevalence = unique(prevalence)) %>%
  filter(important) %>% select(taxon, lifestyle)

(single_ls_prev <- imp_prev_ls %>%
    filter(measure_type == "taxon", taxon %in% single_taxa) %>% 
    group_by(taxon, lifestyle) %>%
    summarize(important = any(!is.na(importance)),
              prevalence = unique(prevalence)) %>%
    ggplot(.,aes(x = important, y = prevalence)) +
    stat_eye() +
    # geom_line(aes(group = taxon)) +
    xlab("Important feature") +
    ylab("Prevalence") +
    add_pvalue(singe_ls_prev_p_val,
               label = "{p.signif}",
               label.size = 4.5,
               tip.length = 0.01,
               xmin = "xmin",
               xmax = "xmax",
               show.legend = FALSE) +
    # facet_wrap(~lifestyle.y) +
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          plot.title = element_text(size=18),
          # legend.position = c(0.8, 0.2),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_blank(),
          strip.background = element_blank()))

## B: Prevalence vs importance ###################
(prev_vs_imp <- imp_prev_ls %>% 
    filter(measure_type == "taxon") %>% 
    group_by(taxon, Lifestyle) %>%
    summarize(importance = mean(importance, na.rm = T),
              prevalence = unique(prevalence)) %>%
    ggplot(.,aes(y = importance, x = prevalence, color = Lifestyle, fill = Lifestyle)) +
    geom_point() +
    geom_smooth(method = "lm") +
    ggpmisc::stat_poly_eq(ggpmisc::use_label("R2", "P"),
                          vstep = 0.1,
                          size = 4) +
    scale_color_manual(values = fixed_colors, 
                       guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1),
                                            direction = "vertical",
                                            title.position = "top")) +
    scale_fill_manual(values = fixed_colors,
                      guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
   labs(title = NULL,
        x = "Prevalence",
        y = "Mean importance") +
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.position = c(0.8, 0.9),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_blank(),
          strip.background = element_blank()))  


## C: KS-test SHAP distributions #####################################
(ks_plot <- ks_res_all %>%
   mutate(Group = case_when(comparison == "between_models" ~ "Between models, within taxa",
                            comparison == "between_taxa" & model == "industrialized" ~ "Between taxa, industrialized",
                            comparison == "between_taxa" & model == "non_industrialized" ~ "Between taxa, non-industrialized")) %>%
   group_by(Group) %>% #summarize(mean_D = mean(D_ks)) 
   ggplot(., aes(x = D_ks, color = Group, fill = Group)) +
   # geom_density(alpha = 0.5) +
   ggdist::stat_slabinterval(.width = 0, alpha = 0, shape = 21) +
   ggdist::stat_slabinterval(.width = 0, slab_alpha = .5, shape = 15, show.legend = F, size = 15) +
   scale_color_manual(values = c("Between models, within taxa" = "#DB3A07FF",
                                 "Between taxa, industrialized" = "#737125",
                                 "Between taxa, non-industrialized" = "#1A97C8")) +
   scale_fill_manual(values = c("Between models, within taxa" = "#DB3A07FF",
                                 "Between taxa, industrialized" = "#737125",
                                 "Between taxa, non-industrialized" = "#1A97C8"), 
                      guide = guide_legend(override.aes = list(size = 5, alpha = 1),
                                           direction = "vertical",
                                           title.position = "top")) +
   labs(title = NULL,
        x = "Kolmogorov–Smirnov D",
        y = "Density") +
   theme(axis.text.x = element_text(size = 12),
         axis.title.x = element_text(size = 16),
         axis.text.y = element_text(size = 12),
         axis.title.y = element_text(size = 16),
         legend.title = element_text(size = 15),
         legend.text = element_text(size = 15),
         # legend.position = "bottom",
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         panel.grid.major.y = element_blank()))
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/ks_figure.png", width = 11, height = 7,
       plot = ks_plot)



(fig_5_combined <- cowplot::plot_grid(single_ls_prev, 
                                      prev_vs_imp, 
                                      get_legend(ks_plot), 
                                      ks_plot + theme(legend.position = "none"), 
                                      labels = c("A)", "B)", "C)", NA),
                                      label_size = 18,
                                      label_fontface = "plain",
                                      rel_widths = c(1, 2),
                                      # align = "v",
                                      # axis = "b",
                                      label_x = -0.01))
                                                       

cowplot::save_plot("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/figure_S5_combined_supplement.png", 
                   fig_5_combined, dpi = 900,
                   bg = "white", base_height = 10, base_width = 14)

# Figure S6: #######################################################################
load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/phyloseq_supplement.RData")
load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/prevalences.RData")

selected_taxa <- c("Staphylococcus", "Faecalibacterium", "Prevotella",
                   "Lactobacillus", "Bifidobacterium")
## A: Prevalences rarefied ##############################
abundance_data_rarefied <- ps_object_genus_raref_comp %>% 
  psmelt() %>%
  select(lifestyle, Abundance, OTU, age) %>%
  filter(OTU %in% selected_taxa) %>%
  mutate(Lifestyle = ifelse(lifestyle == "industrialized", yes = "Industrialized", no = "non-Industrialized"),
         month = ceiling(age / 30)) %>%
  left_join(., short_taxa_names, by = c("OTU" = "taxon")) %>% #pull(display_taxon) %>% unique %>% sort
  filter(display_taxon %in% combined_importance_prev_lifestyles_genus$display_taxon)

(prevalence_plot_raref <- prev_per_month_raref %>%
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
    scale_fill_manual(values = fixed_colors) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    labs(title = NULL,
         x = "Age [Months]",
         y = "Prevalence") +
    facet_wrap(~display_taxon) +
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_blank(),
          strip.text = element_text(size = 14, face = "italic"),
          strip.background = element_blank()))


## B: Abundance over time ######################
abundance_data <- ps_object_genus_comp %>% 
  psmelt() %>%
  select(lifestyle, Abundance, OTU, age) %>%
  filter(OTU %in% selected_taxa) %>%
  mutate(Lifestyle = ifelse(lifestyle == "industrialized", yes = "Industrialized", no = "non-Industrialized"),
         month = ceiling(age / 30)) %>%
  left_join(., short_taxa_names, by = c("OTU" = "taxon")) %>% #pull(display_taxon) %>% unique %>% sort
  filter(display_taxon %in% combined_importance_prev_lifestyles_genus$display_taxon)
  
(abund_plot_suppl <-  abundance_data %>%
    ggplot(., aes(x = age, y = Abundance, color = Lifestyle, fill = Lifestyle)) +
    geom_point(alpha = 0.1, size = 0.3) +
    geom_smooth(method = "loess", alpha = 0) +
    scale_color_manual(values = fixed_colors, 
                       guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1),
                                            direction = "vertical",
                                            # position = "bottom",
                                            title.position = "top")) +
    scale_fill_manual(values = fixed_colors) +
    labs(title = NULL,
         x = "Age [Days]",
         y = "Relative abundance") +
    facet_wrap(~display_taxon) +
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          plot.title = element_text(size=18),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_blank(),
          strip.text = element_text(size = 14, face = "italic"),
          strip.background = element_blank()))


## C: Abundance over time rarefied ######################
(abund_plot_suppl_raref <-  abundance_data_rarefied %>%
    ggplot(., aes(x = age, y = Abundance, color = Lifestyle, fill = Lifestyle)) +
    geom_point(alpha = 0.1, size = 0.3) +
    geom_smooth(method = "loess", alpha = 0) +
    scale_color_manual(values = fixed_colors, 
                       guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1),
                                            direction = "vertical",
                                            # position = "bottom",
                                            title.position = "top")) +
    scale_fill_manual(values = fixed_colors) +
    facet_wrap(~display_taxon) +
   labs(title = NULL,
        x = "Age [Days]",
        y = "Relative abundance")+
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.position = c(0.8, 0.2),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_blank(),
          strip.text = element_text(size = 14, face = "italic"),
          strip.background = element_blank()))



fig_S6_combined <- cowplot::plot_grid(prevalence_plot_raref,
                                      abund_plot_suppl,
                                      abund_plot_suppl_raref,
                                      # rel_widths = c(1.5, 1),
                                      nrow = 3,
                                      # nrow = 3, rel_heights = c(1, 0.2, 2)
                                      labels = c("A)", "B)", "C)"),
                                      label_size = 18,
                                      label_fontface = "plain")

cowplot::save_plot("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/figure_S6_combined_supplement.png", 
                   fig_S6_combined, dpi = 900,
                   bg = "white", base_height = 21.5, base_width = 14)


# Figure S7 ################################
## A: Malnourished PCA ################
fixed_colors <- c(Industrialized = "#737125", `non-Industrialized` = "#1A97C8", Combined = "#2A2359")
fixed_colors_sick <- fixed_colors
names(fixed_colors_sick) <- paste0(names(fixed_colors_sick), " healthy")
# fixed_colors_sick <- c(fixed_colors_sick, `non-Industrialized SAM` = "#25ECA7FF", `Industrialized preterm` = "#4147ADFF")
fixed_colors_sick <- c(fixed_colors_sick, `non-Industrialized SAM` = "#812B8C", `Industrialized preterm` = "#D9731A")

load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/malnurished_pcoa.RData")
explained_var_mal <- ((pcoa_aitch_genus_mal$sdev^2/sum(pcoa_aitch_genus_mal$sdev^2)) * 100) %>%
  round(.,digits = 1)
mdat_malnurished_pca <- mdat_malnurished_all[rownames(pcoa_aitch_genus_mal$x),] %>%
  bind_cols(pcoa_aitch_genus_mal$x[,c(1, 2)])
malnurished_permanova_label <- c("R[Lifestyle]^2 == 0.09*';'~p < 0.001", "~R[Health]^2 == 0.03*';'~p < 0.001")
(malnourished_pca <- mdat_malnurished_pca %>%
  mutate(health_lifestyle = paste(lifestyle, health),
         health_lifestyle = gsub("industrialized", "Industrialized",health_lifestyle) %>%
           gsub("non_", "non-",.)) %>%
  ggplot(., aes(x = PC1, y = PC2, color = health_lifestyle)) +
  geom_point(alpha = 0.8, size = 0.3) +
    annotate("text",
             x = min(mdat_malnurished_pca$PC1) * 1.2,  # Adjust position
             y = c(max(mdat_malnurished_pca$PC2), max(mdat_malnurished_pca$PC2) - 3),
             hjust = 0,
             size = 4,
             label = malnurished_permanova_label, parse = T) +
  scale_size_manual(values = c(0.3, 5)) +
    scale_color_manual(values = fixed_colors_sick,
                       name = "Lifestyle and health",
                       guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
    scale_fill_manual(values = fixed_colors_sick,
                      name = "Lifestyle and health") +
    coord_equal() +
    labs(title = NULL,
         x = paste0("PC 1 (", explained_var_mal[1], "%)"),
         y = paste0("PC 2 (", explained_var_mal[2], "%)")) +
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.text = element_text(size = 13, margin = margin(r = 13)),
          legend.title = element_text(size = 13)))
## B: Preterm PCA ############################
load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/preterm_pcoa.RData")

explained_var_pre <- ((pcoa_aitch_genus_preterm$sdev^2/sum(pcoa_aitch_genus_preterm$sdev^2)) * 100) %>%
  round(.,digits = 1)

mdat_preterm_pca <- mdat_preterm_all[rownames(pcoa_aitch_genus_preterm$x),] %>%
  bind_cols(pcoa_aitch_genus_preterm$x[,c(1, 2)]) %>%
  filter(!is.na(health))
ggplot(mdat_preterm_pca %>% group_by(study, lifestyle, health) %>% summarize(sample_sum = median(sample_sum, na.rm = T)),
       aes( x = health, color = lifestyle, y = sample_sum)) +
  geom_violin() +
  scale_y_log10()
  

preterm_permanova_label <- c("R[Lifestyle]^2 == 0.03*';'~p < 0.001", "~R[Health]^2 == 0.09*';'~p < 0.001")
(preterm_pca <- mdat_preterm_pca %>%
    mutate(health_lifestyle = paste(lifestyle, health),
           health_lifestyle = gsub("industrialized", "Industrialized",health_lifestyle) %>%
             gsub("non_", "non-",.)) %>%
    ggplot(., aes(x = PC1, y = PC2, color = health_lifestyle)) +
  geom_point(alpha = 0.8, size = 0.3) +
    annotate("text",
             x = min(mdat_preterm_pca$PC1) * -0.7,  # Adjust position
             y = c(max(mdat_preterm_pca$PC2), max(mdat_preterm_pca$PC2) - 2.4),
             hjust = 0,
             size = 4,
             label = preterm_permanova_label, parse = T) +
    
  scale_size_manual(values = c(0.3, 5)) +
    scale_color_manual(values = fixed_colors_sick,
                       name = "Lifestyle and health",
                       guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
    scale_fill_manual(values = fixed_colors_sick,
                      name = "Lifestyle and health") +
    coord_equal() +
    labs(title = NULL,
         x = paste0("PC 1 (", explained_var_pre[1], "%)"),
         y = paste0("PC 2 (", explained_var_pre[2], "%)")) +
    theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size = 13, margin = margin(r = 13)),
        legend.title = element_text(size = 13)))


## C: Malnourished SHAP complete ############################
load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/sick_shap_data.RData")
pvals_svz_mal_age <- shap_mdat_mal %>%
  mutate(age_bin = floor(ceiling(age) / 30)) %>%
  group_by(display_taxon, age_bin) %>%
  rstatix::wilcox_test(shap_value ~ health) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>%
  filter(p.adj.signif != "ns")

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
# taxa reducing age:
mean_shap_diff_mal %>%
  filter(display_taxon %in% unique(pvals_svz_mal_age$display_taxon)) %>%
  filter(!is.na(significance), diff_shap < 0) %>% pull(display_taxon) %>% unique %>% length

(heatmap_mal_complete <- mean_shap_diff_mal %>%
    filter(display_taxon %in% unique(pvals_svz_mal_age$display_taxon)) %>%
    mutate(display_taxon = if_else(grepl("_| ", display_taxon) |
                                     display_taxon == "Shannon" |
                                     display_taxon == "Observed",
                                   true = glue("{display_taxon}"),
                                   false = glue("<i>{display_taxon}</i>")) %>%
             gsub("_", ".",.),
           display_taxon = factor(display_taxon, levels = unique(c("Observed", "Shannon", display_taxon)))) %>%
    ggplot(., aes(x = age_bin, y = display_taxon, fill = diff_shap)) +
    geom_tile(color = "white") +                   # White borders between tiles
    geom_text(aes(label = significance), color = "black", size = 2) +  # significance stars
    scale_x_continuous(breaks = c(6, 12, 18)) +
    scale_fill_gradient2(name = "\U0394 mean SHAP",
                         mid = "white",
                         low = c("#4662D7FF", "#36AAF9FF", "#1AE4B6FF"),
                         high = c("#FABA39FF", "#F66B19FF", "#CB2A04FF")) +
    labs(title = "SAM",
         x = "Chronological age [Months]",
         y = NULL) +
    theme(axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.text.y = element_markdown(size = 12),
          plot.title = element_text(size=18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          strip.text = element_text(size = 14)))

## D: Preterm SHAP complete ################################
pvals_svz_preterm_age <- shap_mdat_preterm %>%
  mutate(age_bin = floor(age / 7),
         age_bin = ifelse(age_bin == 12, yes = 11, no = age_bin)) %>% # last bin contains too few values for tests
  group_by(display_taxon, age_bin) %>%
  rstatix::wilcox_test(shap_value ~ health) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>%
  filter(p.adj.signif != "ns")

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
mean_shap_diff_preterm %>%
  filter(display_taxon %in% unique(pvals_svz_preterm_age$display_taxon)) %>%
  filter(!is.na(significance), diff_shap < 0) %>% pull(display_taxon) %>% unique %>% length
mean_shap_diff_preterm %>%
  filter(display_taxon %in% unique(pvals_svz_preterm_age$display_taxon)) %>%
  filter(!is.na(significance), diff_shap > 0) %>% pull(display_taxon) %>% unique %>% length


(heatmap_preterm_complete <- mean_shap_diff_preterm %>%
    filter(display_taxon %in% unique(pvals_svz_preterm_age$display_taxon)) %>%
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
    geom_text(aes(label = significance), color = "black", size = 2, angle = 90) +  # significance stars
    scale_x_continuous(breaks = c(0, 4, 8, 11)) +
    scale_fill_gradient2(name = "\U0394 mean SHAP",
                         mid = "white",
                         low = c("#4662D7FF", "#36AAF9FF", "#1AE4B6FF"),
                         high = c("#FABA39FF", "#F66B19FF", "#CB2A04FF")) +
    labs(title = "Preterm",
         x = "Chronological age [Weeks]",
         y = NULL) +
    theme(axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.text.y = element_markdown(size = 12),
          axis.title.y = element_blank(),
          plot.title = element_text(size=18),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
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
cowplot::set_null_device("agg")

fig_S7_combined <- cowplot::plot_grid(cowplot::plot_grid(malnourished_pca, preterm_pca,
                                                         rel_widths = c(1, 1.25),
                                                         labels = c("A)", "B)"),
                                                         label_size = 18,
                                                         label_fontface = "plain"),
                   get_legend(legend_plot),
                   cowplot::plot_grid(heatmap_mal_complete, heatmap_preterm_complete,
                                      rel_widths = c(1.5, 1),
                                      labels = c("C)", "D)"),
                                      label_size = 18,
                                      label_fontface = "plain"),
                   nrow = 3, rel_heights = c(1, 0.2, 2))

cowplot::save_plot("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/figure_S7_combined_supplement.png", 
                   fig_S7_combined, dpi = 900,
                   bg = "white", base_height = 20.5, base_width = 14)



# Fig. 8 ##############################
load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/sick_age_predictions.RData")
# fixed_colors_sick <- fixed_colors
# names(fixed_colors_sick) <- paste0(names(fixed_colors_sick), " healthy")
# fixed_colors_sick <- c(fixed_colors_sick, `non-Industrialized SAM` = "#25ECA7FF", `Industrialized preterm` = "#4147ADFF")
## A: Malnourished only non-industrialized #####################################
all_preds_combined_only_bangladesh <- all_preds_combined %>% 
  filter(study %in% c("subramanian_2014", "gehrig_2019"))

(mal_pred_age <-  all_preds_combined_only_bangladesh %>%
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
         x = "Chronological Age [days]",
         y = "Microbial Age [days]") +
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

means_maz_mal_only_bangladesh <- all_preds_combined_only_bangladesh %>% 
  group_by(training_set, lifestyle_health) %>%
  dplyr::summarize(mean = mean(MAZ))

pvals_maz_mal_only_bangladesh <- all_preds_combined_only_bangladesh %>%
  group_by(training_set) %>%
  rstatix::wilcox_test(MAZ ~ lifestyle_health) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj", cutpoints = c(0, 1e-03, 0.01, 0.05, 0.1, 1)) %>%
  rstatix::add_xy_position(x = "training_set", dodge = 1) %>%
  left_join(.,means_maz_mal_only_bangladesh, by = c("training_set" = "training_set",
                                    "group1" = "lifestyle_health")) %>%
  left_join(.,means_maz_mal_only_bangladesh, by = c("training_set" = "training_set",
                                    "group2" = "lifestyle_health")) %>%
  mutate(MAZ_diff = abs(mean.x - mean.y)) %>%
  select(MAZ_diff, everything())

## B: Malnurished MAZ only non-industrialized ##################################
(mal_maz_diff <- all_preds_combined_only_bangladesh %>%
    mutate(lifestyle_health = gsub("_", " ", lifestyle_health) %>%
             gsub("industrialized", "Industrialized", .) %>%
             gsub("non ", "non-", .),
           training_set = case_when(training_set == "combined" ~ "Combined",
                                    training_set == "non_industrialized" ~ "non-Industrialized",
                                    training_set == "industrialized" ~ "Industrialized")) %>%
    ggplot(., aes(x = training_set, y = MAZ)) +
    # geom_boxplot(aes(fill = lifestyle_health)) +
   stat_eye(aes(fill = lifestyle_health), position = position_dodge(1), scale = 0.9,
                .width = c(0, 0.5, 0.95), adjust = 1, shape = 23, point_size = 3,
                side = "both") +
   
    scale_fill_manual(values = fixed_colors_sick,
                      name = "Lifestyle and health",
                      guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
    add_pvalue(pvals_maz_mal_only_bangladesh,
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
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)))

## C: Malnourished Features over time ################
(mal_taxa_trends <- shap_mdat_mal %>% filter(lifestyle == "non_industrialized",
                             display_taxon %in% c("Faecalibacterium", "Observed")) %>%
   mutate(lifestyle_health = gsub("_", " ", lifestyle_health) %>%
            gsub("industrialized", "Industrialized", .) %>%
            gsub("non ", "non-", .),
          display_taxon = ifelse(display_taxon == "Observed", yes = display_taxon,
                                  no = glue("<i>{display_taxon}</i>"))) %>% 
  ggplot(., aes(x = age, y = ab_value, color = lifestyle_health, fill = lifestyle_health)) +
  geom_point(size = 0.2, alpha = 0.3) +
   geom_smooth(method = "loess", alpha = 0) +
   scale_color_manual(values = fixed_colors_sick,
                      name = "Lifestyle and health",
                      guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
   scale_fill_manual(values = fixed_colors_sick,
                     name = "Lifestyle and health") +
   facet_wrap(~display_taxon, scales = "free_y") +
   labs(title = NULL,
        y = "Feature value",
        x = "Age [Days]") +
   theme(axis.text.x = element_text(size = 14),
         axis.title.x = element_text(size = 16),
         axis.text.y = element_text(size = 14),
         axis.title.y = element_text(size = 16),
         legend.position = "none",
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         panel.grid.major.y = element_blank(),
         strip.text = element_markdown(size = 14),
         strip.background = element_blank()))



## D: Preterm features over time ########################
(preterm_taxa_trends <- shap_mdat_preterm %>% filter(lifestyle == "industrialized",
                             display_taxon %in% c("Staphylococcus", "Observed")) %>%
   mutate(lifestyle_health = gsub("_", " ", lifestyle_health) %>%
            gsub("industrialized", "Industrialized", .) %>%
            gsub("non ", "non-", .),
          display_taxon = ifelse(display_taxon == "Observed", yes = display_taxon,
                                 no = glue("<i>{display_taxon}</i>"))) %>% 
  ggplot(., aes(x = age, y = ab_value, color = lifestyle_health, fill = lifestyle_health)) +
  geom_point(size = 0.2, alpha = 0.3) +
   geom_smooth(method = "loess", alpha = 0) +
   scale_color_manual(values = fixed_colors_sick,
                      name = "Lifestyle and health",
                      guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
   scale_fill_manual(values = fixed_colors_sick,
                     name = "Lifestyle and health") +
   labs(title = NULL,
        x = "Age [Days]",
        y = "Feature value") +
   facet_wrap(~display_taxon, scales = "free_y") +
   theme(axis.text.x = element_text(size = 14),
         axis.title.x = element_text(size = 16),
         axis.text.y = element_text(size = 14),
         axis.title.y = element_text(size = 16),
         legend.position = "none",
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         panel.grid.major.y = element_blank(),
         strip.text = element_markdown(size = 14),
         strip.background = element_blank()))

(legend_plot_2  <- data.frame(lifestyle_health = names(fixed_colors_sick[-3])) %>%
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


(fig_S8_combined <- cowplot::plot_grid(mal_pred_age + theme(margin = margin(t = 1)),
                                       mal_maz_diff,
                                       mal_taxa_trends, 
                                       preterm_taxa_trends,
                                       get_legend(legend_plot_2),
                                       ncol = 1,
                                      # rel_widths = c(1.5, 1),
                                      labels = c("A)", "B)", "C)", "D)", NA),
                                      label_size = 18,
                                      label_fontface = "plain",
                                      rel_heights = c(1, 1, 1, 1, 0.2)))

cowplot::save_plot("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/figure_S8_combined_supplement.png", 
                   fig_S8_combined, dpi = 900,
                   bg = "white", base_height = 14, base_width = 14)


# Supplementary Tables #################################
# supplementary table 1:
meta_df %>%
  group_by(study, lifestyle) %>%
  summarize(country = paste(unique(country), collapse = ", "),
            n_samples = n(),
            n_individuals = length(unique(subject_ID)),
            cutoff = unique(cutoff)) %>% 
  write.csv("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/supplementary_tables/supplementary_table_S1.csv")

# supplementary table 2 ####################################
# supplementary table with corrs/p-vals from shap, prevalence and mean importance:
shap_cors_filt_xy_plot <- shap_cors_all %>%
  select(display_taxon, display_taxon_filter, model, R2) %>%
  pivot_wider(names_from = "model", values_from = R2) %>%
  filter(!(Industrialized == 0 & `non-Industrialized` == 0), 
         !display_taxon_filter %in% c("Shannon", "Observed")) %>%
  mutate(Classification = case_when(`non-Industrialized` > 0 & Industrialized > 0 ~ "Late colonizer",
                                    `non-Industrialized` < 0 & Industrialized < 0 ~ "Early colonizer",
                                    Industrialized == 0 | `non-Industrialized` == 0 ~ "Single lifestyle",
                                    .default = "Mixed"),
         display_taxon_filter = droplevels(display_taxon_filter))

shap_cors_all %>% mutate(taxon = display_taxon_filter,
                         mean_importance = importance,
                         .keep = "unused") %>%
  select(-display_taxon, -cor_test) %>%
  pivot_wider(names_from = "model", values_from = c("R2", "p_value", "p_adj", "mean_importance", "prevalence")) %>%
  left_join(shap_cors_filt_xy_plot %>% mutate(taxon = display_taxon_filter,
                                              .keep = "unused") %>%
              select(taxon, Classification)) %>% 
  write.csv("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/supplementary_tables/supplementary_table_S2.csv")



# Prevalences ########################
load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/prevalences.RData")
short_taxa_names <- short_taxa_names %>%
  mutate(taxon = taxon %>%
           gsub("X[.]", "[",.) %>%
           gsub("[.][.]", "] ",.) %>%
           gsub("UCG[.]", "UCG-",.) %>%
           gsub("[.]", " ",.))

(prevalence_plot_suppl <- prev_per_month %>%
    mutate(Lifestyle = ifelse(lifestyle == "industrialized", yes = "Industrialized", no = "non-Industrialized")) %>%
    left_join(., short_taxa_names, by = c("OTU" = "taxon")) %>% #pull(display_taxon) %>% unique %>% sort
    filter(display_taxon %in% combined_importance_prev_lifestyles_genus$display_taxon) %>%
    ggplot(., aes(x = month, y = prev, color = Lifestyle, fill = Lifestyle)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", alpha = 0) +
    scale_color_manual(values = fixed_colors, 
                       guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1),
                                            direction = "vertical",
                                            position = "bottom",
                                            title.position = "left")) +
    scale_fill_manual(values = fixed_colors) +
    facet_wrap(~display_taxon) +
    xlab("Age [Months]") +
    ylab("Prevalence")+
    ggtitle("B)")+
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          plot.title = element_text(size=18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_blank(),
          strip.text = element_text(size = 14, face = "italic"),
          strip.background = element_blank()))

ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/prevalences.png", 
       width = 10, plot = prevalence_plot_suppl)


# Figure SX ####################################################################
load(file = "/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/rarefaction_studies.RData")

rf_rf_df_genus %>%
  group_by(study, lifestyle) %>%
  rstatix::cor_test(., vars = c("n_studies", "R2"), method = "spearman")


(rarefaction_studies <- ggplot(rf_rf_df_genus %>% mutate(Lifestyle = ifelse(lifestyle == "industrialized",
                                                                            yes = "Industrialized",
                                                                            no = "non-Industrialized")),
                               aes(x = n_studies, y = R2, color = Lifestyle, fill = Lifestyle)) +
    geom_point(alpha = 0.8, size = 0.5) +
    # geom_smooth(method = "lm", aes(fill = Lifestyle)) +
    # stat_cor(method = "spearman",
    #          cor.coef.name = "spearman",
    #          # size = 5,
    #          show.legend = F) +
    geom_hline(data = plyr::ddply(rf_rf_df_genus, c("study"),
                                  summarize, final_model = R2[which.max(n_studies)]),
               aes(yintercept=final_model), size = 0.4, linetype = "dashed") +
    scale_color_manual(values = fixed_colors, 
                       guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
    scale_fill_manual(values = fixed_colors) +
    # scale_color_manual(values = fixed_colors,
    #                    guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
    facet_wrap(~ study) +
    ylab(bquote("Performance ["~R^2~"]"))+
    xlab("Number of studies") +
    ylim(0, 1.1) +
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size=18),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          strip.text = element_text(size = 15),
          strip.background = element_blank()))
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/rarefaction_studies.png",
       width = 14, height = 7, plot = rarefaction_studies)


## X: importance vs prevalence ###########################
# For downsampling stuff, but keep it out for now
load("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/imp_prevs.RData")

imp_prev_original <- combined_importance_prev_lifestyles_genus %>%
  group_by(taxon, lifestyle) %>%
  summarize(Prevalence = mean(prevalence),
            Importance = mean(importance)) %>%
  mutate(Lifestyle = ifelse(lifestyle == "industrialized", yes = "Industrialized", no = "non-Industrialized"),
         Prevalence = Prevalence / 100,
         figure = "Complete")

imp_prev_downsampling <- imp_prevs %>%
  filter(counts > 1) %>%
  group_by(taxon, lifestyle) %>%
  summarize(Prevalence = mean(prev),
            Importance = mean(importance)) %>%
  mutate(Lifestyle = ifelse(lifestyle == "industrialized", yes = "Industrialized", no = "non-Industrialized"),
         figure = "Downsampling")

imp_prev_original %>% filter(Prevalence > 0.75, Importance > 50) %>% pull(taxon) %>% unique
imp_prev_downsampling %>% filter(Prevalence > 0.75, Importance > 50) %>% pull(taxon) %>% unique

(prev_imp_combined <- bind_rows(imp_prev_original, imp_prev_downsampling) %>%
    filter(figure == "Complete") %>% #, Prevalence > 0.15) %>%
    ggplot(., aes(x = Prevalence, y = Importance, color = Lifestyle, fill = Lifestyle)) +
    geom_point(alpha = 0.8, size = 0.8) +
    geom_smooth(method = "lm", aes(fill = Lifestyle)) +
    stat_cor(method = "spearman",
             cor.coef.name = "spearman",
             size = 5,
             show.legend = F) +
    scale_color_manual(values = fixed_colors, 
                       guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
    scale_fill_manual(values = fixed_colors) +
    xlim(0,1)+
    # ggtitle("B") +
    # facet_wrap(~figure) +
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          panel.grid.majors = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size=18),
          strip.text = element_text(size = 15, hjust = 0),
          strip.background = element_blank(),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)))
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/supplement/prev_imp_combined.png",
       width = 14, height = 7, plot = prev_imp_combined)
