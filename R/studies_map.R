# load library tidyverse
library(tidyverse)
# library(maptools)
library(maps)
library(ggmap)
library(purrr)
library(scatterpie)
library(viridis)
maindir <- "/fast/AG_Forslund/rob/studies/"


# register_google(key = "mQkzTpiaLYjPqXQBotesgif3EfGL2dbrNVOrogg")
# create data for world coordinates using 
# map_data() function
world_coordinates <- map_data("world")

# amplicon_location <- c("Sweden", "Bangladesh", "Netherlands", "Russia", "Finnland", "Estonia",
#              "USA", "Luxembourg", "Ireland", "Bolivia", "Brazil", "India", "Peru",
#              "South Africa", "Malawi", "Denmark", "Norway", "El Salvador", "Gambia",
#              "Malawi", "Puerto Rico", "Chile", "Spain")
# shotgun_location <- c("Italy", "Sweden", "Finnland", "USA", "Russia", "Estonia", 
#                       "Netherlands", "Luxembourg", "UK", "New Zealand", "Bangladesh", 
#                       "El Salvador", "Peru", "Tanzania", "Malawi")
# world_coordinates <- mutate(world_coordinates, fill = case_when(
#   (region %in% amplicon_location) & !(region %in% shotgun_location) ~ "red", 
#   (region %in% shotgun_location) & !(region %in% amplicon_location) ~ "blue",
#   region %in% c(amplicon_location, shotgun_location) ~ "green", 
#   TRUE ~ "white"))
# 
# # create world map using ggplot() function
# ggplot(world_coordinates, aes(long, lat, fill = fill, group = group)) +
#   
#   # geom_map() function takes world coordinates 
#   # as input to plot world map
#   geom_map(
#     data = world_coordinates, map = world_coordinates,
#     aes(long, lat, map_id = region)
#   ) +
#   geom_polygon(colour="gray") +
#   scale_fill_identity(breaks = c("blue", "red", "green"),
#                       labels = c("Shotgun", "16S", "both"),
#                       guide = "legend")

# load 16S data:
setwd("/fast/AG_Forslund/rob/studies/16S/")
source("/fast/AG_Forslund/rob/mm_index/R_scripts/setlists.R")
study_names <- set_list$all#, "song_2021", "iszatt_2019", "goffau_2022")

read_metad <- function(s_name) {
  print(s_name)
  f_name <- list.files(s_name, pattern = "metadata_.*\\_healthy.csv", full.names = T)
  df <- read.table(f_name, header = T, sep = ",")
  df[,c("geographic_location_.latitude.", "geographic_location_.longitude.", "study", "lifestyle", "country")]
}

df_16S <- lapply(study_names, read_metad) %>%
  purrr::reduce(rbind) %>%
  mutate(latitude = round(geographic_location_.latitude., 2),
         longitude = round(geographic_location_.longitude., 2),
         .keep = "unused") %>%
  mutate(lat_lon = paste(latitude, longitude, study, lifestyle, country))

data_loc_16S <- table(df_16S$lat_lon) %>% 
  data.frame %>%
  separate(., col = Var1, into = c("latitude", "longitude", "study", "lifestyle", "country"), sep = " ", convert = T) %>%
  mutate(sample_count = Freq, Type = "16S", .keep = "unused")

# get shotgun data:
setwd("/fast/AG_Forslund/rob/studies/shotgun/")
study_names_shotgun <- list.dirs(recursive = FALSE, full.names = F)
study_names_shotgun <- study_names_shotgun[study_names_shotgun != "failed_ngless"]
study_names_shotgun <- study_names_shotgun[study_names_shotgun != "failed_ngless_2"]
study_names_shotgun <- study_names_shotgun[study_names_shotgun != "test_failed"]
study_names_shotgun <- study_names_shotgun[study_names_shotgun != "vatanen_2022_2"]
# study_names_shotgun <- study_names_shotgun[study_names_shotgun != "bargheet_2023"]
study_names_shotgun <- study_names_shotgun[study_names_shotgun != "gasparrini_2019"]
study_names_shotgun <- study_names_shotgun[study_names_shotgun != "hoskinson_2023"]
study_names_shotgun <- study_names_shotgun[study_names_shotgun != "metaphlan_combined"]
study_names_shotgun <- study_names_shotgun[study_names_shotgun != "bins_for_michal"]
study_names_shotgun <- study_names_shotgun[study_names_shotgun != "selma_royo_2024"]
study_names_shotgun <- study_names_shotgun[study_names_shotgun != "tett_2019"]
study_names_shotgun <- study_names_shotgun[study_names_shotgun != "bottino_2025"]
study_names_shotgun <- study_names_shotgun[study_names_shotgun != "portlock_2024"]
study_names_shotgun <- study_names_shotgun[study_names_shotgun != "deng_2025"]
# study_names_shotgun <- study_names_shotgun[study_names_shotgun != "parnanen_2018"]



df_shotgun <- lapply(study_names_shotgun, read_metad) %>%
  purrr::reduce(rbind) %>%
  mutate(latitude = round(geographic_location_.latitude., 2),
         longitude = round(geographic_location_.longitude., 2),
         .keep = "unused") %>%
  mutate(lat_lon = paste(latitude, longitude, study, lifestyle, country))

data_loc_shotgun <- table(df_shotgun$lat_lon) %>% 
  data.frame %>%
  separate(., col = Var1, into = c("latitude", "longitude", "study", "lifestyle", "country"), sep = " ", convert = T) %>%
  mutate(sample_count = Freq, Type = "shotgun", .keep = "unused")

data_loc <- rbind(data_loc_16S, data_loc_shotgun) %>% 
  # mutate(value = Freq) %>%
  arrange(-sample_count) %>%
  filter(Type == "16S")#%>%
  # pivot_wider(names_from = sequencing, values_from = Freq, values_fill = 0) %>%
  # mutate(total = `16S` + shotgun)
# filter some regions out
world_coordinates <- world_coordinates %>% filter(!(region %in% c("Antarctica", "South Sandwich Islands", "Fiji")))

# add gdp data
gdp_data <- read.table("/fast/AG_Forslund/rob/mm_index/R_scripts/fc190663-b0af-4e9a-a3f7-f2e235923a8f_Data.csv",
                       sep = ",", header = T) %>%
  select(Country.Name, X2019..YR2019.) %>%
  mutate(GDP = as.numeric(X2019..YR2019.)) %>%
  mutate(Country.Name = toupper(Country.Name) %>%
           gsub(" ", "_",.) %>%
           gsub("FINLAND", "FINNLAND",.))


hdi_data <- read.table("/fast/AG_Forslund/rob/mm_index/R_scripts/hdi_data.txt",
                       sep = "\t", col.names = c("region", "hdi")) %>%
  mutate(region = toupper(region) %>%
           gsub(" ", "_",.) %>%
           gsub("FINLAND", "FINNLAND",.))


# 
# world_coordinates <- left_join(world_coordinates, gdp_data, by = c("region" = "Country.Name"))
# world_coordinates <- left_join(world_coordinates, hdi_data, by = "region")
# world_coordinates %>% filter(is.na(hdi)) %>% select(region) %>% unlist %>% unique() %>% sort
  
all_studies <- set_list$all
studies_here <- unique(data_loc$study) %in% all_studies
ggplot() +
  geom_map(
    data = world_coordinates, map = world_coordinates,
    aes(long, lat, map_id = region),
    color = "grey", fill = "white", size = 0.15
  )+
  geom_point(data = data_loc_shotgun,
             aes(x = longitude, y = latitude, size = sample_count, color = study, fill = study),
             alpha = 0.7,
             position = position_jitter(width=0.7),
             shape=21) + #, color = "#9266CC") +

  # scale_color_manual(values = ggplot2::alpha(turbo(length(all_studies))[studies_here], 1), breaks = all_studies[studies_here] )+
  # scale_fill_manual(values = ggplot2::alpha(turbo(length(all_studies))[studies_here], 0.4), breaks = all_studies[studies_here])+
  # ylim(min(world_coordinates$lat), max(world_coordinates$lat)) +
  # coord_quickmap(#xlim = c(min(world_coordinates$long), 180),
  #                ylim = c(min(world_coordinates$lat), max(world_coordinates$lat))) +
  ylim(-50, 78) +
  xlim(-155, 175)+
  # coord_map(projection = "gilbert", xlim=c(-180,180), 
  #           ylim = c(-54, 77)) +
  labs(size = "Sample count") +
  guides(colour = "none", fill = "none") +
  scale_size(breaks = c(500, 1000, 2000)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        # panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=15),
        legend.position = c(0.12,0.15),
        panel.background = element_rect(fill = "#BDE3FF"),
        legend.text = element_text(size=15),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))
ggsave("study_locations.pdf", device = "pdf", path = "/fast/AG_Forslund/rob/studies/", width = 19.2, height = 12.8, units = "cm")


# combined 16S and shotgun:
data_loc_comb <- rbind(data_loc_16S, data_loc_shotgun) %>% 
  # mutate(value = Freq) %>%
  arrange(-sample_count) #%>%
# filter(Sequencing == "16S")#%>%

ggplot() +
  geom_map(
    data = world_coordinates, map = world_coordinates,
    aes(long, lat, map_id = region), # , alpha = hdi
    color = "white", fill = "grey", size = 0.15
  )+
  geom_point(data = data_loc_comb,
             aes(x = longitude, y = latitude, size = sample_count, color = Type, fill = Type),
             alpha = 0.6,
             position = position_jitter(width=0.7),
             shape=21) +
  # scale_fill_manual(values = c("#E1DD37FF", "red")) +
  # scale_color_manual(values = c("#E1DD37FF", "red")) +
  ylim(-50, 78) +
  xlim(-155, 175)+
  # coord_map(projection = "gilbert", xlim=c(-180,180), 
  #           ylim = c(-54, 77)) +
  labs(size = "Samples") +
  guides(colour = "none")+
  scale_size(breaks = c(500, 1000, 2000)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        # panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=15),
        legend.position = c(0.12,0.25),
        panel.background = element_rect(fill = "#BDE3FF"),
        legend.text = element_text(size=15))
ggsave("all_study_locations.pdf", device = "pdf", path = "/fast/AG_Forslund/rob/studies/", width = 19.2, height = 12.8, units = "cm")
ggsave("all_study_locations.png", path = "/fast/AG_Forslund/rob/studies/", width = 19.2, height = 12.8, units = "cm")
# ggsave("all_study_locations_poster.pdf", device = "pdf", path = "/fast/AG_Forslund/rob/studies/", width = 19.2, height = 12.8, units = "cm")

# only 16S, with lifestyle
map_16S <- ggplot() +
  geom_map(
    data = world_coordinates, map = world_coordinates,
    aes(long, lat, map_id = region), # , alpha = hdi
    color = "white", fill = "grey", size = 0.15
  )+
  geom_point(data = data_loc_comb %>% filter (Type == "16S") %>%
               mutate(Lifestyle = ifelse(lifestyle == "industrialized", yes = "Industrialized", no = "Non-industrialized")),
             aes(x = longitude, y = latitude, size = sample_count, color = lifestyle, fill = Lifestyle),
             alpha = 0.6,
             position = position_jitter(width=0.7),
             shape=21) +
  # scale_fill_manual(values = c("#E1DD37FF", "red")) +
  # scale_color_manual(values = c("#E1DD37FF", "red")) +
  ylim(-50, 78) +
  xlim(-155, 175)+
  # coord_map(projection = "gilbert", xlim=c(-180,180), 
  #           ylim = c(-54, 77)) +
  labs(size = "Samples") +
  guides(colour = "none")+
  scale_size(breaks = c(500, 1000, 2000)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        # panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=12),
        legend.position = c(0.12,0.25),
        panel.background = element_rect(fill = "#BDE3FF"),
        legend.text = element_text(size=10))
ggsave("/fast/AG_Forslund/rob/studies/16S_study_locations_ls.pdf",
       width = 23.04, height = 13, units = "cm", plot = map_16S)
ggsave("/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/16S_study_locations_ls.png",
       width = 19.2, height = 12.8, units = "cm", plot = map_16S)
save(world_coordinates, data_loc_comb, map_16S, 
     file = "/fast/AG_Forslund/rob/mm_index/R_scripts/paper_figures/16S_study_locations_ls.RData")
save(world_coordinates, data_loc_comb, 
     file = "/fast/AG_Forslund/rob/mm_index/shotgun_R_scripts/paper_figures/shotgun_study_locations_ls.RData")

# compare hdi and gdp between industrialized and non-industrialized:
countries_16S <- data_loc_comb %>% filter (Type == "16S") %>% select(country, lifestyle) %>%
  distinct() %>%
  left_join(.,hdi_data, by = c("country" = "region")) %>%
  left_join(.,gdp_data, by = c("country" = "Country.Name"))

countries_16S %>%
  pivot_longer(cols = c("hdi", "GDP"), names_to = "measurement", values_to = "value") %>%
  ggplot(., aes(y = value, color = lifestyle, x = lifestyle)) +
  geom_boxplot() + 
  geom_jitter() +
  facet_wrap(~measurement, scales = "free")


# only shotgun, with lifestyle
ggplot() +
  geom_map(
    data = world_coordinates, map = world_coordinates,
    aes(long, lat, map_id = region), # , alpha = hdi
    color = "white", fill = "grey", size = 0.15
  )+
  geom_point(data = data_loc_comb %>% filter (Type == "shotgun"),
             aes(x = longitude, y = latitude, size = sample_count, color = lifestyle, fill = lifestyle),
             alpha = 0.6,
             position = position_jitter(width=0.7),
             shape=21) +
  # scale_fill_manual(values = c("#E1DD37FF", "red")) +
  # scale_color_manual(values = c("#E1DD37FF", "red")) +
  ylim(-50, 78) +
  xlim(-155, 175)+
  # coord_map(projection = "gilbert", xlim=c(-180,180), 
  #           ylim = c(-54, 77)) +
  labs(size = "Samples") +
  guides(colour = "none")+
  scale_size(breaks = c(500, 1000, 2000)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        # panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=15),
        legend.position = c(0.12,0.2),
        panel.background = element_rect(fill = "#BDE3FF"),
        legend.text = element_text(size=15))
ggsave("shotgun_study_locations_ls.pdf", device = "pdf", path = "/fast/AG_Forslund/rob/studies/", width = 23.04, height = 15.36, units = "cm")
ggsave("shotgun_study_locations_ls.png", path = "/fast/AG_Forslund/rob/studies/", width = 19.2, height = 12.8, units = "cm")



# Metadata/studies summary ################################################################################
library(readODS)
library(data.table)
setwd(maindir)
studies_metadata <- readODS::read_ods("/fast/AG_Forslund/rob/Documents/interesting_studies.ods", sheet = "metadata_available") %>%
  tibble::column_to_rownames(var = "metadata_colnames:")
# studies_metadata <- read.table("studies_metadata.txt", sep = "\t", header = T, row.names = 1)
studies_metadata <- studies_metadata[,-c(1, 2)]

studies_metadata[studies_metadata == "x"] <- T
studies_metadata[is.na(studies_metadata)] <- F
studies_metadata <- type.convert(studies_metadata)

studies_metadata["gestational_age",] <- studies_metadata["gestational_age",] | studies_metadata["gestational_age(estimated)",]
# combine antibiotics information
studies_metadata["antibiotics",] <- studies_metadata["antibiotics_before",] | studies_metadata["antibiotics_any",] |
  studies_metadata["antibiotics_one_week_before",]
studies_metadata <- studies_metadata[!rownames(studies_metadata) %in% 
                                       c("gestational_age(estimated)",
                                         "run_accession", "sample_ID", "region",
                                         "study", "family_ID", "batch", "antibiotics_before",
                                         "antibiotics_any", "antibiotics_one_week_before",
                                         "birth_height, firt_solid_food", "height",
                                         "birth_height"),
                                     !colnames(studies_metadata) %in% c(#"Pärnänen_shotgun",
                                       "Vatanen22_shotgun",
                                       "cosmic_shotgun",
                                       # "Bargheet_shotgun",
                                       "Gasparrini_shotgun")]



#studies_metadata <- studies_metadata[order(rowSums(studies_metadata), decreasing = T),]
rownames_ordered <- rownames(studies_metadata)[order(rowSums(studies_metadata, na.rm = T))]
colnames_ordered <- colnames(studies_metadata)[order(colSums(studies_metadata, na.rm = T))]

studies_metadata <- studies_metadata[,unlist(studies_metadata[rownames(studies_metadata) == "subject_ID",])]
studies_metadata$meta_var <- rownames(studies_metadata)
# remove studies without subject_ID information
studies_metadata_long <- gather(studies_metadata, key = "study", value = "availability", -meta_var) %>%
  mutate(meta_var=factor(meta_var, levels = rownames_ordered)) %>%
  mutate(study=factor(study, levels = colnames(studies_metadata)[-length(colnames(studies_metadata))]))
studies_metadata_long <- separate(studies_metadata_long, col = "study", sep = "_", into = c("study", "sequencing"))
studies_metadata_long$study <- factor(studies_metadata_long$study, levels = unique(gsub("_.*", "", colnames_ordered)))
# shorten study names:
short_names <- gsub("_.*", "", colnames_ordered) %>%
  unique %>%
  strtrim(., 8)
studies_metadata_long$study <- strtrim(studies_metadata_long$study, 8) %>% 
  factor(., levels = short_names)
ggplot(studies_metadata_long, aes(x = study, y = meta_var, fill = availability))+
  geom_tile(colour="grey", size=0.5, width = 1, height = 1)+
  #remove x and y axis labels
  labs(x="", y="")+
  #remove extra space
  # scale_y_discrete(expand=c(0, 0))+
  #set a base size for all fonts
  # theme_grey(base_size=15)+
  scale_fill_manual(values=c("white", "black"))+
  facet_wrap(~sequencing, strip.position = "bottom", scales = "free_x") +
  #theme options
  theme(
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank(),
    axis.text.x = element_text(angle = 50, vjust = 1, hjust=1, size = 12),
    axis.text.y = element_text(size = 12),#, face=c("plain","bold","plain", "plain","plain", "plain","plain",
    # "plain","plain","plain", "bold","plain","bold", "plain",
    # "bold","plain","plain")),
    legend.position="none"
  )
ggsave("studies_x_metadata.pdf", device = "pdf", path = maindir)


asv_data <- readODS::read_ods("/fast/AG_Forslund/rob/Documents/interesting_studies.ods", sheet = "complete_16S") %>%
  filter(!is.na(N_ASVs)) %>%
  tibble::column_to_rownames(var = "author") %>%
  select(n_individuals_final, n_samples, N_ASVs) %>%
  mutate(n_samples = as.numeric(n_samples))
# asv_data <- read.table("n_ASVs.txt", header = T, row.names = 1)
ggplot(asv_data, aes(x = n_samples, y = N_ASVs, color = n_individuals_final)) +
  geom_point() +
  scale_color_gradientn(name = "n_individuals", colors = topo.colors(7)) +
  geom_text(
    label=rownames(asv_data), 
    vjust = -0.25,
    hjust = 0.4
  ) +
  scale_x_continuous(n.breaks = 8) +
  scale_y_continuous(n.breaks = 10, limits=c(0,17000)) +
  theme_minimal()

ggsave("dada_2_asvs.pdf", device = "pdf", path = maindir)


# reads lost in each sample over steps for 16S
setwd(paste0(maindir, "16S/"))
stat_files <- list.files(paste0(list.files(), "/dada2"), pattern = "read_loss_per_step.csv", full.names = T)
loss_table <- data.frame(X = integer(),
                         Row.names = character(),
                         read_counts_raw = double(),
                         decontaminated = double(),
                         filtered = double(),
                         denoised = double(),
                         removed_chimeras = double(),
                         study = character())
for(f in stat_files){
  add_table <- read.csv(f)
  add_table$study <- gsub("/dada2.*", "", f)
  loss_table <- rbind(loss_table, add_table)
}
raw_counts <- loss_table$read_counts_raw
for(c in colnames(loss_table)[3:7]){
  loss_table[,c] <- loss_table[,c]/raw_counts
}
loss_table_long <- pivot_longer(loss_table, -c(X, Row.names, study), values_to = "reads", names_to = "step")

loss_table[,c("read_counts_raw", "decontaminated",  "filtered",  "denoised", "removed_chimeras")] <- 
  1-loss_table[,c("read_counts_raw", "decontaminated",  "filtered",  "denoised", "removed_chimeras")]

loss_table$read_counts_raw <- NULL


loss_table_long <- pivot_longer(loss_table, -c(X, Row.names, study), values_to = "reads", names_to = "step")
ggplot(loss_table_long, aes(x=reads, col = step)) + 
  geom_density() +
  ylab("density") +
  xlim(0,1) + 
  ylim(0,11) +
  xlab("lost reads per sample") 



ggplot(loss_table_long) +
  geom_line(aes(factor(step, levels = colnames(loss_table)[3:7]), reads, group = Row.names, color = study)) +
  scale_y_continuous(trans = 'log2')


# overlapping important features:
important_features <- readODS::read_ods("/fast/AG_Forslund/rob/Documents/interesting_studies.ods", sheet = "compare important features") %>%
  select(-`all taxa`)
all_taxa <- unique(unlist(important_features))
rank_df <- data.frame(entry = factor(all_taxa[!is.na(all_taxa)], levels = all_taxa[!is.na(all_taxa)]))

for (col_name in colnames(important_features)) {
  # Get the rank of each entry in the current column (using rank function)
  # If the entry is not present in the column, it will be assigned NA
  rank_df[[col_name]] <- sapply(rank_df$entry, function(x) {
    if (x %in% important_features[[col_name]]) {
      return(which(important_features[[col_name]] == x))
    } else {
      return(NA)
    }
  })
}

rank_df_long <- rank_df %>%
  pivot_longer(cols = -entry, names_to = "column", values_to = "rank") %>%
  mutate(column = factor(column, levels = unique(c("my model", column))))

# Step 5: Plot the heatmap
ggplot(rank_df_long, aes(x = column, y = entry, fill = rank)) +
  geom_tile(color = "white", size = 0.1) +  # Create tiles with borders
  scale_fill_gradient(low = "red", high = "blue", na.value = "gray") +  # Color scale for ranks
  theme_minimal() +  # Clean theme
  labs(title = "Rank Heatmap", x = "Columns", y = "Entries") +  # Labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
ggsave("taxa_importance_ranks.pdf", device = "pdf", path = maindir)

