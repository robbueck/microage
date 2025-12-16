rm(list = ls())
library(filesstrings)
library(tidyverse)
library(tidyr)
library(dplyr)
library(sitar)
library("magrittr")
# convert all male/female versions to m/f
# convert all yes/no versions to y/n
# c-section/vaginal: c/v
# create ages in days and months
# feeding: breast-feeding/formula-feeding/solid-food/unknown/nothing breast/formula/solid/unknown/nothing
# mixed feeding: breast/formula/solid  mixed_b_f_s

# 16S studies
maindir <- "/fast/AG_Forslund/rob/studies/"

# Hesla 2014 (waiting for metadata)#############################################
# only SraRunTable.txt available
setwd(paste(maindir, "16S/hesla_2014/", sep = ""))
ena_file_hesla <- read.table("filereport_read_run_PRJNA263853_tsv.txt", sep = "\t", header = T)
ena_file_hesla <- ena_file_hesla[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]

metadata_hesla <- read.table("SraRunTable.txt", sep = ",", header = T)
# remove columns
metadata_hesla[,c("Assay.Type", "BioSampleModel", "Bytes", "Center.Name", "Collection_Date",
                  "Consent", "DATASTORE.filetype", "DATASTORE.provider", "DATASTORE.region",
                  "env_biome", "env_feature", "env_material", "geo_loc_name_country",
                  "geo_loc_name_country_continent", "geo_loc_name", "Host",
                  "Instrument", "LibraryLayout", "LibrarySelection", "Lat_Lon",
                  "LibrarySource", "Organism", "ReleaseDate", "samp_store_loc",
                  "samp_store_temp", "Library.Name", "Sample_code")] <- NULL

# no duplicates
sum(duplicated(metadata_hesla$Sample.Name))

# remove mother samples
sort_out_samples_hesla <- metadata_hesla %>% filter(Host_Age == "Adult") %>%
   select(Run)

metadata_hesla %<>% filter(Host_Age != "Adult")
# individuals with only one sample
metadata_hesla$subject_ID <- substr(metadata_hesla$Sample.Name, 1, 3)

multiple_samples <- metadata_hesla %>%
   filter(duplicated(subject_ID)) %>%
   select(subject_ID) %>%
   unique %>%
   deframe

sort_out_samples_hesla <- metadata_hesla %>%
   filter(!subject_ID %in% multiple_samples) %>%
   select(Run) %>%
   rbind(sort_out_samples_hesla)

metadata_hesla %<>% filter(subject_ID %in% multiple_samples)

# age column:
metadata_hesla %<>% mutate(age = case_when(
   Host_Age == "3 days old" ~ 3,
   Host_Age == "3 weeks old" ~ 21,
   Host_Age == "2 months old" ~ 60,
   Host_Age == "6 months old" ~ 180
), .keep = "unused")


# unify colnames
metadata_hesla %<>% mutate(run_accession = Run,
                           sample_ID = paste0(Sample.Name, "_hesla"),
                           subject_ID = paste0(subject_ID, "_hesla"),
                           country = "SWEDEN",
                           region = "STOCKHOLM",
                           study = "hesla_2014",
                           lifestyle = "industrialized",
                           geographic_location_.latitude. = 59.090000,
                           geographic_location_.longitude. = 17.560000,
                           .keep = "unused")

metadata_hesla$family_ID <- metadata_hesla$subject_ID
sort(colnames(metadata_hesla))


# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_hesla$Bases, n_lines$read_counts[match(metadata_hesla$run_accession, n_lines$sample_ID)])  # 


write.csv(metadata_hesla, "metadata_hesla_2014_healthy.csv")

# move sorted out samples
sort_out_samples_hesla <- unique(sort_out_samples_hesla)
sort_out_files_hesla <- c(paste(sort_out_samples_hesla, "_1.fastq.gz", sep = ""),
                         paste(sort_out_samples_hesla, "_2.fastq.gz", sep = ""),
                         paste(sort_out_samples_hesla, ".fastq.gz", sep = ""))
sort_out_files_hesla <- sort_out_files_hesla[file.exists(sort_out_files_hesla)]
file.move(sort_out_files_hesla, "sorted_out/")

# move good files
final_samples_hesla <- metadata_hesla$run_accession
final_files_hesla <- c(paste(final_samples_hesla, "_1.fastq.gz", sep = ""),
                      paste(final_samples_hesla, "_2.fastq.gz", sep = ""),
                      paste(final_samples_hesla, ".fastq.gz", sep = ""))
final_files_hesla <- final_files_hesla[file.exists(final_files_hesla)]
file.move(final_files_hesla, "fastq_files/")

# clean up
rm(ena_file_hesla, metadata_hesla, final_files_hesla, final_samples_hesla,
   multiple_samples, sort_out_files_hesla, sort_out_samples_hesla, n_lines)


# Subramanian 2014 #############################################################

setwd(paste(maindir, "16S/subramanian_2014/", sep = ""))
ena_file_subr <- read.table("filereport_read_run_PRJEB5482_tsv.txt", sep = "\t", header = T)
ena_file_subr <- ena_file_subr[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]


metadata_sub_1 <- read.table("suppl_1.txt", sep = "\t", header = T)
metadata_sub_1[c("WHZ", "X", "X.1", "WAZ", "X.2", "X.3", "HAZ", "X.4", "X.5",
                 "Sampling.interval..days..mean...SD", "X.6", "X.7", "Family.ID")] <- NULL
metadata_sub_2 <- read.table("suppl_2.txt", sep = "\t", header = T)
metadata_sub_3 <- read.table("suppl_3.txt", sep = "\t", header = T) # parent samples
metadata_sub_10 <- read.table("suppl_10_treat.txt", sep = "\t", header = T) # SAM samples
metadata_sub_10[c("Food.Intervention.Assignment.x")] <- NULL
metadata_sub_11 <- read.table("suppl_11_treat.txt", sep = "\t", header = T, quote = "")  # SAM samples
metadata_sub_17 <- read.table("suppl_17_treat.txt", sep = "\t", header = T)  # information is already contained in tables above, also contains MAM children

sum(!metadata_sub_2$Fecal.Sample.ID %in% ena_file_subr$sample_alias) # check if everything fits

metadata_sub_combined <- left_join(metadata_sub_1, metadata_sub_2 , by = "Child.ID") %>% 
   right_join(., ena_file_subr, by = c("Fecal.Sample.ID" = "sample_alias"))
# for some samples, more than one library run is available, analyze them separately
# dupls <- metadata_sub_combined[metadata_sub_combined$Fecal.Sample.ID %in% ena_file_subr$sample_title[duplicated(ena_file_subr$sample_title)],]

# duplicated samples seem to have lower read counts, so maybe merge them
dupls_ena <- ena_file_subr[duplicated(ena_file_subr$sample_title),]
hist(ena_file_subr$read_count[ena_file_subr$sample_title %in% dupls_ena$sample_title], breaks = seq(0, 65000, 1000))
hist(ena_file_subr$read_count[!ena_file_subr$sample_title %in% dupls_ena$sample_title], xlim = c(0,65000), breaks = 1000)
hist(ena_file_subr$read_count, xlim = c(0,65000), breaks = 1000)


length(unique(metadata_sub_combined$Fecal.Sample.ID))  # number of unique samples
metadata_sub_combined <- metadata_sub_combined %>% replace(. == "Yes", T) %>% 
  replace(. == "No", F) %>%
  replace(. == "Female", "f") %>%
  replace(. == "Male", "m") %>%
  type.convert()

#convert from months to days
metadata_sub_combined %<>% mutate(Months.of.exclusive.breastfeeding = Months.of.exclusive.breastfeeding *30,
                                 Age.at.first.introduction.of.solid.food..months. = Age.at.first.introduction.of.solid.food..months. * 30)


# remove unused columns
metadata_sub_combined[c("Age.at.first.fecal.sample.collection..days.",
                        "Age.at.last.fecal.sample.collection..days.",
                        "Number.of.fecal.samples.collected", "Age..months")] <- NULL

# unify colnames
metadata_sub_combined %<>% mutate(sample_ID = paste0(Fecal.Sample.ID, "_subramanian"),
                                 subject_ID = paste0(Child.ID, "_subramanian"),
                                 sex = Gender,
                                 exclusive_breastfeeding = Months.of.exclusive.breastfeeding,
                                 first_solid_food = Age.at.first.introduction.of.solid.food..months.,
                                 family_ID = Family.ID,
                                 age = Age..days,
                                 breast = Breast.Milk,
                                 formula = Formula1,
                                 solid = Solid.Foods2,
                                 antibiotics_one_week_before = Antibiotics.within.7.days.prior.to.sample.collection,
                                 country = "BANGLADESH",
                                 region = "DHAKA",
                                 Instrument = "Illumina MiSeq",
                                 study = "subramanian_2014",
                                 lifestyle = "non_industrialized",
                                 .keep = "unused")

# food column
metadata_sub_combined %<>% mutate(food = case_when(
   breast &  !formula & !solid ~ "breast",
   !breast & formula & !solid ~ "formula",
   !breast & !formula & solid ~ "solid",
   breast &  formula & !solid ~ "mixed_b_f",
   breast &  !formula & solid ~ "mixed_b_s",
   !breast &  formula & solid ~ "mixed_f_s",
   breast &  formula & solid ~ "mixed_f_b_s"
   ))

# sort out samples
# sort_out_samples <- metadata_sub_combined %>%
#    filter(sample_ID %in% paste0(metadata_sub_3$Fecal.Sample.ID, "_subramanian") | 
#              sample_ID %in% paste0(metadata_sub_11$Fecal.Sample.ID, "_subramanian")) %>%
#    select(run_accession)
# sort_out_samples <- metadata_sub_combined %>%
#    filter(grepl("Mock", sample_ID) | 
#              is.na(subject_ID)) %>%
#    select(run_accession) %>%
#    rbind(sort_out_samples)
# 
# sort_out_samples <- metadata_sub_combined %>%
#   filter(diarrhoea.at.the.time.of.sample.collection3 == "y") %>%
#   select(run_accession) %>%
#   rbind(sort_out_samples)

samples_raman <- read.csv("../raman_2019/metadata_ram_2019_healthy.csv")
samples_raman %<>% mutate(sample_ID = gsub("_raman_2019", "", sample_ID))
# sort_out_samples <- metadata_sub_combined %>%
#    filter(sample_ID %in% paste0(samples_raman$sample_ID, "_subramanian")) %>%
#    select(run_accession) %>% 
#    rbind(sort_out_samples)


metadata_sub_combined %<>% filter(!sample_ID %in% paste0(metadata_sub_3$Fecal.Sample.ID, "_subramanian"),
                                  !sample_ID %in% paste0(metadata_sub_11$Fecal.Sample.ID, "_subramanian"),
                                  !grepl("Mock", sample_ID),
                                  !is.na(subject_ID),
                                  is.na(diarrhoea.at.the.time.of.sample.collection3) | 
                                    diarrhoea.at.the.time.of.sample.collection3 == F,
                                  subject_ID != "NA_subramanian",
                                  !sample_ID %in% paste0(samples_raman$sample_ID, "_subramanian"))



# individuals with only one sample
multiple_samples <- unique(metadata_sub_combined$subject_ID[duplicated(metadata_sub_combined$subject_ID)])
length(metadata_sub_combined$run_accession[!metadata_sub_combined$subject_ID %in% multiple_samples])  # all individuals have multiple samples


# add height and weight data:
metadata_sub_combined$height <- sitar::LMS2z(metadata_sub_combined$age/365, y = metadata_sub_combined$HAZ, sex = metadata_sub_combined$sex,
                                      measure = "ht", ref = who06, toz = F)
metadata_sub_combined$weight <- sitar::LMS2z(metadata_sub_combined$age/365, y = metadata_sub_combined$WAZ, sex = metadata_sub_combined$sex,
      measure = "wt", ref = who06, toz = F) 

# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_sub_combined$read_count, n_lines$read_counts[match(metadata_sub_combined$run_accession, n_lines$sample_ID)])  # 


# remove redundant cols"
metadata_sub_combined[,c("Cohort", "Training.Validation.Set.....Subject.Allocation")] <- NULL
# add location information
sra_sub <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_sub <- sra_sub[,c("Run", "geographic_location_.latitude.", "geographic_location_.longitude.")]
metadata_sub_combined <- left_join(metadata_sub_combined, sra_sub, by = join_by(run_accession == Run))
# for now, keep duplicated samples separate
write.csv(metadata_sub_combined, "metadata_subramanian_2014_healthy.csv")


metadata_sub_mal_combined <- merge(metadata_sub_10, metadata_sub_11, by = "Child.ID", all.x = T)
sum(!metadata_sub_mal_combined$Fecal.Sample.ID %in% ena_file_subr$sample_alias)  # check if everything fits
nrow(metadata_sub_mal_combined)
metadata_sub_mal_combined <- merge(metadata_sub_mal_combined, ena_file_subr, by.x = "Fecal.Sample.ID", by.y = "sample_alias")
nrow(metadata_sub_mal_combined)
length(unique(metadata_sub_mal_combined$Fecal.Sample.ID))
metadata_sub_mal_combined[metadata_sub_mal_combined == "Male"] <- "m"
metadata_sub_mal_combined[metadata_sub_mal_combined == "Female"] <- "f"
metadata_sub_mal_combined[metadata_sub_mal_combined == "Yes"] <- "y"
metadata_sub_mal_combined[metadata_sub_mal_combined == "No"] <- "n"
metadata_sub_mal_combined$Age.of.child.at.time.of.fecal.sample.collection..months. <- metadata_sub_mal_combined$Age.of.child.at.time.of.fecal.sample.collection..months. * 30


# unify colnames
colnames(metadata_sub_mal_combined)[match(c("Fecal.Sample.ID", "Child.ID", "Gender",
                                            "Weight..kg", "Height..cm", "Weight.for.Height.Z.score..WHZ.",
                                            "Height.for.Age.Z.score..HAZ.", "Weight.for.age.Z.score..WAZ.",
                                            "Age.of.child.at.time.of.fecal.sample.collection..months.", "Antibiotics"),
                                          colnames(metadata_sub_mal_combined))] <- c("sample_ID", "subject_ID", "sex",
                                                                              "weight", "height", "WHZ", "HAZ", "WAZ",
                                                                              "age", "antibiotics_any")
metadata_sub_mal_combined <- metadata_sub_mal_combined %>%
  mutate(country = "BANGLADESH",
         region = "DHAKA",
         study = "subramanian_2014",
         health = "SAM",
         WHZ = as.numeric(WHZ),
         WAZ = as.numeric(WAZ),
         HAZ = as.numeric(HAZ))
sort(colnames(metadata_sub_combined))


write.csv(metadata_sub_mal_combined, "metadata_subramanian_2014_condition.csv")
# move malnurished samples:
malnurished_samples <- metadata_sub_mal_combined$run_accession
malnurished_files <- c(paste(malnurished_samples, "_1.fastq.gz", sep = ""),
                    paste(malnurished_samples, "_2.fastq.gz", sep = ""),
                    paste(malnurished_samples, ".fastq.gz", sep = ""))
malnurished_files <- malnurished_files[file.exists(malnurished_files)]
file.move(malnurished_files, "malnurished_fastq_files/")


# 
# # move good files
# final_samples_sub <- metadata_sub_combined$run_accession
# final_files_sub <- c(paste(final_samples_sub, "_1.fastq.gz", sep = ""),
#                      paste(final_samples_sub, "_2.fastq.gz", sep = ""),
#                      paste(final_samples_sub, ".fastq.gz", sep = ""))
# final_files_sub <- final_files_sub[file.exists(final_files_sub)]
# file.move(final_files_sub, "fastq_files/")


# clean up
rm(ena_file_subr, metadata_sub_1, metadata_sub_2, metadata_sub_3, metadata_sub_10, metadata_sub_11,
   metadata_sub_17, metadata_sub_combined, dupls_ena, malnurished_files, malnurished_samples,
   metadata_sub_mal_combined, sort_out_samples, sort_out_files, multiple_samples,
   final_samples_sub, final_files_sub, n_lines, samples_raman, sra_sub)



# Kristensen	2020 #############################################################

setwd(paste(maindir, "16S/kristensen_2020/", sep = ""))
ena_file_kristensen <- read.table("filereport_read_run_PRJNA554232_tsv.txt", sep = "\t", header = T)
ena_file_kristensen <- ena_file_kristensen[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
metadata_kristensen <- read.table("SraRunTable.txt", sep = ",", header = T)
# remove columns
metadata_kristensen[,c("Assay.Type", "BioSampleModel", "broad.scale_environmental_context",
                       "Bytes", "Center.Name", "Collection_Date", "Consent",
                       "DATASTORE.filetype", "DATASTORE.provider", "DATASTORE.region",
                       "environmental_medium", "Forward_sample", "geo_loc_name_country",
                       "geo_loc_name_country_continent", "Geographic_location",
                       "Host", "Isolate", "latitude_and_longitude",
                       "Isolation_Source", "LibraryLayout", "LibrarySelection",
                       "LibrarySource", "Organism", "Platform", "ReleaseDate",
                       "Reverse_sample")] <- NULL

# no duplicate samples
sum(duplicated(metadata_kristensen$Sample.Name))

# remove mother samples
sort_out_samples_kris <- metadata_kristensen$Run[metadata_kristensen$Sample_Group == "CF"]
metadata_kristensen <- metadata_kristensen[metadata_kristensen$Sample_Group != "CF",]
# individuals with only one sample
multiple_samples <- unique(metadata_kristensen$SUBJECT_ID[duplicated(metadata_kristensen$SUBJECT_ID)])
length(metadata_kristensen$Run[!metadata_kristensen$SUBJECT_ID %in% multiple_samples])  # all individuals have multiple samples



metadata_kristensen %<>% mutate(run_accession = Run,
                                sample_ID = paste0(Sample.Name, "_kristensen"),
                                subject_ID = paste0(SUBJECT_ID, "_kristensen"),
                                age = Sample_Timepoint * 30,
                                country = "NETHERLANDS",
                                region = "NETHERLANDS",
                                study = "kristensen_2020",
                                lifestyle = "industrialized",
                                family_ID = subject_ID,
                                gestational_age = 39.86,
                                geographic_location_.latitude. = 52.092876,
                                geographic_location_.longitude. = 5.104480,
                                .keep = "unused")

# abx information visually taken from "Development of the Nasopharyngeal Microbiota in Infants with Cystic Fibrosis"
first_abx <- data.frame(individual = c("2007_kristensen", "2013_kristensen", "2017_kristensen", "2035_kristensen",
                                       "2036_kristensen", "2038_kristensen"),
                        first_abx = c(150, 180, 150, 180, 180, 150))
metadata_kristensen <- left_join(metadata_kristensen, first_abx, by = join_by(subject_ID == individual))
metadata_kristensen <- metadata_kristensen %>%
  mutate(antibiotics_before = case_when(first_abx < age ~ T,
                                        .default = F))


# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_kristensen$Bases, n_lines$read_counts[match(metadata_kristensen$run_accession, n_lines$sample_ID)])  # 


write.csv(metadata_kristensen, "metadata_kristensen_2020_healthy.csv")

# move sorted out samples
sort_out_samples_kris <- unique(sort_out_samples_kris)
sort_out_files_kris <- c(paste(sort_out_samples_kris, "_1.fastq.gz", sep = ""),
                        paste(sort_out_samples_kris, "_2.fastq.gz", sep = ""),
                        paste(sort_out_samples_kris, ".fastq.gz", sep = ""))
sort_out_files_kris <- sort_out_files_kris[file.exists(sort_out_files_kris)]
file.move(sort_out_files_kris, "sorted_out/")

# move good files
final_samples_kris <- metadata_kristensen$run_accession
final_files_kris <- c(paste(final_samples_kris, "_1.fastq.gz", sep = ""),
                     paste(final_samples_kris, "_2.fastq.gz", sep = ""),
                     paste(final_samples_kris, ".fastq.gz", sep = ""))
final_files_kris <- final_files_kris[file.exists(final_files_kris)]
file.move(final_files_kris, "fastq_files/")

# clean up
rm(ena_file_kristensen, metadata_kristensen, final_files_kris, final_samples_kris,
   multiple_samples, sort_out_files_kris, sort_out_samples_kris, n_lines, first_abx)


# Roswall 2021 (connected with backhed 2015) ###################################
# gestational age only estimated

setwd(paste(maindir, "16S/roswall_2021/", sep = ""))
ena_file_roswall <- read.table("filereport_read_run_PRJEB38986_tsv.txt", sep = "\t", header = T)
ena_file_roswall <- ena_file_roswall[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]

# read and reformat data
metadata_ros <- read.table("suppl_data.txt", sep = "\t", header = T)
metadata_ros$child_id <- gsub("[-].*", "", metadata_ros$sample_alias)
metadata_ros$sample_time <- gsub(".*[-]", "", metadata_ros$sample_alias)
read_counts_ros <- read.table("sample_id_reads.txt", sep = "\t", header = T)
read_counts_ros$Sample.ID <- gsub("[_]", "-", read_counts_ros$Sample.ID)
read_counts_ros$Sample.ID <- gsub("[-]B", "-NB", read_counts_ros$Sample.ID)
dict_ros_back <- read.table("dict.txt", sep = "\t", header = T)
metadata_ros_back <- read.table("suppl_data_2015.txt", sep = "\t", header = T)
metadata_ros_back$child_id_new <- dict_ros_back$STUDY.ID.Growth.study..H2GS[match(metadata_ros_back$STUDY.ID, dict_ros_back$STUDY.ID.metagenomics)]

# check read counts from study
sum(!metadata_ros$sample_alias %in% read_counts_ros$Sample.ID)
sum(!read_counts_ros$Sample.ID %in% metadata_ros$sample_alias)
metadata_ros_combined <- merge(metadata_ros, read_counts_ros, by.x = "sample_alias", by.y = "Sample.ID", all.x = T)
metadata_ros_combined <- merge(metadata_ros_combined, ena_file_roswall, by = "sample_alias")
plot(metadata_ros_combined$read_count, metadata_ros_combined$Number.of.reads)

# merge with backhed 2015 data
sum(!metadata_ros_back$child_id_new %in% metadata_ros_combined$child_id)
sum(!metadata_ros_combined$child_id %in% metadata_ros_back$child_id_new)
nrow(metadata_ros_combined)
metadata_ros_combined <- merge(metadata_ros_combined, metadata_ros_back, by.x = "child_id", by.y = "child_id_new", all.x = T)
nrow(metadata_ros_combined)
sum(is.na(metadata_ros_combined))
metadata_ros_combined[metadata_ros_combined == "vaginal"] <- "v"
metadata_ros_combined[metadata_ros_combined == "sectio" | metadata_ros_combined == "section"] <- "c"
metadata_ros_combined[metadata_ros_combined == "yes"] <- T
metadata_ros_combined[metadata_ros_combined == "no"] <- F
metadata_ros_combined[metadata_ros_combined == "boy"] <- "m"
metadata_ros_combined[metadata_ros_combined == "girl"] <- "f"
metadata_ros_combined[metadata_ros_combined == "exclusively breastfeeding"] <- "breast"
metadata_ros_combined[metadata_ros_combined == "mixed feeding"] <- "mixed_b_f"
metadata_ros_combined[metadata_ros_combined == "exclusively formula feeding" | metadata_ros_combined == "Formula feeding"] <- "formula"
metadata_ros_combined[metadata_ros_combined == "no breastfeeding"] <- "solid"
metadata_ros_combined[metadata_ros_combined == "any breastfeeding"] <- "mixed_b_s"
metadata_ros_combined$sex[metadata_ros_combined$sex == "1"] <- "m"
metadata_ros_combined$sex[metadata_ros_combined$sex == "2"] <- "f"
sum(is.na(metadata_ros_combined))

# remove mother samples
sort_out_samples_ros <- metadata_ros_combined$run_accession[metadata_ros_combined$sex == "mother"]
metadata_ros_combined <- metadata_ros_combined[metadata_ros_combined$sex != "mother",]
# individuals with only one sample
multiple_samples <- unique(metadata_ros_combined$child_id[duplicated(metadata_ros_combined$child_id)])
# sort_out_samples_ros <- c(sort_out_samples_ros, 
#                           metadata_ros_combined$run_accession[!metadata_ros_combined$child_id %in% multiple_samples])
# metadata_ros_combined <- metadata_ros_combined[metadata_ros_combined$child_id %in% multiple_samples,]
# retreive sorted out samples
metadata_ros_combined$run_accession[!metadata_ros_combined$child_id %in% multiple_samples] %>%
  paste("fastq-dump --split-e --gzip ",.) %>%
  data.frame(`#!/bin/bash` = .) %>%
  write.table(., file = "/fast/AG_Forslund/rob/studies/16S/roswall_2021/pull_single.sh",
              sep = " ", quote = F, row.names = F, col.names = "#!/bin/bash")


# fill age column
metadata_ros_combined$age <- NA
metadata_ros_combined$age[metadata_ros_combined$sample_time == "NB"] <- metadata_ros_combined$Age.at.sample.Newborn..days.[metadata_ros_combined$sample_time == "NB"]
metadata_ros_combined$age[metadata_ros_combined$sample_time == "NB" & 
                                          is.na(metadata_ros_combined$age)] <- mean(metadata_ros_combined$Age.at.sample.Newborn..days., na.rm = T)
metadata_ros_combined$age[metadata_ros_combined$sample_time == "4M"] <- metadata_ros_combined$Age.at.Sampling..4.M[metadata_ros_combined$sample_time == "4M"]
metadata_ros_combined$age[metadata_ros_combined$sample_time == "4M" & 
                                          is.na(metadata_ros_combined$age)] <- mean(metadata_ros_combined$Age.at.Sampling..4.M, na.rm = T)
metadata_ros_combined$age[metadata_ros_combined$sample_time == "12M"] <- metadata_ros_combined$Age.at.Sampling.12M[metadata_ros_combined$sample_time == "12M"]
metadata_ros_combined$age[metadata_ros_combined$sample_time == "12M" & 
                                          is.na(metadata_ros_combined$age)] <- mean(metadata_ros_combined$Age.at.Sampling.12M, na.rm = T)
metadata_ros_combined$age[metadata_ros_combined$sample_time == "3Y"] <- 3*365
metadata_ros_combined$age[metadata_ros_combined$sample_time == "5Y"] <- 5*365


# add feeding column
metadata_ros_combined$food <- NA
metadata_ros_combined$food[metadata_ros_combined$timepoint == "NB"] <- metadata_ros_combined$feeding.practice.first.week[metadata_ros_combined$timepoint == "NB"]
metadata_ros_combined$food[metadata_ros_combined$timepoint == "4M"] <- metadata_ros_combined$feeding.practice.4M[metadata_ros_combined$timepoint == "4M"]
metadata_ros_combined$food[metadata_ros_combined$timepoint == "12M"] <- metadata_ros_combined$Any.breastfeeding.12.M[metadata_ros_combined$timepoint == "12M"]
metadata_ros_combined$food[metadata_ros_combined$timepoint == "3Y"] <- "solid"
metadata_ros_combined$food[metadata_ros_combined$timepoint == "5Y"] <- "solid"

# antibiotics
metadata_ros_combined$antibiotics_before <- F
metadata_ros_combined$antibiotics_before[metadata_ros_combined$age >= 100 & 
                                            metadata_ros_combined$Antibiotic.treatment.to.infant.0.4M...........times. > 0] <- T
metadata_ros_combined$antibiotics_before[metadata_ros_combined$age >= 340 & 
                                            metadata_ros_combined$Antibiotic.treatment.to.infant.4.12M.....times. > 0] <- T


# clean up metadata-table
colnames(metadata_ros_combined)
metadata_ros_combined[c("sample_time", "instrument_model", "library_name.x", "library_source",
                        "library_selection", "library_strategy", "forward_file_name",
                        "forward_file_md5", "reverse_file_name", "reverse_file_md5",
                        "X", "X.1", "library_name.y", "sample_title", "Age.at.sample.Newborn..days.",
                        "Age.at.Sampling..4.M", "Age.at.Sampling.12M", "feeding.practice.first.week",
                        "feeding.practice.4M", "Any.breastfeeding.12.M", "GENDER", "Delivery.mode")] <- NULL

sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file[,c("Run", "geographic_location_.latitude.", "geographic_location_.longitude.")]
metadata_ros_combined <- left_join(metadata_ros_combined, sra_file, by = join_by(run_accession == Run))

metadata_ros_combined %<>% mutate(sample_ID = paste0(sample_alias, "_roswall"),
                                subject_ID = paste0(child_id, "_roswall"),
                                birthmode = mode_of_birth,
                                breast = breastfed,
                                country = "SWEDEN",
                                region = "SWEDEN",
                                study = "roswall_2021",
                                lifestyle = "industrialized",
                                family_ID = subject_ID,
                                gestational_age = 40,
                                Instrument = "Illumina MiSeq",
                                .keep = "unused")

# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_ros_combined$read_count, n_lines$read_counts[match(metadata_ros_combined$run_accession, n_lines$sample_ID)])  # 


write.csv(metadata_ros_combined, "metadata_roswall_2021_healthy.csv")

# move sorted out samples
sort_out_samples_ros <- unique(sort_out_samples_ros)
sort_out_files_ros <- c(paste(sort_out_samples_ros, "_1.fastq.gz", sep = ""),
                    paste(sort_out_samples_ros, "_2.fastq.gz", sep = ""),
                    paste(sort_out_samples_ros, ".fastq.gz", sep = ""))
sort_out_files_ros <- sort_out_files_ros[file.exists(sort_out_files_ros)]
file.move(sort_out_files_ros, "sorted_out/")

# move good files
final_samples_ros <- metadata_ros_combined$run_accession
final_files_ros <- c(paste(final_samples_ros, "_1.fastq.gz", sep = ""),
                     paste(final_samples_ros, "_2.fastq.gz", sep = ""),
                     paste(final_samples_ros, ".fastq.gz", sep = ""))
final_files_ros <- final_files_ros[file.exists(final_files_ros)]
file.move(final_files_ros, "fastq_files/")

# clean up
rm(ena_file_roswall, metadata_ros, read_counts_ros, dict_ros_back, metadata_ros_back,
   metadata_ros_combined, sort_out_samples_ros, sort_out_files_ros, multiple_samples,
   final_files_ros, final_samples_ros, n_lines, sra_file)



# Vatanen 2018 #################################################################

setwd(paste(maindir, "16S/vatanen_2018/", sep = ""))
ena_file_vatanen <- read.table("filereport_read_run_PRJNA497734_tsv_AMPLICON.txt", sep = "\t", header = T)
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file[,c("Instrument","Run")]
all(ena_file_vatanen$run_accession %in% sra_file$Run)
ena_file_vatanen <- ena_file_vatanen[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
ena_file_vatanen <- merge(ena_file_vatanen, sra_file, by.x = "run_accession", by.y = "Run")
ena_file_vatanen$sample_title <- NULL

# read and reformat data
birth_vat <- read.table("birth.txt", sep = "\t", header = T)
colnames(birth_vat)[5] <- "birth_mode"
birth_vat$birth_mode[birth_vat$birth_mode == 0] <- "v"
birth_vat$birth_mode[birth_vat$birth_mode == 1] <- "c"
birth_vat[birth_vat == "Male"] <- "m"
birth_vat[birth_vat == "Female"] <- "f"
diabetes_vat <- read.table("diabetes.txt", sep = "\t", header = T)
feeding_vat <- read.table("feeding.txt", sep = "\t", header = T)
feeding_vat$milk_first_three_days[feeding_vat$milk_first_three_days == "mothers_breast_milk"] <- "breast"
feeding_vat$milk_first_three_days[feeding_vat$milk_first_three_days == "multiple_types_or_not_reported"] <- "unknown"
feeding_2_vat <- read.table("feeding_2.txt", sep = "\t", header = T)
feeding_2_vat$start_month <- feeding_2_vat$start_month * 30 
feeding_2_vat <- reshape(feeding_2_vat, timevar = "dietary_compound", idvar = "subjectID", direction = "wide")
growth_vat <- read.table("growth.txt", sep = "\t", header = T)
samples_vat <- read.table("samples.txt", sep = "\t", header = T)
samples_vat <- samples_vat[!samples_vat$gid_16s == "",]   # only samples with 16S-data
samples_vat$gid_wgs<- NULL

# duplicates seem not to have lower read counts, leave them for now
dupls_ena <- ena_file_vatanen[duplicated(ena_file_vatanen$sample_alias),]
hist(ena_file_vatanen$read_count[ena_file_vatanen$sample_alias %in% dupls_ena$sample_alias], breaks = 100)
hist(ena_file_vatanen$read_count[!ena_file_vatanen$sample_alias %in% dupls_ena$sample_alias], xlim = c(0,350000), breaks = 1000)
hist(ena_file_vatanen$read_count, xlim = c(0,350000), breaks = 1000)

#combine data.frames
sum(!ena_file_vatanen$sample_alias %in% samples_vat$sampleID)
sum(!samples_vat$sampleID %in% ena_file_vatanen$sample_alias)
nrow(ena_file_vatanen)
metadata_vat_combined <- merge(samples_vat, ena_file_vatanen, by.x = "gid_16s", by.y = "library_name")
nrow(ena_file_vatanen)

sum(!metadata_vat_combined$subjectID %in% birth_vat$subjectID)  # all samples have birth data
nrow(metadata_vat_combined)
metadata_vat_combined <- merge(metadata_vat_combined, birth_vat, by = "subjectID")
nrow(metadata_vat_combined)

metadata_vat_combined[which(!metadata_vat_combined$subjectID %in% diabetes_vat$subjectID),]  # information missing for one individual only
nrow(metadata_vat_combined)
metadata_vat_combined <- merge(metadata_vat_combined, diabetes_vat, all.x = T,by = "subjectID")
nrow(metadata_vat_combined)

metadata_vat_combined[(!metadata_vat_combined$subjectID %in% feeding_vat$subjectID),] # some samples lack feeding metadata
metadata_vat_combined[(!metadata_vat_combined$subjectID %in% feeding_2_vat$subjectID),]
nrow(metadata_vat_combined)
metadata_vat_combined <- merge(metadata_vat_combined, feeding_vat, all.x = T, by = "subjectID")
metadata_vat_combined <- merge(metadata_vat_combined, feeding_2_vat, all.x = T, by = "subjectID")
nrow(metadata_vat_combined)

sum(!unique(metadata_vat_combined$subjectID) %in% growth_vat$subjectID)  # some samples lack growth data
sum(!unique(growth_vat$subjectID) %in% metadata_vat_combined$subjectID)  
nrow(metadata_vat_combined)
metadata_vat_combined <- merge(metadata_vat_combined, growth_vat, all.x = T,by = "subjectID")
nrow(metadata_vat_combined)


metadata_vat_combined$food <- "unknown"
metadata_vat_combined$food[metadata_vat_combined$age_at_collection >=
                                 metadata_vat_combined$start_month.solid_start_approx] <- "solid"
metadata_vat_combined$food[metadata_vat_combined$age_at_collection <= 
                                 metadata_vat_combined$bf_length] <- "mixed"
metadata_vat_combined$food[metadata_vat_combined$age_at_collection >= 
                                 metadata_vat_combined$start_month.solid_start_approx &
                                 metadata_vat_combined$age_at_collection <= 
                                 metadata_vat_combined$bf_length] <- "mixed_b_s"
metadata_vat_combined$food[metadata_vat_combined$age_at_collection <= 
                                 metadata_vat_combined$start_month.solid_start_approx &
                                 metadata_vat_combined$age_at_collection > 
                                 metadata_vat_combined$bf_length_exclusive] <- "mixed_b_f"
metadata_vat_combined$food[metadata_vat_combined$age_at_collection <= 
                                 metadata_vat_combined$start_month.solid_start_approx &
                                 metadata_vat_combined$age_at_collection > 
                                 metadata_vat_combined$bf_length] <- "formula"
metadata_vat_combined$food[metadata_vat_combined$age_at_collection <= 
                                 metadata_vat_combined$bf_length_exclusive] <- "breast"

# remove some columns
metadata_vat_combined[c("start_month.apple", "start_month.banana", "start_month.barley",
                        "start_month.beef", "start_month.berries",  "start_month.cabbs",
                        "start_month.carrot", "start_month.corn", "start_month.cowsmilk",
                        "start_month.egg", "start_month.fish",
                        "start_month.icecream", "start_month.milkprod", "start_month.oat",                    
                        "start_month.pear", "start_month.peas", "start_month.plum",                     
                        "start_month.pork", "start_month.potato", "start_month.poultry",                  
                        "start_month.rice", "start_month.rye", "start_month.spinach",                  
                        "start_month.tomato", "start_month.wheat", "start_month.mango",                    
                        "start_month.beetroot", "start_month.sweetpota", "start_month.soymilk",                  
                        "start_month.turnip")] <- NULL

# add info on antibiotic usage:
abx_data <- read.table("abx_data_1.txt", sep = "\t", header = T) %>%
  mutate(first_antibiotics = as.numeric(Age..months.) * 30)
abx_data <- aggregate(first_antibiotics ~ Subject, abx_data, function(x) min(x))

metadata_vat_combined <- left_join(metadata_vat_combined, abx_data, by = join_by(subjectID == Subject))

abx_data_2 <- readxl::read_xlsx("1-s2.0-S1931312815000219-mmc2.xlsx")
abx_data_2[,2:86] <- NULL
metadata_vat_combined <- left_join(metadata_vat_combined, abx_data_2,
                                   by = join_by(subjectID == "ID, E=Espoo, Finland, T=Tartu, Estonia"))
abx_data_3 <- readxl::read_xlsx("1-s2.0-S0092867416303981-mmc2.xlsx", sheet = "Basics and medications") %>%
  tibble::column_to_rownames(var = "Participant")
abx_data_3 <- abx_data_3[,c(grep("_Age", colnames(abx_data_3)))]
abx_data_3 <- data.frame(first_abx_2 = apply(abx_data_3, 1, function(x) min(x, na.rm = T))) %>%
  mutate(first_abx_2 = as.numeric(first_abx_2) * 30)

abx_data_3$subject <- rownames(abx_data_3)
metadata_vat_combined <- left_join(metadata_vat_combined, abx_data_3,
                                   by = join_by(subjectID == subject))


metadata_vat_combined <- metadata_vat_combined %>%
  mutate(antibiotics_before = case_when(#first_antibiotics < age_at_collection ~ T, # 467
                                        !is.na(`First antibiotic drug between 3 and 6 months`) & (age_at_collection > 90) ~ T,
                                        !is.na(`First antibiotic drug between 6 and 12 months`) & age_at_collection > 180 ~ T,
                                        !is.na(`First antibiotic drug between 12 and 18 months`) & age_at_collection > 365 ~ T,
                                        !is.na(`First antibiotic drug between 18 and 24 months`) & age_at_collection > 540 ~ T,
                                        !is.na(`First antibiotic drug between 24 and 36 months`) & age_at_collection > 730 ~ T,
                                        first_abx_2 < age_at_collection ~ T, # 529
                                        .default = F)
  )

# remove individuals with one sample only
multiple_samples <- unique(metadata_vat_combined$subjectID[duplicated(metadata_vat_combined$subjectID)])
# sort_out_samples_vat <- metadata_vat_combined$run_accession[!metadata_vat_combined$subjectID %in% multiple_samples]
# metadata_vat_combined <- metadata_vat_combined[metadata_vat_combined$subjectID %in% multiple_samples,]

# retreive sorted out samples
metadata_vat_combined$run_accession[!metadata_vat_combined$subjectID %in% multiple_samples] %>%
  paste("fastq-dump --split-e --gzip ",.) %>%
  data.frame(`#!/bin/bash` = .) %>%
  write.table(., file = "/fast/AG_Forslund/rob/studies/16S/vatanen_2018/pull_single.sh",
              sep = " ", quote = F, row.names = F, col.names = "#!/bin/bash")


# unify colnames
metadata_vat_combined %<>% mutate(sample_ID = paste0(sampleID, "_vatanen_2018"),
                                  subject_ID = paste0(subjectID, "_vatanen_2018"),
                                  age = age_at_collection,
                                  birthmode = birth_mode,
                                  sex = gender,
                                  exclusive_breastfeeding = bf_length_exclusive,
                                  study = "vatanen_2018",
                                  lifestyle = "industrialized",
                                  family_ID = subject_ID,
                                  country = case_when(country == "FIN" ~ "FINNLAND",
                                                      country == "RUS" ~ "RUSSIA",
                                                      country == "EST" ~ "ESTONIA"),
                                  region = case_when(country == "FINNLAND" ~ "FINNLAND",
                                                      country == "RUSSIA" ~ "KARELIA",
                                                      country == "ESTONIA" ~ "ESTONIA"),
                                  geographic_location_.latitude. = case_when(country == "FINNLAND" ~ 60.192059,
                                                     country == "RUSSIA" ~ 61.78491,
                                                     country == "ESTONIA" ~ 58.378025),
                                  geographic_location_.longitude. = case_when(country == "FINNLAND" ~ 24.945831,
                                                     country == "RUSSIA" ~ 34.34691,
                                                     country == "ESTONIA" ~ 26.728493),
                                  .keep = "unused")

metadata_vat_combined$region <- metadata_vat_combined$country
metadata_vat_combined$region[metadata_vat_combined$country == "RUSSIA"] <- "KARELIA"

sort(colnames(metadata_vat_combined))
write.csv(metadata_vat_combined, "metadata_vat_2018_healthy.csv")

# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_vat_combined$read_count, n_lines$read_counts[match(metadata_vat_combined$run_accession, n_lines$sample_ID)])  # 


# sort_out_samples_vat <- unique(sort_out_samples_vat)
# sort_out_files_vat <- c(paste(sort_out_samples_vat, "_1.fastq.gz", sep = ""),
#                         paste(sort_out_samples_vat, "_2.fastq.gz", sep = ""),
#                         paste(sort_out_samples_vat, ".fastq.gz", sep = ""))
# sort_out_files_vat <- sort_out_files_vat[file.exists(sort_out_files_vat)]
# file.move(sort_out_files_vat, "sorted_out/")

# move good files
final_samples_vat <- metadata_vat_combined$run_accession
final_files_vat <- c(paste(final_samples_vat, "_1.fastq.gz", sep = ""),
                     paste(final_samples_vat, "_2.fastq.gz", sep = ""),
                     paste(final_samples_vat, ".fastq.gz", sep = ""))
final_files_vat <- final_files_vat[file.exists(final_files_vat)]
file.move(final_files_vat, "fastq_files/")

rm(ena_file_vatanen, ena_file_vat, birth_vat, diabetes_vat, feeding_2_vat, feeding_vat, growth_vat,
   samples_vat, metadata_vat_combined, dupls_ena, multiple_samples, sort_out_samples_vat,
   sort_out_files_vat, final_samples_vat, final_files_vat, n_lines, abx_data, abx_data_2,
   abx_data_3, sra_file)




# Chu 2017 #####################################################################
# maybe also integrate _pools information? need to know which timepoint is which.
# metadata not good, better leave it out
# 
# setwd(paste(maindir, "16S/chu_2017/", sep = ""))
# ena_file_chu <- read.table("filereport_read_run_PRJNA322188_tsv_stool_only_16S.txt", sep = "\t", header = T)
# ena_file_chu <- ena_file_chu[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
# ena_file_chu$child_ID <- gsub(".Stool.*", "", ena_file_chu$sample_title)
# ena_file_chu$child_ID <- gsub("...STOOL.*", "", ena_file_chu$child_ID)
# ena_file_chu$child_ID <- gsub("S.", "S", ena_file_chu$child_ID)
# ena_file_chu$time_point <- NA
# ena_file_chu$time_point[grep("STOOL1", ena_file_chu$sample_title)] <- 1
# ena_file_chu$time_point[grep("STOOL2", ena_file_chu$sample_title)] <- 36
# ena_file_chu$sample_ID <- paste(ena_file_chu$child_ID, ena_file_chu$time_point, sep = "_")
# # remove unclear time points
# sort_out_samples_chu <- ena_file_chu$run_accession[grep("STOOL3", ena_file_chu$sample_title)]
# ena_file_chu <- ena_file_chu[-grep("STOOL3", ena_file_chu$sample_title),]
# sort_out_samples_chu <- c(sort_out_samples_chu, ena_file_chu$run_accession[grep("STOOL4", ena_file_chu$sample_title)])
# ena_file_chu <- ena_file_chu[-grep("STOOL4", ena_file_chu$sample_title),]
# # remove mother samples and all individuals with one sample only
# sort_out_samples_chu <- c(sort_out_samples_chu, ena_file_chu$run_accession[grep("maternal", ena_file_chu$sample_title)])
# ena_file_chu <- ena_file_chu[-grep("maternal", ena_file_chu$sample_title),]
# sort_out_samples_chu <- c(sort_out_samples_chu, ena_file_chu$run_accession[grep("[.]M[.]STOOL", ena_file_chu$sample_title)])
# ena_file_chu <- ena_file_chu[-grep("[.]M[.]STOOL", ena_file_chu$sample_title),]
# # only keep samples with at least 2 time points
# multiple_samples <- unique(ena_file_chu$child_ID[duplicated(ena_file_chu$child_ID)])
# sort_out_samples_chu <- c(sort_out_samples_chu, ena_file_chu$run_accession[!ena_file_chu$child_ID %in% multiple_samples])
# ena_file_chu <- ena_file_chu[ena_file_chu$child_ID %in% multiple_samples,]
# 
# 
# 
# metadata_chu_child <- read.table("infant_metadata.txt", sep = "\t", header = T)
# metadata_chu_child$Mode.of.Delivery[metadata_chu_child$Mode.of.Delivery == "Cesarean"] <- "c"
# metadata_chu_child$Mode.of.Delivery[metadata_chu_child$Mode.of.Delivery == "Vaginal"] <- "v"
# metadata_chu_mother <- read.table("mother_metadata.txt", sep = "\t", header = T)
# metadata_chu_mother$Mode.of.Delivery <- NULL
# metadata_chu_mother$Term.v..Preterm <- NULL
# metadata_chu_mother$Study <- NULL
# metadata_chu_stats <- read.table("mapping_stats.txt", sep = "\t", header = T)
# metadata_chu_stats <- metadata_chu_stats[metadata_chu_stats$Site == "Stool",]
# metadata_chu_stats <- metadata_chu_stats[grep("maternal", metadata_chu_stats$ID ),]
# metadata_chu_stats$Subject.ID <- gsub("[.]", "", metadata_chu_stats$Subject.ID)
# metadata_chu_stats$Visit[metadata_chu_stats$Visit == "Delivery"] <- "STOOL1"
# metadata_chu_stats$Visit[metadata_chu_stats$Visit == "Postpartum"] <- "STOOL2"
# metadata_chu_stats$sample_ID <- paste(metadata_chu_stats$Subject.ID, metadata_chu_stats$Visit, sep = "_")
# metadata_chu_stats$ID <- NULL
# metadata_chu_pools <- read.table("sequencing_pools.txt", sep = "\t", header = T, comment.char = "")
# metadata_chu_pools <- metadata_chu_pools[metadata_chu_pools$Body.Site == "Stool",]
# metadata_chu_pools <- metadata_chu_pools[metadata_chu_pools$Mother == "Infant",]
# metadata_chu_pools <- metadata_chu_pools[metadata_chu_pools$Cohort == "Longitudinal",]
# 
# sum(duplicated(ena_file_chu$sample_alias)) # no duplicates!
# 
# #combine data.frames
# sum(!unique(ena_file_chu$child_ID) %in% metadata_chu_child$ID)
# sum(!unique(metadata_chu_child$ID) %in% ena_file_chu$child_ID)
# nrow(ena_file_chu)
# metadata_chu_combined <- merge(metadata_chu_child, ena_file_chu, by.x = "ID", by.y = "child_ID")
# nrow(metadata_chu_combined)
# 
# sum(!unique(ena_file_chu$child_ID) %in% metadata_chu_mother$ID)
# nrow(ena_file_chu)
# metadata_chu_combined <- merge(metadata_chu_combined, metadata_chu_mother, by = "ID")
# nrow(metadata_chu_combined)
# 
# sum(!unique(ena_file_chu$child_ID) %in% metadata_chu_mother$ID)
# nrow(ena_file_chu)
# metadata_chu_combined <- merge(metadata_chu_combined, metadata_chu_stats, by = "sample_ID", all.x = T)
# nrow(metadata_chu_combined)
# 
# metadata_chu_combined[metadata_chu_combined == "Yes"] <- "y"
# metadata_chu_combined[metadata_chu_combined == "No"] <- "n"
# metadata_chu_combined$Gender[metadata_chu_combined$Gender == "Female"] <- "f"
# metadata_chu_combined$Gender[metadata_chu_combined$Gender == "Male"] <- "m"
# metadata_chu_combined$antibiotics_any[metadata_chu_combined$Postnatal.Antibiotics == "None"] <- "n"
# metadata_chu_combined$antibiotics_any[metadata_chu_combined$Postnatal.Antibiotics != "None"] <- "y"
# metadata_chu_combined$Postnatal.Antibiotics <- NULL
# # missing timepoint information is from time point 2
# metadata_chu_combined$time_point[is.na(metadata_chu_combined$time_point)] <- 36
# 
# 
# # filter out Preterm samples:
# sort_out_samples_chu <- c(sort_out_samples_chu, metadata_chu_combined$run_accession[metadata_chu_combined$Term.v..Preterm == "Preterm"])
# metadata_chu_combined <- metadata_chu_combined[!metadata_chu_combined$Term.v..Preterm == "Preterm",]
# 
# # move sorted out samples
# sort_out_samples_chu <- unique(sort_out_samples_chu)
# sort_out_files_chu <- c(paste(sort_out_samples_chu, "_1.fastq.gz", sep = ""),
#                         paste(sort_out_samples_chu, "_2.fastq.gz", sep = ""),
#                         paste(sort_out_samples_chu, ".fastq.gz", sep = ""))
# sort_out_files_chu <- sort_out_files_chu[file.exists(sort_out_files_chu)]
# file.move(sort_out_files_chu, "sorted_out/")
# 
# 
# # clean up cols
# metadata_chu_combined[c("Study", "Term.v..Preterm", "Gestational.Age...Delivery",
#                         "Multiples", "Subject.ID", "Site", "Visit", "Total.Paired.Reads",
#                         "Paired.Reads.after.Human.Filtering", "X...", "HUMAnN.Mapped",
#                         "X....1")] <- NULL
# 
# # unify colnames
# colnames(metadata_chu_combined)[match(c("ID", "Gender", "Gestational.Age..Weeks.Days.",
#                                         "Birth.Weight..g.", "Mode.of.Delivery"),
#                                       colnames(metadata_chu_combined))] <- c("subject_ID", "sex", "gestational_age",
#                                                                              "birth_weight", "birthmode")
# 
# metadata_chu_combined$country <- "USA"
# metadata_chu_combined$region <- "TEXAS"
# metadata_chu_combined$study <- "chu_2017"
# 
# 
# 
# # write file
# write.csv(metadata_chu_combined, "metadata_chu_2017_healthy.csv")
# 
# rm(ena_file_chu, multiple_samples, sort_out_samples_chu, metadata_chu_child,
#    metadata_chu_mother, metadata_chu_pools, metadata_chu_stats, metadata_chu_combined, sort_out_files_chu)
# 



# Reyman 2019 (waiting for metadata) ###########################################

setwd(paste(maindir, "16S/reyman_2019/", sep = ""))
ena_file_reyman <- read.table("filereport_read_run_PRJNA481243_tsv.txt", sep = "\t", header = T)
ena_file_reyman <- ena_file_reyman[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]

metadata_rey <- read.table("supplement_1.txt", sep = '\t', header = T)
metadata_rey_sra <- read.table("SraRunTable.txt", sep = ',', header = T) %>%
  mutate(run_accession = Run,
         subject_ID = PID,
         sample_ID = as.character(Sample.Name),
         age = case_when(Time_point == 40 ~ 0,
                         Time_point == 41 ~ 1,
                         Time_point == 42 ~ 7,
                         Time_point == 43 ~ 14,
                         Time_point == 44 ~ 30,
                         Time_point == 45 ~ 60,
                         Time_point == 47 ~ 120,
                         Time_point == 48 ~ 180,
                         Time_point == 49 ~ 270,
                         Time_point == 50 ~ 365,
                         Time_point == 60 ~ NA),
         country = "NETHERLANDS",
         region = "UTRECHT",
         study = "reyman_2019",
         lifestyle = "industrialized",
         family_ID = subject_ID,
         geographic_location_.latitude. = 52.306100,
         geographic_location_.longitude. = 4.690700,
         .keep = "unused")

metadata_rey_sra[,c("Assay.Type", "BioSampleModel", "Bytes", "Center.Name", "Collection_Date",
                    "Consent", "DATASTORE.filetype", "DATASTORE.provider", "DATASTORE.region",
                    "env_biome", "geo_loc_name_country", "geo_loc_name_country_continent",
                    "geo_loc_name", "HOST", "LibraryLayout", "Lat_Lon",
                    "LibrarySelection", "LibrarySource", "Organism", "Platform", "sample_names")] <- NULL
# remove blank samples
sort_out_samples_reyman <- metadata_rey_sra$Run[metadata_rey_sra$Isolation_Source == "Blank"]
metadata_rey_sra <- metadata_rey_sra[metadata_rey_sra$Isolation_Source != "Blank",]


sra_file_bog <- read.table("SraRunTable_bogaert_2023.txt", sep = ",", header = T, comment.char = "", colClasses = "character")
sra_file_bog[,c("Assay.Type", "AvgSpotLen", "BioSampleModel", "Bytes", "Center.Name",
                "Collection_Date", "Consent", "DATASTORE.filetype", "DATASTORE.provider",
                "DATASTORE.region", "geo_loc_name_country", "geo_loc_name_country_continent",
                "geo_loc_name", "HOST", "Instrument", "Lat_Lon", "Library.Name", "LibraryLayout",
                "LibrarySelection", "Organism", "Platform", "ReleaseDate", "create_date",
                "version", "sample_id..run.", "SRA.Study", "LibrarySource",
                "niche..run.", "niche", "Bases", "BioProject", "BioSample", "Experiment",
                "miseq_run_nr..run.", "miseq_run_nr", "Sample.Name", "run_accession")] <- NULL
sra_file_bog <- sra_file_bog %>%
  filter(Isolation_Source == "fae",
         !is.na(TIME)) %>%
  mutate(run_accession = Run,
         birthmode = ifelse(Birth_mode == "vag", yes = "v", no = "c"),
         sex = ifelse(gender == "female", yes = "f", no = "m"),
         subject_ID = gsub("M", "", Subject) %>% as.numeric(),
         sample_ID = Sample_id,
         # age = case_when(TIME == "d0" ~ 0,
         #                 TIME == "d1" ~ 1,
         #                 TIME == "w1" ~ 7,
         #                 TIME == "w2" ~ 14,
         #                 TIME == "m1" ~ 30),
         .keep = "unused") %>% filter(!duplicated(subject_ID))
sra_file_bog[,c("Isolation_Source", "run_accession", "sample_ID", "TIME")] <- NULL


# no duplicates
sum(duplicated(metadata_rey_sra$Sample.Name))

metadata_rey_sra <- left_join(metadata_rey_sra, sra_file_bog, by = "subject_ID")

# remove unknown age samples:
metadata_rey_sra %<>% filter(!is.na(age))

# no individuals with only one sample
multiple_samples <- unique(metadata_rey_sra$PID[duplicated(metadata_rey_sra$PID)])
length(metadata_rey_sra$Run[!metadata_rey_sra$PID %in% multiple_samples])

metadata_rey_sra %<>% mutate(sample_ID = paste0(sample_ID, "_reyman_2019"),
                             subject_ID = paste0(subject_ID, "_reyman_2019"),
                             .keep = "unused")


# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_rey_sra$Bases, n_lines$read_counts[match(metadata_rey_sra$run_accession, n_lines$sample_ID)])  # 


write.csv(metadata_rey_sra, "metadata_rey_2019_healthy.csv")

# move sorted out samples
sort_out_samples_reyman <- unique(ena_file_reyman$run_accession[!ena_file_reyman$run_accession %in% metadata_rey_sra$run_accession])
sort_out_files_reyman <- c(paste("fastq_files/", sort_out_samples_reyman, "_1.fastq.gz", sep = ""),
                           paste("fastq_files/", sort_out_samples_reyman, "_2.fastq.gz", sep = ""),
                           paste("fastq_files/", sort_out_samples_reyman, ".fastq.gz", sep = ""))
sort_out_files_reyman <- sort_out_files_reyman[file.exists(sort_out_files_reyman)]
file.remove(sort_out_files_reyman)


# move good files
final_samples_reyman <- metadata_rey_sra$run_accession
final_files_reyman <- c(paste(final_samples_reyman, "_1.fastq.gz", sep = ""),
                       paste(final_samples_reyman, "_2.fastq.gz", sep = ""),
                       paste(final_samples_reyman, ".fastq.gz", sep = ""))
final_files_reyman <- final_files_reyman[file.exists(final_files_reyman)]
file.move(final_files_reyman, "fastq_files/")


# clean up
rm(ena_file_reyman, metadata_rey, metadata_rey_sra, final_files_reyman, final_samples_reyman,
   sort_out_files_reyman, sort_out_samples_reyman, multiple_samples, n_lines, sra_file_bog, sra_file)



# Wampach 2018 #################################################################

setwd(paste(maindir, "16S/wampach_2018/", sep = ""))
ena_file_wampach <- read.table("COSMIC_metadata_NCBI_16.txt", sep = ",", header = T)
ena_file_wampach <- ena_file_wampach[c("Run", "AvgSpotLen", "Bases", "BioSample",
                                       "Collection_date", "env_biome",
                                       "env_material", "Library.Name", "Sample.Name",
                                       "Host_Age")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file[,c("Instrument","Run")]
all(ena_file_wampach$run_accession %in% sra_file$Run)
ena_file_wampach <- merge(ena_file_wampach, sra_file,by = "Run")

# get time-points and individual ID:
ena_file_wampach <- tidyr::separate(ena_file_wampach, col = "Sample.Name", into = c("time_point", "subject_ID", NA), sep = "[_]", remove = F)
ena_file_wampach$age <- NA
ena_file_wampach$age[ena_file_wampach$time_point == "V1"] <- 1
ena_file_wampach$age[ena_file_wampach$time_point == "V2"] <- 3
ena_file_wampach$age[ena_file_wampach$time_point == "V3"] <- 5
ena_file_wampach$age[ena_file_wampach$time_point == "V4"] <- 28
ena_file_wampach$age[ena_file_wampach$time_point == "V5"] <- 150
ena_file_wampach$age[ena_file_wampach$time_point == "V6"] <- 365
# remove vaginal and control samples and 18S-seq samples
sort_out_samples_wam <- ena_file_wampach$Run[!ena_file_wampach$env_biome == "human gut"]
ena_file_wampach <- ena_file_wampach[ena_file_wampach$env_biome == "human gut",]
sort_out_samples_wam <- c(sort_out_samples_wam, ena_file_wampach$Run[ena_file_wampach$Sample.Name == "Standard_sample"])
ena_file_wampach <- ena_file_wampach[!ena_file_wampach$Sample.Name == "Standard_sample",]
sort_out_samples_wam <- c(sort_out_samples_wam, ena_file_wampach$Run[ena_file_wampach$time_point == "M"])
ena_file_wampach <- ena_file_wampach[!ena_file_wampach$time_point == "M",]
sort_out_samples_wam <- c(sort_out_samples_wam, ena_file_wampach$Run[grep("18S", ena_file_wampach$Library.Name)])
ena_file_wampach <- ena_file_wampach[-grep("18S", ena_file_wampach$Library.Name),]

# remove all individuals with only one sample
multiple_samples <- unique(ena_file_wampach$subject_ID[duplicated(ena_file_wampach$subject_ID)])
# sort_out_samples_wam <- c(sort_out_samples_wam, ena_file_wampach$Run[!ena_file_wampach$subject_ID %in% multiple_samples])
# ena_file_wampach <- ena_file_wampach[ena_file_wampach$subject_ID %in% multiple_samples,]
# retreive sorted out samples
ena_file_wampach$Run[!ena_file_wampach$subject_ID %in% multiple_samples] %>%
  paste("fastq-dump --split-e --gzip ",.) %>%
  # data.frame(`#!/bin/bash` = .) %>%
  write.table(., file = "/fast/AG_Forslund/rob/studies/16S/wampach_2018//pull_single.sh",
              sep = " ", quote = F, row.names = F, col.names = "#!/bin/bash")



metadata_wam_1_init <- read.table("suppl_1_init.txt", sep = '\t', header = T)
metadata_wam_1_init$Delivery.mode[metadata_wam_1_init$Delivery.mode == "Vaginal"] <- "v"
metadata_wam_1_init$Delivery.mode[grep("section", metadata_wam_1_init$Delivery.mode)] <- "c"
metadata_wam_1_init$Gender[metadata_wam_1_init$Gender == "Male"] <- "m"
metadata_wam_1_init$Gender[metadata_wam_1_init$Gender == "Female"] <- "f"
metadata_wam_1_init <- separate(metadata_wam_1_init, col = "Gestational.age....weeks.",
                                into = c("gestational_week", "gestational_day"),
                                sep = "[+]", remove = F, convert = T)
metadata_wam_1_init$gestational_day <- gsub("day.*", "", metadata_wam_1_init$gestational_day) %>% as.numeric() /7
metadata_wam_1_init$gestational_age <- metadata_wam_1_init$gestational_week + metadata_wam_1_init$gestational_day
metadata_wam_1_init[c("gestational_day", "gestational_week")] <- NULL


metadata_wam_1_d1 <- read.table("suppl_1_d1.txt", sep = '\t', header = T)
metadata_wam_1_d1$Type.of.Milk <- gsub("Breast", "breast", metadata_wam_1_d1$Type.of.Milk)
metadata_wam_1_d1$Type.of.Milk <- gsub("Formula", "formula", metadata_wam_1_d1$Type.of.Milk)
metadata_wam_1_d1$Type.of.Milk <- gsub("Combined", "mixed_b_f", metadata_wam_1_d1$Type.of.Milk)
metadata_wam_1_d1$Type.of.Milk <- gsub("not yet fed", "nothing", metadata_wam_1_d1$Type.of.Milk)

metadata_wam_1_d3 <- read.table("suppl_1_d3.txt", sep = '\t', header = T)
metadata_wam_1_d3$Type.of.Milk <- gsub("Breast", "breast", metadata_wam_1_d3$Type.of.Milk)
metadata_wam_1_d3$Type.of.Milk <- gsub("Formula", "formula", metadata_wam_1_d3$Type.of.Milk)
metadata_wam_1_d3$Type.of.Milk <- gsub("Combined", "mixed_b_f", metadata_wam_1_d3$Type.of.Milk)

metadata_wam_1_d5 <- read.table("suppl_1_d5.txt", sep = '\t', header = T)
metadata_wam_1_d5$Type.of.Milk <- gsub("Breast", "breast", metadata_wam_1_d5$Type.of.Milk)
metadata_wam_1_d5$Type.of.Milk <- gsub("Formula", "formula", metadata_wam_1_d5$Type.of.Milk)
metadata_wam_1_d5$Type.of.Milk <- gsub("Combined", "mixed_b_f", metadata_wam_1_d5$Type.of.Milk)

metadata_wam_1_mother <- read.table("suppl_1_init_mother.txt", sep = '\t', header = T)
metadata_wam_1_mother[c("Collected.samples", "Data.available", "Mother.age.range", "Mother.Ethnicity")] <- NULL

metadata_wam_4_d5 <- read.table("suppl_4_d5.txt", sep = '\t', header = T)

metadata_wam_4_m1 <- read.table("suppl_4_1m.txt", sep = '\t', header = T)
metadata_wam_4_m1$Feeding..milk. <- gsub("Breast", "breast", metadata_wam_4_m1$Feeding..milk.)
metadata_wam_4_m1$Feeding..milk. <- gsub("Formula", "formula", metadata_wam_4_m1$Feeding..milk.)
metadata_wam_4_m1$Feeding..milk. <- gsub("Mixed", "mixed_b_f", metadata_wam_4_m1$Feeding..milk.)

metadata_wam_4_m6 <- read.table("suppl_4_6m.txt", sep = '\t', header = T)
metadata_wam_4_m6$Feeding..milk. <- gsub("Breast, Mixed", "mixed_b_s", metadata_wam_4_m6$Feeding..milk.)
metadata_wam_4_m6$Feeding..milk. <- gsub("Breast", "breast", metadata_wam_4_m6$Feeding..milk.)
metadata_wam_4_m6$Feeding..milk. <- gsub("Formula", "formula", metadata_wam_4_m6$Feeding..milk.)
metadata_wam_4_m6$Feeding..milk. <- gsub("Mixed", "mixed_b_f", metadata_wam_4_m6$Feeding..milk.)

metadata_wam_4_y1 <- read.table("suppl_4_1y.txt", sep = '\t', header = T)
metadata_wam_4_y1$Feeding..milk. <- gsub("Breast, Formula, Mixed", "mixed_b_f_s", metadata_wam_4_y1$Feeding..milk.)
metadata_wam_4_y1$Feeding..milk. <- gsub("Breast, Mixed", "mixed_b_s", metadata_wam_4_y1$Feeding..milk.)
metadata_wam_4_y1$Feeding..milk. <- gsub("Formula, Mixed", "mixed_f_s", metadata_wam_4_y1$Feeding..milk.)
metadata_wam_4_y1$Feeding..milk. <- gsub("Breast", "breast", metadata_wam_4_y1$Feeding..milk.)
metadata_wam_4_y1$Feeding..milk. <- gsub("Formula", "formula", metadata_wam_4_y1$Feeding..milk.)


# no duplicates
sum(duplicated(ena_file_wampach$Sample.Name))

#combine data.frames
sum(!unique(ena_file_wampach$child_id) %in% metadata_wam_1_init$Study.ID)
sum(!unique(metadata_wam_1_init$Study.ID) %in% ena_file_wampach$child_id)
nrow(ena_file_wampach)
metadata_wam_combined <- merge(metadata_wam_1_init, ena_file_wampach, by.x = "Study.ID", by.y = "subject_ID")
nrow(metadata_wam_combined)
metadata_wam_combined$weight[metadata_wam_combined$time_point == "V1"] <- 
   metadata_wam_combined$Birth.weight...g.[metadata_wam_combined$time_point == "V1"]



sum(!unique(metadata_wam_combined$Study.ID) %in% metadata_wam_1_mother$Study.ID)
sum(!unique(metadata_wam_1_mother$Study.ID) %in% metadata_wam_combined$Study.ID)
nrow(metadata_wam_combined)
metadata_wam_combined <- merge(metadata_wam_combined, metadata_wam_1_mother, by = "Study.ID")
nrow(metadata_wam_combined)

# add feeding and weight for each time point
metadata_wam_combined$food <- NA
metadata_wam_combined$recent_diarrhoea <- NA
metadata_wam_combined$antibiotics_before <- NA
metadata_wam_combined$antibiotics_before[metadata_wam_combined$time_point %in% c("V1", "V2", "V3")] <- "n"
   
metadata_wam_combined$food[metadata_wam_combined$time_point == "V1"] <- 
   metadata_wam_1_d1$Type.of.Milk[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V1"], metadata_wam_1_d1$Study.ID)]

metadata_wam_combined$food[metadata_wam_combined$time_point == "V2"] <- 
   metadata_wam_1_d3$Type.of.Milk[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V2"], metadata_wam_1_d3$Study.ID)]
metadata_wam_combined$weight[metadata_wam_combined$time_point == "V2"] <- 
   metadata_wam_1_d3$Weight...g.[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V2"], metadata_wam_1_d3$Study.ID)]

metadata_wam_combined$food[metadata_wam_combined$time_point == "V3"] <- 
   metadata_wam_1_d5$Type.of.Milk[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V3"], metadata_wam_1_d5$Study.ID)]
metadata_wam_combined$weight[metadata_wam_combined$time_point == "V3"] <- 
   metadata_wam_1_d5$Weight...g.[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V3"], metadata_wam_1_d5$Study.ID)]
metadata_wam_combined$recent_diarrhoea[metadata_wam_combined$time_point == "V3"] <- 
   metadata_wam_4_d5$Recent.Diarrhoea[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V3"], metadata_wam_4_d5$Study.ID)]

metadata_wam_combined$food[metadata_wam_combined$time_point == "V4"] <- 
   metadata_wam_4_m1$Feeding..milk.[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V4"], metadata_wam_4_m1$Study.ID)]
metadata_wam_combined$weight[metadata_wam_combined$time_point == "V4"] <- 
   metadata_wam_4_m1$Weight..g.[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V4"], metadata_wam_4_m1$Study.ID)]
metadata_wam_combined$recent_diarrhoea[metadata_wam_combined$time_point == "V4"] <- 
   metadata_wam_4_m1$Baby.diarrhoea[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V4"], metadata_wam_4_m1$Study.ID)]
metadata_wam_combined$antibiotics_before[metadata_wam_combined$time_point == "V4"] <- 
   metadata_wam_4_m1$Antibiotics.since.last.visit.[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V4"], metadata_wam_4_m1$Study.ID)]


metadata_wam_combined$food[metadata_wam_combined$time_point == "V5"] <- 
   metadata_wam_4_m6$Feeding..milk.[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V5"], metadata_wam_4_m6$Study.ID)]
metadata_wam_combined$weight[metadata_wam_combined$time_point == "V5"] <- 
   metadata_wam_4_m6$Weight..g.[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V5"], metadata_wam_4_m6$Study.ID)]
metadata_wam_combined$recent_diarrhoea[metadata_wam_combined$time_point == "V5"] <- 
   metadata_wam_4_m6$Recent.baby.Diarrhoea[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V5"], metadata_wam_4_m6$Study.ID)]
metadata_wam_combined$antibiotics_before[metadata_wam_combined$time_point == "V5"] <- 
   metadata_wam_4_m6$Antibiotics.since.last.visit.[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V5"], metadata_wam_4_m1$Study.ID)]


metadata_wam_combined$food[metadata_wam_combined$time_point == "V6"] <- 
   metadata_wam_4_y1$Feeding..milk.[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V6"], metadata_wam_4_y1$Study.ID)]
metadata_wam_combined$weight[metadata_wam_combined$time_point == "V6"] <- 
   metadata_wam_4_y1$Weight..g.[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V6"], metadata_wam_4_y1$Study.ID)]
metadata_wam_combined$recent_diarrhoea[metadata_wam_combined$time_point == "V6"] <- 
   metadata_wam_4_y1$Recent.baby.Diarrhoea[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V6"], metadata_wam_4_y1$Study.ID)]
metadata_wam_combined$antibiotics_before[metadata_wam_combined$time_point == "V6"] <- 
   metadata_wam_4_y1$Antibiotics.since.last.visit.[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V6"], metadata_wam_4_m1$Study.ID)]

metadata_wam_combined$antibiotics_before[grep("Yes", metadata_wam_combined$antibiotics_before)] <- "y"


metadata_wam_combined[metadata_wam_combined == "Yes"] <- T
metadata_wam_combined[metadata_wam_combined == "y"] <- T
metadata_wam_combined[metadata_wam_combined == "No"] <- F
metadata_wam_combined[metadata_wam_combined == "No "] <- F
metadata_wam_combined[metadata_wam_combined == "n"] <- F
metadata_wam_combined <- type.convert(metadata_wam_combined)
# remove unneeded columns
metadata_wam_combined[,c("Study.Allocation", "Cohort", "Gestational.age.category", "Small.for..gestational.age..SGA.",
                         "AvgSpotLen", "Collection_date", "Strep.positive..antenatal.screen",
                         "Gestational.age....weeks.", "Host_Age")] <- NULL


# unify colnames
metadata_wam_combined %<>% mutate(run_accession = Run,
                             sample_ID = paste0(Sample.Name, "_wampach_2018"),
                             subject_ID = paste0(Study.ID, "_wampach_2018"),
                             birthmode = Delivery.mode,
                             birth_weight = Birth.weight...g. / 1000,
                             sex = Gender,
                             birth_height = Length...cm.,
                             country = "LUXEMBOURG",
                             region = "LUXEMBOURG",
                             study = "wampach_2018",
                             lifestyle = "industrialized",
                             family_ID = subject_ID,
                             weight = weight / 1000,
                             geographic_location_.latitude. = 49.611622,
                             geographic_location_.longitude. = 6.131935,
                             .keep = "unused")


# move sorted out samples
sort_out_samples_wam <- unique(sort_out_samples_wam)
sort_out_files_wam <- c(paste(sort_out_samples_wam, "_1.fastq.gz", sep = ""),
                        paste(sort_out_samples_wam, "_2.fastq.gz", sep = ""),
                        paste(sort_out_samples_wam, ".fastq.gz", sep = ""))
sort_out_files_wam <- sort_out_files_wam[file.exists(sort_out_files_wam)]
file.move(sort_out_files_wam, "sorted_out/")

# move good files
final_samples_wam <- metadata_wam_combined$run_accession
final_files_wam <- c(paste(final_samples_wam, "_1.fastq.gz", sep = ""),
                     paste(final_samples_wam, "_2.fastq.gz", sep = ""),
                     paste(final_samples_wam, ".fastq.gz", sep = ""))
final_files_wam <- final_files_wam[file.exists(final_files_wam)]
file.move(final_files_wam, "fastq_files/")


# check read numbers
# n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
# plot(metadata_wam_combined$Bases, n_lines$read_counts[match(metadata_wam_combined$run_accession, n_lines$sample_ID)])  # 
# metadata_wam_combined$read_count <- n_lines$read_counts


# write file
write.csv(metadata_wam_combined, "metadata_wam_2018_healthy.csv")


# clean up
rm(ena_file_wampach, metadata_wam_1_d1, metadata_wam_1_d3, metadata_wam_1_d5, metadata_wam_1_init,
   metadata_wam_1_mother, metadata_wam_4_d5, metadata_wam_4_m1, metadata_wam_4_m6, metadata_wam_4_y1,
   metadata_wam_combined, sort_out_files_wam, sort_out_samples_wam, multiple_samples,
   final_files_wam, final_samples_wam, n_lines, sra_file)





# 
# # Bockulich 2016 #############################################################
# # metadata confusing, need to wait
# # better to use information from SRA, also aligns better with sample ids
# # R-rectal swab, V-vaginal swab, SS-stool sample, SD-dry stool
# # not 100% sure for some metadata
# # what does sampletype "Repeats" mean? maybe sort them out?
# 
# setwd(paste(maindir, "16S/bockulich_2016_single_end/", sep = ""))
# ena_file_bockulich <- read.table("filereport_read_run_PRJEB14529_tsv.txt", sep = "\t", header = T)
# ena_file_bockulich <- ena_file_bockulich[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
# 
# metadata_boc <- read.table("suppl.txt", sep = '\t', header = T)
# metadata_boc_sra <- read.table("SraRunTable.txt", header = T, sep = ",")
# metadata_boc_sra[,c("Assay.Type", "Barcode..exp.", "barcodesequence..exp.",
#                     "Bytes", "center_name..exp.", "Center.Name",
#                     "center_project_name..exp.", "Consent", "DATASTORE.filetype",
#                     "DATASTORE.provider", "DATASTORE.region", "dna_extracted",
#                     "Elevation", "emp_status..exp.", "ENA.FIRST.PUBLIC..run.",
#                     "ENA.FIRST.PUBLIC", "ENA.LAST.UPDATE..run.", "ENA.LAST.UPDATE",
#                     "env_biome", "env_feature", "geo_loc_name_country_continent",
#                     "geo_loc_name", "INSDC_center_alias", "INSDC_center_name",
#                     "INSDC_first_public", "INSDC_last_update", "INSDC_status",
#                     "Instrument", "Latitude", "LibraryLayout", "LibrarySelection",
#                     "LibrarySource", "linker..exp.", "linkerprimersequence..exp.",
#                     "LONGITUDE", "Organism", "physical_specimen_location",
#                     "physical_specimen_remaining", "Platform", "primer..exp.",
#                     "ReleaseDate", "run_center..exp.", "sequencing_meth..exp.",
#                     "target_gene..exp.", "target_subfragment..exp.", "COLLECTION_TIMESTAMP",
#                     "sample_summary", "Bases", "host_subject_id", "run_date..exp.", "Sample_Name",
#                     "SAMPLE_TYPE", "StudyID", "SRA.Study", "Submitter_Id", "month_of_life",
#                     "month")] <- NULL
# 
# 
# # no duplicates
# sum(duplicated(ena_file_bockulich$sample_title))
# 
# #combine data.frames
# sum(!unique(ena_file_bockulich$sample_title) %in% metadata_boc$sample_name)
# sum(!unique(metadata_boc$sample_name) %in% ena_file_bockulich$sample_title)
# nrow(ena_file_bockulich)
# metadata_boc_combined <- merge(metadata_boc, ena_file_bockulich, by.x = "sample_name", by.y = "sample_title")
# nrow(metadata_boc_combined)
# 
# metadata_boc_combined[,c("description", "dna_extracted", "elevation", "empo_1", "empo_2", "empo_3",
#                          "env_biome", "env_feature", "host_age", "host_age_normalized_years",
#                          "host_age_units", "host_body_mass_index", "host_common_name",
#                          "host_height", "host_height_units", "host_scientific_name",
#                          "host_taxid", "host_weight", "host_weight_units", "latitude",
#                          "longitude", "physical_specimen_location", "physical_specimen_remaining",
#                          "qiita_study_id", "sample_alias", "collection_timestamp", "geo_loc_name",
#                          "country", "month_of_life", "month", "lifestage", "mom_child", "delivery",
#                          "sample_summary", "abx1_pmp_all_bymonth", "abx_all_sources", "abx_name",
#                          "abx_pmp_all", "antiexposedall", "course_num", "diet_2", "diet_2_month",
#                          "diet_3", "host_subject_id", "sex", "sampletype", "mom_courses_ofdeliveryabx",
#                          "mom_ld_abx", "mom_pre_post_ld", "mom_prenatal_abx", "mom_prenatal_abx_class",
#                          "mom_prenatal_abx_trimester", "env_material", "baby_sex")] <- NULL
# 
# metadata_boc_combined <- merge(metadata_boc_combined, metadata_boc_sra, by.x = "sample_name", by.y = "Library.Name")
# 
# metadata_boc_combined[metadata_boc_combined == "na"] <- NA
# 
# # sort out mother samples
# sort_out_samples_boc <- metadata_boc_combined$run_accession[metadata_boc_combined$mom_child == "M"]
# metadata_boc_combined <- metadata_boc_combined[!metadata_boc_combined$mom_child == "M",]
# # no individuals with one sample only
# multiple_samples <- unique(metadata_boc_combined$studyid[duplicated(metadata_boc_combined$studyid)])
# sum(!metadata_boc_combined$studyid %in% multiple_samples)
# 
# 
# # remove duplicate columns and unify entries
# metadata_boc_combined[,c("Day_of_life", "env_package", "scientific_name",
#                          "taxon_id", "mom_courses_ofdeliveryabx", "mom_pre_post_ld", "host_body_habitat",
#                          "host_body_product", "host_body_site", "sample_type", "Run", "run_prefix..exp.",
#                          "barcodewell..exp.", "External_Id", "diet_2", "Diet", "diet_2_month", "diet_3")] <- NULL
# 
# # unify colnames
# colnames(metadata_boc_combined)[match(c("day_of_life", "delivery", "geo_loc_name_country",
#                                         "sample_name", "studyid", "delivery"),
#                                       colnames(metadata_boc_combined))] <- c("age", "delivery", "country",
#                                                                              "sample_ID", "subject_ID", "birthmode")
# 
# (colnames(metadata_boc_combined))
# metadata_boc_combined$sex[metadata_boc_combined$sex == "female"] <- "f"
# metadata_boc_combined$sex[metadata_boc_combined$sex == "male"] <- "m"
# 
# metadata_boc_combined$birthmode[metadata_boc_combined$birthmode == "Vaginal"] <- "v"
# metadata_boc_combined$birthmode[metadata_boc_combined$birthmode == "Cesarean"] <- "c"
# 
# metadata_boc_combined$sampletype[metadata_boc_combined$sampletype == "Stool_Stabilizer"] <- "SS"
# metadata_boc_combined$sampletype[metadata_boc_combined$sampletype == "Dry_Stool"] <- "SD"
# 
# metadata_boc_combined$region <- "NYC"
# metadata_boc_combined$study <- "bockulich_2016"
# 
# metadata_boc_combined$age <- as.numeric(metadata_boc_combined$age)
# 
# 
# 
# # check new files:
# paired_1 <- read.table("../bockulich_2016_paired_end/download_per_library_1/read_counts_raw_1.txt", sep = "\t")
# paired_2 <- read.table("../bockulich_2016_paired_end/download_per_library_2/read_counts_raw_2.txt", sep = "\t")
# paired_3 <- read.table("../bockulich_2016_paired_end/download_per_library_3/read_counts_raw_3.txt", sep = "\t")
# paired_4 <- read.table("../bockulich_2016_paired_end/download_per_library_4/read_counts_raw_4.txt", sep = "\t")
# paired_5 <- read.table("../bockulich_2016_paired_end/download_per_library_5/read_counts_raw_5.txt", sep = "\t")
# paired <- rbind(paired_2, paired_3, paired_4, paired_5)
# 
# paired$V1 <- gsub("[_].*", "", paired$V1)
# merged_rc <- merge(paired, metadata_boc_combined[,c("sample_ID", "read_count")], by.x = "V1", by.y = "sample_ID")
# plot(merged_rc$V2, merged_rc$read_count)
# 
# # move sorted out samples
# sort_out_samples_boc <- unique(sort_out_samples_boc)
# sort_out_files_boc <- c(paste(sort_out_samples_boc, "_1.fastq.gz", sep = ""),
#                         paste(sort_out_samples_boc, "_2.fastq.gz", sep = ""),
#                         paste(sort_out_samples_boc, ".fastq.gz", sep = ""))
# sort_out_files_boc <- sort_out_files_boc[file.exists(sort_out_files_boc)]
# file.move(sort_out_files_boc, "sorted_out/")
# 
# 
# # move good files
# final_samples_bock <- metadata_boc_combined$run_accession
# final_files_bock <- c(paste(final_samples_bock, "_1.fastq.gz", sep = ""),
#                      paste(final_samples_bock, "_2.fastq.gz", sep = ""),
#                      paste(final_samples_bock, ".fastq.gz", sep = ""))
# final_files_bock <- final_files_bock[file.exists(final_files_bock)]
# file.move(final_files_bock, "fastq_files/")
# 
# 
# # check read numbers
# n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
# plot(metadata_boc_combined$read_count, n_lines$read_counts[match(metadata_boc_combined$run_accession, n_lines$sample_ID)])  # 
# 
# # write file
# write.csv(metadata_boc_combined, "metadata_boc_2016_healthy.csv")
# 
# 
# # clean up
# rm(ena_file_bockulich, metadata_boc, metadata_boc_combined, metadata_boc_sra,
#    sort_out_files_boc, sort_out_samples_boc, final_files_bock, final_samples_bock,
#    multiple_samples, n_lines)





# Bockulich 2016 ###############################################################
# better to use information from SRA, also aligns better with sample ids
# R-rectal swab, V-vaginal swab, SS-stool sample, SD-dry stool
# not 100% sure for some metadata
# what does sampletype "Repeats" mean? maybe sort them out?

setwd(paste(maindir, "16S/bockulich_2016/", sep = ""))
ena_file_bockulich <- read.table("filereport_read_run_PRJEB14529_tsv.txt", sep = "\t", header = T)
ena_file_bockulich <- ena_file_bockulich[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]

metadata_boc <- read.table("suppl.txt", sep = '\t', header = T)
metadata_boc_sra <- read.table("SraRunTable.txt", header = T, sep = ",")
metadata_boc_sra[,c("Assay.Type", "Barcode..exp.", "barcodesequence..exp.",
                    "Bytes", "center_name..exp.", "Center.Name",
                    "center_project_name..exp.", "Consent", "DATASTORE.filetype",
                    "DATASTORE.provider", "DATASTORE.region", "dna_extracted",
                    "Elevation", "emp_status..exp.", "ENA.FIRST.PUBLIC..run.",
                    "ENA.FIRST.PUBLIC", "ENA.LAST.UPDATE..run.", "ENA.LAST.UPDATE",
                    "env_biome", "env_feature", "geo_loc_name_country_continent",
                    "geo_loc_name", "INSDC_center_alias", "INSDC_center_name",
                    "INSDC_first_public", "INSDC_last_update", "INSDC_status", "LibraryLayout", "LibrarySelection",
                    "LibrarySource", "linker..exp.", "linkerprimersequence..exp.",
                    "Organism", "physical_specimen_location",
                    "physical_specimen_remaining", "Platform", "primer..exp.",
                    "ReleaseDate", "run_center..exp.", "sequencing_meth..exp.",
                    "target_gene..exp.", "target_subfragment..exp.", "COLLECTION_TIMESTAMP",
                    "sample_summary", "Bases", "host_subject_id", "run_date..exp.", "Sample_Name",
                    "SAMPLE_TYPE", "StudyID", "SRA.Study", "Submitter_Id", "month_of_life",
                    "month")] <- NULL


# no duplicates
sum(duplicated(ena_file_bockulich$sample_title))

#combine data.frames
sum(!unique(ena_file_bockulich$sample_title) %in% metadata_boc$sample_name)
sum(!unique(metadata_boc$sample_name) %in% ena_file_bockulich$sample_title)
nrow(ena_file_bockulich)
metadata_boc_combined <- merge(metadata_boc, ena_file_bockulich, by.x = "sample_name", by.y = "sample_title")
nrow(metadata_boc_combined)

metadata_boc_combined[,c("description", "dna_extracted", "elevation", "empo_1", "empo_2", "empo_3",
                         "env_biome", "env_feature", "host_age", "host_age_normalized_years",
                         "host_age_units", "host_body_mass_index", "host_common_name",
                         "host_height", "host_height_units", "host_scientific_name",
                         "host_taxid", "host_weight", "host_weight_units", "latitude",
                         "longitude", "physical_specimen_location", "physical_specimen_remaining",
                         "qiita_study_id", "sample_alias", "collection_timestamp", "geo_loc_name",
                         "country", "month_of_life", "month", "lifestage", "mom_child", "delivery",
                         "sample_summary", "abx1_pmp_all_bymonth", "abx_all_sources", "abx_name",
                         "abx_pmp_all", "antiexposedall", "course_num", "diet_2", "diet_2_month",
                         "diet_3", "host_subject_id", "sex", "sampletype", "mom_courses_ofdeliveryabx",
                         "mom_ld_abx", "mom_pre_post_ld", "mom_prenatal_abx", "mom_prenatal_abx_class",
                         "mom_prenatal_abx_trimester", "env_material", "baby_sex")] <- NULL

metadata_boc_combined <- merge(metadata_boc_combined, metadata_boc_sra, by.x = "sample_name", by.y = "Library.Name")

metadata_boc_combined[metadata_boc_combined == "na"] <- NA

# sort out mother samples
sort_out_samples_boc <- metadata_boc_combined$sample_name[metadata_boc_combined$mom_child == "M"]
metadata_boc_combined <- metadata_boc_combined[!metadata_boc_combined$mom_child == "M",]
# no individuals with one sample only
multiple_samples <- unique(metadata_boc_combined$studyid[duplicated(metadata_boc_combined$studyid)])
sum(!metadata_boc_combined$studyid %in% multiple_samples)


# remove duplicate columns and unify entries
metadata_boc_combined[,c("Day_of_life", "env_package", "scientific_name",
                         "taxon_id", "mom_courses_ofdeliveryabx", "mom_pre_post_ld", "host_body_habitat",
                         "host_body_product", "host_body_site", "sample_type", "Run", "run_prefix..exp.",
                         "barcodewell..exp.", "External_Id", "diet_2", "Diet", "diet_2_month", "diet_3")] <- NULL

# unify colnames
metadata_boc_combined %<>% mutate(sample_ID = paste0(sample_name, "_bockulich_2016"),
                                  subject_ID = paste0(studyid, "_bockulich_2016"),
                                  birthmode = delivery,
                                  age = as.numeric(day_of_life),
                                  country = geo_loc_name_country,
                                  region = "NYC",
                                  study = "bockulich_2016",
                                  lifestyle = "industrialized",
                                  family_ID = subject_ID,
                                  geographic_location_.latitude. = Latitude,
                                  geographic_location_.longitude. = LONGITUDE,
                                  antibiotics_before = ifelse(antiexposedall == "y", yes = T, no = F),
                                  .keep = "unused")


metadata_boc_combined$sex[metadata_boc_combined$sex == "female"] <- "f"
metadata_boc_combined$sex[metadata_boc_combined$sex == "male"] <- "m"

metadata_boc_combined$birthmode[metadata_boc_combined$birthmode == "Vaginal"] <- "v"
metadata_boc_combined$birthmode[metadata_boc_combined$birthmode == "Cesarean"] <- "c"

metadata_boc_combined$sampletype[metadata_boc_combined$sampletype == "Stool_Stabilizer"] <- "SS"
metadata_boc_combined$sampletype[metadata_boc_combined$sampletype == "Dry_Stool"] <- "SD"


# check new files:
paired_1 <- read.table("download_per_library_1/read_counts_raw_1.txt", sep = "\t")
paired_2 <- read.table("download_per_library_2/read_counts_raw_2.txt", sep = "\t")
paired_3 <- read.table("download_per_library_3/read_counts_raw_3.txt", sep = "\t")
paired_4 <- read.table("download_per_library_4/read_counts_raw_4.txt", sep = "\t")
paired_5 <- read.table("download_per_library_5/read_counts_raw_5.txt", sep = "\t")
paired <- rbind(paired_2, paired_3, paired_4, paired_5)

paired$V1 <- gsub("[_].*", "", paired$V1)
merged_rc <- merge(paired %>% mutate(V1 = paste0(V1, "_bockulich_2016")),
                   metadata_boc_combined[,c("sample_ID", "read_count")], by.x = "V1", by.y = "sample_ID")
plot(merged_rc$V2, merged_rc$read_count)

# move sorted out samples
sort_out_samples_boc <- unique(sort_out_samples_boc)
sort_out_files_boc <- list.files(pattern = paste0(sort_out_samples_boc, collapse = "|"))
file.move(sort_out_files_boc, "sorted_out/")

# 
# # move good files
# final_samples_bock <- metadata_boc_combined$sample_ID
# final_files_bock <- list.files(pattern = paste0(final_samples_bock, collapse = "|"))
# file.move(final_files_bock, "fastq_files/")


# check read numbers
# n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
# plot(metadata_boc_combined$read_count, n_lines$read_counts[match(metadata_boc_combined$sample_ID, n_lines$sample_ID)])  # 

# write file
write.csv(metadata_boc_combined, "metadata_boc_2016_healthy.csv")


# clean up
rm(ena_file_bockulich, metadata_boc, metadata_boc_combined, metadata_boc_sra,
   sort_out_files_boc, sort_out_samples_boc, final_files_bock, final_samples_bock,
   multiple_samples, n_lines, merged_rc, paired, paired_1, paired_2, paired_3,
   paired_4, paired_5)




# Hill #########################################################################
# Issues: Not clear, when which food is introduced. For some samples feeding information is lacking

setwd(paste(maindir, "16S/hill_2017/", sep = ""))
ena_file_hill <- read.table("filereport_read_run_PRJNA339264_tsv.txt", sep = "\t", header = T)
ena_file_hill <- ena_file_hill[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias")]
ena_file_hill <- separate(ena_file_hill, "sample_alias", into = c(NA,"ID", NA, "age"),
                          sep = c(2,5,6), convert = T, remove = F)
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file[,c("Instrument","Run")]
all(ena_file_hill$run_accession %in% sra_file$Run)
ena_file_hill <- merge(ena_file_hill, sra_file,by.x = "run_accession", by.y = "Run")

ena_file_hill$age <- gsub("a", "", ena_file_hill$age) %>% 
   gsub(pattern = "b", replacement = "") %>%
   gsub(pattern = "c", replacement = "") %>%
   gsub(pattern = "ue", replacement = 0) %>%
   gsub(pattern = "stori.*", replacement = "") %>%
   as.numeric()
ena_file_hill$age <- ena_file_hill$age * 7


metadata_hil <- read.table("suppl_1.txt", sep = '\t', header = T)
metadata_hil <- separate(metadata_hil, "Gestational.Age...weeks...days.", into = c("Gestational_age", "Gestational_days"),
                          sep = "[+]", convert = T, remove = T)
metadata_hil$Gestational_days[is.na(metadata_hil$Gestational_days)] <- 0
metadata_hil$Gestational_age <- metadata_hil$Gestational_age + metadata_hil$Gestational_days/7
metadata_hil$Actual.Breastfeeding.Duration[metadata_hil$Actual.Breastfeeding.Duration == ""] <- 
   metadata_hil$Breast.Feeding.Status[metadata_hil$Actual.Breastfeeding.Duration == ""]
metadata_hil[,c("Gestational_days", "Mother.breast.fed.as.child.", "Mother.Smoker",
                "Other.persons.in.home.smoker", "Mother.BF.as.a.child", "Breast.Feeding.Status")] <- NULL



metadata_hill_preterms <- read.table("suppl_2_preterm.txt", sep = '\t', header = T)
metadata_hill_preterms[c("Delivery.Mode", "Gestational.Age..days.", "Gestational.Age..weeks.", "Total.Days")] <- NULL

# no duplicates
sum(duplicated(ena_file_hill$sample_alias))

#combine data.frames
sum(!unique(ena_file_hill$ID) %in% metadata_hil$Subject.Number)
sum(!unique(metadata_hil$Subject.Number) %in% ena_file_hill$ID)
nrow(ena_file_hill)
metadata_hil_combined <- merge(metadata_hil, ena_file_hill, by.x = "Subject.Number", by.y = "ID")
nrow(metadata_hil_combined)

nrow(metadata_hil_combined)
metadata_hil_combined <- merge(metadata_hil_combined, metadata_hill_preterms, by = "Subject.Number", all.x = T)
nrow(metadata_hil_combined)

# remove all individuals with only one sample
multiple_samples <- unique(metadata_hil_combined$Subject.Number[duplicated(metadata_hil_combined$Subject.Number)])
sort_out_samples_hil <- metadata_hil_combined$run_accession[!metadata_hil_combined$Subject.Number %in% multiple_samples]
metadata_hil_combined <- metadata_hil_combined[metadata_hil_combined$Subject.Number %in% multiple_samples,]

# change column entries
metadata_hil_combined$preterm <- F
metadata_hil_combined$preterm[metadata_hil_combined$Subject.Number %in% metadata_hill_preterms$Subject.Number] <- T
sort_out_samples_hil <- c(sort_out_samples_hil, metadata_hil_combined$run_accession[metadata_hil_combined$preterm])
metadata_hil_combined <- metadata_hil_combined[!metadata_hil_combined$preterm,]
metadata_hil_combined$Gender[metadata_hil_combined$Gender == "Female"] <- "f"
metadata_hil_combined$Gender[metadata_hil_combined$Gender == "Male"] <- "m"
metadata_hil_combined$Delivery.Mode[metadata_hil_combined$Delivery.Mode == "CS"] <- "c"
metadata_hil_combined$Delivery.Mode[metadata_hil_combined$Delivery.Mode %in% c("SVD", "Kiwi Vacum", "NBFD", "Vacum", "forceps")] <- "v"
metadata_hil_combined$Delivery.Mode[metadata_hil_combined$Delivery.Mode == ""] <- NA
metadata_hil_combined$Formula.Used[metadata_hil_combined$Formula.Used %in% c("None", "none", "no formula", "No Formula", "No formula")] <- "n"
metadata_hil_combined$Formula.Used[!metadata_hil_combined$Formula.Used %in% c("n", "")] <- "y"
metadata_hil_combined[metadata_hil_combined == "Yes"] <- T
metadata_hil_combined[metadata_hil_combined == "No"] <- F
metadata_hil_combined[metadata_hil_combined == "yes"] <- T
metadata_hil_combined[metadata_hil_combined == "no"] <- F

# feeding
metadata_hil_combined$food <- NA
t_5_hill <- c("18 months", "7 months", "7", "6 months", "10 months", "12 months +",
              "8 months", "11 months", "12 months", "12 months ", "7.5 months",
              "10 Months", "9 months", "9.5 months", "12months +")
t_4_hill <- c("3 months", "3.5 months", "8 weeks", "5 months", "4 months",
              "2 months", "5.5 months", "4.5 months", "10 weeks", "12 weeks",
              "Less than 4 months", "Greater than 4 months", "greater than 4 months")
t_3_hill <- c("6 weeks", "5 weeks", "1 month", "1 months")
metadata_hil_combined$food[metadata_hil_combined$Actual.Breastfeeding.Duration == "No"] <- "formula"
metadata_hil_combined$food[metadata_hil_combined$Actual.Breastfeeding.Duration %in% t_5_hill] <- "breast"
metadata_hil_combined$food[metadata_hil_combined$age <= 56 &
                                 metadata_hil_combined$Actual.Breastfeeding.Duration %in% t_4_hill] <- "breast"
metadata_hil_combined$food[metadata_hil_combined$age <= 28 &
                                 metadata_hil_combined$Actual.Breastfeeding.Duration %in% t_3_hill] <- "breast"
metadata_hil_combined$food[metadata_hil_combined$age <= 7 &
                                 metadata_hil_combined$Actual.Breastfeeding.Duration == "2 weeks"] <- "breast"
# fill up all NAs with formula, if formula is used at any time point
metadata_hil_combined$food[is.na(metadata_hil_combined$food) &
                                 metadata_hil_combined$Formula.Used == "y"] <- "formula"
# every infant younger than 112 days is not going to receive any solid food
metadata_hil_combined$food[metadata_hil_combined$age <=56 & 
                                 is.na(metadata_hil_combined$food) &
                                 metadata_hil_combined$Actual.Breastfeeding.Duration != ""] <- "formula"
# assume that non breast or formula fed infants at this age receive solid food
metadata_hil_combined$food[metadata_hil_combined$age == 168 & is.na(metadata_hil_combined$food)] <- "solid"
# every child receiving formula and breast, but no time point for formula introduction can be inferred is assumed to be mixed_b_f
subj_assigned_form <- unique(metadata_hil_combined$Subject.Number[metadata_hil_combined$food == "formula"])
subj_orig_form <- unique(metadata_hil_combined$Subject.Number[metadata_hil_combined$Formula.Used == "y"])
no_intersect_ids <- subj_orig_form[!subj_orig_form %in% subj_assigned_form]
no_intersect_index <- metadata_hil_combined$Subject.Number %in% no_intersect_ids
metadata_hil_combined$food[no_intersect_index] <- "mixed_b_f"

# delet cols
metadata_hil_combined[c("Days.in.neonatal.unit", "Benzylpenicillin",
                         "Gentamicin", "Benzylpenicilin...Gentamicin",
                         "Teicoplanin","flucloxicillin","Meropenem", "Second.course",
                         "preterm", "CS.Type", "Weight..g.")] <- NULL

# unify colnames
metadata_hil_combined %<>% mutate(sample_ID = paste0(sample_alias, "_hill_2017"),
                                  subject_ID = paste0(Subject.Number, "_hill_2017"),
                                  birthmode = Delivery.Mode,
                                  gestational_age = Gestational_age,
                                  birth_weight = Birth.Weight..g. / 1000,
                                  antibiotics_any = case_when(Infant.Antibiotics.Taken.after.initial.discharge.from.hospital == "" ~ NA,
                                                              Infant.Antibiotics.Taken.after.initial.discharge.from.hospital == "FALSE" ~ F,
                                                              Infant.Antibiotics.Taken.after.initial.discharge.from.hospital == "augmentin" ~ T,
                                                              Infant.Antibiotics.Taken.after.initial.discharge.from.hospital == "TRUE" ~ T),
                                  formula = Formula.Used,
                                  sex = Gender,
                                  country = "IRELAND",
                                  region = "CORK",
                                  study = "hill_2017",
                                  lifestyle = "industrialized",
                                  family_ID = subject_ID,
                                  geographic_location_.latitude. = 51.903614,
                                  geographic_location_.longitude. = -8.468399,
                                  .keep = "unused")

# move sorted out samples
sort_out_samples_hil <- unique(sort_out_samples_hil)
sort_out_files_hil <- c(paste(sort_out_samples_hil, "_1.fastq.gz", sep = ""),
                        paste(sort_out_samples_hil, "_2.fastq.gz", sep = ""),
                        paste(sort_out_samples_hil, ".fastq.gz", sep = ""))
sort_out_files_hil <- sort_out_files_hil[file.exists(sort_out_files_hil)]
file.move(sort_out_files_hil, "sorted_out/")

# move good files
final_samples_hil <- metadata_hil_combined$run_accession
final_files_hil <- c(paste(final_samples_hil, "_1.fastq.gz", sep = ""),
                     paste(final_samples_hil, "_2.fastq.gz", sep = ""),
                     paste(final_samples_hil, ".fastq.gz", sep = ""))
final_files_hil <- final_files_hil[file.exists(final_files_hil)]
file.move(final_files_hil, "fastq_files/")

# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_hil_combined$read_count, n_lines$read_counts[match(metadata_hil_combined$run_accession, n_lines$sample_ID)])  # 


# write file
write.csv(metadata_hil_combined, "metadata_hil_2017_healthy.csv")


# clean up
rm(ena_file_hill, metadata_hil, metadata_hil_combined, metadata_hill_preterms,
   sort_out_files_hil, sort_out_samples_hil, multiple_samples,
   no_intersect_ids, no_intersect_index, subj_assigned_form, subj_orig_form,
   t_3_hill, t_4_hill, t_5_hill, final_files_hil, final_samples_hil, n_lines, sra_file)



# Sprockett 2020 ###############################################################
# assume no access to antibiotics

setwd(paste(maindir, "16S/sprockett_2020/", sep = ""))
ena_file_sprockett <- read.table("filereport_read_run_PRJNA574920_tsv.txt", sep = "\t", header = T)
ena_file_sprockett <- ena_file_sprockett[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file[,c("Instrument","Run")]
all(ena_file_sprockett$run_accession %in% sra_file$Run)
ena_file_sprockett <- merge(ena_file_sprockett, sra_file,by.x = "run_accession", by.y = "Run")


metadata_spr <- read.table("suppl.txt", sep = '\t', header = T)
metadata_spr$Age_Days[is.na(metadata_spr$Age_Days)] <- metadata_spr$Age_Months[is.na(metadata_spr$Age_Days)]*30
metadata_spr[c("Age_Months", "Age_Years", "Days_After_Birth", "Months_After_Birth",
               "Age_Group", "Age_Class", "AgeGroup_By_SampleType", "Complete_Dyads_Saliva",
               "Study_Code")] <- NULL

metadata_spr_neopt <- read.table("Tsimane_Neopterin_Data.txt", sep = '\t', header = T)

# no duplicates
sum(duplicated(ena_file_sprockett$sample_alias))

#combine data.frames
sum(!unique(ena_file_sprockett$sample_alias) %in% metadata_spr$Sample_ID)
sum(!unique(metadata_spr$Sample_ID) %in% ena_file_sprockett$sample_alias)
nrow(ena_file_sprockett)
metadata_spr_combined <- merge(metadata_spr, ena_file_sprockett, by.x = "Sample_ID", by.y = "sample_alias")
nrow(metadata_spr_combined)

# remove saliva and mother samples:
sort_out_samples_spr <- metadata_spr_combined$run_accession[metadata_spr_combined$Sample_Type == "Saliva"]
metadata_spr_combined <- metadata_spr_combined[metadata_spr_combined$Sample_Type != "Saliva", ]
# only keep children below 3 years
sort_out_samples_spr <- c(sort_out_samples_spr, metadata_spr_combined$run_accession[metadata_spr_combined$Age_Days > 1095])
metadata_spr_combined <- metadata_spr_combined[metadata_spr_combined$Age_Days <= 1095,]

# remove all individuals with only one sample
multiple_samples <- unique(metadata_spr_combined$Subject_ID[duplicated(metadata_spr_combined$Subject_ID)])
# sort_out_samples_spr <- c(sort_out_samples_spr, 
#                           metadata_spr_combined$run_accession[!metadata_spr_combined$Subject_ID %in% multiple_samples])
# metadata_spr_combined <- metadata_spr_combined[metadata_spr_combined$Subject_ID %in% multiple_samples,]
# retreive sorted out samples
metadata_spr_combined$run_accession[!metadata_spr_combined$Subject_ID %in% multiple_samples] %>%
  paste("fastq-dump --split-e --gzip ",.) %>%
  data.frame(`#!/bin/bash` = .) %>%
  write.table(., file = "/fast/AG_Forslund/rob/studies/16S/sprockett_2020//pull_single.sh",
              sep = " ", quote = F, row.names = F, col.names = "#!/bin/bash")

# rename entries
metadata_spr_combined[metadata_spr_combined == "Yes"] <- T
metadata_spr_combined[metadata_spr_combined == "No"] <- F
metadata_spr_combined$Sex[metadata_spr_combined$Sex == "Female"] <- "f"
metadata_spr_combined$Sex[metadata_spr_combined$Sex == "Male"] <- "m"
metadata_spr_combined$Feeding_Status_Corrected[metadata_spr_combined$Feeding_Status_Corrected == "Exclusively_Breastfed"] <- "breast"
metadata_spr_combined$Feeding_Status_Corrected[metadata_spr_combined$Feeding_Status_Corrected == "Complementary_Foods"] <- "mixed_b_f"
metadata_spr_combined$Delivery_Mode <- "v"

# clean cols
metadata_spr_combined[c("Sample_Label", "Sample_Type", "Country", "Complete_Dyads_Feces")] <- NULL

# unify colnames
metadata_spr_combined %<>% mutate(sample_ID = paste0(Sample_ID, "_sprockett_2020"),
                                  subject_ID = paste0(Subject_ID, "_sprockett_2020"),
                                  antibiotics_any = F,
                                  antibiotics_before = F,
                                  antibiotics_one_week_before = F,
                                  age = Age_Days,
                                  food = Feeding_Status_Corrected,
                                  birthmode = Delivery_Mode,
                                  sex = Sex,
                                  country = "BOLIVIA",
                                  region = "TSIMANE",
                                  study = "sprockett_2020",
                                  lifestyle = "non_industrialized",
                                  family_ID = Dyad,
                                  geographic_location_.latitude. = -14.86537907781345,
                                  geographic_location_.longitude. = -66.72422768406338,
                                  .keep = "unused")

# move sorted out samples
sort_out_samples_spr <- unique(sort_out_samples_spr)
sort_out_files_spr <- c(paste(sort_out_samples_spr, "_1.fastq.gz", sep = ""),
                        paste(sort_out_samples_spr, "_2.fastq.gz", sep = ""),
                        paste(sort_out_samples_spr, ".fastq.gz", sep = ""))
sort_out_files_spr <- sort_out_files_spr[file.exists(sort_out_files_spr)]
file.move(sort_out_files_spr, "sorted_out/")

# move good files
final_samples_spr <- metadata_spr_combined$run_accession
final_files_spr <- c(paste(final_samples_spr, "_1.fastq.gz", sep = ""),
                     paste(final_samples_spr, "_2.fastq.gz", sep = ""),
                     paste(final_samples_spr, ".fastq.gz", sep = ""))
final_files_spr <- final_files_spr[file.exists(final_files_spr)]
file.move(final_files_spr, "fastq_files/")


# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_spr_combined$read_count, n_lines$read_counts[match(metadata_spr_combined$run_accession, n_lines$sample_ID)])  # 


# write file
write.csv(metadata_spr_combined, "metadata_spr_2020_healthy.csv")


# clean up
rm(ena_file_sprockett, metadata_spr, metadata_spr_combined, metadata_spr_neopt,
   multiple_samples, sort_out_files_spr, sort_out_samples_spr, final_samples_spr,
   final_files_spr, n_lines, sra_file)




# Raman 2019 ###################################################################

setwd(paste(maindir, "16S/raman_2019/", sep = ""))
ena_file_raman <- read.table("filereport_read_run_PRJEB27068_tsv_human.txt", sep = "\t", header = T)
ena_file_raman <- ena_file_raman[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file[,c("Instrument","Run")]
all(ena_file_raman$run_accession %in% sra_file$Run)
ena_file_raman <- merge(ena_file_raman, sra_file,by.x = "run_accession", by.y = "Run")


metadata_ram_1 <- read.table("supp_t1.txt", sep = '\t', header = T)
metadata_ram_1$Cohort <- "Bangladesh"
metadata_ram_1$Age <- metadata_ram_1$Age..months. * 30
metadata_ram_1$Age..months. <- NULL
metadata_ram_1$Monthly.Timepoint <- NULL
metadata_ram_1 <- separate(metadata_ram_1, "SampleID", into = c("Child_ID", NA),
                           sep = "[.]", convert = T, remove = F)
# expand the sex information to all samples
sex_dataframe <- metadata_ram_1[metadata_ram_1$Sex != "", c("Sex", "Child_ID")]
metadata_ram_1$Sex <- NULL
metadata_ram_1 <- inner_join(metadata_ram_1, sex_dataframe, by = "Child_ID")


metadata_ram_4 <- read.table("supp_t4.txt", sep = '\t', header = T)
colnames(metadata_ram_4)[3] <- "SampleID"
metadata_ram_4$Age <- metadata_ram_4$Age.of.donor.at.time.of.fecal.sample.collection..months. * 30
metadata_ram_4$Age.of.donor.at.time.of.fecal.sample.collection..months. <- NULL
metadata_ram_4$Child.ID <- NULL
metadata_ram_4 <- separate(metadata_ram_4, "SampleID", into = c("Child_ID", NA, NA),
                         sep = "[.]", convert = T, remove = F)

metadata_ram <- plyr::rbind.fill(metadata_ram_1, metadata_ram_4)

# no duplicates
sum(duplicated(ena_file_raman$sample_alias))

#combine data.frames
sum(!unique(ena_file_raman$sample_alias) %in% metadata_ram$SampleID) # some samples are in the second file
sum(!unique(metadata_ram$SampleID) %in% ena_file_raman$sample_alias)
nrow(ena_file_raman)
metadata_ram_combined <- merge(metadata_ram, ena_file_raman, by.x = "SampleID", by.y = "sample_alias")
nrow(metadata_ram_combined)
# 
# # which samples are not in raman, but in subr.
# samples_sub$Birth.Cohort[!samples_sub$sample_ID %in% metadata_ram_combined$SampleID] %>% unique()
# 
# # check read counts in the two datasets
# n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
# common_samples <- unique(samples_sub$sample_ID[samples_sub$sample_ID %in% metadata_ram_combined$SampleID])
# common_samples_ram <- metadata_ram_combined[metadata_ram_combined$SampleID %in% common_samples,]
# common_samples_subr <- samples_sub[samples_sub$sample_ID %in% common_samples,]
# sum(duplicated(common_samples_subr$sample_ID))
# common_samples_subr <- common_samples_subr[!duplicated(common_samples_subr$sample_ID),]
# plot(common_samples_ram$read_count, common_samples_subr$read_count[match(common_samples_ram$SampleID, common_samples_subr$sample_ID)])  # 
# (common_samples_ram$read_count != 
#       common_samples_subr$read_count[match(common_samples_ram$SampleID,
#                                            common_samples_subr$sample_ID)])  %>% sum


# # remove samples from subramanian
# samples_sub <- read.table("../subramanian_2014/metadata_subramanian_2014_healthy.csv", sep = ",", header = T)
# sort_out_samples_ram <- metadata_ram_combined$run_accession[metadata_ram_combined$SampleID %in% samples_sub$sample_ID]
# metadata_ram_combined <- metadata_ram_combined[!metadata_ram_combined$SampleID %in% samples_sub$sample_ID,]



# combine barcode columns
metadata_ram_combined$Barcode.Sequence[is.na(metadata_ram_combined$Barcode.Sequence)] <- 
   metadata_ram_combined$Barcode.sequence[is.na(metadata_ram_combined$Barcode.Sequence)]
metadata_ram_combined$Barcode.sequence <- NULL

# add antibiotics information from subramanian
abx_subr_data <- read.table("suppl_2_subr.txt", sep = "\t", header = T) %>%
  select(Fecal.Sample.ID, Breast.Milk, Formula1, Solid.Foods2, diarrhoea.at.the.time.of.sample.collection3,
         Antibiotics.within.7.days.prior.to.sample.collection)
metadata_ram_combined <- left_join(metadata_ram_combined, abx_subr_data, by = join_by(SampleID == Fecal.Sample.ID))

# remove diarrhoea samples
metadata_ram_combined <- metadata_ram_combined %>%
  filter(is.na(diarrhoea.at.the.time.of.sample.collection3) | diarrhoea.at.the.time.of.sample.collection3 == "No")
# remove all individuals with only one sample
multiple_samples <- unique(metadata_ram_combined$Child_ID[duplicated(metadata_ram_combined$Child_ID)])
sort_out_samples_ram <- metadata_ram_combined$run_accession[!metadata_ram_combined$Child_ID %in% multiple_samples]
metadata_ram_combined <- metadata_ram_combined[metadata_ram_combined$Child_ID %in% multiple_samples,]

# unify colnames
metadata_ram_combined %<>% mutate(sample_ID = paste0(SampleID, "_raman_2019"),
                                  subject_ID = paste0(Child_ID, "_raman_2019"),
                                  country = gsub(" ", "_", toupper(Cohort)),
                                  region = country,
                                  age = Age,
                                  sex = Sex,
                                  antibiotics_one_week_before = ifelse(Antibiotics.within.7.days.prior.to.sample.collection == "Yes", yes = T, no = F),
                                  food = case_when(Breast.Milk == "Yes" & Formula1 == "No" & Solid.Foods2 == "No" ~ "breast",
                                                   Breast.Milk == "Yes" & Formula1 == "Yes" & Solid.Foods2 == "No" ~ "mixed_b_f",
                                                   Breast.Milk == "Yes" & Formula1 == "No" & Solid.Foods2 == "Yes" ~ "mixed_b_s",
                                                   Breast.Milk == "Yes" & Formula1 == "Yes" & Solid.Foods2 == "Yes" ~ "mixed_b_f_s",
                                                   Breast.Milk == "No" & Formula1 == "Yes" & Solid.Foods2 == "No" ~ "formula",
                                                   Breast.Milk == "No" & Formula1 == "Yes" & Solid.Foods2 == "Yes" ~ "mixed_f_s",
                                                   Breast.Milk == "No" & Formula1 == "No" & Solid.Foods2 == "Yes" ~ "solid",
                                                   .default = NA),
                                  study = "raman_2019",
                                  lifestyle = "non_industrialized",
                                  family_ID = subject_ID,
                                  geographic_location_.latitude. = case_when(country == "BANGLADESH" ~ 23.933,
                                                                             country == "BRAZIL" ~ -3.731862,
                                                                             country == "INDIA" ~ 12.934968,
                                                                             country == "PERU" ~ -4.0000,
                                                                             country == "SOUTH_AFRICA" ~ -22.9456),
                                  geographic_location_.longitude. = case_when(country == "BANGLADESH" ~ 88.983,
                                                                              country == "BRAZIL" ~ -38.526669,
                                                                              country == "INDIA" ~ 79.146881,
                                                                              country == "PERU" ~ -74.5000,
                                                                              country == "SOUTH_AFRICA" ~ 30.4850),
                                  .keep = "unused")

metadata_ram_combined$region[metadata_ram_combined$country == "BANGLADESH"] <- "DHAKA"
metadata_ram_combined$region[metadata_ram_combined$country == "BRAZIL"] <- "FORTALEZA"
metadata_ram_combined$lifestyle[metadata_ram_combined$country == "BRAZIL"] <- "industrialized"
metadata_ram_combined$region[metadata_ram_combined$country == "INDIA"] <- "VELLORE"
metadata_ram_combined$region[metadata_ram_combined$country == "PERU"] <- "LORETO"
metadata_ram_combined$region[metadata_ram_combined$country == "SOUTH_AFRICA"] <- "VENDA"

metadata_ram_combined$sex[metadata_ram_combined$sex == "female"] <- "f"
metadata_ram_combined$sex[metadata_ram_combined$sex == "male"] <- "m"
sort(colnames(metadata_ram_combined))

# add height and weight data:
metadata_ram_combined$height <- sitar::LMS2z(metadata_ram_combined$age/365, y = metadata_ram_combined$HAZ, sex = metadata_ram_combined$sex,
                                      measure = "ht", ref = who06, toz = F)
metadata_ram_combined$weight <- sitar::LMS2z(metadata_ram_combined$age/365, y = metadata_ram_combined$WAZ, sex = metadata_ram_combined$sex,
                                      measure = "wt", ref = who06, toz = F) 

# move sorted out samples
# sort_out_samples_ram <- ena_file_raman$run_accession[!ena_file_raman$run_accession %in% metadata_ram_combined$run_accession]
# sort_out_samples_ram <- unique(sort_out_samples_ram)
# sort_out_files_ram <- c(paste(sort_out_samples_ram, "_1.fastq.gz", sep = ""),
#                         paste(sort_out_samples_ram, "_2.fastq.gz", sep = ""),
#                         paste(sort_out_samples_ram, ".fastq.gz", sep = ""))
# sort_out_files_ram <- sort_out_files_ram[file.exists(sort_out_files_ram)]
# file.move(sort_out_files_ram, "sorted_out/")

# move good files
# final_samples_ram <- metadata_ram_combined$run_accession
# final_files_ram <- c(paste(final_samples_ram, "_1.fastq.gz", sep = ""),
#                      paste(final_samples_ram, "_2.fastq.gz", sep = ""),
#                      paste(final_samples_ram, ".fastq.gz", sep = ""))
# final_files_ram <- final_files_ram[file.exists(final_files_ram)]
# file.move(final_files_ram, "fastq_files/")


# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_ram_combined$read_count, n_lines$read_counts[match(metadata_ram_combined$run_accession, n_lines$sample_ID)])  # 
 # some samples have many unpaired reads, think about what to do with that

# write file
write.csv(metadata_ram_combined, "metadata_ram_2019_healthy.csv")
# write.csv(metadata_ram_combined_malnurished, "metadata_ram_2019_condition.csv")

# clean up
rm(ena_file_raman, metadata_ram, metadata_ram_1, metadata_ram_4, metadata_ram_combined,
   sex_dataframe, multiple_samples, sort_out_files_ram, sort_out_samples_ram,
   final_samples_ram, final_files_ram, n_lines, sra_file, abx_subr_data)


# Kamng’ona	2019 ###############################################################
# only SraRunTable.txt available
setwd(paste(maindir, "16S/kamngona_2019/", sep = ""))
ena_file_kamngona <- read.table("filereport_read_run_PRJEB29433_tsv.txt", sep = "\t", header = T)
ena_file_kamngona <- ena_file_kamngona[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]

metadata_kamngona <- read.table("SraRunTable.txt", sep = ",", header = T)
# remove columns
metadata_kamngona[,c("Assay.Type", "Bytes", "Center.Name", "Consent",
                     "DATASTORE.filetype", "DATASTORE.provider", "DATASTORE.region",
                     "ENA.FIRST.PUBLIC..run.", "ENA.FIRST.PUBLIC",
                     "ENA.LAST.UPDATE..run.", "ENA.LAST.UPDATE", "INSDC_center_alias",
                     "INSDC_center_name", "INSDC_first_public", "INSDC_last_update",
                     "INSDC_status", "Library.Name", "LibraryLayout",
                     "LibrarySelection", "LibrarySource", "Platform", "project_title",
                     "ReleaseDate", "SAMPLE_TYPE", "Sequencing_method", "SRA.Study",
                     "Submitter_Id", "Organism", "common_name")] <- NULL

# no duplicates
sum(duplicated(metadata_kamngona$Sample_Name))

# subject_ID
metadata_kamngona$subject_ID <- paste(metadata_kamngona$Family_Id, metadata_kamngona$Family_member, sep = "_")

#remove controls
sort_out_samples_kamn <- metadata_kamngona$Run[grep("InternalControl", metadata_kamngona$Sample_Name)]
metadata_kamngona <- metadata_kamngona[-grep("InternalControl", metadata_kamngona$Sample_Name),]
# remove mother samples
sort_out_samples_kamn <- c(sort_out_samples_kamn, metadata_kamngona$Run[metadata_kamngona$Family_member == "Mother"])
metadata_kamngona <- metadata_kamngona[metadata_kamngona$Family_member != "Mother",]
# individuals with only one sample
multiple_samples <- unique(metadata_kamngona$subject_ID[duplicated(metadata_kamngona$subject_ID)])
sort_out_samples_kamn <- c(sort_out_samples_kamn,
                           metadata_kamngona$Run[!metadata_kamngona$subject_ID %in% multiple_samples])
metadata_kamngona <- metadata_kamngona[metadata_kamngona$subject_ID %in% multiple_samples,]

# age column:
metadata_kamngona$age <- NA
metadata_kamngona$age[metadata_kamngona$age_in_months == "1M"] <- 30
metadata_kamngona$age[metadata_kamngona$age_in_months == "6M"] <- 182
metadata_kamngona$age[metadata_kamngona$age_in_months == "12M"] <- 365
metadata_kamngona$age[metadata_kamngona$age_in_months == "18M"] <- 547
metadata_kamngona$age[metadata_kamngona$age_in_months == "30M"] <- 912
metadata_kamngona$age_in_months <- NULL


# unify colnames
metadata_kamngona %<>% mutate(run_accession = Run,
                              sample_ID = paste0(Sample_Name, "_kamngona_2019"),
                              subject_ID = paste0(subject_ID, "_kamngona_2019"),
                              country = "MALAWI",
                              region = "MALAWI",
                              study = "kamngona_2019",
                              lifestyle = "non_industrialized",
                              family_ID = Family_Id,
                              geographic_location_.latitude. = -14.486173,
                              geographic_location_.longitude. = 35.253304,
                              .keep = "unused")

# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_kamngona$Bases, n_lines$read_counts[match(metadata_kamngona$run_accession, n_lines$sample_ID)])  # 


write.csv(metadata_kamngona, "metadata_kamngona_2019_healthy.csv")

# move sorted out samples
sort_out_samples_kamn <- unique(sort_out_samples_kamn)
sort_out_files_kamn <- c(paste(sort_out_samples_kamn, "_1.fastq.gz", sep = ""),
                          paste(sort_out_samples_kamn, "_2.fastq.gz", sep = ""),
                          paste(sort_out_samples_kamn, ".fastq.gz", sep = ""))
sort_out_files_kamn <- sort_out_files_kamn[file.exists(sort_out_files_kamn)]
file.move(sort_out_files_kamn, "sorted_out/")

# move good files
final_samples_kamn <- metadata_kamngona$run_accession
final_files_kamn <- c(paste(final_samples_kamn, "_1.fastq.gz", sep = ""),
                       paste(final_samples_kamn, "_2.fastq.gz", sep = ""),
                       paste(final_samples_kamn, ".fastq.gz", sep = ""))
final_files_kamn <- final_files_kamn[file.exists(final_files_kamn)]
file.move(final_files_kamn, "fastq_files/")

# clean up
rm(ena_file_kamngona, metadata_kamngona, final_files_kamn, final_samples_kamn,
   multiple_samples, sort_out_files_kamn, sort_out_samples_kamn, n_lines)


# Pannaraj	2017 ###############################################################

setwd(paste(maindir, "16S/pannaraj_2017/", sep = ""))
library(xml2)
library(XML)

ncbi_xml_file_pan <- read_xml("biosample_result.xml")
ncbi_xml_pan <- xmlParse(ncbi_xml_file_pan)
ncbi_dict_pan <- xmlToDataFrame(nodes=getNodeSet(ncbi_xml_pan, "//Ids"))
colnames(ncbi_dict_pan) <- c("run_accession", "sample_ID", "other_sample_ID")

metadata_pan <- read.table("metadata.txt", sep = '\t', header = T, comment.char = "")
metadata_pan[c("Exclusion", "ChosenDuplicate", "HasDouble", "Pair", "Visit", "Treatment",
               "CSectionPlanned", "BinnedAge", "BinnedAge2", "BinnedAge3", "BinnedAge4",
               "BinnedAge5", "PercentBFBinned", "PercentBFBinned2", "PercentBFBinned3",
               "PercentBFBinned4", "PercentBFBinned2WhileYoung", "SampleType", "Sample")] <- NULL

# no duplicates
sum(duplicated(ncbi_dict_pan$run_accession))

#combine data.frames
sum(!unique(ncbi_dict_pan$sample_ID) %in% metadata_pan$X.SampleID) # no metadata for some samples
sum(!unique(metadata_pan$X.SampleID) %in% ncbi_dict_pan$sample_ID)
metadata_pan_combined <- merge(metadata_pan, ncbi_dict_pan, by.x = "X.SampleID", by.y = "sample_ID")


# remove mother samples 
sort_out_samples_pan <- metadata_pan_combined$run_accession[metadata_pan_combined$MomBB == "Mom"]
metadata_pan_combined <- metadata_pan_combined[metadata_pan_combined$MomBB != "Mom",]
metadata_pan_combined$MomBB <- NULL
# remove all individuals with only one sample
multiple_samples <- unique(metadata_pan_combined$PairID[duplicated(metadata_pan_combined$PairID)])
sort_out_samples_pan <- c(sort_out_samples_pan, metadata_pan_combined$run_accession[!metadata_pan_combined$PairID %in% multiple_samples])
metadata_pan_combined <- metadata_pan_combined[metadata_pan_combined$PairID %in% multiple_samples,]

# unify colnames
metadata_pan_combined %<>% mutate(sample_ID = paste0(X.SampleID, "_pannaraj_2017"),
                                  subject_ID = paste0(PairID, "_pannaraj_2017"),
                                  age = BBAge,
                                  birthmode = Delivery,
                                  sex = BabyGender,
                                  antibiotics_any = case_when(AnyInfantAbx == "" ~ NA,
                                                              AnyInfantAbx == "N" ~ F,
                                                              AnyInfantAbx == "Y" ~ T),
                                  antibiotics_one_week_before = case_when(BabyAbxSinceLast14days == "" ~ NA,
                                                                          BabyAbxSinceLast14days == "Y" ~ T,
                                                                          BabyAbxSinceLast14days == "N" ~ F),
                                  exclusive_breastfeeding = ExclusiveBF,
                                  read_count = NumReads,
                                  country = "USA",
                                  region = CollectionSite,
                                  study = "pannaraj_2017",
                                  lifestyle = "industrialized",
                                  family_ID = subject_ID,
                                  geographic_location_.latitude. = 25.430000,
                                  geographic_location_.longitude. = -80.160000,
                                  Instrument = "Illumina MiSeq",
                                  .keep = "unused")

metadata_pan_combined$region[metadata_pan_combined$region == "FL"] <- "FLORIDA"
metadata_pan_combined$region[metadata_pan_combined$region == "CA"] <- "CALIFORNIA"

# rename entries
metadata_pan_combined$birthmode[metadata_pan_combined$birthmode == "Vaginal"] <- "v"
metadata_pan_combined$birthmode[metadata_pan_combined$birthmode == "C-section"] <- "c"
metadata_pan_combined$birthmode[metadata_pan_combined$birthmode == ""] <- NA
metadata_pan_combined$sex[metadata_pan_combined$sex == "F"] <- "f"
metadata_pan_combined$sex[metadata_pan_combined$sex == "M"] <- "m"
metadata_pan_combined$sex[metadata_pan_combined$sex == ""] <- NA
metadata_pan_combined[metadata_pan_combined == "N"] <- F
metadata_pan_combined[metadata_pan_combined == "Y"] <- T

# feeding practice
metadata_pan_combined$AgeSolidIntro[metadata_pan_combined$AgeSolidIntro %in% c("NO SOL", "Unknown", "", " ")] <- "9999"  # set no solid, avoid NAs
metadata_pan_combined$AgeSolidIntro <- as.numeric(metadata_pan_combined$AgeSolidIntro)
metadata_pan_combined$AgeFormulaIntro[metadata_pan_combined$AgeFormulaIntro %in% c("NO F", "Unknown", "", " ")] <- "9999"
metadata_pan_combined$AgeFormulaIntro <- gsub("<", "", metadata_pan_combined$AgeFormulaIntro) %>% as.numeric()
metadata_pan_combined$PercentBF[metadata_pan_combined$PercentBF %in% c("", "Unknown", "unknown")] <- 0
metadata_pan_combined$PercentBF <- as.numeric(metadata_pan_combined$PercentBF)
# metadata_pan_combined$formula_stopped <- NA
# metadata_pan_combined$formula_stopped[grep("formula stopped", metadata_pan_combined$Notes)] <- 
#    metadata_pan_combined$Notes[grep("formula stopped", metadata_pan_combined$Notes)]
# metadata_pan_combined$formula_stopped <- gsub("formula stopped @ ", "", metadata_pan_combined$formula_stopped)
# metadata_pan_combined$formula_stopped <- gsub("days.*", "", metadata_pan_combined$formula_stopped) %>% as.numeric()
metadata_pan_combined$food <- NA
metadata_pan_combined$food[metadata_pan_combined$BFMixFM == "FM"] <- "formula"
metadata_pan_combined$food[metadata_pan_combined$PercentBF == 100] <- "breast"
metadata_pan_combined$food[is.na(metadata_pan_combined$PercentBF)] <- "mixed_b_f"
metadata_pan_combined$food[metadata_pan_combined$age >= metadata_pan_combined$AgeFormulaIntro &
                              metadata_pan_combined$SolidsIntroduced %in% c("n", "") &
                              metadata_pan_combined$PercentBF > 0] <- "mixed_b_f"
metadata_pan_combined$food[metadata_pan_combined$age < metadata_pan_combined$AgeFormulaIntro &
                              metadata_pan_combined$SolidsIntroduced == "y" &
                              metadata_pan_combined$PercentBF > 0] <- "mixed_b_s"
metadata_pan_combined$food[metadata_pan_combined$age >= metadata_pan_combined$AgeFormulaIntro &
                              metadata_pan_combined$SolidsIntroduced == "y" &
                              metadata_pan_combined$PercentBF > 0] <- "mixed_b_f_s"
metadata_pan_combined$food[metadata_pan_combined$age >= metadata_pan_combined$AgeFormulaIntro &
                              metadata_pan_combined$SolidsIntroduced == "y" &
                              metadata_pan_combined$PercentBF == 0] <- "mixed_f_s"
metadata_pan_combined$food[metadata_pan_combined$age < metadata_pan_combined$AgeFormulaIntro &
                              metadata_pan_combined$SolidsIntroduced == "y" &
                              metadata_pan_combined$PercentBF == 0] <- "solid"
metadata_pan_combined$food[metadata_pan_combined$BFMixFM == "SOLID"] <- "solid"
metadata_pan_combined$food[is.na(metadata_pan_combined$food) &
                              metadata_pan_combined$PercentBF > 0] <- "breast"
metadata_pan_combined$food[is.na(metadata_pan_combined$food) &
                              metadata_pan_combined$age < 5] <- "nothing"
sort(colnames(metadata_pan_combined))

# move sorted out samples
sort_out_samples_pan <- unique(sort_out_samples_pan)
sort_out_files_pan <- c(paste(sort_out_samples_pan, "_1.fastq.gz", sep = ""),
                        paste(sort_out_samples_pan, "_2.fastq.gz", sep = ""),
                        paste(sort_out_samples_pan, ".fastq.gz", sep = ""))
sort_out_files_pan <- sort_out_files_pan[file.exists(sort_out_files_pan)]
file.move(sort_out_files_pan, "sorted_out/")

# move good files
final_samples_pan <- metadata_pan_combined$run_accession
final_files_pan <- c(paste(final_samples_pan, "_1.fastq.gz", sep = ""),
                     paste(final_samples_pan, "_2.fastq.gz", sep = ""),
                     paste(final_samples_pan, ".fastq.gz", sep = ""))
final_files_pan <- final_files_pan[file.exists(final_files_pan)]
file.move(final_files_pan, "fastq_files/")

# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_pan_combined$read_count, n_lines$read_counts[match(metadata_pan_combined$run_accession, n_lines$sample_ID)])  # 


# write file
write.csv(metadata_pan_combined, "metadata_pan_2017_healthy.csv")

# clean up
rm(metadata_pan, metadata_pan_combined, ncbi_dict_pan, ncbi_xml_file_pan,
   ncbi_xml_pan, multiple_samples, sort_out_files_pan, sort_out_samples_pan,
   final_samples_pan, final_files_pan, n_lines, sra_file)


# Stokholm	2018 ###############################################################
# only SraRunTable.txt available
setwd(paste(maindir, "16S/stokholm_2018/", sep = ""))
ena_file_stokholm <- read.table("filereport_read_run_PRJNA417357_tsv.txt", sep = "\t", header = T)
ena_file_stokholm <- ena_file_stokholm[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]

metadata_stokholm <- read.table("SraRunTable.txt", sep = ",", header = T)
# remove columns
metadata_stokholm[,c("BioSampleModel", "Bytes", "Center.Name", "Collection_Date",
                     "Consent", "DATASTORE.filetype", "DATASTORE.provider", "DATASTORE.region",
                     "env_biome", "env_feature", "env_material", "Env_package",
                     "geo_loc_name_country", "geo_loc_name_country_continent",
                     "geo_loc_name", "Host", "LibraryLayout", "Lat_Lon",
                     "LibrarySelection", "LibrarySource", "Organism", "Platform",
                     "ReleaseDate", "reverse_read_length..run.", "Library.Name")] <- NULL

# no duplicates
sum(duplicated(metadata_stokholm$Sample.Name))

# individuals with only one sample
multiple_samples <- unique(metadata_stokholm$Individual[duplicated(metadata_stokholm$Individual)])
sort_out_samples_stok <- metadata_stokholm$Run[!metadata_stokholm$Individual %in% multiple_samples]
metadata_stokholm <- metadata_stokholm[metadata_stokholm$Individual %in% multiple_samples,]

# age column:
metadata_stokholm$age <- NA
metadata_stokholm$age[metadata_stokholm$Timepoint == "One.Week"] <- 7
metadata_stokholm$age[metadata_stokholm$Timepoint == "One.Month"] <- 30
metadata_stokholm$age[metadata_stokholm$Timepoint == "One.Year"] <- 365
metadata_stokholm$Timepoint <- NULL

# unify colnames
metadata_stokholm %<>% mutate(run_accession = Run,
                              sample_ID = paste0(Sample.Name, "_stokholm_2018"),
                              subject_ID = paste0(Individual, "_stokholm_2018"),
                              country = "DENMARK",
                              region = "COPENHAGEN",
                              study = "stokholm_2018",
                              lifestyle = "industrialized",
                              family_ID = subject_ID,
                              geographic_location_.latitude. = 55.740813,
                              geographic_location_.longitude. = 12.544197,
                              .keep = "unused")

# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_stokholm$Bases, n_lines$read_counts[match(metadata_stokholm$run_accession, n_lines$sample_ID)])  # 


write.csv(metadata_stokholm, "metadata_stokholm_2018_healthy.csv")

# move sorted out samples
sort_out_samples_stok <- unique(sort_out_samples_stok)
sort_out_files_stok <- c(paste(sort_out_samples_stok, "_1.fastq.gz", sep = ""),
                          paste(sort_out_samples_stok, "_2.fastq.gz", sep = ""),
                          paste(sort_out_samples_stok, ".fastq.gz", sep = ""))
sort_out_files_stok <- sort_out_files_stok[file.exists(sort_out_files_stok)]
file.move(sort_out_files_stok, "sorted_out/")

# move good files
final_samples_stok <- metadata_stokholm$run_accession
final_files_stok <- c(paste(final_samples_stok, "_1.fastq.gz", sep = ""),
                       paste(final_samples_stok, "_2.fastq.gz", sep = ""),
                       paste(final_samples_stok, ".fastq.gz", sep = ""))
final_files_stok <- final_files_stok[file.exists(final_files_stok)]
file.move(final_files_stok, "fastq_files/")

# clean up
rm(ena_file_stokholm, metadata_stokholm, final_files_stok, final_samples_stok,
   multiple_samples, sort_out_files_stok, sort_out_samples_stok, n_lines)





# Gupta	2019 ###################################################################

# only SraRunTable.txt available
setwd(paste(maindir, "16S/gupta_2019/", sep = ""))
ena_file_gupta <- read.table("filereport_read_run_PRJNA543007_tsv.txt", sep = "\t", header = T)
ena_file_gupta <- ena_file_gupta[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]

metadata_gupta <- read.table("SraRunTable.txt", sep = ",", header = T)
# remove columns
metadata_gupta[,c("Assay.Type", "BioSampleModel", "Bytes", "Center.Name",
                  "Collection_Date", "Consent", "DATASTORE.filetype", "DATASTORE.provider",
                  "DATASTORE.region", "geo_loc_name_country", "geo_loc_name_country_continent",
                  "geo_loc_name", "Host", "Lat_Lon", "LibraryLayout",
                  "LibrarySelection", "LibrarySource", "Organism", "Platform", "ReleaseDate",
                  "isolation", "Library.Name", "Isolation_Source")] <- NULL

# no duplicates
sum(duplicated(metadata_gupta$Sample.Name))

# remove non fecal
sort_out_samples_gupta <- metadata_gupta$Run[metadata_gupta$Type == "Hypopharyngeal"]
metadata_gupta <- metadata_gupta[metadata_gupta$Type != "Hypopharyngeal",]
# individuals with only one sample
# no individual information, not possible

# age column:
metadata_gupta$age <- NA
metadata_gupta$age[metadata_gupta$Time == "1 week"] <- 7
metadata_gupta$age[metadata_gupta$Time == "1 month"] <- 30
metadata_gupta$age[metadata_gupta$Time == "1 year"] <- 365
metadata_gupta$Time <- NULL

# unify colnames
metadata_gupta %<>% mutate(run_accession = Run,
                              sample_ID = paste0(Sample.Name, "_gupta_2019"),
                              country = "DENMARK",
                              region = "COPENHAGEN",
                              study = "gupta_2019",
                              lifestyle = "industrialized",
                              .keep = "unused")

# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_gupta$Bases, n_lines$read_counts[match(metadata_gupta$run_accession, n_lines$sample_ID)])  # 


write.csv(metadata_gupta, "metadata_gupta_2019_healthy.csv")

# move sorted out samples
sort_out_samples_gupta <- unique(sort_out_samples_gupta)
sort_out_files_gupta <- c(paste(sort_out_samples_gupta, "_1.fastq.gz", sep = ""),
                          paste(sort_out_samples_gupta, "_2.fastq.gz", sep = ""),
                          paste(sort_out_samples_gupta, ".fastq.gz", sep = ""))
sort_out_files_gupta <- sort_out_files_gupta[file.exists(sort_out_files_gupta)]
file.move(sort_out_files_gupta, "sorted_out/")

# move good files
final_samples_gupta <- metadata_gupta$run_accession
final_files_gupta <- c(paste(final_samples_gupta, "_1.fastq.gz", sep = ""),
                       paste(final_samples_gupta, "_2.fastq.gz", sep = ""),
                       paste(final_samples_gupta, ".fastq.gz", sep = ""))
final_files_gupta <- final_files_gupta[file.exists(final_files_gupta)]
file.move(final_files_gupta, "fastq_files/")

# clean up
rm(ena_file_gupta, metadata_gupta, final_files_gupta, final_samples_gupta,
   sort_out_files_gupta, sort_out_samples_gupta, n_lines)



# Iszatt	2019 #################################################################

setwd(paste(maindir, "16S/iszatt_2019/", sep = ""))
ena_file_iszatt <- read.table("filereport_read_run_PRJEB29078_tsv.txt", sep = "\t", header = T)
ena_file_iszatt <- ena_file_iszatt[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file[,c("Instrument","Run")]
all(ena_file_iszatt$run_accession %in% sra_file$Run)
ena_file_iszatt <- merge(ena_file_iszatt, sra_file,by.x = "run_accession", by.y = "Run")


metadata_isz <- read.table("11937_20220510-095154.txt", sep = '\t', header = T)


# no duplicates
sum(duplicated(ena_file_iszatt$sample_alias))

#combine data.frames
sum(!unique(ena_file_iszatt$library_name) %in% metadata_isz$sample_name)
sum(!unique(metadata_isz$sample_name) %in% ena_file_iszatt$library_name)
nrow(ena_file_iszatt)
metadata_isz_combined <- merge(metadata_isz, ena_file_iszatt, by.x = "sample_name", by.y = "library_name")
nrow(metadata_isz_combined)

# remove mothers, preterms and control samples
sort_out_samples_isz <- metadata_isz_combined$sample_name[metadata_isz_combined$day == "not applicable"]
metadata_isz_combined <- metadata_isz_combined[metadata_isz_combined$day != "not applicable",]
sort_out_samples_isz <- c(sort_out_samples_isz, 
                          metadata_isz_combined$sample_name[metadata_isz_combined$host_life_stage %in% c("adult", "not provided")])
metadata_isz_combined <- metadata_isz_combined[!metadata_isz_combined$host_life_stage %in% c("adult", "not provided"),]
sort_out_samples_isz <- c(sort_out_samples_isz, 
                          metadata_isz_combined$sample_name[metadata_isz_combined$pretermsm %in% c("1", "")])
metadata_isz_combined <- metadata_isz_combined[!metadata_isz_combined$pretermsm %in% c("1", ""),]

# # remove low read counts
# sort_out_samples_isz <- c(sort_out_samples_isz,
#                           metadata_isz_combined$sample_name[metadata_isz_combined$read_count < 100])
# metadata_isz_combined <- metadata_isz_combined[metadata_isz_combined$read_count >= 100,]


# remove all individuals with only one sample
multiple_samples <- unique(metadata_isz_combined$host_subject_id_original[duplicated(metadata_isz_combined$host_subject_id_original)])
sort_out_samples_isz <- c(sort_out_samples_isz, 
                          metadata_isz_combined$sample_name[!metadata_isz_combined$host_subject_id_original %in% multiple_samples])
metadata_isz_combined <- metadata_isz_combined[metadata_isz_combined$host_subject_id_original %in% multiple_samples,]


# remove cols
metadata_isz_combined[,c("sample_name", "collection_timestamp", "description",
                         "dna_extracted", "elevation", "empo_1", "empo_2", "empo_3",
                         "env_biome", "env_feature", "env_material", "env_package",
                         "host_age_units", "host_body_habitat",
                         "host_body_mass_index", "host_body_product", "host_body_site",
                         "host_common_name", "host_height", "host_height_units",
                         "host_life_stage", "host_scientific_name", "host_subject_id",
                         "host_taxid", "host_weight", "host_weight_units",
                         "m1smok4gi", "m2smok4g", "physical_specimen_location",
                         "physical_specimen_remaining", "qiita_study_id", "sample_type",
                         "scientific_name", "sex_numeric", "taxon_id", "title",
                         "pretermsm", "day")] <- NULL

# unify colnames
metadata_isz_combined %<>% mutate(sample_ID = anonymized_name,
                                  subject_ID = paste0(host_subject_id_original, "_iszatt_2019"),
                                  country = toupper(geo_loc_name),
                                  age = host_age,
                                  antibiotics_one_week_before = F,
                                  region = country,
                                  study = "iszatt_2019",
                                  lifestyle = "industrialized",
                                  family_ID = mother_baby_pair,
                                  geographic_location_.latitude. = latitude,
                                  geographic_location_.longitude. = longitude, 
                                  .keep = "unused")

metadata_isz_combined$sex[metadata_isz_combined$sex == "female"] <- "f"
metadata_isz_combined$sex[metadata_isz_combined$sex == "male"] <- "m"
metadata_isz_combined$age <- as.numeric(metadata_isz_combined$age)

metadata_isz_combined$sample_ID <- paste0("11937.", metadata_isz_combined$sample_ID)

metadata_isz_combined[metadata_isz_combined == ""] <- NA

# move sorted out samples
# sort_out_samples_isz <- unique(sort_out_samples_isz)
sort_out_files_isz <- list.files(pattern = paste0(sort_out_samples_isz[1:trunc(length(sort_out_samples_isz)/2)], collapse = "_|"))
sort_out_files_isz <- c(sort_out_files_isz, 
                        list.files(pattern = paste0(sort_out_samples_isz[trunc(length(sort_out_samples_isz)/2):length(sort_out_samples_isz)], collapse = "_|")))
sort_out_files_isz <- unique(sort_out_files_isz)
file.move(sort_out_files_isz, "sorted_out/")

# move good files
final_samples_isz <- metadata_isz_combined$sample_ID
final_files_isz <- list.files(pattern = paste0(final_samples_isz[1:trunc(length(final_samples_isz)/3)], collapse = "_|"))
final_files_isz <- c(final_files_isz, 
                     list.files(pattern = paste0(final_samples_isz[trunc(length(final_samples_isz)/3+1):(length(final_samples_isz)/(3/2))], collapse = "_|")))
final_files_isz <- c(final_files_isz, 
                     list.files(pattern = paste0(final_samples_isz[trunc(length(final_samples_isz)/(3/2)+1):length(final_samples_isz)], collapse = "_|")))
final_files_isz <- unique(final_files_isz)
file.move(final_files_isz, "fastq_files/")


# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_isz_combined$read_count, n_lines$read_counts[match(metadata_isz_combined$sample_ID,
                                                                 n_lines$sample_ID)])
metadata_isz_combined$read_count <- n_lines$read_counts[match(metadata_isz_combined$sample_ID,
                                                              n_lines$sample_ID)]

# write file
write.csv(metadata_isz_combined, "metadata_iszatt_2019_healthy.csv")


# clean up
rm(ena_file_iszatt, metadata_isz, metadata_isz_combined, multiple_samples, sort_out_files_isz,
   sort_out_samples_isz, final_samples_isz, final_files_isz, n_lines, sra_file)




# Gehrig	2019 #################################################################

setwd(paste(maindir, "16S/gehrig_2019/", sep = ""))
ena_file_gehr <- read.table("filereport_read_run_PRJEB26419_tsv_human.txt", sep = "\t", header = T)
ena_file_gehr <- ena_file_gehr[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]

# 54 malnurished Children with SAM + treatment
# Table S1 B per subject info 54 individuals
metadata_gehr_mal_sam_1 <- read.table("metadata_malnurished_s1b.txt", sep = "\t", header = T, na.strings = c("NA", "-", "")) %>% 
  mutate(subject_ID = gsub("-", ".", PID),
         health = "SAM",
         .keep = "unused")
# table S7 minus samples from subr, samples from the first tp, 29 samples + 9 from subr
metadata_gehr_age_mal <- read.table("metadata_malnurished_age_s7.txt", sep = "\t", header = T, na.strings = c("NA", "-", "")) %>% 
  select(Child.ID, Sample.ID, Age.of.donor.at.time.of.fecal.sample.collection..months., WHZ, WAZ, HAZ) %>%
  mutate(subject_ID = gsub("-", ".", Child.ID),
         sample_ID = Sample.ID,
         age = Age.of.donor.at.time.of.fecal.sample.collection..months.,
         health = "SAM",
         .keep = "unused")
all(metadata_gehr_age_mal$subject_ID %in% metadata_gehr_mal_sam_1$subject_ID)
# table 22 SAM, first time point only, 30 samples, should be the same as in table S1
# unclear about the samples not in S1 or S7, leave them out for now
# metadata_gehr_age_mal_4 <- read.table("metadata_malnurished_age_s22_sam.txt", sep = "\t", header = T, na.strings = c("NA", "-", "")) %>%
#   select(-c("Trial.week", "X", "MAZ", "Plasma.SampleID", "MUAC", "Weight..kg.", "Length..cm.")) %>%
#   mutate(age = Child.age..months.,
#          subject_ID = PID,
#          sample_ID = Fecal.SampleID,
#          sex = tolower(Gender),
#          health = "SAM",
#          .keep = "unused")

# samples from the first tp, table S21f, only sam, only first time points, 74, not related to sequencing
# metadata_gehr_age_mal_2 <- read.table("metadata_malnurished_s21f_sam.txt", sep = "\t", header = T, na.strings = c("NA", "-", "")) %>%
#   mutate(age = Age..months.,
#          subject_ID = paste0("PS.", sprintf("%03d", PID)),
#          sample_ID = SampleID,
#          health = "SAM",
#          .keep = "unused")


# table S22 MAM
metadata_gehr_age_mal_3 <- read.table("metadata_malnurished_age_s22_mam.txt", sep = "\t", header = T, na.strings = c("NA", "-", "")) %>%
  select(-Plasma.SampleID, -MUAC, -Timepoint) %>%
  filter(Trial.week <= 5) %>%   # only before treatment
  mutate(weight = Weight..kg.,
         length = Length..cm.,
         age = Child.age..months.,
         sample_ID = Fecal.SampleID,
         subject_ID = gsub("..G", "", sample_ID),
         health = "MAM",
         .keep = "unused")

# metadata age healthy table S4 for RF model
metadata_gehr_age <- read.table("metadata_age_s4.txt", sep = "\t", header = T, na.strings = c("NA", "-", "")) %>%
  select(Sample.ID, Age.of.donor.at.time.of.fecal.sample.collection..months., WHZ, WAZ, HAZ) %>%
  mutate(sample_ID = Sample.ID,
         subject_ID = gsub("\\.m.*", "", sample_ID),
         age = Age.of.donor.at.time.of.fecal.sample.collection..months.,
         health = "healthy",
         .keep = "unused")

# healthy and sick shotgun sequencing results
metadata_gehr_age_2 <- read.table("metadata_age_shotgun_s5b.txt", sep = "\t", header = T, na.strings = c("NA", "-", "")) %>%
  select(Sample.ID, Age..months.) %>%
  mutate(sample_ID = Sample.ID,
         subject_ID = ifelse(grepl("Bgsng|BG1C", sample_ID), yes = gsub("\\.m.*", "", sample_ID),
                          no = gsub("\\.S[014].*", "", sample_ID)),
         age = Age..months.,
         health = case_when(grepl("BG1C", sample_ID) ~ "healthy",
                            grepl("Bgsng", sample_ID) ~ "healthy",
                            grepl("PS.", sample_ID) ~ "SAM"),
         .keep = "unused") %>%
  filter(!grepl("S06|S07|S08|S09|S10|S11|S12|S13|S14", sample_ID))


metadata_gehr_age_combined <- bind_rows(metadata_gehr_age, metadata_gehr_age_mal,
                                        metadata_gehr_age_mal_3) %>%
  left_join(metadata_gehr_mal_sam_1 %>% mutate(sex = tolower(Gender)) %>% select(sex, subject_ID)) %>%
  full_join(.,metadata_gehr_age_2, by = "sample_ID") %>%
  mutate(age = ifelse(is.na(age.x), yes = age.y, no = age.x),
         health = ifelse(is.na(health.x), yes = health.y, no = health.x),
         subject_ID = ifelse(is.na(subject_ID.x), yes = subject_ID.y, no = subject_ID.x),
         .keep = "unused")

metadata_gehr_sra <- read.table("SraRunTable.txt", sep = ",", header = T)
# remove some cols:
metadata_gehr_sra[,c("Assay.Type", "AvgSpotLen", "Bases", "BioProject", "BioSample",
                     "Bytes", "Center.Name", "Consent", "DATASTORE.filetype",
                     "DATASTORE.provider", "DATASTORE.region", "ENA.FIRST.PUBLIC..run.",
                     "ENA.FIRST.PUBLIC", "ENA.LAST.UPDATE..run.", "ENA.LAST.UPDATE",
                     "INSDC_center_alias", "INSDC_center_name", "INSDC_first_public",
                     "INSDC_last_update", "INSDC_status", "LibrarySelection",
                     "Platform", "ReleaseDate", "LibrarySource", "LibraryLayout")] <- NULL
dim(ena_file_gehr)
dim(metadata_gehr_age_combined)
sum(ena_file_gehr$sample_alias %in% metadata_gehr_age_combined$sample_ID)
metadata_gehr_combined <- left_join(ena_file_gehr, metadata_gehr_age_combined, by = c("sample_alias" = "sample_ID"))
metadata_gehr_combined <- left_join(metadata_gehr_combined, metadata_gehr_sra, by = c("run_accession" = "Run"))

table(metadata_gehr_combined$Library.Name, useNA = "always")
# remove non-human samples:
metadata_gehr_combined <- metadata_gehr_combined %>% 
  filter(!common_name %in% c("house mouse", "pig"), 
         !Organism %in% c("bacterium", "mouse gut metagenome"),
         !Library.Name %in% c("table S16C, mouse sample",
                             "table S15B, mouse fecal sample")) %>%
  mutate(health = case_when(!is.na(health) ~ health,
                            is.na(health) & Library.Name == "table S23, MDCF trial" ~ "MAM",
                            is.na(health) & Library.Name == "table S5A, SAM trial" ~ "SAM"))
  
# differentiate between shotgun and 16S:
metadata_gehr_combined <- metadata_gehr_combined[-grep("metagenomic", metadata_gehr_combined$sample_alias),]
# remove raman samples:
# samples_sub <- read.table("../subramanian_2014/metadata_subramanian_2014_healthy.csv", sep = ",", header = T)
samples_ram <- read.table("../raman_2019/metadata_ram_2019_healthy.csv", sep = ",", header = T)
# sort_out_samples <- c(sort_out_samples,
#                       metadata_gehr_combined$run_accession[metadata_gehr_combined$sample_alias %in% samples_sub$sample_ID])
# metadata_gehr_combined <- metadata_gehr_combined[!metadata_gehr_combined$sample_alias %in% samples_sub$sample_ID,]
metadata_gehr_combined <- metadata_gehr_combined[!paste0(metadata_gehr_combined$sample_alias, "_raman_2019") %in% samples_ram$sample_ID,]

# clean up cols
metadata_gehr_combined[,c("sample_title", "library_name", "Barcode.sequence", "Submitter_Id", 
                          "Sample_Name", "Organism", "common_name", "Number.of.V4.16S.rDNA.reads")] <- NULL

# no duplicated samples
sum(duplicated(metadata_gehr_combined$sample_alias))

# add antibiotics information from subramanian
abx_subr_data <- read.table("suppl_2_subr.txt", sep = "\t", header = T) %>%
  select(Fecal.Sample.ID, Breast.Milk, Formula1, Solid.Foods2, diarrhoea.at.the.time.of.sample.collection3,
         Antibiotics.within.7.days.prior.to.sample.collection)
sum(metadata_gehr_combined$sample_alias %in% abx_subr_data$Fecal.Sample.ID)
metadata_gehr_combined <- left_join(metadata_gehr_combined, abx_subr_data, by = join_by(sample_alias == Fecal.Sample.ID))

# remove diarrhoea samples
metadata_gehr_combined <- metadata_gehr_combined %>%
  filter(is.na(diarrhoea.at.the.time.of.sample.collection3) | diarrhoea.at.the.time.of.sample.collection3 == "No")

# unify colnames
metadata_gehr_combined %<>% mutate(sample_ID = paste0(sample_alias, "_gehrig_2019"),
                                  subject_ID = gsub("\\.m.*", "", sample_alias) %>%
                                    gsub("\\.S.*", "", .) %>%
                                    paste0(., "_gehrig_2019"),
                                  country = "BANGLADESH",
                                  age = age * 30,
                                  antibiotics_one_week_before = ifelse(Antibiotics.within.7.days.prior.to.sample.collection == "Yes", yes = T, no = F),
                                  food = case_when(Breast.Milk == "Yes" & Formula1 == "No" & Solid.Foods2 == "No" ~ "breast",
                                                   Breast.Milk == "Yes" & Formula1 == "Yes" & Solid.Foods2 == "No" ~ "mixed_b_f",
                                                   Breast.Milk == "Yes" & Formula1 == "No" & Solid.Foods2 == "Yes" ~ "mixed_b_s",
                                                   Breast.Milk == "Yes" & Formula1 == "Yes" & Solid.Foods2 == "Yes" ~ "mixed_b_f_s",
                                                   Breast.Milk == "No" & Formula1 == "Yes" & Solid.Foods2 == "No" ~ "formula",
                                                   Breast.Milk == "No" & Formula1 == "Yes" & Solid.Foods2 == "Yes" ~ "mixed_f_s",
                                                   Breast.Milk == "No" & Formula1 == "No" & Solid.Foods2 == "Yes" ~ "solid",
                                                   .default = NA),
                                  
                                  region = "DHAKA",
                                  study = "gehrig_2019",
                                  lifestyle = "non_industrialized",
                                  family_ID = subject_ID,
                                  geographic_location_.latitude. = 23.933,
                                  geographic_location_.longitude. = 88.983,
                                  .keep = "unused") %>%
  filter(!is.na(age), !is.na(subject_ID), !is.na(sample_ID))

# # add height and weight data:
# metadata_gehr_combined$height <- sitar::LMS2z(metadata_gehr_combined$age/365, y = metadata_gehr_combined$HAZ, sex = metadata_gehr_combined$sex,
#                                       measure = "ht", ref = who06, toz = F)
# metadata_gehr_combined$weight <- sitar::LMS2z(metadata_gehr_combined$age/365, y = metadata_gehr_combined$WAZ, sex = metadata_gehr_combined$sex,
#                                       measure = "wt", ref = who06, toz = F) 


# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_gehr_combined$read_count, n_lines$read_counts[match(metadata_gehr_combined$run_accession, n_lines$sample_ID)])  # 

write.csv(metadata_gehr_combined %>% filter(health == "healthy"), "metadata_gehrig_2019_healthy.csv")
write.csv(metadata_gehr_combined %>% filter(health %in% c("SAM", "MAM")), "metadata_gehrig_2019_condition.csv")

# # move sorted out samples
malnurished_samples <- metadata_gehr_combined %>% filter(health == "malnurished") %>% pull(run_accession)
malnurished_files <- c(paste(malnurished_samples, "_1.fastq.gz", sep = ""),
                    paste(malnurished_samples, "_2.fastq.gz", sep = ""),
                    paste(malnurished_samples, ".fastq.gz", sep = ""))
malnurished_files <- malnurished_files[file.exists(malnurished_files)]
file.move(malnurished_files, "malnurished_fastq_files/")

# move good files from shotgun folder
final_samples_gehr <- metadata_gehr_combined %>% filter(is.na(health)) %>% pull(run_accession)
final_files_gehr <- c(paste0(final_samples_gehr, "_1.fastq.gz"),
                     paste0(final_samples_gehr, "_2.fastq.gz"),
                     paste0(final_samples_gehr, ".fastq.gz"))
final_files_gehr <- final_files_gehr[file.exists(final_files_gehr)]
file.move(final_files_gehr, "fastq_files/")


# clean up
rm(ena_file_gehr, metadata_gehr_age, metadata_gehr_combined, metadata_gehr_sra,
   samples_ram, final_files_gehr, final_samples_gehr, multiple_samples,
   n_lines, sra_file, sort_out_files, sort_out_samples, abx_subr_data,
   malnurished_files, malnurished_samples, metadata_gehr_age_2, metadata_gehr_age_combined,
   metadata_gehr_age_mal, metadata_gehr_age_mal_2, metadata_gehr_age_mal_3, metadata_gehr_age_mal_4,
   all_samples, dupl_samples, metadata_gehr_mal_sam_1)




# Lim	2015 #####################################################################
setwd(paste(maindir, "16S/lim_2015/", sep = ""))
ena_file_lim <- read.table("filereport_read_run_PRJNA284162_tsv.txt", sep = "\t", header = T)
ena_file_lim <- ena_file_lim[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file[,c("Instrument","Run")]
all(ena_file_lim$run_accession %in% sra_file$Run)
ena_file_lim <- merge(ena_file_lim, sra_file,by.x = "run_accession", by.y = "Run")


metadata_lim <- read.table("metadata.txt", header = T, sep = "\t")
metadata_combined_lim <- merge(ena_file_lim, metadata_lim, by.x = "sample_alias", by.y = "sample_ID", all = T)
metadata_combined_lim$sample_ID <- metadata_combined_lim$sample_alias
metadata_combined_lim <- separate(metadata_combined_lim, col = "sample_alias", into = c("subject_ID", NA, NA), sep = "[_]")
metadata_combined_lim <- metadata_combined_lim[metadata_combined_lim$library_name == "",]
metadata_combined_lim$library_name <- NULL


metadata_combined_lim %<>% mutate(sample_ID = paste0(sample_ID, "_lim_2015"),
                                  subject_ID = paste0(subject_ID, "_lim_2015"),
                                  antibiotics_before = ifelse(antibiotics_before == "y", yes = T, no = F),
                                  country = "USA",
                                  region = "MISSOURI",
                                  study = "lim_2015",
                                  lifestyle = "industrialized",
                                  family_ID = subject_ID,
                                  geographic_location_.latitude. = 38.627000,
                                  geographic_location_.longitude. = -90.199389,
                                  .keep = "unused")

# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_combined_lim$read_count, n_lines$read_counts[match(metadata_combined_lim$run_accession, n_lines$sample_ID)])  # 


write.csv(metadata_combined_lim, "metadata_lim_2015_healthy.csv")

sort_out_samples_lim <- ena_file_lim$run_accession[!ena_file_lim$run_accession %in% metadata_combined_lim$run_accession]
# move sorted out samples
sort_out_samples_lim <- unique(sort_out_samples_lim)
sort_out_files_lim <- c(paste(sort_out_samples_lim, "_1.fastq.gz", sep = ""),
                        paste(sort_out_samples_lim, "_2.fastq.gz", sep = ""),
                        paste(sort_out_samples_lim, ".fastq.gz", sep = ""))
sort_out_files_lim <- sort_out_files_lim[file.exists(sort_out_files_lim)]
file.move(sort_out_files_lim, "sorted_out/")



final_samples_lim <- metadata_combined_lim$run_accession
final_files_lim <- c(paste0(final_samples_lim, "_1.fastq.gz"),
                      paste0(final_samples_lim, "_2.fastq.gz"),
                      paste0(final_samples_lim, ".fastq.gz"))
final_files_lim <- final_files_lim[file.exists(final_files_lim)]
file.move(final_files_lim, "fastq_files/")


rm(ena_file_lim, metadata_lim, metadata_combined_lim, final_files_lim,
   final_samples_lim, sort_out_files_lim, sort_out_samples_lim, n_lines)

# Goffau	2022 #################################################################
setwd(paste(maindir, "16S/goffau_2022/", sep = ""))
ena_file_gof <- read.table("filereport_read_run_PRJEB28671_tsv_human.txt", sep = "\t", header = T, comment.char = "")
ena_file_gof <- ena_file_gof[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file[,c("Instrument","Run")]
all(ena_file_gof$run_accession %in% sra_file$Run)
ena_file_gof <- merge(ena_file_gof, sra_file,by.x = "run_accession", by.y = "Run")


metadata_goffau <- read.table("metadata.txt", sep = "\t", header = T)
metadata_goffau$sample_ID <- paste(metadata_goffau$patient.ID, metadata_goffau$timepoints, sep = "_")

# only use samples that have metadata and a fastq-file
metadata_gof_combined <- merge(metadata_goffau, ena_file_gof, by.x = "sample_ID", by.y = "sample_title") 

# unify colnames
colnames(metadata_gof_combined)[match(c("patient.ID", "gender", "age.at.sampling"),
                                      colnames(metadata_gof_combined))] <- c("subject_ID", "sex", "age")

metadata_gof_combined$age <- metadata_gof_combined$age * 30
metadata_gof_combined$sex[metadata_gof_combined$sex == "Female"] <- "f"
metadata_gof_combined$sex[metadata_gof_combined$sex == "Male"] <- "m"

metadata_gof_combined <- metadata_gof_combined[metadata_gof_combined$group == "placebo", ]
metadata_gof_combined <- metadata_gof_combined[metadata_gof_combined$read_count >= 100, ]
multiple_samples <- unique(metadata_gof_combined$subject_ID[duplicated(metadata_gof_combined$subject_ID)])
metadata_gof_combined <- metadata_gof_combined[metadata_gof_combined$subject_ID %in% multiple_samples,]

metadata_gof_combined[,c("age.at.enrolment", "X3agegroups.based.on.age.at.enroment",
                         "X11agegroups.based.on.age.at.enrolment", "group",
                         "X11agegroups.based.on.age.at.sampling", "X3agegroups.based.on.age.at.sampling")] <- NULL

metadata_gof_combined %<>% mutate(sample_ID = paste0(sample_ID, "_goffau_2022"),
                                  subject_ID = paste0(subject_ID, "_goffau_2022"),
                                  country = "GAMBIA",
                                  region = "GAMBIA",
                                  study = "goffau_2022",
                                  lifestyle = "non_industrialized",
                                  family_ID = subject_ID,
                                  geographic_location_.latitude. = 13.4063945,
                                  geographic_location_.longitude. = -14.273027,
                                  .keep = "unused")

duplicated(metadata_gof_combined$sample_ID) %>% sum

# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_gof_combined$read_count, n_lines$read_counts[match(metadata_gof_combined$run_accession, n_lines$sample_ID)])  # 


write.csv(metadata_gof_combined, "metadata_goffau_2022_healthy.csv")

sort_out_samples_gof <- ena_file_gof$run_accession[!ena_file_gof$run_accession %in% metadata_gof_combined$run_accession]
# move sorted out samples
sort_out_samples_gof <- unique(sort_out_samples_gof)
sort_out_files_gof <- c(paste(sort_out_samples_gof, "_1.fastq.gz", sep = ""),
                        paste(sort_out_samples_gof, "_2.fastq.gz", sep = ""),
                        paste(sort_out_samples_gof, ".fastq.gz", sep = ""))
sort_out_files_gof <- sort_out_files_gof[file.exists(sort_out_files_gof)]
file.move(sort_out_files_gof, "sorted_out/")



final_samples_gof <- metadata_gof_combined$run_accession
final_files_gof <- c(paste0(final_samples_gof, "_1.fastq.gz"),
                     paste0(final_samples_gof, "_2.fastq.gz"),
                     paste0(final_samples_gof, ".fastq.gz"))
final_files_gof <- final_files_gof[file.exists(final_files_gof)]
file.move(final_files_gof, "fastq_files/")


rm(ena_file_gof, metadata_gof_combined, metadata_goffau, final_files_gof, final_samples_gof,
   sort_out_files_gof, sort_out_samples_gof, multiple_samples, n_lines)



# Blanton	2016 #################################################################
setwd(paste(maindir, "16S/blanton_2016/", sep = ""))
ena_file_blant <- read.table("filereport_read_run_PRJEB9853_tsv.txt", sep = "\t", header = T, comment.char = "")
ena_file_blant <- ena_file_blant[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
# duplicated samples?
any(duplicated(ena_file_blant$sample_alias))
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file[,c("Instrument","Run")]
all(ena_file_blant$run_accession %in% sra_file$Run)
ena_file_blant <- merge(ena_file_blant, sra_file,by.x = "run_accession", by.y = "Run")


metadata_subj_blant <- read.table("metadata_subj.txt", sep = "\t", header = T)
metadata_sample_blant <- read.table("metadata_sample.txt", sep = "\t", header = T)
metadata_sample_blant <- metadata_sample_blant %>% dplyr::mutate(MUAC = Mid.Upper.Arm.Circumference..MUAC. , .keep = "unused")
metadata_sample_blant_2 <- read.table("metadata_smith.txt", sep = "\t", header = T)
metadata_sample_blant_3 <- read.table("metadata_twins.txt", sep = "\t", header = T) %>%
  select(-Family.ID,
         -Gender)

location_metadata <- metadata_sample_blant_2[,c("Site", "Health.category", "Twin.pair.phenotype",
                                                "Zygosity", "Person.ID", "RUTF", "Month", "Year")]
location_metadata <- location_metadata[!duplicated(location_metadata$Person.ID),]
# metadata_sample_blant_2 <- metadata_sample_blant_2 %>% 
#   mutate(Sample.ID = paste(Person.ID, gsub("h.*[.]", "", Sample.ID), sep = "."),
#          Child.ID = Person.ID,
#          Age..months. = Age..months,
#          .keep = "unused") %>%
#   filter(!Sample.ID %in% metadata_sample_blant$Sample.ID,
#           Health.category == "Healthy") %>%
#   select(-c(Site, Family.ID, Health.category, Twin.pair.phenotype,
#             Zygosity, Gender, RUTF, Month, Year))
# metadata_sample_combined <- bind_rows(metadata_sample_blant, metadata_sample_blant_2)

sum(metadata_sample_blant$Sample.ID %in% ena_file_blant$sample_title)
metadata_blant_combined <- merge(metadata_sample_blant, ena_file_blant, by.x = "Sample.ID", by.y = "sample_title") 
metadata_blant_combined <- merge(metadata_blant_combined, metadata_subj_blant, by = "Child.ID") 


metadata_blant_combined$Gender[metadata_blant_combined$Gender == "Female"] <- "f"
metadata_blant_combined$Gender[metadata_blant_combined$Gender == "Male"] <- "m"
metadata_blant_combined$Age..months. <- metadata_blant_combined$Age..months. * 30


# add additional info from smith et al
metadata_blant_combined <- merge(metadata_blant_combined, location_metadata, by.x = "Child.ID", by.y = "Person.ID", all.x = T)
metadata_from_smith <- metadata_sample_blant_2[c("Fever", "Diarrhea",
                                                "Vomiting", "Sample.ID")]
metadata_from_smith$Sample.ID <- sub('.', "", metadata_from_smith$Sample.ID) %>%
  gsub("Post", "", .) %>%
  gsub("Pre", "", .) %>%
  gsub("RUTF", "", .) %>%
  gsub("Mod", "", .) %>%
  gsub("kw", "", .) %>%
  gsub("Kw", "", .)
metadata_blant_combined$sample_alias %in% metadata_from_smith$Sample.ID 
metadata_blant_combined <- left_join(metadata_blant_combined, metadata_from_smith, by = join_by(sample_alias == Sample.ID))
metadata_blant_combined$Child.ID %in% metadata_sample_blant_3$Child.ID 
metadata_blant_combined <- left_join(metadata_blant_combined, metadata_sample_blant_3, by = "Child.ID")

# remove sick children
metadata_blant_combined <- metadata_blant_combined %>%
  filter(is.na(Fever) | Fever == 0,
         is.na(Diarrhea) | Diarrhea == 0,
         is.na(Vomiting) | Vomiting == 0)
metadata_blant_combined[,c("Fever", "Diarrhea", "Vomiting", "percent.samples.with.diarrhea.within.7.days.prior.to.sampling",
                           "Training.Validation.Set.....Subject.Allocation")] <- NULL
# remove a low read count sample:
metadata_blant_combined <- metadata_blant_combined %>%
  filter(run_accession != "ERR1256933")
# remove single sample children
multiple_samples <- unique(metadata_blant_combined$Child.ID[duplicated(metadata_blant_combined$Child.ID)])
# metadata_blant_combined <- metadata_blant_combined[metadata_blant_combined$Child.ID %in% multiple_samples,]
# retreive sorted out samples
metadata_blant_combined$run_accession[!metadata_blant_combined$Child.ID %in% multiple_samples] %>%
  paste("fastq-dump --split-e --gzip ",.) %>%
  data.frame(`#!/bin/bash` = .) %>%
  write.table(., file = "/fast/AG_Forslund/rob/studies/16S/blanton_2016///pull_single.sh",
              sep = " ", quote = F, row.names = F, col.names = "#!/bin/bash")





plot(metadata_blant_combined$read_count, metadata_blant_combined$Reads.per.Sample)

# unify colnames
metadata_blant_combined %<>% mutate(sample_ID = paste0(Sample.ID, "_blanton_2016"),
                                    subject_ID = paste0(Child.ID, "_blanton_2016"),
                                    age = Age..months.,
                                    antibiotics_any = ifelse(Fraction.of.samples.collected.where.antibiotics.had.been.consumed.within.prior.7.days == 0,
                                                             yes = F, no = T),
                                    weight = Weight..kg.,
                                    height = Height..cm.,
                                    family_ID = Family.ID,
                                    sex = Gender,
                                    region = Site,
                                    country = "MALAWI",
                                    region = "MALAWI",
                                    study = "blanton_2016",
                                    lifestyle = "non_industrialized",
                                    family_ID = subject_ID,
                                    geographic_location_.latitude. = -14.486173,
                                    geographic_location_.longitude. = 35.253304,
                                    .keep = "unused")
# # remove samples from kamngona:
# metadata_kamngona <- read.csv("/fast/AG_Forslund/rob/studies/16S/kamngona_2019/metadata_kamngona_2019_healthy.csv")
# metadata_blant_combined$run_accession %in% metadata_kamngona$run_accession


# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_blant_combined$read_count, n_lines$read_counts[match(metadata_blant_combined$run_accession, n_lines$sample_ID)])  # 


sort_out_samples_blant <- ena_file_blant$run_accession[!ena_file_blant$run_accession %in% metadata_blant_combined$run_accession]
# move sorted out samples
sort_out_samples_blant <- unique(sort_out_samples_blant)
sort_out_files_blant <- c(paste(sort_out_samples_blant, "_1.fastq.gz", sep = ""),
                        paste(sort_out_samples_blant, "_2.fastq.gz", sep = ""),
                        paste(sort_out_samples_blant, ".fastq.gz", sep = ""))
sort_out_files_blant <- sort_out_files_blant[file.exists(sort_out_files_blant)]
# file.move(sort_out_files_blant, "sorted_out/")


final_samples_blant <- metadata_blant_combined$run_accession
final_files_blant <- c(paste0(final_samples_blant, "_1.fastq.gz"),
                     paste0(final_samples_blant, "_2.fastq.gz"),
                     paste0(final_samples_blant, ".fastq.gz"))
final_files_blant <- final_files_blant[file.exists(final_files_blant)]
final_files_blant <- unique(final_files_blant)
# file.move(final_files_blant, "fastq_files/")

write.csv(metadata_blant_combined, "metadata_blanton_2016_healthy.csv")

rm(ena_file_blant, final_files_blant, final_samples_blant,
   metadata_blant_combined, metadata_sample_blant, metadata_subj_blant, multiple_samples,
   sort_out_files_blant, sort_out_samples_blant, location_metadata, metadata_sample_blant_2,
   n_lines, sra_file, metadata_from_smith, metadata_sample_blant_3)



# Song	2018 ###################################################################
setwd(paste(maindir, "16S/song_2021", sep = ""))
ena_file_song1 <- read.table("filereport_read_run_PRJEB14529_tsv.txt", sep = "\t", header = T, comment.char = "")
ena_file_song2 <- read.table("filereport_read_run_PRJEB36857_tsv_filter.txt", sep = "\t", header = T, comment.char = "")
ena_file_song3 <- read.table("filereport_read_run_PRJEB36860_tsv_filter.txt", sep = "\t", header = T, comment.char = "")
ena_file_song <- rbind(ena_file_song1, ena_file_song2, ena_file_song3)
ena_file_song <- ena_file_song[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
# duplicated samples?
any(duplicated(ena_file_song$sample_alias))

sra_file_1 <- read.table("SraRunTable_1.txt", sep = ",", header = T, comment.char = "", colClasses = "character")
sra_file_2 <- read.table("SraRunTable_2.txt", sep = ",", header = T, comment.char = "", colClasses = "character")
sra_file_3 <- read.table("SraRunTable_3.txt", sep = ",", header = T, comment.char = "", colClasses = "character")
sra_file <- bind_rows(sra_file_1, sra_file_2, sra_file_3) %>%
  mutate_all(type.convert, as.is=TRUE)

metadata_song_combined <- left_join(ena_file_song, sra_file, by = join_by(run_accession == Run))
metadata_song_combined[,c("sample_title", "Assay.Type", "Barcode..exp.", "Bytes",
                          "Center.Name", "COLLECTION_TIMESTAMP", "Consent", "DATASTORE.filetype",
                          "DATASTORE.provider", "DATASTORE.region", "dna_extracted",
                          "empo_1", "empo_2", "empo_3", "ENA.FIRST.PUBLIC..run.",
                          "ENA_first_public", "ENA.LAST.UPDATE..run.", "ENA.LAST.UPDATE",
                          "env_biome", "env_material", "geo_loc_name_country_continent",
                          "host_height_units", "host_weight_units",
                          "Library.Name", "LibraryLayout", "LibrarySelection",
                          "LibrarySource", "linker..exp.", "Organism",
                          "pcr_primers..exp.", "physical_specimen_remaining", "Platform",
                          "primer..exp.", "ReleaseDate", "create_date", "physical_specimen_location",
                          "sequencing_meth..exp.", "target_gene..exp.", "target_subfragment..exp.",
                          "tube_id", "Elevation", "External_Id", "INSDC_center_alias",
                          "INSDC_center_name", "INSDC_first_public", "INSDC_last_update",
                          "INSDC_status", "Sample_name", "Submitter_Id", "birth_mode_ms",
                          "env_feature", "Env_package", "body_site_corrected",
                          "body_site_type", "date_sampling_category_days", 
                          "date_sampling_category", "familyid", "body_site_orig",
                          "host_body_habitat", "host_body_mass_index", "host_body_product",
                          "host_height", "host_scientific_name",
                          "host_taxid", "host_weight", "center_project_name..exp.",
                          "sequencing_id..exp.", "mother_abx_2nd_trimester", "mother_abx_1st_trimester",
                          "mother_abx_3rd_trimester", "extraction_robot..exp.",
                          "extractionkit_lot..exp.", "mastermix_lot..exp.",
                          "orig_name..exp.", "plating..exp.", "primer_date..exp.",
                          "primer_plate..exp.", "processing_robot..exp.", "project_name..exp.",
                          "tm300_8_tool..exp.", "tm50_8_tool..exp.", "water_lot..exp.",
                          "well_id..exp.", "mother_abx_1st_trimester_name", "mother_abx_2nd_trimester_name",
                          "mother_abx_3rd_trimester_name", "TITLE", "barcodesequence..exp.",
                          "barcodewell..exp.",
                          "diet_2_month", "emp_status..exp.", "linkerprimersequence..exp.",
                          "mom_child", "mom_ld_abx", "mom_prenatal_abx", "run_prefix_prep",
                          "sample_summary", "sampletype", "samplewell..exp.", "StudyID",
                          "Day_of_life", "delivery", "month_of_life","Month", "abx_all_sources",
                          "antiexposedall", "diet_2", "diet_3", "Diet", "mom_courses_ofdeliveryabx",
                          "mom_prenatal_abx_class", "mom_prenatal_abx_trimester", "abx_name",
                          "abx_pmp_all", "course_num", "mom_pre_post_ld", "abx1_pmp_all_bymonth",
                          "experiment_center..exp.", "experiment_title..exp.", "public",
                          "qiita..exp.", "qiita_option_processing..exp.", "sample_center..exp.",
                          "abx_past_6_mos", "manuscript_use", "Bases", "date_sampling_category_days_continuous",
                          "run_date..exp.", "family_relationship", "runid..exp.",
                          "seeding_method", "anonymized_name", "plate_renamed..exp.",
                          "tm1000_8_tool..exp.", "mother_prenatal_gbs", "COUNTRY",
                          "geo_loc_name", "exclusive_breastfeed", "baby_sex", "library_name",
                          "sample_alias", "center_name..exp.", "version", "SRA.Study",
                          "run_lane..exp.", "sample_plate..exp.", "run_center..exp.",
                          "well_description..exp.", "primer_plate")] <- NULL

metadata_song_combined %<>% filter(Sample_Type == "feces",
                                   host_common_name == "human",
                                   host_body_site == "UBERON:feces",
                                   is.na(mom_baby) | mom_baby == "Baby",
                                   read_count > 100,
                                   Birth_mode != "CSseed",
                                   host_subject_id != "dad01",
                                   host_subject_id != "mum01")
metadata_song_combined[,c("Sample_Type", "host_common_name", "host_body_site", "mom_baby")] <- NULL

# multiple samples for each individual exist
multiple_samples <- unique(metadata_song_combined$host_subject_id[duplicated(metadata_song_combined$host_subject_id)])
metadata_song_combined <- metadata_song_combined[metadata_song_combined$host_subject_id %in% multiple_samples,]

metadata_song_combined[metadata_song_combined == "Yes"] <- T
metadata_song_combined[metadata_song_combined == "No"] <- F

#age
metadata_song_combined$age <- metadata_song_combined$real_sampling_time
metadata_song_combined$age <- metadata_song_combined$age %>%
  gsub("Birth", "0", .) %>%
  gsub("D", "", .)
metadata_song_combined$age[is.na(metadata_song_combined$age)] <- metadata_song_combined$Host_age[is.na(metadata_song_combined$age)]
metadata_song_combined$age[metadata_song_combined$age == ""] <- 
  metadata_song_combined$date_sampling[metadata_song_combined$age == ""]
month_data <- grep("M", metadata_song_combined$age)
metadata_song_combined$age[month_data] <- metadata_song_combined$age[month_data] %>%
  gsub("M", "", .) %>%
  gsub("onth_", "", .) %>%
  as.numeric() %>%
  multiply_by(30) %>%
  as.character()
metadata_song_combined$age <- as.numeric(metadata_song_combined$age)
metadata_song_combined[,c("age_in_years", "real_sampling_time", "Host_age", "host_age_units", "date_sampling")] <- NULL


metadata_song_combined <- metadata_song_combined %>%
  mutate(antibiotics_one_week_before = current_abx,
         birthmode = case_when(Birth_mode %in% c("CS", "CSseed", "CSself") ~ "c",
                               Birth_mode == "Vag" ~ "v"),
         country = geo_loc_name_country, 
         faimly_ID = case_when(!is.na(familyid_unique) ~ familyid_unique,
                               is.na(familyid_unique) ~ host_subject_id),
         subject_ID = host_subject_id,
         food = case_when(current_breast_feeding == "y" & current_formula == "y" & current_solids == "y" ~ "mixed_b_f_s",
                          current_breast_feeding == "y" & current_formula == "y" & current_solids == "n" ~ "mixed_b_f",
                          current_breast_feeding == "y" & current_formula == "n" & current_solids == "y" ~ "mixed_b_s",
                          current_breast_feeding == "y" & current_formula == "n" & current_solids == "n" ~ "breast",
                          current_breast_feeding == "n" & current_formula == "y" & current_solids == "y" ~ "mixed_f_s",
                          current_breast_feeding == "n" & current_formula == "y" & current_solids == "n" ~ "formula",
                          current_breast_feeding == "n" & current_formula == "n" & current_solids == "y" ~ "solid"),
         sex = ifelse(sex == "", yes = NA, no = sex),
         region = case_when(state == "Santiago" ~ "SANTIAGO",
                            state == "NY" ~ "NEW YORK",
                            state == "CO" ~ "COLORADO",
                            state == "PuertoRico" ~ "PUERTO RICO",
                            state == "WA" ~ "WASHINGTON",
                            state == "WashingtonDC" ~ "WASHINGTONDC",
                            state == "FL" ~ "FLORIDA",
                            state == "" ~ country,
                            is.na(state) ~ country),
         geographic_location_.latitude. = LATITUDE,
         geographic_location_.longitude. = LONGITUDE,
         batch = run_prefix..exp.,
         lifestyle = "industrialized",
         study = "song_2021", 
         .keep = "unused") %>%
  mutate(sample_ID = paste("song", subject_ID, age, sep = "_"), .keep = "all") %>%
  type.convert()


# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_song_combined$read_count, n_lines$read_counts[match(metadata_song_combined$run_accession, n_lines$sample_ID)])  # 


sort_out_samples_song <- ena_file_song$run_accession[!ena_file_song$run_accession %in% metadata_song_combined$run_accession]
# move sorted out samples
sort_out_samples_song <- unique(sort_out_samples_song)
sort_out_files_song <- c(paste(sort_out_samples_song, "_1.fastq.gz", sep = ""),
                          paste(sort_out_samples_song, "_2.fastq.gz", sep = ""),
                          paste(sort_out_samples_song, ".fastq.gz", sep = ""))
sort_out_files_song <- sort_out_files_song[file.exists(sort_out_files_song)]
file.move(sort_out_files_song, "sorted_out/")


final_samples_song <- metadata_song_combined$run_accession
final_files_song <- c(paste0(final_samples_song, "_1.fastq.gz"),
                       paste0(final_samples_song, "_2.fastq.gz"),
                       paste0(final_samples_song, ".fastq.gz"))
final_files_song <- final_files_song[file.exists(final_files_song)]
file.move(final_files_song, "fastq_files/")

write.csv(metadata_song_combined, "metadata_song_2021_healthy.csv")

rm(ena_file_song, ena_file_song1, ena_file_song2, ena_file_song3, metadata_song_combined,
   sra_file, sra_file_1, sra_file_2, sra_file_3, final_files_song, final_samples_song,
   sort_out_samples_song, sort_out_files_song, month_data, multiple_samples, n_lines)



# Combellick	2018 #############################################################



# Bogaert	2023 #################################################################
# setwd(paste(maindir, "16S/bogaert_2023", sep = ""))
# sra_file <- read.table("SraRunTable.txt", sep = ",", header = T, comment.char = "", colClasses = "character")
# sra_file[,c("Assay.Type", "AvgSpotLen", "BioSampleModel", "Bytes", "Center.Name",
#             "Collection_Date", "Consent", "DATASTORE.filetype", "DATASTORE.provider",
#             "DATASTORE.region", "geo_loc_name_country", "geo_loc_name_country_continent",
#             "geo_loc_name", "HOST", "Lat_Lon", "Library.Name", "LibraryLayout",
#             "LibrarySelection", "Organism", "Platform", "ReleaseDate", "create_date",
#             "version", "sample_id..run.", "SRA.Study", "LibrarySource")] <- NULL
# metadata_combined_bog <- sra_file %>%
#   filter(Isolation_Source == "fae",
#          !is.na(TIME)) %>%
#   mutate(run_accession = Run,
#          birthmode = ifelse(Birth_mode == "vag", yes = "v", no = "n"),
#          sex = ifelse(gender == "female", yes = "f", no = "m"),
#          age = case_when(TIME == "d0" ~ 0,
#                          TIME == "d1" ~ 1,
#                          TIME == "w1" ~ 7,
#                          TIME == "w2" ~ 14,
#                          TIME == "m1" ~ 30),
#          subject_ID = Subject,
#          sample_ID = Sample_id,
#          family_ID = Subject,
#          country = "NETHERLANDS",
#          region = "UTRECHT",
#          study = "bogaert_2023",
#          lifestyle = "industrialized",
#          .keep = "unused")
# metadata_combined_bog[,c("Isolation_Source", "niche..run.", "niche")] <- NULL
# 
# # remove samples already present in reyman 2019
# reyman_samples <- read.table("../reyman_2019/metadata_rey_2019_healthy.csv", sep = ",", header = T) %>%
#   mutate(sample_ID = gsub("_reyman_2019", "", sample_ID))
# metadata_combined_bog <- metadata_combined_bog %>%
#   filter(!(sample_ID %in% reyman_samples$sample_ID))
# 
# multiple_samples <- unique(metadata_combined_bog$subject_ID[duplicated(metadata_combined_bog$subject_ID)])
# metadata_combined_bog <- metadata_combined_bog[metadata_combined_bog$subject_ID %in% multiple_samples,]
# 
# 
# n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
# plot(metadata_combined_bog$Bases, n_lines$read_counts[match(metadata_combined_bog$run_accession, n_lines$sample_ID)])  # 
# 
# 
# sort_out_samples_bog <- sra_file$Run[!sra_file$Run %in% metadata_combined_bog$run_accession]
# # move sorted out samples
# sort_out_samples_bog <- unique(sort_out_samples_bog)
# sort_out_files_bog <- c(paste(sort_out_samples_bog, "_1.fastq.gz", sep = ""),
#                          paste(sort_out_samples_bog, "_2.fastq.gz", sep = ""),
#                          paste(sort_out_samples_bog, ".fastq.gz", sep = ""))
# sort_out_files_bog <- sort_out_files_bog[file.exists(sort_out_files_bog)]
# file.move(sort_out_files_bog, "sorted_out/")
# 
# 
# final_samples_bog <- metadata_combined_bog$run_accession
# final_files_bog <- c(paste0(final_samples_bog, "_1.fastq.gz"),
#                       paste0(final_samples_bog, "_2.fastq.gz"),
#                       paste0(final_samples_bog, ".fastq.gz"))
# final_files_bog <- final_files_bog[file.exists(final_files_bog)]
# file.move(final_files_bog, "fastq_files/")
# 
# write.csv(metadata_combined_bog, "metadata_bogaert_2023_healthy.csv")
# 
# rm(metadata_combined_bog, sra_file, n_lines, final_samples_bog, final_files_bog,
# plot(metadata_combined_bog$Bases, n_lines$read_counts[match(metadata_combined_bog sort_out_files_bog, sort_out_samples_bog, multiple_samples, reyman_samples)




# Muinck 2018 ##################################################################
setwd(paste(maindir, "16S/muinck_2018", sep = ""))
gender_data <- readxl::read_xlsx("metadata/GenderAndDOB.xlsx") %>%
  mutate(subject_ID = `Child ID`,
         sex = ifelse(Gender == "Female", "f", "m"),
         birthdate = `DOB (DD.MM.YYYY)`,
         .keep = "unused")
# read in sra files:
files <- list.files("metadata/", pattern = "SraRunTable", full.names = TRUE)
tables <- list()
for (file in files) {
  table <- read.csv(file, header = TRUE)
  tables[[file]] <- table
}
combined_table <- do.call(rbind, tables)
# clean up columns
combined_table[,c("Assay.Type", "BioSampleModel", "Bytes", "Center.Name", "Consent",
                  "DATASTORE.filetype", "DATASTORE.provider", "DATASTORE.region",
                  "geo_loc_name_country", "geo_loc_name_country_continent", "geo_loc_name",
                  "HOST", "Isolation_Source", "Library.Name", "LibraryLayout", 
                  "LibrarySelection", "LibrarySource", "Organism", "Platform", "ReleaseDate",
                  "create_date", "version")] <- NULL
rownames(combined_table) <- combined_table$Run
combined_table$subject_ID <- as.numeric(gsub("^ID(\\d+)_.*", "\\1", combined_table$Sample.Name))
metadata_muinck_combined <- left_join(combined_table, gender_data, by = "subject_ID")

metadata_muinck_combined <- metadata_muinck_combined %>%
  mutate(run_accession = Run,
         geographic_location_.latitude. = 59.91,
         geographic_location_.longitude. = 10.75,
         age = (as.Date(Collection_Date) - as.Date(birthdate)) %>% as.numeric(), 
         country = "NORWAY",
         region = "Oslo",
         lifestyle = "industrialized",
         study = "muinck_2018",
         antibiotics_one_week_before = ifelse(subject_ID == 10 & age > 22 & age < 40, T, F),
         antibiotics_before = ifelse(subject_ID == 10 & age > 22, T, F),
         antibiotics_any = ifelse(subject_ID == 10, T, F),
         food = case_when(age < first_solid_food ~ formula_or_breast,
                          age >= first_solid_food & formula_or_breast == "breast" ~ "mixed_b_s",
                          age >= first_solid_food & formula_or_breast == "formula" ~ "mixed_f_s",
                          age >= first_solid_food & formula_or_breast == "mixed_b_f" ~ "mixed_b_f_s"),
         sample_ID = paste(Sample.Name, study, sep = "_"),
         subject_ID = paste(subject_ID, study, sep = "_"))

metadata_muinck_combined[,c("Sample.Name", "Run", "Collection_Date", "birthdate",
                            "first_solid_food", "formula_or_breast")] <- NULL
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_muinck_combined$Bases, n_lines$read_counts[match(metadata_muinck_combined$run_accession, n_lines$sample_ID)])  # 

write.csv(metadata_muinck_combined, "metadata_muinck_2018_healthy.csv")

rm(metadata_muinck_combined, combined_table, gender_data, table, tables, file, files, n_lines)


# Davis ########################################################################
setwd(paste(maindir, "16S/davis_2017", sep = ""))
ena_file_davis <- read.table("filereport_read_run_PRJEB15633_tsv.txt", sep = "\t", header = T, comment.char = "",
                             colClasses = c(library_name = "character", sample_title = "character"))
ena_file_davis <- ena_file_davis[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T, comment.char = "", colClasses = "character")
sra_file[,c("Assay.Type", "AvgSpotLen", "body_habitat", "body_product", "body_site",
            "Bytes", "center_name..exp.", "Center.Name", "center_project_name..exp.",
            "COLLECTION_TIMESTAMP", "conception_season", "Consent", "DATASTORE.filetype",
            "DATASTORE.provider", "DATASTORE.region", "dna_extracted", "Elevation",
            "emp_status..exp.", "ENA.FIRST.PUBLIC..run.", "ENA.FIRST.PUBLIC", "ENA.LAST.UPDATE..run.",
            "ENA.LAST.UPDATE", "env_biome", "env_feature", "env_material", "Env_package",
            "experiment_center..exp.", "experiment_title..exp.", "geo_loc_name_country",
            "geo_loc_name_country_continent", "geo_loc_name", "host_common_name",
            "host_scientific_name", "host_taxid", "illumina_technology..exp.", "INSDC_center_alias",
            "INSDC_center_name", "INSDC_first_public", "INSDC_last_update", "INSDC_status",
            "LibraryLayout", "LibrarySelection", "LibrarySource",
            "linker..exp.", "Organism", "pcr_primers..exp.", "physical_specimen_location",
            "physical_specimen_remaining", "Platform", "primer..exp.", "public",
            "ReleaseDate", "required_sample_info_status", "run_center..exp.", "run_date..exp.",
            "create_date", "version", "run_prefix..exp.", "sample_center..exp.",
            "SAMPLE_TYPE", "sequencing_meth..exp.", "study_center..exp.",
            "target_gene..exp.", "target_subfragment..exp."
            )] <- NULL
qiita_file <- read.table("10297_20230201-070349.txt", header = T, sep = "\t", colClasses = c(sample_name = "character"))
qiita_file[,c("collection_timestamp", "country", "description", "diet", "dna_extracted",
              "elevation", "empo_1", "empo_2", "empo_3", "env_biome", "env_feature",
              "env_material", "env_package", "geo_loc_name", "host_age", "host_age_normalized_years",
              "host_age_units", "host_body_habitat", "host_body_product", "host_body_site",
              "host_common_name", "host_height", "host_height_units", "host_scientific_name",
              "host_taxid", "host_weight", "host_weight_units", "lifestage",
              "physical_specimen_location", "physical_specimen_remaining", "public",
              "qiita_study_id", "sample_type", "scientific_name", "baby_status",
              "birth_season", "calprotectin_units", "calprotectin", "host_subject_id",
              "season_after_lactation_time", "secretor_status", "sex", "sick_v_not_sick",
              "latitude", "longitude", "week")] <- NULL

metadata_combined_davis <- merge(ena_file_davis, sra_file, by.x = "run_accession", by.y = "Run")
metadata_combined_davis <- merge(metadata_combined_davis, qiita_file, by.x = "library_name", by.y = "sample_name")
metadata_combined_davis <- metadata_combined_davis %>%
  filter(sick_v_not_sick == "Not_Sick")

# keep only individuals with more than 1 sample:
multiple_samples <- unique(metadata_combined_davis$host_subject_id[duplicated(metadata_combined_davis$host_subject_id)])
# metadata_combined_davis <- metadata_combined_davis[metadata_combined_davis$host_subject_id %in% multiple_samples,]
metadata_combined_davis$run_accession[!metadata_combined_davis$host_subject_id %in% multiple_samples] %>%
  paste("fastq-dump --split-e --gzip ",.) %>%
  data.frame(`#!/bin/bash` = .) %>%
  write.table(., file = "/fast/AG_Forslund/rob/studies/16S/davis_2017/pull_single.sh",
              sep = " ", quote = F, row.names = F, col.names = "#!/bin/bash")


metadata_combined_davis <- metadata_combined_davis %>%
  mutate(age = as.numeric(Week) * 7,
         study = "davis_2017",
         subject_ID = paste(host_subject_id, study, sep = "_"),
         sample_ID = paste(sample_title, study, sep = "_"),
         family_ID = paste(host_subject_id, study, sep = "_"),
         geographic_location_.latitude. = LATITUDE,
         geographic_location_.longitude. = LONGITUDE,
         country = "GAMBIA",
         region = "GAMBIA",
         lifestyle = "non_industrialized",
         food = "breast",
         sex = ifelse(sex == "female", "f", "m"),
         .keep = "unused")



n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_combined_davis$Bases, n_lines$read_counts[match(metadata_combined_davis$run_accession, n_lines$sample_ID)])  # 


sort_out_samples_davis <- ena_file_davis$run_accession[!ena_file_davis$run_accession %in% metadata_combined_davis$run_accession]
# move sorted out samples
sort_out_samples_davis <- unique(sort_out_samples_davis)
sort_out_files_davis <- c(paste(sort_out_samples_davis, "_1.fastq.gz", sep = ""),
                        paste(sort_out_samples_davis, "_2.fastq.gz", sep = ""),
                        paste(sort_out_samples_davis, ".fastq.gz", sep = ""))
sort_out_files_davis <- sort_out_files_davis[file.exists(sort_out_files_davis)]
file.move(sort_out_files_davis, "sorted_out/")


final_samples_davis <- metadata_combined_davis$run_accession
final_files_davis <- c(paste0(final_samples_davis, "_1.fastq.gz"),
                     paste0(final_samples_davis, "_2.fastq.gz"),
                     paste0(final_samples_davis, ".fastq.gz"))
final_files_davis <- final_files_davis[file.exists(final_files_davis)]
file.move(final_files_davis, "fastq_files/")

write.csv(metadata_combined_davis, "metadata_davis_2017_healthy.csv")

rm(ena_file_davis, final_files_davis, final_samples_davis, metadata_combined_davis,
multiple_samples, n_lines, qiita_file, sort_out_files_davis, sort_out_samples_davis, sra_file)


# Jokela 2023 ##################################################################
setwd(paste(maindir, "16S/jokela_2023", sep = ""))
ena_file_jokela <- read.table("filereport_read_run_PRJEB55243_tsv.txt", sep = "\t", header = T, comment.char = "")
ena_file_jokela <- ena_file_jokela[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T, comment.char = "")
sra_file[,c("Assay.Type", "broad.scale_environmental_context", "Bytes", "Center.Name", 
            "Collection_Date", "Consent", "DATASTORE.filetype", "DATASTORE.provider",
            "DATASTORE.region", "ENA.FIRST.PUBLIC..run.", "ENA.FIRST.PUBLIC",
            "ENA.LAST.UPDATE..run.", "ENA.LAST.UPDATE", "environmental_medium",
            "geo_loc_name_country", "geo_loc_name_country_continent", "geographic_location_.country_and.or_sea.",
            "host_body_product", "Library.Name", "LibraryLayout", "LibrarySelection", 
            "LibrarySource", "local_environmental_context", "Organism", "Platform",
            "ReleaseDate", "create_date", "version", "Sequencing_method", "External_Id",
            "INSDC_center_alias", "INSDC_center_name", "INSDC_first_public",
            "INSDC_last_update", "INSDC_status", "Scientific_Name", "organism"
)] <- NULL

metadata_combined_jokela <- merge(ena_file_jokela, sra_file, by.x = "run_accession", by.y = "Run")
# only children, remove duplicate samples, keep the one with the higher read count
metadata_combined_jokela <- metadata_combined_jokela[grep("I.*", metadata_combined_jokela$host_subject_id),]
metadata_combined_jokela <- metadata_combined_jokela %>%
  group_by(sample_alias) %>%
  slice(which.max(read_count)) %>%
  ungroup()

metadata_combined_jokela <- metadata_combined_jokela %>%
  separate(sample_alias, into = c("individual", NA, "timep"), sep = "_", convert = T, remove = F)

metadata_combined_jokela %<>% mutate(
  age_unclear = Host_Age * 7,
  age = case_when(
    timep == "3w" ~ 21,
    timep == "6w" ~ 42,
    timep == "3m" ~ 90,
    timep == "6m" ~ 180,
    timep == "9m" ~ 273,
    timep == "12m" ~ 365,
    timep == "18m" ~ 546,
    timep == "24m" ~ 730
  ),
  sex = ifelse(host_sex == "female", yes = "f", no = "m"),
  birthmode = case_when(
    delivery_mode == "infant born by C-Section" ~ "c",
    delivery_mode == "infant born by vaginal delivery" ~ "v",
    delivery_mode == "infant born by vaginal delivery with intrapartum antibiotics" ~ "v"
  ),
  country = "FINNLAND",
  region = "Helsinki",
  lifestyle = "industrialized",
  study = "jokela_2023",
  subject_ID = individual,
  family_ID = individual,
  sample_ID = sample_alias,
  .keep = "unused")


# keep only individuals with more than 1 sample:
multiple_samples <- unique(metadata_combined_jokela$subject_ID[duplicated(metadata_combined_jokela$subject_ID)])
metadata_combined_jokela <- metadata_combined_jokela[metadata_combined_jokela$subject_ID %in% multiple_samples,]


n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_combined_jokela$Bases, n_lines$read_counts[match(metadata_combined_jokela$run_accession, n_lines$sample_ID)])  # 


sort_out_samples_jokela <- ena_file_jokela$run_accession[!ena_file_jokela$run_accession %in% metadata_combined_jokela$run_accession]
# move sorted out samples
sort_out_samples_jokela <- unique(sort_out_samples_jokela)
sort_out_files_jokela <- c(paste(sort_out_samples_jokela, "_1.fastq.gz", sep = ""),
                          paste(sort_out_samples_jokela, "_2.fastq.gz", sep = ""),
                          paste(sort_out_samples_jokela, ".fastq.gz", sep = ""))
sort_out_files_jokela <- sort_out_files_jokela[file.exists(sort_out_files_jokela)]
file.move(sort_out_files_jokela, "sorted_out/")


final_samples_jokela <- metadata_combined_jokela$run_accession
final_files_jokela <- c(paste0(final_samples_jokela, "_1.fastq.gz"),
                       paste0(final_samples_jokela, "_2.fastq.gz"),
                       paste0(final_samples_jokela, ".fastq.gz"))
final_files_jokela <- final_files_jokela[file.exists(final_files_jokela)]
file.move(final_files_jokela, "fastq_files/")

write.csv(metadata_combined_jokela, "metadata_jokela_2023_healthy.csv")

rm("ena_file_jokela", "final_files_jokela", "final_samples_jokela", "metadata_combined_jokela",
   "multiple_samples", "sort_out_files_jokela", "sort_out_samples_jokela", "sra_file")


# Jia	2019 #####################################################################
setwd(paste(maindir, "16S/jia_2019", sep = ""))
ena_file_jia <- read.table("filereport_read_run_PRJNA517572_tsv.txt", sep = "\t", header = T, comment.char = "")
ena_file_jia <- ena_file_jia[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T, comment.char = "")
sra_file[,c("Assay.Type", "AvgSpotLen", "BioSampleModel", "Bytes", "Center.Name",
            "collection_date", "Consent", "DATASTORE.filetype", "DATASTORE.provider",
            "DATASTORE.region", "env_broad_scale", "env_local_scale", "env_medium",
            "geo_loc_name_country", "geo_loc_name_country_continent", "geo_loc_name",
            "HOST", "lat_lon", "LibraryLayout", "LibrarySelection", "LibrarySource", 
            "Organism", "Platform", "ReleaseDate", "create_date", "version"
)] <- NULL

metadata_combined_jia <- merge(ena_file_jia, sra_file, by.x = "run_accession", by.y = "Run")
# remove preterms
metadata_combined_jia <- metadata_combined_jia %>%
  separate(library_name, into = c("subject_ID", "timepoint"), sep = "_", convert = T) %>%
  separate(subject_ID, into = c("term", NA), sep = 1, convert = T, remove = F) %>%
  filter(term == "F")


metadata_combined_jia %<>% mutate(
  # sex = ifelse(host_sex == "female", yes = "f", no = "m"),
  birthmode = "v",
  country = "CHINA",
  region = "Shanghai",
  lifestyle = "industrialized",
  study = "jia_2019",
  subject_ID = subject_ID,
  family_ID = subject_ID,
  sample_ID = sample_alias,
  geographic_location_.latitude. = 31.224361,
  geographic_location_.longitude. = 121.469170,
  .keep = "unused")

# keep only individuals with more than 1 sample:
multiple_samples <- unique(metadata_combined_jia$subject_ID[duplicated(metadata_combined_jia$subject_ID)])
metadata_combined_jia <- metadata_combined_jia[metadata_combined_jia$subject_ID %in% multiple_samples,]


n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_combined_jia$Bases, n_lines$read_counts[match(metadata_combined_jia$run_accession, n_lines$sample_ID)])  # 


sort_out_samples_jia <- ena_file_jia$run_accession[!ena_file_jia$run_accession %in% metadata_combined_jia$run_accession]
# move sorted out samples
sort_out_samples_jia <- unique(sort_out_samples_jia)
sort_out_files_jia <- c(paste(sort_out_samples_jia, "_1.fastq.gz", sep = ""),
                           paste(sort_out_samples_jia, "_2.fastq.gz", sep = ""),
                           paste(sort_out_samples_jia, ".fastq.gz", sep = ""))
sort_out_files_jia <- sort_out_files_jia[file.exists(sort_out_files_jia)]
file.move(sort_out_files_jia, "sorted_out/")


final_samples_jia <- metadata_combined_jia$run_accession
final_files_jia <- c(paste0(final_samples_jia, "_1.fastq.gz"),
                        paste0(final_samples_jia, "_2.fastq.gz"),
                        paste0(final_samples_jia, ".fastq.gz"))
final_files_jia <- final_files_jia[file.exists(final_files_jia)]
file.move(final_files_jia, "fastq_files/")

write.csv(metadata_combined_jia, "metadata_jia_2019_healthy.csv")

rm("ena_file_jia", "metadata_combined_jia", "multiple_samples", "sort_out_files_jia",
   "sort_out_samples_jia", "sra_file", final_files_jia, final_samples_jia)



# Beller 2021 ##################################################################
setwd(paste(maindir, "16S/beller_2021", sep = ""))
ena_file_beller <- read.table("filereport_read_run_PRJEB40751_tsv.txt", sep = "\t", header = T, comment.char = "")
ena_file_beller <- ena_file_beller[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T, comment.char = "")
sra_file[,c("Assay.Type", "AvgSpotLen", "Bytes", "Center.Name", "common_name", "Consent",
            "DATASTORE.filetype", "DATASTORE.provider", "DATASTORE.region", "ENA.FIRST.PUBLIC..run.",
            "ena_first_public", "ENA.LAST.UPDATE..run.", "ena_last_update", "geo_loc_name_country",
            "geo_loc_name_country_continent", "geographic_location_.country_and.or_sea.",
            "INSDC_center_name", "INSDC_first_public", "INSDC_last_update", "INSDC_status",
            "Library.Name", "LibraryLayout", "LibrarySelection", "LibrarySource",
            "Organism", "Platform", "ReleaseDate", "create_date", "version", "INSDC_center_alias"
            
)] <- NULL
metadata_combined_beller <- merge(ena_file_beller, sra_file, by.x = "run_accession", by.y = "Run")
metadata_1 <- read.table("BabyGut16S_metadata.csv", sep = ",", header = T)
metadata_1[,c("Holiday", "InfantID", "X.days")] <- NULL
metadata_1$sickness <- rowSums(metadata_1[,c("ICD_simple_respiratory", "ICD_simple_urinary", "ICD_simple_ear",
              "ICD_simple_gastrointestinal", "ICD_simple_skin", "ICD_simple_fungal", "ICD_simple_chickenpox",
              "ICD_simple_teething", "ICD_simple_reflux", "ICD_simple_constipation",
              "ICD_simple_NoroVirusInfection", "ICD_simple_CryptosporidiumInfection",
              "ICD_simple_UnspecifiedVirusInfection", "ICD_simple_Fever", "ICD_simple_GeneralSicknessSymptoms",
              "ICD_simple_Diarrea")] == "Y") > 0
metadata_combined_beller <- merge(metadata_combined_beller, metadata_1, by.x = "Sample_Name", by.y = "Sample.ID")

metadata_combined_beller %<>% mutate(
  gestational_age = case_when(
    InfantID == "S003" ~ 40,
    InfantID == "S004" ~ 40.5,
    InfantID == "S005" ~ 37,
    InfantID == "S006" ~ 40,
    InfantID == "S007" ~ 40,
    InfantID == "S009" ~ 41.5,
    InfantID == "S010" ~ 40,
    InfantID == "S011" ~ 39
    ),
  birth_height = case_when(
    InfantID == "S003" ~ 53,
    InfantID == "S004" ~ 51,
    InfantID == "S005" ~ 47,
    InfantID == "S006" ~ 50.5,
    InfantID == "S007" ~ 52.5,
    InfantID == "S009" ~ 54.5,
    InfantID == "S010" ~ 52,
    InfantID == "S011" ~ 52
  ),
  birth_weight = case_when(
    InfantID == "S003" ~ 3.19,
    InfantID == "S004" ~ 3.62,
    InfantID == "S005" ~ 2.86,
    InfantID == "S006" ~ 3.66,
    InfantID == "S007" ~ 3.85,
    InfantID == "S009" ~ 3.82,
    InfantID == "S010" ~ 4.08,
    InfantID == "S011" ~ 3.87
  ),
  age = Age,
  sex = tolower(Gender),
  birthmode = "v",
  food = case_when(
    FOOD == "BreastOnly" ~ "breast",
    FOOD == "NoSolidFood" ~ "mixed_b_f",
    FOOD == "SolidFood" ~ "mixed_b_f_s"
  ),
  antibiotics_before = case_when(
    (InfantID == "S003" & Age >= 354) ~ T,
    (InfantID == "S004" & Age >= 155) ~ T,
    (InfantID == "S010" & Age >= 214) ~ T,
    .default = F
  ),
  antibiotics_one_week_before = case_when(
    (InfantID == "S003" & Age >= 354 & Age <=367) ~ T,
    (InfantID == "S004" & Age >= 155 & Age <=168) ~ T,
    (InfantID == "S010" & Age >= 214 & Age <=228) ~ T,
    .default = F
  ),
  country = "BELGIUM",
  region = "BELGIUM",
  lifestyle = "industrialized",
  study = "beller_2021",
  sample_ID = paste(Sample_Name, study, sep = "_"),
  subject_ID = paste(InfantID, study, sep = "_"),
  family_ID = paste(InfantID, study, sep = "_"),
  geographic_location_.latitude. = 50.8823,
  geographic_location_.longitude. = 4.7138,
  .keep = "unused")

# remove sick samples:
metadata_combined_beller <- metadata_combined_beller %>%
  filter(ICD_simple_respiratory == "N",
         ICD_simple_urinary == "N",
         ICD_simple_gastrointestinal == "N",
         ICD_simple_fungal == "N",
         ICD_simple_chickenpox == "N",
         ICD_simple_constipation == "N",
         ICD_simple_NoroVirusInfection == "N",
         ICD_simple_CryptosporidiumInfection == "N",
         ICD_simple_UnspecifiedVirusInfection == "N",
         ICD_simple_Fever == "N",
         ICD_simple_GeneralSicknessSymptoms == "N",
         ICD_simple_Diarrea == "N"
         # ICD_simple_ear == "N",
         # ICD_simple_skin == "N",
         # ICD_simple_teething == "N",
         # ICD_simple_reflux == "N"
         )


metadata_combined_beller[,c("library_name", "X", "ICD_simple_respiratory", "ICD_simple_urinary",
                            "ICD_simple_gastrointestinal", "ICD_simple_fungal",
                            "ICD_simple_chickenpox", "ICD_simple_constipation",
                            "ICD_simple_NoroVirusInfection", "ICD_simple_CryptosporidiumInfection",
                            "ICD_simple_UnspecifiedVirusInfection", "ICD_simple_Fever",
                            "ICD_simple_GeneralSicknessSymptoms", "ICD_simple_Diarrea")] <- NULL

# keep only individuals with more than 1 sample:
multiple_samples <- unique(metadata_combined_beller$subject_ID[duplicated(metadata_combined_beller$subject_ID)])
metadata_combined_beller <- metadata_combined_beller[metadata_combined_beller$subject_ID %in% multiple_samples,]


n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_combined_beller$Bases, n_lines$read_counts[match(metadata_combined_beller$run_accession, n_lines$sample_ID)])  # 


sort_out_samples_beller <- ena_file_beller$run_accession[!ena_file_beller$run_accession %in% metadata_combined_beller$run_accession]
# move sorted out samples
sort_out_samples_beller <- unique(sort_out_samples_beller)
sort_out_files_beller <- c(paste(sort_out_samples_beller, "_1.fastq.gz", sep = ""),
                        paste(sort_out_samples_beller, "_2.fastq.gz", sep = ""),
                        paste(sort_out_samples_beller, ".fastq.gz", sep = ""))
sort_out_files_beller <- sort_out_files_beller[file.exists(sort_out_files_beller)]
file.move(sort_out_files_beller, "sorted_out/")


final_samples_beller <- metadata_combined_beller$run_accession
final_files_beller <- c(paste0(final_samples_beller, "_1.fastq.gz"),
                     paste0(final_samples_beller, "_2.fastq.gz"),
                     paste0(final_samples_beller, ".fastq.gz"))
final_files_beller <- final_files_beller[file.exists(final_files_beller)]
file.move(final_files_beller, "fastq_files/")

write.csv(metadata_combined_beller, "metadata_beller_2021_healthy.csv")

rm("ena_file_beller", "final_files_beller", "final_samples_beller", "metadata_1",
   "metadata_combined_beller", "sort_out_files_beller", "sort_out_samples_beller",
   "sra_file", n_lines, multiple_samples)



# Korpela	2020 #################################################################
setwd(paste(maindir, "16S/korpela_2020", sep = ""))
metadata_korpela <- readxl::read_xlsx("RSV_metadata.xlsx")

metadata_korpela <- metadata_korpela %>% 
  mutate(country = "FINNLAND",
         region = "FINNLAND",
         family_ID = subject,
         subject_ID = subject,
         sample_ID = sampleID,
         # age = ?,
         # antibiotics_any = ifelse(is.na(timeSinceLastAb.days), yes = F, no = T),
         Instrument = "Illumina MiSeq",
         geographic_location_.latitude. = 60.192059,
         geographic_location_.longitude. = 24.945831,
         run_accession = sampleID,
         study = "korpela_2020",
         lifestyle = "industrialized")
sort_out_indiv <- metadata_korpela %>%
  filter(Group == "AB") %>%
  select(subject_ID)
metadata_korpela <- metadata_korpela %>%
  filter(Group == "Ctrl")


sort_out_files_korpela <- list.files(pattern = paste0(sort_out_indiv$subject_ID, ".*.fastq.gz", collapse = "|"))
file.move(sort_out_files_korpela, "sorted_out/")

keep_files_korpela <- list.files(pattern = paste0(metadata_korpela$subject_ID, ".*.fastq.gz", collapse = "|"))
file.move(keep_files_korpela, "fastq_files/")


write.csv(metadata_korpela, "metadata_korpela_2020_healthy.csv")

rm("keep_files_korpela", "metadata_korpela", "sort_out_files_korpela", "sort_out_indiv",
   "sort_out_samples_korpela")





# Kortekangas	2020 #############################################################
setwd(paste(maindir, "16S/kortekangas_2020", sep = ""))
ena_file_kort <- read.table("filereport_read_run_PRJEB29433_tsv.txt", sep = "\t", header = T, comment.char = "")
ena_file_kort <- ena_file_kort[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T, comment.char = "")
sra_file[,c("Assay.Type", "AvgSpotLen", "Bytes", "Center.Name", "Consent", 
            "DATASTORE.filetype", "DATASTORE.provider", "DATASTORE.region",
            "ENA.FIRST.PUBLIC..run.", "ENA.LAST.UPDATE..run.", "INSDC_center_name",
            "INSDC_first_public", "INSDC_last_update", "INSDC_status", "Library.Name",
            "LibraryLayout", "LibrarySelection", "LibrarySource", "Organism",
            "Platform", "project_title", "ReleaseDate", "create_date", "version",
            "Sequencing_method", "ENA_first_public", "ENA.LAST.UPDATE", "Scientific_Name",
            "ena_first_public", "ENA_last_update", "INSDC_center_alias"
)] <- NULL

metadata_combined_kort <- merge(ena_file_kort, sra_file, by.x = "run_accession", by.y = "Run")
# remove moms and samples from blanton
metadata_combined_kort <- metadata_combined_kort %>%
  filter(previous_accession != "PRJEB9853") %>%
  filter(Family_member %in% c("Child1", "Child2")) %>%
  select(-previous_accession, -sample_type, -library_name)
  


metadata_combined_kort %<>% mutate(
  age = case_when(age_in_months == "30M" ~ 900,
                  age_in_months == "1M" ~ 30,
                  age_in_months == "6M" ~ 180,
                  age_in_months == "12M" ~ 365,
                  age_in_months == "18M" ~ 545),
  country = "MALAWI",
  region = "Mangochi",
  lifestyle = "non_industrialized",
  study = "kortekangas_2020",
  subject_ID = paste(Family_Id, Family_member, study, sep = "_"),
  family_ID = Family_Id,
  sample_ID = paste(Sample_name, sep = "_"),
  geographic_location_.latitude. = -14.46,
  geographic_location_.longitude. = 35.27,
  .keep = "unused")

# keep only individuals with more than 1 sample:
multiple_samples <- unique(metadata_combined_kort$subject_ID[duplicated(metadata_combined_kort$subject_ID)])
# metadata_combined_kort <- metadata_combined_kort[metadata_combined_kort$subject_ID %in% multiple_samples,]

metadata_combined_kort$run_accession[!metadata_combined_kort$subject_ID %in% multiple_samples] %>%
  paste("fastq-dump --split-e --gzip ",.) %>%
  data.frame(`#!/bin/bash` = .) %>%
  write.table(., file = "/fast/AG_Forslund/rob/studies/16S/kortekangas_2020//pull_single.sh",
              sep = " ", quote = F, row.names = F, col.names = "#!/bin/bash")



n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_combined_kort$Bases, n_lines$read_counts[match(metadata_combined_kort$run_accession, n_lines$sample_ID)])  # 


sort_out_samples_kort <- ena_file_kort$run_accession[!ena_file_kort$run_accession %in% metadata_combined_kort$run_accession]
# move sorted out samples
sort_out_samples_kort <- unique(sort_out_samples_kort)
sort_out_files_kort <- c(paste(sort_out_samples_kort, "_1.fastq.gz", sep = ""),
                        paste(sort_out_samples_kort, "_2.fastq.gz", sep = ""),
                        paste(sort_out_samples_kort, ".fastq.gz", sep = ""))
sort_out_files_kort <- sort_out_files_kort[file.exists(sort_out_files_kort)]
file.move(sort_out_files_kort, "sorted_out/")


final_samples_kort <- metadata_combined_kort$run_accession
final_files_kort <- c(paste0(final_samples_kort, "_1.fastq.gz"),
                     paste0(final_samples_kort, "_2.fastq.gz"),
                     paste0(final_samples_kort, ".fastq.gz"))
final_files_kort <- final_files_kort[file.exists(final_files_kort)]
file.move(final_files_kort, "fastq_files/")

write.csv(metadata_combined_kort, "metadata_kortekangas_2020_healthy.csv")
rm(ena_file_kort, metadata_combined_kort, n_lines, sra_file, final_files_kort,
   final_samples_kort, multiple_samples, sort_out_files_kort, sort_out_samples_kort)

# Bender	2016 #############################################################
setwd(paste(maindir, "16S/bender_2016", sep = ""))
ena_file_bend <- read.table("fastq-run-info.tsv", sep = "\t", header = T, comment.char = "")
ena_file_bend <- ena_file_bend[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias")]
sra_file <- read.table("SraRunTable.csv", sep = ",", header = T, comment.char = "")
sra_file[,c("alignment_software..exp.", "Assay.Type", "AvgSpotLen", "BioSampleModel", 
            "Bytes", "Center.Name", "Collection_Date", "Consent", "DATASTORE.filetype",
            "DATASTORE.provider", "DATASTORE.region", "Experiment", "geo_loc_name_country",
            "geo_loc_name_country_continent", "geo_loc_name", "HOST", "LibraryLayout",
            "LibrarySelection", "LibrarySource", "Organism", "Platform", "ReleaseDate",
            "create_date", "version")] <- NULL

metadata_combined_bend <- merge(ena_file_bend, sra_file, by.x = "run_accession", by.y = "Run")
# remove moms and samples from blanton
metadata_combined_bend <- metadata_combined_bend %>%
  separate(sample_alias, into = c("sample_ID", "isolation_source_2")) %>%
  filter(isolation_source == "stool")

metadata_combined_bend %<>% mutate(
  age = BBAge,
  antibiotics_any = F,
  antibiotics_before = F,
  antibiotics_one_week_before = F,
  country = "Haiti",
  region = "Port_au_prince",
  lifestyle = "non_industrialized",
  study = "bender_2016",
  geographic_location_.latitude. = 18.32,
  geographic_location_.longitude. = -72.20,
  subject_ID = paste(PairID, study, sep = "_"),
  family_ID = PairID,
  sample_ID = paste(Sample.Name, study, sep = "_"),
  .keep = "unused")

# keep only individuals with more than 1 sample:
# multiple_samples <- unique(metadata_combined_kort$subject_ID[duplicated(metadata_combined_kort$subject_ID)])
# metadata_combined_kort <- metadata_combined_kort[metadata_combined_kort$subject_ID %in% multiple_samples,]


# n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
# plot(metadata_combined_bend$Bases, n_lines$read_counts[match(metadata_combined_bend$run_accession, n_lines$sample_ID)])  # 


# samples to keep
final_samples_bend <- metadata_combined_bend$run_accession
final_files_bend <- c(paste0(final_samples_bend, "_1.fastq.gz"),
                      paste0(final_samples_bend, "_2.fastq.gz"),
                      paste0(final_samples_bend, ".fastq.gz"))
final_files_bend <- final_files_bend[file.exists(final_files_bend)]
file.move(final_files_bend, "fastq_files/")

write.csv(metadata_combined_bend, "metadata_bender_2016_healthy.csv")
rm(ena_file_bend, metadata_combined_bend, sra_file, final_samples_bend, final_files_bend)

# Morandini	2023 #############################################################
setwd(paste(maindir, "16S/morandini_2023", sep = ""))
ena_file_moran <- read.table("fastq-run-info.tsv", sep = "\t", header = T, comment.char = "")
ena_file_moran <- ena_file_moran[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias")]
sra_file <- read.table("SraRunTable.csv", sep = ",", header = T, comment.char = "")
sra_file[,c("Assay.Type", "AvgSpotLen", "BioSampleModel", 
            "Bytes", "Center.Name", "Collection_Date", "Consent", "DATASTORE.filetype",
            "DATASTORE.provider", "DATASTORE.region", "env_broad_scale",
            "env_local_scale", "env_medium", "Experiment", "geo_loc_name_country",
            "geo_loc_name_country_continent", "HOST", "lat_lon", "LibraryLayout",
            "LibrarySelection", "LibrarySource", "Organism", "Platform", "ReleaseDate",
            "create_date", "version")] <- NULL
metadata_morandini <- readxl::read_xlsx("1-s2.0-S2589004223022137-mmc2.xlsx")

metadata_combined_mor <- merge(metadata_morandini, sra_file, by.x = "SampleID", by.y = "Sample.Name") %>%
  left_join(., ena_file_moran, by = c("Run" = "run_accession"), keep = T) %>%
  filter(AdultChild == "Child") %>%
  filter(is.na(HealthStatusC) | HealthStatusC %in% c("Delayed cicatrization of belly button",
                                                     "Formula fed in the first 48h of life, breastfed after",
                                                     "Minor sebhorroic dermatitis",
                                                     "Occasionally fed with cow milk"))


metadata_combined_mor %<>% mutate(
  age = AgeC * 30,
  sex = ifelse(Gender == "F", yes = "f", no = "m"),
  birthmode = ifelse(cSection, yes = "c", no = "v"),
  height = Height,
  weight = Weight,
  antibiotics_any = F,
  antibiotics_one_week_before = F,
  country = "SENEGAL",
  region = ifelse(Location == "Rural", yes = "Widou Thiengoli", no = "DAKAR"),
  lifestyle = "non_industrialized",
  study = "morandini_2023",
  geographic_location_.latitude. = ifelse(Location == "Rural", yes = 15.522, no = 14.693),
  geographic_location_.longitude. = ifelse(Location == "Rural", yes = -14.01, no = -17.448),
  subject_ID = paste(host_subject_id, study, sep = "_"),
  family_ID = host_subject_id,
  sample_ID = paste(SampleID, study, sep = "_"),
  .keep = "unused")


# samples to keep
final_samples_moran <- metadata_combined_mor$Run
final_files_moran <- c(paste0(final_samples_moran, "_1.fastq.gz"),
                      paste0(final_samples_moran, "_2.fastq.gz"),
                      paste0(final_samples_moran, ".fastq.gz"))
final_files_moran <- final_files_moran[file.exists(final_files_moran)]
file.move(final_files_moran, "fastq_files/")

write.csv(metadata_combined_mor, "metadata_morandini_2023_healthy.csv")
rm(ena_file_moran, metadata_combined_mor, metadata_morandini, sra_file, final_files_moran,
   final_samples_moran)


# preterm studies #########################################

## Alcon-Giner	2020 #############
setwd(paste(maindir, "16S/alcon_giner_2020/", sep = ""))
sra_file <- read.table("SraRunTable.csv", sep = ",", header = T)
sra_file <- sra_file %>% select(Run, Assay.Type, geographic_location_.latitude.,
                                geographic_location_.longitude., Instrument, Organism,
                                Bases, Sample_name) %>%
  mutate(sample_ID = Sample_name, 
         run_accession = Run,
         country = "UNITED_KINGDOM",
         study = "alcon_giner_2020",
         lifestyle = "industrialized",
         health = "preterm",
         .keep = "unused") %>%
  filter(Organism == "metagenome")
excel_info <- readxl::read_xlsx("1-s2.0-S2666379120300987-mmc3.xlsx", na = c("", "NA"))
excel_info <- excel_info %>% select(ENA_unique_name, suspectedNEC, sex, delivery, 
                                    birthweight,	gestageweeks,	gestagedays, antibiotics,
                                    diettype,	daysofprobiotics,	NICUstay, daysfrombirth, 
                                    Treatment, BAMBIinfant) %>%
  mutate(gestagedays = ifelse(is.na(gestagedays) , yes = 0, no = as.numeric(gestagedays)), 
         sex = tolower(sex),
         sample_ID = ENA_unique_name,
         suspectedNEC = suspectedNEC == "y",
         birthmode = tolower(delivery),
         birth_weight = birthweight,
         antibiotics_any = antibiotics == "y",
         age = daysfrombirth,
         gestational_age = gestageweeks + (gestagedays / 7),
         subject_ID = BAMBIinfant,
         family_ID = BAMBIinfant,
         .keep = "unused")

metadata_combined <- left_join(sra_file, excel_info, by = "sample_ID") %>%
  mutate(sample_ID = paste(sample_ID, study, sep = "_"))

final_samples <- metadata_combined$run_accession
final_files <- c(paste0(final_samples, "_1.fastq.gz"),
                 paste0(final_samples, "_2.fastq.gz"),
                 paste0(final_samples, ".fastq.gz"))
final_files <- final_files[file.exists(final_files)]
file.move(final_files, "fastq_files_preterm/")

write.csv(metadata_combined, "metadata_alcon_giner_2020_preterm.csv")
rm(excel_info, filereport, metadata_combined, sra_file, final_files, final_samples)


## Gibson	2016 #############
setwd(paste(maindir, "16S/gibson_2016/", sep = ""))
sra_file <- read.table("SraRunTable.csv", sep = ",", header = T)
sra_file <- sra_file %>% select(Run, Assay.Type, lat_lon, Instrument, LibraryLayout,
                                Bases, BirthWeight, host_body_product, Library.Name,
                                host_day_of_life, host_sex, host_subject_id) %>%
  separate(., col = "Library.Name", into = c("type", "sample_ID", "sample_ID_2"), sep = "[_]") %>%
  mutate(age = host_day_of_life,
         birth_weight = BirthWeight,
         run_accession = Run,
         country = "USA",
         study = "gibson_2016",
         subject_ID = paste(host_subject_id, study, sep = "_"), 
         lifestyle = "industrialized",
         health = "preterm",
         .keep = "unused") %>%
  filter(Assay.Type == "AMPLICON", host_body_product == "stool")
excel_individual_info <- readxl::read_xlsx("41564_2016_BFnmicrobiol201624_MOESM249_ESM.xlsx", skip = 1,
                                na = c("", "NA"), sheet = "Table 1b - Individual Metadata")
excel_individual_info <- excel_individual_info %>% select(-c("CRIB II Score", "Number of Samples")) %>%
  mutate(sex = ifelse(Gender == "Female", yes = "f", no = "m"),
         birthmode = ifelse(`Delivery Method` == "Cesarean Section", yes = "c", no = "v"),
         antibiotics_any = `Antibiotics at Birth` == "YES",
         gestational_age = `Gestational Age (Weeks)`,
         subject_ID = paste0(Individual, "gibson_2016"),
         family_ID = Individual,
         .keep = "unused")

metadata_combined <- left_join(sra_file, excel_individual_info, by = "subject_ID") %>%
  mutate(sample_ID = paste(sample_ID, study, sep = "_"))

final_samples <- metadata_combined$run_accession
final_files <- c(paste0(final_samples, "_1.fastq.gz"),
                 paste0(final_samples, "_2.fastq.gz"),
                 paste0(final_samples, ".fastq.gz"))
final_files <- final_files[file.exists(final_files)]
file.move(final_files, "fastq_files_preterm/")

write.csv(metadata_combined, "metadata_gibson_2016_preterm.csv")
rm(excel_individual_info, metadata_combined, sra_file, final_files, final_samples)


## Rao	2021 #############
setwd(paste(maindir, "16S/rao_2021/", sep = ""))
sra_file <- read.table("SraRunTable.csv", sep = ",", header = T, quote = "\"")
metadata_combined <- sra_file %>% select(Run, Bases, common_name, Instrument, Submitter_Id) %>%
  separate(., col = "Submitter_Id", into = c("type", "sample_ID", "time_point", "kingdom1", "kingdom2"),
           sep = "[_]") %>%
  mutate(run_accession = Run,
         family_ID = sample_ID,
         country = "USA",
         study = "rao_2021",
         subject_ID = paste(sample_ID, study, sep = "_"),
         lifestyle = "industrialized",
         health = "preterm",
         sample_ID = paste(sample_ID, study, sep = "_"),
         .keep = "unused") %>%
  filter(common_name == "human gut metagenome", kingdom1 %in% c("bac16S", "bac16SV3V4"))

final_samples <- metadata_combined$run_accession
final_files <- c(paste0(final_samples, "_1.fastq.gz"),
                 paste0(final_samples, "_2.fastq.gz"),
                 paste0(final_samples, ".fastq.gz"))
final_files <- final_files[file.exists(final_files)]
file.move(final_files, "fastq_files_preterm/")

write.csv(metadata_combined, "metadata_rao_2021_preterm.csv")
rm(metadata_combined, sra_file, final_files, final_samples)


## Seki	2021 #####################
# no age information
# setwd(paste(maindir, "16S/seki_2021/", sep = ""))
# sra_file <- read.table("SraRunTable.csv", sep = ",", header = T)
# metadata_combined <- sra_file %>% select(Run, Assay.Type, Bases, env_medium,
#                                          lat_lon, Instrument, Library.Name, host_sex, 
#                                          host_subject_id, GA, is_healthy) %>%
#   separate(., col = "Submitter_Id", into = c("type", "sample_ID", "time_point", "kingdom1", "kingdom2"),
#            sep = "[_]") %>%
#   mutate(run_accession = Run,
#          sex = ifelse(host_sex == "female", yes = "f", no = "m"),
#          subjsect_ID = sample_ID,
#          family_ID = sample_ID,
#          country = "AUSTRIA",
#          study = "seki_2021",
#          lifestyle = "industrialized",
#          health = "preterm",
#          sample_ID = paste(sample_ID, study, sep = "_"),
#          .keep = "unused") %>%
#   filter(env_medium == "stool", kingdom1 %in% c("bac16S", "bac16SV3V4"))


## Ryan	2019 ############################
setwd(paste(maindir, "16S/ryan_2019/", sep = ""))
sra_file <- read.table("SraRunTable.csv", sep = ",", header = T)
sra_file <- sra_file %>% select(Run, Instrument,
                                Bases, Sample.Name) %>%
  separate(., col = "Sample.Name", into = c(NA, "subject_ID", NA, "age"),
           sep = "[_]", remove = F) %>%
  mutate(sample_ID = Sample.Name, 
         subject_ID = paste("ryan_2019", subject_ID, sep = "_"),
         family_ID = subject_ID,
         run_accession = Run,
         country = "AUSTRALIA",
         study = "ryan_2019",
         lifestyle = "industrialized",
         health = "preterm",
         age = gsub("day", "", age) %>% as.numeric(),
         .keep = "unused")
final_samples <- sra_file$run_accession
final_files <- c(paste0(final_samples, "_1.fastq.gz"),
                 paste0(final_samples, "_2.fastq.gz"),
                 paste0(final_samples, ".fastq.gz"))
final_files <- final_files[file.exists(final_files)]
file.move(final_files, "fastq_files_preterm/")

write.csv(sra_file, "metadata_ryan_2019_preterm.csv")
rm(sra_file, final_files, final_samples)


## Kamdar	2020 ###########################
setwd(paste(maindir, "16S/kamdar_2020/", sep = ""))
sra_file <- read.table("SraRunTable.csv", sep = ",", header = T)
sra_file <- sra_file %>% select(Run, seq_methods, Description, lat_lon,
                                Bases, Sample.Name) %>%
  filter(Description == "Preterm faecal sample") %>%
  mutate(Instrument = seq_methods,
         sample_ID = Sample.Name,
         run_accession = Run,
         country = "UNITED_KINGDOM",
         study = "kamdar_2020",
         lifestyle = "industrialized",
         health = "preterm",
         .keep = "unused")

excel_individual_info <- readxl::read_xlsx("41467_2020_14923_MOESM4_ESM.xlsx", skip = 2,
                                             na = c("", "NA"), sheet = "Fig 4A") %>%
  select(`Sample ID`, Participant, Day) %>%
  mutate(age = Day,
         subject_ID = Participant,
         family_ID = Participant,
         .keep = "unused")
metadata_combined <- left_join(sra_file, excel_individual_info, by = c("sample_ID" = "Sample ID")) %>%
  mutate(sample_ID = paste(sample_ID, "kamdar_2020", sep = "_"),
         subject_ID = paste(subject_ID, "kamdar_2020", sep = "_")) %>%
  filter(!is.na(age))

write.csv(metadata_combined, "metadata_kamdar_2020_preterm.csv")

final_samples <- metadata_combined$run_accession
final_files <- c(paste0(final_samples, "_1.fastq.gz"),
                 paste0(final_samples, "_2.fastq.gz"),
                 paste0(final_samples, ".fastq.gz"))
final_files <- final_files[file.exists(final_files)]
file.move(final_files, "fastq_files_preterm/")
rm(excel_individual_info, sra_file, metadata_combined, final_files, final_samples)

# Yuan	2019 ###############################
# setwd(paste(maindir, "16S/yuan_2019/", sep = ""))
# sra_file <- read.table("SraRunTable.csv", sep = ",", header = T)
# sra_file_2 <- read.table("SraRunTable_PRJNA493729.csv", sep = ",", header = T)
# 
# sra_file <- bind_rows(sra_file, sra_file_2) %>% select(Run, Instrument, Library.Name,
#                                 Bases, Sample.Name) %>%
#   mutate(sample_ID = Sample.Name, 
#          run_accession = Run,
#          country = "CHINA",
#          study = "yuan_2019",
#          lifestyle = "industrialized",
#          health = "preterm",
#          .keep = "unused")
# 
# excel_individual_info <- readxl::read_xlsx("41467_2020_14923_MOESM4_ESM.xlsx", skip = 2,
#                                            na = c("", "NA"), sheet = "Fig 4A") %>%
#   select(`Sample ID`, Participant, Day) %>%
#   mutate(age = Day,
#          subject_ID = Participant,
#          family_ID = Participant,
#          .keep = "unused")
# metadata_combined <- left_join(sra_file, excel_individual_info, by = c("sample_ID" = "Sample ID")) %>%
#   filter(!is.na(age))
# 
# final_samples <- metadata_combined$run_accession
# final_files <- c(paste0(final_samples, "_1.fastq.gz"),
#                  paste0(final_samples, "_2.fastq.gz"),
#                  paste0(final_samples, ".fastq.gz"))
# final_files <- final_files[file.exists(final_files)]
# file.move(final_files, "fastq_files_preterm/")
# rm(excel_individual_info, sra_file, metadata_combined, final_files, final_samples)


# Gregory	2016 ############################
setwd(paste(maindir, "16S/gregory_2016/", sep = ""))
sra_file <- read.table("SraRunTable.csv", sep = ",", header = T) %>%
  select(Run, Instrument, Host_age, Bases, Sample.Name, host_birthweight_.grams.,
         host_diet, host_sex, lat_lon, Library.Name) %>%
  separate(., Library.Name, into = c("subject_ID", "time_point"), sep = "_", remove = F) %>%
  mutate(birthweight = host_birthweight_.grams.,
         sex = ifelse(host_sex == "female", yes = "f", no = "m"),
         family_ID = subject_ID,
         run_accession = Run,
         country = "USA",
         study = "gregory_2016",
         sample_ID = paste(Library.Name, study, sep = "_"),
         subject_ID = paste(sample_ID, study, sep = "_"),
         lifestyle = "industrialized",
         health = "preterm",
         .keep = "unused")

write.csv(sra_file, "metadata_gregory_2016_preterm.csv")

final_samples <- sra_file$run_accession
final_files <- c(paste0(final_samples, "_1.fastq.gz"),
                 paste0(final_samples, "_2.fastq.gz"),
                 paste0(final_samples, ".fastq.gz"))
final_files <- final_files[file.exists(final_files)]
file.move(final_files, "fastq_files_preterm/")
rm(sra_file, final_files, final_samples)



# shotgun studies ##############################################################
################################################################################



# Ferretti	2018 ###############################################################
# C10039 with two times t0
# C10019 with probably one t1 mislableled as mother
# no c-section birth

setwd(paste(maindir, "shotgun/ferretti_2018/", sep = ""))
ena_file_fer <- read.table("filereport_read_run_PRJNA352475_tsv.txt", sep = "\t", header = T)
ena_file_fer <- ena_file_fer[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file[,c("Instrument","Run")]
all(ena_file_fer$run_accession %in% sra_file$Run)
ena_file_fer <- merge(ena_file_fer, sra_file,by.x = "run_accession", by.y = "Run")


metadata_fer <- read.table("r_metadata.csv", sep = ',', header = T)
metadata_fer <- separate(metadata_fer, col = "subjectID", into = c("subjectID", "m"), sep = "[_]")
# for one sample two mother samples are available, keep both to check, if any of them is a child
metadata_fer$m[metadata_fer$subjectID == "C10039"] <- NA
# for one sample, m and family_role disagree, keep this sample (child in m-column)
metadata_fer$family_role[is.na(metadata_fer$m)] <- "child"
metadata_fer$age[is.na(metadata_fer$m)] <- 0
metadata_fer$BMI[is.na(metadata_fer$m)] <- NA
metadata_fer$feeding_practice[metadata_fer$X == "CA_C10019IS2318FE_t1M15"] <- "exclusively_breastfeeding"
# for one sample two mother samples are available, keep both to check, if any of them is a child


metadata_fer_2 <- read.table("suppl_2.txt", sep = "\t", header = T)
metadata_fer_2[metadata_fer_2 == "na"] <- NA
metadata_fer_2$Couple.ID <- paste("C", metadata_fer_2$Couple.ID, sep = "")
metadata_fer_2$Weight.at.birth..kg. <- gsub(",", ".", metadata_fer_2$Weight.at.birth..kg.) %>% as.numeric()
metadata_fer_2$Weight.at.3.days..kg. <- gsub(",", ".", metadata_fer_2$Weight.at.3.days..kg.) %>% as.numeric()
metadata_fer_2[c("X", "X.1", "X.2", "X.3")] <- NULL

# no duplicates
sum(duplicated(ena_file_fer$sample_alias))

#combine data.frames
sum(!unique(ena_file_fer$sample_alias) %in% metadata_fer$X) # non fecal samples are not in the metadata
sum(!unique(metadata_fer$X) %in% ena_file_fer$sample_alias) # files for all metadata exist
nrow(ena_file_fer)
metadata_fer_combined <- merge(metadata_fer, ena_file_fer, by.x = "X", by.y = "sample_alias", all.y = T)
nrow(metadata_fer_combined)

sum(!unique(metadata_fer_combined$subjectID) %in% metadata_fer_2$Couple.ID) # 1 caused by NAs
sum(!unique(metadata_fer_2$Couple.ID) %in% metadata_fer_combined$subjectID) # files for all metadata exist
nrow(metadata_fer_combined)
metadata_fer_combined <- merge(metadata_fer_combined, metadata_fer_2, by.x = "subjectID", by.y = "Couple.ID", all.x = T)
nrow(metadata_fer_combined)


# multiple samples for each individual exist
multiple_samples <- unique(metadata_fer_combined$subjectID[duplicated(metadata_fer_combined$subjectID)])
sum(!metadata_fer_combined$subjectID %in% multiple_samples)


# remove mother samples and non-fecal samples
sort_out_samples_fer <- metadata_fer_combined$run_accession[is.na(metadata_fer_combined$body_site)]
sort_out_samples_fer <- c(sort_out_samples_fer, metadata_fer_combined$run_accession[metadata_fer_combined$family_role == "mother"])
metadata_fer_combined <- metadata_fer_combined %>% 
  filter(!is.na(body_site)) %>%
  filter(family_role != "mother") %>%
  filter(!is.na(infant_age))

# test <- list.dirs("/fast/AG_Forslund/rob/studies/shotgun/ferretti_2018/fastq_files/preprocessed/", recursive = F) %>%
#   gsub("/fast/AG_Forslund/rob/studies/shotgun/ferretti_2018/fastq_files/preprocessed//", "",.) %>%
#   grep("SRR",., value = T)
# sort_out_samples_fer[sort_out_samples_fer %in% test]
# test[!test %in% metadata_fer_combined$run_accession]
# metadata_fer_combined$run_accession[!test %in% test]

# rename entries 
metadata_fer_combined[metadata_fer_combined == "NO"] <- F
metadata_fer_combined[metadata_fer_combined == "no"] <- F
metadata_fer_combined[metadata_fer_combined == "YES"] <- T
# feeding
metadata_fer_combined$food <- NA
unique(metadata_fer_combined$Exclusively.breastfed.at.1.day) # only breastfed
metadata_fer_combined$food[metadata_fer_combined$infant_age == 1] <- "breast"
unique(metadata_fer_combined$Exclusively.breastfed.at.3.days) # only breastfed
metadata_fer_combined$food[metadata_fer_combined$infant_age == 3] <- "breast"
metadata_fer_combined$food[metadata_fer_combined$infant_age == 7] <- "breast" # no data, but as most children are still breastfed at 1 month, assume here, too
metadata_fer_combined$food[metadata_fer_combined$infant_age == 30 & 
                                metadata_fer_combined$Exclusively.breastfed.at.1.month... == "TRUE"] <- "breast"
metadata_fer_combined$food[metadata_fer_combined$infant_age == 30 & 
                                metadata_fer_combined$Exclusively.breastfed.at.1.month... == "FALSE"] <- "mixed_b_f"
metadata_fer_combined$food[metadata_fer_combined$infant_age == 120 & 
                                metadata_fer_combined$Exclusively.breastfed.at.1.month... == "TRUE"] <- "breast"
metadata_fer_combined$food[metadata_fer_combined$infant_age == 120 & 
                                metadata_fer_combined$Exclusively.bottle.fed.at.4.months == "TRUE"] <- "formula"
metadata_fer_combined$food[metadata_fer_combined$infant_age == 120 & 
                                metadata_fer_combined$Mixed.breastfed.at.4.months == "TRUE"] <- "mixed_b_f"
# no solid food here at 4 months
# weight
metadata_fer_combined$weigth <- NA
metadata_fer_combined$weigth[metadata_fer_combined$infant_age == 1 & !is.na(metadata_fer_combined$infant_age)] <- 
   metadata_fer_combined$Weight.at.birth..kg.[metadata_fer_combined$infant_age == 1 & !is.na(metadata_fer_combined$infant_age)]
metadata_fer_combined$weigth[metadata_fer_combined$infant_age == 3 & !is.na(metadata_fer_combined$infant_age)] <- 
   metadata_fer_combined$Weight.at.3.days..kg.[metadata_fer_combined$infant_age == 3 & !is.na(metadata_fer_combined$infant_age)]
# convert other columns
metadata_fer_combined$mother_age <- as.numeric(metadata_fer_combined$Age)
metadata_fer_combined$Age <- NULL
colnames(metadata_fer_combined)[2] <- "sampleID"

# antibiotics one week before
metadata_fer_combined$antibiotics_one_week_before <- NA
metadata_fer_combined$antibiotics_before <- NA
logical_cols <- sapply(metadata_fer_combined, function(x) all(x %in% c("TRUE", "FALSE", NA)))
metadata_fer_combined[logical_cols] <- lapply(metadata_fer_combined[logical_cols], as.logical)
metadata_fer_combined <- metadata_fer_combined %>%
  mutate(antibiotics_one_week_before = case_when(infant_age == 1 ~ F,
                                                 (infant_age == 3) & Antibiotic.therapy.in.the.first.three.days.of.life ~ T,
                                                 (infant_age == 3) & !Antibiotic.therapy.in.the.first.three.days.of.life ~ F,
                                                 (infant_age == 7) & Antibiotic.therapy.in.the.first.three.days.of.life ~ T,
                                                 (infant_age == 7) & !Antibiotic.therapy.in.the.first.three.days.of.life ~ F,
                                                 (infant_age == 30) & Antibiotic.treatment.in.the.week.before.sampling.at.1.month ~ T,
                                                 (infant_age == 30) & !Antibiotic.treatment.in.the.week.before.sampling.at.1.month ~ F,
                                                 (infant_age == 120) & Antibiotic.treatment.in.the.week.before.sampling.at.4.months ~ T,
                                                 (infant_age == 120) & !Antibiotic.treatment.in.the.week.before.sampling.at.4.months ~ F),
         antibiotics_before = case_when(antibiotics_one_week_before ~ T,
                                        infant_age == 1 ~ F,
                                        (infant_age == 3) & !Antibiotic.therapy.in.the.first.three.days.of.life ~ F,
                                        (infant_age == 7) & !Antibiotic.therapy.in.the.first.three.days.of.life ~ F,
                                        infant_age == 30 & (Antibiotic.treatment.between.three.days.and.1.month |
                                                              Antibiotic.therapy.in.the.first.three.days.of.life) ~ T,
                                        infant_age == 30 & !Antibiotic.treatment.between.three.days.and.1.month &
                                          !Antibiotic.therapy.in.the.first.three.days.of.life ~ F,
                                        infant_age == 120 & Antibiotic.treatment.since.birth ~ T,
                                        infant_age == 120 & !Antibiotic.treatment.since.birth ~ F))



# delete columns
metadata_fer_combined[c("antibiotics_current_use", "NCBI_accession", "BMI",
                        "m", "body_site", "body_subside", "study_condition",
                        "disease", "age", "non_industrialized", "minimum_read_length",
                        "median_read_length", "family_role", "curator", "BMI..before.pregnancy.",
                        "Physical.exercise", "Exclusively.breastfed.at.1.day",
                        "Exclusively.breastfed.at.3.days", "Exclusively.breastfed.at.1.month...",
                        "Exclusively.breastfed.at.4.months", "Mixed.breastfed.at.4.months",
                        "Exclusively.bottle.fed.at.4.months", "Introduction.of.solid.food.at.4.months",
                        "age_category", "body_subsite", "Weight.at.3.days..kg."
                        # "Antibiotic.therapy.in.the.first.three.days.of.life",
                        # "Antibiotic.treatment.in.the.week.before.sampling.at.1.month",
                        # "Antibiotic.treatment.between.three.days.and.1.month",
                        # "Antibiotic.treatment.in.the.week.before.sampling.at.4.months",
                        # "feeding_practice", "Antibiotic.treatment.since.birth"
                        )] <- NULL


# unify colnames
metadata_fer_combined %<>% mutate(subject_ID = subjectID,
                                  sample_ID = sampleID,
                                  age = infant_age,
                                  family_ID = family, 
                                  birth_weight = Weight.at.birth..kg.,
                                  antibiotics_any = as.logical(Antibiotic.treatment.since.birth),
                                  geographic_location_.latitude. = 46.066666,
                                  geographic_location_.longitude. = 11.116667,
                                  .keep = "unused")



metadata_fer_combined$country <- "ITALY"
metadata_fer_combined$region <- "TRENTO"
metadata_fer_combined$study <- "ferretti_2018"
metadata_fer_combined$lifestyle <- "industrialized"
metadata_fer_combined$birthmode <-"v"
sort(colnames(metadata_fer_combined))

# move sorted out samples
sort_out_samples_fer <- unique(sort_out_samples_fer)
sort_out_files_fer <- c(paste(sort_out_samples_fer, "_1.fastq.gz", sep = ""),
                        paste(sort_out_samples_fer, "_2.fastq.gz", sep = ""),
                        paste(sort_out_samples_fer, ".fastq.gz", sep = ""))
sort_out_files_fer <- sort_out_files_fer[file.exists(sort_out_files_fer)]
file.move(sort_out_files_fer, "sorted_out/")


# move good files
final_samples_fer <- metadata_fer_combined$run_accession
final_files_fer <- c(paste(final_samples_fer, "_1.fastq.gz", sep = ""),
                     paste(final_samples_fer, "_2.fastq.gz", sep = ""),
                     paste(final_samples_fer, ".fastq.gz", sep = ""))
final_files_fer <- final_files_fer[file.exists(final_files_fer)]
file.move(final_files_fer, "fastq_files/")

# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_fer_combined$read_count, n_lines$read_counts[match(metadata_fer_combined$run_accession, n_lines$sample_ID)])  # 

# write file
write.csv(metadata_fer_combined, "metadata_fer_2018_healthy.csv")


# clean up
rm(ena_file_fer, metadata_fer, metadata_fer_2, metadata_fer_combined,
   multiple_samples, sort_out_samples_fer, sort_out_files_fer, final_samples_fer,
   final_files_fer, n_lines, sra_file, logical_cols)


# Bäckhed	2015 #################################################################
# gestational age estimated

setwd(paste(maindir, "shotgun/backhed_2015/", sep = ""))
ena_file_backhed <- read.table("filereport_read_run_PRJEB6456_tsv.txt", sep = "\t", header = T)
ena_file_backhed <- ena_file_backhed[c("study_accession", "run_accession", "read_count", "sample_alias")]
ena_file_backhed$sample_alias <- gsub("[_]", "-", ena_file_backhed$sample_alias)
ena_file_backhed$sample_alias <- gsub("[-]B", "-NB", ena_file_backhed$sample_alias)
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file[,c("Instrument","Run")]
all(ena_file_backhed$run_accession %in% sra_file$Run)
ena_file_backhed <- merge(ena_file_backhed, sra_file,by.x = "run_accession", by.y = "Run")


# read and reformat data
metadata_bac <- read.table("metadata_1.txt", sep = "\t", header = T)
metadata_bac <- metadata_bac[-grep("-3Y", metadata_bac$sample_alias),]
metadata_bac <- metadata_bac[-grep("-5Y", metadata_bac$sample_alias),]
metadata_bac[c("instrument_model", "library_name", "library_source", "library_selection",
               "library_strategy", "design_description", "library_construction_protocol", 
               "insert_size", "forward_file_name", "forward_file_md5", "reverse_file_name",
               "reverse_file_md5", "X", "X.1")] <- NULL
dict_ros_back <- read.table("dict.txt", sep = "\t", header = T)
metadata_ros_back <- read.table("metadata_2.txt", sep = "\t", header = T)
metadata_ros_back$child_id_new <- dict_ros_back$STUDY.ID.Growth.study..H2GS[match(metadata_ros_back$STUDY.ID, dict_ros_back$STUDY.ID.metagenomics)]


# merge dataframes
sum(!metadata_bac$sample_alias %in% ena_file_backhed$sample_alias)
metadata_bac$sample_alias[(!metadata_bac$sample_alias %in% ena_file_backhed$sample_alias)]
sum(!ena_file_backhed$sample_alias %in% metadata_bac$sample_alias)  # some samples lack metadata
metadata_bac_combined <- merge(metadata_bac, ena_file_backhed, by = "sample_alias", all.y = T)


# merge with roswall 2021 data
metadata_bac_combined$child_id <- as.numeric(gsub("[-].*", "", metadata_bac_combined$sample_alias))
metadata_bac_combined$timepoint <- gsub(".*[-]", "", metadata_bac_combined$sample_alias)

sum(!unique(metadata_ros_back$child_id_new) %in% metadata_bac_combined$child_id)
sum(!unique(metadata_bac_combined$child_id) %in% metadata_ros_back$child_id_new)
unique(metadata_bac_combined$child_id)[(!unique(metadata_bac_combined$child_id) %in% metadata_ros_back$child_id_new)]

nrow(metadata_bac_combined)
metadata_bac_combined <- merge(metadata_bac_combined, metadata_ros_back, by.x = "child_id", by.y = "child_id_new", all.x = T)
nrow(metadata_bac_combined)
metadata_bac_combined[metadata_bac_combined == "vaginal"] <- "v"
metadata_bac_combined[metadata_bac_combined == "sectio" | metadata_bac_combined == "section"] <- "c"
metadata_bac_combined[metadata_bac_combined == "yes"] <- "y"
metadata_bac_combined[metadata_bac_combined == "no"] <- "n"
metadata_bac_combined[metadata_bac_combined == "boy"] <- "m"
metadata_bac_combined[metadata_bac_combined == "girl"] <- "f"
metadata_bac_combined[metadata_bac_combined == "exclusively breastfeeding"] <- "breast"
metadata_bac_combined[metadata_bac_combined == "mixed feeding"] <- "mixed_b_f"
metadata_bac_combined[metadata_bac_combined == "exclusively formula feeding" | metadata_bac_combined == "Formula feeding"] <- "formula"
metadata_bac_combined[metadata_bac_combined == "no breastfeeding"] <- "solid"
metadata_bac_combined[metadata_bac_combined == "any breastfeeding"] <- "mixed_b_s"
metadata_bac_combined$sex[metadata_bac_combined$sex == "1"] <- "m"
metadata_bac_combined$sex[metadata_bac_combined$sex == "2"] <- "f"
sum(is.na(metadata_bac_combined))

# remove mother samples
sort_out_samples_bac <- metadata_bac_combined$run_accession[metadata_bac_combined$timepoint == "M"]
metadata_bac_combined <- metadata_bac_combined[metadata_bac_combined$timepoint != "M",]


# sex-column overlaps with GENDER
any(!metadata_bac_combined$sex == metadata_bac_combined$GENDER)
# combine sex columns
metadata_bac_combined$sex[is.na(metadata_bac_combined$sex)] <- metadata_bac_combined$GENDER[is.na(metadata_bac_combined$sex)]
metadata_bac_combined$GENDER <- NULL
# expand the sex information to all samples
sex_dataframe <- metadata_bac_combined[!is.na(metadata_bac_combined$sex) , c("sex", "child_id")] %>% unique
metadata_bac_combined$sex <- NULL
metadata_bac_combined <- inner_join(metadata_bac_combined, sex_dataframe, by = "child_id")

# mode_of_birth-column overlaps with Delivery.mode
any(!metadata_bac_combined$mode_of_birth == metadata_bac_combined$Delivery.mode, na.rm = T)
any(! metadata_bac_combined$Delivery.mode == metadata_bac_combined$mode_of_birth, na.rm = T)
# combine birth mode columns
metadata_bac_combined$mode_of_birth[is.na(metadata_bac_combined$mode_of_birth)] <- metadata_bac_combined$Delivery.mode[is.na(metadata_bac_combined$mode_of_birth)]
metadata_bac_combined$Delivery.mode <- NULL
# expand the sex information to all samples
birth_mode_dataframe <- metadata_bac_combined[!is.na(metadata_bac_combined$mode_of_birth) , c("mode_of_birth", "child_id")] %>% unique
metadata_bac_combined$mode_of_birth <- NULL
metadata_bac_combined <- inner_join(metadata_bac_combined, birth_mode_dataframe, by = "child_id")


# individuals with only one sample
multiple_samples <- unique(metadata_bac_combined$child_id[duplicated(metadata_bac_combined$child_id)])
sort_out_samples_bac <- c(sort_out_samples_bac, 
                          metadata_bac_combined$run_accession[!metadata_bac_combined$child_id %in% multiple_samples])
metadata_bac_combined <- metadata_bac_combined[metadata_bac_combined$child_id %in% multiple_samples,]


# fill age column
metadata_bac_combined$age <- NA
metadata_bac_combined$age[metadata_bac_combined$timepoint == "NB"] <- metadata_bac_combined$Age.at.sample.Newborn..days.[metadata_bac_combined$timepoint == "NB"]
metadata_bac_combined$age[metadata_bac_combined$timepoint == "NB" & 
                             is.na(metadata_bac_combined$age)] <- mean(metadata_bac_combined$Age.at.sample.Newborn..days., na.rm = T)
metadata_bac_combined$age[metadata_bac_combined$timepoint == "4M"] <- metadata_bac_combined$Age.at.Sampling..4.M[metadata_bac_combined$timepoint == "4M"]
metadata_bac_combined$age[metadata_bac_combined$timepoint == "4M" & 
                             is.na(metadata_bac_combined$age)] <- mean(metadata_bac_combined$Age.at.Sampling..4.M, na.rm = T)
metadata_bac_combined$age[metadata_bac_combined$timepoint == "12M"] <- metadata_bac_combined$Age.at.Sampling.12M[metadata_bac_combined$timepoint == "12M"]
metadata_bac_combined$age[metadata_bac_combined$timepoint == "12M" & 
                             is.na(metadata_bac_combined$age)] <- mean(metadata_bac_combined$Age.at.Sampling.12M, na.rm = T)

# add feeding column
metadata_bac_combined$food <- NA
metadata_bac_combined$food[metadata_bac_combined$timepoint == "NB"] <- metadata_bac_combined$feeding.practice.first.week[metadata_bac_combined$timepoint == "NB"]
metadata_bac_combined$food[metadata_bac_combined$timepoint == "4M"] <- metadata_bac_combined$feeding.practice.4M[metadata_bac_combined$timepoint == "4M"]
metadata_bac_combined$food[metadata_bac_combined$timepoint == "12M"] <- metadata_bac_combined$Any.breastfeeding.12.M[metadata_bac_combined$timepoint == "12M"]

# antibiotics
metadata_bac_combined$antibiotics_before <- F
metadata_bac_combined$antibiotics_before[metadata_bac_combined$age >= 100 & 
                                            metadata_bac_combined$Antibiotic.treatment.to.infant.0.4M...........times. > 0] <- T
metadata_bac_combined$antibiotics_before[metadata_bac_combined$age >= 340 & 
                                            metadata_bac_combined$Antibiotic.treatment.to.infant.4.12M.....times. > 0] <- T


# clean up metadata-table
metadata_bac_combined[c("timepoint", "Age.at.Sampling..4.M", "Age.at.Sampling.12M",
                        "feeding.practice.first.week", "feeding.practice.4M",
                        "Any.breastfeeding.12.M")] <- NULL

# unify colnames
colnames(metadata_bac_combined)[match(c("child_id", "sample_alias", "breastfed",
                                        "mode_of_birth"),
                                      colnames(metadata_bac_combined))] <- c("subject_ID", "sample_ID", "breast",
                                                                             "birthmode")
metadata_bac_combined$country <- "SWEDEN"
metadata_bac_combined$region <- "SWEDEN"
metadata_bac_combined$study <- "backhed_2015"
metadata_bac_combined$lifestyle <- "industrialized"
metadata_bac_combined$family_ID <- metadata_bac_combined$subject_ID
metadata_bac_combined$gestational_age <- 40
sort(colnames(metadata_bac_combined))

# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_bac_combined$read_count, n_lines$read_counts[match(metadata_bac_combined$run_accession, n_lines$sample_ID)])  # 

metadata_bac_combined$geographic_location_.latitude. <- 56.67446
metadata_bac_combined$geographic_location_.longitude. <- 12.85676

write.csv(metadata_bac_combined, "metadata_backhed_2015_healthy.csv")

# move sorted out samples
sort_out_samples_bac <- unique(sort_out_samples_bac)
sort_out_files_bac <- c(paste(sort_out_samples_bac, "_1.fastq.gz", sep = ""),
                        paste(sort_out_samples_bac, "_2.fastq.gz", sep = ""),
                        paste(sort_out_samples_bac, ".fastq.gz", sep = ""))
sort_out_files_bac <- sort_out_files_bac[file.exists(sort_out_files_bac)]
file.move(sort_out_files_bac, "sorted_out/")


# move good files
final_samples_bac <- metadata_bac_combined$run_accession
final_files_bac <- c(paste(final_samples_bac, "_1.fastq.gz", sep = ""),
                     paste(final_samples_bac, "_2.fastq.gz", sep = ""),
                     paste(final_samples_bac, ".fastq.gz", sep = ""))
final_files_bac <- final_files_bac[file.exists(final_files_bac)]
file.move(final_files_bac, "fastq_files/")

# clean up
rm(ena_file_backhed, metadata_bac, metadata_ros_back, dict_ros_back,
   metadata_bac_combined, sort_out_samples_bac, sort_out_files_bac, 
   multiple_samples, sex_dataframe, birth_mode_dataframe, final_samples_bac,
   final_files_bac, n_lines)



# Yassour	2018 #################################################################

setwd(paste(maindir, "shotgun/yassour_2018/", sep = ""))
ena_file_yas <- read.table("filereport_read_run_PRJNA475246_tsv.txt", sep = "\t", header = T)
ena_file_yas <- ena_file_yas[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file[,c("Instrument","Run")]
all(ena_file_yas$run_accession %in% sra_file$Run)
ena_file_yas <- merge(ena_file_yas, sra_file,by.x = "run_accession", by.y = "Run")


read_counts_yas <- read.table("metadata_read_counts.txt", sep = "\t", header = T)
read_counts_yas$Total.number.of.reads <- as.numeric(gsub(",", "", read_counts_yas$Total.number.of.reads))
read_counts_yas$sample_title <- paste(read_counts_yas$Family, "-", 
                                      read_counts_yas$Type, ":", read_counts_yas$Sample.timing, sep = "")

metadata_yas <- read.table("metadata_1.txt", sep = "\t", header = T)
metadata_yas[c("Number.of.child.samples", "Number.of.mother.samples", "Total.number.of.family.samples")] <- NULL

# no duplicates
sum(duplicated(ena_file_yas$sample_alias))

#combine data.frames
sum(!unique(ena_file_yas$sample_title) %in% read_counts_yas$sample_title) # some samples are in the second file
sum(!unique(read_counts_yas$sample_title) %in% ena_file_yas$sample_title)
nrow(ena_file_yas)
metadata_yas_combined <- merge(read_counts_yas, ena_file_yas, by = "sample_title")
nrow(metadata_yas_combined)

# read counts correlate well
plot(metadata_yas_combined$read_count, metadata_yas_combined$Total.number.of.reads)
metadata_yas_combined$Total.number.of.reads <- NULL

sum(!unique(metadata_yas_combined$Family) %in% metadata_yas$Family) # some samples are in the second file
sum(!unique(metadata_yas$Family) %in% metadata_yas_combined$Family)
nrow(metadata_yas_combined)
metadata_yas_combined <- merge(metadata_yas_combined, metadata_yas, by = "Family")
nrow(metadata_yas_combined)

# remove mother samples
sort_out_samples_yas <- metadata_yas_combined$run_accession[metadata_yas_combined$Type == "Mother"]
metadata_yas_combined <- metadata_yas_combined[metadata_yas_combined$Type != "Mother",]
# individuals with only one sample
multiple_samples <- unique(metadata_yas_combined$Family[duplicated(metadata_yas_combined$Family)])
sort_out_samples_yas <- c(sort_out_samples_yas, 
                          metadata_yas_combined$run_accession[!metadata_yas_combined$Family %in% multiple_samples])
metadata_yas_combined <- metadata_yas_combined[metadata_yas_combined$Family %in% multiple_samples,]

# age
metadata_yas_combined$age <- NA
metadata_yas_combined$age[metadata_yas_combined$Sample.timing == "Birth"] <- 0
metadata_yas_combined$age[metadata_yas_combined$Sample.timing == "14 days"] <- 14
metadata_yas_combined$age[metadata_yas_combined$Sample.timing == "1 month"] <- 30
metadata_yas_combined$age[metadata_yas_combined$Sample.timing == "2 months"] <- 60
metadata_yas_combined$age[metadata_yas_combined$Sample.timing == "3 months"] <- 90
metadata_yas_combined$Sample.timing <- NULL

# rename stuff
metadata_yas_combined$Delivery.mode[metadata_yas_combined$Delivery.mode == "Vaginal"] <- "v"
metadata_yas_combined$Delivery.mode[metadata_yas_combined$Delivery.mode == "C-section"] <- "c"
metadata_yas_combined$Gender[metadata_yas_combined$Gender == "Girl"] <- "f"
metadata_yas_combined$Gender[metadata_yas_combined$Gender == "Boy"] <- "m"
metadata_yas_combined$Feeding.state.at.2.weeks[metadata_yas_combined$Feeding.state.at.2.weeks == 
                                                  "Exclusive breastmilk"] <- "breast"
metadata_yas_combined$Feeding.state.at.2.weeks[metadata_yas_combined$Feeding.state.at.2.weeks == 
                                                  "Breastmilk & formula"] <- "mixed_b_f"
metadata_yas_combined$Feeding.state.at.1mo[metadata_yas_combined$Feeding.state.at.1mo == 
                                              "Exclusive breastmilk"] <- "breast"
metadata_yas_combined$Feeding.state.at.1mo[metadata_yas_combined$Feeding.state.at.1mo == 
                                              "Breastmilk & formula"] <- "mixed_b_f"
metadata_yas_combined$Feeding.state.at.2mo[metadata_yas_combined$Feeding.state.at.2mo == 
                                              "Exclusive breastmilk"] <- "breast"
metadata_yas_combined$Feeding.state.at.2mo[metadata_yas_combined$Feeding.state.at.2mo == 
                                              "Breastmilk & formula"] <- "mixed_b_f"
metadata_yas_combined$Feeding.state.at.3mo[metadata_yas_combined$Feeding.state.at.3mo == 
                                              "Exclusive breastmilk"] <- "breast"
metadata_yas_combined$Feeding.state.at.3mo[metadata_yas_combined$Feeding.state.at.3mo == 
                                              "Breastmilk & formula"] <- "mixed_b_f"
metadata_yas_combined$food <- NA
metadata_yas_combined$food[metadata_yas_combined$age == 0] <- "breast"
metadata_yas_combined$food[metadata_yas_combined$age == 14] <- metadata_yas_combined$Feeding.state.at.2.weeks[metadata_yas_combined$age == 14]
metadata_yas_combined$food[metadata_yas_combined$age == 30] <- metadata_yas_combined$Feeding.state.at.1mo[metadata_yas_combined$age == 30]
metadata_yas_combined$food[metadata_yas_combined$age == 60] <- metadata_yas_combined$Feeding.state.at.2mo[metadata_yas_combined$age == 60]
metadata_yas_combined$food[metadata_yas_combined$age == 90] <- metadata_yas_combined$Feeding.state.at.3mo[metadata_yas_combined$age == 90]

# clean cols
metadata_yas_combined[c("Type", "Formula.first.day", "Amount.of.breastmilk.at.age.3months..grams.",
                        "Feeding.state.at.2.weeks", "Formula.used.at.2.weeks",
                        "Feeding.state.at.1mo", "Formula.used.at.1mo", "Feeding.state.at.2mo",
                        "Formula.used.at.2mo", "Feeding.state.at.3mo", "Formula.used.at.3mo")] <-NULL

# unify colnames
metadata_yas_combined %<>% mutate(family_ID = Family,
                                  sample_ID = sample_title,
                                  birthmode = Delivery.mode,
                                  sex = "Gender",
                                  birth_weight =Birth.Weight/1000,
                                  gestational_age = Gestational.age,
                                  geographic_location_.latitude. = 60.192059,
                                  geographic_location_.longitude. = 24.945831,
                                  .keep = "unused")

metadata_yas_combined$subject_ID <- metadata_yas_combined$family_ID
metadata_yas_combined$country <- "FINLAND"
metadata_yas_combined$region <- "FINLAND"
metadata_yas_combined$study <- "yassour_2018"
metadata_yas_combined$lifestyle <- "industrialized"
sort(colnames(metadata_yas_combined))

# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_yas_combined$read_count, n_lines$read_counts[match(metadata_yas_combined$run_accession, n_lines$sample_ID)])  # 


write.csv(metadata_yas_combined, "metadata_yassour_2018_healthy.csv")

# move sorted out samples
sort_out_samples_yas <- unique(sort_out_samples_yas)
sort_out_files_yas <- c(paste(sort_out_samples_yas, "_1.fastq.gz", sep = ""),
                        paste(sort_out_samples_yas, "_2.fastq.gz", sep = ""),
                        paste(sort_out_samples_yas, ".fastq.gz", sep = ""))
sort_out_files_yas <- sort_out_files_yas[file.exists(sort_out_files_yas)]
file.move(sort_out_files_yas, "sorted_out/")

# move good files
final_samples_yas <- metadata_yas_combined$run_accession
final_files_yas <- c(paste(final_samples_yas, "_1.fastq.gz", sep = ""),
                     paste(final_samples_yas, "_2.fastq.gz", sep = ""),
                     paste(final_samples_yas, ".fastq.gz", sep = ""))
final_files_yas <- final_files_yas[file.exists(final_files_yas)]
file.move(final_files_yas, "fastq_files/")

# clean up
rm(ena_file_yas, metadata_yas, metadata_yas_combined, read_counts_yas,
   multiple_samples, sort_out_samples_yas, sort_out_files_yas, final_samples_yas,
   final_files_yas, n_lines, sra_file)


# Lou	2021 #####################################################################

setwd(paste(maindir, "shotgun/lou_2021/", sep = ""))
ena_file_lou <- read.table("filereport_read_run_PRJNA698986_tsv.txt", sep = "\t", header = T)
ena_file_lou <- ena_file_lou[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T, comment.char = "")
sra_file <- sra_file[,c("Instrument","Run")]
all(ena_file_lou$run_accession %in% sra_file$Run)
ena_file_lou <- merge(ena_file_lou, sra_file,by.x = "run_accession", by.y = "Run")

# remove controls and mother samples
sort_out_samples_lou <- ena_file_lou$run_accession[ena_file_lou$sample_title != "Infant "]
ena_file_lou <- ena_file_lou[ena_file_lou$sample_title == "Infant ",]
ena_file_lou <- separate(ena_file_lou, col = "sample_alias", 
                         into = c("L", "infant_ID", "age"), sep = "[_]", remove = F,
                         convert = T)
ena_file_lou <- separate(ena_file_lou, col = "age", 
                         into = c("age", "G"), sep = "G", remove = F,
                         convert = T)


read_counts_lou <- read.table("metadata_reads.txt", sep = "\t", header = T, comment.char = "")
metadata_lou_1 <- read.table("github_metadata.txt", sep = ",", header = T)
metadata_lou_1 <- separate(metadata_lou_1, col = "infant", into = c(NA, "infant_ID"),
                           sep = "[_]", remove = T, convert = T)
metadata_lou_2 <- read.table("metadata_1.txt", sep = "\t", header = T, comment.char = "")
metadata_lou_2$infant_id <- gsub("[#]", "", metadata_lou_2$Infant) %>%
   type.convert()
metadata_lou_2[c("Infant", "POD.in.NICU.if.born.prematurely", 'Hisanic', "Discharge.DOL",
                 "NEC.LOS", "LOS.DOL", "Blood.isolate", "first.NEC.DOL",
                 "recurred.NEC.DOL", "Maternal.related.information", "Maternal.Sample",
                 "Maternal.Group.B.Streptococcal..GBS..colonization.status",
                 "Term", "Gender", "Race", "Delivery", "Weight..g.", "Discharge_DOL",
                 "Disease")] <- NULL

# no duplicates
sum(duplicated(ena_file_lou$sample_alias))

# combine data.frames
sum(!unique(ena_file_lou$sample_alias) %in% read_counts_lou$public_code)
sum(!unique(read_counts_lou$public_code) %in% ena_file_lou$sample_alias)
nrow(ena_file_lou)
check_read_counts <- merge(read_counts_lou, ena_file_lou, by.x = "public_code", by.y = "sample_alias", all.y = T)
nrow(check_read_counts)
plot(check_read_counts$reads, check_read_counts$read_count)

sum(!unique(ena_file_lou$infant_ID) %in% metadata_lou_1$infant_ID)  # some samples lack metadata
sum(!unique(metadata_lou_1$infant_ID) %in% ena_file_lou$infant_ID)  # all metadata has samples
nrow(ena_file_lou)
metadata_lou_combined <- merge(metadata_lou_1, ena_file_lou, by = "infant_ID", all.y = T)
nrow(metadata_lou_combined)
# remove samples without metadata (contaminated according to paper)
sort_out_samples_lou <- c(sort_out_samples_lou, 
                          metadata_lou_combined$run_accession[is.na(metadata_lou_combined$num_persisters)])
metadata_lou_combined <- metadata_lou_combined[!is.na(metadata_lou_combined$num_persisters),]

# take feeding information from this table instead of the other one, sounds more realistic
metadata_lou_combined[c("Feeds", "Prolacta", "Cessation_BRM", "Solid_Food_DOL" )] <- NULL
sum(!unique(metadata_lou_combined$infant_ID) %in% metadata_lou_2$infant_id)  # some samples lack metadata
sum(!unique(metadata_lou_2$infant_id) %in% metadata_lou_combined$infant_ID)  # all metadata has samples
nrow(metadata_lou_combined)
metadata_lou_combined <- merge(metadata_lou_combined, metadata_lou_2, by.x = "infant_ID", by.y = "infant_id", all.x = T)
nrow(metadata_lou_combined)

# each child has multiple samples
multiple_samples <- unique(metadata_lou_combined$infant_ID[duplicated(metadata_lou_combined$infant_ID)])
metadata_lou_combined$run_accession[!metadata_lou_combined$infant_ID %in% multiple_samples]

# remove preterms
sort_out_samples_lou <- c(sort_out_samples_lou, 
                          metadata_lou_combined$run_accession[metadata_lou_combined$Term == "Preterm"])
metadata_lou_combined <- metadata_lou_combined[metadata_lou_combined$Term != "Preterm" ,]

# rename stuff
metadata_lou_combined$Gender[metadata_lou_combined$Gender == "Female"] <- "f"
metadata_lou_combined$Gender[metadata_lou_combined$Gender == "Female"] <- "f"
metadata_lou_combined$Delivery[metadata_lou_combined$Delivery == "C-section"] <- "c"
metadata_lou_combined$Delivery[metadata_lou_combined$Delivery == "Vaginal"] <- "v"
metadata_lou_combined[metadata_lou_combined == "No"] <- "n"

# feeding
# metadata_lou_combined$Feeds[metadata_lou_combined$Feeds == "BRM"] <- "breast"
# metadata_lou_combined$Feeds[metadata_lou_combined$Feeds == "Mix"] <- "mixed_b_f"
metadata_lou_combined$food <- NA
metadata_lou_combined$food[metadata_lou_combined$age <= metadata_lou_combined$Breastfeeding.Cessation.DOL & 
                         metadata_lou_combined$age <= metadata_lou_combined$Solid.Food.DOL &
                         metadata_lou_combined$Feeds == "BRM"] <- "breast"
metadata_lou_combined$food[metadata_lou_combined$age <= metadata_lou_combined$Breastfeeding.Cessation.DOL & 
                              metadata_lou_combined$age <= metadata_lou_combined$Solid.Food.DOL &
                              metadata_lou_combined$Feeds == "Mix"] <- "mixed_b_f"
metadata_lou_combined$food[metadata_lou_combined$age <= metadata_lou_combined$Breastfeeding.Cessation.DOL & 
                              metadata_lou_combined$age > metadata_lou_combined$Solid.Food.DOL &
                              metadata_lou_combined$Feeds == "BRM"] <- "mixed_b_s"
metadata_lou_combined$food[metadata_lou_combined$age > metadata_lou_combined$Breastfeeding.Cessation.DOL & 
                              metadata_lou_combined$age <= metadata_lou_combined$Solid.Food.DOL] <- "formula"
metadata_lou_combined$food[metadata_lou_combined$age > metadata_lou_combined$Breastfeeding.Cessation.DOL & 
                              metadata_lou_combined$age > metadata_lou_combined$Solid.Food.DOL] <- "solid"
metadata_lou_combined$food[is.na(metadata_lou_combined$Breastfeeding.Cessation.DOL) &
                              metadata_lou_combined$age > metadata_lou_combined$Solid.Food.DOL] <- "solid"
metadata_lou_combined$Feeds <- NULL

metadata_lou_combined[metadata_lou_combined == ""] <- NA
metadata_lou_combined <- metadata_lou_combined[,!unlist(lapply(metadata_lou_combined, function(x) all(is.na(x))))]

# Antibiotics one week/some time before sample or at all
metadata_lou_combined$antibiotics_one_week_before <- F
for (abx in c("AbxA.DOL", "AbxB.DOL", "AbxC.DOL", "AbxD.DOL", "AbxE.DOL", "AbxF.DOL", "AbxG.DOL")) {
   abx_start <- paste(abx, "_s", sep = "")
   abx_end <- paste(abx, "_e", sep = "")
   metadata_lou_combined <- separate(metadata_lou_combined, into = c(abx_start, abx_end),
                                     col = abx, sep = "[-]", remove = T, convert = T)
   metadata_lou_combined$antibiotics_one_week_before[metadata_lou_combined$age >= metadata_lou_combined[c(abx_start)] &
                                                metadata_lou_combined$age -7 <= metadata_lou_combined[c(abx_end)]] <- T
}

metadata_lou_combined$antibiotics_before <- F
for (abx in c("AbxA.DOL", "AbxB.DOL", "AbxC.DOL", "AbxD.DOL", "AbxE.DOL", "AbxF.DOL", "AbxG.DOL")) {
   abx_start <- paste(abx, "_s", sep = "")
   metadata_lou_combined$antibiotics_before[metadata_lou_combined$age >= metadata_lou_combined[c(abx_start)]] <- T
}

metadata_lou_combined$antibiotics_any <- F
for (abx in c("AbxA.DOL", "AbxB.DOL", "AbxC.DOL", "AbxD.DOL", "AbxE.DOL", "AbxF.DOL", "AbxG.DOL")) {
   abx_start <- paste(abx, "_s", sep = "")
   metadata_lou_combined$antibiotics_any[!is.na(metadata_lou_combined[c(abx_start)])] <- T
}

# delete columns
metadata_lou_combined[c("Term", "num_persisters",  "perc_persisters", "num_early_colonizers",
                        "Disease", "Discharge_DOL", "Initial.Abx.Admission",
                        "Maternal.Chorioamnionitis", "G", "Prolacta", "sample_title",
                        "Maternal.Anepartum.Antibiotics", "Maternal.Antibiotics.Usage")] <- NULL

# unify colnames
colnames(metadata_lou_combined)[match(c("infant_ID", "Gender", "Delivery", "Weight",
                                        "sample_alias", "Gestational.Age"),
                                      colnames(metadata_lou_combined))] <- c("subject_ID", "sex", "birthmode", "weight",
                                                                             "sample_ID", "gestational_age")

metadata_lou_combined$geographic_location_.latitude. <- 40.440624
metadata_lou_combined$geographic_location_.longitude. <- -79.995888
metadata_lou_combined$family_ID <- metadata_lou_combined$subject_ID
metadata_lou_combined$country <- "USA"
metadata_lou_combined$region <- "PENNSYLVANIA"
metadata_lou_combined$study <- "lou_2021"
metadata_lou_combined$lifestyle <- "industrialized"
sort(colnames(metadata_lou_combined))

# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_lou_combined$read_count, n_lines$read_counts[match(metadata_lou_combined$run_accession, n_lines$sample_ID)])  # 


write.csv(metadata_lou_combined, "metadata_lou_2021_healthy.csv")

# move sorted out samples
sort_out_samples_lou <- unique(sort_out_samples_lou)
sort_out_files_lou <- c(paste(sort_out_samples_lou, "_1.fastq.gz", sep = ""),
                        paste(sort_out_samples_lou, "_2.fastq.gz", sep = ""),
                        paste(sort_out_samples_lou, ".fastq.gz", sep = ""))
sort_out_files_lou <- sort_out_files_lou[file.exists(sort_out_files_lou)]
file.move(sort_out_files_lou, "sorted_out/")

# move good files
final_samples_lou <- metadata_lou_combined$run_accession
final_files_lou <- c(paste(final_samples_lou, "_1.fastq.gz", sep = ""),
                     paste(final_samples_lou, "_2.fastq.gz", sep = ""),
                     paste(final_samples_lou, ".fastq.gz", sep = ""))
final_files_lou <- final_files_lou[file.exists(final_files_lou)]
file.move(final_files_lou, "fastq_files/")


# clean up
rm(check_read_counts, ena_file_lou, metadata_lou_1, metadata_lou_2,
   metadata_lou_combined, read_counts_lou, abx, abx_start, abx_end,
   multiple_samples, sort_out_samples_lou, sort_out_files_lou, final_samples_lou,
   final_files_lou, n_lines, sra_file)


# Vatanen	2018 #################################################################

setwd(paste(maindir, "shotgun/vatanen_2018/", sep = ""))
ena_file_vatanen <- read.table("filereport_read_run_PRJNA497734_tsv_WGS.txt", sep = "\t", header = T)
ena_file_vatanen <- ena_file_vatanen[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
ena_file_vatanen$sample_title <- NULL
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file[,c("Instrument","Run")]
all(ena_file_vatanen$run_accession %in% sra_file$Run)
ena_file_vatanen <- merge(ena_file_vatanen, sra_file,by.x = "run_accession", by.y = "Run")

# read and reformat data
birth_vat <- read.table("birth.txt", sep = "\t", header = T)
colnames(birth_vat)[5] <- "birth_mode"
birth_vat$birth_mode[birth_vat$birth_mode == 0] <- "v"
birth_vat$birth_mode[birth_vat$birth_mode == 1] <- "c"
birth_vat[birth_vat == "Male"] <- "m"
birth_vat[birth_vat == "Female"] <- "f"
diabetes_vat <- read.table("diabetes.txt", sep = "\t", header = T)
feeding_vat <- read.table("feeding.txt", sep = "\t", header = T)
feeding_vat$milk_first_three_days[feeding_vat$milk_first_three_days == "mothers_breast_milk"] <- "breast"
feeding_vat$milk_first_three_days[feeding_vat$milk_first_three_days == "multiple_types_or_not_reported"] <- "unknown"
feeding_2_vat <- read.table("feeding_2.txt", sep = "\t", header = T)
feeding_2_vat$start_month <- feeding_2_vat$start_month * 30 
feeding_2_vat <- reshape(feeding_2_vat, timevar = "dietary_compound", idvar = "subjectID", direction = "wide")
feeding_2_vat[c("milk_first_three_days", "start_month.apple", "start_month.banana",
                "start_month.barley", "start_month.beef", "start_month.berries",
                "start_month.cabbs", "start_month.carrot", "start_month.corn",
                "start_month.cowsmilk", "start_month.egg", "start_month.fish",
                "start_month.icecream", "start_month.milkprod", "start_month.oat",
                "start_month.pear", "start_month.peas", "start_month.plum",
                "start_month.pork", "start_month.potato", "start_month.poultry",
                "start_month.rice", "start_month.rye", "start_month.spinach",
                "start_month.tomato", "start_month.wheat", "start_month.mango",
                "start_month.beetroot", "start_month.sweetpota", "start_month.soymilk",
                "start_month.turnip")] <- NULL
growth_vat <- read.table("growth.txt", sep = "\t", header = T)
samples_vat <- read.table("samples.txt", sep = "\t", header = T)
samples_vat <- samples_vat[!samples_vat$gid_wgs == "",]   # only samples with 16S-data
samples_vat$gid_16s <- NULL

# no duplicated samples
sum(duplicated(ena_file_vatanen$sample_alias))

#combine data.frames
sum(!ena_file_vatanen$library_name %in% samples_vat$gid_wgs)
sum(!samples_vat$gid_wgs %in% ena_file_vatanen$library_name)
nrow(ena_file_vatanen)
metadata_vat_combined <- merge(samples_vat, ena_file_vatanen, by.x = "gid_wgs", by.y = "library_name")
nrow(metadata_vat_combined)

sum(!metadata_vat_combined$subjectID %in% birth_vat$subjectID)  # all samples have birth data
nrow(metadata_vat_combined)
metadata_vat_combined <- merge(metadata_vat_combined, birth_vat, by = "subjectID")
nrow(metadata_vat_combined)

metadata_vat_combined[which(!metadata_vat_combined$subjectID %in% diabetes_vat$subjectID),]  # information missing for one individual only
nrow(metadata_vat_combined)
metadata_vat_combined <- merge(metadata_vat_combined, diabetes_vat, all.x = T,by = "subjectID")
nrow(metadata_vat_combined)

metadata_vat_combined[(!metadata_vat_combined$subjectID %in% feeding_vat$subjectID),] # some samples lack feeding metadata
metadata_vat_combined[(!metadata_vat_combined$subjectID %in% feeding_2_vat$subjectID),]
nrow(metadata_vat_combined)
metadata_vat_combined <- merge(metadata_vat_combined, feeding_vat, all.x = T, by = "subjectID")
metadata_vat_combined <- merge(metadata_vat_combined, feeding_2_vat, all.x = T, by = "subjectID")
nrow(metadata_vat_combined)

sum(!unique(metadata_vat_combined$subjectID) %in% growth_vat$subjectID)  # some samples lack growth data
sum(!unique(growth_vat$subjectID) %in% metadata_vat_combined$subjectID)  
nrow(metadata_vat_combined)
metadata_vat_combined <- merge(metadata_vat_combined, growth_vat, all.x = T,by = "subjectID")
nrow(metadata_vat_combined)

# food column
metadata_vat_combined$food <- "unknown"
metadata_vat_combined$food[metadata_vat_combined$age_at_collection >=
                                 metadata_vat_combined$start_month.solid_start_approx] <- "solid"
metadata_vat_combined$food[metadata_vat_combined$age_at_collection <= 
                                 metadata_vat_combined$bf_length] <- "mixed"
metadata_vat_combined$food[metadata_vat_combined$age_at_collection >= 
                                 metadata_vat_combined$start_month.solid_start_approx &
                                 metadata_vat_combined$age_at_collection <= 
                                 metadata_vat_combined$bf_length] <- "mixed_b_s"
metadata_vat_combined$food[metadata_vat_combined$age_at_collection <= 
                                 metadata_vat_combined$start_month.solid_start_approx &
                                 metadata_vat_combined$age_at_collection > 
                                 metadata_vat_combined$bf_length_exclusive] <- "mixed_b_f"
metadata_vat_combined$food[metadata_vat_combined$age_at_collection <= 
                                 metadata_vat_combined$start_month.solid_start_approx &
                                 metadata_vat_combined$age_at_collection > 
                                 metadata_vat_combined$bf_length] <- "formula"
metadata_vat_combined$food[metadata_vat_combined$age_at_collection <= 
                                 metadata_vat_combined$bf_length_exclusive] <- "breast"


# remove some columns
metadata_vat_combined[c("start_month.apple", "start_month.banana", "start_month.barley",
                        "start_month.beef", "start_month.berries",  "start_month.cabbs",
                        "start_month.carrot", "start_month.corn", "start_month.cowsmilk",
                        "start_month.egg", "start_month.fish",
                        "start_month.icecream", "start_month.milkprod", "start_month.oat",                    
                        "start_month.pear", "start_month.peas", "start_month.plum",                     
                        "start_month.pork", "start_month.potato", "start_month.poultry",                  
                        "start_month.rice", "start_month.rye", "start_month.spinach",                  
                        "start_month.tomato", "start_month.wheat", "start_month.mango",                    
                        "start_month.beetroot", "start_month.sweetpota", "start_month.soymilk",                  
                        "start_month.turnip")] <- NULL

# add info on antibiotic usage:
abx_data <- read.table("abx_data_1.txt", sep = "\t", header = T) %>%
  mutate(age = as.numeric(Age..months.),
         .keep = "unused")
abx_data <- abx_data %>%
  group_by(Subject) %>%
  summarise(first_antibiotics = sum(age, na.rm = T)) %>%
  mutate(first_antibiotics = as.numeric(first_antibiotics) * 30)
metadata_vat_combined <- left_join(metadata_vat_combined, abx_data, by = join_by(subjectID == Subject))

abx_data_2 <- readxl::read_xlsx("1-s2.0-S1931312815000219-mmc2.xlsx")
abx_data_2[,2:86] <- NULL
metadata_vat_combined <- left_join(metadata_vat_combined, abx_data_2,
                                   by = join_by(subjectID == "ID, E=Espoo, Finland, T=Tartu, Estonia"))
abx_data_3 <- readxl::read_xlsx("1-s2.0-S0092867416303981-mmc2.xlsx", sheet = "Basics and medications") %>%
  tibble::column_to_rownames(var = "Participant")
abx_data_3 <- abx_data_3[,c(grep("_Age", colnames(abx_data_3)))]
abx_data_3 <- data.frame(first_abx_2 = apply(abx_data_3, 1, function(x) min(x, na.rm = T))) %>%
  mutate(first_abx_2 = as.numeric(first_abx_2) * 30)

abx_data_3$subject <- rownames(abx_data_3)
metadata_vat_combined <- left_join(metadata_vat_combined, abx_data_3,
                                   by = join_by(subjectID == subject))

metadata_vat_combined <- metadata_vat_combined %>%
  mutate(antibiotics_before = case_when(first_antibiotics < age_at_collection ~ T,
                                        !is.na("First antibiotic drug between 3 and 6 months") & age_at_collection > 90 ~ T,
                                        !is.na("First antibiotic drug between 6 and 12 months") & age_at_collection > 180 ~ T,
                                        !is.na("First antibiotic drug between 12 and 18 months") & age_at_collection > 365 ~ T,
                                        !is.na("First antibiotic drug between 18 and 24 months") & age_at_collection > 540 ~ T,
                                        !is.na("First antibiotic drug between 24 and 36 months") & age_at_collection > 730 ~ T,
                                        first_abx_2 < age_at_collection ~ T,
                                        .default = F))
metadata_vat_combined[,c("First antibiotic drug between 3 and 6 months", "Second antibiotic drug between 3 and 6 months",  
                         "Third antibiotic drug between 3 and 6 months", "First antibiotic drug between 6 and 12 months", 
                         "Second antibiotic drug between 6 and 12 months", "Third antibiotic drug between 6 and 12 months",  
                         "Fourth antibiotic drug between 6 and 12 months", "Fifth antibiotic drug between 6 and 12 months",  
                         "First antibiotic drug between 12 and 18 months",  "Second antibiotic drug between 12 and 18 months",
                         "Third antibiotic drug between 12 and 18 months",  "Fourth antibiotic drug between 12 and 18 months",
                         "First antibiotic drug between 18 and 24 months",  "Second antibiotic drug between 18 and 24 months",
                         "Third antibiotic drug between 18 and 24 months",  "First antibiotic drug between 24 and 36 months",
                         "Second antibiotic drug between 24 and 36 months", "Third antibiotic drug between 24 and 36 months",
                         "Fourth antibiotic drug between 24 and 36 months")] <- NULL
         


# remove individuals with one sample only
# no, keep them, for comparison with 16S-data!
# multiple_samples <- unique(metadata_vat_combined$subjectID[duplicated(metadata_vat_combined$subjectID)])
# sort_out_samples_vat <- metadata_vat_combined$run_accession[!metadata_vat_combined$subjectID %in% multiple_samples]
# metadata_vat_combined <- metadata_vat_combined[metadata_vat_combined$subjectID %in% multiple_samples,]

# unify colnames
colnames(metadata_vat_combined)[match(c("subjectID", "sampleID", "age_at_collection",
                                        "birth_mode", "gender", "bf_length_exclusive"),
                                      colnames(metadata_vat_combined))] <- c("subject_ID", "sample_ID", "age",
                                                                             "birthmode", "sex", "exclusive_breastfeeding")

metadata_vat_combined$country[metadata_vat_combined$country == "FIN"] <- "FINNLAND"
metadata_vat_combined$country[metadata_vat_combined$country == "RUS"] <- "RUSSIA"
metadata_vat_combined$country[metadata_vat_combined$country == "EST"] <- "ESTONIA"
metadata_vat_combined$region <- metadata_vat_combined$country
metadata_vat_combined$region[metadata_vat_combined$country == "RUSSIA"] <- "KARELIA"
metadata_vat_combined$study <- "vatanen_2018"
metadata_vat_combined$lifestyle <- "industrialized"
metadata_vat_combined$family_ID <- metadata_vat_combined$subject_ID
sort(colnames(metadata_vat_combined))

metadata_vat_combined %<>% mutate(geographic_location_.latitude. = case_when(country == "FINNLAND" ~ 60.192059,
                                                                             country == "RUSSIA" ~ 61.78491,
                                                                             country == "ESTONIA" ~ 58.378025),
                                  geographic_location_.longitude. = case_when(country == "FINNLAND" ~ 24.945831,
                                                                              country == "RUSSIA" ~ 34.34691,
                                                                              country == "ESTONIA" ~ 26.728493))


# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_vat_combined$read_count, n_lines$read_counts[match(metadata_vat_combined$run_accession, n_lines$sample_ID)])  # 


write.csv(metadata_vat_combined, "metadata_vat_2018_healthy.csv")

# sort_out_samples_vat <- unique(sort_out_samples_vat)
# sort_out_files_vat <- c(paste(sort_out_samples_vat, "_1.fastq.gz", sep = ""),
#                         paste(sort_out_samples_vat, "_2.fastq.gz", sep = ""),
#                         paste(sort_out_samples_vat, ".fastq.gz", sep = ""))
# sort_out_files_vat <- sort_out_files_vat[file.exists(sort_out_files_vat)]
# file.move(sort_out_files_vat, "sorted_out/")

# move good files
final_samples_vat <- metadata_vat_combined$run_accession
final_files_vat <- c(paste(final_samples_vat, "_1.fastq.gz", sep = ""),
                     paste(final_samples_vat, "_2.fastq.gz", sep = ""),
                     paste(final_samples_vat, ".fastq.gz", sep = ""))
final_files_vat <- final_files_vat[file.exists(final_files_vat)]
file.move(final_files_vat, "fastq_files/")


rm(ena_file_vatanen, birth_vat, diabetes_vat, feeding_2_vat, feeding_vat, growth_vat,
   samples_vat, metadata_vat_combined, multiple_samples, sort_out_samples_vat,
   sort_out_files_vat, final_samples_vat, final_files_vat, n_lines, abx_data, abx_data_2,
   abx_data_3, sra_file)




# # Chu	2017 ###################################################################
# # maybe also integrate _pools information? need to know which timepoint is which.
# # metadata not good, better leave it out
# 
# setwd(paste(maindir, "shotgun/chu_2017/", sep = ""))
# ena_file_chu <- read.table("filereport_read_run_PRJNA322188_tsv_stool_only_wgs.txt", sep = "\t", header = T)
# ena_file_chu <- ena_file_chu[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
# ena_file_chu$child_ID <- gsub(".Stool.*", "", ena_file_chu$sample_title)
# ena_file_chu$sample_ID <- ena_file_chu$sample_title
# # remove unclear time points
# sort_out_samples_chu <- ena_file_chu$run_accession[grep("STOOL3", ena_file_chu$sample_title)]
# ena_file_chu <- ena_file_chu[-grep("STOOL3", ena_file_chu$sample_title),]
# sort_out_samples_chu <- c(sort_out_samples_chu, ena_file_chu$run_accession[grep("STOOL4", ena_file_chu$sample_title)])
# ena_file_chu <- ena_file_chu[-grep("STOOL4", ena_file_chu$sample_title),]
# # remove mother samples and all individuals with one sample only
# sort_out_samples_chu <- c(sort_out_samples_chu, ena_file_chu$run_accession[grep("maternal", ena_file_chu$sample_title)])
# ena_file_chu <- ena_file_chu[-grep("maternal", ena_file_chu$sample_title),]
# sort_out_samples_chu <- c(sort_out_samples_chu, ena_file_chu$run_accession[grep("[.]M[.]STOOL", ena_file_chu$sample_title)])
# ena_file_chu <- ena_file_chu[-grep("[.]M[.]STOOL", ena_file_chu$sample_title),]
# # only keep samples with at least 2 time points
# multiple_samples <- unique(ena_file_chu$child_ID[duplicated(ena_file_chu$child_ID)])
# sort_out_samples_chu <- c(sort_out_samples_chu, ena_file_chu$run_accession[!ena_file_chu$child_ID %in% multiple_samples])
# ena_file_chu <- ena_file_chu[ena_file_chu$child_ID %in% multiple_samples,]
# 
# 
# 
# metadata_chu_child <- read.table("infant_metadata.txt", sep = "\t", header = T)
# metadata_chu_child$Mode.of.Delivery[metadata_chu_child$Mode.of.Delivery == "Cesarean"] <- "c"
# metadata_chu_child$Mode.of.Delivery[metadata_chu_child$Mode.of.Delivery == "Vaginal"] <- "v"
# metadata_chu_mother <- read.table("mother_metadata.txt", sep = "\t", header = T)
# metadata_chu_mother$Mode.of.Delivery <- NULL
# metadata_chu_mother$Term.v..Preterm <- NULL
# metadata_chu_mother$Study <- NULL
# metadata_chu_stats <- read.table("mapping_stats.txt", sep = "\t", header = T)
# metadata_chu_stats <- metadata_chu_stats[metadata_chu_stats$Site == "Stool",]
# metadata_chu_stats <- metadata_chu_stats[grep("maternal", metadata_chu_stats$ID ),]
# metadata_chu_stats$Subject.ID <- gsub("[.]", "", metadata_chu_stats$Subject.ID)
# metadata_chu_stats$Visit[metadata_chu_stats$Visit == "Delivery"] <- "STOOL1"
# metadata_chu_stats$Visit[metadata_chu_stats$Visit == "Postpartum"] <- "STOOL2"
# metadata_chu_stats$sample_ID <- paste(metadata_chu_stats$Subject.ID, metadata_chu_stats$Visit, sep = "_")
# metadata_chu_stats$ID <- NULL
# metadata_chu_pools <- read.table("sequencing_pools.txt", sep = "\t", header = T, comment.char = "")
# metadata_chu_pools <- metadata_chu_pools[metadata_chu_pools$Body.Site == "Stool",]
# metadata_chu_pools <- metadata_chu_pools[metadata_chu_pools$Mother == "Infant",]
# metadata_chu_pools <- metadata_chu_pools[metadata_chu_pools$Cohort == "Longitudinal",]
# 
# sum(duplicated(ena_file_chu$sample_alias)) # no duplicates!
# 
# #combine data.frames
# sum(!unique(ena_file_chu$child_ID) %in% metadata_chu_child$ID)
# sum(!unique(metadata_chu_child$ID) %in% ena_file_chu$child_ID)
# nrow(ena_file_chu)
# metadata_chu_combined <- merge(metadata_chu_child, ena_file_chu, by.x = "ID", by.y = "child_ID")
# nrow(metadata_chu_combined)
# 
# sum(!unique(ena_file_chu$child_ID) %in% metadata_chu_mother$ID)
# nrow(ena_file_chu)
# metadata_chu_combined <- merge(metadata_chu_combined, metadata_chu_mother, by = "ID")
# nrow(metadata_chu_combined)
# 
# sum(!unique(ena_file_chu$child_ID) %in% metadata_chu_mother$ID)
# nrow(ena_file_chu)
# metadata_chu_combined <- merge(metadata_chu_combined, metadata_chu_stats, by = "sample_ID", all.x = T)
# nrow(metadata_chu_combined)
# 
# metadata_chu_combined[metadata_chu_combined == "Yes"] <- "y"
# metadata_chu_combined[metadata_chu_combined == "No"] <- "n"
# metadata_chu_combined$Gender[metadata_chu_combined$Gender == "Female"] <- "f"
# metadata_chu_combined$Gender[metadata_chu_combined$Gender == "Male"] <- "m"
# metadata_chu_combined$antibiotics_any[metadata_chu_combined$Postnatal.Antibiotics == "None"] <- "n"
# metadata_chu_combined$antibiotics_any[metadata_chu_combined$Postnatal.Antibiotics != "None"] <- "y"
# metadata_chu_combined$Postnatal.Antibiotics <- NULL
# # missing timepoint information is from time point 2
# metadata_chu_combined$time_point[is.na(metadata_chu_combined$time_point)] <- 36
# 
# 
# # filter out Preterm samples:
# sort_out_samples_chu <- c(sort_out_samples_chu, metadata_chu_combined$run_accession[metadata_chu_combined$Term.v..Preterm == "Preterm"])
# metadata_chu_combined <- metadata_chu_combined[!metadata_chu_combined$Term.v..Preterm == "Preterm",]
# 
# # move sorted out samples
# sort_out_samples_chu <- unique(sort_out_samples_chu)
# sort_out_files_chu <- c(paste(sort_out_samples_chu, "_1.fastq.gz", sep = ""),
#                         paste(sort_out_samples_chu, "_2.fastq.gz", sep = ""),
#                         paste(sort_out_samples_chu, ".fastq.gz", sep = ""))
# sort_out_files_chu <- sort_out_files_chu[file.exists(sort_out_files_chu)]
# file.move(sort_out_files_chu, "sorted_out/")
# 
# 
# # clean up cols
# metadata_chu_combined[c("Study", "Term.v..Preterm", "Gestational.Age...Delivery",
#                         "Multiples", "Subject.ID", "Site", "Visit", "Total.Paired.Reads",
#                         "Paired.Reads.after.Human.Filtering", "X...", "HUMAnN.Mapped",
#                         "X....1")] <- NULL
# 
# # unify colnames
# colnames(metadata_chu_combined)[match(c("ID", "Gender", "Gestational.Age..Weeks.Days.",
#                                         "Birth.Weight..g.", "Mode.of.Delivery"),
#                                       colnames(metadata_chu_combined))] <- c("subject_ID", "sex", "gestational_age",
#                                                                              "birth_weight", "birthmode")
# 
# metadata_chu_combined$country <- "USA"
# metadata_chu_combined$region <- "TEXAS"
# metadata_chu_combined$study <- "chu_2017"
# 
# 
# 
# # write file
# write.csv(metadata_chu_combined, "metadata_chu_2017_healthy.csv")
# 
# rm(ena_file_chu, multiple_samples, sort_out_samples_chu, metadata_chu_child,
#    metadata_chu_mother, metadata_chu_pools, metadata_chu_stats, metadata_chu_combined, sort_out_files_chu)



# Reyman 2019 (waiting for metadata) ###########################################

setwd(paste(maindir, "shotgun/reyman_2019/", sep = ""))
ena_file_reyman <- read.table("filereport_read_run_PRJNA555020_tsv.txt", sep = "\t", header = T)
ena_file_reyman <- ena_file_reyman[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]

metadata_rey <- read.table("supplement_1.txt", sep = '\t', header = T)
metadata_rey_sra <- read.table("SraRunTable.txt", sep = ',', header = T) %>%
  mutate(run_accession = Run,
         subject_ID = PID,
         sample_ID = as.character(Sample.Name),
         age = Time_point,
         .keep = "unused")



metadata_rey_sra[,c("Assay.Type", "BioSampleModel", "Bytes", "Center.Name", "Collection_Date",
                    "Consent", "DATASTORE.filetype", "DATASTORE.provider", "DATASTORE.region",
                    "env_biome", "geo_loc_name_country", "geo_loc_name_country_continent",
                    "geo_loc_name", "HOST", "Lat_Lon", "LibraryLayout",
                    "LibrarySelection", "LibrarySource", "Organism", "Platform", "sample_names")] <- NULL

sra_file_bog <- read.table("SraRunTable_bogaert_2023.txt", sep = ",", header = T, comment.char = "", colClasses = "character")
sra_file_bog[,c("Assay.Type", "AvgSpotLen", "BioSampleModel", "Bytes", "Center.Name",
            "Collection_Date", "Consent", "DATASTORE.filetype", "DATASTORE.provider",
            "DATASTORE.region", "geo_loc_name_country", "geo_loc_name_country_continent",
            "geo_loc_name", "HOST", "Instrument", "Library.Name", "LibraryLayout",
            "LibrarySelection", "Organism", "Platform", "ReleaseDate", "create_date",
            "version", "sample_id..run.", "SRA.Study", "LibrarySource",
            "niche..run.", "niche", "Bases", "BioProject", "BioSample", "Experiment",
            "miseq_run_nr..run.", "miseq_run_nr", "Sample.Name", "run_accession")] <- NULL
sra_file_bog <- sra_file_bog %>%
  filter(Isolation_Source == "fae",
         !is.na(TIME)) %>%
  mutate(run_accession = Run,
         birthmode = ifelse(Birth_mode == "vag", yes = "v", no = "c"),
         sex = ifelse(gender == "female", yes = "f", no = "m"),
         subject_ID = gsub("M", "", Subject) %>% as.numeric(),
         sample_ID = Sample_id,
         family_ID = Subject,
         age = 7,
         .keep = "unused") %>% filter(!duplicated(subject_ID))
sra_file_bog[,c("Isolation_Source", "TIME", "sample_ID", "run_accession")] <- NULL


# no duplicates
sum(duplicated(metadata_rey_sra$Sample.Name))

metadata_rey_sra <- left_join(metadata_rey_sra, sra_file_bog, by = "subject_ID")
metadata_rey_sra$age <- 42
metadata_rey_sra$age.x <- NULL
metadata_rey_sra$age.y <- NULL


metadata_rey_sra$country <- "NETHERLANDS"
metadata_rey_sra$region <- "UTRECHT"
metadata_rey_sra$study <- "reyman_2019"
metadata_rey_sra$lifestyle <- "industrialized"
metadata_rey_sra$family_ID <- metadata_rey_sra$subject_ID
sort(colnames(metadata_rey_sra))
metadata_rey_sra$geographic_location_.latitude. <- 52.306100
metadata_rey_sra$geographic_location_.longitude. <- 4.690700



# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_rey_sra$Bases, n_lines$read_counts[match(metadata_rey_sra$run_accession, n_lines$sample_ID)])  # 


write.csv(metadata_rey_sra, "metadata_rey_2019_healthy.csv")

# do not sort out samples with only one time point so I can use them still for functional profiling comparison with picrust from 16S
# # move sorted out samples
# sort_out_samples_reyman <- unique(sort_out_samples_reyman)
# sort_out_files_reyman <- c(paste(sort_out_samples_reyman, "_1.fastq.gz", sep = ""),
#                            paste(sort_out_samples_reyman, "_2.fastq.gz", sep = ""),
#                            paste(sort_out_samples_reyman, ".fastq.gz", sep = ""))
# sort_out_files_reyman <- sort_out_files_reyman[file.exists(sort_out_files_reyman)]
# file.move(sort_out_files_reyman, "sorted_out/")

# move good files
final_samples_reyman <- metadata_rey_sra$run_accession
final_files_reyman <- c(paste(final_samples_reyman, "_1.fastq.gz", sep = ""),
                        paste(final_samples_reyman, "_2.fastq.gz", sep = ""),
                        paste(final_samples_reyman, ".fastq.gz", sep = ""))
final_files_reyman <- final_files_reyman[file.exists(final_files_reyman)]
file.move(final_files_reyman, "fastq_files/")


# clean up
rm(ena_file_reyman, metadata_rey, metadata_rey_sra, final_files_reyman, final_samples_reyman,
   n_lines, sra_file_bog)




# Wampach	2018 #################################################################

setwd(paste(maindir, "shotgun/wampach_2018/", sep = ""))
ena_file_wampach <- read.table("COSMIC_metadata_NCBI_WGS.txt", sep = ",", header = T)
ena_file_wampach <- ena_file_wampach[c("Run", "AvgSpotLen", "Bases", "BioSample",
                                       "Collection_date", "env_biome",
                                       "env_material", "Library.Name", "Sample.Name",
                                       "Host_Age")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file[,c("Instrument","Run")]
all(ena_file_wampach$run_accession %in% sra_file$Run)
ena_file_wampach <- merge(ena_file_wampach, sra_file,by = "Run")

# get time-points and individual ID:
ena_file_wampach <- separate(ena_file_wampach, col = "Sample.Name", into = c("time_point", "subject_ID"), sep = "[_]", remove = F)
ena_file_wampach$age <- NA
ena_file_wampach$age[ena_file_wampach$time_point == "V1"] <- 1
ena_file_wampach$age[ena_file_wampach$time_point == "V2"] <- 3
ena_file_wampach$age[ena_file_wampach$time_point == "V3"] <- 5
# remove vaginal and control samples
sort_out_samples_wam <- ena_file_wampach$Run[!ena_file_wampach$env_biome == "human gut"]
ena_file_wampach <- ena_file_wampach[ena_file_wampach$env_biome == "human gut",]
sort_out_samples_wam <- c(sort_out_samples_wam, ena_file_wampach$Run[ena_file_wampach$time_point == "M"])
ena_file_wampach <- ena_file_wampach[!ena_file_wampach$time_point == "M",]
# all individuals have several samples
multiple_samples <- unique(ena_file_wampach$subject_ID[duplicated(ena_file_wampach$subject_ID)])
length(ena_file_wampach$Run[!ena_file_wampach$subject_ID %in% multiple_samples])



metadata_wam_1_init <- read.table("suppl_1_init.txt", sep = '\t', header = T)
metadata_wam_1_init$Delivery.mode[metadata_wam_1_init$Delivery.mode == "Vaginal"] <- "v"
metadata_wam_1_init$Delivery.mode[grep("section", metadata_wam_1_init$Delivery.mode)] <- "c"
metadata_wam_1_init$Gender[metadata_wam_1_init$Gender == "Male"] <- "m"
metadata_wam_1_init$Gender[metadata_wam_1_init$Gender == "Female"] <- "f"
metadata_wam_1_init <- separate(metadata_wam_1_init, col = "Gestational.age....weeks.",
                                into = c("gestational_week", "gestational_day"),
                                sep = "[+]", remove = T, convert = T)
metadata_wam_1_init$gestational_day <- gsub("day.*", "", metadata_wam_1_init$gestational_day) %>% as.numeric() /7
metadata_wam_1_init$gestational_age <- metadata_wam_1_init$gestational_week + metadata_wam_1_init$gestational_day
metadata_wam_1_init[c("gestational_day", "gestational_week")] <- NULL

metadata_wam_1_d1 <- read.table("suppl_1_d1.txt", sep = '\t', header = T)
metadata_wam_1_d1$Type.of.Milk <- gsub("Breast", "breast", metadata_wam_1_d1$Type.of.Milk)
metadata_wam_1_d1$Type.of.Milk <- gsub("Formula", "formula", metadata_wam_1_d1$Type.of.Milk)
metadata_wam_1_d1$Type.of.Milk <- gsub("Combined", "mixed_b_f", metadata_wam_1_d1$Type.of.Milk)
metadata_wam_1_d1$Type.of.Milk <- gsub("not yet fed", "nothing", metadata_wam_1_d1$Type.of.Milk)

metadata_wam_1_d3 <- read.table("suppl_1_d3.txt", sep = '\t', header = T)
metadata_wam_1_d3$Type.of.Milk <- gsub("Breast", "breast", metadata_wam_1_d3$Type.of.Milk)
metadata_wam_1_d3$Type.of.Milk <- gsub("Formula", "formula", metadata_wam_1_d3$Type.of.Milk)
metadata_wam_1_d3$Type.of.Milk <- gsub("Combined", "mixed_b_f", metadata_wam_1_d3$Type.of.Milk)

metadata_wam_1_d5 <- read.table("suppl_1_d5.txt", sep = '\t', header = T)
metadata_wam_1_d5$Type.of.Milk <- gsub("Breast", "breast", metadata_wam_1_d5$Type.of.Milk)
metadata_wam_1_d5$Type.of.Milk <- gsub("Formula", "formula", metadata_wam_1_d5$Type.of.Milk)
metadata_wam_1_d5$Type.of.Milk <- gsub("Combined", "mixed_b_f", metadata_wam_1_d5$Type.of.Milk)

metadata_wam_1_mother <- read.table("suppl_1_init_mother.txt", sep = '\t', header = T)
metadata_wam_1_mother[c("Collected.samples", "Data.available", "Mother.age.range", "Mother.Ethnicity")] <- NULL

metadata_wam_4_d5 <- read.table("suppl_4_d5.txt", sep = '\t', header = T)

metadata_wam_4_m1 <- read.table("suppl_4_1m.txt", sep = '\t', header = T)
metadata_wam_4_m1$Feeding..milk. <- gsub("Breast", "breast", metadata_wam_4_m1$Feeding..milk.)
metadata_wam_4_m1$Feeding..milk. <- gsub("Formula", "formula", metadata_wam_4_m1$Feeding..milk.)
metadata_wam_4_m1$Feeding..milk. <- gsub("Mixed", "mixed_b_f", metadata_wam_4_m1$Feeding..milk.)

metadata_wam_4_m6 <- read.table("suppl_4_6m.txt", sep = '\t', header = T)
metadata_wam_4_m6$Feeding..milk. <- gsub("Breast, Mixed", "mixed_b_s", metadata_wam_4_m6$Feeding..milk.)
metadata_wam_4_m6$Feeding..milk. <- gsub("Breast", "breast", metadata_wam_4_m6$Feeding..milk.)
metadata_wam_4_m6$Feeding..milk. <- gsub("Formula", "formula", metadata_wam_4_m6$Feeding..milk.)
metadata_wam_4_m6$Feeding..milk. <- gsub("Mixed", "mixed_b_f", metadata_wam_4_m6$Feeding..milk.)

metadata_wam_4_y1 <- read.table("suppl_4_1y.txt", sep = '\t', header = T)
metadata_wam_4_y1$Feeding..milk. <- gsub("Breast, Formula, Mixed", "mixed_b_f_s", metadata_wam_4_y1$Feeding..milk.)
metadata_wam_4_y1$Feeding..milk. <- gsub("Breast, Mixed", "mixed_b_s", metadata_wam_4_y1$Feeding..milk.)
metadata_wam_4_y1$Feeding..milk. <- gsub("Formula, Mixed", "mixed_f_s", metadata_wam_4_y1$Feeding..milk.)
metadata_wam_4_y1$Feeding..milk. <- gsub("Breast", "breast", metadata_wam_4_y1$Feeding..milk.)
metadata_wam_4_y1$Feeding..milk. <- gsub("Formula", "formula", metadata_wam_4_y1$Feeding..milk.)


# no duplicates
sum(duplicated(ena_file_wampach$Sample.Name))

#combine data.frames
sum(!unique(ena_file_wampach$subject_ID) %in% metadata_wam_1_init$Study.ID)
sum(!unique(metadata_wam_1_init$Study.ID) %in% ena_file_wampach$subject_ID)
nrow(ena_file_wampach)
metadata_wam_combined <- merge(metadata_wam_1_init, ena_file_wampach, by.x = "Study.ID", by.y = "subject_ID")
nrow(metadata_wam_combined)
metadata_wam_combined$weight <- NA
metadata_wam_combined$weight[metadata_wam_combined$time_point == "V1"] <- 
   metadata_wam_combined$Birth.weight...g.[metadata_wam_combined$time_point == "V1"]



sum(!unique(metadata_wam_combined$Study.ID) %in% metadata_wam_1_mother$Study.ID)
sum(!unique(metadata_wam_1_mother$Study.ID) %in% metadata_wam_combined$Study.ID)
nrow(metadata_wam_combined)
metadata_wam_combined <- merge(metadata_wam_combined, metadata_wam_1_mother, by = "Study.ID")
nrow(metadata_wam_combined)

# add feeding and weight for each time point
metadata_wam_combined$food <- NA
metadata_wam_combined$recent_diarrhoea <- NA
metadata_wam_combined$food[metadata_wam_combined$time_point == "V1"] <- 
   metadata_wam_1_d1$Type.of.Milk[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V1"], metadata_wam_1_d1$Study.ID)]

metadata_wam_combined$food[metadata_wam_combined$time_point == "V2"] <- 
   metadata_wam_1_d3$Type.of.Milk[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V2"], metadata_wam_1_d3$Study.ID)]
metadata_wam_combined$weight[metadata_wam_combined$time_point == "V2"] <- 
   metadata_wam_1_d3$Weight...g.[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V2"], metadata_wam_1_d3$Study.ID)]

metadata_wam_combined$food[metadata_wam_combined$time_point == "V3"] <- 
   metadata_wam_1_d5$Type.of.Milk[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V3"], metadata_wam_1_d5$Study.ID)]
metadata_wam_combined$weight[metadata_wam_combined$time_point == "V3"] <- 
   metadata_wam_1_d5$Weight...g.[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V3"], metadata_wam_1_d5$Study.ID)]
metadata_wam_combined$recent_diarrhoea[metadata_wam_combined$time_point == "V3"] <- 
   metadata_wam_4_d5$Recent.Diarrhoea[match(metadata_wam_combined$Study.ID[metadata_wam_combined$time_point == "V3"], metadata_wam_4_d5$Study.ID)]

metadata_wam_combined[metadata_wam_combined == "Yes"] <- T
metadata_wam_combined[metadata_wam_combined == "No"] <- F
metadata_wam_combined[metadata_wam_combined == "No "] <- F

# remove unneeded columns
metadata_wam_combined[,c("Study.Allocation", "Cohort", "Gestational.age.category", "Small.for..gestational.age..SGA.",
                         "AvgSpotLen", "Collection_date", "Strep.positive..antenatal.screen", "Host_Age")] <- NULL


# unify colnames
metadata_wam_combined %<>% mutate(run_accession = Run,
                                  sample_ID = paste0(Sample.Name, "_wampach_2018"),
                                  subject_ID = paste0(Study.ID, "_wampach_2018"),
                                  birthmode = Delivery.mode,
                                  birth_weight = Birth.weight...g. / 1000,
                                  sex = Gender,
                                  birth_height = Length...cm.,
                                  country = "LUXEMBOURG",
                                  region = "LUXEMBOURG",
                                  study = "wampach_2018",
                                  lifestyle = "industrialized",
                                  family_ID = subject_ID,
                                  weight = weight / 1000,
                                  geographic_location_.latitude. = 49.611622,
                                  geographic_location_.longitude. = 6.131935,
                                  .keep = "unused")


metadata_wam_combined[c("antibiotics_any",
                        "antibiotics_before",
                        "antibiotics_one_week_before")] <- F
sort(colnames(metadata_wam_combined))

# move sorted out samples
sort_out_samples_wam <- unique(sort_out_samples_wam)
sort_out_files_wam <- c(paste(sort_out_samples_wam, "_1.fastq.gz", sep = ""),
                        paste(sort_out_samples_wam, "_2.fastq.gz", sep = ""),
                        paste(sort_out_samples_wam, ".fastq.gz", sep = ""))
sort_out_files_wam <- sort_out_files_wam[file.exists(sort_out_files_wam)]
file.move(sort_out_files_wam, "sorted_out/")

# move good files
final_samples_wam <- metadata_wam_combined$run_accession
final_files_wam <- c(paste(final_samples_wam, "_1.fastq.gz", sep = ""),
                     paste(final_samples_wam, "_2.fastq.gz", sep = ""),
                     paste(final_samples_wam, ".fastq.gz", sep = ""))
final_files_wam <- final_files_wam[file.exists(final_files_wam)]
file.move(final_files_wam, "fastq_files/")

# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_wam_combined$Bases, n_lines$read_counts[match(metadata_wam_combined$run_accession, n_lines$sample_ID)])  # 
metadata_wam_combined$read_count <- n_lines$read_counts[match(metadata_wam_combined$run_accession, n_lines$sample_ID)]

# write file
write.csv(metadata_wam_combined, "metadata_wam_2018_healthy.csv")


# clean up
rm(ena_file_wampach, metadata_wam_1_d1, metadata_wam_1_d3, metadata_wam_1_d5, metadata_wam_1_init,
   metadata_wam_1_mother, metadata_wam_4_d5, metadata_wam_4_m1, metadata_wam_4_m6, metadata_wam_4_y1,
   metadata_wam_combined, sort_out_files_wam, sort_out_samples_wam, multiple_samples,
   final_samples_wam, final_files_wam, n_lines, sra_file)




# Gasparrini	2019 #############################################################

setwd(paste(maindir, "shotgun/gasparrini_2019/", sep = ""))
ena_file_gasparrini <- read.table("filereport_read_run_PRJNA489090_tsv.txt", sep = "\t", header = T)
ena_file_gasparrini <- ena_file_gasparrini[c("study_accession", "run_accession",
                                             "library_name", "read_count", "sample_alias")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file[,c("Instrument","Run")]
all(ena_file_gasparrini$run_accession %in% sra_file$Run)
ena_file_gasparrini <- merge(ena_file_gasparrini, sra_file,by.x = "run_accession", by.y = "Run")


metadata_gas <- read.table("41564_2019_550_MOESM3_ESM.txt", sep = "\t", header = T)


# no duplicated samples
sum(duplicated(ena_file_gasparrini$sample_alias))

#combine data.frames
sum(!ena_file_gasparrini$sample_alias %in% metadata_gas$Sample)
sum(!metadata_gas$Sample %in% ena_file_gasparrini$sample_alias)  # 8 missing samples for 4 individuals
nrow(ena_file_gasparrini)
metadata_gas_combined <- merge(metadata_gas, ena_file_gasparrini, by.x = "Sample", by.y = "sample_alias")
nrow(metadata_gas_combined)
# sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
# sra_file <- sra_file[,c("Run", "geographic_location_.latitude.", "geographic_location_.longitude.")]

#unify colnames
metadata_gas_combined %<>% mutate(sample_ID = Sample,
                                  subject_ID = Individual,
                                  age = DOL,
                                  gestational_age = GestationalAgeWeeks,
                                  birth_weight = Birthweight/1000,
                                  sex = Gender,
                                  birthmode = DeliveryMode,
                                  geographic_location_.latitude. = 38.627000,
                                  geographic_location_.longitude. = -90.199389,
                                  .keep = "unused")

metadata_gas_combined$country <- "USA"
metadata_gas_combined$region <- "MISSOURI"
metadata_gas_combined$lifestyle <- "industrialized"
metadata_gas_combined$study <- "gasparrini_2019"

# rename entries
metadata_gas_combined$age <- floor(metadata_gas_combined$age)
metadata_gas_combined$sex[metadata_gas_combined$sex == "Female"] <- "f"
metadata_gas_combined$sex[metadata_gas_combined$sex == "Male"] <- "m"
metadata_gas_combined$birthmode[metadata_gas_combined$birthmode == "Cesarean"] <- "c"
metadata_gas_combined$birthmode[metadata_gas_combined$birthmode == "Vaginal"] <- "v"
metadata_gas_combined$NecrotizingEnterocolitis[metadata_gas_combined$NecrotizingEnterocolitis == "No"] <- 0
metadata_gas_combined$NecrotizingEnterocolitis[metadata_gas_combined$NecrotizingEnterocolitis == "Yes"] <- 1
metadata_gas_combined$NecrotizingEnterocolitis <- as.numeric(metadata_gas_combined$NecrotizingEnterocolitis)
metadata_gas_combined$Chorioamnionitis[metadata_gas_combined$Chorioamnionitis == "No"] <- 0
metadata_gas_combined$Chorioamnionitis[metadata_gas_combined$Chorioamnionitis == "Yes"] <- 1
metadata_gas_combined$Chorioamnionitis <- as.numeric(metadata_gas_combined$Chorioamnionitis)
metadata_gas_combined$GBS[metadata_gas_combined$GBS == "Negative"] <- 0
metadata_gas_combined$GBS[metadata_gas_combined$GBS == "Positive"] <- 1
metadata_gas_combined$GBS[metadata_gas_combined$GBS == "Unknown"] <- NA
metadata_gas_combined$GBS <- as.numeric(metadata_gas_combined$GBS)

metadata_gas_combined$food <- NA
metadata_gas_combined$food[metadata_gas_combined$age <= 75 &
                              metadata_gas_combined$EnteralFeeds_2mo == "Formula "] <- "formula"
metadata_gas_combined$food[metadata_gas_combined$age <= 75 &
                              metadata_gas_combined$EnteralFeeds_2mo == "Breastmilk and formula"] <- "mixed_b_f"
metadata_gas_combined$food[metadata_gas_combined$age <= 75 &
                              metadata_gas_combined$EnteralFeeds_2mo %in% c("Breastmilk ",
                                                                            "Breastmilk",
                                                                            "Breastmilk  ")] <- "mixed_b_f"
# delete cols
metadata_gas_combined[c("PostmenstrualAgeDays", "PostmenstrualAgeWeeks",
                        "GestationalAgeDays", "EnteralFeeds_2mo")] <- NULL
metadata_gas_combined$family_ID <- metadata_gas_combined$subject_ID
sort(colnames(metadata_gas_combined))

# split to preterms and terms
metadata_gas_combined_preterms <- metadata_gas_combined[metadata_gas_combined$Group != "Term",]
sort_out_samples_gas <- metadata_gas_combined_preterms$run_accession
metadata_gas_combined <- metadata_gas_combined[metadata_gas_combined$Group == "Term",]

# no antibiotics in terms:
metadata_gas_combined[c("antibiotics_any", "antibiotics_before", "antibiotics_one_week_before")] <- "n"


#remove cols with only zero
metadata_gas_combined <- metadata_gas_combined[, colSums(metadata_gas_combined != 0, na.rm = T) > 0]
metadata_gas_combined_preterms <- metadata_gas_combined_preterms[, colSums(metadata_gas_combined_preterms != 0, na.rm = T) > 0]


# all individuals have multiple samples
multiple_samples <- unique(metadata_gas_combined_preterms$subject_ID[duplicated(metadata_gas_combined_preterms$subject_ID)])
sum(!metadata_gas_combined_preterms$subject_ID %in% multiple_samples)
multiple_samples <- unique(metadata_gas_combined$subject_ID[duplicated(metadata_gas_combined$subject_ID)])
sum(!metadata_gas_combined$subject_ID %in% multiple_samples)

# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_gas_combined$read_count, n_lines$read_counts[match(metadata_gas_combined$run_accession, n_lines$sample_ID)])  # 


write.csv(metadata_gas_combined, "metadata_gas_2018_healthy.csv")
write.csv(metadata_gas_combined_preterms, "metadata_gas_2018_preterm.csv")

sort_out_samples_gas <- unique(sort_out_samples_gas)
sort_out_files_gas <- c(paste(sort_out_samples_gas, "_1.fastq.gz", sep = ""),
                        paste(sort_out_samples_gas, "_2.fastq.gz", sep = ""),
                        paste(sort_out_samples_gas, ".fastq.gz", sep = ""))
sort_out_files_gas <- sort_out_files_gas[file.exists(sort_out_files_gas)]
file.move(sort_out_files_gas, "sorted_out/")


# move good files
final_samples_gas <- metadata_gas_combined$run_accession
final_files_gas <- c(paste(final_samples_gas, "_1.fastq.gz", sep = ""),
                     paste(final_samples_gas, "_2.fastq.gz", sep = ""),
                     paste(final_samples_gas, ".fastq.gz", sep = ""))
final_files_gas <- final_files_gas[file.exists(final_files_gas)]
file.move(final_files_gas, "fastq_files/")

rm(ena_file_gasparrini, metadata_gas, metadata_gas_combined,
   metadata_gas_combined_preterms, multiple_samples, sort_out_files_gas,
   sort_out_samples_gas, final_samples_gas, final_files_gas, n_lines)





# Bhanu Busi	2021 #############################################################

setwd(paste(maindir, "shotgun/bhanu_busi_2021/", sep = ""))
ena_file_bhanu <- read.table("filereport_read_run_PRJNA595749_tsv.txt", sep = "\t", header = T)
ena_file_bhanu <- ena_file_bhanu[c("study_accession", "run_accession", "library_name",
                                   "read_count", "sample_alias", "sample_title")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file[,c("Instrument","Run")]
all(ena_file_bhanu$run_accession %in% sra_file$Run)
ena_file_bhanu <- merge(ena_file_bhanu, sra_file,by.x = "run_accession", by.y = "Run")

# get only fastq-files
sort_out_samples_bha <- ena_file_bhanu$run_accession[-grep("COSMIC", ena_file_bhanu$sample_alias)]
ena_file_bhanu <- ena_file_bhanu[grep("COSMIC", ena_file_bhanu$sample_alias),]
ena_file_bhanu <- separate(ena_file_bhanu, "sample_title",
                           into = c("subject_ID", "time_point"), sep = "V", remove = F)
# get time-points and individual ID:
ena_file_bhanu$age <- NA
ena_file_bhanu$age[ena_file_bhanu$time_point == "3"] <- 5
ena_file_bhanu$age[ena_file_bhanu$time_point == "4"] <- 30
ena_file_bhanu$age[ena_file_bhanu$time_point == "5"] <- 180
ena_file_bhanu$age[ena_file_bhanu$time_point == "6"] <- 365

# all individuals have more thanb one sample
multiple_samples <- unique(ena_file_bhanu$subject_ID[duplicated(ena_file_bhanu$subject_ID)])
length(ena_file_bhanu$run_accession[!ena_file_bhanu$subject_ID %in% multiple_samples])


metadata_bha_1_init <- read.table("suppl_1_init.txt", sep = '\t', header = T)
metadata_bha_1_init$Delivery.mode[metadata_bha_1_init$Delivery.mode == "Vaginal"] <- "v"
metadata_bha_1_init$Delivery.mode[grep("section", metadata_bha_1_init$Delivery.mode)] <- "c"
metadata_bha_1_init$Gender[metadata_bha_1_init$Gender == "Male"] <- "m"
metadata_bha_1_init$Gender[metadata_bha_1_init$Gender == "Female"] <- "f"
metadata_bha_1_init <- separate(metadata_bha_1_init, col = "Gestational.age....weeks.",
                                into = c("gestational_week", "gestational_day"),
                                sep = "[+]", remove = T, convert = T)
metadata_bha_1_init$gestational_day <- gsub("day.*", "", metadata_bha_1_init$gestational_day) %>% as.numeric() /7
metadata_bha_1_init$gestational_age <- metadata_bha_1_init$gestational_week + metadata_bha_1_init$gestational_day
metadata_bha_1_init[c("gestational_day", "gestational_week")] <- NULL

metadata_bha_1_d1 <- read.table("suppl_1_d1.txt", sep = '\t', header = T)
metadata_bha_1_d1$Type.of.Milk <- gsub("Breast", "breast", metadata_bha_1_d1$Type.of.Milk)
metadata_bha_1_d1$Type.of.Milk <- gsub("Formula", "formula", metadata_bha_1_d1$Type.of.Milk)
metadata_bha_1_d1$Type.of.Milk <- gsub("Combined", "mixed_b_f", metadata_bha_1_d1$Type.of.Milk)
metadata_bha_1_d1$Type.of.Milk <- gsub("not yet fed", "nothing", metadata_bha_1_d1$Type.of.Milk)

metadata_bha_1_d3 <- read.table("suppl_1_d3.txt", sep = '\t', header = T)
metadata_bha_1_d3$Type.of.Milk <- gsub("Breast", "breast", metadata_bha_1_d3$Type.of.Milk)
metadata_bha_1_d3$Type.of.Milk <- gsub("Formula", "formula", metadata_bha_1_d3$Type.of.Milk)
metadata_bha_1_d3$Type.of.Milk <- gsub("Combined", "mixed_b_f", metadata_bha_1_d3$Type.of.Milk)

metadata_bha_1_d5 <- read.table("suppl_1_d5.txt", sep = '\t', header = T)
metadata_bha_1_d5$Type.of.Milk <- gsub("Breast", "breast", metadata_bha_1_d5$Type.of.Milk)
metadata_bha_1_d5$Type.of.Milk <- gsub("Formula", "formula", metadata_bha_1_d5$Type.of.Milk)
metadata_bha_1_d5$Type.of.Milk <- gsub("Combined", "mixed_b_f", metadata_bha_1_d5$Type.of.Milk)

metadata_bha_1_mother <- read.table("suppl_1_init_mother.txt", sep = '\t', header = T)
metadata_bha_1_mother[c("Collected.samples", "Data.available", "Mother.age.range", "Mother.Ethnicity")] <- NULL

metadata_bha_4_d5 <- read.table("suppl_4_d5.txt", sep = '\t', header = T)

metadata_bha_4_m1 <- read.table("suppl_4_1m.txt", sep = '\t', header = T)
metadata_bha_4_m1$Feeding..milk. <- gsub("Breast", "breast", metadata_bha_4_m1$Feeding..milk.)
metadata_bha_4_m1$Feeding..milk. <- gsub("Formula", "formula", metadata_bha_4_m1$Feeding..milk.)
metadata_bha_4_m1$Feeding..milk. <- gsub("Mixed", "mixed_b_f", metadata_bha_4_m1$Feeding..milk.)

metadata_bha_4_m6 <- read.table("suppl_4_6m.txt", sep = '\t', header = T)
metadata_bha_4_m6$Feeding..milk. <- gsub("Breast, Mixed", "mixed_b_s", metadata_bha_4_m6$Feeding..milk.)
metadata_bha_4_m6$Feeding..milk. <- gsub("Breast", "breast", metadata_bha_4_m6$Feeding..milk.)
metadata_bha_4_m6$Feeding..milk. <- gsub("Formula", "formula", metadata_bha_4_m6$Feeding..milk.)
metadata_bha_4_m6$Feeding..milk. <- gsub("Mixed", "mixed_b_f", metadata_bha_4_m6$Feeding..milk.)

metadata_bha_4_y1 <- read.table("suppl_4_1y.txt", sep = '\t', header = T)
metadata_bha_4_y1$Feeding..milk. <- gsub("Breast, Formula, Mixed", "mixed_b_f_s", metadata_bha_4_y1$Feeding..milk.)
metadata_bha_4_y1$Feeding..milk. <- gsub("Breast, Mixed", "mixed_b_s", metadata_bha_4_y1$Feeding..milk.)
metadata_bha_4_y1$Feeding..milk. <- gsub("Formula, Mixed", "mixed_f_s", metadata_bha_4_y1$Feeding..milk.)
metadata_bha_4_y1$Feeding..milk. <- gsub("Breast", "breast", metadata_bha_4_y1$Feeding..milk.)
metadata_bha_4_y1$Feeding..milk. <- gsub("Formula", "formula", metadata_bha_4_y1$Feeding..milk.)


# no duplicates
sum(duplicated(ena_file_bhanu$sample_alias))

#combine data.frames
sum(!unique(ena_file_bhanu$subject_ID) %in% metadata_bha_1_init$Study.ID)
sum(!unique(metadata_bha_1_init$Study.ID) %in% ena_file_bhanu$subject_ID)
nrow(ena_file_bhanu)
metadata_bha_combined <- merge(metadata_bha_1_init, ena_file_bhanu, by.x = "Study.ID", by.y = "subject_ID")
nrow(metadata_bha_combined)
metadata_bha_combined$weight <- NA
metadata_bha_combined$weight[metadata_bha_combined$time_point == "1"] <- 
   metadata_bha_combined$Birth.weight...g.[metadata_bha_combined$time_point == "1"]


# not mother data for these samples
sum(!unique(metadata_bha_combined$Study.ID) %in% metadata_bha_1_mother$Study.ID)
sum(!unique(metadata_bha_1_mother$Study.ID) %in% metadata_bha_combined$Study.ID)


# add feeding and weight for each time point
metadata_bha_combined$food <- NA
metadata_bha_combined$recent_diarrhoea <- NA
metadata_bha_combined$antibiotics_before <- NA
metadata_bha_combined$antibiotics_before[metadata_bha_combined$time_point %in% c("1", "2", "3")] <- F
metadata_bha_combined$food[metadata_bha_combined$time_point == "1"] <- 
   metadata_bha_1_d1$Type.of.Milk[match(metadata_bha_combined$Study.ID[metadata_bha_combined$time_point == "1"], metadata_bha_1_d1$Study.ID)]

metadata_bha_combined$food[metadata_bha_combined$time_point == "2"] <- 
   metadata_bha_1_d3$Type.of.Milk[match(metadata_bha_combined$Study.ID[metadata_bha_combined$time_point == "2"], metadata_bha_1_d3$Study.ID)]
metadata_bha_combined$weight[metadata_bha_combined$time_point == "2"] <- 
   metadata_bha_1_d3$Weight...g.[match(metadata_bha_combined$Study.ID[metadata_bha_combined$time_point == "2"], metadata_bha_1_d3$Study.ID)]

metadata_bha_combined$food[metadata_bha_combined$time_point == "3"] <- 
   metadata_bha_1_d5$Type.of.Milk[match(metadata_bha_combined$Study.ID[metadata_bha_combined$time_point == "3"], metadata_bha_1_d5$Study.ID)]
metadata_bha_combined$weight[metadata_bha_combined$time_point == "3"] <- 
   metadata_bha_1_d5$Weight...g.[match(metadata_bha_combined$Study.ID[metadata_bha_combined$time_point == "3"], metadata_bha_1_d5$Study.ID)]
metadata_bha_combined$recent_diarrhoea[metadata_bha_combined$time_point == "3"] <- 
   metadata_bha_4_d5$Recent.Diarrhoea[match(metadata_bha_combined$Study.ID[metadata_bha_combined$time_point == "3"], metadata_bha_4_d5$Study.ID)]

metadata_bha_combined$food[metadata_bha_combined$time_point == "4"] <- 
   metadata_bha_4_m1$Feeding..milk.[match(metadata_bha_combined$Study.ID[metadata_bha_combined$time_point == "4"], metadata_bha_4_m1$Study.ID)]
metadata_bha_combined$weight[metadata_bha_combined$time_point == "4"] <- 
   metadata_bha_4_m1$Weight..g.[match(metadata_bha_combined$Study.ID[metadata_bha_combined$time_point == "4"], metadata_bha_4_m1$Study.ID)]
metadata_bha_combined$recent_diarrhoea[metadata_bha_combined$time_point == "4"] <- 
   metadata_bha_4_m1$Baby.diarrhoea[match(metadata_bha_combined$Study.ID[metadata_bha_combined$time_point == "4"], metadata_bha_4_m1$Study.ID)]
metadata_bha_combined$antibiotics_before[metadata_bha_combined$time_point == "4"] <- 
   metadata_bha_4_m1$Antibiotics.since.last.visit.[match(metadata_bha_combined$Study.ID[metadata_bha_combined$time_point == "4"], metadata_bha_4_m1$Study.ID)]

metadata_bha_combined$food[metadata_bha_combined$time_point == "5"] <- 
   metadata_bha_4_m6$Feeding..milk.[match(metadata_bha_combined$Study.ID[metadata_bha_combined$time_point == "5"], metadata_bha_4_m6$Study.ID)]
metadata_bha_combined$weight[metadata_bha_combined$time_point == "5"] <- 
   metadata_bha_4_m6$Weight..g.[match(metadata_bha_combined$Study.ID[metadata_bha_combined$time_point == "5"], metadata_bha_4_m6$Study.ID)]
metadata_bha_combined$recent_diarrhoea[metadata_bha_combined$time_point == "5"] <- 
   metadata_bha_4_m6$Recent.baby.Diarrhoea[match(metadata_bha_combined$Study.ID[metadata_bha_combined$time_point == "5"], metadata_bha_4_m6$Study.ID)]
metadata_bha_combined$antibiotics_before[metadata_bha_combined$time_point == "5"] <- 
   metadata_bha_4_m6$Antibiotics.since.last.visit.[match(metadata_bha_combined$Study.ID[metadata_bha_combined$time_point == "5"], metadata_bha_4_m6$Study.ID)]

metadata_bha_combined$food[metadata_bha_combined$time_point == "6"] <- 
   metadata_bha_4_y1$Feeding..milk.[match(metadata_bha_combined$Study.ID[metadata_bha_combined$time_point == "6"], metadata_bha_4_y1$Study.ID)]
metadata_bha_combined$weight[metadata_bha_combined$time_point == "6"] <- 
   metadata_bha_4_y1$Weight..g.[match(metadata_bha_combined$Study.ID[metadata_bha_combined$time_point == "6"], metadata_bha_4_y1$Study.ID)]
metadata_bha_combined$recent_diarrhoea[metadata_bha_combined$time_point == "6"] <- 
   metadata_bha_4_y1$Recent.baby.Diarrhoea[match(metadata_bha_combined$Study.ID[metadata_bha_combined$time_point == "6"], metadata_bha_4_y1$Study.ID)]
metadata_bha_combined$antibiotics_before[metadata_bha_combined$time_point == "6"] <- 
   metadata_bha_4_y1$Antibiotics.since.last.visit.[match(metadata_bha_combined$Study.ID[metadata_bha_combined$time_point == "6"], metadata_bha_4_y1$Study.ID)]

metadata_bha_combined$antibiotics_before[grep("Yes", metadata_bha_combined$antibiotics_before)] <- T

metadata_bha_combined[metadata_bha_combined == "Yes"] <- T
metadata_bha_combined[metadata_bha_combined == "No"] <- F
metadata_bha_combined[metadata_bha_combined == "No "] <- F

# remove unneeded columns
metadata_bha_combined[,c("Study.Allocation", "Cohort", "Gestational.age.category", "Small.for..gestational.age..SGA.",
                         "AvgSpotLen", "Collection_date", "Strep.positive..antenatal.screen", "Host_Age")] <- NULL


# unify colnames
metadata_bha_combined %<>% mutate(subject_ID = Study.ID,
                                  birthmode = Delivery.mode,
                                  birth_weight = Birth.weight...g. /1000,
                                  birth_height = Length...cm.,
                                  sex = Gender,
                                  sample_ID = sample_title,
                                  weight = weight / 1000,
                                  geographic_location_.latitude. = 49.611622,
                                  geographic_location_.longitude. = 6.131935,
                                  antibiotics_before = as.logical(antibiotics_before),
                                  .keep = "unused")

metadata_bha_combined$country <- "LUXEMBOURG"
metadata_bha_combined$region <- "LUXEMBOURG"
metadata_bha_combined$study <- "bhanu_busi_2021"
metadata_bha_combined$lifestyle <- "industrialized"
metadata_bha_combined$family_ID <- metadata_bha_combined$subject_ID
sort(colnames(metadata_bha_combined))


# move sorted out samples
sort_out_samples_bha <- unique(sort_out_samples_bha)
sort_out_files_bha <- c(paste(sort_out_samples_bha, "_1.fastq.gz", sep = ""),
                        paste(sort_out_samples_bha, "_2.fastq.gz", sep = ""),
                        paste(sort_out_samples_bha, ".fastq.gz", sep = ""))
sort_out_files_bha <- sort_out_files_bha[file.exists(sort_out_files_bha)]
file.move(sort_out_files_bha, "sorted_out/")

# move good files
final_samples_bha <- metadata_bha_combined$run_accession
final_files_bha <- c(paste(final_samples_bha, "_1.fastq.gz", sep = ""),
                     paste(final_samples_bha, "_2.fastq.gz", sep = ""),
                     paste(final_samples_bha, ".fastq.gz", sep = ""))
final_files_bha <- final_files_bha[file.exists(final_files_bha)]
file.move(final_files_bha, "fastq_files/")

# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_bha_combined$read_count, n_lines$read_counts[match(metadata_bha_combined$run_accession, n_lines$sample_ID)])  # 

# write file
write.csv(metadata_bha_combined, "metadata_bha_2021_healthy.csv")


# clean up
rm(ena_file_bhanu, metadata_bha_1_d1, metadata_bha_1_d3, metadata_bha_1_d5,
   metadata_bha_1_init, metadata_bha_1_mother, metadata_bha_4_d5, metadata_bha_4_m1,
   metadata_bha_4_m6, metadata_bha_4_y1, metadata_bha_combined, multiple_samples,
   sort_out_files_bha, sort_out_samples_bha, final_samples_bha, final_files_bha,
   n_lines)



# Wampach 2018 and Bhanu Busi 2021 combined ####################################


# setwd(paste(maindir, "shotgun/wampach_bhanu_busi_combined/", sep = ""))
# # load two ena files separately, combine then
# ena_file_cosmic_1 <- read.table("filereport_read_run_PRJNA595749_tsv.txt", sep = "\t", header = T)
# sort_out_samples_cosmic <- ena_file_cosmic_1$run_accession[-grep("COSMIC", ena_file_cosmic_1$sample_alias)]
# ena_file_cosmic_1 <- ena_file_cosmic_1[grep("COSMIC", ena_file_cosmic_1$sample_alias),]
# ena_file_cosmic_1 <- separate(ena_file_cosmic_1, "sample_title",
#                               into = c("subject_ID", "time_point"), sep = "V", remove = F)
# 
# ena_file_cosmic_2 <- read.table("filereport_read_run_PRJNA379120_tsv.txt", sep = "\t", header = T)
# sort_out_samples_cosmic <- c(sort_out_samples_cosmic, 
#                              ena_file_cosmic_2$run_accession[-grep("human stool", ena_file_cosmic_2$sample_title)])
# ena_file_cosmic_2 <- ena_file_cosmic_2[grep("human stool", ena_file_cosmic_2$sample_title),]
# sort_out_samples_cosmic <- c(sort_out_samples_cosmic, 
#                              ena_file_cosmic_2$run_accession[grep("16S", ena_file_cosmic_2$library_name)])
# ena_file_cosmic_2 <- ena_file_cosmic_2[-grep("16S", ena_file_cosmic_2$library_name),]
# sort_out_samples_cosmic <- c(sort_out_samples_cosmic, 
#                              ena_file_cosmic_2$run_accession[grep("18S", ena_file_cosmic_2$library_name)])
# ena_file_cosmic_2 <- ena_file_cosmic_2[-grep("18S", ena_file_cosmic_2$library_name),]
# sort_out_samples_cosmic <- c(sort_out_samples_cosmic, 
#                              ena_file_cosmic_2$run_accession[grep("M[_]", ena_file_cosmic_2$sample_alias)])
# ena_file_cosmic_2 <- ena_file_cosmic_2[-grep("M[_]", ena_file_cosmic_2$sample_alias),]
# ena_file_cosmic_2 <- separate(ena_file_cosmic_2, "sample_alias",
#                               into = c("time_point", "subject_ID"), sep = "[_]", remove = F)
# ena_file_cosmic_2$time_point <- gsub("V", "", ena_file_cosmic_2$time_point)
# 
# ena_file_cosmic <- rbind(ena_file_cosmic_1, ena_file_cosmic_2)
# ena_file_cosmic <- ena_file_cosmic[c("study_accession", "run_accession", "library_name",
#                                    "read_count", "sample_alias", "sample_title", "time_point", "subject_ID")]
# 
# # get time-points and individual ID:
# ena_file_cosmic$age <- NA
# ena_file_cosmic$age[ena_file_cosmic$time_point == "1"] <- 1
# ena_file_cosmic$age[ena_file_cosmic$time_point == "2"] <- 3
# ena_file_cosmic$age[ena_file_cosmic$time_point == "3"] <- 5
# ena_file_cosmic$age[ena_file_cosmic$time_point == "4"] <- 30
# ena_file_cosmic$age[ena_file_cosmic$time_point == "5"] <- 180
# ena_file_cosmic$age[ena_file_cosmic$time_point == "6"] <- 365
# 
# # all individuals have more than one sample
# multiple_samples <- unique(ena_file_cosmic$subject_ID[duplicated(ena_file_cosmic$subject_ID)])
# length(ena_file_cosmic$run_accession[!ena_file_cosmic$subject_ID %in% multiple_samples])
# 
# 
# metadata_cosmic_1_init <- read.table("suppl_1_init.txt", sep = '\t', header = T)
# metadata_cosmic_1_init$Delivery.mode[metadata_cosmic_1_init$Delivery.mode == "Vaginal"] <- "v"
# metadata_cosmic_1_init$Delivery.mode[grep("section", metadata_cosmic_1_init$Delivery.mode)] <- "c"
# metadata_cosmic_1_init$Gender[metadata_cosmic_1_init$Gender == "Male"] <- "m"
# metadata_cosmic_1_init$Gender[metadata_cosmic_1_init$Gender == "Female"] <- "f"
# metadata_cosmic_1_init <- separate(metadata_cosmic_1_init, col = "Gestational.age....weeks.",
#                                 into = c("gestational_week", "gestational_day"),
#                                 sep = "[+]", remove = T, convert = T)
# metadata_cosmic_1_init$gestational_day <- gsub("day.*", "", metadata_cosmic_1_init$gestational_day) %>% as.numeric() /7
# metadata_cosmic_1_init$gestational_age <- metadata_cosmic_1_init$gestational_week + metadata_cosmic_1_init$gestational_day
# metadata_cosmic_1_init[c("gestational_day", "gestational_week", "Gestational.age....weeks.")] <- NULL
# 
# metadata_cosmic_1_d1 <- read.table("suppl_1_d1.txt", sep = '\t', header = T)
# metadata_cosmic_1_d3 <- read.table("suppl_1_d3.txt", sep = '\t', header = T)
# metadata_cosmic_1_d5 <- read.table("suppl_1_d5.txt", sep = '\t', header = T)
# 
# metadata_cosmic_1_mother <- read.table("suppl_1_init_mother.txt", sep = '\t', header = T)
# metadata_cosmic_1_mother[c("Collected.samples", "Data.available", "Mother.age.range", "Mother.Ethnicity")] <- NULL
# 
# metadata_cosmic_4_d5 <- read.table("suppl_4_d5.txt", sep = '\t', header = T)
# metadata_cosmic_4_m1 <- read.table("suppl_4_1m.txt", sep = '\t', header = T)
# metadata_cosmic_4_m6 <- read.table("suppl_4_6m.txt", sep = '\t', header = T)
# metadata_cosmic_4_y1 <- read.table("suppl_4_1y.txt", sep = '\t', header = T)
# 
# 
# # no duplicates
# sum(duplicated(ena_file_cosmic$sample_alias))
# 
# #combine data.frames
# sum(!unique(ena_file_cosmic$subject_ID) %in% metadata_cosmic_1_init$Study.ID)
# sum(!unique(metadata_cosmic_1_init$Study.ID) %in% ena_file_cosmic$subject_ID)  # no data for some individuals
# nrow(ena_file_cosmic)
# metadata_cosmic_combined <- merge(metadata_cosmic_1_init, ena_file_cosmic, by.x = "Study.ID", by.y = "subject_ID")
# nrow(metadata_cosmic_combined)
# metadata_cosmic_combined$weight <- NA
# metadata_cosmic_combined$weight[metadata_cosmic_combined$time_point == "1"] <- 
#    metadata_cosmic_combined$Birth.weight...g.[metadata_cosmic_combined$time_point == "1"]
# 
# 
# # not mother data for some samples
# sum(!unique(metadata_cosmic_combined$Study.ID) %in% metadata_cosmic_1_mother$Study.ID)
# sum(!unique(metadata_cosmic_1_mother$Study.ID) %in% metadata_cosmic_combined$Study.ID)  # two mothers without samples
# nrow(metadata_cosmic_combined)
# metadata_cosmic_combined <- merge(metadata_cosmic_combined, metadata_cosmic_1_mother, by = "Study.ID", all.x = T)
# nrow(metadata_cosmic_combined)
# 
# # add feeding and weight for each time point
# metadata_cosmic_combined$food <- NA
# metadata_cosmic_combined$recent_diarrhoea <- NA
# metadata_cosmic_combined$antibiotics_before <- NA
# metadata_cosmic_combined$antibiotics_before[metadata_cosmic_combined$time_point %in% c("1", "2", "3")] <- "n"
# metadata_cosmic_combined$food[metadata_cosmic_combined$time_point == "1"] <- 
#    metadata_cosmic_1_d1$Type.of.Milk[match(metadata_cosmic_combined$Study.ID[metadata_cosmic_combined$time_point == "1"], metadata_cosmic_1_d1$Study.ID)]
# 
# metadata_cosmic_combined$food[metadata_cosmic_combined$time_point == "2"] <- 
#    metadata_cosmic_1_d3$Type.of.Milk[match(metadata_cosmic_combined$Study.ID[metadata_cosmic_combined$time_point == "2"], metadata_cosmic_1_d3$Study.ID)]
# metadata_cosmic_combined$weight[metadata_cosmic_combined$time_point == "2"] <- 
#    metadata_cosmic_1_d3$Weight...g.[match(metadata_cosmic_combined$Study.ID[metadata_cosmic_combined$time_point == "2"], metadata_cosmic_1_d3$Study.ID)]
# 
# metadata_cosmic_combined$food[metadata_cosmic_combined$time_point == "3"] <- 
#    metadata_cosmic_1_d5$Type.of.Milk[match(metadata_cosmic_combined$Study.ID[metadata_cosmic_combined$time_point == "3"], metadata_cosmic_1_d5$Study.ID)]
# metadata_cosmic_combined$weight[metadata_cosmic_combined$time_point == "3"] <- 
#    metadata_cosmic_1_d5$Weight...g.[match(metadata_cosmic_combined$Study.ID[metadata_cosmic_combined$time_point == "3"], metadata_cosmic_1_d5$Study.ID)]
# metadata_cosmic_combined$recent_diarrhoea[metadata_cosmic_combined$time_point == "3"] <- 
#    metadata_cosmic_4_d5$Recent.Diarrhoea[match(metadata_cosmic_combined$Study.ID[metadata_cosmic_combined$time_point == "3"], metadata_cosmic_4_d5$Study.ID)]
# 
# metadata_cosmic_combined$food[metadata_cosmic_combined$time_point == "4"] <- 
#    metadata_cosmic_4_m1$Feeding..milk.[match(metadata_cosmic_combined$Study.ID[metadata_cosmic_combined$time_point == "4"], metadata_cosmic_4_m1$Study.ID)]
# metadata_cosmic_combined$weight[metadata_cosmic_combined$time_point == "4"] <- 
#    metadata_cosmic_4_m1$Weight..g.[match(metadata_cosmic_combined$Study.ID[metadata_cosmic_combined$time_point == "4"], metadata_cosmic_4_m1$Study.ID)]
# metadata_cosmic_combined$recent_diarrhoea[metadata_cosmic_combined$time_point == "4"] <- 
#    metadata_cosmic_4_m1$Baby.diarrhoea[match(metadata_cosmic_combined$Study.ID[metadata_cosmic_combined$time_point == "4"], metadata_cosmic_4_m1$Study.ID)]
# metadata_cosmic_combined$antibiotics_before[metadata_cosmic_combined$time_point == "4"] <- 
#    metadata_cosmic_4_m1$Antibiotics.since.last.visit.[match(metadata_cosmic_combined$Study.ID[metadata_cosmic_combined$time_point == "4"], metadata_cosmic_4_m1$Study.ID)]
# 
# metadata_cosmic_combined$food[metadata_cosmic_combined$time_point == "5"] <- 
#    metadata_cosmic_4_m6$Feeding..milk.[match(metadata_cosmic_combined$Study.ID[metadata_cosmic_combined$time_point == "5"], metadata_cosmic_4_m6$Study.ID)]
# metadata_cosmic_combined$weight[metadata_cosmic_combined$time_point == "5"] <- 
#    metadata_cosmic_4_m6$Weight..g.[match(metadata_cosmic_combined$Study.ID[metadata_cosmic_combined$time_point == "5"], metadata_cosmic_4_m6$Study.ID)]
# metadata_cosmic_combined$recent_diarrhoea[metadata_cosmic_combined$time_point == "5"] <- 
#    metadata_cosmic_4_m6$Recent.baby.Diarrhoea[match(metadata_cosmic_combined$Study.ID[metadata_cosmic_combined$time_point == "5"], metadata_cosmic_4_m6$Study.ID)]
# metadata_cosmic_combined$antibiotics_before[metadata_cosmic_combined$time_point == "5"] <- 
#    metadata_cosmic_4_m6$Antibiotics.since.last.visit.[match(metadata_cosmic_combined$Study.ID[metadata_cosmic_combined$time_point == "5"], metadata_cosmic_4_m6$Study.ID)]
# 
# metadata_cosmic_combined$food[metadata_cosmic_combined$time_point == "6"] <- 
#    metadata_cosmic_4_y1$Feeding..milk.[match(metadata_cosmic_combined$Study.ID[metadata_cosmic_combined$time_point == "6"], metadata_cosmic_4_y1$Study.ID)]
# metadata_cosmic_combined$weight[metadata_cosmic_combined$time_point == "6"] <- 
#    metadata_cosmic_4_y1$Weight..g.[match(metadata_cosmic_combined$Study.ID[metadata_cosmic_combined$time_point == "6"], metadata_cosmic_4_y1$Study.ID)]
# metadata_cosmic_combined$recent_diarrhoea[metadata_cosmic_combined$time_point == "6"] <- 
#    metadata_cosmic_4_y1$Recent.baby.Diarrhoea[match(metadata_cosmic_combined$Study.ID[metadata_cosmic_combined$time_point == "6"], metadata_cosmic_4_y1$Study.ID)]
# metadata_cosmic_combined$antibiotics_before[metadata_cosmic_combined$time_point == "6"] <- 
#    metadata_cosmic_4_y1$Antibiotics.since.last.visit.[match(metadata_cosmic_combined$Study.ID[metadata_cosmic_combined$time_point == "6"], metadata_cosmic_4_y1$Study.ID)]
# 
# metadata_cosmic_combined$antibiotics_before[grep("Yes", metadata_cosmic_combined$antibiotics_before)] <- "y"
# 
# 
# metadata_cosmic_combined[metadata_cosmic_combined == "Yes"] <- "y"
# metadata_cosmic_combined[metadata_cosmic_combined == "No"] <- "n"
# metadata_cosmic_combined[metadata_cosmic_combined == "No "] <- "n"
# metadata_cosmic_combined$food <- gsub("Breast, Formula, Mixed", "mixed_b_f_s", metadata_cosmic_combined$food)
# metadata_cosmic_combined$food <- gsub("Formula, Mixed", "mixed_f_s", metadata_cosmic_combined$food)
# metadata_cosmic_combined$food <- gsub("Mixed", "mixed_b_f", metadata_cosmic_combined$food)
# metadata_cosmic_combined$food <- gsub("Breast", "breast", metadata_cosmic_combined$food)
# metadata_cosmic_combined$food <- gsub("Formula", "formula", metadata_cosmic_combined$food)
# metadata_cosmic_combined$food <- gsub("Combined", "mixed_b_f", metadata_cosmic_combined$food)
# 
# # remove unneeded columns
# metadata_cosmic_combined[c("Study.Allocation", "Cohort", "Gestational.age.category",
#                            "Small.for..gestational.age..SGA.", "Strep.positive..antenatal.screen")] <- NULL
# 
# # unify colnames
# metadata_bha_combined %<>% mutate(subject_ID = Study.ID,
#                                   birthmode = Delivery.mode,
#                                   birth_weight = Birth.weight...g. /1000,
#                                   birth_height = Length...cm.,
#                                   sex = Gender,
#                                   sample_ID = sample_title,
#                                   weight = weight / 1000,
#                                   geographic_location_.latitude. = 49.611622,
#                                   geographic_location_.longitude. = 6.131935,
#                                   .keep = "unused")
# 
# metadata_cosmic_combined$country <- "LUXEMBOURG"
# metadata_cosmic_combined$region <- "LUXEMBOURG"
# metadata_cosmic_combined$study <- "cosmic"
# metadata_cosmic_combined$lifestyle <- "industrialized"
# metadata_cosmic_combined$family_ID <- metadata_cosmic_combined$subject_ID
# sort(colnames(metadata_cosmic_combined))
# 
# 
# # move sorted out samples
# sort_out_samples_cosmic <- unique(sort_out_samples_cosmic)
# sort_out_files_cosmic <- c(paste(sort_out_samples_cosmic, "_1.fastq.gz", sep = ""),
#                         paste(sort_out_samples_cosmic, "_2.fastq.gz", sep = ""),
#                         paste(sort_out_samples_cosmic, ".fastq.gz", sep = ""))
# sort_out_files_cosmic <- sort_out_files_cosmic[file.exists(sort_out_files_cosmic)]
# file.move(sort_out_files_cosmic, "sorted_out/")
# 
# # move good files
# final_samples_cosmic <- metadata_cosmic_combined$run_accession
# final_files_cosmic <- c(paste(final_samples_cosmic, "_1.fastq.gz", sep = ""),
#                      paste(final_samples_cosmic, "_2.fastq.gz", sep = ""),
#                      paste(final_samples_cosmic, ".fastq.gz", sep = ""))
# final_files_cosmic <- final_files_cosmic[file.exists(final_files_cosmic)]
# file.move(final_files_cosmic, "fastq_files/")
# 
# # check read numbers
# n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
# plot(metadata_cosmic_combined$read_count, n_lines$read_counts[match(metadata_cosmic_combined$run_accession, n_lines$sample_ID)])  # 
# 
# 
# # write file
# write.csv(metadata_cosmic_combined, "metadata_cosmic_healthy.csv")
# 
# 
# # clean up
# rm(ena_file_cosmic, ena_file_cosmic_1, ena_file_cosmic_2, metadata_cosmic_1_d1,
#    metadata_cosmic_1_d3, metadata_cosmic_1_d5, metadata_cosmic_1_init, metadata_cosmic_1_mother,
#    metadata_cosmic_4_d5, metadata_cosmic_4_m1, metadata_cosmic_4_m6, metadata_cosmic_4_y1,
#    metadata_cosmic_combined, multiple_samples, sort_out_files_cosmic, sort_out_samples_cosmic,
#    final_samples_cosmic, final_files_cosmic, n_lines)



# Baumann-Dudenhoeffer	2018 ###################################################
# waiting for response
# only SraRunTable available
setwd(paste(maindir, "shotgun/baumann_dudenhoeffer_2018", sep = ""))
ena_file_baumann <- read.table("filereport_read_run_PRJNA473126_tsv.txt", sep = "\t", header = T)
ena_file_baumann <- ena_file_baumann[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]

metadata_bau <- read.table("SraRunTable.txt", sep = ",", header = T)
# remove columns
metadata_bau[,c("Assay.Type", "BioSampleModel", "Bytes", "Center.Name", "Collection_Date",
                "Consent", "DATASTORE.filetype", "DATASTORE.provider", "DATASTORE.region",
                "env_biome", "env_feature", "env_material", 
                "geo_loc_name_country_continent", "Host", "host_tissue_sampled",
                "LibraryLayout", "LibrarySelection",
                "LibrarySource", "medic_hist_perform", "misc_param", "Organism",
                "Platform", "ReleaseDate", "samp_store_temp", "host_age_days",
                "Host_Age", "host_body_product", "geo_loc_name")] <- NULL

# 45 duplicated samples
sum(duplicated(metadata_bau$Sample.Name))

# no individuals with only one sample
multiple_samples <- unique(metadata_bau$host_subject_id[duplicated(metadata_bau$host_subject_id)])
metadata_bau$Run[!metadata_bau$host_subject_id %in% multiple_samples]


# unify colnames
metadata_bau %<>% mutate(run_accession = Run,
                         country = geo_loc_name_country,
                         age = host_age_days_at_time_of_survey,
                         antibiotics_one_week_before = host_antimicrobial_in_last_7days,
                         birth_weight = host_birthweight_g / 1000,
                         birthmode = host_del_route,
                         family_ID = host_family,
                         gestational_age = host_gestationalage_atdelivery,
                         antibiotics_any = host_lifetime_antibiotics_exposure,
                         sex = host_sex,
                         subject_ID = host_subject_id,
                         sample_ID = Sample.Name,
                         geographic_location_.latitude. = 38.627000,
                         geographic_location_.longitude. = -90.199389,
                         food = Host_diet,
                         .keep = "unused")

metadata_bau$region <- "MISSOURI"
metadata_bau$study <- "baumann_dudenhoeffer_2018"
metadata_bau$lifestyle <- "industrialized"
sort(colnames(metadata_bau))

# feeding column
metadata_bau$food[metadata_bau$food %in% c("Exclusively Breastfed\\, N/A",
                                           "Exclusively Breastfed\\, N/A\\, Lactose"
                                           )] <- "breast"
metadata_bau$food[metadata_bau$food %in% c("Mostly Breastfed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Maltodextrin\\, Polydextrose",
                                           "Mostly Breastfed\\, Hydrolyzed Cow''s Milk Formula\\, Corn Syrup",
                                           "Mostly Breastfed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose",
                                           "Exclusively Breastfed\\, Hydrolyzed Cow''s Milk Formula\\, Corn Syrup",
                                           "Mostly Breastfed\\, Cow''s Milk Formula\\, Lactose",
                                           "Mostly Breastfed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Polydextrose",
                                           "Mostly Breastfed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Fructo-oligosaccharides\\, Lactose\\, Polydextrose"
                                           )] <- "mixed_b_f"
metadata_bau$food[metadata_bau$food %in% c("Exclusively Breastfed\\, N/A\\, Cereal/Starch\\, Fruits/Vegetables",
                                           "Exclusively Breastfed\\, N/A\\, Cereal/Starch",
                                           "Exclusively Breastfed\\, N/A\\, Cereal/Starch\\, Fruits/Vegetables\\, Yogurt"
                                           )] <- "mixed_b_s"
metadata_bau$food[metadata_bau$food %in% c("Mostly Breastfed\\, Hydrolyzed Cow''s Milk Formula\\, Corn Syrup\\, Cereal/Starch\\, Fruits/Vegetables",
                                           "Mostly Breastfed\\, Cow''s Milk Formula\\, Lactose\\, Juice/Sweetened Drinks",
                                           "Mostly Breastfed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Cereal/Starch\\, Fruits/Vegetables",
                                           "Mostly Breastfed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Maltodextrin\\, Polydextrose\\, Juice/Sweetened Drinks",
                                           "Mostly Breastfed\\, Cow''s Milk Formula\\, Lactose\\, Cereal/Starch\\, Fruits/Vegetables",
                                           "Mostly Breastfed\\, Cow''s Milk Formula\\, Lactose\\, Fruits/Vegetables",
                                           "Exclusively Breastfed\\, N/A\\, Galacto-oligosaccharides\\, Lactose\\, Polydextrose\\, Cereal/Starch\\, Fruits/Vegetables",
                                           "Mostly Breastfed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Fructo-oligosaccharides\\, Lactose\\, Maltodextrin",
                                           "Mostly Breastfed\\, Cow''s Milk Formula\\, Lactose\\, Cereal/Starch",
                                           "Mostly Breastfed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Cereal/Starch",
                                           "Exclusively Breastfed\\, N/A\\, Juice/Sweetened Drinks"
                                           )] <- "mixed_b_f_s"
metadata_bau$food[metadata_bau$food %in% c("Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose",
                                           "Exclusively Formula-Fed\\, Soy Formula\\, Fructo-oligosaccharides\\, Corn Syrup",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Fructo-oligosaccharides\\, Lactose",
                                           "Exclusively Formula-Fed\\, Hydrolyzed Cow''s Milk Formula\\, Corn Syrup",
                                           "Mostly Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Lactose",
                                           "Mostly Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Polydextrose",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Polydextrose",
                                           "Mostly Formula-Fed\\, Hydrolyzed Cow''s Milk Formula\\, Corn Syrup",
                                           "Mostly Formula-Fed\\, Hydrolyzed Cow''s Milk Formula\\, Lactose\\, Maltodextrin",
                                           "Mostly Formula-Fed\\, Hydrolyzed Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Corn Syrup",
                                           "Mostly Formula-Fed\\, Hydrolyzed Cow''s Milk Formula\\, Corn Syrup\\, Lactose",
                                           "Exclusively Formula-Fed\\, Soy Formula\\, Corn Syrup",
                                           "Exclusively Formula-Fed\\, Soy Formula\\, Fructo-oligosaccharides",
                                           "Exclusively Formula-Fed\\, Hydrolyzed Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Corn Syrup\\, Lactose\\, Polydextrose",
                                           "Exclusively Formula-Fed\\, Hydrolyzed Cow''s Milk Formula\\, Corn Syrup\\, Lactose",
                                           "Mostly Formula-Fed\\, Hydrolyzed Cow''s Milk Formula",
                                           "Mostly Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Fructo-oligosaccharides\\, Lactose",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Maltodextrin\\, Polydextrose",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Fructo-oligosaccharides\\, Corn Syrup\\, Lactose\\, Maltodextrin\\, Polydextrose",
                                           "Mostly Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Maltodextrin\\, Polydextrose",
                                           "Exclusively Formula-Fed\\, Hydrolyzed Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Corn Syrup\\, Lactose",
                                           "Mostly Formula-Fed\\, Cow''s Milk Formula\\, Lactose"
                                           )] <- "formula"
metadata_bau$food[metadata_bau$food %in% c("Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Lactose\\, Cereal/Starch\\, Fruits/Vegetables\\, Meat/Fish/Eggs",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Lactose\\, Cereal/Starch\\, Fruits/Vegetables\\, Yogurt",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Lactose\\, Cereal/Starch",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Lactose\\, Cereal/Starch\\, Fruits/Vegetables",
                                           "Exclusively Formula-Fed\\, Soy Formula\\, Fructo-oligosaccharides\\, Corn Syrup\\, Cereal/Starch",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Fructo-oligosaccharides\\, Lactose\\, Cereal/Starch\\, Fruits/Vegetables\\, Meat/Fish/Eggs",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Fructo-oligosaccharides\\, Lactose\\, Cereal/Starch\\, Fruits/Vegetables",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Fructo-oligosaccharides\\, Lactose\\, Cereal/Starch",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Lactose\\, Cereal/Starch\\, Fruits/Vegetables\\, Sweets",
                                           "Exclusively Formula-Fed\\, Hydrolyzed Cow''s Milk Formula\\, Corn Syrup\\, Cereal/Starch\\, Fruits/Vegetables",
                                           "Mostly Formula-Fed\\, Hydrolyzed Cow''s Milk Formula\\, Juice/Sweetened Drinks",
                                           "Exclusively Formula-Fed\\, Hydrolyzed Cow''s Milk Formula\\, Corn Syrup\\, Cereal/Starch\\, Fruits/Vegetables\\, Juice/Sweetened Drinks",
                                           "Mostly Formula-Fed\\, Hydrolyzed Cow''s Milk Formula\\, Corn Syrup\\, Juice/Sweetened Drinks",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Corn Syrup\\, Lactose\\, Cereal/Starch\\, Fruits/Vegetables",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Cereal/Starch\\, Fruits/Vegetables",
                                           "Exclusively Formula-Fed\\, Hydrolyzed Cow''s Milk Formula\\, Corn Syrup\\, Fruits/Vegetables\\, Meat/Fish/Eggs\\, Juice/Sweetened Drinks",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Lactose\\, Cereal/Starch\\, Juice/Sweetened Drinks",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Cereal/Starch",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Lactose\\, Cereal/Starch\\, Fruits/Vegetables\\, Juice/Sweetened Drinks",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Fructo-oligosaccharides\\, Lactose\\, Juice/Sweetened Drinks",
                                           "Whole Milk Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Polydextrose\\, Cereal/Starch\\, Fruits/Vegetables\\, Whole Milk",
                                           "Exclusively Formula-Fed\\, Hydrolyzed Cow''s Milk Formula\\, Lactose\\, Maltodextrin\\, Fruits/Vegetables",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Lactose\\, Cereal/Starch\\, Fruits/Vegetables\\, Meat/Fish/Eggs\\, Juice/Sweetened Drinks",
                                           "Mostly Formula-Fed\\, Soy Formula\\, Fructo-oligosaccharides\\, Corn Syrup\\, Cereal/Starch\\, Fruits/Vegetables",
                                           "Exclusively Formula-Fed\\, Soy Formula\\, Fructo-oligosaccharides\\, Corn Syrup\\, Cereal/Starch\\, Fruits/Vegetables",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Polydextrose\\, Fruits/Vegetables\\, Juice/Sweetened Drinks",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Polydextrose\\, Cereal/Starch\\, Fruits/Vegetables\\, Juice/Sweetened Drinks",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Polydextrose\\, Cereal/Starch\\, Fruits/Vegetables",
                                           "Exclusively Formula-Fed\\, Soy Formula\\, Fructo-oligosaccharides\\, Corn Syrup\\, Cereal/Starch\\, Fruits/Vegetables\\, Meat/Fish/Eggs\\, Juice/Sweetened Drinks",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Polydextrose\\, Cereal/Starch",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Lactose\\, Cereal/Starch\\, Yogurt (Inactive)",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Lactose\\, Cereal/Starch\\, Yogurt (Inactive)\\, Juice/Sweetened Drinks",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Lactose\\, Cereal/Starch\\, Fruits/Vegetables\\, Yogurt (Inactive)\\, Juice/Sweetened Drinks",
                                           "Exclusively Formula-Fed\\, Soy Formula\\, Fructo-oligosaccharides\\, Corn Syrup\\, Fruits/Vegetables",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Lactose\\, Cereal/Starch\\, Fruits/Vegetables\\, Meat/Fish/Eggs\\, Yogurt (Inactive)\\, Juice/Sweetened Drinks",
                                           "Exclusively Formula-Fed\\, Hydrolyzed Cow''s Milk Formula\\, Corn Syrup\\, Fruits/Vegetables\\, Yogurt (Inactive)\\, Juice/Sweetened Drinks",
                                           "Exclusively Formula-Fed\\, Hydrolyzed Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Corn Syrup\\, Lactose\\, Polydextrose\\, Cereal/Starch\\, Fruits/Vegetables",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Cereal/Starch\\, Fruits/Vegetables\\, Meat/Fish/Eggs\\, Whole Milk",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Corn Syrup\\, Lactose\\, Maltodextrin\\, Polydextrose\\, Cereal/Starch\\, Fruits/Vegetables\\, Meat/Fish/Eggs",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Cereal/Starch\\, Fruits/Vegetables\\, Meat/Fish/Eggs\\, Yogurt\\, Whole Milk",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Corn Syrup\\, Lactose\\, Maltodextrin\\, Polydextrose\\, Cereal/Starch",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Maltodextrin\\, Polydextrose\\, Cereal/Starch\\, Fruits/Vegetables",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Polydextrose\\, Cereal/Starch\\, Fruits/Vegetables\\, Meat/Fish/Eggs\\, Juice/Sweetened Drinks",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Lactose\\, Cereal/Starch\\, Fruits/Vegetables\\, Meat/Fish/Eggs\\, Yogurt",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Juice/Sweetened Drinks",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Polydextrose\\, Cereal/Starch\\, Fruits/Vegetables\\, Meat/Fish/Eggs",
                                           "Mostly Formula-Fed\\, Cow''s Milk Formula\\, Lactose\\, Cereal/Starch",
                                           "Whole Milk Fed\\, Cow''s Milk Formula\\, Cereal/Starch\\, Fruits/Vegetables\\, Meat/Fish/Eggs\\, Yogurt\\, Whole Milk\\, Juice/Sweetened Drinks\\, Sweets",
                                           "Mostly Formula-Fed\\, Cow''s Milk Formula\\, Galacto-oligosaccharides\\, Lactose\\, Cereal/Starch\\, Fruits/Vegetables",
                                           "Exclusively Formula-Fed\\, Cow''s Milk Formula\\, Lactose\\, Fruits/Vegetables",
                                           "Mostly Formula-Fed\\, Cow''s Milk Formula\\, Lactose\\, Cereal/Starch\\, Fruits/Vegetables\\, Juice/Sweetened Drinks"
                                           )] <- "mixed_f_s"

# unify col entries columns
metadata_bau[metadata_bau == "True"] <- T
metadata_bau[metadata_bau == "False"] <- F
metadata_bau$antibiotics_any[metadata_bau$antibiotics_any == 1] <- T
metadata_bau$antibiotics_any[metadata_bau$antibiotics_any == 0] <- F
metadata_bau$birthmode[metadata_bau$birthmode == "Vaginal"] <- "v"
metadata_bau$birthmode[metadata_bau$birthmode == "Cesarean"] <- "c"
metadata_bau$sex[metadata_bau$sex == "female"] <- "f"
metadata_bau$sex[metadata_bau$sex == "male"] <- "m"
metadata_bau <- metadata_bau %>% mutate(
  antibiotics_any = as.logical(antibiotics_any),
  antibiotics_one_week_before = as.logical(antibiotics_one_week_before)
)

# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_bau$Bases, n_lines$read_counts[match(metadata_bau$run_accession, n_lines$sample_ID)])  # 


write.csv(metadata_bau, "metadata_baumann_dudenhoeffer_2018_healthy.csv")

# nothing to sort out

# move good files
final_samples_bau <- metadata_bau$run_accession
final_files_bau <- c(paste(final_samples_bau, "_1.fastq.gz", sep = ""),
                       paste(final_samples_bau, "_2.fastq.gz", sep = ""),
                       paste(final_samples_bau, ".fastq.gz", sep = ""))
final_files_bau <- final_files_bau[file.exists(final_files_bau)]
file.move(final_files_bau, "fastq_files/")

# clean up
rm(ena_file_baumann, metadata_bau, final_files_bau, final_samples_bau, multiple_samples, n_lines, sra_file)


# Shao	2019 ###################################################################

setwd(paste(maindir, "shotgun/shao_2019/", sep = ""))
ena_file_shao <- read.table("filereport_read_run_PRJEB32631_tsv.txt", sep = "\t", header = T, comment.char = "")
ena_file_shao <- ena_file_shao[c("study_accession", "run_accession",
                                   "read_count", "sample_alias", "sample_title", "secondary_sample_accession")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file[,c("Instrument","Run")]
all(ena_file_shao$run_accession %in% sra_file$Run)
ena_file_shao <- merge(ena_file_shao, sra_file,by.x = "run_accession", by.y = "Run")


metadata_shao <- read.table("metadata.txt", sep = "\t", header = T)

# no duplicated samples
sum(duplicated(ena_file_shao$secondary_sample_accession))

#combine data.frames
sum(!ena_file_shao$secondary_sample_accession %in% metadata_shao$Accession)
nrow(ena_file_shao)
metadata_sha_combined <- merge(metadata_shao, ena_file_shao, by.x = "Accession", by.y = "secondary_sample_accession")
nrow(metadata_sha_combined)

# remove mother samples
sort_out_samples_shao <- metadata_sha_combined$run_accession[metadata_sha_combined$Time_point == "Mother"]
metadata_sha_combined <- metadata_sha_combined[metadata_sha_combined$Time_point != "Mother",]
# remove NAs in age column
sort_out_nas <- metadata_sha_combined$Time_point == "Infancy" &is.na(metadata_sha_combined$Infancy_sampling_age_months)
sort_out_samples_shao <- c(sort_out_samples_shao,
                           metadata_sha_combined$run_accession[sort_out_nas])
metadata_sha_combined <- metadata_sha_combined[!sort_out_nas,]

# unify colnames
metadata_sha_combined %<>% mutate(subject_ID = Individual,
                                  birthmode = Delivery_mode,
                                  sex = Gender,
                                  birth_weight = Birth_weight / 1000,
                                  food = Feeding_method,
                                  sample_ID = sample_title,
                                  geographic_location_.latitude. = 51.509865,
                                  geographic_location_.longitude. = -0.118092,
                                  .keep = "unused")


metadata_sha_combined$family_ID <- metadata_sha_combined$subject_ID
metadata_sha_combined$country <- "UK"
metadata_sha_combined$region <- "UK"
metadata_sha_combined$lifestyle <- "industrialized"
metadata_sha_combined$study <- "shao_2019"

# rename entries
metadata_sha_combined[metadata_sha_combined == "Yes"] <- T
metadata_sha_combined[metadata_sha_combined == "No"] <- F
metadata_sha_combined[metadata_sha_combined == "no"] <- F
metadata_sha_combined$sex[metadata_sha_combined$sex == "Female"] <- "f"
metadata_sha_combined$sex[metadata_sha_combined$sex == "Male"] <- "m"
metadata_sha_combined$birthmode[metadata_sha_combined$birthmode == "Caesarean"] <- "c"
metadata_sha_combined$birthmode[metadata_sha_combined$birthmode == "Vaginal"] <- "v"

# age column
metadata_sha_combined$age <- NA
metadata_sha_combined$age[metadata_sha_combined$Time_point == "Infancy"] <- 
   metadata_sha_combined$Infancy_sampling_age_months[metadata_sha_combined$Time_point == "Infancy"] %>% 
   as.numeric() * 30
metadata_sha_combined$age[metadata_sha_combined$Infancy_sampling_age_months == "Neonatal"] <- 
   metadata_sha_combined$Time_point[metadata_sha_combined$Infancy_sampling_age_months == "Neonatal"] %>%
   as.numeric()
metadata_sha_combined[,c("Time_point", "Infancy_sampling_age_months")] <- NULL

# food column, assume that children consume solid food from infancy period on (90 days)
metadata_sha_combined$food[metadata_sha_combined$food == "BF" & metadata_sha_combined$age < 90] <- "breast"
metadata_sha_combined$food[metadata_sha_combined$food == "BF" & metadata_sha_combined$age >= 90] <- "mixed_b_s"
metadata_sha_combined$food[metadata_sha_combined$food == "NoBF" & metadata_sha_combined$age < 90] <- "formula"
metadata_sha_combined$food[metadata_sha_combined$food == "NoBF" & metadata_sha_combined$age >= 90] <- "mixed_f_s"
metadata_sha_combined$food[metadata_sha_combined$food == "Mixed" & metadata_sha_combined$age < 90] <- "mixed_b_f"
metadata_sha_combined$food[metadata_sha_combined$food == "Mixed" & metadata_sha_combined$age >= 90] <- "mixed_b_f_s"

# antibiotics
metadata_sha_combined$antibiotics_any <- NA
metadata_sha_combined$antibiotics_any[metadata_sha_combined$Abx_Baby_in_hospital == "TRUE" | 
                                         metadata_sha_combined$Abx_Baby_after_hospital == "TRUE"] <- T
metadata_sha_combined$antibiotics_any[metadata_sha_combined$Abx_Baby_in_hospital == "FALSE" & 
                                         metadata_sha_combined$Abx_Baby_after_hospital == "FALSE"] <- F


# remove individuals with only one sample
multiple_samples <- unique(metadata_sha_combined$subject_ID[duplicated(metadata_sha_combined$subject_ID)])
sort_out_samples_shao <- c(sort_out_samples_shao, metadata_sha_combined$run_accession[!metadata_sha_combined$subject_ID %in% multiple_samples])
metadata_sha_combined <- metadata_sha_combined[metadata_sha_combined$subject_ID %in% multiple_samples,]

# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_sha_combined$read_count, n_lines$read_counts[match(metadata_sha_combined$run_accession, n_lines$sample_ID)])  # 


write.csv(metadata_sha_combined, "metadata_shao_2019_healthy.csv")

sort_out_samples_shao <- unique(sort_out_samples_shao)
sort_out_files_shao <- c(paste(sort_out_samples_shao, "_1.fastq.gz", sep = ""),
                        paste(sort_out_samples_shao, "_2.fastq.gz", sep = ""),
                        paste(sort_out_samples_shao, ".fastq.gz", sep = ""))
sort_out_files_shao <- sort_out_files_shao[file.exists(sort_out_files_shao)]
file.move(sort_out_files_shao, "sorted_out/")


# move good files
final_samples_shao <- metadata_sha_combined$run_accession
final_files_shao <- c(paste(final_samples_shao, "_1.fastq.gz", sep = ""),
                     paste(final_samples_shao, "_2.fastq.gz", sep = ""),
                     paste(final_samples_shao, ".fastq.gz", sep = ""))
final_files_shao <- final_files_shao[file.exists(final_files_shao)]
file.move(final_files_shao, "fastq_files/")

rm(ena_file_shao, metadata_sha_combined, metadata_shao, final_files_shao, final_samples_shao,
   multiple_samples, sort_out_samples_shao, sort_out_files_shao, sort_out_nas, n_lines, sra_file)




# Coker	2021 ###################################################################
# waiting for response

# Murphy	2019 #################################################################
# only SraRunTable.txt available
setwd(paste(maindir, "shotgun/murphy_2019/", sep = ""))
ena_file_murphy <- read.table("filereport_read_run_PRJNA345144_tsv.txt", sep = "\t", header = T)
ena_file_murphy <- ena_file_murphy[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]

metadata_murphy <- read.table("SraRunTable.txt", sep = ",", header = T)
# remove columns
metadata_murphy[,c("Assay.Type", "BioSampleModel", "Bytes", "Center.Name", "Collection_Date",
                  "complete_set_faecal_samples", "Consent", "DATASTORE.filetype",
                  "DATASTORE.provider", "DATASTORE.region",
                  "env_biome", "env_feature", "env_material", "eth",
                  "faecal_sample_taken_at_12_months", "faecal_sample_taken_at_24_mths",
                  "faecal_sample_taken_at_3_months", "faecal_sample_taken_at_4_yrs",
                  "faecal_sample_taken_at_6_yrs", "faecal_sample_taken_at_birth",
                  "geo_loc_name_country", "geo_loc_name_country_continent",
                  "geo_loc_name", "host_body_product", "Host", "Lat_Lon", 
                  "LibraryLayout", "LibrarySelection", "LibrarySource",
                  "Organism", "Platform", "ReleaseDate")] <- NULL

# no duplicates
sum(duplicated(metadata_murphy$Sample.Name))

# remove intervention samples
sort_out_samples_murphy <- metadata_murphy$Run[metadata_murphy$studygroup != "placeb"]
metadata_murphy <- metadata_murphy[metadata_murphy$studygroup == "placeb",]
# individuals with only one sample
multiple_samples <- unique(metadata_murphy$Patient_ID[duplicated(metadata_murphy$Patient_ID)])
sort_out_samples_murphy <- c(sort_out_samples_murphy,
                             metadata_murphy$Run[!metadata_murphy$Patient_ID %in% multiple_samples])
metadata_murphy <- metadata_murphy[metadata_murphy$Patient_ID %in% multiple_samples,]

# age column:
metadata_murphy$age <- NA
metadata_murphy$age[metadata_murphy$Host_Age == "At_Birth"] <- 3
metadata_murphy$age[metadata_murphy$Host_Age == "3_months"] <- 90
metadata_murphy$age[metadata_murphy$Host_Age == "12_months"] <- 365
metadata_murphy$age[metadata_murphy$Host_Age == "24_months"] <- 730
metadata_murphy$Host_Age <- NULL

# antibiotics column
metadata_murphy$antibiotics_before <- F
metadata_murphy$antibiotics_before[metadata_murphy$age >= 90 & 
                                      metadata_murphy$Antibiotics_before_3_months == 1] <- T
metadata_murphy$antibiotics_before[metadata_murphy$age >= 180 & 
                                      metadata_murphy$Antibiotics_before_6_months == 1] <- T
metadata_murphy$Antibiotics_before_3_months <- NULL
metadata_murphy$Antibiotics_before_6_months <- NULL

metadata_murphy$caesar[metadata_murphy$caesar == 1] <- "c"
metadata_murphy$caesar[metadata_murphy$caesar == 0] <- "v"

# unify colnames
colnames(metadata_murphy)[match(c("Run", "caesar", "host_sex", "Patient_ID", "Sample.Name"),
                               colnames(metadata_murphy))] <- c("run_accession", "birthmode",
                                                                "sex", "subject_ID",
                                                                "sample_ID")



metadata_murphy$country <- "NEW_ZEALAND"
metadata_murphy$region <- "NEW_ZEALAND"
metadata_murphy$study <- "murphy_2019"
metadata_murphy$lifestyle <- "industrialized"
metadata_murphy$family_ID <- metadata_murphy$subject_ID
metadata_murphy$sex[metadata_murphy$sex == "female"] <- "f"
metadata_murphy$sex[metadata_murphy$sex == "male"] <- "m"
sort(colnames(metadata_murphy))
metadata_murphy$geographic_location_.latitude. <- -45.8593 
metadata_murphy$geographic_location_.longitude. <- 170.5083


# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_murphy$Bases, n_lines$read_counts[match(metadata_murphy$run_accession, n_lines$sample_ID)])  # 


write.csv(metadata_murphy, "metadata_murphy_2019_healthy.csv")

# move sorted out samples
sort_out_samples_murphy <- unique(sort_out_samples_murphy)
sort_out_files_murphy <- c(paste(sort_out_samples_murphy, "_1.fastq.gz", sep = ""),
                          paste(sort_out_samples_murphy, "_2.fastq.gz", sep = ""),
                          paste(sort_out_samples_murphy, ".fastq.gz", sep = ""))
sort_out_files_murphy <- sort_out_files_murphy[file.exists(sort_out_files_murphy)]
file.move(sort_out_files_murphy, "sorted_out/")

# move good files
final_samples_murphy <- metadata_murphy$run_accession
final_files_murphy <- c(paste(final_samples_murphy, "_1.fastq.gz", sep = ""),
                       paste(final_samples_murphy, "_2.fastq.gz", sep = ""),
                       paste(final_samples_murphy, ".fastq.gz", sep = ""))
final_files_murphy <- final_files_murphy[file.exists(final_files_murphy)]
file.move(final_files_murphy, "fastq_files/")

# clean up
rm(ena_file_murphy, metadata_murphy, final_files_murphy, final_samples_murphy,
   multiple_samples, sort_out_files_murphy, sort_out_samples_murphy, n_lines)




# Tapiainen	2019 ###############################################################
# waiting for response


# Gehrig	2019 #################################################################

setwd(paste(maindir, "shotgun/gehrig_2019/", sep = ""))
ena_file_gehr <- read.table("filereport_read_run_PRJEB26419_tsv_human.txt", sep = "\t", header = T)
ena_file_gehr <- ena_file_gehr[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
ena_file_gehr <- separate(ena_file_gehr, col = "sample_alias", into = c("sample_ID", "meta"), sep = ".meta", remove = T)

metadata_gehr_age <- read.table("metadata_age_shotgun.txt", sep = "\t", header = T)
metadata_gehr_sra <- read.table("SraRunTable.txt", sep = ",", header = T)
# remove some cols:
metadata_gehr_sra[,c("Assay.Type", "AvgSpotLen", "Bases", "BioProject", "BioSample",
                     "Bytes", "Center.Name", "Consent", "DATASTORE.filetype",
                     "DATASTORE.provider", "DATASTORE.region", "ENA.FIRST.PUBLIC..run.",
                     "ENA.FIRST.PUBLIC", "ENA.LAST.UPDATE..run.", "ENA.LAST.UPDATE",
                     "INSDC_center_alias", "INSDC_center_name", "INSDC_first_public",
                     "INSDC_last_update", "INSDC_status", "LibrarySelection",
                     "Platform", "ReleaseDate", "LibrarySource", "LibraryLayout")] <- NULL

metadata_gehr_combined <- merge(ena_file_gehr, metadata_gehr_age, by.x = "sample_ID", by.y = "Sample.ID", all.x = T)
metadata_gehr_combined <- merge(metadata_gehr_combined, metadata_gehr_sra, by.x = "run_accession", by.y = "Run", all = T)

unique(metadata_gehr_combined$library_name)
# remove non-human and SAM/MAM samples:
sort_out_samples <- metadata_gehr_combined$run_accession[metadata_gehr_combined$common_name %in% c("house mouse", "pig")]
metadata_gehr_combined <- metadata_gehr_combined[!metadata_gehr_combined$common_name %in% c("house mouse", "pig"),]
sort_out_samples <- c(sort_out_samples,
                      metadata_gehr_combined$run_accession[metadata_gehr_combined$Organism %in% c("bacterium", "mouse gut metagenome")])
metadata_gehr_combined <- metadata_gehr_combined[!metadata_gehr_combined$Organism %in% c("bacterium", "mouse gut metagenome"),]
sort_out_samples <- c(sort_out_samples,
                      metadata_gehr_combined$run_accession[metadata_gehr_combined$Library.Name %in% c("table S23, MDCF trial",
                                                                                                      "table S5A, SAM trial",
                                                                                                      "table S16C, mouse sample",
                                                                                                      "table S15B, mouse fecal sample")])
metadata_gehr_combined <- metadata_gehr_combined[!metadata_gehr_combined$Library.Name %in% c("table S23, MDCF trial",
                                                                                             "table S5A, SAM trial",
                                                                                             "table S16C, mouse sample",
                                                                                             "table S15B, mouse fecal sample"),]
# differentiate between shotgun and 16S:
sort_out_samples <- c(sort_out_samples,
                      metadata_gehr_combined$run_accession[is.na(metadata_gehr_combined$meta)])
metadata_gehr_combined <- metadata_gehr_combined[!is.na(metadata_gehr_combined$meta),]

# keep only healthy children
sort_out_samples <- c(sort_out_samples,
                      metadata_gehr_combined$run_accession[is.na(metadata_gehr_combined$Age..months.)])
metadata_gehr_combined <- metadata_gehr_combined[!is.na(metadata_gehr_combined$Age..months.),]

# add antibiotics information from subramanian
abx_subr_data <- read.table("suppl_2_subr.txt", sep = "\t", header = T) %>%
  select(Fecal.Sample.ID, Breast.Milk, Formula1, Solid.Foods2, diarrhoea.at.the.time.of.sample.collection3,
         Antibiotics.within.7.days.prior.to.sample.collection)
sum(metadata_gehr_combined$sample_ID %in% abx_subr_data$Fecal.Sample.ID)
metadata_gehr_combined <- left_join(metadata_gehr_combined, abx_subr_data, by = join_by(sample_ID == Fecal.Sample.ID))

# no diarrhoea samples
metadata_gehr_combined$diarrhoea.at.the.time.of.sample.collection3 <- NULL

# unify colnames
metadata_gehr_combined %<>% mutate(antibiotics_one_week_before = ifelse(Antibiotics.within.7.days.prior.to.sample.collection == "Yes", yes = T, no = F),
                                   food = case_when(Breast.Milk == "Yes" & Formula1 == "No" & Solid.Foods2 == "No" ~ "breast",
                                                    Breast.Milk == "Yes" & Formula1 == "Yes" & Solid.Foods2 == "No" ~ "mixed_b_f",
                                                    Breast.Milk == "Yes" & Formula1 == "No" & Solid.Foods2 == "Yes" ~ "mixed_b_s",
                                                    Breast.Milk == "Yes" & Formula1 == "Yes" & Solid.Foods2 == "Yes" ~ "mixed_b_f_s",
                                                    Breast.Milk == "No" & Formula1 == "Yes" & Solid.Foods2 == "No" ~ "formula",
                                                    Breast.Milk == "No" & Formula1 == "Yes" & Solid.Foods2 == "Yes" ~ "mixed_f_s",
                                                    Breast.Milk == "No" & Formula1 == "No" & Solid.Foods2 == "Yes" ~ "solid",
                                                    .default = NA),
                                   .keep = "unused")


# no individuals with only one sample
metadata_gehr_combined <- separate(metadata_gehr_combined, "Submitter_Id",
                                   into = c("subject_ID", NA), sep = ".m", remove = T)
multiple_samples <- unique(metadata_gehr_combined$subject_ID[duplicated(metadata_gehr_combined$subject_ID)])
sum(!metadata_gehr_combined$subject_ID %in% multiple_samples)


# clean up cols
metadata_gehr_combined[,c("meta", "sample_title", "Index.1", "Index.2",
                          "Library.Name", "Organism", "common_name",
                          "library_name")] <- NULL

# no duplicated samples
sum(duplicated(metadata_gehr_combined$sample_alias))

# unify colnames
colnames(metadata_gehr_combined)[match(c("Age..months."),
                                       names(metadata_gehr_combined))] <- c("age")
metadata_gehr_combined$age <- metadata_gehr_combined$age * 30

metadata_gehr_combined$country <- "BANGLADESH"
metadata_gehr_combined$region <- "DHAKA"
metadata_gehr_combined$study <- "gehrig_2019"
metadata_gehr_combined$lifestyle <- "non_industrialized"
metadata_gehr_combined$family_ID <- metadata_gehr_combined$subject_ID
metadata_gehr_combined$geographic_location_.latitude. = 23.933
metadata_gehr_combined$geographic_location_.longitude. = 88.983


# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_gehr_combined$read_count, n_lines$read_counts[match(metadata_gehr_combined$run_accession, n_lines$sample_ID)])  # 

metadata_gehr_combined$run_accession[metadata_gehr_combined$read_count != n_lines$read_counts[match(metadata_gehr_combined$run_accession, n_lines$sample_ID)]]

write.csv(metadata_gehr_combined, "metadata_gehrig_2019_healthy.csv")

# # move sorted out samples
sort_out_samples <- unique(sort_out_samples)
sort_out_files <- c(paste(sort_out_samples, "_1.fastq.gz", sep = ""),
                    paste(sort_out_samples, "_2.fastq.gz", sep = ""),
                    paste(sort_out_samples, ".fastq.gz", sep = ""))
sort_out_files <- sort_out_files[file.exists(sort_out_files)]
file.move(sort_out_files, "sorted_out/")

# move good files from shotgun folder
final_samples_gehr <- metadata_gehr_combined$run_accession
final_files_gehr <- c(paste0(final_samples_gehr, "_1.fastq.gz"),
                      paste0(final_samples_gehr, "_2.fastq.gz"),
                      paste0(final_samples_gehr, ".fastq.gz"))
final_files_gehr <- final_files_gehr[file.exists(final_files_gehr)]
file.move(final_files_gehr, "fastq_files/")


# clean up
rm(ena_file_gehr, metadata_gehr_age, metadata_gehr_combined, metadata_gehr_sra,
   final_files_gehr, final_samples_gehr, multiple_samples, n_lines, sort_out_files, sort_out_samples,
   abx_subr_data)




# Olm	2022 #####################################################################
setwd(paste(maindir, "shotgun/olm_2022/", sep = ""))
ena_file_olm <- read.table("filereport_read_run_PRJEB49206_tsv.txt", sep = "\t", header = T)
ena_file_olm <- ena_file_olm[,c("study_accession", "run_accession", "library_name",
                               "read_count", "sample_alias", "sample_title")]


sra_file_olm <- read.table("SraRunTable.txt", header = T, sep = ",")
sra_file_olm[,c("Assay.Type", "Center.Name", "Consent", "ENA.FIRST.PUBLIC..run.",
                "ENA.LAST.UPDATE..run.", "geo_loc_name_country",
                "geo_loc_name_country_continent", "LibraryLayout",
                "LibrarySelection", "LibrarySource", "Platform", "ReleaseDate",
                "AvgSpotLen", "Bases", "Bytes", "DATASTORE.filetype", "DATASTORE.provider",
                "DATASTORE.region", "ENA.FIRST.PUBLIC", "ENA.LAST.UPDATE", "INSDC_center_name",
                "INSDC_first_public", "INSDC_last_update", "INSDC_status", "dna_extracted",
                "physical_specimen_location", "physical_specimen_remaining", "Elevation",
                "env_feature", "env_material", "gender", "host_common_name", "host_scientific_name",
                "host_taxid", "host_weight",  "age_units",
                "height_units", "maize_observation", "sample_source", "sample_type_expanded",
                "sample_type_field_notes", "seasonality", "seasonality_unique_subjects",
                "weight_units", "Condition", "COUNTRY", "host_age_units", "host_height_units",
                "host_weight_units", "alco", "brush_freq", "brush", "drk", "exer_freq", "exer",
                "ferm_freq", "ferm", "Fish", "food", "fuel", "gra", "hle", "hlth_tim", "kit_loc",
                "lit", "Meat", "med", "milk", "mlkt", "noc", "scar", "sick2",
                "sick_loc2", "sisnu", "smok", "soda", "toba", "toi", "trad_med",
                "ygrt_freq", "ygrt", "Collection_Date", "environment_.biome.",
                "environment_.feature.", "geographic_location_.country_and.or_sea.",
                "geographic_location_.latitude.", "geographic_location_.longitude.",
                "human_gut_environmental_package", "Investigation_type", "project_name",
                "Sequencing_method", "annotation_source", "extrachromosomal_elements",
                "isolation_and_growth_condition", "number_of_replicons", "ploidy",
                "reference_for_biomaterial", "sub_species", "BioSampleModel",
                "env_broad_scale", "env_local_scale", "env_medium", "Host_diet",
                "Host", "isol_growth_condt", "Isolate", "Lat_Lon", "misc_param",
                "host_body_mass_index", "host_height", "parasite", "Age",
                "environment_.material.", "Popoulation")] <- NULL
metadata_olm_1 <- read.table("metadata_1.txt", header = T, sep = "\t")
metadata_olm_1 <- metadata_olm_1[!duplicated(metadata_olm_1$run_acc),]
metadata_olm_1 <- metadata_olm_1[metadata_olm_1$age <= 3,]
metadata_olm_1[,c("location", "Mom_Infant", "unique_bases", "age", "rarefied_pairs")] <- NULL
metadata_olm_2 <- read.table("metadata_2.txt", header = T, sep = "\t")
metadata_olm_2[,c("Study", "country", "Gbp", "lifestyle", "breastfeeding_status")] <- NULL
metadata_olm_3 <- read.table("metadata_3.txt", header = T, sep = "\t")
metadata_olm_3[,c("Study", "lifestyle", "study_accession")] <- NULL
metadata_olm_qiita <- read.table("11358_20221025-080819.txt", header = T, sep = "\t")
# clean up qiita record before merging
metadata_olm_qiita[c("dna_extracted", "elevation", "empo_1", "empo_2", "geo_loc_name",
                     "host_age_units", "host_body_mass_index","host_body_product",
                     "host_body_site", "host_common_name", "host_height", "host_height_units",
                     "host_life_stage", "host_scientific_name", "host_taxid", "host_weight",
                     "host_weight_units", "latitude", "longitude", "maize_observation",
                     "physical_specimen_location", "physical_specimen_remaining",
                     "qiita_empo_2", "qiita_empo_3", "qiita_study_id", "sample_type_expanded",
                     "sample_type_field_notes", "scientific_name", "season", "seasonality",
                     "seasonality_unique_subjects", "story", "taxon_id", "title",
                     "water_source", "qiita_empo_1")] <- NULL
metadata_olm_qiita <- metadata_olm_qiita[metadata_olm_qiita$env_biome == "anthropogenic terrestrial biome",]
metadata_olm_qiita <- metadata_olm_qiita[metadata_olm_qiita$env_material == "feces",]
metadata_olm_qiita <- metadata_olm_qiita[metadata_olm_qiita$sample_source == "HUMAN",]
metadata_olm_qiita$host_age[metadata_olm_qiita$host_age == "not collected"] <- NA
metadata_olm_qiita$host_age <- as.numeric(metadata_olm_qiita$host_age)
metadata_olm_qiita <- metadata_olm_qiita[metadata_olm_qiita$host_age <= 3 | is.na(metadata_olm_qiita$host_age),]
metadata_olm_qiita[, c("country", "empo_3", "env_biome", "env_feature", "env_material",
                       "env_package", "host_body_habitat", "qiita_empo_1", "sample_source",
                       "sample_type")] <- NULL
metadata_olm_qiita$description <- gsub("feces sample ", "", metadata_olm_qiita$description)


metadata_olm_combined <- merge(ena_file_olm, sra_file_olm, by.x = "run_accession", by.y = "Run")
# remove non-human and adult samples:
metadata_olm_combined <- metadata_olm_combined %>%
  filter(Organism == "human gut metagenome",
         host_life_stage != "adult",
         Host_age != "adults over eighteen") %>%
  separate(., col = "sample_title", into = c(NA, "qiita_id"), sep = "[.]", remove = F) %>%
  left_join(., metadata_olm_qiita, by = c("qiita_id" = "description")) %>%
  left_join(., metadata_olm_1, by = c("External_Id" = "BioSample")) %>%
  left_join(., metadata_olm_2, by = c("run_accession" = "run_acc")) %>%
  left_join(., metadata_olm_3 %>% 
              mutate(SampleName = as.character(SampleName)), by = c("qiita_id" = "SampleName")) %>%
  # merge columns
  mutate(SRA = ifelse(is.na(SRA), yes = ifelse(is.na(SRA.x), yes = SRA.y, no = SRA.x), no = SRA),
         age_months = ifelse(is.na(age_months), yes = ifelse(is.na(age_months.x), yes = age_months.y, no = age_months.x), no = age_months),
         age_months = ifelse(is.na(age_months), yes = as.numeric(Host_age) * 12, no = age_months),
         bushcamp = ifelse(is.na(bushcamp), yes = ifelse(is.na(bush_camp.x), yes = bush_camp.y, no = bush_camp.x), no = bushcamp),
         subject = ifelse(is.na(subject), yes = ifelse(is.na(host_subject_id.x), yes = host_subject_id.y, no = host_subject_id.x), no = subject),
         sample_name = ifelse(is.na(sample_name), yes = NCBI_SampleName, no = sample_name),
         run_accession = ifelse(is.na(run_accession), yes = run_acc, no = run_accession),
         study_accession_1 = "PRJEB49206",
         study_accession_2 = "PRJEB27517",
         .keep = "unused") %>%
  filter(!is.na(age_months), age_months <= 24, !is.na(subject), !is.na(run_accession)) %>%
  select(-c("Organism", "empo_1", "empo_2", "empo_3", "host_body_habitat",
            "host_body_site", "host_life_stage", "env_biome", "Env_package",
            "SAMPLE_TYPE", "location", "host_body_product", "SRA.Study", 
            "INSDC_center_alias", "assembly_acc.y", "BioSample.y", "COLLECTION_TIMESTAMP", 
            "experimentId.x", "experimentId.y", "host_age", "read_count", "Sample_Name",
            "sex.x", "gender", "study_accession.x", "study_accession.y", "unique_reads", 
            "BioProject", "BioSample.x", "External_Id", "Sample.Name", "DYAD", 
            "sample_alias", "sample_title", "FecalSample_ID", "Submitter_Id",
            "run_acc", "Library.Name", "SampleName.x", "SampleName.y", "Population",
            "geo_loc_name"))



# unify colnames
metadata_olm_combined <- metadata_olm_combined %>%
  mutate(sample_ID = library_name,
         sex = sex.y,
         subject_ID = subject,
         family_ID = family,
         age = age_months * 30,
         geographic_location_.latitude. = LATITUDE,
         geographic_location_.longitude. = LONGITUDE,
         family_ID = ifelse(is.na(family_ID), yes = subject_ID, no = family_ID),
         country = "TANZANIA",
         region = "HADZA",
         lifestyle = "non_industrialized",
         study = "olm_2022",
         birthmode = "v",
         sex = ifelse(sex == "female", yes = "f", no = "m")
         )

# remove single sample individuals
# multiple_samples <- unique(metadata_olm_combined$subject_ID[duplicated(metadata_olm_combined$subject_ID)])
# metadata_olm_combined <- metadata_olm_combined[metadata_olm_combined$subject_ID %in% multiple_samples,]

metadata_olm_combined$run_accession %>%
  paste("fastq-dump --split-e --gzip ",.) %>%
  data.frame(`#!/bin/bash` = .) %>%
  write.table(., file = "/fast/AG_Forslund/rob/studies/shotgun/olm_2022/pull_single.sh",
              sep = " ", quote = F, row.names = F, col.names = "#!/bin/bash")


write.csv(metadata_olm_combined, "metadata_olm_2022_healthy.csv")


# move good files from shotgun folder
final_samples_olm <- metadata_olm_combined$run_accession
final_files_olm <- c(paste0(final_samples_olm, "_1.fastq.gz"),
                      paste0(final_samples_olm, "_2.fastq.gz"),
                      paste0(final_samples_olm, ".fastq.gz"))
final_files_olm <- final_files_olm[file.exists(final_files_olm)]
file.move(final_files_olm, "fastq_files/")


# clean up
rm(ena_file_olm, metadata_olm_1, metadata_olm_2, metadata_olm_3, metadata_olm_combined,
   metadata_olm_qiita, sra_file_olm, final_samples_olm, final_files_olm, n_lines,
   multiple_samples)




# Smith	2013 ###################################################################
# setwd(paste(maindir, "shotgun/smith_2013", sep = ""))
# ena_file_smith <- read.table("filereport_read_run_PRJEB3350_tsv.txt", sep = "\t", header = T, comment.char = "")
# ena_file_smith <- ena_file_smith %>% 
#    select(c(study_accession, run_accession, library_name, read_count, sample_alias, sample_title)) %>%
#    separate(sample_title, c("subject_ID", "time_point"), sep = "[.]", remove = F, convert = T) %>%
#    mutate(subject_ID = str_sub(subject_ID, start = 2), time_point = as.numeric(gsub("\\D", "", time_point))) %>%
#    mutate(sample_ID = paste(subject_ID, time_point, sep = "."))
# 
# # duplicated samples?
# any(duplicated(ena_file_smith$sample_alias))
# 
# # metadata_subj_smith <- read.table("metadata_subj.txt", sep = "\t", header = T)
# # metadata_sample_smith <- read.table("metadata_sample.txt", sep = "\t", header = T)
# # metadata_sample_smith <- metadata_sample_smith %>% dplyr::mutate(MUAC = Mid.Upper.Arm.Circumference..MUAC. , .keep = "unused")
# metadata_sample_smith_2 <- read.table("metadata_smith.txt", sep = "\t", header = T)
# metadata_sample_smith_2 <- metadata_sample_smith_2 %>%
#   dplyr::mutate(Sample.ID = paste(str_remove(Person.ID, "^0+"), gsub("h.*[.]", "", Sample.ID), sep = "."),
#          Child.ID = Person.ID,
#          .keep = "unused") %>%
#   filter(Health.category == "Healthy")
# # metadata_sample_combined <- bind_rows(metadata_sample_smith, metadata_sample_smith_2)
# 
# sum(metadata_sample_smith_2$Sample.ID %in% ena_file_smith$sample_ID)
# metadata_smith_combined <- merge(metadata_sample_smith_2, ena_file_smith, by.x = "Sample.ID", by.y = "sample_ID") 
# 
# 
# metadata_smith_combined$Gender[metadata_smith_combined$Gender == "Female"] <- "f"
# metadata_smith_combined$Gender[metadata_smith_combined$Gender == "Male"] <- "m"
# metadata_smith_combined$Age..months <- metadata_smith_combined$Age..months * 30
# 
# # no single sample children
# multiple_samples <- unique(metadata_smith_combined$Child.ID[duplicated(metadata_smith_combined$Child.ID)])
# sum(!metadata_smith_combined$Child.ID %in% multiple_samples)
# 
# 
# # # add additional info from smith et al
# # metadata_blant_combined <- merge(metadata_blant_combined, location_metadata, by.x = "Child.ID", by.y = "Person.ID", all.x = T)
# # unify colnames
# metadata_smith_combined <- metadata_smith_combined %>% dplyr::mutate(samples_ID = Sample.ID, family_ID = Family.ID, sex = Gender,
#                                                               age = Age..months, region = Site, .keep = "unused")
# 
# metadata_smith_combined$country <- "MALAWI"
# metadata_smith_combined$study <- "smith_2013"
# metadata_smith_combined$lifestyle <- "non_industrialized"
# 
# # add height and weight data:
# metadata_smith_combined$height <- sitar::LMS2z(metadata_smith_combined$age/365, y = metadata_smith_combined$HAZ, sex = metadata_smith_combined$sex,
#                                       measure = "ht", ref = who06, toz = F)
# metadata_smith_combined$weight <- sitar::LMS2z(metadata_smith_combined$age/365, y = metadata_smith_combined$WAZ, sex = metadata_smith_combined$sex,
#                                       measure = "wt", ref = who06, toz = F) 
# 
# 
# sort_out_samples_smith <- ena_file_smith$run_accession[!ena_file_smith$run_accession %in% metadata_smith_combined$run_accession]
# # move sorted out samples
# sort_out_samples_smith <- unique(sort_out_samples_smith)
# sort_out_files_smith <- c(paste(sort_out_samples_smith, "_1.fastq.gz", sep = ""),
#                           paste(sort_out_samples_smith, "_2.fastq.gz", sep = ""),
#                           paste(sort_out_samples_smith, ".fastq.gz", sep = ""))
# sort_out_files_smith <- sort_out_files_smith[file.exists(sort_out_files_smith)]
# file.move(sort_out_files_smith, "sorted_out/")
# 
# 
# final_samples_smith <- metadata_smith_combined$run_accession
# final_files_smith <- c(paste0(final_samples_smith, "_1.fastq.gz"),
#                        paste0(final_samples_smith, "_2.fastq.gz"),
#                        paste0(final_samples_smith, ".fastq.gz"))
# final_files_smith <- final_files_smith[file.exists(final_files_smith)]
# file.move(final_files_smith, "fastq_files/")
# 
# write.csv(metadata_smith_combined, "metadata_smith_2013_healthy.csv")
# 
# rm(ena_file_smith, metadata_smith_combined, metadata_sample_smith_2, final_samples_smith,
#    final_files_smith, sort_out_samples_smith, sort_out_files_smith, multiple_samples)






# Vatanen 2022 #################################################################

setwd(paste(maindir, "shotgun/vatanen_2022", sep = ""))
ena_file_vat <- read.table("filereport_read_run_PRJNA806984_tsv.txt", sep = "\t", header = T)
ena_file_vat <- ena_file_vat[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file[,c("Instrument","Run")]
all(ena_file_vat$run_accession %in% sra_file$Run)
ena_file_vat <- merge(ena_file_vat, sra_file,by.x = "run_accession", by.y = "Run")



metadata_vat_combined <- ena_file_vat %>% separate(sample_alias, c("subject_ID", "age"), sep = "[_]", remove = F) %>%
  mutate(country = "BANGLADESH",
         study = "vatanen_2022",
         lifestyle = "non_industrialized",
         region = "DHAKA",
         sample_ID = sample_alias,
         family_ID = subject_ID,
         geographic_location_.latitude. = 23.933,
         geographic_location_.longitude. = 88.983,
         diarrhea = ifelse(grepl("DS", age), yes = "y", no = "n")) %>%
  filter(diarrhea == "n") %>%
  mutate(age = as.numeric(gsub("m", "", 
                    gsub("DS","",
                         gsub("birth", 0, age)))) * 30)
any(duplicated(metadata_vat_combined$sample_ID))

# remove diarrhea samples


# no individuals with only one sample
multiple_samples <- unique(metadata_vat_combined$subject_ID[duplicated(metadata_vat_combined$subject_ID)])
metadata_vat_combined <- metadata_vat_combined %>%
  filter(subject_ID %in% multiple_samples)
sum(!metadata_vat_combined$subject_ID %in% multiple_samples)


# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_vat_combined$read_count, n_lines$read_counts[match(metadata_vat_combined$run_accession, n_lines$sample_ID)])  # 

write.csv(metadata_vat_combined, "metadata_vatanen_2022_healthy.csv")

# # move sorted out samples
sort_out_samples <- ena_file_vat$run_accession[!ena_file_vat$run_accession %in% metadata_vat_combined$run_accession]
sort_out_files <- c(paste(sort_out_samples, "_1.fastq.gz", sep = ""),
                    paste(sort_out_samples, "_2.fastq.gz", sep = ""),
                    paste(sort_out_samples, ".fastq.gz", sep = ""))
sort_out_files <- sort_out_files[file.exists(sort_out_files)]
file.move(sort_out_files, "sorted_out/")

# move good files from shotgun folder
final_samples <- metadata_vat_combined$run_accession
final_files <- c(paste0(final_samples, "_1.fastq.gz"),
                      paste0(final_samples, "_2.fastq.gz"),
                      paste0(final_samples, ".fastq.gz"))
final_files <- final_files[file.exists(final_files)]
file.move(final_files, "fastq_files/")


# clean up
rm(ena_file_vat, metadata_vat_combined, final_files, final_samples,
   multiple_samples, sort_out_samples, sort_out_files, n_lines, sra_file)



# Robertson	2023 ###############################################################
setwd(paste(maindir, "shotgun/robertson_2023", sep = ""))
ena_file_robertson <- read.table("filereport_read_run_PRJEB51728_tsv.txt", sep = "\t", header = T, comment.char = "")
sra_file <- read.table("SraRunTable_2.txt", sep = ",", header = T, comment.char = "", colClasses = "character")
metadata_robertson_combined <- left_join(ena_file_robertson, sra_file, by = join_by(run_accession == Run))

metadata_robertson_combined[,c("tax_id", "scientific_name", "instrument_platform",
                               "library_name", "library_layout", "library_strategy",
                               "study_alias", "experiment_alias", "run_alias", "fastq_bytes",
                               "fastq_md5", "fastq_ftp", "submitted_md5", "submitted_ftp",
                               "sra_md5", "sra_ftp", "broker_name", "sample_title",
                               "Assay.Type", "BioProject", "BioSample",
                               "broad.scale_environmental_context", "Center.Name",
                               "Collection_Date", "Consent", "ENA.FIRST.PUBLIC..run.",
                               "ENA.LAST.UPDATE..run.", "environmental_medium",
                               "geo_loc_name_country_continent", "geographic_location_.country_and.or_sea.",
                               "host_body_product", "LibraryLayout", "LibrarySelection",
                               "LibrarySource", "local_environmental_context", "Organism",
                               "organism", "Platform", "project_name", "ReleaseDate",
                               "Scientific_Name", "SRA.Study", "Sample.Name")] <- NULL

metadata_robertson_combined <- metadata_robertson_combined %>%
  mutate(age = Host_Age,
         sex = case_when(host_sex == "female" ~ "f",
                         host_sex == "male" ~ "m"),
         country = toupper(geo_loc_name_country),
         region = toupper(geo_loc_name_country),
         study = "robertson_2023",
         sample_ID = sample_alias,
         lifestyle = "non_industrialized",
         .keep = "unused") %>%
  filter(!is.na(age))


metadata_1 <- readRDS("metaphlan.rds.gz")
metadata_1[,c("split", "name", "RelAb", "pct_human", "visitlabel", "pct_unknown",
              "evenness", "invsimpson", "richness", "shannon", "simpson", "antihelminth_2w",
              "anylatrine", "apostolic", "age", "sex", "benzylpenicillin", "birth_month",
              "carbonated_1d", "carbonated_7d", "carbonated_ever", "cheese_1d",
              "cheese_7d", "depression", "eatsoil", "ebf_init", "ebf24hr", "ebfever",
              "edyrs", "employed", "epds", "fecesobserved", "female", "gripe_1d",
              "gripe_7d", "gripe_ever", "handwash_water", "handwashing_b0", "hc",
              "hcz", "hh_mdds", "hhsize", "hwstation_items_b0", "improve_latrine",
              "improvedfloor", "iycf", "laz", "livestockinhse", "loosestool",
              "marital_status", "mom_age", "mom_cd4", "mom_edu", "mom_hb",
              "mom_height", "neonatal_death_elig", "opendefecation", "otherliq_1d",
              "otherliq_7d", "otherliq_ever", "othersoft_1d", "othersoft_7d", "othersoft_ever",
              "ownschickens", "parity", "pm12", "pmtct_artuse", "pmtct_contri",
              "prop_opendefecation", "religion", "sam_dnk", "sam_none", "sam_nosuckle",
              "sam_noweightgain", "sam_oedema", "sam_referral", "sam_wasting", "sample_month",
              "sample_season", "timetowater", "treatwater", "vel_laz", "vel_whz",
              "wash", "water_1d", "water_7d", "water_ever", "waz", "wealth_index_score",
              "whz", "age_group", "log_age", "Visit", "SubjectID_mom", "csi", "erythromycin",
              "ironfortified_1d", "ironfortified_7d", "ironfortified_ever", "mom_muac",
              "muac", "muacz", "oralmulti_1d", "oralmulti_7d", "oralmulti_ever",
              "perinealskininfection", "antibiotic_dnk", "cheese_ever", "drink_breastfeed",
              "formula_1d", "fruitveg_1d", "fruitveg_7d", "fruitveg_ever", "pwdrmilk_ever",
              "redmeat_7d", "redmeat_ever", "sugarwater_ever")] <- NULL
metadata_1 <- distinct(metadata_1) %>%
  select_if(~ !all(is.na(.)))
metadata_robertson_combined <- left_join(metadata_robertson_combined, metadata_1,
                                         by = join_by(sample_ID == Sample)) %>%
  filter(diarrhea_bloodstool == "no" | is.na(diarrhea_bloodstool),
         inf_m18_hivstatus == "HUU" | inf_m18_hivstatus == "(Missing)" | is.na(inf_m18_hivstatus),
         infant_hiv_status == "HUU" | is.na(infant_hiv_status),
         preg_hiv_status == "Negative",
         age != "",
         healthy_cohort,
         cc_outcome_elig != "AGA/preterm",
         !(diarrhea_1d == "yes" & age < 7), # remove if diarrhea one week before
         !(diarrhea_2w == "yes" & age < 14 & age > 7))
metadata_robertson_combined[,c("diarrhea_bloodstool", "inf_m18_hivstatus", "infant_hiv_status",
                               "preg_hiv_status", "healthy_cohort", "cc_outcome_elig",
                               "diarrhea_1d", "diarrhea_2w")] <- NULL
metadata_robertson_combined[,c("Bytes", "DATASTORE.filetype", "DATASTORE.provider", "DATASTORE.region", "ENA_first_public",
                               "ENA.LAST.UPDATE", "INSDC_center_alias", "INSDC_center_name", "INSDC_first_public",
                               "INSDC_last_update", "INSDC_status", "create_date", "version", "Submitter_Id",
                               "beans_1d", "beans_7d", "beans_ever", "cookoil_7d", "cookoil_ever", "cookoils_1d",
                               "dairy_1d", "dairy_7d", "dairy_ever", "darkgreenveg_1d", "darkgreenveg_7d", "darkgreenveg_ever",
                               "dessert_1d", "dessert_7d", "dessert_ever", "eggfishpoultry_1d", "eggfishpoultry_7d",
                               "eggfishpoultry_ever", "eggs_1d", "eggs_7d", "eggs_ever", "fish_1d", "fish_7d", "fish_ever",
                               "fleshfoods_1d", "fleshfoods_7d", "fleshfoods_ever", "freshmilk_1d", "freshmilk_7d",
                               "freshmilk_ever", "grainsroots_1d", "grainsroots_7d", "grainsroots_ever", "insects_1d",
                               "insects_7d", "insects_ever", "juice_1d", "juice_7d", "juice_ever", "legumesnuts_1d",
                               "legumesnuts_7d", "legumesnuts_ever", "mahewu_1d", "mahewu_7d", "mahewu_ever",
                               "mangopaw_1d", "mangopaw_7d", "mangopaw_ever", "meats_1d", "meats_7d", "meats_ever",
                               "nuts_1d", "nuts_7d", "nuts_ever", "oilfats_1d", "oilfats_7d", "oilfats_ever", "oralmed_1d",
                               "oralmed_7d", "oralmed_ever", "organmeats_1d", "organmeats_7d", "organmeats_ever", "ors_1d",
                               "ors_7d", "ors_ever", "othereggs_1d", "othereggs_7d", "othereggs_ever", "otherfruit_1d",
                               "otherfruit_7d", "otherfruit_ever", "otherfruitveg_1d", "otherfruitveg_7d",
                               "otherfruitveg_ever", "otherveg_1d", "otherveg_7d", "otherveg_ever", "porridge_1d",
                               "porridge_7d", "porridge_ever", "potatoyam_1d", "potatoyam_7d", "potatoyam_ever",
                               "poultry_1d", "poultry_7d", "poultry_ever", "pwdrmilk_1d", "pwdrmilk_7d", "redmeat_1d",
                               "rutf_1d", "rutf_7d", "rutf_ever", "sadza_1d", "sadza_7d", "sadza_ever", "sugar_1d",
                               "sugar_7d", "sugar_ever", "sugarwater_1d", "sugarwater_7d", "teasugar_1d", "teasugar_7d",
                               "teasugar_ever", "teasugmilk_1d", "teasugmilk_7d", "teasugmilk_ever", "timessolidfeed_1d",
                               "vitafoods_1d", "vitafoods_7d", "vitafoods_ever", "yelloworgeveg_1d", "yelloworgeveg_7d",
                               "yelloworgeveg_ever", "yogurt_1d", "yogurt_7d", "yogurt_ever")] <- NULL
metadata_robertson_combined <- metadata_robertson_combined %>%
  mutate(family_ID = HouseholdID,
         subject_ID = SubjectID_child,
         birthmode = case_when(AC38 == "normal vaginal" ~ "v",
                               AC38 == "forceps_vacuum" ~ "v",
                               AC38 == "Caesarean section" ~ "c"),
         birth_height = birthlength,
         birth_weight = birthweight,
         gestational_age = gaw_final,
         height = length,
         antibiotics_any = case_when(amoxicillin == "yes" |
                                    cotrimoxazole_adj == "yes" |
                                    otherantibiotic == "yes" |
                                    metronidazole == "yes" |
                                    penicillin == "yes" |
                                    antibiotics_exposed ~ T,
                                    amoxicillin == "no" &
                                    cotrimoxazole_adj == "no" &
                                    otherantibiotic == "no" &
                                    penicillin == "no" &
                                    gentamicin == "no" &
                                    nitrofurantoin == "no" & 
                                    metronidazole == "no" ~ F, .default = NA),
         food = ifelse(breastfeednow == "yes", yes = "breast", no = NA),
         geographic_location_.latitude. = -19.453424,
         geographic_location_.longitude. = 30.154827,
         .keep = "unused") %>%
  select_if(~ !all(is.na(.)))

# multiple samples for each individual exist
multiple_samples <- unique(metadata_robertson_combined$subject_ID[duplicated(metadata_robertson_combined$subject_ID)])
metadata_robertson_combined <- metadata_robertson_combined[metadata_robertson_combined$subject_ID %in% multiple_samples,]


# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_robertson_combined$read_count, n_lines$read_counts[match(metadata_robertson_combined$run_accession, n_lines$sample_ID)])  # 


sort_out_samples_rob <- ena_file_robertson$run_accession[!ena_file_robertson$run_accession %in% metadata_robertson_combined$run_accession]
# move sorted out samples
sort_out_samples_rob <- unique(sort_out_samples_rob)
sort_out_files_rob <- c(paste(sort_out_samples_rob, "_1.fastq.gz", sep = ""),
                         paste(sort_out_samples_rob, "_2.fastq.gz", sep = ""),
                         paste(sort_out_samples_rob, ".fastq.gz", sep = ""))
sort_out_files_rob <- sort_out_files_rob[file.exists(sort_out_files_rob)]
file.move(sort_out_files_rob, "sorted_out/")


final_samples_rob <- metadata_robertson_combined$run_accession
final_files_rob <- c(paste0(final_samples_rob, "_1.fastq.gz"),
                      paste0(final_samples_rob, "_2.fastq.gz"),
                      paste0(final_samples_rob, ".fastq.gz"))
final_files_rob <- final_files_rob[file.exists(final_files_rob)]
file.move(final_files_rob, "fastq_files/")

write.csv(metadata_robertson_combined, "metadata_robertson_2023_healthy.csv")
rm(ena_file_robertson, metadata_1, metadata_robertson_combined, sra_file, sort_out_samples_rob,
   sort_out_files_rob, final_samples_rob, final_files_rob, multiple_samples, n_lines)


# Valles-Colomer	2023 #########################################################
setwd(paste(maindir, "shotgun/valles_colomer_2023", sep = ""))
ena_file_valles1 <- read.table("filereport_read_run_PRJEB45799_tsv.txt", sep = "\t", header = T, comment.char = "")
ena_file_valles2 <- read.table("filereport_read_run_PRJNA613947_tsv.txt", sep = "\t", header = T, comment.char = "")
ena_file_valles3 <- read.table("filereport_read_run_PRJNA716780_tsv.txt", sep = "\t", header = T, comment.char = "")
ena_file_valles <- rbind(ena_file_valles1, ena_file_valles2, ena_file_valles3)
ena_file_valles <- ena_file_valles[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
# duplicated samples?
any(duplicated(ena_file_valles$sample_alias))

sra_file_1 <- read.table("SraRunTable_PRJEB45799.txt", sep = ",", header = T, comment.char = "", colClasses = "character", quote = "\"")
sra_file_2 <- read.table("SraRunTable_PRJNA613947.txt", sep = ",", header = T, comment.char = "", colClasses = "character", quote = "\"")
sra_file_3 <- read.table("SraRunTable_PRJNA716780.txt", sep = ",", header = T, comment.char = "", colClasses = "character", quote = "\"")
sra_file <- bind_rows(sra_file_1, sra_file_2, sra_file_3) %>%
  mutate_all(type.convert, as.is=TRUE)

metadata <- read.table("metadata_1.txt", sep = "\t", header = T)

metadata_valles_combined <- left_join(ena_file_valles, sra_file, by = join_by(run_accession == Run))
metadata_valles_combined <- left_join(metadata_valles_combined, metadata, by = join_by(Sample_name == sampleID))
metadata_valles_combined[c("library_name", "Assay.Type", "Consent", "ENA.FIRST.PUBLIC..run.",
                           "ENA_first_public", "ENA.LAST.UPDATE..run.", "ENA_last_update",
                           "INSDC_center_alias", "INSDC_center_name", "INSDC_first_public",
                           "INSDC_last_update", "INSDC_status", "LibraryLayout",
                           "LibrarySource", "Organism", "Platform",
                           "ReleaseDate", "Bytes", "DATASTORE.filetype", "DATASTORE.provider",
                           "DATASTORE.region", "create_date", "version", "BioSampleModel",
                           "Collection_Date", "env_broad_scale", "env_local_scale",
                           "env_medium", "geo_loc_name_country_continent", "HOST",
                           "Isolation_source", "sequencing_platform", "number_bases",
                           "minimum_read_length", "median_read_length", "X",
                           "age_category", "Antibiotics", "Dataset", "mother_subjectID",
                           "subjectID", "Timepoint")] <- NULL
metadata_valles_combined <- metadata_valles_combined %>%
  separate(Phenotype, into = c(NA, NA, NA, "age", NA), sep = "_", convert = T) %>%
  filter(is.na(Host_disease) | Host_disease == "healthy control",
         is.na(body_site) | body_site == "stool",
         is.na(Host_age) | Host_age <= 3,
         is.na(age) | age <= 3,
         !is.na(age) | !is.na(Host_age)) %>%  # remove missing age information
  dplyr::mutate(birthmode = case_when(is.na(mode_delivery) ~ NA,
                               mode_delivery == "c_section" ~ "c",
                               mode_delivery == "vaginal" ~ "v"),
         sex = case_when(gender == "female" ~ "f",
                         gender == "male" ~ "m",
                         host_sex == "female" ~ "f",
                         host_sex == "male" ~ "m"),
         age = ifelse(is.na(age), yes = (Host_age * 365), no = (age * 365)),
         country = case_when(country == "ARG" ~ "ARGENTINIA",
                             country == "ITA" ~ "ITALY",
                             country == "CHN" ~ "CHINA",
                             country == "COL" ~ "COLOMBIA",
                             country == "GNB" ~ "GUINEA_BISSAU",
                             geo_loc_name_country == "Italy" ~ "ITALY",
                             geo_loc_name_country == "China" ~ "CHINA"),
         subject_ID = ifelse(is.na(Subject_ID), yes = host_subject_id, no = Subject_ID),
         family_ID = ifelse(is.na(familyID), yes = subject_ID, no = familyID),
         study = "valles_colomer_2023",
         region = country,
         lifestyle = case_when(country == "ARGENTINIA" ~ "non_industrialized",
                               country == "ITALY" ~ "industrialized",
                               country == "CHINA" ~ "industrialized",
                               country == "COLOMBIA" ~ "non_industrialized",
                               country == "GUINEA_BISSAU" ~ "non_industrialized"),
         sample_ID = paste0("valles_colomer_2023_", sample_alias),
         Instrument = ifelse(Instrument == "unspecified", "Illumina NovaSeq 6000", Instrument),
         geographic_location_.latitude.  = case_when(country == "ITALY" ~ 46.066666,
                                                     country == "CHINA" ~ 31.2304,
                                                     country == "ARGENTINIA" ~ -29.202528,
                                                     country == "COLOMBIA" ~ 11.4113,
                                                     country == "GUINEA_BISSAU" ~ 11.2489),
         geographic_location_.longitude. = case_when(country == "ITALY" ~ 11.116667,
                                                     country == "CHINA" ~ 121.4737,
                                                     country == "ARGENTINIA" ~ -61.7260,
                                                     country == "COLOMBIA" ~ -72.9062,
                                                     country == "GUINEA_BISSAU" ~ -15.8616),
         .keep = "unused") %>%
  select_if(~ !all(is.na(.))) %>%
  filter(!is.na(subject_ID))
metadata_valles_combined[,c("Host_disease", "body_site", "LibrarySelection",
                            "Submitter_Id", "AvgSpotLen", "Bases", "zygosity",
                            "non_Westernized", "SRA.Study")] <- NULL

# only keep subjects with multiple samples:
multiple_samples <- unique(metadata_valles_combined$subject_ID[duplicated(metadata_valles_combined$subject_ID)])
# metadata_valles_combined <- metadata_valles_combined[metadata_valles_combined$subject_ID %in% multiple_samples,]
metadata_valles_combined$run_accession[!metadata_valles_combined$subject_ID %in% multiple_samples] %>%
  paste("fastq-dump --split-e --gzip ",.) %>%
  data.frame(`#!/bin/bash` = .) %>%
  write.table(., file = "/fast/AG_Forslund/rob/studies/shotgun/valles_colomer_2023/pull_single.sh",
              sep = " ", quote = F, row.names = F, col.names = "#!/bin/bash")


# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_valles_combined$read_count, n_lines$read_counts[match(metadata_valles_combined$run_accession, n_lines$sample_ID)])  # 


sort_out_samples_val <- ena_file_valles$run_accession[!ena_file_valles$run_accession %in% metadata_valles_combined$run_accession]
# move sorted out samples
sort_out_files_val <- c(paste(sort_out_samples_val, "_1.fastq.gz", sep = ""),
                        paste(sort_out_samples_val, "_2.fastq.gz", sep = ""),
                        paste(sort_out_samples_val, ".fastq.gz", sep = ""))
sort_out_files_val <- sort_out_files_val[file.exists(sort_out_files_val)]
file.move(sort_out_files_val, "sorted_out/")


final_samples_val <- metadata_valles_combined$run_accession
final_files_val <- c(paste0(final_samples_val, "_1.fastq.gz"),
                     paste0(final_samples_val, "_2.fastq.gz"),
                     paste0(final_samples_val, ".fastq.gz"))
final_files_val <- final_files_val[file.exists(final_files_val)]
file.move(final_files_val, "fastq_files/")

write.csv(metadata_valles_combined, "metadata_valles_colomer_2023_healthy.csv")
rm(metadata_valles_combined, sra_file, sra_file_1, sra_file_2, sra_file_3, final_files_val,
   final_samples_val, sort_out_files_val, sort_out_samples_val, metadata, ena_file_valles,
   ena_file_valles1, ena_file_valles2, ena_file_valles3, multiple_samples, n_lines)



# Pärnänen	2018 ###############################################################
setwd(paste(maindir, "shotgun/parnanen_2018", sep = ""))
ena_file_par <- read.table("filereport_read_run_PRJNA384716_tsv.txt", sep = "\t", header = T, comment.char = "")
ena_file_par <- ena_file_par[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T, comment.char = "")
sra_file[,c("Assay.Type", "AvgSpotLen", "BioSampleModel", "Bytes", "Center.Name",
            "collection_date", "Consent", "DATASTORE.filetype", "DATASTORE.provider",
            "DATASTORE.region", "geo_loc_name_country", "geo_loc_name_country_continent",
            "geo_loc_name", "HOST", "LibraryLayout", "LibrarySelection", "LibrarySource",
            "Organism", "Platform", "ReleaseDate", "create_date", "version"
)] <- NULL

metadata_combined_par <- merge(ena_file_par, sra_file, by.x = "run_accession", by.y = "Run")
metadata_1 <- read.table("metadata_raw.tsv", sep = "\t", header = T)
metadata_combined_par <- merge(metadata_combined_par, metadata_1, by.x = "library_name", by.y = "SAMPLE_ID")
all(metadata_combined_par$library_name == metadata_combined_par$sample_alias)
all(metadata_combined_par$library_name == metadata_combined_par$Sample.Name)
all(metadata_combined_par$library_name == metadata_combined_par$Library.Name)
all(metadata_combined_par$library_name == metadata_combined_par$unique_identifier)
metadata_combined_par[,c("sample_alias", "Sample.Name", "Library.Name", "unique_identifier",
                         "BFAT6MO")] <- NULL
# only children, remove duplicate samples, keep the one with the higher read count
metadata_combined_par <- metadata_combined_par %>% 
  filter(isolation_source == "feces",
         TIME != "32WK",
         TYPE != "AU2")

metadata_combined_par %<>% mutate(
  food = case_when(
    (TIME == "1MONTH") & (BFMONTHS >= 1) ~ "breast",
    (TIME == "1MONTH") & (BFMONTHS < 1) ~ "formula",
    (TIME == "6MONTH") & (BFMONTHS >= 6) ~ "breast",
    (TIME == "6MONTH") & (BFMONTHS < 6) ~ "formula"
  ),
  age = case_when(
    TIME == "1MONTH" ~ 30,
    TIME == "6MONTH" ~ 180),
  gestational_age = GE.Weeks.,
  birthmode = "v",
  country = "FINNLAND",
  region = "Turku",
  lifestyle = "industrialized",
  study = "parnanen_2018",
  sample_ID = library_name,
  subject_ID = PAIR,
  family_ID = PAIR,
  geographic_location_.latitude. = 60.45,
  geographic_location_.longitude. = 22.27,
  .keep = "unused")

# keep only individuals with more than 1 sample:
multiple_samples <- unique(metadata_combined_par$subject_ID[duplicated(metadata_combined_par$subject_ID)])
metadata_combined_par <- metadata_combined_par[metadata_combined_par$subject_ID %in% multiple_samples,]


n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_combined_par$Bases, n_lines$read_counts[match(metadata_combined_par$run_accession, n_lines$sample_ID)])  # 


sort_out_samples_par <- ena_file_par$run_accession[!ena_file_par$run_accession %in% metadata_combined_par$run_accession]
# move sorted out samples
sort_out_samples_par <- unique(sort_out_samples_par)
sort_out_files_par <- c(paste(sort_out_samples_par, "_1.fastq.gz", sep = ""),
                           paste(sort_out_samples_par, "_2.fastq.gz", sep = ""),
                           paste(sort_out_samples_par, ".fastq.gz", sep = ""))
sort_out_files_par <- sort_out_files_par[file.exists(sort_out_files_par)]
file.move(sort_out_files_par, "sorted_out/")


final_samples_par <- metadata_combined_par$run_accession
final_files_par <- c(paste0(final_samples_par, "_1.fastq.gz"),
                        paste0(final_samples_par, "_2.fastq.gz"),
                        paste0(final_samples_par, ".fastq.gz"))
final_files_par <- final_files_par[file.exists(final_files_par)]
file.move(final_files_par, "fastq_files/")

write.csv(metadata_combined_par, "metadata_parnanen_2018_healthy.csv")
rm(ena_file_par, metadata_1, metadata_combined_par, n_lines, sra_file, 
   final_samples_par, final_files_par, sort_out_files_par, sort_out_samples_par)



# Vatanen 2022 2 ###############################################################

setwd(paste(maindir, "shotgun/vatanen_2022_2", sep = ""))
ena_file_vat <- read.table("filereport_read_run_PRJNA821542_tsv_filtered.txt", sep = "\t", header = T)
ena_file_vat <- ena_file_vat[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file[,c( "Assay.Type", "AvgSpotLen", "BioSampleModel", "Bytes", "Center.Name",
             "Collection_Date", "Consent", "DATASTORE.filetype", "DATASTORE.provider",
             "DATASTORE.region", "env_biome", "env_feature", "env_material", "Host",
             "geo_loc_name_country", "geo_loc_name_country_continent", "geo_loc_name",
             "LibraryLayout", "LibrarySelection", "LibrarySource", "Organism", "Platform",
             "ReleaseDate", "create_date", "version")] <- NULL
metadata_vat_combined <- merge(ena_file_vat, sra_file,by.x = "run_accession", by.y = "Run")



metadata_vat_combined <- metadata_vat_combined %>% separate(misc_param, c(NA, NA, "family_ID"), sep = " ", remove = F) %>%
  mutate(country = "FINNLAND",
         age = case_when(Host_age == "0.5 months" ~ 15,
                         Host_age == "1 months" ~ 30,
                         Host_age == "10 months" ~ 300,
                         Host_age == "11 months" ~ 330,
                         Host_age == "12 months" ~ 365,
                         Host_age == "2 months" ~ 60,
                         Host_age == "3 months" ~ 90,
                         Host_age == "4 months" ~ 120,
                         Host_age == "5 months" ~ 150,
                         Host_age == "6 months" ~ 180,
                         Host_age == "7 months" ~ 210,
                         Host_age == "8 months" ~ 240,
                         Host_age == "9 months" ~ 270,
                         Host_age == "infant age 12 months" ~ 365,
                         Host_age == "infant age 18 months" ~ 545,
                         Host_age == "infant age 3 months" ~ 90,
                         Host_age == "infant birth" ~ 0,
                         Host_age == "missing" ~ NA,
                         Host_age == "near birth" ~ 1,
                         Host_age == "third trimester of pregnancy" ~ NA),
         study = "vatanen_2022_2",
         lifestyle = "industrialized",
         region = "TAMPERE",
         sample_ID = Sample.Name,
         subject_ID = family_ID,
         family_ID = family_ID,
         geographic_location_.latitude. = 61.4978,
         geographic_location_.longitude. = 23.7610,
         diarrhea = ifelse(grepl("DS", age), yes = "y", no = "n"),
         .keep = "unused") %>%
  filter(!is.na(age))
any(duplicated(metadata_vat_combined$sample_ID))


# remove individuals with only one sample
multiple_samples <- unique(metadata_vat_combined$subject_ID[duplicated(metadata_vat_combined$subject_ID)])
metadata_vat_combined <- metadata_vat_combined %>%
  filter(subject_ID %in% multiple_samples)


# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_vat_combined$read_count, n_lines$read_counts[match(metadata_vat_combined$run_accession, n_lines$sample_ID)])  # 

write.csv(metadata_vat_combined, "metadata_vatanen_2022_2_healthy.csv")

# # move sorted out samples
sort_out_samples <- ena_file_vat$run_accession[!ena_file_vat$run_accession %in% metadata_vat_combined$run_accession]
sort_out_files <- c(paste(sort_out_samples, "_1.fastq.gz", sep = ""),
                    paste(sort_out_samples, "_2.fastq.gz", sep = ""),
                    paste(sort_out_samples, ".fastq.gz", sep = ""))
sort_out_files <- sort_out_files[file.exists(sort_out_files)]
file.move(sort_out_files, "sorted_out/")

# move good files from shotgun folder
final_samples <- metadata_vat_combined$run_accession
final_files <- c(paste0(final_samples, "_1.fastq.gz"),
                 paste0(final_samples, "_2.fastq.gz"),
                 paste0(final_samples, ".fastq.gz"))
final_files <- final_files[file.exists(final_files)]
file.move(final_files, "fastq_files/")


# clean up
rm(ena_file_vat, metadata_vat_combined, final_files, final_samples,
   multiple_samples, sort_out_samples, sort_out_files, n_lines, sra_file)


# Hoskinson	2023 ###############################################################

# setwd(paste(maindir, "shotgun/hoskinson_2023", sep = ""))
# ena_file_hos <- read.table("filereport_read_run_PRJNA838575_tsv.txt", sep = "\t", header = T)
# ena_file_hos <- ena_file_hos[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
# sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
# sra_file[,c( "Assay.Type", "AvgSpotLen", "BioSampleModel", "Bytes", "Center.Name",
#              "Collection_Date", "DATASTORE.filetype", "consent", "DATASTORE.provider", "DATASTORE.region",
#              "env_broad_scale", "env_local_scale", "env_medium", "geo_loc_name_country",
#              "geo_loc_name_country_continent", "geo_loc_name", "Host", "lat_lon",
#              "LibraryLayout", "LibrarySelection", "LibrarySource", "Organism",
#              "Platform", "ReleaseDate", "create_date")] <- NULL
# metadata_hos_combined <- merge(ena_file_hos, sra_file,by.x = "run_accession", by.y = "Run")
# 
# 
# metadata_hos_combined <- metadata_hos_combined %>%
#   mutate(country = "CANNADA",
#          age = case_when(Host_age == "1 year" ~ 365,
#                          Host_age == "13 month" ~ 90),
#          study = "hoskinson_2023",
#          lifestyle = "industrialized",
#          region = "CANNADA",
#          sample_ID = Sample.Name,
#          subject_ID = family_ID,
#          family_ID = family_ID,
#          geographic_location_.latitude. = 43.651070,
#          geographic_location_.longitude. = -79.347015) %>%
#   filter(!is.na(age))
# any(duplicated(metadata_hos_combined$sample_ID))
# 
# 
# # remove individuals with only one sample
# multiple_samples <- unique(metadata_vat_combined$subject_ID[duplicated(metadata_vat_combined$subject_ID)])
# metadata_vat_combined <- metadata_vat_combined %>%
#   filter(subject_ID %in% multiple_samples)
# 
# 
# # check read numbers
# n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
# plot(metadata_vat_combined$read_count, n_lines$read_counts[match(metadata_vat_combined$run_accession, n_lines$sample_ID)])  # 
# 
# write.csv(metadata_vat_combined, "metadata_vatanen_2022_2_healthy.csv")
# 
# # # move sorted out samples
# sort_out_samples <- ena_file_vat$run_accession[!ena_file_vat$run_accession %in% metadata_vat_combined$run_accession]
# sort_out_files <- c(paste(sort_out_samples, "_1.fastq.gz", sep = ""),
#                     paste(sort_out_samples, "_2.fastq.gz", sep = ""),
#                     paste(sort_out_samples, ".fastq.gz", sep = ""))
# sort_out_files <- sort_out_files[file.exists(sort_out_files)]
# file.move(sort_out_files, "sorted_out/")
# 
# # move good files from shotgun folder
# final_samples <- metadata_vat_combined$run_accession
# final_files <- c(paste0(final_samples, "_1.fastq.gz"),
#                  paste0(final_samples, "_2.fastq.gz"),
#                  paste0(final_samples, ".fastq.gz"))
# final_files <- final_files[file.exists(final_files)]
# file.move(final_files, "fastq_files/")
# 
# 
# # clean up
# rm(ena_file_vat, metadata_vat_combined, final_files, final_samples,
#    multiple_samples, sort_out_samples, sort_out_files, n_lines, sra_file)



# Bargheet	2023 ###############################################################

setwd(paste(maindir, "shotgun/bargheet_2023", sep = ""))
ena_file <- read.table("filereport_read_run_PRJNA898628_tsv_filtered.txt", sep = "\t", header = T)
ena_file <- ena_file[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file[,c( "Assay.Type", "AvgSpotLen", "BioSampleModel", "Bytes", "Center.Name",
             "Collection_Date", "Consent", "DATASTORE.filetype", "DATASTORE.provider",
             "DATASTORE.region", "env_broad_scale", "env_local_scale",
             "env_medium", "geo_loc_name_country", "geo_loc_name_country_continent",
             "geo_loc_name", "Host", "lat_lon",
             "LibraryLayout", "LibrarySelection", "LibrarySource", "Organism", "Platform",
             "ReleaseDate", "create_date", "version")] <- NULL
metadata_bargh_combined <- merge(ena_file, sra_file,by.x = "run_accession", by.y = "Run")


metadata_bargh_combined <- metadata_bargh_combined %>% separate(sample_alias,
                                                                c(NA, "subject_ID", "sample", "time_point"),
                                                                sep = "-",
                                                                remove = F,
                                                                convert = F) %>%
   mutate(country = "NORWAY",
          age = case_when(time_point == "Day120_R1" ~ 120,
                          time_point == "Day28_R1" ~ 28,
                          time_point == "Day365_R1" ~ 365,
                          time_point == "Day7_R1" ~ 7),
          study = "bargheet_2023",
          antibiotics_any = F,
          antibiotics_before = F,
          antibiotics_one_week_before = F,
          birthmode = "v",
          lifestyle = "industrialized",
          region = "TROMSO",
          sample_ID = sample_alias,
          family_ID = subject_ID,
          geographic_location_.latitude. = 60.4720 ,
          geographic_location_.longitude. = 8.4689)
any(duplicated(metadata_bargh_combined$sample_ID))


# remove individuals with only one sample
multiple_samples <- unique(metadata_bargh_combined$subject_ID[duplicated(metadata_bargh_combined$subject_ID)])
metadata_bargh_combined <- metadata_bargh_combined %>%
   filter(subject_ID %in% multiple_samples)


# check read numbers
n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
plot(metadata_bargh_combined$read_count, n_lines$read_counts[match(metadata_bargh_combined$run_accession, n_lines$sample_ID)])  # 

write.csv(metadata_bargh_combined, "metadata_bargheet_2023_healthy.csv")

# # move sorted out samples
sort_out_samples <- ena_file$run_accession[!ena_file$run_accession %in% metadata_bargh_combined$run_accession]
sort_out_files <- c(paste(sort_out_samples, "_1.fastq.gz", sep = ""),
                    paste(sort_out_samples, "_2.fastq.gz", sep = ""),
                    paste(sort_out_samples, ".fastq.gz", sep = ""))
sort_out_files <- sort_out_files[file.exists(sort_out_files)]
file.move(sort_out_files, "sorted_out/")

# move good files from shotgun folder
final_samples <- metadata_bargh_combined$run_accession
final_files <- c(paste0(final_samples, "_1.fastq.gz"),
                 paste0(final_samples, "_2.fastq.gz"),
                 paste0(final_samples, ".fastq.gz"))
final_files <- final_files[file.exists(final_files)]
file.move(final_files, "fastq_files/")

# clean up
rm( "ena_file", "final_files", "final_samples", "metadata_bargh_combined",
    "multiple_samples", "sort_out_files", "sort_out_samples", "sra_file", n_lines )


# Patangia	2024 ###############################################################

setwd(paste(maindir, "shotgun/patangia_2024", sep = ""))
ena_file <- read.table("filereport_read_run_PRJNA971895_tsv.txt", sep = "\t", header = T)
ena_file <- ena_file[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file[,c( "Assay.Type", "AvgSpotLen", "BioSampleModel", "Bytes", "Center.Name",
             "Collection_Date", "Consent", "DATASTORE.filetype", "DATASTORE.provider",
             "DATASTORE.region", "forward", "geo_loc_name_country", "geo_loc_name_country_continent",
             "geo_loc_name", "host_developmental_stage", "Host", "isolation_source", "lat_lon",
             "LibraryLayout", "LibrarySelection", "LibrarySource", "Organism", "Platform",
             "ReleaseDate","Reverse", "create_date", "version", "Library_ID", "Library.Name")] <- NULL
metadata_pat_combined <- merge(ena_file, sra_file,by.x = "run_accession", by.y = "Run")

metadata_pat_combined <- metadata_pat_combined %>%
  filter(Group..run. != "blank_r") %>%
  separate(sample_alias,
           c("subject_ID", NA),
           sep = "[wy]",
           remove = F,
           convert = F) %>%
  mutate(subject_ID = gsub("[IM]", "", subject_ID),
         age = case_when(Time..run. == "w1" ~ 7,
                         Time..run. == "w4" ~ 28,
                         Time..run. == "w8" ~ 56,
                         Time..run. == "w24" ~ 168,
                         Time..run. == "w104" ~ 728),
         study = "bargheet_2023",
         antibiotics_before = case_when(Group..run. == "CS_noab" ~ F,
                                        Group..run. == "CS_ab" ~ T,
                                        Group..run. == "VD_noab" ~ F),
         antibiotics_any = antibiotics_before,
         birthmode = case_when(Group..run. == "CS_noab" ~ "c",
                               Group..run. == "CS_ab" ~ "c",
                               Group..run. == "VD_noab" ~ "v"),
         country = "IRELAND",
         region = "CORK",
         study = "patangia_2024",
         lifestyle = "industrialized",
         geographic_location_.latitude. = 51.903614,
         geographic_location_.longitude. = -8.468399,
         .keep = "unused") %>%
  mutate(family_ID = subject_ID,
         sample_ID = gsub("[IM]","", sample_alias) %>% paste(.,"patangia_2024", sep = "_"))

# check metadata from hill:
metadata_hil_combined <- read.csv("/fast/AG_Forslund/rob/studies/16S/hill_2017/metadata_hil_2017_healthy.csv")
metadata_hil_combined <- metadata_hil_combined %>% mutate(subject_ID = gsub("_hill_2017", "", subject_ID)) %>%
  select(subject_ID, gestational_age, birth_weight, sex) %>%
  distinct
metadata_pat_combined$subject_ID %in% metadata_hil_combined$subject_ID

metadata_pat_combined <- left_join(metadata_pat_combined, metadata_hil_combined, by = "subject_ID")


# remove individuals with only one sample
multiple_samples <- unique(metadata_pat_combined$subject_ID[duplicated(metadata_pat_combined$subject_ID)])
metadata_pat_combined <- metadata_pat_combined %>%
  filter(subject_ID %in% multiple_samples)


# check read numbers
# n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
# plot(metadata_pat_combined$read_count, n_lines$read_counts[match(metadata_pat_combined$run_accession, n_lines$sample_ID)])  # 

write.csv(metadata_pat_combined, "metadata_patangia_2024_healthy.csv")

# # move sorted out samples
sort_out_samples <- ena_file$run_accession[!ena_file$run_accession %in% metadata_pat_combined$run_accession]
sort_out_files <- c(paste(sort_out_samples, "_1.fastq.gz", sep = ""),
                    paste(sort_out_samples, "_2.fastq.gz", sep = ""),
                    paste(sort_out_samples, ".fastq.gz", sep = ""))
sort_out_files <- sort_out_files[file.exists(sort_out_files)]
file.move(sort_out_files, "sorted_out/")

# move good files from shotgun folder
final_samples <- metadata_pat_combined$run_accession
final_files <- c(paste0(final_samples, "_1.fastq.gz"),
                 paste0(final_samples, "_2.fastq.gz"),
                 paste0(final_samples, ".fastq.gz"))
final_files <- final_files[file.exists(final_files)]
file.move(final_files, "fastq_files/")

# clean up
rm("ena_file", "final_files", "final_samples", "metadata_hil_combined", "metadata_pat_combined", 
   "multiple_samples", "sort_out_files", "sort_out_samples", "sra_file", n_lines)


# Selma-Royo	2024 ###############################################################

setwd(paste(maindir, "shotgun/selma_royo_2024", sep = ""))
ena_file <- read.table("filereport_read_run_PRJEB74322_tsv.txt", sep = "\t", header = T)
ena_file <- ena_file[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
# sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
# sra_file[,c( "Assay.Type", "AvgSpotLen", "BioSampleModel", "Bytes", "Center.Name",
#              "Collection_Date", "Consent", "DATASTORE.filetype", "DATASTORE.provider",
#              "DATASTORE.region", "forward", "geo_loc_name_country", "geo_loc_name_country_continent",
#              "geo_loc_name", "host_developmental_stage", "Host", "isolation_source", "lat_lon",
#              "LibraryLayout", "LibrarySelection", "LibrarySource", "Organism", "Platform",
#              "ReleaseDate","Reverse", "create_date", "version", "Library_ID", "Library.Name")] <- NULL
# metadata_sel_combined <- merge(ena_file, sra_file,by.x = "run_accession", by.y = "Run")
# 
# 
# write.csv(metadata_sel_combined, "metadata_selma_royo_2024_healthy.csv")
# 
# # # move sorted out samples
# sort_out_samples <- ena_file$run_accession[!ena_file$run_accession %in% metadata_sel_combined$run_accession]
# sort_out_files <- c(paste(sort_out_samples, "_1.fastq.gz", sep = ""),
#                     paste(sort_out_samples, "_2.fastq.gz", sep = ""),
#                     paste(sort_out_samples, ".fastq.gz", sep = ""))
# sort_out_files <- sort_out_files[file.exists(sort_out_files)]
# file.move(sort_out_files, "sorted_out/")
# 
# # move good files from shotgun folder
# final_samples <- metadata_sel_combined$run_accession
# final_files <- c(paste0(final_samples, "_1.fastq.gz"),
#                  paste0(final_samples, "_2.fastq.gz"),
#                  paste0(final_samples, ".fastq.gz"))
# final_files <- final_files[file.exists(final_files)]
# file.move(final_files, "fastq_files/")
# 
# # clean up
# rm("ena_file", "final_files", "final_samples", "metadata_hil_combined", "metadata_pat_combined", 
#    "multiple_samples", "sort_out_files", "sort_out_samples", "sra_file")

# Dubois	2024 ###############################################################

setwd(paste(maindir, "shotgun/dubois_2024", sep = ""))
ena_file <- read.table("filereport_read_run_PRJEB52774_tsv.txt", sep = "\t", header = T)
ena_file <- ena_file[c("study_accession", "run_accession", "read_count", "sample_alias", "sample_title")]

sra_file <- read.table("SraRunTable_PRJEB52774.txt", sep = ",", header = T)
sra_file[,c( "Assay.Type", "AvgSpotLen", "Bytes", "Center.Name",
             "Collection_Date", "Consent", "DATASTORE.filetype", "DATASTORE.provider",
             "DATASTORE.region", "ENA.FIRST.PUBLIC..run.", "ENA.LAST.UPDATE..run.",
             "geo_loc_name_country", "geo_loc_name_country_continent",
             "geographic_location_.country_and.or_sea.", "host_body_product",
             "human_gut_environmental_package", "INSDC_center_name",
             "Investigation_type", "Library.Name", "LibraryLayout", "LibrarySelection",
             "LibrarySource", "Organism", "Platform", "ReleaseDate", "create_date",
             "version", "Sequencing_method", "INSDC_first_public", "INSDC_last_update",
             "INSDC_status", "INSDC_center_alias", "environment_.biome.", "environment_.feature.",
             "environment_.material.", "ena_first_public", "ENA_last_update",
             "Scientific_Name", "broad.scale_environmental_context", "environmental_medium",
             "local_environmental_context", "ENA.FIRST.PUBLIC", "ENA.LAST.UPDATE",
             "Cohort")] <- NULL

metadata_dub_combined <- merge(ena_file, sra_file,by.x = "run_accession", by.y = "Run")
xls_file <- readxl::read_xlsx("mmc2.xlsx", sheet = "Table S1", skip = 1)
metadata_dub_combined <- left_join(metadata_dub_combined, xls_file,by = c("sample_alias" = "ENA_SampleAlias"))

metadata_dub_combined <- metadata_dub_combined %>%
  filter(sample_title == "infant stool sample", infant_HadProbiotics == "FALSE") %>%
  mutate(age = Host_age * 7,
         subject_ID = host_subject_id,
         sample_ID = Sample_Name %>% paste(.,"dubois_2024", sep = "_"),
         study = "dubois_2024",
         birthmode = case_when(infant_DeliveryMode == "Vaginal" ~ "v",
                               infant_DeliveryMode == "C-section" ~ "c"),
         country = "FINNLAND",
         region = "FINNLAND",
         study = "dubois_2024",
         lifestyle = "industrialized",
         sex = ifelse(host_sex == "male", yes = "m", no = "f"),
         .keep = "unused") %>%
  mutate(family_ID = subject_ID)

# check read numbers
# n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
# plot(metadata_dub_combined$read_count, n_lines$read_counts[match(metadata_pat_combined$run_accession, n_lines$sample_ID)])  # 

write.csv(metadata_dub_combined, "metadata_dubois_2024_healthy.csv")

# move sorted out samples
sort_out_samples <- ena_file$run_accession[!ena_file$run_accession %in% metadata_dub_combined$run_accession]
sort_out_files <- c(paste(sort_out_samples, "_1.fastq.gz", sep = ""),
                    paste(sort_out_samples, "_2.fastq.gz", sep = ""),
                    paste(sort_out_samples, ".fastq.gz", sep = ""))
sort_out_files <- sort_out_files[file.exists(sort_out_files)]
file.move(sort_out_files, "sorted_out/")

# move good files from shotgun folder
final_samples <- metadata_dub_combined$run_accession
final_files <- c(paste0(final_samples, "_1.fastq.gz"),
                 paste0(final_samples, "_2.fastq.gz"),
                 paste0(final_samples, ".fastq.gz"))
final_files <- final_files[file.exists(final_files)]
file.move(final_files, "fastq_files/")

# clean up
rm("ena_file", "final_files", "final_samples", "metadata_dub_combined", "xls_file",
   "multiple_samples", "sort_out_files", "sort_out_samples", "sra_file")

# Manara	2023 ###############################################################

setwd(paste(maindir, "shotgun/manara_2023", sep = ""))
ena_file <- read.table("filereport_read_run_PRJNA504891_tsv.txt", sep = "\t", header = T)
ena_file <- ena_file[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file %>%
  select(Run, Experiment, Instrument,Library.Name, Organism, 
         c("Bases", "BioProject", "BioSample", "Sample.Name", "SRA.Study", "Host",
           "SubjectID"))
metadata_man_combined <- merge(ena_file, sra_file,by.x = "run_accession", by.y = "Run")
xls_file <- readxl::read_xlsx("1-s2.0-S0960982223004591-mmc2.xlsx") %>%
  column_to_rownames("sampleID") %>% t() %>%
  data.frame
xls_file$sample_ID <- rownames(xls_file)
metadata_man_combined <- left_join(metadata_man_combined, xls_file, by = c("library_name" = "sample_ID"))


metadata_man_combined <- metadata_man_combined %>%
  filter(Host == "Homo sapiens", as.numeric(age) <= 3) %>%
  mutate(age = as.numeric(age) * 365,
         subject_ID = SubjectID,
         study = "manara_2023",
         family_ID = family,
         sample_ID = library_name,
         birthmode = "v",
         country = "ETHIOPIA",
         region = "ETHIOPIA",
         lifestyle = "non_industrialized",
         geographic_location_.latitude. = 8.826290727911537,
         geographic_location_.longitude. = 39.174184044654815,
         .keep = "unused") %>%
  select(-westernization, -mother_infant, -age_category, -body_site, -Organism, -delivery_mode)

# check read numbers
# n_lines <- read.table("read_counts_raw.txt", header = F, col.names = c("sample_ID", "read_counts"), sep = "\t")
# plot(metadata_pat_combined$read_count, n_lines$read_counts[match(metadata_pat_combined$run_accession, n_lines$sample_ID)])  # 

write.csv(metadata_man_combined, "metadata_manara_2023_healthy.csv")

# # move sorted out samples
sort_out_samples <- ena_file$run_accession[!ena_file$run_accession %in% metadata_man_combined$run_accession]
sort_out_files <- c(paste(sort_out_samples, "_1.fastq.gz", sep = ""),
                    paste(sort_out_samples, "_2.fastq.gz", sep = ""),
                    paste(sort_out_samples, ".fastq.gz", sep = ""))
sort_out_files <- sort_out_files[file.exists(sort_out_files)]
file.move(sort_out_files, "sorted_out/")

# move good files from shotgun folder
final_samples <- metadata_man_combined$run_accession
final_files <- c(paste0(final_samples, "_1.fastq.gz"),
                 paste0(final_samples, "_2.fastq.gz"),
                 paste0(final_samples, ".fastq.gz"))
final_files <- final_files[file.exists(final_files)]
file.move(final_files, "fastq_files/")

# clean up
rm("ena_file", "final_files", "final_samples", "metadata_hil_combined", "metadata_man_combined", 
   "multiple_samples", "sort_out_files", "sort_out_samples", "sra_file", "xls_file")



# Xiao	2021 ###########################
setwd(paste(maindir, "shotgun/xiao_2021", sep = ""))
ena_file <- read.table("filereport_read_run_PRJNA695070_tsv.txt", sep = "\t", header = T)
ena_file <- ena_file[c("study_accession", "run_accession", "library_name", "read_count", "sample_alias", "sample_title")]
sra_file <- read.table("SraRunTable.txt", sep = ",", header = T)
sra_file <- sra_file %>%
  select(c("Run", "Bases", "BioProject", "BioSample", "Experiment", "Instrument",
           "Library.Name", "Sample.Name", "SRA.Study"))
xls_file <- readxl::read_xlsx("ChinaMetadata.xlsx") 
metadata_xiao_combined <- merge(ena_file, sra_file,by.x = "run_accession", by.y = "Run") %>%
  left_join(xls_file, ., by = c("SampleID" = "run_accession"))

metadata_xiao_combined <- metadata_xiao_combined %>%
  mutate(run_accession = SampleID,
         age = days,
         subject_ID = SubjectID,
         study = "xiao_2021",
         family_ID = subject_ID,
         sample_ID = Library.Name,
         birthmode = ifelse(deliveryMode == "Vaginal", yes = "v", no = "c"),
         country = "CHINA",
         region = "WENZHOU",
         lifestyle = "industrialized",
         antibiotics_any = F,
         antibiotics_before = F,
         antibiotics_one_week_before = F,
         geographic_location_.latitude. = 27.994267,
         geographic_location_.longitude. = 120.699364,
         .keep = "unused")
write.csv(metadata_xiao_combined, "metadata_xiao_2021_healthy.csv")

# # move sorted out samples
sort_out_samples <- ena_file$run_accession[!ena_file$run_accession %in% metadata_xiao_combined$run_accession]
sort_out_files <- c(paste(sort_out_samples, "_1.fastq.gz", sep = ""),
                    paste(sort_out_samples, "_2.fastq.gz", sep = ""),
                    paste(sort_out_samples, ".fastq.gz", sep = ""))
sort_out_files <- sort_out_files[file.exists(sort_out_files)]
file.move(sort_out_files, "sorted_out/")

# move good files from shotgun folder
final_samples <- metadata_xiao_combined$run_accession
final_files <- c(paste0(final_samples, "_1.fastq.gz"),
                 paste0(final_samples, "_2.fastq.gz"),
                 paste0(final_samples, ".fastq.gz"))
final_files <- final_files[file.exists(final_files)]
file.move(final_files, "fastq_files/")

# clean up
rm("ena_file", "final_files", "final_samples", "metadata_xiao_combined", 
   "sort_out_files", "sort_out_samples", "sra_file", "xls_file")

# Portlock	2025   ##########################
setwd(paste(maindir, "shotgun/portlock_2024", sep = ""))
sra_file <- read.table("SraRunTable.csv", sep = ",", header = T) # only one timepoint
rm(sra_file)

# Bottino	2025   ##########################
setwd(paste(maindir, "shotgun/bottino_2025", sep = ""))
ena_file <- read.table("filereport_read_run_PRJNA1128723_tsv.txt", sep = "\t", header = T)
sra_file <- read.table("SraRunTable.csv", sep = ",", header = T)
metadat_bottino <- readxl::read_xlsx("41467_2025_56072_MOESM5_ESM.xlsx",
                                     sheet = "Figure1B", skip = 2, col_names = T) %>%
  filter(datasource == "1kDLEAP-KHULA")
# Khula study is this one

metadat_bottino_combined <- sra_file %>% left_join(., metadat_bottino, by = c("Library.Name" = "sample")) %>%
  select(-c("datacolor", "Assay.Type", "AvgSpotLen", "BioSampleModel", "Bytes", "Center.Name",
            "Collection_Date", "Consent", "DATASTORE.filetype", "DATASTORE.provider", "DATASTORE.region",
            "geo_loc_name_country", "geo_loc_name_country_continent", "geo_loc_name", "HOST",
            "study_name", "Organism", "Platform", "ReleaseDate", "create_date", "version",
            "samp_collect_device", "LibraryLayout", "LibrarySelection", "LibrarySource", "isolation_source")) %>%
  mutate(study = "bottino_2025",
         run_accession = Run,
         age = ageMonths * 30,
         study = "bottino_2025",
         geographic_location_.latitude. = -33.95,
         geographic_location_.longitude. = 18.46,
         country = "SOUTH_AFRICA",
         region = "CAPE_TOWN",
         sample_ID = Sample.Name,
         subject_ID = subject_id,
         family_ID = subject_ID,
         lifestyle = "non_industrialized")

write.csv(metadat_bottino_combined, "metadata_bottino_2025_healthy.csv")

# move good files from shotgun folder
final_samples <- metadat_bottino_combined$run_accession
final_files <- c(paste0(final_samples, "_1.fastq.gz"),
                 paste0(final_samples, "_2.fastq.gz"),
                 paste0(final_samples, ".fastq.gz"))
final_files <- final_files[file.exists(final_files)]
file.move(final_files, "fastq_files/")

# clean up
rm("ena_file", "final_files", "final_samples", metadat_bottino_combined, 
   "sra_file", metadat_bottino)


# Chatzigiannidou	2025   ##########################
setwd(paste(maindir, "shotgun/chatzigiannidou_2025", sep = ""))
ena_file <- read.table("fastq-run-info.tsv", sep = "\t", header = T)
sra_file <- read.table("SraRunTable.csv", sep = ",", header = T) %>%
  select(-c("Assay.Type", "AvgSpotLen", "broad.scale_environmental_context", "Bytes",
            "Center.Name", "Collection_Date", "Consent", "DATASTORE.filetype",
            "DATASTORE.provider", "DATASTORE.region", "ENA.FIRST.PUBLIC..run.",
            "ENA_first_public", "ENA.LAST.UPDATE..run.", "ENA.STATUS.ID..exp.",
            "ENA.STATUS.ID..run.", "environmental_medium", "geo_loc_name_country", 
            "geo_loc_name_country_continent", "geographic_location_.country_and.or_sea.",
            "INSDC_center_name", "INSDC_status", "Library.Name","LibraryLayout",
            "LibrarySelection", "LibrarySource", "local_environmental_context",
            "Organism", "project_name", "ReleaseDate", "create_date", "version",
            "Scientific_Name"))
mdat_chatz <- sra_file %>%
  mutate(run_accession = Run,
         subject_ID = host_subject_id,
         infant = grepl("B", host_subject_id),
         sample_ID = Submitter_Id,
         family_ID = host_subject_id,
         lifestyle = "industrialized",
         study = "chatzigiannidou_2025",
         country = "SWEDEN",
         age = case_when(timepoint == "12 months" ~ 365,
                         timepoint == "18 months" ~ 540,
                         timepoint == "2 months" ~ 60,
                         timepoint == "24 months" ~ 730,
                         timepoint == "6 months" ~ 180,
                         timepoint == "60 months" ~ 1800,
                         timepoint == "3 weeks" ~ 21,
                         timepoint == "3-6 days" ~ 5),
         region = "SWEDEN",
         .keep = "unused")
mdat_chatz_infants <- mdat_chatz %>% filter(infant)
mdat_chatz_mothers <- mdat_chatz %>% filter(!infant)
write.csv(mdat_chatz_infants, "metadata_chatzigiannidou_2025_healthy.csv")
write.csv(mdat_chatz_mothers, "metadata_chatzigiannidou_2025_mothers.csv")

# move good files from shotgun folder
final_samples <- mdat_chatz_infants$run_accession
final_files <- c(paste0(final_samples, "_1.fastq.gz"),
                 paste0(final_samples, "_2.fastq.gz"),
                 paste0(final_samples, ".fastq.gz"))
final_files <- final_files[file.exists(final_files)]
file.move(final_files, "fastq_files/")

final_samples_m <- mdat_chatz_mothers$run_accession
final_files_m <- c(paste0(final_samples_m, "_1.fastq.gz"),
                 paste0(final_samples_m, "_2.fastq.gz"),
                 paste0(final_samples_m, ".fastq.gz"))
final_files_m <- final_files_m[file.exists(final_files_m)]
file.move(final_files_m, "mother_fastq/")

# clean up
rm("ena_file", "final_files", "final_samples", "sra_file", mdat_chatz, 
   mdat_chatz_infants, mdat_chatz_mothers, final_files_m, final_samples_m)



# Metadata/studies summary #####################################################
library(ggplot2)
library(readODS)
library(data.table)
setwd(maindir)
studies_metadata <- readODS::read_ods("/fast/AG_Forslund/rob/Documents/interesting_studies.ods", sheet = "metadata_available") %>%
   tibble::column_to_rownames(var = "metadata_colnames.")
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
   select(n_individuals_final, n_samples, N_ASVs)
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


rm(studies_metadata, studies_metadata_long, maindir, rownames_ordered, asv_data,
   colnames_ordered, add_table, loss_table, loss_table_long, c, f, raw_counts, stat_files)


# check raman, subr, vatanen, gehrig and raman samples:
setwd(paste0(maindir, "16S"))
subr_samples <- read.csv("subramanian_2014/metadata_subramanian_2014_healthy.csv")
raman_samples <- read.csv("raman_2019/metadata_ram_2019_healthy.csv")
gehrig_samples <- read.csv("gehrig_2019/metadata_gehrig_2019_healthy.csv")
gehrig_shotgun_samples <- read.csv("../shotgun/gehrig_2019/metadata_gehrig_2019_healthy.csv")

unique(gehrig_shotgun_samples$subject_ID)[!unique(gehrig_shotgun_samples$subject_ID) %in% raman_samples$subject_ID]
unique(gehrig_shotgun_samples$subject_ID)[!unique(gehrig_shotgun_samples$subject_ID) %in% gehrig_samples$subject_ID]
unique(gehrig_shotgun_samples$subject_ID)[!unique(gehrig_shotgun_samples$subject_ID) %in% subr_samples$subject_ID]

unique(gehrig_samples$subject_ID) %in% raman_samples$subject_ID

