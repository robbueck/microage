rm(list = ls())
maindir <- "/fast/AG_Forslund/rob/studies/16S/hickman_2024/"
setwd(maindir)
library(dada2)
library(tidyverse)
library("reshape2")
library("plyr")
library("dplyr")
library("ggplot2")
library("phangorn")
library("picante")
library("btools")
library("microbiome")
library("gridExtra")


# load additional custom functions
source("/fast/AG_Forslund/rob/mm_index/R_scripts/dada_2_functions.R")

# switch parts of the pipeline on/off
filter_step <- F
error_and_asv_step <- F
remove_chimeras_step <- F
id_taxa_step <- F
tree_step <- F
phyloseq_step <- T
ordination_step <- F
rarecurve_step <- T



# load metadata
metadata_hick <- read.table("metadata_hickman_2024_healthy.csv", sep = ",", header = T)

# remove samples with too high adapter content:
fastqc_forward <- read.table("/fast/AG_Forslund/rob/studies/16S/hickman_2024/fastqc/decontam/forward/multiqc_data/multiqc_fastqc.txt",
                             header = T, sep = "\t") %>%
  mutate(sample_ID = gsub("_dec_1", "", Sample))
fastqc_reverse <- read.table("/fast/AG_Forslund/rob/studies/16S/hickman_2024/fastqc/decontam/reverse/multiqc_data/multiqc_fastqc.txt",
                             header = T, sep = "\t") %>%
  mutate(sample_ID = gsub("_dec_2", "", Sample))
fastqc_merge <- left_join(fastqc_forward, fastqc_reverse, by = "sample_ID")
good_files <- fastqc_merge %>%
  filter(adapter_content.x == "pass", adapter_content.y == "pass") %>%
  pull(sample_ID)

# remove samples with low read counts
read_counts_raw <- read.table("read_counts_raw.txt", header = F,
                              col.names = c("sample_ID", "read_counts_raw"), row.names = 1)
read_counts_raw <- read_counts_raw %>% 
  rownames_to_column() %>%
  filter(rowname %in% good_files, read_counts_raw > 2000)

metadata_hick <- metadata_hick %>%
  filter(run_accession %in% read_counts_raw$rowname)

# filter, learn errors and denoise
print(paste("analyzing", length(metadata_hick$run_accession), sep = " "))
print("Filter")
if(filter_step) {
  out <- filtering_fastq(sample_names = metadata_hick$run_accession,
                     trimLeft = 40,
                     truncLen=245,
                     maxN=0,
                     maxEE=c(2,2),
                     truncQ=2,
                     cores = 20)
  saveRDS(out, "dada2/hickman_out_table.rds")
} else {
  out <- readRDS("dada2/hickman_out_table.rds")
}

if(error_and_asv_step){
  total_bases <- sum(out[,2]) * 205 # because error  estimation is done for forward and reverse
  st.all <- error_and_asv(sample_names = metadata_hick$run_accession, nbases = 0.25* total_bases)
  saveRDS(st.all, "dada2/hickman_seqtab.rds")
} else {
  st.all <- readRDS("dada2/hickman_seqtab.rds")
}
  


# how much got filtered out?
head(out)
cdf_data_filter <- melt(loss_ratio(out[,2], out[,1]),
                        value.name = "filter_ratio")
ggplot(cdf_data_filter, aes(filter_ratio)) +
  stat_ecdf() +
  ylab("samples") +
  xlim(0,1)
ggsave("dada2/cfd_filter_ratios.pdf", device = "pdf")

# assess filter criteria for samples from each library
# fnFs <- sort(file.path(maindir, "fastq_files/decontaminated", paste0(metadata_hick$run_accession, "_dec_1.fastq.gz")))
# fnRs <- sort(file.path(maindir, "fastq_files/decontaminated", paste0(metadata_hick$run_accession, "_dec_2.fastq.gz")))
# filtFs <- sort(file.path(maindir, "fastq_files/filtered", paste0(metadata_hick$run_accession, "_dec_filt_1.fastq.gz")))
# filtRs <- sort(file.path(maindir, "fastq_files/filtered", paste0(metadata_hick$run_accession, "_dec_filt_2.fastq.gz")))
# plotQualityProfile(c(fnFs[1], filtFs[1], fnRs[1], filtRs[1]))
# ggsave("dada2/quality_profile.pdf", device = "pdf")


# Inspect distribution of sequence lengths
cdf_data_seqlen <- melt(nchar(getSequences(st.all)),
                        value.name = "seq_length")
ggplot(cdf_data_seqlen, aes(seq_length)) +
  stat_ecdf() +
  ylab("samples")
ggsave("dada2/cfd_seq_lengths_per_library.pdf", device = "pdf")

# maybe split dataset for chimera removal:

# filter out some low prevalence/abundance ASVs:
# ASV_prev <- colSums(st.all > 9)
# # retain ASVs with more than 9 reads in at least 2 samples
# st.all_filt <- st.all[,ASV_prev > 1]
# sum(st.all_filt)/sum(st.all) # how much is kept after this filtering
# print("Number of cores detected:")
# print(RcppParallel::defaultNumThreads())
# remove chimeras
print("remove chimeras")
if(remove_chimeras_step){
  st.all.nochim <- removeBimeraDenovo(st.all_filt, method="consensus", multithread=20, verbose=TRUE)
  saveRDS(st.all.nochim, "dada2/st_all_nochim.rds")
} else{
  st.all.nochim <- readRDS("dada2/st_all_nochim.rds")
}
dim(st.all.nochim)
sum(st.all.nochim)/sum(st.all)


# where do most reads get lost?
# load read counts for decontamination
read_counts_raw <- read.table("read_counts_raw.txt", header = F,
                              col.names = c("sample_ID", "read_counts_raw"), row.names = 1)
rownames(out) <- gsub("_dec_1.fastq.gz", "", rownames(out))
loss_table <- merge(out, read_counts_raw, by = 0)
head(loss_table)
loss_table <- merge(loss_table,
                    cbind(denoised = rowSums(st.all), removed_chimeras = rowSums(st.all.nochim)),
                    by.x = "Row.names", by.y = 0)

track <- cbind(loss_ratio(loss_table$reads.in, loss_table$read_counts_raw),  # decontamination
               loss_ratio(loss_table$reads.out, loss_table$reads.in),        # filtering
               loss_ratio(loss_table$denoised, loss_table$reads.out),        # denoising
               loss_ratio(loss_table$removed_chimeras, loss_table$denoised), # remove chimeras
               loss_ratio(loss_table$removed_chimeras, loss_table$read_counts_raw))        # total
colnames(track) <- c("decontamination", "filtering", "denoising", "remove_chimeras", "total_loss")
rownames(track) <- loss_table$Row.names
head(track)

cdf_data_step_loss <- melt(track,
                 value.name = "lost_reads")
head(cdf_data_step_loss)
colnames(cdf_data_step_loss)[1:2] <- c("sample", "step")
# violin plot:
ggplot(cdf_data_step_loss, aes(x=step, y = lost_reads, color = step)) + 
  geom_violin() +
  ylim(0,1) +
  geom_boxplot(width=0.1) +
  coord_flip() +
  theme(legend.position="none")
ggsave("dada2/violin_filter_ratios_per_step.pdf", device = "pdf")

# cdf_data_step_loss <- ddply(cdf_data_step_loss, .(step), transform, ecd=ecdf(lost_reads)(lost_reads))
ggplot(cdf_data_step_loss, aes(x=lost_reads, col = step)) +
  stat_ecdf() +
  ylab("samples") +
  xlim(0,1)
cdf_data_step_loss[is.na(cdf_data_step_loss$lost_reads),]
ggsave("dada2/cfd_filter_ratios_per_step.pdf", device = "pdf")

# save loss table to be used for further stats
loss_table <- dplyr::transmute(loss_table, Row.names = Row.names, read_counts_raw = read_counts_raw,
                               decontaminated = reads.in, filtered = reads.out,
                               denoised = denoised, removed_chimeras = removed_chimeras)
write.csv(loss_table, "dada2/read_loss_per_step.csv")


# assign taxonomy
library(DECIPHER)
print("assign taxonomy")
dna <- DNAStringSet(getSequences(st.all.nochim)) # Create a DNAStringSet from the ASVs
if(id_taxa_step) {
  load("/fast/AG_Forslund/rob/references/silva/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
  ids <- IdTaxa(dna, trainingSet, processors=10, verbose=FALSE)
  saveRDS(ids, "dada2/taxa_ids.rds")
} else {
  ids <- readRDS("dada2/taxa_ids.rds")
}
taxa <- taxa_to_matrix(ids)

if (tree_step) {
  tree <- build_tree(st.all.nochim)
  saveRDS(tree, "dada2/tree.rds")
} else {
  tree <- readRDS("dada2/tree.rds")
}

# handoff to phyloseq
library(phyloseq)
library(Biostrings)
theme_set(theme_bw())
if(phyloseq_step) {
  metadata_hick <- read.table("metadata_hickman_2024_healthy.csv", sep = ",", header = T)
  #add metadata
  head(metadata_hick)
  rownames(metadata_hick) <- metadata_hick$run_accession
  all(rownames(metadata_hick) %in% rownames(st.all.nochim))
  all(rownames(st.all.nochim) %in% rownames(metadata_hick))
  metadata_hick <- metadata_hick[match(rownames(st.all.nochim), rownames(metadata_hick)),]
  all(rownames(metadata_hick) == rownames(st.all.nochim))
  metadata_hick$lifestyle <- "westernized"
  ps <- phyloseq(otu_table(st.all.nochim, taxa_are_rows=FALSE), 
                 sample_data(metadata_hick), 
                 tax_table(taxa), tree)
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, dna)
  taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
  saveRDS(ps, file = "/fast/AG_Forslund/rob/studies/16S/hickman_2024/dada2/phyloseq_hickman_2024.rds")
} else {
  ps <- readRDS("/fast/AG_Forslund/rob/studies/16S/hickman_2024/dada2/phyloseq_hickman_2024.rds")
}

# compare a-diversity between different groups
# diversity_plots(ps_object = ps, x = "age", colour = "age")
# ggsave("dada2/a_div_age.pdf", device = "pdf")

# check metadata in PCoA
# remove zero count samples
ps <- prune_samples(sample_sums(ps) >= 10, ps)
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
if(ordination_step) {
  dist.bray <- phyloseq::distance(ps.prop, method = "bray")
  # ordinate
  pcoa_bray <- pcoa(dist.bray)
  save(dist.bray, pcoa_bray, file = "dada2/ordinations.RData")
} else {
  load("dada2/ordinations.RData")
}

rownames(metadata_hick) <- metadata_hick$run_accession
cor(pcoa_bray$vectors,
    metadata_hick[rownames(pcoa_bray$vectors),c("age")], use = "complete.obs") %>% head


#PCoA
# bray-curtis
p1 <- plot_ordination(ps.prop, pcoa_bray, color="age", title="Bray NMDS") +
  scale_color_gradientn(name = "age", colors = topo.colors(7))
ggsave("dada2/PCoA_age.pdf", device = "pdf", plot = p1)

# rarefaction filtering of samples:
ps@sam_data$sample_sum <- sample_sums(ps)
ps@sam_data$study <- "hickman_2024"
ps_genus <- aggregate_taxa(ps, level = "genus")
otu_matrix <- t(otu_table(ps_genus))
class(otu_matrix) <- "matrix"

if(rarecurve_step) {
  raredat <- rarecurve(otu_matrix, step = 200, tidy = T)
  raredat <- left_join(raredat, sample_data(ps), by = c("Site" = "run_accession")) %>%
    select(Site, Sample, Species, X, subject_ID, sample_ID, age,
           lifestyle, sample_sum)
  saveRDS(raredat, "/fast/AG_Forslund/rob/studies/16S/hickman_2024/dada2/raredat.rds")
} else {
  raredat <- readRDS("/fast/AG_Forslund/rob/studies/16S/hickman_2024/dada2/raredat.rds")
}
n_reads <- data.frame(n_reads = sample_sums(ps_genus), 
                      names = sample_names(ps_genus))


cutoff <- 1000
rareplot <- ggplot(raredat) +
  theme_bw() +
  scale_color_discrete(guide = "none") +
  geom_line(aes(x = Sample, y = Species, group = Site, color = subject_ID)) +
  xlim(0, 25000) +
  geom_vline(xintercept = cutoff) +
  # ylim(0,500) +
  labs(title="Rarefaction curves") + xlab("Sequenced Reads") + ylab('ASVs Detected')

n_reads_plot <- ggplot(n_reads, aes(x=n_reads)) +
  theme_bw() +
  geom_freqpoly(binwidth = 250) +
  xlim(0, 25000) +
  geom_vline(xintercept = cutoff) +
  scale_y_log10()+
  ylab("N_samples")
grid.arrange(rareplot, n_reads_plot, heights = c(1, 0.3))
ps_genus_filter <- ps_genus %>% subset_samples(sample_sum >= cutoff)


# predict age ######################
library(caret)
source("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/regression_functions.R")

load("/fast/AG_Forslund/rob/mm_index/R_scripts/regression_models/data/final_no_ls_genus.RData")
final_model_genus_no_ls$rf1
for_caret_list_test <- create_caret_df(ps_object = ps_genus_filter, transformation = "compositional",
                                       additional_cols = c("Observed", "Shannon"),
                                       only_multiple_samples = F,
                                       filter_features = final_model_genus_no_ls$rf1$finalModel$xNames)
predictions <- cbind(for_caret_list_test$metadata,
                     predict(final_model_genus_no_ls, newdata = for_caret_list_test$features)) %>%
  mutate(age_cat = as.factor(age))
ggplot(predictions, aes(x = age, y = rf1)) +
  geom_point() +
  geom_smooth(method = "loess")
performance_lm <- lm(rank(rf1) ~ age, data = predictions)
summary(performance_lm)$r.squared
# check MAZ-score:
maz <- predictions %>%
  mutate(age_cat = as.factor(age)) %>% 
  group_by(age_cat) %>%
  dplyr::summarise(group_median = median(rf1),
                   group_sd = sd(rf1)) %>%
  ungroup() %>%
  left_join(., predictions, by = "age_cat") %>%
  mutate(MAZ = (rf1 - group_median) / group_sd) %>%
  ungroup() %>%
  mutate(birthmode_detail = case_when(grepl("C-Section", delivery_mode) ~ "c",
                                      grepl("antibiotics", delivery_mode) ~ "v_abx",
                                      grepl("vaginal", delivery_mode) ~ "v"))


# check connections with MAZ:
ggplot(maz, aes(x = birthmode_detail, y = MAZ)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.3) +
  facet_wrap(~age_cat)

ggplot(maz, aes(x = sex, y = MAZ)) +
  geom_boxplot()  +
  # geom_jitter() +
  facet_wrap(~age_cat)

(res <- maz %>%
  rstatix::group_by(age_cat) %>%
  rstatix::wilcox_test(MAZ ~ birthmode_detail))  

# export to be used with picrust:
writeXStringSet(ps@refseq, "dada2/ASVs.fasta", format = "fasta")
otu_table_out <- t(as.data.frame(ps@otu_table))
write.table(data.frame(OTU = rownames(otu_table_out), otu_table_out), file = "dada2/otu_table.tsv", sep = "\t", row.names = F, quote = F)
