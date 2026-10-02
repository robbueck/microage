rm(list = ls())
maindir <- "/fast/AG_Forslund/rob/studies/16S/vatanen_2018"
setwd(maindir)
library(dada2)
library("reshape2")
library("plyr")
library("dplyr")
library("ggplot2")
library("phangorn")
library("picante")
library("btools")
library("microbiome")
library("umap")
library("gridExtra")
library("Rtsne")
library("RColorBrewer")


# load additional custom functions
source("/fast/AG_Forslund/rob/mm_index/R_scripts/dada_2_functions.R")

# switch parts of the pipeline on/off
filter_step <- F
error_and_asv_step <- F
remove_chimeras_step <- F
id_taxa_step <- F
tree_step <- F
phyloseq_step <- T


# load metadata
metadata_vat <- read.table("metadata_vat_2018_healthy.csv", sep = ",", header = T)
library_names <- unique(metadata_vat$cohort)
if(any(is.na(library_names))){
  stop("Library names contain NAs")
}
out <- vector("list", length(library_names))
seqtab <- vector("list", length(library_names))
names(out) <- library_names
names(seqtab) <- library_names

# filter, learn errors and denoise seperately

for(library in library_names) {
  library_samples <- metadata_vat$run_accession[metadata_vat$cohort == library]
  print(paste("analyzing", length(library_samples), "samples in",
              library, which(library_names == library), "of",
              length(library_names), "libraries",
              sep = " "))
  print("Filter")
  if(filter_step) {
    out[[library]] <- filtering_fastq(sample_names = library_samples,
                                  trimLeft = 10,
                                  truncLen=150,
                                  maxN=0,
                                  maxEE=c(2,2),
                                  truncQ=2)
    saveRDS(out[[library]], paste0(maindir, "/dada2/lib_", library, "_out_table.rds")) # save for each run separately
  } else {
    out[[library]] <- readRDS(paste0(maindir, "/dada2/lib_", library, "_out_table.rds"))
  }
}

for(library in library_names) {
  print("error estimation and asv step")
  library_samples <- metadata_vat$run_accession[metadata_vat$cohort == library]
  print(paste("analyzing", length(library_samples), "samples in",
              library, which(library_names == library), "of",
              length(library_names), "libraries",
              sep = " "))
  if(error_and_asv_step){
    total_bases <- sum(out[[library]][,2]) * 140 # because error estimation is done for forward and reverse
    seqtab[[library]] <- error_and_asv(sample_names = library_samples, libname = library, nbases = 0.25*total_bases)
    saveRDS(seqtab[[library]], paste0(maindir, "/dada2/lib_", library, "_seqtab.rds")) # save for each run separately
  } else {
    seqtab[[library]] <- readRDS(paste0(maindir, "/dada2/lib_", library, "_seqtab.rds"))
  }
}

st.all <- mergeSequenceTables(tables = seqtab)

# how much got filtered out in each library?
head(out[[1]])
out[["combined"]] <-  do.call("rbind", out)
cdf_data_filter <- melt(lapply(out, function(x){loss_ratio(x[,2], x[,1])}),
                 value.name = "filter_ratio")
colnames(cdf_data_filter)[2] <- "Library"
# violin plot:
sample_size = cdf_data_filter %>% group_by(Library) %>% dplyr::summarize(num=n())
cdf_data_filter <- cdf_data_filter %>% 
  left_join(sample_size) %>%
  mutate(Library = paste0(Library, "\n", "n=", num))
ggplot(cdf_data_filter, aes(x=Library, y = filter_ratio, color = Library)) + 
  geom_violin() +
  ylim(0,1) +
  geom_boxplot(width=0.1) +
  coord_flip() +
  theme(legend.position="none")
ggsave("dada2/violin_filter_ratios_per_library.pdf", device = "pdf")
# cdf_data_filter <- ddply(cdf_data_filter, .(Library), transform, ecd=ecdf(filter_ratio)(filter_ratio))
ggplot(cdf_data_filter, aes(x=filter_ratio, col = Library)) + 
  stat_ecdf() +
  ylab("samples") +
  xlim(0,1)
ggsave("dada2/cfd_filter_ratios_per_library.pdf", device = "pdf")


# assess filter criteria for samples from each library
library_samples <- metadata_vat$run_accession[metadata_vat$cohort == "karelia"]
fnFs <- sort(file.path(maindir, "fastq_files/decontaminated", paste0(library_samples, "_dec_1.fastq.gz")))
fnRs <- sort(file.path(maindir, "fastq_files/decontaminated", paste0(library_samples, "_dec_2.fastq.gz")))
filtFs <- sort(file.path(maindir, "fastq_files/filtered", paste0(library_samples, "_dec_filt_1.fastq.gz")))
filtRs <- sort(file.path(maindir, "fastq_files/filtered", paste0(library_samples, "_dec_filt_2.fastq.gz")))
plotQualityProfile(c(fnFs[1], filtFs[1], fnRs[1], filtRs[1]))
ggsave("dada2/quality_profile.pdf", device = "pdf")

# Inspect distribution of sequence lengths
cdf_data_seqlen <- melt(lapply(seqtab, function(x){nchar(getSequences(x))}),
                 value.name = "seq_length")
colnames(cdf_data_seqlen)[2] <- "Library"
# violin plot:
sample_size = cdf_data_seqlen %>% group_by(Library) %>% dplyr::summarize(num=n())
cdf_data_seqlen <- cdf_data_seqlen %>% 
  left_join(sample_size) %>%
  mutate(Library = paste0(Library, "\n", "n=", num))
ggplot(cdf_data_seqlen, aes(x=Library, y = seq_length, color = Library)) + 
  geom_violin() +
  geom_boxplot(width=0.1) +
  coord_flip() +
  theme(legend.position="none")
ggsave("dada2/violin_seq_lengths_per_library.pdf", device = "pdf")

# cdf_data_seqlen <- ddply(cdf_data_seqlen, .(Library), transform, ecd=ecdf(seq_length)(seq_length))
ggplot(cdf_data_seqlen, aes(x=seq_length, col = Library)) + 
  stat_ecdf() +
  ylab("samples") 
ggsave("dada2/cfd_seq_lengths_per_library.pdf", device = "pdf")



# remove chimeras
print("remove chimeras")
if(remove_chimeras_step){
  st.all.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=10, verbose=TRUE)
  saveRDS(st.all.nochim, "dada2/st_all_nochim.rds")
} else {
  st.all.nochim <- readRDS("dada2/st_all_nochim.rds")
}
dim(st.all.nochim)
sum(st.all.nochim)/sum(st.all)


# where do most reads get lost?
# load read counts for decontamination
read_counts_raw <- read.table("read_counts_raw.txt", header = F,
                              col.names = c("sample_ID", "read_counts_raw"), row.names = 1)
rownames(out$combined) <- gsub("_dec_1.fastq.gz", "", rownames(out$combined))
loss_table <- merge(out$combined, read_counts_raw, by = 0)
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
  ylab("samples")
cdf_data_step_loss[is.na(cdf_data_step_loss$lost_reads),]
ggsave("dada2/cfd_filter_ratios_per_step.pdf", device = "pdf")

# save loss table to be used for further stats
loss_table <- dplyr::transmute(loss_table, Row.names = Row.names, read_counts_raw = read_counts_raw,
                               decontaminated = reads.in, filtered = reads.out,
                               denoised = denoised, removed_chimeras = removed_chimeras)
write.csv(loss_table, "dada2/read_loss_per_step.csv")


# assign taxonomy
library(DECIPHER)
print("id taxa")
dna <- DNAStringSet(getSequences(st.all.nochim)) # Create a DNAStringSet from the ASVs
start.time_silva <- Sys.time()
if(id_taxa_step) {
  load("/fast/AG_Forslund/rob/references/silva/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
  ids <- IdTaxa(dna, trainingSet, processors=10, verbose=FALSE)
  saveRDS(ids, "dada2/taxa_ids.rds")
} else {
  ids <- readRDS("dada2/taxa_ids.rds")
}
end.time_silva <- Sys.time()
time_silva <- end.time_silva - start.time_silva

taxa <- taxa_to_matrix(ids)

start.time_gtdb <- Sys.time()
if(id_taxa_step) {
  load("/fast/AG_Forslund/rob/references/GTDB_r207-mod_April2022/GTDB_r207-mod_April2022.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
  ids_gtdb <- IdTaxa(dna, trainingSet, processors=10, verbose=FALSE)
  saveRDS(ids_gtdb, "dada2/taxa_ids_GTDB.rds")
} else {
  ids_gtdb <- readRDS("dada2/taxa_ids_GTDB.rds")
}
end.time_gtdb <- Sys.time()
time_gtdb <- end.time_gtdb - start.time_gtdb

taxa_gtdb <- taxa_to_matrix(ids_gtdb)

taxdb_table <- t(data.frame(silva = (!is.na(taxa)) %>% colSums() / nrow(taxa),
                            gtdb = (!is.na(taxa_gtdb)) %>% colSums() / nrow(taxa_gtdb)))



print("tree step")
if (tree_step) {
  tree <- build_tree(st.all.nochim)
  saveRDS(tree, "dada2/tree.rds")
} else {
  tree <- readRDS("dada2/tree.rds")
}


# handoff to phyloseq
library(phyloseq)
library(Biostrings)
library(ggplot2)
theme_set(theme_bw())
if(phyloseq_step) {
  metadata_vat <- read.table("metadata_vat_2018_healthy.csv", sep = ",", header = T)
  #add metadata
  head(metadata_vat)
  rownames(metadata_vat) <- metadata_vat$run_accession
  all(rownames(metadata_vat) %in% rownames(st.all.nochim))
  all(rownames(st.all.nochim) %in% rownames(metadata_vat))
  metadata_vat <- metadata_vat[match(rownames(st.all.nochim), rownames(metadata_vat)),]
  all(rownames(metadata_vat) == rownames(st.all.nochim))
  # number of ASVs unassigned in whole study
  metadata_vat$n_unassigned_asvs_study <-   sum(grepl("unclassified", taxa[,6]))
  metadata_vat$n_total_asvs_study <- nrow(taxa)
  # number of ASVs not assigned to genus level per sample
  matching_taxa <- rownames(taxa)[grepl("unclassified", taxa[,6])]
  unassigned_taxa_present_per_sample <- apply(st.all.nochim[,matching_taxa], 1, function(x) sum(x > 0))
  total_taxa_present_per_sample <- apply(st.all.nochim, 1, function(x) sum(x > 0))
  metadata$n_unassigned_asvs_sample <- unassigned_taxa_present_per_sample
  metadata$n_total_asvs_sample <- total_taxa_present_per_sample

  ps <- phyloseq(otu_table(st.all.nochim, taxa_are_rows=FALSE), 
                 sample_data(metadata_vat), 
                 tax_table(taxa), tree)
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, dna)
  taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
  # rarefaction for alpha diversity
  if(min(sample_sums(ps)) < 2000) {
    ps_low_counts <- ps %>% prune_samples(samples = (sample_sums(.) < 2000))
    ps_raref <- ps %>% rarefy_even_depth(., sample.size = 2000, trimOTUs = F) %>%
      merge_phyloseq(.,ps_low_counts)
  } else {
    ps_raref <- ps %>% rarefy_even_depth(., sample.size = 2000)
  }
  # remove low abundant taxa:
  remove_taxa_asv <- ps_raref %>%
    prevalence(., detection = 10)
  remove_taxa_asv <- names(remove_taxa_asv[remove_taxa_asv < 5e-5])
  ps_raref <- ps_raref %>%
    merge_taxa2_fixed(.,taxa = remove_taxa_asv, name = "Other")
  # add alpha diversity on ASV level
  alpha.div_gen <- microbiome::alpha(ps_raref, index = c("Observed", "Shannon"))
  alpha.div_gen <- alpha.div_gen[sample_names(ps),]
  if (all(sample_names(ps) == rownames(alpha.div_gen))) {
    ps@sam_data$shannon_asv <- alpha.div_gen$observed
    ps@sam_data$richness_asv <- alpha.div_gen$diversity_shannon
  } else {stop(" Names don't match")}
  
  # save object
  saveRDS(ps, file = "/fast/AG_Forslund/rob/studies/16S/vatanen_2018/dada2/phyloseq_vatanen_2018.rds")
} else {
  ps <- readRDS("/fast/AG_Forslund/rob/studies/16S/vatanen_2018/dada2/phyloseq_vatanen_2018.rds")
}